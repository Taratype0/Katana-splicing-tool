from __future__ import annotations

import gzip
import re
import json
from pathlib import Path

import pandas as pd
import numpy as np

from katana_next.services.project_service import ProjectContext
from katana_next.science.gene_filter import apply_gene_filter


ATTRIBUTE_RE = re.compile(r'(\S+)\s+"([^"]*)"')


def transcript_gene_map(gtf: Path, cache_path: Path) -> pd.DataFrame:
    from katana_next.infrastructure.fingerprints import fingerprint_file
    identity = fingerprint_file(gtf).model_dump(mode="json")
    marker = cache_path.with_suffix(".source.json")
    if cache_path.is_file() and marker.is_file() and json.loads(marker.read_text(encoding="utf-8")) == identity:
        return pd.read_parquet(cache_path)
    opener = gzip.open if gtf.suffix.lower() == ".gz" else open
    rows: dict[str, tuple[str, str]] = {}
    with opener(gtf, "rt", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] not in {"transcript", "exon"}:
                continue
            attributes = dict(ATTRIBUTE_RE.findall(parts[8]))
            transcript = attributes.get("transcript_id", "").split(".", 1)[0]
            gene = attributes.get("gene_id", "").split(".", 1)[0]
            symbol = attributes.get("gene_name", gene)
            if transcript and transcript not in rows:
                rows[transcript] = (gene, symbol)
    frame = pd.DataFrame(
        [
            {"transcript_id": transcript, "gene_id": gene, "gene_symbol": symbol}
            for transcript, (gene, symbol) in rows.items()
        ]
    )
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(cache_path, index=False)
    marker.write_text(json.dumps(identity), encoding="utf-8")
    return frame


def load_quant_support(context: ProjectContext, cache_path: Path, *, sample_ids: set[str] | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    gtf = context.manifest.annotation.gtf
    if gtf is None or not gtf.is_file():
        raise ValueError("Isoform support requires the project GTF.")
    mapping = transcript_gene_map(gtf, cache_path)
    frames = []
    reference_ids = None
    for sample in context.manifest.samples:
        if sample_ids is not None and sample.sample_id not in sample_ids:
            continue
        quant = sample.assets.get("quant_sf")
        if quant is None or not quant.is_file():
            raise ValueError(f"Missing quant.sf for required sample {sample.sample_id}; unmeasured is not zero.")
        source = pd.read_csv(quant, sep="\t", low_memory=False)
        required = {"Name", "TPM", "NumReads"}
        if not required.issubset(source.columns):
            raise ValueError(f"Invalid quant.sf schema for {sample.sample_id}: missing {sorted(required - set(source.columns))}")
        if source.empty or source.Name.isna().any():
            raise ValueError(f"Empty quant.sf or missing transcript identity: {sample.sample_id}")
        ids = source.Name.astype(str).str.split("|").str[0].str.replace(r"\.\d+$", "", regex=True)
        if ids.eq("").any() or ids.duplicated().any():
            raise ValueError(f"Duplicate/empty normalized transcript IDs: {sample.sample_id}")
        for column in ("TPM", "NumReads"):
            source[column] = pd.to_numeric(source[column], errors="coerce")
            if not np.isfinite(source[column]).all() or source[column].lt(0).any():
                raise ValueError(f"Invalid non-finite/negative {column} in quant.sf: {sample.sample_id}")
        if reference_ids is not None and set(ids) != reference_ids:
            raise ValueError(f"quant.sf transcript catalogs differ for {sample.sample_id}; use a common Salmon reference, not zero filling.")
        reference_ids = set(ids)
        frame = pd.DataFrame(
            {
                "transcript_id": ids,
                "tpm": source["TPM"],
                "estimated_reads": source["NumReads"],
            }
        ).merge(mapping, on="transcript_id", how="left")
        frame["sample_id"] = sample.sample_id
        frame["annotation_status"] = np.where(frame.gene_id.notna(), "mapped", "not in declared GTF")
        frame["condition"] = sample.condition
        frame["row_id"] = sample.sample_id + "::" + frame["transcript_id"]
        frame = apply_gene_filter(context, frame)
        frames.append(frame)
    if not frames:
        raise ValueError("No valid quant.sf sample assets are available.")
    transcript = pd.concat(frames, ignore_index=True)
    gene = (
        transcript[transcript.gene_id.notna()].groupby(["sample_id", "condition", "gene_id", "gene_symbol"], dropna=False)
        .agg(total_tpm=("tpm", "sum"), transcript_count=("transcript_id", "nunique"))
        .reset_index()
    )
    gene["row_id"] = gene["sample_id"] + "::" + gene["gene_id"].fillna("").astype(str)
    return transcript, gene


def isoform_contrast_summary(context: ProjectContext, transcript: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for contrast in context.manifest.contrasts:
        numerator = transcript[transcript["sample_id"].isin(contrast.numerator.sample_ids)]
        denominator = transcript[transcript["sample_id"].isin(contrast.denominator.sample_ids)]
        available = set(transcript.sample_id)
        required = set((*contrast.numerator.sample_ids, *contrast.denominator.sample_ids))
        if not required.issubset(available):
            continue
        numerator_mean = numerator.groupby(["transcript_id", "gene_id", "gene_symbol"], dropna=False)["tpm"].mean().rename("numerator_mean_tpm")
        denominator_mean = denominator.groupby(["transcript_id", "gene_id", "gene_symbol"], dropna=False)["tpm"].mean().rename("denominator_mean_tpm")
        joined = pd.concat([numerator_mean, denominator_mean], axis=1).reset_index()
        joined["delta_tpm"] = joined["numerator_mean_tpm"] - joined["denominator_mean_tpm"]
        joined["comparison_id"] = contrast.contrast_id
        joined["numerator"] = contrast.numerator.label
        joined["denominator"] = contrast.denominator.label
        joined["row_id"] = contrast.contrast_id + "::" + joined["transcript_id"]
        rows.append(joined)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def dominant_isoform_switches(context: ProjectContext, transcript: pd.DataFrame) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    gene_columns = ["gene_id", "gene_symbol"]
    for contrast in context.manifest.contrasts:
        if not set((*contrast.numerator.sample_ids, *contrast.denominator.sample_ids)).issubset(set(transcript.sample_id)):
            continue
        selected = transcript[
            transcript["sample_id"].isin(
                (*contrast.numerator.sample_ids, *contrast.denominator.sample_ids)
            )
        ].copy()
        selected = selected[selected.gene_id.notna()]
        if selected.empty:
            continue
        selected["contrast_arm"] = ""
        selected.loc[
            selected["sample_id"].isin(contrast.numerator.sample_ids), "contrast_arm"
        ] = "numerator"
        selected.loc[
            selected["sample_id"].isin(contrast.denominator.sample_ids), "contrast_arm"
        ] = "denominator"
        arm_means = (
            selected.groupby(
                ["contrast_arm", "gene_id", "gene_symbol", "transcript_id"],
                dropna=False,
            )["tpm"]
            .mean()
            .reset_index()
        )
        top_by_arm: dict[str, pd.DataFrame] = {}
        for arm in ("numerator", "denominator"):
            top = (
                arm_means.loc[arm_means["contrast_arm"].eq(arm)]
                .sort_values(
                    [*gene_columns, "tpm", "transcript_id"],
                    ascending=[True, True, False, True],
                    na_position="last",
                )
                .drop_duplicates(gene_columns, keep="first")
                .rename(
                    columns={
                        "transcript_id": f"{arm}_dominant_transcript",
                        "tpm": f"{arm}_dominant_tpm",
                    }
                )
            )
            top_by_arm[arm] = top[
                [
                    *gene_columns,
                    f"{arm}_dominant_transcript",
                    f"{arm}_dominant_tpm",
                ]
            ]
        joined = top_by_arm["numerator"].merge(
            top_by_arm["denominator"], on=gene_columns, how="inner"
        )
        joined = joined.loc[
            joined["numerator_dominant_transcript"].ne(
                joined["denominator_dominant_transcript"]
            ) & joined["numerator_dominant_tpm"].gt(0) & joined["denominator_dominant_tpm"].gt(0)
        ].copy()
        if joined.empty:
            continue
        joined.insert(0, "comparison_id", contrast.contrast_id)
        joined.insert(
            0,
            "row_id",
            contrast.contrast_id + "::" + joined["gene_id"].fillna("").astype(str),
        )
        joined["support_semantics"] = (
            "descriptive mean-TPM dominant transcript change; "
            "DTU significance is reported separately"
        )
        frames.append(joined)
    if frames:
        return pd.concat(frames, ignore_index=True)
    return pd.DataFrame(
        columns=[
            "row_id",
            "comparison_id",
            *gene_columns,
            "numerator_dominant_transcript",
            "numerator_dominant_tpm",
            "denominator_dominant_transcript",
            "denominator_dominant_tpm",
            "support_semantics",
        ]
    )
