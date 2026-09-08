from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from katana_next.science.statistics import beta_binomial_two_group, bh_adjust
from katana_next.services.project_service import ProjectContext


REQUIRED_CELL_METADATA = {"cell_id", "sample_id", "donor_id", "cell_type", "condition"}
REQUIRED_JUNCTION_COUNTS = {"cell_id", "event_id", "inclusion", "exclusion"}


def load_single_cell_metadata(context: ProjectContext) -> pd.DataFrame:
    asset = next(
        (
            item
            for item in context.manifest.assets
            if item.kind == "single_cell_metadata" and item.path.is_file()
        ),
        None,
    )
    if asset is None:
        raise ValueError("Single-cell analysis requires a single_cell_metadata asset.")
    frame = pd.read_csv(asset.path, sep=None, engine="python")
    missing = sorted(REQUIRED_CELL_METADATA - set(frame.columns))
    if missing:
        raise ValueError(f"Single-cell metadata is missing columns: {missing}")
    frame = frame.copy()
    for column in REQUIRED_CELL_METADATA:
        frame[column] = frame[column].astype(str).str.strip()
    if frame["cell_id"].duplicated().any():
        duplicates = frame.loc[frame["cell_id"].duplicated(), "cell_id"].head(5).tolist()
        raise ValueError(f"Single-cell metadata contains duplicate cell IDs: {duplicates}")
    if "batch" not in frame:
        frame["batch"] = ""
    if "assay" not in frame:
        frame["assay"] = "unspecified"
    frame["batch"] = frame["batch"].fillna("").astype(str)
    frame["assay"] = frame["assay"].fillna("unspecified").astype(str)
    frame["row_id"] = frame["cell_id"]
    return frame


def load_single_cell_counts(context: ProjectContext) -> pd.DataFrame:
    asset = next(
        (
            item
            for item in context.manifest.assets
            if item.kind == "single_cell_junction_counts" and item.path.is_file()
        ),
        None,
    )
    if asset is None:
        raise ValueError(
            "Quantification requires a single_cell_junction_counts asset or a BAM extraction run."
        )
    frame = pd.read_csv(asset.path, sep=None, engine="python")
    missing = sorted(REQUIRED_JUNCTION_COUNTS - set(frame.columns))
    if missing:
        raise ValueError(f"Single-cell junction counts are missing columns: {missing}")
    frame = frame.copy()
    frame["cell_id"] = frame["cell_id"].astype(str).str.strip()
    frame["event_id"] = frame["event_id"].astype(str).str.strip()
    frame["inclusion"] = pd.to_numeric(frame["inclusion"], errors="coerce")
    frame["exclusion"] = pd.to_numeric(frame["exclusion"], errors="coerce")
    invalid = (
        frame[["inclusion", "exclusion"]].isna().any(axis=1)
        | frame["inclusion"].lt(0)
        | frame["exclusion"].lt(0)
    )
    if invalid.any():
        raise ValueError(f"Single-cell junction counts contain {int(invalid.sum())} invalid rows.")
    defaults = {
        "gene_id": "",
        "gene_symbol": "",
        "event_type": "junction_choice",
        "annotated": False,
    }
    for column, default in defaults.items():
        if column not in frame:
            frame[column] = default
    frame["gene_id"] = frame["gene_id"].fillna("").astype(str)
    frame["gene_symbol"] = frame["gene_symbol"].fillna("").astype(str)
    frame["event_type"] = frame["event_type"].fillna("junction_choice").astype(str)
    frame["annotated"] = frame["annotated"].fillna(False).astype(bool)
    frame["total"] = frame["inclusion"] + frame["exclusion"]
    frame["psi_posterior"] = (frame["inclusion"] + 0.5) / (frame["total"] + 1.0)
    frame["row_id"] = frame["cell_id"] + "::" + frame["event_id"]
    return frame


def single_cell_audit(
    metadata: pd.DataFrame,
    counts: pd.DataFrame | None,
    *,
    min_cells_per_donor_cell_type: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    design = (
        metadata.groupby(["condition", "cell_type", "donor_id", "sample_id", "assay"], dropna=False)
        .size()
        .rename("cells")
        .reset_index()
    )
    design["passes_cell_count"] = design["cells"].ge(min_cells_per_donor_cell_type)
    design["row_id"] = design["condition"] + "::" + design["cell_type"] + "::" + design["donor_id"]
    rows = [
        {
            "row_id": "cell_identity",
            "gate": "Cell identity",
            "status": "pass",
            "observed": str(metadata["cell_id"].nunique()),
            "required": "unique cell_id",
            "interpretation": "Every cell maps to one sample and donor.",
        },
        {
            "row_id": "biological_units",
            "gate": "Independent biological units",
            "status": "pass" if metadata["donor_id"].nunique() >= 4 else "warning",
            "observed": str(metadata["donor_id"].nunique()),
            "required": "at least 2 donors per tested condition; more recommended",
            "interpretation": "Cells are not treated as independent biological replicates.",
        },
        {
            "row_id": "assay",
            "gate": "Assay declaration",
            "status": "warning" if metadata["assay"].eq("unspecified").any() else "pass",
            "observed": ", ".join(sorted(metadata["assay"].unique())),
            "required": "droplet or full-length chemistry declaration",
            "interpretation": "Quantification and claim boundaries depend on assay chemistry.",
        },
    ]
    if counts is not None:
        matched = counts["cell_id"].isin(set(metadata["cell_id"]))
        rows.append(
            {
                "row_id": "count_cell_mapping",
                "gate": "Junction-count cell mapping",
                "status": "pass" if bool(matched.all()) else "fail",
                "observed": f"{int(matched.sum())}/{len(matched)} count rows",
                "required": "all count cells present in metadata",
                "interpretation": "Unmapped cells invalidate donor-aware aggregation.",
            }
        )
    return pd.DataFrame(rows), design


def donor_celltype_pseudobulk(
    metadata: pd.DataFrame,
    counts: pd.DataFrame,
    *,
    min_cell_event_reads: float = 1.0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    merged = counts.merge(
        metadata[["cell_id", "sample_id", "donor_id", "cell_type", "condition", "batch", "assay"]],
        on="cell_id",
        how="inner",
        validate="many_to_one",
    )
    merged["detected"] = merged["total"].ge(min_cell_event_reads)
    grouping = [
        "event_id",
        "gene_id",
        "gene_symbol",
        "event_type",
        "annotated",
        "sample_id",
        "donor_id",
        "cell_type",
        "condition",
        "batch",
        "assay",
    ]
    aggregate = (
        merged.groupby(grouping, dropna=False)
        .agg(
            inclusion=("inclusion", "sum"),
            exclusion=("exclusion", "sum"),
            contributing_cells=("cell_id", "nunique"),
            detected_cells=("detected", "sum"),
        )
        .reset_index()
    )
    aggregate["total"] = aggregate["inclusion"] + aggregate["exclusion"]
    aggregate["psi_posterior"] = (aggregate["inclusion"] + 0.5) / (aggregate["total"] + 1.0)
    aggregate["measurement_se"] = np.sqrt(
        aggregate["psi_posterior"] * (1 - aggregate["psi_posterior"])
        / (aggregate["total"] + 1.0)
    )
    aggregate["row_id"] = (
        aggregate["event_id"]
        + "::"
        + aggregate["cell_type"]
        + "::"
        + aggregate["donor_id"]
    )
    coverage = (
        aggregate.groupby(["cell_type", "condition", "donor_id"], dropna=False)
        .agg(
            cells=("contributing_cells", "max"),
            measured_events=("total", lambda values: int(pd.Series(values).gt(0).sum())),
            median_event_reads=("total", "median"),
            median_measurement_se=("measurement_se", "median"),
        )
        .reset_index()
    )
    coverage["row_id"] = coverage["cell_type"] + "::" + coverage["condition"] + "::" + coverage["donor_id"]
    return aggregate, coverage


def differential_pseudobulk(
    pseudobulk: pd.DataFrame,
    *,
    numerator: str,
    denominator: str,
    min_total_reads: float,
    min_donors_per_condition: int,
    fdr_cutoff: float,
    min_abs_dpsi: float,
) -> pd.DataFrame:
    selected = pseudobulk[pseudobulk["condition"].isin([numerator, denominator])].copy()
    rows: list[dict[str, Any]] = []
    for (cell_type, event_id), group in selected.groupby(["cell_type", "event_id"], sort=False):
        group = group[group["total"].ge(min_total_reads)].copy()
        donors = group.groupby("condition")["donor_id"].nunique()
        if donors.get(numerator, 0) < min_donors_per_condition or donors.get(denominator, 0) < min_donors_per_condition:
            continue
        groups = group["condition"].eq(numerator).astype(int)
        fit = beta_binomial_two_group(group["inclusion"], group["total"], groups)
        first = group.iloc[0]
        rows.append(
            {
                "row_id": f"{cell_type}::{event_id}",
                "cell_type": cell_type,
                "event_id": event_id,
                "gene_id": first.get("gene_id", ""),
                "gene_symbol": first.get("gene_symbol", ""),
                "event_type": first.get("event_type", "junction_choice"),
                "numerator": numerator,
                "denominator": denominator,
                "numerator_donors": int(donors.get(numerator, 0)),
                "denominator_donors": int(donors.get(denominator, 0)),
                "numerator_cells": int(group.loc[group["condition"].eq(numerator), "contributing_cells"].sum()),
                "denominator_cells": int(group.loc[group["condition"].eq(denominator), "contributing_cells"].sum()),
                "dpsi": fit["dpsi"],
                "log_odds": fit["log_odds"],
                "overdispersion_rho": fit["rho"],
                "pvalue": fit["pvalue"],
                "model_status": fit["status"],
            }
        )
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    result["fdr"] = result.groupby("cell_type")["pvalue"].transform(
        lambda values: bh_adjust(values)
    )
    result["significant"] = result["fdr"].le(fdr_cutoff) & result["dpsi"].abs().ge(min_abs_dpsi)
    result["direction"] = np.where(result["dpsi"].ge(0), "higher inclusion in numerator", "lower inclusion in numerator")
    return result.sort_values(["significant", "fdr", "dpsi"], ascending=[False, True, False], na_position="last")


def extract_junction_choices_from_bams(
    context: ProjectContext,
    output_path: Path,
    *,
    cell_barcode_tag: str = "CB",
    umi_tag: str = "UB",
    min_mapping_quality: int = 20,
    max_reads_per_bam: int = 0,
) -> pd.DataFrame:
    backend = "pysam"
    try:
        import pysam
    except ImportError:
        pysam = None
        backend = "bamnostic"
        try:
            import bamnostic  # noqa: F401
        except ImportError as exc:
            raise RuntimeError("BAM extraction requires pysam or bamnostic.") from exc

    metadata_cells: set[str] = set()
    try:
        metadata_cells = set(load_single_cell_metadata(context)["cell_id"].astype(str))
    except (ValueError, FileNotFoundError):
        pass

    molecules: set[tuple[str, str, str, int, int, str, str]] = set()
    for sample in context.manifest.samples:
        bam = sample.assets.get("single_cell_bam") or sample.assets.get("bam")
        if bam is None or not bam.is_file():
            continue
        if backend == "pysam":
            handle = pysam.AlignmentFile(str(bam), "rb")
            records = handle.fetch(until_eof=True)
        else:
            import bamnostic

            index_candidates = (Path(str(bam) + ".bai"), bam.with_suffix(".bai"))
            index_path = next((path for path in index_candidates if path.is_file()), None)
            kwargs = {"filepath_index": str(index_path)} if index_path else {}
            handle = bamnostic.AlignmentFile(str(bam), "rb", **kwargs)
            records = iter(handle)
        try:
            for read_index, read in enumerate(records, start=1):
                if max_reads_per_bam and read_index > max_reads_per_bam:
                    break
                if any(bool(getattr(read, flag, False)) for flag in ("is_unmapped", "is_secondary", "is_supplementary")):
                    continue
                if int(getattr(read, "mapping_quality", getattr(read, "mapq", 0)) or 0) < min_mapping_quality:
                    continue
                try:
                    raw_cell = str(read.get_tag(cell_barcode_tag))
                except (KeyError, AttributeError):
                    continue
                prefixed_cell = f"{sample.sample_id}:{raw_cell}"
                cell = raw_cell if raw_cell in metadata_cells else prefixed_cell
                try:
                    umi = str(read.get_tag(umi_tag))
                except (KeyError, AttributeError):
                    umi = str(getattr(read, "query_name", read_index))
                chrom = str(getattr(read, "reference_name", "") or "")
                if not chrom and backend == "pysam":
                    chrom = str(handle.get_reference_name(read.reference_id))
                reference = int(read.reference_start)
                try:
                    strand = str(read.get_tag("XS"))
                except (KeyError, AttributeError):
                    strand = "."
                cigar = list(getattr(read, "cigartuples", None) or getattr(read, "cigar", []) or [])
                for item in cigar:
                    operation, length = _normalize_cigar_item(item)
                    if operation == 3:
                        start, end = reference, reference + int(length)
                        molecules.add((sample.sample_id, cell, chrom, start, end, strand, umi))
                    if operation in {0, 2, 3, 7, 8}:
                        reference += int(length)
        finally:
            handle.close()
    if not molecules:
        raise ValueError("No UMI-supported splice junctions were detected in configured BAM files.")
    raw = pd.DataFrame(
        molecules,
        columns=["sample_id", "cell_id", "chrom", "start", "end", "strand", "molecule"],
    )
    junctions = (
        raw.groupby(["sample_id", "cell_id", "chrom", "start", "end", "strand"])
        .size()
        .rename("junction_count")
        .reset_index()
    )
    choices = junctions[["chrom", "start", "end", "strand"]].drop_duplicates().reset_index(drop=True)
    choices["cluster_id"] = _splice_site_components(choices)
    choices["event_id"] = (
        choices["chrom"] + ":" + choices["start"].astype(str) + "-" + choices["end"].astype(str) + ":" + choices["strand"]
    )
    choices = _annotate_junction_choices(choices, context.manifest.annotation.gtf)
    junctions = junctions.merge(choices, on=["chrom", "start", "end", "strand"], how="left", validate="many_to_one")

    observed_clusters = (
        junctions.groupby(["sample_id", "cell_id", "cluster_id"], as_index=False)["junction_count"]
        .sum()
        .rename(columns={"junction_count": "cluster_total"})
    )
    expanded = observed_clusters.merge(
        choices[
            ["cluster_id", "event_id", "gene_id", "gene_symbol", "annotated"]
        ],
        on="cluster_id",
        how="left",
        validate="many_to_many",
    )
    inclusion = junctions[["sample_id", "cell_id", "cluster_id", "event_id", "junction_count"]]
    expanded = expanded.merge(
        inclusion,
        on=["sample_id", "cell_id", "cluster_id", "event_id"],
        how="left",
        validate="one_to_one",
    )
    expanded["junction_count"] = expanded["junction_count"].fillna(0).astype(int)
    output = pd.DataFrame({
        "cell_id": expanded["cell_id"],
        "event_id": expanded["event_id"],
        "cluster_id": expanded["cluster_id"],
        "gene_id": expanded["gene_id"],
        "gene_symbol": expanded["gene_symbol"],
        "event_type": "junction_choice",
        "annotated": expanded["annotated"],
        "inclusion": expanded["junction_count"],
        "exclusion": expanded["cluster_total"] - expanded["junction_count"],
        "count_backend": backend,
    })
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(output_path, sep="\t", index=False)
    return output


def _annotate_junction_choices(
    choices: pd.DataFrame, gtf_path: Path | None
) -> pd.DataFrame:
    annotated = choices.copy()
    annotated["gene_id"] = ""
    annotated["gene_symbol"] = ""
    annotated["annotated"] = False
    if gtf_path is None or not gtf_path.is_file():
        return annotated

    transcripts: dict[
        tuple[str, str, str], list[tuple[int, int, str, str]]
    ] = {}
    gene_spans: dict[tuple[str, str], list[object]] = {}
    with gtf_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            attributes = _gtf_attributes(fields[8])
            gene_id = attributes.get("gene_id", "")
            transcript_id = attributes.get("transcript_id", "")
            if not gene_id or not transcript_id:
                continue
            chrom, strand = fields[0], fields[6]
            start, end = int(fields[3]), int(fields[4])
            symbol = attributes.get("gene_name", gene_id)
            transcripts.setdefault((chrom, strand, transcript_id), []).append(
                (start, end, gene_id, symbol)
            )
            key = (chrom, gene_id)
            if key not in gene_spans:
                gene_spans[key] = [start - 1, end, symbol]
            else:
                gene_spans[key][0] = min(int(gene_spans[key][0]), start - 1)
                gene_spans[key][1] = max(int(gene_spans[key][1]), end)

    known: dict[tuple[str, int, int], list[tuple[str, str, str]]] = {}
    for (chrom, strand, _transcript), exons in transcripts.items():
        ordered = sorted(exons, key=lambda item: (item[0], item[1]))
        for left, right in zip(ordered, ordered[1:], strict=False):
            key = (chrom, int(left[1]), int(right[0]) - 1)
            known.setdefault(key, []).append((left[2], left[3], strand))

    bins: dict[tuple[str, int], list[tuple[int, int, str, str]]] = {}
    bin_size = 1_000_000
    for (chrom, gene_id), (start, end, symbol) in gene_spans.items():
        record = (int(start), int(end), gene_id, str(symbol))
        for index in range(int(start) // bin_size, int(end) // bin_size + 1):
            bins.setdefault((chrom, index), []).append(record)

    gene_ids: list[str] = []
    symbols: list[str] = []
    exact_flags: list[bool] = []
    for row in annotated.itertuples(index=False):
        chrom, start, end, strand = (
            str(row.chrom),
            int(row.start),
            int(row.end),
            str(row.strand),
        )
        exact = [
            item
            for item in known.get((chrom, start, end), [])
            if strand == "." or item[2] == strand
        ]
        candidates = [(item[0], item[1]) for item in exact]
        if not candidates:
            nearby: dict[str, tuple[int, int, str, str]] = {}
            for index in range(start // bin_size, end // bin_size + 1):
                for item in bins.get((chrom, index), []):
                    nearby[item[2]] = item
            spans = [
                item for item in nearby.values() if item[0] <= start < end <= item[1]
            ]
            spans.sort(key=lambda item: (item[1] - item[0], item[2]))
            candidates = [(item[2], item[3]) for item in spans[:1]]
        unique = list(dict.fromkeys(candidates))
        gene_ids.append(";".join(item[0] for item in unique))
        symbols.append(";".join(item[1] for item in unique))
        exact_flags.append(bool(exact))
    annotated["gene_id"] = gene_ids
    annotated["gene_symbol"] = symbols
    annotated["annotated"] = exact_flags
    return annotated


def _gtf_attributes(value: str) -> dict[str, str]:
    return {
        key: content
        for key, content in re.findall(r'(\S+)\s+"([^"]*)"', value)
    }


def _normalize_cigar_item(item: tuple[int, int] | tuple[str, str]) -> tuple[int, int]:
    first, second = item
    try:
        first_int, second_int = int(first), int(second)
    except (TypeError, ValueError):
        return 0, 0
    if first_int in {0, 1, 2, 3, 4, 5, 6, 7, 8}:
        return first_int, second_int
    return second_int, first_int


def _splice_site_components(choices: pd.DataFrame) -> pd.Series:
    """Cluster junctions connected by a shared donor or acceptor splice site."""
    parent = list(range(len(choices)))

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left: int, right: int) -> None:
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    for _, group in choices.groupby(["chrom", "strand"], sort=False):
        by_start: dict[int, int] = {}
        by_end: dict[int, int] = {}
        for index, row in group.iterrows():
            start, end = int(row["start"]), int(row["end"])
            if start in by_start:
                union(index, by_start[start])
            else:
                by_start[start] = index
            if end in by_end:
                union(index, by_end[end])
            else:
                by_end[end] = index
    roots = [find(index) for index in range(len(choices))]
    labels = {root: position + 1 for position, root in enumerate(dict.fromkeys(roots))}
    return pd.Series(
        [f"{choices.iloc[index]['chrom']}:{choices.iloc[index]['strand']}:JC{labels[root]:06d}" for index, root in enumerate(roots)],
        index=choices.index,
    )
