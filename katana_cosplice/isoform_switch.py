from __future__ import annotations

from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from katana_cosplice.utils import bh_adjust, normalize_gene_id


def _attributes(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in text.split(";"):
        item = item.strip()
        if item and " " in item:
            key, value = item.split(" ", 1)
            result[key] = value.strip().strip('"')
    return result


def transcript_annotation(gtf: Path) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    with gtf.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            attrs = _attributes(fields[8])
            transcript = attrs.get("transcript_id", "")
            gene_id = normalize_gene_id(attrs.get("gene_id", ""))
            if transcript and gene_id:
                rows.append({
                    "transcript_id": transcript.split(".")[0],
                    "gene_id": gene_id,
                    "gene_symbol": attrs.get("gene_name", ""),
                    "transcript_name": attrs.get("transcript_name", ""),
                    "transcript_type": attrs.get("transcript_type", attrs.get("transcript_biotype", "")),
                })
    return pd.DataFrame(rows).drop_duplicates("transcript_id")


def quantify_isoform_switches(
    quant_root: Path,
    gtf: Path,
    sample_metadata: pd.DataFrame,
    module_genes: pd.DataFrame,
    *,
    condition_column: str,
    min_gene_tpm: float,
    min_abs_delta_if: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    annotation = transcript_annotation(gtf)
    if annotation.empty or "sample_id" not in sample_metadata or condition_column not in sample_metadata:
        return pd.DataFrame(), pd.DataFrame()
    wanted_genes = set(module_genes["gene_id"].map(normalize_gene_id)) - {""}
    annotation = annotation.loc[annotation["gene_id"].isin(wanted_genes)]
    tx_to_gene = annotation.set_index("transcript_id")["gene_id"]
    quant_parts: list[pd.DataFrame] = []
    valid_samples = set(sample_metadata["sample_id"].astype(str))
    for path in sorted(quant_root.glob("**/quant.sf")):
        sample = path.parent.name
        if sample not in valid_samples:
            continue
        frame = pd.read_csv(path, sep="\t", usecols=["Name", "TPM"])
        frame["transcript_id"] = frame["Name"].astype(str).str.split(".").str[0]
        frame["gene_id"] = frame["transcript_id"].map(tx_to_gene)
        frame["TPM"] = pd.to_numeric(frame["TPM"], errors="coerce")
        frame = frame.dropna(subset=["gene_id", "TPM"])
        if frame.empty:
            continue
        frame["sample_id"] = sample
        gene_total = frame.groupby("gene_id", observed=True)["TPM"].transform("sum")
        frame["isoform_fraction"] = np.where(gene_total > 0, frame["TPM"] / gene_total, np.nan)
        quant_parts.append(frame[["sample_id", "gene_id", "transcript_id", "TPM", "isoform_fraction"]])
    if not quant_parts:
        return pd.DataFrame(), pd.DataFrame()
    quant = pd.concat(quant_parts, ignore_index=True).merge(
        sample_metadata[["sample_id", condition_column]].astype({"sample_id": str}),
        on="sample_id", how="inner", validate="many_to_one",
    )
    means = quant.groupby(["gene_id", "transcript_id", condition_column], as_index=False, observed=True).agg(
        mean_TPM=("TPM", "mean"), mean_isoform_fraction=("isoform_fraction", "mean"), n_samples=("sample_id", "nunique")
    )
    gene_tpm = means.groupby(["gene_id", condition_column], observed=True)["mean_TPM"].sum()
    rows: list[dict[str, object]] = []
    conditions = sorted(means[condition_column].astype(str).unique())
    annotation_index = annotation.set_index("transcript_id")
    for condition_a, condition_b in combinations(conditions, 2):
        left = means.loc[means[condition_column].astype(str).eq(condition_a)]
        right = means.loc[means[condition_column].astype(str).eq(condition_b)]
        for gene_id in sorted(set(left["gene_id"]) & set(right["gene_id"])):
            if min(gene_tpm.get((gene_id, condition_a), 0), gene_tpm.get((gene_id, condition_b), 0)) < min_gene_tpm:
                continue
            a = left.loc[left["gene_id"].eq(gene_id)].set_index("transcript_id")
            b = right.loc[right["gene_id"].eq(gene_id)].set_index("transcript_id")
            transcripts = a.index.union(b.index)
            frac_a = a["mean_isoform_fraction"].reindex(transcripts).fillna(0.0)
            frac_b = b["mean_isoform_fraction"].reindex(transcripts).fillna(0.0)
            dominant_a = str(frac_a.idxmax())
            dominant_b = str(frac_b.idxmax())
            delta = frac_b - frac_a
            max_abs_delta = float(delta.abs().max())
            if dominant_a == dominant_b or max_abs_delta < min_abs_delta_if:
                continue
            symbol = ""
            if dominant_a in annotation_index.index:
                symbol = str(annotation_index.loc[dominant_a, "gene_symbol"])
            rows.append({
                "gene_id": gene_id, "gene_symbol": symbol,
                "condition_A": condition_a, "condition_B": condition_b,
                "dominant_transcript_A": dominant_a, "dominant_transcript_B": dominant_b,
                "dominant_IF_A": float(frac_a[dominant_a]), "dominant_IF_B": float(frac_b[dominant_b]),
                "delta_IF_new_dominant": float(delta[dominant_b]),
                "max_abs_delta_IF": max_abs_delta,
                "gene_mean_TPM_A": float(gene_tpm[(gene_id, condition_a)]),
                "gene_mean_TPM_B": float(gene_tpm[(gene_id, condition_b)]),
                "switch_definition": "dominant transcript changes and max abs(delta isoform fraction) passes configured cutoff",
                "statistical_scope": "descriptive switch; statistical support is inherited from configured DTU evidence",
            })
    switches = pd.DataFrame(rows)
    if switches.empty:
        return switches, pd.DataFrame()
    support_columns = [c for c in ("gene_id", "DTU_significant", "best_DTU_qvalue") if c in module_genes]
    if support_columns:
        support = module_genes[support_columns].drop_duplicates("gene_id")
        switches = switches.merge(support, on="gene_id", how="left", validate="many_to_one")
    if "DTU_significant" not in switches:
        switches["DTU_significant"] = False
    else:
        switches["DTU_significant"] = switches["DTU_significant"].fillna(False).astype(bool)
    enrichment = enrich_isoform_switches(switches, module_genes)
    return switches.sort_values("max_abs_delta_IF", ascending=False), enrichment


def enrich_isoform_switches(
    switches: pd.DataFrame,
    module_genes: pd.DataFrame,
) -> pd.DataFrame:
    """Test module enrichment for descriptive and DTU-supported switches separately."""
    if switches.empty or not {"module", "gene_id"}.issubset(module_genes):
        return pd.DataFrame()
    module_map = module_genes[["module", "gene_id"]].drop_duplicates()
    universe = set(module_map["gene_id"])
    enrichment_rows = []
    scopes = {"descriptive_switch": switches}
    if "DTU_significant" in switches:
        dtu = switches.loc[switches["DTU_significant"].fillna(False).astype(bool)]
        scopes["DTU_supported_switch"] = dtu
    for switch_scope, scope_frame in scopes.items():
        selected = set(scope_frame["gene_id"]) & universe
        for module, part in module_map.groupby("module", observed=True):
            members = set(part["gene_id"]) & universe
            overlap = members & selected
            pvalue = (
                float(hypergeom.sf(len(overlap) - 1, len(universe), len(selected), len(members)))
                if overlap and selected else 1.0
            )
            enrichment_rows.append({
                "module": module,
                "switch_scope": switch_scope,
                "module_genes": len(members),
                "switch_genes_in_universe": len(selected),
                "overlap_switch_genes": len(overlap),
                "overlap_genes": ";".join(sorted(overlap)),
                "pvalue": pvalue,
                "method": (
                    "one-sided hypergeometric enrichment of dominant-isoform switches; "
                    "DTU-supported scope requires inherited DTU_significant evidence"
                ),
            })
    enrichment = pd.DataFrame(enrichment_rows)
    enrichment["FDR"] = bh_adjust(enrichment["pvalue"])
    return enrichment.sort_values(["switch_scope", "FDR", "pvalue"])
