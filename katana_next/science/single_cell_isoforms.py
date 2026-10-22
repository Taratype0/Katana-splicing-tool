"""Donor-aware transcript usage, separate from bulk co-splicing networks.

Motivated by Mantri et al., bioRxiv 10.64898/2026.09.11.748700.
This is not an implementation of their long-read preprocessing or cell-level LRT.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
import pandas as pd

from katana_next.science.experimental_design import donor_covariates
from katana_next.science.statistics import bh_adjust, dirichlet_multinomial_two_group


@dataclass
class IsoformUsageResult:
    usage: pd.DataFrame
    diversity: pd.DataFrame
    differential: pd.DataFrame
    effects: pd.DataFrame
    diagnostics: pd.DataFrame


def validate_isoform_inputs(counts: pd.DataFrame, metadata: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    required_counts = {"cell_id", "gene_id", "transcript_id", "count"}
    required_meta = {"cell_id", "donor_id", "cell_type", "condition"}
    for frame, required in ((counts, required_counts), (metadata, required_meta)):
        if not required.issubset(frame):
            raise ValueError(f"Missing required columns: {sorted(required - set(frame))}")
        for col in required - {"count"}:
            if frame[col].isna().any() or frame[col].astype(str).str.strip().eq("").any():
                raise ValueError(f"Empty identity: {col}")
    counts, metadata = counts.copy(), metadata.copy()
    for col in required_counts - {"count"}:
        counts[col] = counts[col].astype(str)
    for col in required_meta:
        metadata[col] = metadata[col].astype(str)
    if metadata.cell_id.duplicated().any():
        raise ValueError("cell_id must be globally unique, including the sample prefix")
    if counts.duplicated(["cell_id", "transcript_id"]).any():
        raise ValueError("Duplicate cell-transcript counts; aggregate technical partitions explicitly")
    if counts.groupby("transcript_id").gene_id.nunique().gt(1).any():
        raise ValueError("A transcript must belong to exactly one gene")
    if not set(counts.cell_id).issubset(set(metadata.cell_id)):
        raise ValueError("Count matrix contains cells absent from metadata")
    values = pd.to_numeric(counts["count"], errors="raise").to_numpy(float)
    if not np.isfinite(values).all() or (values < 0).any() or (values % 1).any():
        raise ValueError("Requires nonnegative integer molecule/uniquely assigned fragment counts, not TPM or fractional estimates")
    if (values > np.iinfo(np.int64).max // 1000).any():
        raise ValueError("Counts exceed safe integer range")
    counts["count"] = values.astype(np.int64)
    if counts.empty or counts["count"].sum() == 0:
        raise ValueError("No detected transcript molecules")
    if "batch" not in metadata:
        metadata["batch"] = ""
    metadata["batch"] = metadata.batch.fillna("").astype(str)
    if metadata.groupby(["donor_id", "condition"]).batch.nunique().gt(1).any():
        raise ValueError("Multiple batches within donor-condition; resolve design before aggregation")
    return counts, metadata


def _entropy(counts: np.ndarray) -> tuple[float, float]:
    total = counts.sum()
    if total <= 0:
        return float("nan"), float("nan")
    p = counts[counts > 0] / total
    return float(-np.log2(np.square(p).sum())), float(-(p * np.log2(p)).sum())


def analyze_isoform_usage(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    numerator: str,
    denominator: str,
    min_donors: int = 3,
    min_gene_count: int = 20,
    min_expressing_cells: int = 5,
    rarefaction_depth: int = 20,
    rarefaction_draws: int = 30,
    seed: int = 17,
    adjust_batch: bool = True,
    checkpoint: Callable | None = None,
    progress: Callable[[float, str], None] | None = None,
) -> IsoformUsageResult:
    """Sum counts within donor/cell-type; never treat cells as replicates.

    Absent genes have unmeasured usage, not zero usage. Rarefaction samples
    without replacement within each measured donor-gene; its spread is Monte
    Carlo sensitivity, not a confidence interval on biological uncertainty.
    """
    if numerator == denominator or not numerator or not denominator:
        raise ValueError("Select two distinct, explicit conditions")
    if min_donors < 2 or min(min_gene_count, min_expressing_cells, rarefaction_depth, rarefaction_draws) < 1:
        raise ValueError("Invalid replication or measurement thresholds")
    counts, metadata = validate_isoform_inputs(counts, metadata)
    if not {numerator, denominator}.issubset(set(metadata.condition)):
        raise ValueError("Selected conditions not found in metadata")
    unit = ["donor_id", "condition", "cell_type", "batch"]
    joined = counts.merge(metadata[list(dict.fromkeys(["cell_id", *unit]))], on="cell_id", validate="many_to_one")
    joined = joined[joined.condition.isin([numerator, denominator])]
    # Retain the complete supplied transcript catalogue, including explicit zeros.
    catalogue = counts.groupby("gene_id", sort=True).transcript_id.unique()
    expressing = joined[joined["count"].gt(0)].groupby([*unit, "gene_id"]).cell_id.nunique()
    summed = joined.groupby([*unit, "gene_id", "transcript_id"], sort=True)["count"].sum()
    cell_totals = joined.groupby(["cell_id", "gene_id"])["count"].transform("sum")
    joined["cell_fraction"] = joined["count"].div(cell_totals.where(cell_totals.gt(0)))
    fraction_sums = joined.groupby([*unit, "gene_id", "transcript_id"]).cell_fraction.sum()
    rng = np.random.default_rng(seed)
    usage_rows, diversity_rows, qc_rows, tests, effects = [], [], [], [], []
    units = metadata[metadata.condition.isin([numerator, denominator])][unit].drop_duplicates()
    for gene_index, (gene, transcripts) in enumerate(catalogue.items()):
        if checkpoint:
            checkpoint()
        if progress and gene_index % max(1, len(catalogue) // 100) == 0:
            progress(.1 + .7 * gene_index / len(catalogue), f"Transcript usage: gene {gene_index + 1}/{len(catalogue)}")
        gene_units = []
        for key in units.itertuples(index=False, name=None):
            key = tuple(key)
            vector = np.array([summed.get((*key, gene, t), 0) for t in transcripts], dtype=np.int64)
            total = int(vector.sum())
            cells = int(expressing.get((*key, gene), 0))
            measured = total >= min_gene_count and cells >= min_expressing_cells
            identity = dict(zip(unit, key)) | {"gene_id": gene}
            qc_rows.append(identity | {"gene_count": total, "expressing_cells": cells, "measured": measured,
                "reason": "measured" if measured else "low_gene_count" if total < min_gene_count else "few_expressing_cells"})
            fraction = vector / total if measured else np.full(len(vector), np.nan)
            usage_rows.extend(identity | {"transcript_id": t, "count": int(v), "gene_count": total,
                "fraction": float(p), "measured": measured} for t, v, p in zip(transcripts, vector, fraction))
            h2, h1 = _entropy(vector) if measured else (np.nan, np.nan)
            cell_mean = np.array([fraction_sums.get((*key, gene, t), 0) for t in transcripts]) / max(cells, 1)
            cell_h2, cell_h1 = _entropy(cell_mean) if measured else (np.nan, np.nan)
            rare = []
            if measured and total >= rarefaction_depth:
                rare = [_entropy(rng.multivariate_hypergeometric(vector, rarefaction_depth))[0] for _ in range(rarefaction_draws)]
            diversity_rows.append(identity | {"renyi2_bits": h2, "shannon_bits": h1,
                "mean_cell_fraction_renyi2_bits": cell_h2, "mean_cell_fraction_shannon_bits": cell_h1,
                "rarefied_renyi2_mean": float(np.mean(rare)) if rare else np.nan,
                "rarefied_mc_sd": float(np.std(rare)) if rare else np.nan,
                "rarefaction_depth": rarefaction_depth, "rarefaction_status": "complete" if rare else "not_measurable",
                "gene_count": total, "expressing_cells": cells, "transcripts_in_catalogue": len(transcripts)})
            if measured:
                gene_units.append((identity, vector))
        for cell_type in units.cell_type.unique():
            members = [(identity, vector) for identity, vector in gene_units if identity["cell_type"] == cell_type]
            rows = pd.DataFrame([identity for identity, _ in members])
            record = {"gene_id": gene, "cell_type": cell_type, "numerator": numerator, "denominator": denominator,
                "pvalue": np.nan, "model_status": "insufficient_donors", "numerator_donors": 0, "denominator_donors": 0}
            if not members:
                tests.append(record)
                continue
            arm = rows.condition.eq(numerator).to_numpy()
            record.update(numerator_donors=int(arm.sum()), denominator_donors=int((~arm).sum()))
            matrix = np.vstack([vector for _, vector in members])
            proportions = matrix / matrix.sum(axis=1, keepdims=True)
            if arm.any() and (~arm).any():
                delta = proportions[arm].mean(axis=0) - proportions[~arm].mean(axis=0)
                loo = []
                for donor in rows.donor_id.unique():
                    keep = rows.donor_id.ne(donor).to_numpy()
                    # Keep biological units together across paired conditions.
                    if min(int((keep & arm).sum()), int((keep & ~arm).sum())) >= 2:
                        loo.append(proportions[keep & arm].mean(axis=0) - proportions[keep & ~arm].mean(axis=0))
                for index, (t, change) in enumerate(zip(transcripts, delta)):
                    consistency = float(np.mean(np.sign(np.asarray(loo)[:, index]) == np.sign(change))) if loo and change != 0 else np.nan
                    effects.append({"gene_id": gene, "cell_type": cell_type, "transcript_id": t,
                        "numerator": numerator, "denominator": denominator, "delta_fraction": float(change),
                        "lodo_direction_consistency": consistency, "lodo_evaluable": len(loo),
                        "lodo_total_donors": int(rows.donor_id.nunique()),
                        "lodo_status": "complete" if len(loo) == rows.donor_id.nunique() else "partial" if loo else "not_estimable"})
            observed_choices = matrix.sum(axis=0) > 0
            model_matrix = matrix[:, observed_choices]
            pooled = []
            if model_matrix.shape[1] > 30:
                # The declared paper-style dimensionality cap pools the tail;
                # it does not delete molecules or select on condition effects.
                order = np.argsort(-model_matrix.sum(axis=0), kind="stable")
                pooled = np.asarray(transcripts)[observed_choices][order[29:]].tolist()
                model_matrix = np.column_stack([model_matrix[:, order[:29]], model_matrix[:, order[29:]].sum(axis=1)])
            record.update(tested_transcript_choices=model_matrix.shape[1], pooled_transcripts=";".join(pooled))
            if observed_choices.sum() < 2:
                record["model_status"] = "single_transcript"
            elif min(int(arm.sum()), int((~arm).sum())) >= min_donors:
                covariates, paired, batch = donor_covariates(rows, adjust_for_batch=adjust_batch)
                fit = dirichlet_multinomial_two_group(model_matrix, arm.astype(int), covariates=covariates)
                record.update(pvalue=fit["pvalue"], model_status=fit["status"], likelihood_ratio=fit["likelihood_ratio"],
                    paired_donor_adjusted=paired, batch_adjusted=batch,
                    inference_note="Small donor counts require independent validation; asymptotic LRT")
            tests.append(record)
    differential = pd.DataFrame(tests)
    # One declared family across tested genes and cell types, never after significance selection.
    differential["fdr"] = bh_adjust(differential.pvalue)
    effect_frame = pd.DataFrame(effects, columns=["gene_id", "cell_type", "transcript_id", "numerator", "denominator", "delta_fraction", "lodo_direction_consistency", "lodo_evaluable", "lodo_total_donors", "lodo_status"])
    effect_frame = effect_frame.merge(differential[["gene_id", "cell_type", "fdr", "model_status"]], on=["gene_id", "cell_type"], how="left")
    return IsoformUsageResult(pd.DataFrame(usage_rows), pd.DataFrame(diversity_rows), differential, effect_frame, pd.DataFrame(qc_rows))
