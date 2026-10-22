from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd

from katana_next.science.statistics import (
    beta_binomial_group_test,
    bh_adjust,
    dirichlet_multinomial_two_group,
)
from katana_next.services.project_service import ProjectContext
from katana_next.science.experimental_design import donor_covariates
from katana_next.infrastructure.model_checkpoint import model_groups


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
        if frame[column].isna().any() or frame[column].astype(str).str.strip().eq("").any():
            raise ValueError(f"Single-cell metadata requires nonempty {column} in every row.")
        frame[column] = frame[column].astype(str).str.strip()
    if frame["cell_id"].duplicated().any():
        duplicates = (
            frame.loc[frame["cell_id"].duplicated(), "cell_id"].head(5).tolist()
        )
        raise ValueError(
            f"Single-cell metadata contains duplicate cell IDs: {duplicates}"
        )
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
    frame = (pd.read_parquet(asset.path) if asset.path.suffix.lower() == ".parquet"
             else pd.read_csv(asset.path, sep="\t") if asset.path.suffix.lower() == ".tsv"
             else pd.read_csv(asset.path, sep=None, engine="python"))
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
        | ~np.isfinite(frame[["inclusion", "exclusion"]]).all(axis=1)
        | frame[["inclusion", "exclusion"]].mod(1).ne(0).any(axis=1)
        | frame["inclusion"].lt(0)
        | frame["exclusion"].lt(0)
    )
    if invalid.any():
        raise ValueError(
            f"Single-cell junction counts contain {int(invalid.sum())} invalid rows."
        )
    if frame.duplicated(["cell_id", "event_id"]).any():
        raise ValueError("Duplicate cell-event count rows: merge disjoint technical counts explicitly before importing.")
    defaults = {
        "gene_id": "",
        "gene_symbol": "",
        "event_type": "junction_choice",
        "cluster_id": "",
        "annotated": False,
    }
    for column, default in defaults.items():
        if column not in frame:
            frame[column] = default
    frame["gene_id"] = frame["gene_id"].fillna("").astype(str)
    frame["gene_symbol"] = frame["gene_symbol"].fillna("").astype(str)
    frame["event_type"] = frame["event_type"].fillna("junction_choice").astype(str)
    frame["cluster_id"] = frame["cluster_id"].fillna("").astype(str)
    frame.loc[frame["cluster_id"].eq(""), "cluster_id"] = frame.loc[
        frame["cluster_id"].eq(""), "event_id"
    ]
    frame["annotated"] = frame["annotated"].fillna(False).astype(str).str.lower().isin({"true", "1", "yes"})
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
        metadata.groupby(
            ["condition", "cell_type", "donor_id", "sample_id", "assay"], dropna=False
        )
        .size()
        .rename("cells")
        .reset_index()
    )
    design["passes_cell_count"] = design["cells"].ge(min_cells_per_donor_cell_type)
    design["row_id"] = (
        design["condition"] + "::" + design["cell_type"] + "::" + design["donor_id"]
    )
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
            "status": "warning"
            if metadata["assay"].eq("unspecified").any()
            else "pass",
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
        metadata[
            [
                "cell_id",
                "sample_id",
                "donor_id",
                "cell_type",
                "condition",
                "batch",
                "assay",
            ]
        ],
        on="cell_id",
        how="inner",
        validate="many_to_one",
    )
    merged["detected"] = merged["total"].ge(min_cell_event_reads)
    grouping = [
        "event_id",
        "cluster_id",
        "gene_id",
        "gene_symbol",
        "event_type",
        "annotated",
        "donor_id",
        "cell_type",
        "condition",
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
    # These attributes belong to the observation, not each of its events.
    unit_keys = ["donor_id", "cell_type", "condition"]
    unit_metadata = metadata.groupby(unit_keys, dropna=False).agg(
        sample_id=("sample_id", _joined_unique),
        batch=("batch", _joined_unique),
        assay=("assay", _joined_unique),
    ).reset_index()
    aggregate = aggregate.merge(unit_metadata, on=unit_keys, validate="many_to_one")
    aggregate["total"] = aggregate["inclusion"] + aggregate["exclusion"]
    aggregate["psi_posterior"] = (aggregate["inclusion"] + 0.5) / (
        aggregate["total"] + 1.0
    )
    aggregate["measurement_se"] = np.sqrt(
        aggregate["psi_posterior"]
        * (1 - aggregate["psi_posterior"])
        / (aggregate["total"] + 1.0)
    )
    aggregate["row_id"] = (
        aggregate["event_id"]
        + "::"
        + aggregate["cell_type"]
        + "::"
        + aggregate["donor_id"]
        + "::"
        + aggregate["condition"]
    )
    coverage = (
        aggregate.groupby(["cell_type", "condition", "donor_id"], dropna=False)
        .agg(
            cells=("contributing_cells", "max"),
            measured_events=(
                "total",
                lambda values: int(pd.Series(values).gt(0).sum()),
            ),
            median_event_reads=("total", "median"),
            median_measurement_se=("measurement_se", "median"),
        )
        .reset_index()
    )
    coverage["row_id"] = (
        coverage["cell_type"]
        + "::"
        + coverage["condition"]
        + "::"
        + coverage["donor_id"]
    )
    return aggregate, coverage


def gene_splicing_deviation_scores(
    metadata: pd.DataFrame,
    counts: pd.DataFrame,
    *,
    min_gene_reads: float = 3.0,
    min_gene_events: int = 2,
) -> pd.DataFrame:
    """Compute a count-standardized gene-level junction deviation per cell.

    This is a transparent sparse-data summary inspired by the role of SpliZ,
    but it is not an implementation of the published SpliZ rank statistic.
    """

    annotated = counts.copy()
    annotated["gene_key"] = annotated["gene_symbol"].where(
        annotated["gene_symbol"].ne(""), annotated["gene_id"]
    )
    annotated = annotated[annotated["gene_key"].ne("")].copy()
    if annotated.empty:
        return _empty_gene_score_table()
    event_counts = annotated.groupby("gene_key")["event_id"].nunique()
    eligible_genes = event_counts[event_counts.ge(min_gene_events)].index
    annotated = annotated[annotated["gene_key"].isin(eligible_genes)].copy()
    if annotated.empty:
        return _empty_gene_score_table()

    global_usage = annotated.groupby(["gene_key", "event_id"], as_index=False).agg(
        inclusion=("inclusion", "sum"), total=("total", "sum")
    )
    global_usage["expected_psi"] = (global_usage["inclusion"] + 0.5) / (
        global_usage["total"] + 1.0
    )
    scored = annotated.merge(
        global_usage[["gene_key", "event_id", "expected_psi"]],
        on=["gene_key", "event_id"],
        how="left",
        validate="many_to_one",
    )
    variance = (
        scored["total"] * scored["expected_psi"] * (1.0 - scored["expected_psi"])
    ).clip(lower=0.25)
    scored["standardized_residual"] = (
        scored["inclusion"] - scored["total"] * scored["expected_psi"]
    ) / np.sqrt(variance)
    output = scored.groupby(["cell_id", "gene_key"], as_index=False).agg(
        junction_deviation_z=("standardized_residual", "sum"),
        total_junction_reads=("total", "sum"),
        measured_events=("event_id", lambda values: int(pd.Series(values).nunique())),
    )
    output = output[
        output["total_junction_reads"].ge(min_gene_reads)
        & output["measured_events"].ge(min_gene_events)
    ].copy()
    if output.empty:
        return _empty_gene_score_table()
    output["junction_deviation_z"] /= np.sqrt(
        output["total_junction_reads"].clip(lower=1.0)
    )
    output = output.merge(
        metadata[["cell_id", "sample_id", "donor_id", "cell_type", "condition"]],
        on="cell_id",
        how="left",
        validate="many_to_one",
    )
    output["row_id"] = output["cell_id"] + "::" + output["gene_key"]
    return output


def differential_pseudobulk(
    pseudobulk: pd.DataFrame,
    *,
    numerator: str,
    denominator: str,
    min_total_reads: float,
    min_donors_per_condition: int,
    fdr_cutoff: float,
    min_abs_dpsi: float,
    bootstrap_iterations: int = 200,
    random_seed: int = 17,
    adjust_for_batch: bool = True,
    checkpoint=None,
    report=None,
) -> pd.DataFrame:
    selected = pseudobulk[pseudobulk["condition"].isin([numerator, denominator])].copy()
    rows: list[dict[str, Any]] = []
    grouped = selected.groupby(["cell_type", "event_id"], sort=False)
    for (cell_type, event_id), group in model_groups(grouped, rows, checkpoint, report):
        group = group[group["total"].ge(min_total_reads)].copy()
        donors = group.groupby("condition")["donor_id"].nunique()
        if (
            donors.get(numerator, 0) < min_donors_per_condition
            or donors.get(denominator, 0) < min_donors_per_condition
        ):
            continue
        groups = group["condition"].eq(numerator).astype(int)
        batch_covariates, paired, batch_adjusted = donor_covariates(group, adjust_for_batch=adjust_for_batch)
        batch_levels = group["batch"].nunique()
        fit = beta_binomial_group_test(
            group["inclusion"],
            group["total"],
            groups,
            covariates=batch_covariates,
            analytic_gradient=True,
        )
        consistency, ci_low, ci_high = _replicate_effect_diagnostics(
            group,
            numerator=numerator,
            denominator=denominator,
            iterations=bootstrap_iterations,
            seed=random_seed + len(rows),
        )
        first = group.iloc[0]
        rows.append(
            {
                "row_id": f"{cell_type}::{event_id}",
                "cell_type": cell_type,
                "event_id": event_id,
                "cluster_id": first.get("cluster_id", event_id),
                "gene_id": first.get("gene_id", ""),
                "gene_symbol": first.get("gene_symbol", ""),
                "event_type": first.get("event_type", "junction_choice"),
                "numerator": numerator,
                "denominator": denominator,
                "numerator_donors": int(donors.get(numerator, 0)),
                "denominator_donors": int(donors.get(denominator, 0)),
                "numerator_cells": int(
                    group.loc[
                        group["condition"].eq(numerator), "contributing_cells"
                    ].sum()
                ),
                "denominator_cells": int(
                    group.loc[
                        group["condition"].eq(denominator), "contributing_cells"
                    ].sum()
                ),
                "dpsi": fit["dpsi"],
                "log_odds": fit["log_odds"],
                "overdispersion_rho": fit["rho"],
                "dpsi_ci_low": ci_low,
                "dpsi_ci_high": ci_high,
                "leave_one_donor_out_consistency": consistency,
                "pvalue": fit["pvalue"],
                "model_status": fit["status"],
                "null_log_likelihood": fit.get("null_log_likelihood", np.nan),
                "full_log_likelihood": fit.get("full_log_likelihood", np.nan),
                "optimizer": fit.get("optimizer", "analytic_gradient_nested_start"),
                "batch_adjusted": batch_adjusted,
                "paired_donor_adjusted": paired,
                "batch_levels": batch_levels,
            }
        )
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    result["fdr"] = result.groupby("cell_type")["pvalue"].transform(
        lambda values: bh_adjust(values)
    )
    result["significant"] = result["fdr"].le(fdr_cutoff) & result["dpsi"].abs().ge(
        min_abs_dpsi
    )
    result["direction"] = np.where(
        result["dpsi"].ge(0),
        "higher inclusion in numerator",
        "lower inclusion in numerator",
    )
    return result.sort_values(
        ["significant", "fdr", "dpsi"],
        ascending=[False, True, False],
        na_position="last",
    )


def differential_multijunction_pseudobulk(
    pseudobulk: pd.DataFrame,
    *,
    numerator: str,
    denominator: str,
    min_total_reads: float,
    min_donors_per_condition: int,
    fdr_cutoff: float,
    min_abs_delta_usage: float,
    adjust_for_batch: bool = True,
    checkpoint=None,
    report=None,
) -> pd.DataFrame:
    """Fit a donor-level Dirichlet-multinomial omnibus test per junction group."""

    required = {
        "cluster_id",
        "event_id",
        "inclusion",
        "donor_id",
        "cell_type",
        "condition",
    }
    if not required.issubset(pseudobulk.columns):
        return _empty_multijunction_table()
    selected = pseudobulk[pseudobulk["condition"].isin([numerator, denominator])].copy()
    rows: list[dict[str, Any]] = []
    grouped = selected.groupby(["cell_type", "cluster_id"], sort=False)
    for (cell_type, cluster_id), group in model_groups(grouped, rows, checkpoint, report):
        if group["event_id"].nunique() < 2:
            continue
        matrix = group.pivot_table(
            index=["condition", "donor_id", "batch"],
            columns="event_id",
            values="inclusion",
            aggfunc="sum",
            fill_value=0,
        )
        totals = matrix.sum(axis=1)
        matrix = matrix.loc[totals.ge(min_total_reads)]
        if matrix.empty:
            continue
        conditions = matrix.index.get_level_values("condition")
        model_units = matrix.index.to_frame(index=False)
        donor_counts = model_units.groupby("condition")["donor_id"].nunique()
        if (
            donor_counts.get(numerator, 0) < min_donors_per_condition
            or donor_counts.get(denominator, 0) < min_donors_per_condition
        ):
            continue
        groups = np.asarray(conditions == numerator, dtype=int)
        batch_covariates, paired, batch_adjusted = donor_covariates(model_units, adjust_for_batch=adjust_for_batch)
        batch_levels = model_units["batch"].nunique()
        fit = dirichlet_multinomial_two_group(
            matrix.to_numpy(float), groups, covariates=batch_covariates
        )
        event_ids = matrix.columns.astype(str).tolist()
        delta = np.asarray(fit.get("delta_proportions", []), dtype=float)
        if len(delta) != len(event_ids):
            delta = np.repeat(np.nan, len(event_ids))
        strongest_index = (
            int(np.nanargmax(np.abs(delta))) if np.isfinite(delta).any() else 0
        )
        maximum_delta = (
            float(np.max(np.abs(delta[np.isfinite(delta)])))
            if np.isfinite(delta).any()
            else float("nan")
        )
        rows.append(
            {
                "row_id": f"{cell_type}::{cluster_id}",
                "cell_type": str(cell_type),
                "cluster_id": str(cluster_id),
                "gene_id": _joined_unique(group["gene_id"]),
                "gene_symbol": _joined_unique(group["gene_symbol"]),
                "junction_choices": len(event_ids),
                "event_ids": ";".join(event_ids),
                "strongest_event_id": event_ids[strongest_index],
                "max_abs_delta_usage": maximum_delta,
                "strongest_delta_usage": float(delta[strongest_index]),
                "numerator_donors": int(donor_counts.get(numerator, 0)),
                "denominator_donors": int(donor_counts.get(denominator, 0)),
                "overdispersion_rho": fit["rho"],
                "likelihood_ratio": fit["likelihood_ratio"],
                "degrees_of_freedom": fit["degrees_of_freedom"],
                "pvalue": fit["pvalue"],
                "model_status": fit["status"],
                "batch_adjusted": batch_adjusted,
                "paired_donor_adjusted": paired,
                "batch_levels": batch_levels,
            }
        )
    result = pd.DataFrame(rows)
    if result.empty:
        return _empty_multijunction_table()
    result["fdr"] = result.groupby("cell_type")["pvalue"].transform(bh_adjust)
    result["significant"] = result["fdr"].le(fdr_cutoff) & result[
        "max_abs_delta_usage"
    ].ge(min_abs_delta_usage)
    return result.sort_values(
        ["significant", "fdr", "max_abs_delta_usage"],
        ascending=[False, True, False],
        na_position="last",
    )


def _replicate_effect_diagnostics(
    group: pd.DataFrame,
    *,
    numerator: str,
    denominator: str,
    iterations: int,
    seed: int,
) -> tuple[float, float, float]:
    donor = group.groupby(["condition", "donor_id"], as_index=False).agg(
        inclusion=("inclusion", "sum"), total=("total", "sum")
    )

    codes, donor_ids = pd.factorize(donor.donor_id, sort=False)
    psi = (donor.inclusion.to_numpy(float)+.5)/(donor.total.to_numpy(float)+1)
    masks = [donor.condition.eq(c).to_numpy() for c in (numerator, denominator)]
    sums = np.column_stack([np.bincount(codes, weights=psi*mask, minlength=len(donor_ids)) for mask in masks])
    sizes = np.column_stack([np.bincount(codes, weights=mask, minlength=len(donor_ids)) for mask in masks])
    totals, counts = sums.sum(axis=0), sizes.sum(axis=0)
    means = np.divide(totals, counts, out=np.full(2,np.nan), where=counts>0)
    full = means[0]-means[1]
    loo_counts = counts-sizes
    loo = np.divide(totals-sums, loo_counts, out=np.full_like(sums,np.nan), where=loo_counts>0)
    leave_one_out = loo[:,0]-loo[:,1]
    valid = leave_one_out[np.isfinite(leave_one_out)]
    consistency = (
        float(np.mean(np.sign(valid) == np.sign(full)))
        if len(valid) and np.isfinite(full) and full != 0
        else float("nan")
    )
    if iterations < 2:
        return consistency, float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    # Freeze donor blocks once; rebuilding pandas groups inside each draw is costly.
    profiles = [tuple(sorted(c for c,n in zip((numerator,denominator), row) if n>0)) for row in sizes]
    strata = [np.asarray([j for j,p in enumerate(profiles) if p == profile])
              for profile in dict.fromkeys(profiles)]
    if len(strata) == 1:
        selected = rng.choice(strata[0],size=(iterations,len(strata[0])),replace=True)
        means = sums[selected].sum(axis=1)/sizes[selected].sum(axis=1)
        bootstrap = means[:,0]-means[:,1]
    else:
        bootstrap = []
        for _ in range(iterations):
            selected = np.concatenate([rng.choice(indices, size=len(indices), replace=True) for indices in strata])
            means = sums[selected].sum(axis=0)/sizes[selected].sum(axis=0)
            bootstrap.append(float(means[0]-means[1]))
    finite = np.asarray(
        [value for value in bootstrap if np.isfinite(value)], dtype=float
    )
    if not len(finite):
        return consistency, float("nan"), float("nan")
    low, high = np.quantile(finite, [0.025, 0.975])
    return consistency, float(low), float(high)


def _joined_unique(values: pd.Series) -> str:
    return ";".join(sorted({str(value) for value in values if str(value).strip()}))


def _categorical_covariates(
    values: pd.Series, *, enabled: bool
) -> tuple[np.ndarray | None, int]:
    labels = values.fillna("").astype(str).str.strip()
    levels = int(labels.nunique())
    if not enabled or levels <= 1:
        return None, levels
    encoded = pd.get_dummies(labels, drop_first=True, dtype=float)
    return encoded.to_numpy(float), levels


def _empty_gene_score_table() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "row_id",
            "cell_id",
            "gene_key",
            "junction_deviation_z",
            "total_junction_reads",
            "measured_events",
            "sample_id",
            "donor_id",
            "cell_type",
            "condition",
        ]
    )


def _empty_multijunction_table() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "row_id",
            "cell_type",
            "cluster_id",
            "gene_id",
            "gene_symbol",
            "junction_choices",
            "event_ids",
            "strongest_event_id",
            "max_abs_delta_usage",
            "strongest_delta_usage",
            "numerator_donors",
            "denominator_donors",
            "overdispersion_rho",
            "likelihood_ratio",
            "degrees_of_freedom",
            "pvalue",
            "model_status",
            "batch_adjusted",
            "batch_levels",
            "fdr",
            "significant",
        ]
    )


def extract_junction_choices_from_bams(
    context: ProjectContext,
    output_path: Path,
    *,
    cell_barcode_tag: str = "CB",
    umi_tag: str = "UB",
    min_mapping_quality: int = 20,
    max_reads_per_bam: int = 0,
    check_cancelled: Callable[[], None] | None = None,
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

    from katana_next.science.molecule_store import MoleculeStore

    metadata = load_single_cell_metadata(context)
    full_length_assays = {"full_length", "smart_seq2", "smart_seq3", "smart-seq2", "smart-seq3"}
    count_modes = {}
    with MoleculeStore(output_path.parent) as molecules:
        for sample in context.manifest.samples:
            bam = sample.assets.get("single_cell_bam") or sample.assets.get("bam")
            if bam is None or not bam.is_file():
                continue
            sample_cells = metadata[metadata["sample_id"].eq(sample.sample_id)]
            if sample_cells.empty:
                raise ValueError(f"BAM sample {sample.sample_id} has no cell metadata.")
            metadata_cells = set(sample_cells["cell_id"])
            raw_lookup = {}
            if "barcode" in sample_cells:
                if sample_cells.barcode.isna().any() or sample_cells.barcode.duplicated().any():
                    raise ValueError(f"Barcode mapping must be unique and nonmissing within {sample.sample_id}")
                raw_lookup = dict(zip(sample_cells.barcode.astype(str), sample_cells.cell_id))
            full_length = sample_cells["assay"].str.lower().isin(full_length_assays).all()
            singleton = str(sample_cells.iloc[0]["cell_id"]) if full_length and len(sample_cells) == 1 else None
            declared_modes = set(sample_cells.get("count_mode", pd.Series(dtype=str)).dropna().astype(str))
            if len(declared_modes) > 1 or declared_modes - {"umi", "read_name"}:
                raise ValueError(f"Choose one explicit count_mode (umi or read_name) for {sample.sample_id}.")
            count_mode = next(iter(declared_modes), "read_name" if full_length else "umi")
            for cell_id in metadata_cells:
                count_modes[cell_id] = count_mode
            if backend == "pysam":
                handle = pysam.AlignmentFile(str(bam), "rb")
                records = handle.fetch(until_eof=True)
            else:
                import bamnostic

                index_candidates = (Path(str(bam) + ".bai"), bam.with_suffix(".bai"))
                index_path = next(
                    (path for path in index_candidates if path.is_file()), None
                )
                kwargs = {"filepath_index": str(index_path)} if index_path else {}
                handle = bamnostic.AlignmentFile(str(bam), "rb", **kwargs)
                records = iter(handle)
            try:
                for read_index, read in enumerate(records, start=1):
                    if check_cancelled and read_index % 20000 == 1:
                        check_cancelled()
                    if max_reads_per_bam and read_index > max_reads_per_bam:
                        break
                    if any(
                        bool(getattr(read, flag, False))
                        for flag in ("is_unmapped", "is_secondary", "is_supplementary", "is_qcfail")
                    ):
                        continue
                    if (
                        int(getattr(read, "mapping_quality", getattr(read, "mapq", 0)) or 0)
                        < min_mapping_quality
                    ):
                        continue
                    try:
                        raw_cell = str(read.get_tag(cell_barcode_tag))
                        prefixed_cell = f"{sample.sample_id}:{raw_cell}"
                        cell = raw_lookup.get(raw_cell, raw_cell if raw_cell in metadata_cells else prefixed_cell)
                    except (KeyError, AttributeError):
                        if singleton is None:
                            continue
                        cell = singleton
                    if cell not in metadata_cells:
                        continue
                    if count_mode == "umi":
                        try:
                            umi = str(read.get_tag(umi_tag))
                        except (KeyError, AttributeError):
                            continue
                    else:
                        umi = str(getattr(read, "query_name", "") or "")
                    if not umi:
                        continue
                    chrom = str(getattr(read, "reference_name", "") or "")
                    if not chrom and backend == "pysam":
                        chrom = str(handle.get_reference_name(read.reference_id))
                    reference = int(read.reference_start)
                    try:
                        strand = str(read.get_tag("XS"))
                    except (KeyError, AttributeError):
                        strand = "."
                    cigar = list(
                        getattr(read, "cigartuples", None)
                        or getattr(read, "cigar", [])
                        or []
                    )
                    for item in cigar:
                        operation, length = _normalize_cigar_item(item)
                        if operation == 3:
                            start, end = reference, reference + int(length)
                            molecules.add(
                                (sample.sample_id, cell, chrom, start, end, strand, umi)
                            )
                        if operation in {0, 2, 3, 7, 8}:
                            reference += int(length)
            finally:
                from katana_cosplice.bam_io import close_alignment
                close_alignment(handle)
        junctions = molecules.counts()
    if junctions.empty:
        raise ValueError("No splice-junction counts matched the assay, barcode, count mode and metadata. Droplet data require UMI tags; full-length single-cell BAMs can use read names.")
    choices = (
        junctions[["chrom", "start", "end", "strand"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    choices["cluster_id"] = _splice_site_components(choices)
    choices["event_id"] = (
        choices["chrom"]
        + ":"
        + choices["start"].astype(str)
        + "-"
        + choices["end"].astype(str)
        + ":"
        + choices["strand"]
    )
    choices = _annotate_junction_choices(choices, context.manifest.annotation.gtf)
    junctions = junctions.merge(
        choices,
        on=["chrom", "start", "end", "strand"],
        how="left",
        validate="many_to_one",
    )

    observed_clusters = (
        junctions.groupby(["sample_id", "cell_id", "cluster_id"], as_index=False)[
            "junction_count"
        ]
        .sum()
        .rename(columns={"junction_count": "cluster_total"})
    )
    expanded = observed_clusters.merge(
        choices[["cluster_id", "event_id", "gene_id", "gene_symbol", "annotated"]],
        on="cluster_id",
        how="left",
        validate="many_to_many",
    )
    inclusion = junctions[
        ["sample_id", "cell_id", "cluster_id", "event_id", "junction_count"]
    ]
    expanded = expanded.merge(
        inclusion,
        on=["sample_id", "cell_id", "cluster_id", "event_id"],
        how="left",
        validate="one_to_one",
    )
    expanded["junction_count"] = expanded["junction_count"].fillna(0).astype(int)
    output = pd.DataFrame(
        {
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
            "count_mode": expanded["cell_id"].map(count_modes),
        }
    )
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

    transcripts: dict[tuple[str, str, str], list[tuple[int, int, str, str]]] = {}
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
    return {key: content for key, content in re.findall(r'(\S+)\s+"([^"]*)"', value)}


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
        [
            f"{choices.iloc[index]['chrom']}:{choices.iloc[index]['strand']}:JC{labels[root]:06d}"
            for index, root in enumerate(roots)
        ],
        index=choices.index,
    )
