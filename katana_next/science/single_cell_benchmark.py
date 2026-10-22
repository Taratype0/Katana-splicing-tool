from __future__ import annotations

from dataclasses import dataclass
from typing import Any
import warnings

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
    silhouette_score,
)

from katana_cosplice.benchmark import fit_exact_wgcna_labels
from katana_cosplice.config import BenchmarkConfig, NetworkConfig
from katana_next.science.single_cell_programs import (
    fit_single_cell_tom_labels,
    smooth_cell_event_matrix,
)
from katana_next.science.statistics import safe_spearman
from katana_next.science.experimental_design import observation_units, bootstrap_blocks


@dataclass(frozen=True)
class SingleCellBenchmarkResult:
    quantification_summary: pd.DataFrame
    program_summary: pd.DataFrame
    program_assignments: pd.DataFrame
    bootstrap_stability: pd.DataFrame
    method_audit: pd.DataFrame


def validate_psi_observations(frame: pd.DataFrame, value: str, keys: list[str]) -> pd.DataFrame:
    """Keep missing PSI, reject malformed values or many-to-many accuracy pairs."""
    result = frame.copy()
    if result[keys].isna().any().any() or result[keys].astype(str).apply(lambda c: c.str.strip().eq("")).any().any():
        raise ValueError("PSI identity columns must not be empty")
    if result.duplicated(keys).any():
        raise ValueError(f"Duplicate PSI observations for {keys}")
    values = pd.to_numeric(result[value], errors="raise")
    if np.isinf(values).any() or ((values < 0) | (values > 1)).any():
        raise ValueError("PSI must be in [0, 1] or missing; gene splicing scores and percentages are not event PSI")
    result[value] = values
    return result


def benchmark_single_cell_methods(
    metadata: pd.DataFrame,
    counts: pd.DataFrame,
    pseudobulk: pd.DataFrame,
    *,
    embedding: pd.DataFrame | None = None,
    reference_psi: pd.DataFrame | None = None,
    reference_modules: pd.DataFrame | None = None,
    competitor_psi: pd.DataFrame | None = None,
    min_cell_event_reads: float = 1.0,
    min_donors: int = 15,
    min_unit_reads: float = 10.0,
    min_unit_fraction: float = 0.70,
    min_psi_variance: float = 0.0025,
    min_module_size: int = 5,
    max_events: int = 300,
    neighbors: int = 15,
    bootstrap_iterations: int = 10,
    soft_power: int = 6,
    random_seed: int = 17,
) -> SingleCellBenchmarkResult:
    """Compare quantification and program methods on explicitly shared universes."""

    quantification, quant_audit = _quantification_benchmark(
        metadata,
        counts,
        embedding=embedding,
        reference_psi=reference_psi,
        competitor_psi=competitor_psi,
        min_reads=min_cell_event_reads,
        neighbors=neighbors,
        random_seed=random_seed,
    )
    program, assignments, bootstrap, program_audit = _program_benchmark(
        pseudobulk,
        reference_modules=reference_modules,
        min_donors=min_donors,
        min_unit_reads=min_unit_reads,
        min_unit_fraction=min_unit_fraction,
        min_psi_variance=min_psi_variance,
        min_module_size=min_module_size,
        max_events=max_events,
        bootstrap_iterations=bootstrap_iterations,
        soft_power=soft_power,
        random_seed=random_seed,
    )
    audit = pd.concat([quant_audit, program_audit], ignore_index=True, sort=False)
    audit["row_id"] = audit["benchmark_track"] + "::" + audit["method"]
    return SingleCellBenchmarkResult(
        quantification_summary=quantification,
        program_summary=program,
        program_assignments=assignments,
        bootstrap_stability=bootstrap,
        method_audit=audit,
    )


def _quantification_benchmark(
    metadata: pd.DataFrame,
    counts: pd.DataFrame,
    *,
    embedding: pd.DataFrame | None,
    reference_psi: pd.DataFrame | None,
    competitor_psi: pd.DataFrame | None,
    min_reads: float,
    neighbors: int,
    random_seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    estimates = counts[
        ["cell_id", "event_id", "inclusion", "total", "psi_posterior"]
    ].copy()
    estimates["raw_psi"] = np.divide(
        estimates["inclusion"],
        estimates["total"],
        out=np.full(len(estimates), np.nan, dtype=float),
        where=estimates["total"].ge(min_reads),
    )
    estimates["measured_posterior_psi"] = estimates["psi_posterior"].where(
        estimates["total"].ge(min_reads)
    )
    long = pd.concat(
        [
            estimates[["cell_id", "event_id", "raw_psi"]]
            .rename(columns={"raw_psi": "estimated_psi"})
            .assign(method="Observed read fraction", uses_external_embedding=False),
            estimates[["cell_id", "event_id", "measured_posterior_psi"]]
            .rename(columns={"measured_posterior_psi": "estimated_psi"})
            .assign(method="Beta posterior mean", uses_external_embedding=False),
        ],
        ignore_index=True,
    )
    audit_rows = [
        {
            "benchmark_track": "PSI quantification",
            "method": "Observed read fraction",
            "status": "available",
            "implementation": "inclusion / total for measured cell-event pairs",
        },
        {
            "benchmark_track": "PSI quantification",
            "method": "Beta posterior mean",
            "status": "available",
            "implementation": "Jeffreys-prior Beta posterior mean",
        },
        {
            "benchmark_track": "Published method map",
            "method": "SICILIAN",
            "status": "reference_only",
            "implementation": "Statistical splice-junction confidence filtering after alignment.",
            "scientific_endpoint": "single-cell junction precision and replicate concordance",
            "published_statistic": "junction-level SICILIAN confidence model",
            "reference": "https://doi.org/10.1186/s13059-021-02434-8",
        },
        {
            "benchmark_track": "Published method map",
            "method": "SpliZ / SpliZVD",
            "status": "reference_only",
            "implementation": "Rank-based gene-cell splicing Z score and SVD-derived junction residual score; not an event-PSI estimator.",
            "scientific_endpoint": "cell-state-associated gene splicing",
            "published_statistic": "SpliZ or SpliZVD with permutation/FDR in the published workflow",
            "reference": "https://doi.org/10.1038/s41592-022-01400-x",
        },
        {
            "benchmark_track": "Published method map",
            "method": "Psix",
            "status": "reference_only",
            "implementation": "Probabilistic local cell-manifold association and co-regulated exon modules.",
            "scientific_endpoint": "exon association with continuous cell state",
            "published_statistic": "manifold autocorrelation-inspired Psix score",
            "reference": "https://doi.org/10.1101/gr.276109.121",
        },
        {
            "benchmark_track": "Published method map",
            "method": "BRIE2",
            "status": "reference_only",
            "implementation": "Bayesian hierarchical PSI inference with covariate-aware feature testing.",
            "scientific_endpoint": "cell-level exon inclusion and feature association",
            "published_statistic": "posterior PSI and ELBO gain",
            "reference": "https://doi.org/10.1186/s13059-021-02461-5",
        },
        {
            "benchmark_track": "Published method map",
            "method": "scQuint",
            "status": "reference_only",
            "implementation": "Alternative-intron groups with Dirichlet-multinomial latent and differential models.",
            "scientific_endpoint": "cell-type embedding and donor-aware differential junction usage",
            "published_statistic": "Dirichlet-multinomial GLM likelihood-ratio test with BH FDR",
            "reference": "https://doi.org/10.7554/eLife.73520",
        },
        {
            "benchmark_track": "Published method map",
            "method": "MARVEL",
            "status": "reference_only",
            "implementation": "Integrated event validation, PSI modality, differential-splicing and functional-annotation workflow.",
            "scientific_endpoint": "end-to-end single-cell alternative-splicing interpretation",
            "published_statistic": "bimodality correction, Anderson-Darling/DTS and pseudobulk permutation tracks",
            "reference": "https://doi.org/10.1093/nar/gkac1260",
        },
        {
            "benchmark_track": "Published method map",
            "method": "SCSES / DOLPHIN",
            "status": "reference_only",
            "implementation": "Recent network-diffusion or exon-junction graph recovery/representation methods.",
            "scientific_endpoint": "sparse signal recovery or splice-aware cell representation",
            "published_statistic": "reference correlation/RMSE or downstream representation metrics",
            "reference": "https://doi.org/10.1038/s41467-025-64517-5; https://doi.org/10.1038/s41467-025-61580-w",
        },
    ]
    if embedding is not None:
        shared_cells = pd.Index(metadata["cell_id"].astype(str)).intersection(
            embedding["cell_id"].astype(str)
        )
        matrix = estimates.pivot_table(
            index="event_id", columns="cell_id", values="measured_posterior_psi", aggfunc="first"
        ).reindex(columns=shared_cells)
        coordinates = (
            embedding.set_index("cell_id")
            .reindex(shared_cells)[[c for c in embedding if c.startswith("state_")]]
            .to_numpy(float)
        )
        smoothed = smooth_cell_event_matrix(matrix, coordinates, neighbors)
        smooth_long = (
            smoothed.rename_axis(index="event_id", columns="cell_id")
            .stack(future_stack=True)
            .rename("estimated_psi")
            .reset_index()
        )
        smooth_long["method"] = "External-manifold neighbour mean"
        smooth_long["uses_external_embedding"] = True
        long = pd.concat([long, smooth_long], ignore_index=True, sort=False)
        audit_rows.append(
            {
                "benchmark_track": "PSI quantification",
                "method": "External-manifold neighbour mean",
                "status": "available",
                "implementation": "count-informed PSI smoothed over the declared expression/state kNN graph",
            }
        )
    if competitor_psi is not None and not competitor_psi.empty:
        required = {"method", "cell_id", "event_id", "estimated_psi"}
        missing = sorted(required.difference(competitor_psi.columns))
        if missing:
            raise ValueError(f"Competitor PSI input is missing columns: {missing}")
        imported = competitor_psi[
            ["method", "cell_id", "event_id", "estimated_psi"]
        ].copy()
        imported = validate_psi_observations(imported, "estimated_psi", ["method", "cell_id", "event_id"])
        imported["method"] = "Imported: " + imported["method"].astype(str)
        imported["uses_external_embedding"] = np.nan
        long = pd.concat([long, imported], ignore_index=True, sort=False)
        for method in sorted(imported["method"].unique()):
            audit_rows.append(
                {
                    "benchmark_track": "PSI quantification",
                    "method": method,
                    "status": "imported",
                    "implementation": "user-declared external method output; Katana did not run the competitor",
                }
            )

    if reference_psi is None or reference_psi.empty:
        audit_rows.append(
            {
                "benchmark_track": "PSI quantification",
                "method": "Reference truth",
                "status": "unavailable",
                "implementation": "Declare single_cell_reference_psi for accuracy metrics.",
            }
        )
        return pd.DataFrame(), pd.DataFrame(audit_rows)
    reference = reference_psi.copy()
    psi_column = next(
        (
            column
            for column in ("reference_psi", "true_psi", "psi")
            if column in reference
        ),
        None,
    )
    if psi_column is None or not {"cell_id", "event_id"}.issubset(reference.columns):
        raise ValueError(
            "Reference PSI requires cell_id, event_id and reference_psi/true_psi/psi."
        )
    reference = reference[["cell_id", "event_id", psi_column]].rename(
        columns={psi_column: "reference_psi"}
    )
    reference = validate_psi_observations(reference, "reference_psi", ["cell_id", "event_id"])
    joined = long.merge(
        reference, on=["cell_id", "event_id"], how="inner", validate="many_to_one"
    )
    total_reference = len(reference.dropna(subset=["reference_psi"]))
    rows: list[dict[str, Any]] = []
    for method, group in joined.groupby("method", sort=True):
        valid = group.dropna(subset=["estimated_psi", "reference_psi"])
        correlation, correlation_pvalue = safe_spearman(
            valid["estimated_psi"], valid["reference_psi"]
        )
        rmse = (
            float(
                np.sqrt(np.mean((valid["estimated_psi"] - valid["reference_psi"]) ** 2))
            )
            if len(valid)
            else np.nan
        )
        clustering = _cell_clustering_metrics(
            group,
            metadata,
            random_seed=random_seed,
        )
        rows.append(
            {
                "row_id": method,
                "method": method,
                "reference_pairs": total_reference,
                "evaluated_pairs": len(valid),
                "coverage_fraction": len(valid) / total_reference
                if total_reference
                else np.nan,
                "rmse": rmse,
                "spearman": correlation,
                "spearman_pvalue": correlation_pvalue,
                "cell_type_ari": clustering["ari"],
                "cell_type_nmi": clustering["nmi"],
                "cell_embedding_silhouette": clustering["silhouette"],
                "uses_external_embedding": bool(
                    group["uses_external_embedding"].fillna(False).any()
                ),
                "interpretation": "Accuracy metrics use only exact shared cell-event pairs.",
            }
        )
    return pd.DataFrame(rows).sort_values(["rmse", "method"]), pd.DataFrame(audit_rows)


def _cell_clustering_metrics(
    values: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    random_seed: int,
) -> dict[str, float]:
    matrix = values.pivot_table(
        index="cell_id", columns="event_id", values="estimated_psi", aggfunc="first"
    )
    if matrix.shape[0] < 4 or matrix.shape[1] < 2:
        return {"ari": np.nan, "nmi": np.nan, "silhouette": np.nan}
    keep = matrix.notna().mean(axis=0).ge(0.20) & matrix.var(axis=0).gt(0)
    matrix = matrix.loc[:, keep]
    if matrix.shape[1] < 2:
        return {"ari": np.nan, "nmi": np.nan, "silhouette": np.nan}
    matrix = matrix.fillna(matrix.median(axis=0))
    labels = metadata.set_index("cell_id").reindex(matrix.index)["cell_type"]
    valid = labels.notna()
    matrix, labels = matrix.loc[valid], labels.loc[valid].astype(str)
    clusters = labels.nunique()
    if len(matrix) <= clusters or clusters < 2:
        return {"ari": np.nan, "nmi": np.nan, "silhouette": np.nan}
    components = min(10, len(matrix) - 1, matrix.shape[1])
    latent = PCA(n_components=components, random_state=random_seed).fit_transform(
        matrix
    )
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="KMeans is known to have a memory leak on Windows with MKL",
            category=UserWarning,
        )
        predicted = KMeans(
            n_clusters=clusters, n_init=20, random_state=random_seed
        ).fit_predict(latent)
    return {
        "ari": float(adjusted_rand_score(labels, predicted)),
        "nmi": float(normalized_mutual_info_score(labels, predicted)),
        "silhouette": float(silhouette_score(latent, labels)),
    }


def _program_benchmark(
    pseudobulk: pd.DataFrame,
    *,
    reference_modules: pd.DataFrame | None,
    min_donors: int,
    min_unit_reads: float,
    min_unit_fraction: float,
    min_psi_variance: float,
    min_module_size: int,
    max_events: int,
    bootstrap_iterations: int,
    soft_power: int,
    random_seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    summaries: list[dict[str, Any]] = []
    assignments: list[pd.DataFrame] = []
    bootstrap_rows: list[dict[str, Any]] = []
    audit_rows: list[dict[str, Any]] = []
    for cell_type, frame in pseudobulk.groupby("cell_type", sort=True):
        frame = observation_units(frame)
        donor_metadata = frame.drop_duplicates("observation_id").set_index("observation_id")
        donors = donor_metadata.index.astype(str)
        if donor_metadata["donor_id"].nunique() < min_donors:
            audit_rows.append(
                {
                    "benchmark_track": f"Replicate programs: {cell_type}",
                    "method": "All methods",
                    "status": "not_run",
                    "implementation": f"{len(donors)} independent donors; {min_donors} required",
                }
            )
            continue
        matrix, reliability = _shared_program_matrix(
            frame,
            donors,
            min_unit_reads=min_unit_reads,
            min_unit_fraction=min_unit_fraction,
            min_psi_variance=min_psi_variance,
            max_events=max_events,
        )
        if len(matrix) < min_module_size:
            continue
        network = NetworkConfig(
            method="bicor",
            signed=True,
            soft_power=soft_power,
            min_module_size=min_module_size,
            max_events=len(matrix),
        )
        benchmark_settings = BenchmarkConfig(
            max_events=len(matrix),
            bootstrap_iterations=bootstrap_iterations,
            wgcna_merge_cut_height=0.25,
            random_seed=random_seed,
        )
        methods: dict[str, tuple[np.ndarray | None, str]] = {}
        weighted, _tom, weighted_engine = fit_single_cell_tom_labels(
            matrix,
            reliability=reliability,
            min_module_size=min_module_size,
            soft_power=soft_power,
        )
        methods["Katana uncertainty-weighted TOM"] = (weighted, weighted_engine)
        unweighted, _tom, unweighted_engine = fit_single_cell_tom_labels(
            matrix,
            reliability=np.ones(len(matrix), dtype=float),
            min_module_size=min_module_size,
            soft_power=soft_power,
        )
        methods["Unweighted same R TOM"] = (unweighted, unweighted_engine)
        try:
            exact, exact_audit = fit_exact_wgcna_labels(
                matrix, network, benchmark_settings
            )
            methods["Exact R WGCNA"] = (exact, exact_audit)
        except (FileNotFoundError, RuntimeError, ValueError) as error:
            methods["Exact R WGCNA"] = (None, str(error))

        rng = np.random.default_rng(random_seed)
        bootstrap_indices = [
            bootstrap_blocks(donor_metadata.reindex(donors), rng)
            for _ in range(bootstrap_iterations)
        ]
        for method, (labels, implementation) in methods.items():
            if labels is None:
                audit_rows.append(
                    {
                        "benchmark_track": f"Replicate programs: {cell_type}",
                        "method": method,
                        "status": "unavailable",
                        "implementation": implementation,
                    }
                )
                continue
            label_series = pd.Series(labels, index=matrix.index)
            metrics = _partition_metrics(matrix, labels)
            reference_metrics = _reference_partition_metrics(
                label_series,
                reference_modules,
                cell_type=str(cell_type),
            )
            run_ari: list[float] = []
            for iteration, indices in enumerate(bootstrap_indices, 1):
                boot = matrix.iloc[:, indices].copy()
                boot.columns = [
                    f"boot_{iteration}_{index}" for index in range(len(indices))
                ]
                try:
                    if method == "Exact R WGCNA":
                        boot_labels, _audit = fit_exact_wgcna_labels(
                            boot, network, benchmark_settings
                        )
                        if boot_labels is None:
                            raise RuntimeError(_audit)
                    else:
                        boot_labels, _boot_tom, _engine = fit_single_cell_tom_labels(
                            boot,
                            reliability=reliability
                            if method == "Katana uncertainty-weighted TOM"
                            else np.ones(len(boot), dtype=float),
                            min_module_size=min_module_size,
                            soft_power=soft_power,
                        )
                    ari = float(adjusted_rand_score(labels, boot_labels))
                    run_ari.append(ari)
                    bootstrap_rows.append(
                        {
                            "row_id": f"{cell_type}::{method}::{iteration}",
                            "cell_type": cell_type,
                            "method": method,
                            "iteration": iteration,
                            "adjusted_rand_index": ari,
                            "status": "complete",
                        }
                    )
                except (RuntimeError, ValueError, FileNotFoundError) as error:
                    bootstrap_rows.append(
                        {
                            "row_id": f"{cell_type}::{method}::{iteration}",
                            "cell_type": cell_type,
                            "method": method,
                            "iteration": iteration,
                            "adjusted_rand_index": np.nan,
                            "status": "failed",
                            "error": str(error),
                        }
                    )
            summaries.append(
                {
                    "row_id": f"{cell_type}::{method}",
                    "cell_type": cell_type,
                    "method": method,
                    "events": len(matrix),
                    "donors": donor_metadata["donor_id"].nunique(),
                    "donor_condition_observations": len(donors),
                    **metrics,
                    **reference_metrics,
                    "mean_bootstrap_ari": float(np.mean(run_ari))
                    if run_ari
                    else np.nan,
                    "min_bootstrap_ari": float(np.min(run_ari)) if run_ari else np.nan,
                    "interpretation": "Methods use the identical event and donor universe; no composite winner is manufactured.",
                }
            )
            table = pd.DataFrame(
                {
                    "row_id": [
                        f"{cell_type}::{method}::{event}" for event in matrix.index
                    ],
                    "cell_type": cell_type,
                    "method": method,
                    "event_id": matrix.index,
                    "module": label_series.to_numpy(),
                }
            )
            assignments.append(table)
            audit_rows.append(
                {
                    "benchmark_track": f"Replicate programs: {cell_type}",
                    "method": method,
                    "status": "complete",
                    "implementation": implementation,
                }
            )
    return (
        pd.DataFrame(summaries),
        pd.concat(assignments, ignore_index=True) if assignments else pd.DataFrame(),
        pd.DataFrame(bootstrap_rows),
        pd.DataFrame(audit_rows),
    )


def _shared_program_matrix(
    frame: pd.DataFrame,
    donors: pd.Index,
    *,
    min_unit_reads: float,
    min_unit_fraction: float,
    min_psi_variance: float,
    max_events: int,
) -> tuple[pd.DataFrame, np.ndarray]:
    psi = frame.pivot_table(
        index="event_id", columns="observation_id", values="psi_posterior", aggfunc="first"
    ).reindex(columns=donors)
    total = frame.pivot_table(
        index="event_id", columns="observation_id", values="total", aggfunc="first"
    ).reindex(index=psi.index, columns=donors)
    se = frame.pivot_table(
        index="event_id", columns="observation_id", values="measurement_se", aggfunc="first"
    ).reindex(index=psi.index, columns=donors)
    observed = psi.where(total.ge(min_unit_reads))
    valid = observed.notna().mean(axis=1)
    variance = observed.var(axis=1, ddof=1)
    median_reads = total.where(total.ge(min_unit_reads)).median(axis=1)
    median_se = se.where(total.ge(min_unit_reads)).median(axis=1)
    reliability = (
        valid
        * median_reads.div(median_reads + 20.0)
        * (1.0 - 2.0 * median_se).clip(0.0, 1.0)
    )
    event_table = pd.DataFrame(
        {"valid": valid, "variance": variance, "reliability": reliability}
    )
    selected = (
        event_table.query(
            "valid >= @min_unit_fraction and variance >= @min_psi_variance"
        )
        .assign(score=lambda table: table["variance"] * table["reliability"])
        .sort_values("score", ascending=False)
        .head(max_events)
    )
    matrix = observed.reindex(selected.index)
    matrix = matrix.T.fillna(matrix.median(axis=1)).T
    return matrix, selected["reliability"].fillna(0.0).to_numpy(float)


def _partition_metrics(matrix: pd.DataFrame, labels: np.ndarray) -> dict[str, Any]:
    label_series = pd.Series(labels, index=matrix.index)
    retained = label_series.ne("grey")
    sizes = label_series.loc[retained].value_counts()
    n_modules = len(sizes)
    silhouette = np.nan
    if n_modules >= 2 and retained.sum() > n_modules:
        silhouette = float(
            silhouette_score(
                matrix.loc[retained], label_series.loc[retained], metric="cosine"
            )
        )
    return {
        "modules": n_modules,
        "grey_fraction": float(label_series.eq("grey").mean()),
        "largest_module_fraction": float(sizes.max() / len(matrix))
        if not sizes.empty
        else 0.0,
        "silhouette_cosine": silhouette,
        "valid_partition": bool(
            n_modules >= 2
            and not sizes.empty
            and float(sizes.max() / len(matrix)) <= 0.80
        ),
    }


def _reference_partition_metrics(
    labels: pd.Series,
    reference: pd.DataFrame | None,
    *,
    cell_type: str,
) -> dict[str, float | int]:
    empty = {
        "reference_events": 0,
        "reference_ari": np.nan,
        "reference_nmi": np.nan,
    }
    if reference is None or reference.empty:
        return empty
    required = {"event_id", "reference_module"}
    if not required.issubset(reference.columns):
        raise ValueError(
            "Reference modules require event_id and reference_module columns."
        )
    selected = reference.copy()
    if "cell_type" in selected:
        selected = selected[
            selected["cell_type"].astype(str).isin([cell_type, "all", "*"])
        ]
    selected = selected.drop_duplicates("event_id").set_index("event_id")
    selected.index = selected.index.astype(str)
    shared = labels.index.astype(str).intersection(selected.index)
    if len(shared) < 2:
        return empty
    observed = labels.reindex(shared).astype(str)
    expected = selected.reindex(shared)["reference_module"].astype(str)
    return {
        "reference_events": len(shared),
        "reference_ari": float(adjusted_rand_score(expected, observed)),
        "reference_nmi": float(normalized_mutual_info_score(expected, observed)),
    }
