from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.stats import norm
from sklearn.metrics import roc_auc_score

from .bayesian_preservation import bayesian_multi_cohort_preservation
from .differential_network import compare_condition_networks
from .multiview import fuse_event_edge_views
from .uncertainty import estimate_beta_binomial_uncertainty


@dataclass
class MethodExtensionBenchmarkResult:
    synthetic_metrics: pd.DataFrame
    real_count_metrics: pd.DataFrame
    public_preservation_metrics: pd.DataFrame
    recommendation: pd.DataFrame


def _pair_scores(matrix: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    correlation = np.corrcoef(matrix.to_numpy(dtype=float))
    upper = np.triu_indices_from(correlation, k=1)
    return upper[0], upper[1], np.abs(correlation[upper])


def _top_edge_set(scores: np.ndarray, top_k: int) -> set[int]:
    top_k = min(max(int(top_k), 1), len(scores))
    return set(np.argpartition(scores, -top_k)[-top_k:].tolist())


def _jaccard(left: set[int], right: set[int]) -> float:
    union = left | right
    return len(left & right) / len(union) if union else 1.0


def _synthetic_uncertainty(seed: int) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    rng = np.random.default_rng(seed)
    n_modules = 3
    events_per_module = 36
    n_samples = 24
    event_modules = np.repeat(np.arange(n_modules), events_per_module)
    sample_ids = [f"S{index + 1:02d}" for index in range(n_samples)]
    event_ids = [f"E{index + 1:03d}" for index in range(len(event_modules))]
    latent = rng.normal(size=(n_modules, n_samples))
    true_psi = np.vstack(
        [
            expit(
                -0.1
                + 1.25 * latent[module]
                + rng.normal(0.0, 0.12, size=n_samples)
            )
            for module in event_modules
        ]
    )
    high_depth = np.tile(np.arange(events_per_module) < 20, n_modules)
    depth_mean = np.where(high_depth, 140.0, 7.0)
    totals = rng.poisson(depth_mean[:, None], size=true_psi.shape) + 1
    concentration = np.where(high_depth, 120.0, 8.0)[:, None]
    observed_probability = rng.beta(
        np.clip(true_psi * concentration, 0.05, None),
        np.clip((1.0 - true_psi) * concentration, 0.05, None),
    )
    inclusion = rng.binomial(totals, observed_probability)
    skipping = totals - inclusion
    ijc = pd.DataFrame(inclusion, index=event_ids, columns=sample_ids)
    sjc = pd.DataFrame(skipping, index=event_ids, columns=sample_ids)
    raw = ijc / (ijc + sjc)
    samples = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "condition": np.repeat(["A", "B", "C", "D"], 6),
        }
    )
    qc = pd.DataFrame({"stability_score": 1.0}, index=event_ids)
    uncertainty = estimate_beta_binomial_uncertainty(ijc, sjc, samples, qc)
    rows, cols, raw_scores = _pair_scores(raw)
    _, _, posterior_scores = _pair_scores(uncertainty.posterior_psi)
    reliability = uncertainty.event_metrics["network_reliability"].to_numpy()
    reliability_factor = np.sqrt(reliability[rows] * reliability[cols])
    raw_weighted_scores = raw_scores * reliability_factor
    weighted_scores = posterior_scores * reliability_factor
    truth = (event_modules[rows] == event_modules[cols]).astype(int)
    recoverable_truth = (
        truth.astype(bool) & high_depth[rows] & high_depth[cols]
    ).astype(int)
    top_k = int(recoverable_truth.sum())
    metrics = pd.DataFrame(
        [
            {
                "algorithm": "uncertainty_weighting",
                "configuration": "raw_PSI",
                "metric": "recoverable_edge_AUROC",
                "value": roc_auc_score(recoverable_truth, raw_scores),
                "higher_is_better": True,
            },
            {
                "algorithm": "uncertainty_weighting",
                "configuration": "raw_PSI_reliability_weighted",
                "metric": "recoverable_edge_AUROC",
                "value": roc_auc_score(recoverable_truth, raw_weighted_scores),
                "higher_is_better": True,
            },
            {
                "algorithm": "uncertainty_weighting",
                "configuration": "beta_binomial_posterior",
                "metric": "recoverable_edge_AUROC",
                "value": roc_auc_score(recoverable_truth, posterior_scores),
                "higher_is_better": True,
            },
            {
                "algorithm": "uncertainty_weighting",
                "configuration": "beta_binomial_posterior_weighted",
                "metric": "recoverable_edge_AUROC",
                "value": roc_auc_score(recoverable_truth, weighted_scores),
                "higher_is_better": True,
            },
            {
                "algorithm": "uncertainty_weighting",
                "configuration": "raw_PSI",
                "metric": "recoverable_top_edge_precision",
                "value": float(
                    recoverable_truth[list(_top_edge_set(raw_scores, top_k))].mean()
                ),
                "higher_is_better": True,
            },
            {
                "algorithm": "uncertainty_weighting",
                "configuration": "raw_PSI_reliability_weighted",
                "metric": "recoverable_top_edge_precision",
                "value": float(
                    recoverable_truth[
                        list(_top_edge_set(raw_weighted_scores, top_k))
                    ].mean()
                ),
                "higher_is_better": True,
            },
            {
                "algorithm": "uncertainty_weighting",
                "configuration": "beta_binomial_posterior",
                "metric": "recoverable_top_edge_precision",
                "value": float(
                    recoverable_truth[
                        list(_top_edge_set(posterior_scores, top_k))
                    ].mean()
                ),
                "higher_is_better": True,
            },
            {
                "algorithm": "uncertainty_weighting",
                "configuration": "beta_binomial_posterior_weighted",
                "metric": "recoverable_top_edge_precision",
                "value": float(
                    recoverable_truth[
                        list(_top_edge_set(weighted_scores, top_k))
                    ].mean()
                ),
                "higher_is_better": True,
            },
        ]
    )
    edge_table = pd.DataFrame(
        {
            "event_a": np.asarray(event_ids)[rows],
            "event_b": np.asarray(event_ids)[cols],
            "adjacency": raw_scores,
            "truth": truth,
        }
    )
    event_features = pd.DataFrame(
        {
            "event_key": event_ids,
            "regulator_axis": event_modules * 6.0 + rng.normal(0.0, 0.2, len(event_ids)),
            "adversarial_noise": rng.normal(size=len(event_ids)),
        }
    )
    return metrics, {
        "edge_table": edge_table,
        "event_features": event_features,
        "raw": raw,
        "posterior": uncertainty.posterior_psi,
        "reliability": uncertainty.event_metrics,
    }


def _synthetic_multiview(
    inputs: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    edges = inputs["edge_table"]
    features = inputs["event_features"]
    informative = fuse_event_edge_views(
        edges,
        features,
        auxiliary_columns=["regulator_axis"],
        as_weight=0.80,
        edge_threshold=0.0,
        min_module_size=10,
    ).fused_edges
    adversarial = fuse_event_edge_views(
        edges,
        features,
        auxiliary_columns=["adversarial_noise"],
        as_weight=0.80,
        edge_threshold=0.0,
        min_module_size=10,
    ).fused_edges
    truth_lookup = edges.set_index(["event_a", "event_b"])["truth"]
    informative_truth = informative.set_index(["event_a", "event_b"]).index.map(
        truth_lookup
    ).to_numpy(dtype=int)
    adversarial_truth = adversarial.set_index(["event_a", "event_b"]).index.map(
        truth_lookup
    ).to_numpy(dtype=int)
    truth = edges["truth"].to_numpy(dtype=int)
    base_auc = roc_auc_score(truth, edges["adjacency"])
    return pd.DataFrame(
        [
            {
                "algorithm": "multiview_AS_primary",
                "configuration": "AS_only",
                "metric": "true_edge_AUROC",
                "value": base_auc,
                "higher_is_better": True,
            },
            {
                "algorithm": "multiview_AS_primary",
                "configuration": "informative_auxiliary_20pct",
                "metric": "true_edge_AUROC",
                "value": roc_auc_score(
                    informative_truth, informative["fused_adjacency"]
                ),
                "higher_is_better": True,
            },
            {
                "algorithm": "multiview_AS_primary",
                "configuration": "adversarial_auxiliary_20pct",
                "metric": "true_edge_AUROC",
                "value": roc_auc_score(
                    adversarial_truth, adversarial["fused_adjacency"]
                ),
                "higher_is_better": True,
            },
            {
                "algorithm": "multiview_AS_primary",
                "configuration": "informative_auxiliary_20pct",
                "metric": "minimum_AS_primary_weight",
                "value": float(informative["AS_primary_weight"].min()),
                "higher_is_better": True,
            },
        ]
    )


def _synthetic_differential(seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed + 1)
    n_per_condition = 12
    sample_ids = [f"A{index + 1}" for index in range(n_per_condition)] + [
        f"B{index + 1}" for index in range(n_per_condition)
    ]
    latent_a = rng.normal(size=n_per_condition)
    latent_b = rng.normal(size=n_per_condition)
    gained_a = rng.normal(size=(10, n_per_condition))
    gained_b = np.vstack(
        [latent_b + rng.normal(0.0, 0.08, n_per_condition) for _ in range(10)]
    )
    stable_a = np.vstack(
        [latent_a + rng.normal(0.0, 0.08, n_per_condition) for _ in range(10)]
    )
    # The null module uses the same sample-level topology on both sides. This
    # makes it a strict false-positive guardrail rather than another noisy
    # alternative.
    stable_b = stable_a.copy()
    values = np.vstack(
        [
            np.hstack([gained_a, gained_b]),
            np.hstack([stable_a, stable_b]),
        ]
    )
    event_keys = [f"D{index + 1:02d}" for index in range(len(values))]
    matrix = pd.DataFrame(values, index=event_keys, columns=sample_ids)
    samples = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "condition": ["A"] * n_per_condition + ["B"] * n_per_condition,
        }
    )
    assignments = pd.DataFrame(
        {
            "event_key": event_keys,
            "consensus_module": ["gained"] * 10 + ["stable"] * 10,
        }
    )
    result = compare_condition_networks(
        matrix,
        samples,
        assignments,
        permutations=199,
        max_output_edges_per_module_pair=100,
        random_seed=seed,
    )
    states = result.module_changes.set_index("module")
    gained_detected = float(
        states.loc["gained", "density_change_direction"] == "gained_in_B"
    )
    stable_false_positive = float(
        states.loc["stable", "differential_network_state"]
        != "stable_or_underpowered"
    )
    return pd.DataFrame(
        [
            {
                "algorithm": "differential_network",
                "configuration": "permutation_density_and_topology",
                "metric": "known_gained_module_detected",
                "value": gained_detected,
                "higher_is_better": True,
            },
            {
                "algorithm": "differential_network",
                "configuration": "permutation_density_and_topology",
                "metric": "stable_module_false_positive",
                "value": stable_false_positive,
                "higher_is_better": False,
            },
        ]
    )


def _synthetic_bayesian(seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed + 2)
    rows: list[dict[str, object]] = []
    profiles = {
        "preserved": (12.0, 11.0, 0.85, 0.82),
        "partial": (4.0, 3.5, 0.55, 0.50),
        "null": (0.0, 0.0, 0.05, 0.02),
    }
    for module, profile in profiles.items():
        for cohort in ("C1", "C2", "C3", "C4"):
            zdensity = profile[0] + rng.normal(0.0, 0.6)
            zconnectivity = profile[1] + rng.normal(0.0, 0.6)
            rows.append(
                {
                    "consensus_module": module,
                    "cohort": cohort,
                    "Zdensity": zdensity,
                    "Zconnectivity": zconnectivity,
                    "density_empirical_pvalue": 2.0 * norm.sf(abs(zdensity)),
                    "connectivity_empirical_pvalue": 2.0
                    * norm.sf(abs(zconnectivity)),
                    "program_score_preservation": np.clip(
                        profile[2] + rng.normal(0.0, 0.03), -0.99, 0.99
                    ),
                    "membership_profile_correlation": np.clip(
                        profile[3] + rng.normal(0.0, 0.03), -0.99, 0.99
                    ),
                }
            )
    posterior = bayesian_multi_cohort_preservation(pd.DataFrame(rows))
    summary = posterior.module_posteriors.set_index("consensus_module")
    return pd.DataFrame(
        [
            {
                "algorithm": "bayesian_preservation",
                "configuration": "hierarchical_multi_cohort",
                "metric": "preserved_module_probability",
                "value": summary.loc[
                    "preserved", "conservative_probability_partial_preservation"
                ],
                "higher_is_better": True,
            },
            {
                "algorithm": "bayesian_preservation",
                "configuration": "hierarchical_multi_cohort",
                "metric": "null_module_probability",
                "value": summary.loc[
                    "null", "conservative_probability_partial_preservation"
                ],
                "higher_is_better": False,
            },
        ]
    )


def run_synthetic_method_ablation(seed: int = 41) -> pd.DataFrame:
    uncertainty, inputs = _synthetic_uncertainty(seed)
    return pd.concat(
        [
            uncertainty,
            _synthetic_multiview(inputs),
            _synthetic_differential(seed),
            _synthetic_bayesian(seed),
        ],
        ignore_index=True,
    )


def benchmark_real_read_uncertainty(
    cosplice_output: Path,
    *,
    max_events: int = 600,
    bootstrap_iterations: int = 20,
    top_edges: int = 3_000,
    seed: int = 91,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    qc_path = cosplice_output / "03_stable_as_qc" / "stable_as_event_qc.tsv"
    psi_path = cosplice_output / "03_stable_as_qc" / "stable_psi_matrix.tsv"
    ijc_path = cosplice_output / "02_matrix" / "sample_level_ijc_matrix.tsv"
    sjc_path = cosplice_output / "02_matrix" / "sample_level_sjc_matrix.tsv"
    samples_path = cosplice_output / "02_matrix" / "sample_metadata_resolved.tsv"
    required = [qc_path, psi_path, ijc_path, sjc_path, samples_path]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing real-count benchmark inputs: {missing}")
    qc = pd.read_csv(qc_path, sep="\t", index_col=0, low_memory=False)
    selected = (
        pd.to_numeric(qc["stability_score"], errors="coerce")
        .sort_values(ascending=False)
        .head(max_events)
        .index
    )
    raw = pd.read_csv(psi_path, sep="\t", index_col=0).reindex(selected)
    ijc = pd.read_csv(ijc_path, sep="\t", index_col=0).reindex(selected)
    sjc = pd.read_csv(sjc_path, sep="\t", index_col=0).reindex(selected)
    samples = pd.read_csv(samples_path, sep="\t", low_memory=False)
    uncertainty = estimate_beta_binomial_uncertainty(
        ijc,
        sjc,
        samples,
        qc.reindex(selected),
    )
    reliability = uncertainty.event_metrics["network_reliability"].to_numpy()
    rows, cols, raw_full = _pair_scores(raw)
    _, _, posterior_full = _pair_scores(uncertainty.posterior_psi)
    reliability_factor = np.sqrt(reliability[rows] * reliability[cols])
    full_scores = {
        "raw_PSI": raw_full,
        "raw_PSI_reliability_weighted": raw_full * reliability_factor,
        "beta_binomial_posterior": posterior_full,
        "beta_binomial_posterior_weighted": posterior_full * reliability_factor,
    }
    references = {
        name: _top_edge_set(score, top_edges)
        for name, score in full_scores.items()
    }
    metadata = samples.set_index("sample_id").reindex(raw.columns)
    condition = metadata["condition"].astype(str)
    groups = {
        name: np.flatnonzero(condition.to_numpy() == name)
        for name in condition.unique()
    }
    rng = np.random.default_rng(seed)
    jaccards: dict[str, list[float]] = {name: [] for name in full_scores}
    for _ in range(bootstrap_iterations):
        sampled = np.concatenate(
            [rng.choice(indices, len(indices), replace=True) for indices in groups.values()]
        )
        _, _, raw_score = _pair_scores(raw.iloc[:, sampled])
        _, _, posterior_score = _pair_scores(
            uncertainty.posterior_psi.iloc[:, sampled]
        )
        bootstrap_scores = {
            "raw_PSI": raw_score,
            "raw_PSI_reliability_weighted": raw_score * reliability_factor,
            "beta_binomial_posterior": posterior_score,
            "beta_binomial_posterior_weighted": posterior_score
            * reliability_factor,
        }
        for name, score in bootstrap_scores.items():
            jaccards[name].append(
                _jaccard(references[name], _top_edge_set(score, top_edges))
            )
    metrics = pd.DataFrame(
        [
            {
                "dataset": cosplice_output.parent.name,
                "algorithm": name,
                "metric": "bootstrap_top_edge_Jaccard",
                "mean": float(np.mean(values)),
                "median": float(np.median(values)),
                "standard_deviation": float(np.std(values, ddof=1)),
                "n_iterations": bootstrap_iterations,
                "n_events": len(selected),
            }
            for name, values in jaccards.items()
        ]
    )
    audit = uncertainty.event_metrics.reset_index()
    return metrics, audit


def benchmark_public_bayesian_preservation(
    preservation_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    preservation = pd.read_csv(preservation_path, sep="\t", low_memory=False)
    module_columns = tuple(
        column
        for column in ("method", "discovery_cohort", "consensus_module")
        if column in preservation
    )
    result = bayesian_multi_cohort_preservation(
        preservation,
        module_columns=module_columns,
    )
    fixed = preservation.assign(
        fixed_preserved=preservation["preservation_class"].eq("preserved"),
        fixed_partial_or_better=preservation["preservation_class"].isin(
            ["preserved", "partially_preserved"]
        ),
    )
    fixed_summary = (
        fixed.groupby("method", observed=True)
        .agg(
            n_projection_tests=("preservation_class", "size"),
            fixed_preserved_fraction=("fixed_preserved", "mean"),
            fixed_partial_or_better_fraction=("fixed_partial_or_better", "mean"),
        )
        .reset_index()
    )
    bayesian_summary = (
        result.module_posteriors.groupby("method", observed=True)
        .agg(
            n_discovery_modules=("bayesian_preservation_class", "size"),
            bayesian_preserved_fraction=(
                "bayesian_preservation_class",
                lambda values: values.eq("preserved").mean(),
            ),
            bayesian_partial_or_better_fraction=(
                "bayesian_preservation_class",
                lambda values: values.isin(["preserved", "partially_preserved"]).mean(),
            ),
            median_partial_posterior=(
                "conservative_probability_partial_preservation",
                "median",
            ),
        )
        .reset_index()
    )
    return fixed_summary.merge(bayesian_summary, on="method", how="outer"), (
        result.module_posteriors
    )


def benchmark_real_differential_network(
    cosplice_output: Path,
    *,
    events_per_module: int = 40,
    permutations: int = 99,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    matrix_path = cosplice_output / "04_network" / "observed" / "selected_event_matrix.tsv"
    samples_path = cosplice_output / "02_matrix" / "sample_metadata_resolved.tsv"
    consensus_path = (
        cosplice_output
        / "04c_consensus_modules"
        / "consensus_event_module_assignments.tsv"
    )
    observed_path = (
        cosplice_output / "04_network" / "observed" / "event_module_assignments.tsv"
    )
    required = [matrix_path, samples_path, consensus_path, observed_path]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing differential-network benchmark input: {missing}")
    matrix = pd.read_csv(matrix_path, sep="\t", index_col=0)
    samples = pd.read_csv(samples_path, sep="\t", low_memory=False)
    if "condition" not in samples:
        raise ValueError("Setd2 differential-network benchmark requires condition metadata.")
    consensus = pd.read_csv(consensus_path, sep="\t", low_memory=False)
    observed = pd.read_csv(observed_path, sep="\t", low_memory=False)
    connectivity = observed.drop_duplicates("event_key").set_index("event_key")[
        "network_connectivity"
    ]
    consensus["network_connectivity"] = consensus["event_key"].map(connectivity)
    selected = (
        consensus.loc[consensus["consensus_module"].ne("grey")]
        .sort_values(
            ["consensus_module", "network_connectivity"],
            ascending=[True, False],
        )
        .groupby("consensus_module", observed=True)
        .head(events_per_module)
    )
    matrix = matrix.reindex(selected["event_key"]).dropna()
    selected = selected.loc[selected["event_key"].isin(matrix.index)]
    result = compare_condition_networks(
        matrix,
        samples,
        selected[["event_key", "consensus_module"]],
        condition_column="condition",
        min_samples_per_condition=3,
        soft_power=6,
        permutations=permutations,
        max_output_edges_per_module_pair=250,
        random_seed=117,
    )
    changes = result.module_changes
    metrics = pd.DataFrame(
        [
            {
                "dataset": cosplice_output.parent.name,
                "algorithm": "differential_network",
                "metric": "module_pair_tests",
                "mean": float(len(changes)),
                "median": np.nan,
                "standard_deviation": np.nan,
                "n_iterations": permutations,
                "n_events": len(matrix),
            },
            {
                "dataset": cosplice_output.parent.name,
                "algorithm": "differential_network",
                "metric": "significant_density_or_topology_modules",
                "mean": float(
                    changes[
                        [
                            "has_significant_density_change",
                            "has_significant_topology_rewiring",
                        ]
                    ]
                    .any(axis=1)
                    .sum()
                ),
                "median": np.nan,
                "standard_deviation": np.nan,
                "n_iterations": permutations,
                "n_events": len(matrix),
            },
        ]
    )
    return metrics, changes, result.event_rewiring, result.condition_edges


def benchmark_real_multiview(
    cosplice_output: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    edges_path = cosplice_output / "04_network" / "observed" / "network_edges.tsv"
    motif_path = (
        cosplice_output
        / "06e_module_interpretation"
        / "event_rbp_motif_hits.tsv"
    )
    isoform_path = (
        cosplice_output
        / "06e_module_interpretation"
        / "dominant_isoform_switches.tsv"
    )
    if not edges_path.exists():
        raise FileNotFoundError(f"Missing AS edge table: {edges_path}")
    edges = pd.read_csv(edges_path, sep="\t", low_memory=False)
    event_keys = sorted(
        set(edges["event_a"].astype(str)).union(edges["event_b"].astype(str))
    )
    features = pd.DataFrame({"event_key": event_keys})
    features["gene_id"] = features["event_key"].str.split("|", regex=False).str[0]
    auxiliary_columns: list[str] = []
    if motif_path.exists():
        motif = pd.read_csv(motif_path, sep="\t", low_memory=False)
        motif_summary = motif.groupby("event_key", observed=True).agg(
            motif_site_count=("site_count", "sum"),
            motif_max_relative_score=("max_relative_score", "max"),
        )
        features = features.join(motif_summary, on="event_key")
        auxiliary_columns.extend(["motif_site_count", "motif_max_relative_score"])
    if isoform_path.exists():
        isoform = pd.read_csv(isoform_path, sep="\t", low_memory=False)
        isoform_summary = isoform.groupby("gene_id", observed=True).agg(
            isoform_max_abs_delta_IF=("max_abs_delta_IF", "max"),
            DTU_supported=("DTU_significant", "max"),
        )
        isoform_summary["DTU_supported"] = (
            isoform_summary["DTU_supported"].astype(bool).astype(float)
        )
        features = features.join(isoform_summary, on="gene_id")
        auxiliary_columns.extend(["isoform_max_abs_delta_IF", "DTU_supported"])
    if not auxiliary_columns:
        raise ValueError("No real mapped event-level auxiliary evidence was found.")
    result = fuse_event_edge_views(
        edges,
        features,
        auxiliary_columns=auxiliary_columns,
        as_weight=0.80,
        edge_threshold=0.10,
        min_module_size=20,
        resolution=1.0,
    )
    fused = result.fused_edges
    metrics = pd.DataFrame(
        [
            {
                "dataset": cosplice_output.parent.name,
                "algorithm": "AS_primary_multiview",
                "metric": "evaluated_AS_edges",
                "mean": float(len(fused)),
                "median": np.nan,
                "standard_deviation": np.nan,
                "n_iterations": 1,
                "n_events": len(event_keys),
            },
            {
                "dataset": cosplice_output.parent.name,
                "algorithm": "AS_primary_multiview",
                "metric": "edges_with_auxiliary_evidence_fraction",
                "mean": float(fused["auxiliary_view_coverage"].gt(0).mean()),
                "median": np.nan,
                "standard_deviation": np.nan,
                "n_iterations": 1,
                "n_events": len(event_keys),
            },
            {
                "dataset": cosplice_output.parent.name,
                "algorithm": "AS_primary_multiview",
                "metric": "median_absolute_edge_reweighting",
                "mean": float(
                    (fused["fused_adjacency"] - fused["adjacency"]).abs().median()
                ),
                "median": np.nan,
                "standard_deviation": np.nan,
                "n_iterations": 1,
                "n_events": len(event_keys),
            },
        ]
    )
    return metrics, fused, result.assignments, result.view_audit


def _recommendations(
    synthetic: pd.DataFrame,
    real: pd.DataFrame,
    public: pd.DataFrame,
) -> pd.DataFrame:
    def metric_value(
        algorithm: str, configuration: str, metric: str
    ) -> float:
        match = synthetic.loc[
            synthetic["algorithm"].eq(algorithm)
            & synthetic["configuration"].eq(configuration)
            & synthetic["metric"].eq(metric),
            "value",
        ]
        return float(match.iloc[0]) if not match.empty else np.nan

    uncertainty_precision = synthetic.loc[
        synthetic["algorithm"].eq("uncertainty_weighting")
        & synthetic["metric"].eq("recoverable_top_edge_precision")
    ].set_index("configuration")["value"]
    raw_precision = float(uncertainty_precision["raw_PSI"])
    real_index = (
        real.set_index("algorithm")["median"].to_dict() if not real.empty else {}
    )
    raw_stability = float(real_index.get("raw_PSI", np.nan))
    best_mode = max(real_index, key=real_index.get) if real_index else "not_evaluated"
    best_stability = float(real_index.get(best_mode, np.nan))
    best_precision = float(uncertainty_precision.get(best_mode, np.nan))
    informative_auc = metric_value(
        "multiview_AS_primary", "informative_auxiliary_20pct", "true_edge_AUROC"
    )
    base_auc = metric_value(
        "multiview_AS_primary", "AS_only", "true_edge_AUROC"
    )
    adversarial_auc = metric_value(
        "multiview_AS_primary", "adversarial_auxiliary_20pct", "true_edge_AUROC"
    )
    gained = metric_value(
        "differential_network",
        "permutation_density_and_topology",
        "known_gained_module_detected",
    )
    false_positive = metric_value(
        "differential_network",
        "permutation_density_and_topology",
        "stable_module_false_positive",
    )
    preserved = metric_value(
        "bayesian_preservation",
        "hierarchical_multi_cohort",
        "preserved_module_probability",
    )
    null_probability = metric_value(
        "bayesian_preservation",
        "hierarchical_multi_cohort",
        "null_module_probability",
    )
    return pd.DataFrame(
        [
            {
                "algorithm": "uncertainty_weighted_network",
                "recommended": bool(
                    best_mode != "raw_PSI"
                    and best_precision >= raw_precision - 0.02
                    and best_stability > raw_stability
                ),
                "evidence": (
                    f"best mode={best_mode}; synthetic top-edge precision "
                    f"{raw_precision:.3f} -> {best_precision:.3f}; real bootstrap "
                    f"median Jaccard {raw_stability:.3f} -> {best_stability:.3f}"
                ),
                "scope": "requires replicate-level IJC/SJC counts",
            },
            {
                "algorithm": "differential_network",
                "recommended": bool(gained >= 1.0 and false_positive <= 0.0),
                "evidence": (
                    f"known gained-module detection={gained:.0f}; stable-module "
                    f"false positive={false_positive:.0f}"
                ),
                "scope": "condition topology follow-up; not a mean dPSI test",
            },
            {
                "algorithm": "AS_primary_multiview",
                "recommended": bool(
                    informative_auc > base_auc and adversarial_auc >= base_auc - 0.05
                ),
                "evidence": (
                    f"edge AUROC AS={base_auc:.3f}, informative={informative_auc:.3f}, "
                    f"adversarial={adversarial_auc:.3f}"
                ),
                "scope": (
                    "optional when mapped event-level auxiliary evidence exists; "
                    "auxiliary evidence cannot create AS edges"
                ),
            },
            {
                "algorithm": "bayesian_cross_cohort_preservation",
                "recommended": bool(preserved >= 0.95 and null_probability <= 0.20),
                "evidence": (
                    f"preserved posterior={preserved:.3f}; null posterior="
                    f"{null_probability:.3f}; public methods evaluated={len(public)}"
                ),
                "scope": "requires at least two external cohort projections",
            },
        ]
    )


def render_benchmark_plot(
    result: MethodExtensionBenchmarkResult, output_path: Path
) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    synthetic = result.synthetic_metrics
    panels = [
        (
            "uncertainty_weighting",
            "recoverable_top_edge_precision",
            "Known-truth edge precision",
        ),
        (
            "multiview_AS_primary",
            "true_edge_AUROC",
            "AS-primary fusion guardrail",
        ),
    ]
    for axis, (algorithm, metric, title) in zip(axes.flat[:2], panels, strict=True):
        part = synthetic.loc[
            synthetic["algorithm"].eq(algorithm)
            & synthetic["metric"].eq(metric)
        ]
        axis.bar(
            np.arange(len(part)),
            part["value"],
            color=["#6b7280", "#2563eb", "#7c3aed", "#0f766e"][: len(part)],
        )
        axis.set_xticks(np.arange(len(part)))
        axis.set_xticklabels(part["configuration"], rotation=25, ha="right")
        axis.set_ylim(0, 1.05)
        axis.set_title(title, fontweight="bold")
        axis.set_ylabel(metric)
    real = result.real_count_metrics.loc[
        result.real_count_metrics["metric"].eq("bootstrap_top_edge_Jaccard")
    ]
    axes[1, 0].bar(
        np.arange(len(real)),
        real["median"],
        color=["#6b7280", "#2563eb", "#7c3aed", "#0f766e"][: len(real)],
    )
    axes[1, 0].set_xticks(np.arange(len(real)))
    axes[1, 0].set_xticklabels(real["algorithm"], rotation=20, ha="right")
    axes[1, 0].set_ylim(0, max(1.0, float(real["median"].max()) * 1.15))
    axes[1, 0].set_title("Setd2 count-aware edge stability", fontweight="bold")
    axes[1, 0].set_ylabel("Median bootstrap top-edge Jaccard")
    public = result.public_preservation_metrics
    positions = np.arange(len(public))
    width = 0.35
    axes[1, 1].bar(
        positions - width / 2,
        public["fixed_partial_or_better_fraction"],
        width,
        label="Fixed threshold",
        color="#6b7280",
    )
    axes[1, 1].bar(
        positions + width / 2,
        public["bayesian_partial_or_better_fraction"],
        width,
        label="Bayesian posterior",
        color="#0f766e",
    )
    axes[1, 1].set_xticks(positions)
    axes[1, 1].set_xticklabels(public["method"], rotation=20, ha="right")
    axes[1, 1].set_ylim(0, 1.05)
    axes[1, 1].set_title("Three-cohort preservation calls", fontweight="bold")
    axes[1, 1].set_ylabel("Fraction partial-or-better")
    axes[1, 1].legend(frameon=False)
    fig.suptitle(
        "CoSplice method-extension ablation and external validation",
        fontsize=16,
        fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output_path


def run_method_extension_benchmark(
    *,
    cosplice_output: Path,
    public_preservation_path: Path,
    output_dir: Path,
) -> MethodExtensionBenchmarkResult:
    output_dir.mkdir(parents=True, exist_ok=True)
    synthetic = run_synthetic_method_ablation()
    real, uncertainty_audit = benchmark_real_read_uncertainty(cosplice_output)
    differential_metrics, differential_modules, differential_events, differential_edges = (
        benchmark_real_differential_network(cosplice_output)
    )
    multiview_metrics, fused_edges, multiview_assignments, view_audit = (
        benchmark_real_multiview(cosplice_output)
    )
    real = pd.concat(
        [real, differential_metrics, multiview_metrics],
        ignore_index=True,
        sort=False,
    )
    public, public_posteriors = benchmark_public_bayesian_preservation(
        public_preservation_path
    )
    recommendation = _recommendations(synthetic, real, public)
    result = MethodExtensionBenchmarkResult(
        synthetic,
        real,
        public,
        recommendation,
    )
    synthetic.to_csv(output_dir / "synthetic_known_truth_ablation.tsv", sep="\t", index=False)
    real.to_csv(output_dir / "setd2_read_count_ablation.tsv", sep="\t", index=False)
    uncertainty_audit.to_csv(
        output_dir / "setd2_uncertainty_event_audit.tsv", sep="\t", index=False
    )
    differential_modules.to_csv(
        output_dir / "setd2_differential_module_network.tsv",
        sep="\t",
        index=False,
    )
    differential_events.to_csv(
        output_dir / "setd2_differential_event_rewiring.tsv",
        sep="\t",
        index=False,
    )
    differential_edges.to_csv(
        output_dir / "setd2_differential_top_edges.tsv",
        sep="\t",
        index=False,
    )
    fused_edges.to_csv(
        output_dir / "setd2_multiview_fused_edges.tsv",
        sep="\t",
        index=False,
    )
    multiview_assignments.to_csv(
        output_dir / "setd2_multiview_module_assignments.tsv",
        sep="\t",
        index=False,
    )
    view_audit.to_csv(
        output_dir / "setd2_multiview_view_audit.tsv",
        sep="\t",
        index=False,
    )
    public.to_csv(
        output_dir / "public_three_cohort_preservation_head_to_head.tsv",
        sep="\t",
        index=False,
    )
    public_posteriors.to_csv(
        output_dir / "public_bayesian_module_posteriors.tsv", sep="\t", index=False
    )
    recommendation.to_csv(
        output_dir / "method_recommendation.tsv", sep="\t", index=False
    )
    render_benchmark_plot(result, output_dir / "method_extension_benchmark.png")
    return result
