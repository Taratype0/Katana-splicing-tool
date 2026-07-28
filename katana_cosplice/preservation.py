from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from katana_cosplice.config import ExternalCohortConfig, PreservationConfig


Progress = Callable[[str], None] | None


@dataclass
class PreservationResult:
    module_preservation: pd.DataFrame
    event_mapping: pd.DataFrame
    null_summary: pd.DataFrame
    target_audit: pd.DataFrame


def _read_matrix(path: Path) -> pd.DataFrame:
    matrix = pd.read_csv(path, sep=None, engine="python", index_col=0)
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)
    matrix = matrix.apply(pd.to_numeric, errors="coerce")
    if not matrix.index.is_unique or not matrix.columns.is_unique:
        raise ValueError(f"External PSI matrix labels must be unique: {path}")
    return matrix


def _validate_samples(path: Path, columns: pd.Index) -> None:
    samples = pd.read_csv(path, sep=None, engine="python")
    sample_column = "sample_id" if "sample_id" in samples.columns else "sample" if "sample" in samples.columns else ""
    if not sample_column:
        raise ValueError(f"External sample metadata requires sample_id or sample: {path}")
    available = set(samples[sample_column].astype(str))
    missing = [sample for sample in columns.astype(str) if sample not in available]
    if missing:
        raise ValueError(f"External sample metadata does not cover matrix columns: {missing[:5]}")


def _standardize_rows(matrix: pd.DataFrame) -> np.ndarray:
    values = matrix.to_numpy(dtype=np.float32)
    means = np.nanmean(values, axis=1, keepdims=True)
    values = np.where(np.isfinite(values), values, means)
    centered = values - values.mean(axis=1, keepdims=True)
    scale = centered.std(axis=1, ddof=1, keepdims=True)
    return centered / np.where(scale < 1e-8, 1.0, scale)


def _correlation_matrix(standardized: np.ndarray) -> np.ndarray:
    denominator = max(1, standardized.shape[1] - 1)
    correlation = (standardized @ standardized.T) / denominator
    correlation = np.clip(correlation, -1.0, 1.0).astype(np.float32, copy=False)
    np.fill_diagonal(correlation, 1.0)
    return correlation


def _network_statistics(
    reference_corr: np.ndarray,
    target_corr: np.ndarray,
    indices: np.ndarray,
) -> dict[str, float]:
    reference = np.abs(reference_corr[np.ix_(indices, indices)])
    target = np.abs(target_corr[np.ix_(indices, indices)])
    upper = np.triu_indices(len(indices), k=1)
    reference_density = float(reference[upper].mean()) if len(upper[0]) else 0.0
    target_density = float(target[upper].mean()) if len(upper[0]) else 0.0
    reference_connectivity = reference.sum(axis=1) - 1.0
    target_connectivity = target.sum(axis=1) - 1.0
    connectivity_correlation = (
        float(np.corrcoef(reference_connectivity, target_connectivity)[0, 1])
        if len(indices) > 2
        and np.std(reference_connectivity) > 1e-8
        and np.std(target_connectivity) > 1e-8
        else 0.0
    )
    return {
        "reference_density": reference_density,
        "target_density": target_density,
        "density_ratio": target_density / reference_density if reference_density > 0 else 0.0,
        "connectivity_correlation": connectivity_correlation,
    }


def _program_statistics(
    reference_values: np.ndarray,
    target_values: np.ndarray,
    indices: np.ndarray,
) -> dict[str, float]:
    reference = reference_values[indices]
    target = target_values[indices]
    _ref_u, ref_s, ref_vt = np.linalg.svd(reference, full_matrices=False)
    _target_u, target_s, target_vt = np.linalg.svd(target, full_matrices=False)
    ref_score = ref_vt[0]
    target_score = target_vt[0]
    ref_membership = np.asarray([np.corrcoef(row, ref_score)[0, 1] for row in reference])
    target_membership = np.asarray([np.corrcoef(row, target_score)[0, 1] for row in target])
    membership_correlation = (
        float(np.corrcoef(ref_membership, target_membership)[0, 1])
        if np.std(ref_membership) > 1e-8 and np.std(target_membership) > 1e-8
        else 0.0
    )
    membership_correlation = abs(membership_correlation)
    ref_total = float(np.square(ref_s).sum())
    target_total = float(np.square(target_s).sum())
    ref_variance = float(ref_s[0] ** 2 / ref_total) if ref_total else 0.0
    target_variance = float(target_s[0] ** 2 / target_total) if target_total else 0.0
    variance_similarity = (
        min(ref_variance / target_variance, target_variance / ref_variance)
        if ref_variance > 0 and target_variance > 0
        else 0.0
    )
    return {
        "membership_profile_correlation": membership_correlation,
        "reference_program_variance_explained": ref_variance,
        "target_program_variance_explained": target_variance,
        "program_variance_similarity": variance_similarity,
        "program_score_preservation": membership_correlation * variance_similarity,
    }


def _zscore(value: float, null: np.ndarray) -> float:
    scale = float(np.std(null, ddof=1))
    return (value - float(np.mean(null))) / scale if scale > 1e-12 else 0.0


def _empirical_p(value: float, null: np.ndarray) -> float:
    return float((1 + np.count_nonzero(null >= value)) / (len(null) + 1))


def _category(zsummary: float, settings: PreservationConfig) -> str:
    if zsummary >= settings.preserved_zsummary:
        return "preserved"
    if zsummary >= settings.partial_zsummary:
        return "partially_preserved"
    return "not_preserved"


def evaluate_external_module_preservation(
    reference_matrix: pd.DataFrame,
    consensus_assignments: pd.DataFrame,
    cohorts: list[ExternalCohortConfig],
    settings: PreservationConfig,
    progress: Progress = None,
) -> PreservationResult:
    reference_matrix = reference_matrix.copy()
    reference_matrix.index = reference_matrix.index.astype(str)
    module_series = consensus_assignments.set_index("event_key")["consensus_module"].astype(str)
    module_series.index = module_series.index.astype(str)
    module_series = module_series.loc[module_series.ne("grey")]
    module_rows: list[dict[str, object]] = []
    mapping_rows: list[dict[str, object]] = []
    null_rows: list[dict[str, object]] = []
    audit_rows: list[dict[str, object]] = []
    rng = np.random.default_rng(settings.random_seed)
    for cohort_index, cohort in enumerate(cohorts, start=1):
        if not cohort.psi_matrix.exists() or not cohort.sample_metadata.exists():
            audit_rows.append(
                {
                    "cohort": cohort.name,
                    "status": "missing_input",
                    "psi_matrix": str(cohort.psi_matrix),
                    "sample_metadata": str(cohort.sample_metadata),
                    "message": "PSI matrix or sample metadata is missing.",
                }
            )
            continue
        target = _read_matrix(cohort.psi_matrix)
        _validate_samples(cohort.sample_metadata, target.columns)
        target_quality = target.notna().sum(axis=1).ge(3) & target.std(axis=1, skipna=True).gt(1e-8)
        target_events = set(target.index[target_quality])
        module_events = set(module_series.index)
        shared = pd.Index(
            [event for event in reference_matrix.index if event in target_events and event in module_events],
            dtype=object,
        )
        if len(shared) < settings.min_shared_events:
            audit_rows.append(
                {
                    "cohort": cohort.name,
                    "status": "insufficient_shared_events",
                    "psi_matrix": str(cohort.psi_matrix),
                    "sample_metadata": str(cohort.sample_metadata),
                    "message": f"Only {len(shared)} shared non-grey events.",
                }
            )
            continue
        if len(shared) > settings.dense_max_events:
            raise ValueError(
                f"External cohort {cohort.name} has {len(shared)} shared events, exceeding "
                f"preservation.dense_max_events={settings.dense_max_events}."
            )
        reference_values = _standardize_rows(reference_matrix.reindex(shared))
        target_values = _standardize_rows(target.reindex(shared))
        reference_corr = _correlation_matrix(reference_values)
        target_corr = _correlation_matrix(target_values)
        positions = pd.Series(np.arange(len(shared)), index=shared)
        shared_modules = module_series.reindex(shared)
        for event_key, module in shared_modules.items():
            mapping_rows.append(
                {
                    "cohort": cohort.name,
                    "event_key": event_key,
                    "consensus_module": module,
                    "mapped_exact_event_key": True,
                }
            )
        for module, event_keys in shared_modules.groupby(shared_modules, observed=True):
            indices = positions.reindex(event_keys.index).to_numpy(dtype=int)
            if len(indices) < settings.min_shared_events:
                continue
            network_stats = _network_statistics(reference_corr, target_corr, indices)
            program_stats = _program_statistics(reference_values, target_values, indices)
            density_null = np.zeros(settings.permutations, dtype=float)
            connectivity_null = np.zeros(settings.permutations, dtype=float)
            for permutation in range(settings.permutations):
                random_indices = rng.choice(len(shared), size=len(indices), replace=False)
                random_stats = _network_statistics(reference_corr, target_corr, random_indices)
                density_null[permutation] = random_stats["target_density"]
                connectivity_null[permutation] = random_stats["connectivity_correlation"]
            z_density = _zscore(network_stats["target_density"], density_null)
            z_connectivity = _zscore(network_stats["connectivity_correlation"], connectivity_null)
            zsummary = (z_density + z_connectivity) / 2.0
            module_rows.append(
                {
                    "cohort": cohort.name,
                    "consensus_module": str(module),
                    "reference_module_events": int(module_series.eq(module).sum()),
                    "shared_module_events": int(len(indices)),
                    "shared_module_fraction": float(len(indices) / max(1, module_series.eq(module).sum())),
                    **network_stats,
                    **program_stats,
                    "Zdensity": z_density,
                    "Zconnectivity": z_connectivity,
                    "Zsummary_wgcna_like": zsummary,
                    "density_empirical_pvalue": _empirical_p(network_stats["target_density"], density_null),
                    "connectivity_empirical_pvalue": _empirical_p(
                        network_stats["connectivity_correlation"], connectivity_null
                    ),
                    "preservation_class": _category(zsummary, settings),
                }
            )
            null_rows.append(
                {
                    "cohort": cohort.name,
                    "consensus_module": str(module),
                    "permutations": settings.permutations,
                    "density_null_mean": float(density_null.mean()),
                    "density_null_sd": float(density_null.std(ddof=1)),
                    "connectivity_null_mean": float(connectivity_null.mean()),
                    "connectivity_null_sd": float(connectivity_null.std(ddof=1)),
                }
            )
        audit_rows.append(
            {
                "cohort": cohort.name,
                "status": "finished",
                "psi_matrix": str(cohort.psi_matrix),
                "sample_metadata": str(cohort.sample_metadata),
                    "message": (
                        f"Mapped {len(shared)} exact event keys with >=3 finite target PSI values "
                        "and non-zero target variance."
                    ),
            }
        )
        if progress:
            progress(f"External preservation {cohort_index}/{len(cohorts)}: {cohort.name}")
    return PreservationResult(
        pd.DataFrame(module_rows, columns=[
            "cohort", "consensus_module", "reference_module_events", "shared_module_events",
            "shared_module_fraction", "reference_density", "target_density", "density_ratio",
            "connectivity_correlation", "membership_profile_correlation",
            "reference_program_variance_explained", "target_program_variance_explained",
            "program_variance_similarity", "program_score_preservation", "Zdensity",
            "Zconnectivity", "Zsummary_wgcna_like", "density_empirical_pvalue",
            "connectivity_empirical_pvalue", "preservation_class",
        ]),
        pd.DataFrame(mapping_rows, columns=[
            "cohort", "event_key", "consensus_module", "mapped_exact_event_key",
        ]),
        pd.DataFrame(null_rows, columns=[
            "cohort", "consensus_module", "permutations", "density_null_mean",
            "density_null_sd", "connectivity_null_mean", "connectivity_null_sd",
        ]),
        pd.DataFrame(audit_rows, columns=[
            "cohort", "status", "psi_matrix", "sample_metadata", "message",
        ]),
    )


def render_preservation_plots(result: PreservationResult, output_dir: Path) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    if result.module_preservation.empty:
        return []
    generated: list[Path] = []
    pivot = result.module_preservation.pivot(
        index="consensus_module", columns="cohort", values="Zsummary_wgcna_like"
    )
    fig, ax = plt.subplots(figsize=(max(8, pivot.shape[1] * 1.4), max(5.5, pivot.shape[0] * 0.65)))
    sns.heatmap(
        pivot,
        annot=True,
        fmt=".1f",
        cmap="RdYlGn",
        center=2.0,
        ax=ax,
        cbar_kws={"label": "WGCNA-like Zsummary"},
    )
    ax.set_title("External consensus-module preservation", fontweight="bold")
    fig.tight_layout()
    path = output_dir / "external_module_preservation_zsummary.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    generated.append(path)

    fig, ax = plt.subplots(figsize=(9, 6))
    sns.scatterplot(
        data=result.module_preservation,
        x="membership_profile_correlation",
        y="density_ratio",
        hue="preservation_class",
        size="shared_module_events",
        sizes=(50, 250),
        ax=ax,
    )
    ax.set_title("External density and membership preservation", fontweight="bold")
    fig.tight_layout()
    path = output_dir / "external_module_preservation_components.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    generated.append(path)
    return generated
