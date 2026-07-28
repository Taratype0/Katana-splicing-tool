from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Callable

import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from katana_cosplice.benchmark import (
    _fit_method,
    _program_score,
    _rscript_with_wgcna,
)
from katana_cosplice.config import (
    BenchmarkConfig,
    HybridProgramConfig,
    NetworkConfig,
)
from katana_cosplice.network import _module_labels, correlation_matrix


Progress = Callable[[str], None] | None


@dataclass
class HybridProgramResult:
    event_assignments: pd.DataFrame
    program_summary: pd.DataFrame
    program_scores: pd.DataFrame
    bootstrap_runs: pd.DataFrame
    engine_audit: pd.DataFrame
    overlap: pd.DataFrame


def _python_soft_tom_labels(
    matrix: pd.DataFrame,
    network: NetworkConfig,
    settings: HybridProgramConfig,
) -> np.ndarray:
    """Signed soft adjacency + TOM fallback when exact R WGCNA is unavailable."""
    if len(matrix) > settings.python_max_events:
        raise ValueError(
            f"Python soft TOM requires {len(matrix)} events, exceeding "
            f"hybrid_programs.python_max_events={settings.python_max_events}."
        )
    correlation = correlation_matrix(matrix, network.method).astype(np.float32, copy=False)
    base = (1.0 + correlation) / 2.0 if network.signed else np.abs(correlation)
    adjacency = np.power(np.clip(base, 0.0, 1.0), network.soft_power).astype(
        np.float32,
        copy=False,
    )
    np.fill_diagonal(adjacency, 0.0)
    connectivity = adjacency.sum(axis=1, dtype=np.float64)
    shared = adjacency @ adjacency
    denominator = np.minimum.outer(connectivity, connectivity) + 1.0 - adjacency
    tom = np.divide(
        shared + adjacency,
        denominator,
        out=np.zeros_like(shared),
        where=denominator > 1e-12,
    )
    tom = np.clip((tom + tom.T) / 2.0, 0.0, 1.0)
    np.fill_diagonal(tom, 1.0)
    distance = np.clip(1.0 - tom, 0.0, 1.0)
    hierarchy = linkage(squareform(distance, checks=False), method="average")
    raw = fcluster(
        hierarchy,
        t=settings.python_tom_cut_height,
        criterion="distance",
    )
    min_size = settings.min_module_size or network.min_module_size
    labels, _ = _module_labels(raw, min_size)
    return labels


def _diffuse_labeler(
    matrix: pd.DataFrame,
    network: NetworkConfig,
    benchmark: BenchmarkConfig,
    settings: HybridProgramConfig,
) -> tuple[np.ndarray, str, str]:
    requested = settings.diffuse_engine
    rscript, r_audit = _rscript_with_wgcna()
    if requested in {"auto", "wgcna_exact"} and rscript:
        labels, _ = _fit_method(
            "wgcna_exact",
            matrix,
            replace(
                network,
                min_module_size=settings.min_module_size or network.min_module_size,
            ),
            replace(benchmark, max_events=len(matrix)),
            rscript,
        )
        return labels, "wgcna_exact", r_audit
    if requested == "wgcna_exact":
        raise RuntimeError(r_audit)
    labels = _python_soft_tom_labels(matrix, network, settings)
    fallback = (
        "Python signed soft-adjacency/TOM fallback; not exact WGCNA blockwiseModules."
        if requested == "python_soft_tom"
        else f"Exact WGCNA unavailable ({r_audit}); used Python signed soft TOM."
    )
    return labels, "python_soft_tom", fallback


def _bootstrap_indices(
    samples: pd.DataFrame,
    columns: pd.Index,
    rng: np.random.Generator,
) -> np.ndarray:
    metadata = samples.set_index("sample_id").reindex(columns)
    if metadata.index.has_duplicates or metadata["condition"].isna().any():
        raise ValueError("Hybrid bootstrap requires unique sample_id and condition metadata.")
    indices: list[int] = []
    positions = pd.Series(np.arange(len(columns)), index=columns)
    for _condition, frame in metadata.groupby("condition", sort=False):
        group_positions = positions.reindex(frame.index).to_numpy(dtype=int)
        indices.extend(rng.choice(group_positions, size=len(group_positions), replace=True))
    return np.asarray(indices, dtype=int)


def _program_scores(matrix: pd.DataFrame, labels: np.ndarray, prefix: str) -> pd.DataFrame:
    output: dict[str, np.ndarray] = {}
    for module in sorted(set(labels).difference({"grey"})):
        indices = np.flatnonzero(labels == module)
        if len(indices) >= 2:
            output[f"{prefix}{module}"] = _program_score(
                matrix.iloc[indices].to_numpy(dtype=float)
            )
    scores = pd.DataFrame(output, index=matrix.columns)
    scores.index.name = "sample_id"
    return scores


def _diffuse_bootstrap_confidence(
    matrix: pd.DataFrame,
    baseline_labels: np.ndarray,
    samples: pd.DataFrame,
    network: NetworkConfig,
    benchmark: BenchmarkConfig,
    settings: HybridProgramConfig,
    engine: str,
    progress: Progress,
) -> tuple[np.ndarray, pd.DataFrame]:
    n_events = len(matrix)
    event_totals = np.zeros(n_events, dtype=float)
    rows: list[dict[str, object]] = []
    rng = np.random.default_rng(settings.random_seed)
    run_settings = replace(settings, diffuse_engine=engine)
    for iteration in range(1, settings.bootstrap_iterations + 1):
        sample_indices = _bootstrap_indices(samples, matrix.columns, rng)
        bootstrap = matrix.iloc[:, sample_indices].copy()
        bootstrap.columns = [
            f"diffuse_{iteration:03d}_{position:03d}"
            for position in range(len(sample_indices))
        ]
        labels, _engine, _audit = _diffuse_labeler(
            bootstrap,
            network,
            benchmark,
            run_settings,
        )
        for module in sorted(set(baseline_labels).difference({"grey"})):
            indices = np.flatnonzero(baseline_labels == module)
            observed = labels[indices]
            counts = pd.Series(observed[observed != "grey"]).value_counts()
            denominator = len(indices) * max(1, len(indices) - 1)
            paired = float(sum(count * (count - 1) for count in counts.to_numpy()))
            module_coassignment = paired / denominator if denominator else 0.0
            rows.append(
                {
                    "bootstrap_iteration": iteration,
                    "diffuse_program": f"D-{module}",
                    "n_events": len(indices),
                    "within_program_coassignment": module_coassignment,
                }
            )
            count_map = counts.to_dict()
            for index, label in zip(indices, observed, strict=True):
                if label != "grey" and len(indices) > 1:
                    event_totals[index] += (count_map.get(label, 0) - 1) / (
                        len(indices) - 1
                    )
        if progress:
            progress(
                f"Diffuse TOM bootstrap {iteration}/{settings.bootstrap_iterations}"
            )
    return event_totals / settings.bootstrap_iterations, pd.DataFrame(rows)


def _event_membership(
    matrix: pd.DataFrame,
    labels: np.ndarray,
    scores: pd.DataFrame,
    prefix: str,
) -> np.ndarray:
    values = np.full(len(matrix), np.nan, dtype=float)
    for module in set(labels).difference({"grey"}):
        program = f"{prefix}{module}"
        score = scores[program]
        for index in np.flatnonzero(labels == module):
            event = pd.to_numeric(matrix.iloc[index], errors="coerce")
            valid = event.notna() & score.notna()
            if valid.sum() >= 3 and event.loc[valid].std() > 1e-12:
                values[index] = abs(
                    float(np.corrcoef(event.loc[valid], score.loc[valid])[0, 1])
                )
    return values


def _confidence_class(value: float, settings: HybridProgramConfig) -> str:
    if value >= settings.robust_confidence:
        return "robust"
    if value >= settings.moderate_confidence:
        return "moderate"
    return "exploratory"


def _program_overlap(
    tight: pd.Series,
    diffuse: pd.Series,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    tight_sets = {
        f"T-{module}": set(tight.index[tight.eq(module)])
        for module in set(tight).difference({"grey"})
    }
    diffuse_sets = {
        f"D-{module}": set(diffuse.index[diffuse.eq(module)])
        for module in set(diffuse).difference({"grey"})
    }
    for diffuse_name, diffuse_events in diffuse_sets.items():
        for tight_name, tight_events in tight_sets.items():
            intersection = len(diffuse_events & tight_events)
            union = len(diffuse_events | tight_events)
            rows.append(
                {
                    "diffuse_program": diffuse_name,
                    "tight_program": tight_name,
                    "shared_events": intersection,
                    "jaccard": intersection / union if union else 0.0,
                    "diffuse_fraction_shared": intersection / len(diffuse_events),
                    "tight_fraction_shared": intersection / len(tight_events),
                }
            )
    return pd.DataFrame(rows)


def build_hybrid_programs(
    matrix: pd.DataFrame,
    qc_metrics: pd.DataFrame,
    samples: pd.DataFrame,
    tight_assignments: pd.DataFrame,
    network: NetworkConfig,
    benchmark: BenchmarkConfig,
    settings: HybridProgramConfig,
    tight_coassignment: pd.Series | None = None,
    progress: Progress = None,
) -> HybridProgramResult:
    if matrix.empty:
        raise ValueError("Hybrid programs require a non-empty stable PSI matrix.")
    tight = (
        tight_assignments.set_index("event_key")
        .reindex(matrix.index)
        .copy()
    )
    tight_labels = tight["module"].fillna("grey").astype(str)
    diffuse_labels, engine, audit = _diffuse_labeler(
        matrix,
        network,
        benchmark,
        settings,
    )
    diffuse_scores = _program_scores(matrix, diffuse_labels, "D-")
    diffuse_bootstrap, bootstrap_runs = _diffuse_bootstrap_confidence(
        matrix,
        diffuse_labels,
        samples,
        network,
        benchmark,
        settings,
        engine,
        progress,
    )
    diffuse_membership = _event_membership(
        matrix,
        diffuse_labels,
        diffuse_scores,
        "D-",
    )
    qc = qc_metrics.reindex(matrix.index)
    if "replicate_repeatability" not in qc:
        qc["replicate_repeatability"] = 0.0
    repeatability = (
        pd.to_numeric(qc["replicate_repeatability"], errors="coerce")
        .clip(0.0, 1.0)
        .fillna(0.0)
        .to_numpy()
    )
    tight_membership_source = (
        tight["module_membership"]
        if "module_membership" in tight
        else pd.Series(0.0, index=tight.index)
    )
    tight_membership = (
        pd.to_numeric(tight_membership_source, errors="coerce")
        .abs()
        .clip(0.0, 1.0)
        .fillna(0.0)
        .to_numpy()
    )
    if tight_coassignment is not None:
        tight_bootstrap_source = tight_coassignment.reindex(tight.index)
        tight_bootstrap_source_name = "CoSplice bootstrap co-assignment"
    elif "mean_incident_bootstrap_support" in tight:
        tight_bootstrap_source = tight["mean_incident_bootstrap_support"]
        tight_bootstrap_source_name = "network edge bootstrap support fallback"
    else:
        tight_bootstrap_source = pd.Series(0.0, index=tight.index)
        tight_bootstrap_source_name = "unavailable; assigned zero"
    tight_bootstrap = (
        pd.to_numeric(tight_bootstrap_source, errors="coerce")
        .clip(0.0, 1.0)
        .fillna(0.0)
        .to_numpy()
    )
    tight_confidence = (
        0.45 * tight_bootstrap + 0.35 * tight_membership + 0.20 * repeatability
    )
    diffuse_confidence = (
        0.45 * diffuse_bootstrap
        + 0.35 * np.nan_to_num(diffuse_membership)
        + 0.20 * repeatability
    )
    assignments = pd.DataFrame(
        {
            "event_key": matrix.index,
            "tight_program": [
                f"T-{module}" if module != "grey" else "grey"
                for module in tight_labels
            ],
            "diffuse_program": [
                f"D-{module}" if module != "grey" else "grey"
                for module in diffuse_labels
            ],
            "tight_membership": tight_membership,
            "tight_bootstrap_support": tight_bootstrap,
            "tight_internal_confidence": tight_confidence,
            "diffuse_membership_kME": diffuse_membership,
            "diffuse_bootstrap_coassignment": diffuse_bootstrap,
            "diffuse_internal_confidence": diffuse_confidence,
            "replicate_repeatability": repeatability,
        }
    )
    assignments["program_scope"] = np.select(
        [
            assignments["tight_program"].ne("grey")
            & assignments["diffuse_program"].ne("grey"),
            assignments["tight_program"].ne("grey"),
            assignments["diffuse_program"].ne("grey"),
        ],
        ["both", "tight_only", "diffuse_only"],
        default="unassigned",
    )
    summaries: list[dict[str, object]] = []
    for program_type, column, confidence_column, membership_column, bootstrap_column in (
        (
            "tight",
            "tight_program",
            "tight_internal_confidence",
            "tight_membership",
            "tight_bootstrap_support",
        ),
        (
            "diffuse",
            "diffuse_program",
            "diffuse_internal_confidence",
            "diffuse_membership_kME",
            "diffuse_bootstrap_coassignment",
        ),
    ):
        for program, frame in assignments.loc[
            assignments[column].ne("grey")
        ].groupby(column, sort=True):
            internal = float(frame[confidence_column].median())
            summaries.append(
                {
                    "program": program,
                    "program_type": program_type,
                    "n_events": len(frame),
                    "median_membership": float(frame[membership_column].median()),
                    "median_bootstrap_coassignment": float(
                        frame[bootstrap_column].median()
                    ),
                    "median_replicate_repeatability": float(
                        frame["replicate_repeatability"].median()
                    ),
                    "internal_confidence": internal,
                    "external_validation_status": "not_configured",
                    "external_preservation_confidence": np.nan,
                    "overall_confidence": internal,
                    "confidence_scope": "internal_only",
                    "confidence_class": _confidence_class(internal, settings),
                }
            )
    summary = pd.DataFrame(summaries)
    tight_scores = _program_scores(matrix, tight_labels.to_numpy(object), "T-")
    scores = pd.concat([tight_scores, diffuse_scores], axis=1)
    overlap = _program_overlap(
        tight_labels,
        pd.Series(diffuse_labels, index=matrix.index),
    )
    if not overlap.empty:
        nearest = (
            overlap.sort_values("jaccard", ascending=False)
            .drop_duplicates("diffuse_program")
            .set_index("diffuse_program")
        )
        summary["nearest_other_program"] = summary["program"].map(
            nearest["tight_program"]
        )
        summary["nearest_program_jaccard"] = summary["program"].map(
            nearest["jaccard"]
        )
    else:
        summary["nearest_other_program"] = ""
        summary["nearest_program_jaccard"] = np.nan
    audit_table = pd.DataFrame(
        [
            {
                "requested_engine": settings.diffuse_engine,
                "resolved_engine": engine,
                "matrix_events": len(matrix),
                "matrix_samples": matrix.shape[1],
                "bootstrap_iterations": settings.bootstrap_iterations,
                "tight_bootstrap_source": tight_bootstrap_source_name,
                "message": audit,
            }
        ]
    )
    return HybridProgramResult(
        assignments,
        summary,
        scores,
        bootstrap_runs,
        audit_table,
        overlap,
    )


def add_external_preservation_confidence(
    summary: pd.DataFrame,
    preservation: pd.DataFrame,
    settings: HybridProgramConfig,
) -> pd.DataFrame:
    output = summary.copy()
    if preservation.empty:
        return output
    class_score = {
        "preserved": 1.0,
        "partially_preserved": 0.5,
        "not_preserved": 0.0,
    }
    values = preservation.copy()
    values["external_score"] = values["preservation_class"].map(class_score)
    external = values.groupby("consensus_module")["external_score"].mean()
    output["external_preservation_confidence"] = output["program"].map(external)
    available = output["external_preservation_confidence"].notna()
    output.loc[available, "external_validation_status"] = "evaluated"
    output.loc[available, "confidence_scope"] = "internal_plus_external"
    output.loc[available, "overall_confidence"] = (
        0.75 * output.loc[available, "internal_confidence"]
        + 0.25 * output.loc[available, "external_preservation_confidence"]
    )
    output["confidence_class"] = output["overall_confidence"].map(
        lambda value: _confidence_class(float(value), settings)
    )
    return output


def render_hybrid_program_plots(
    result: HybridProgramResult,
    output_dir: Path,
) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    generated: list[Path] = []
    summary = result.program_summary
    if not summary.empty:
        fig, ax = plt.subplots(figsize=(11, 6.5))
        palette = {"tight": "#0B7285", "diffuse": "#E67700"}
        sns.scatterplot(
            data=summary,
            x="n_events",
            y="overall_confidence",
            hue="program_type",
            size="median_replicate_repeatability",
            palette=palette,
            sizes=(45, 240),
            ax=ax,
        )
        ax.axhline(0.75, color="#2B8A3E", linestyle="--", linewidth=1)
        ax.axhline(0.55, color="#F08C00", linestyle=":", linewidth=1)
        ax.set_title("Tight and diffuse co-splicing program confidence", fontweight="bold")
        ax.set_xlabel("Stable AS events")
        ax.set_ylabel("Program confidence")
        ax.set_ylim(0, 1.02)
        fig.tight_layout()
        path = output_dir / "tight_diffuse_program_confidence.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        generated.append(path)
    if not result.overlap.empty:
        matrix = result.overlap.pivot(
            index="diffuse_program",
            columns="tight_program",
            values="jaccard",
        )
        fig, ax = plt.subplots(
            figsize=(max(8, matrix.shape[1] * 0.65), max(6, matrix.shape[0] * 0.45))
        )
        sns.heatmap(
            matrix,
            cmap="YlGnBu",
            vmin=0,
            vmax=1,
            ax=ax,
            cbar_kws={"label": "Event-set Jaccard"},
        )
        ax.set_title("Tight versus diffuse program overlap", fontweight="bold")
        fig.tight_layout()
        path = output_dir / "tight_diffuse_program_overlap.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        generated.append(path)
    return generated
