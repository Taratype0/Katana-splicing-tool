from __future__ import annotations

from dataclasses import dataclass, replace
from html import escape
from pathlib import Path
from typing import Callable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.optimize import linear_sum_assignment
from scipy.stats import hypergeom
from sklearn.metrics import adjusted_rand_score
from statsmodels.stats.multitest import multipletests

from katana_cosplice.config import NetworkConfig, ValidationConfig
from katana_cosplice.network import _program_scores, discover_module_labels


Progress = Callable[[str], None] | None


@dataclass
class ValidationResult:
    event_stability: pd.DataFrame
    module_stability: pd.DataFrame
    bootstrap_runs: pd.DataFrame
    parameter_sensitivity: pd.DataFrame
    preservation_pairs: pd.DataFrame
    preservation_summary: pd.DataFrame
    program_score_confidence: pd.DataFrame


def _module_sets(labels: pd.Series) -> dict[str, set[str]]:
    return {
        str(module): set(group.index.astype(str))
        for module, group in labels.groupby(labels, observed=True)
        if str(module) != "grey"
    }


def _match_modules(reference: pd.Series, candidate: pd.Series) -> tuple[dict[str, str], pd.DataFrame]:
    reference_sets = _module_sets(reference)
    candidate_sets = _module_sets(candidate)
    reference_names = sorted(reference_sets)
    candidate_names = sorted(candidate_sets)
    if not reference_names or not candidate_names:
        return {}, pd.DataFrame(
            columns=[
                "reference_module",
                "candidate_module",
                "intersection",
                "union",
                "jaccard",
                "is_optimal_match",
            ]
        )
    jaccard = np.zeros((len(reference_names), len(candidate_names)), dtype=float)
    intersection = np.zeros_like(jaccard, dtype=int)
    for row, reference_name in enumerate(reference_names):
        for column, candidate_name in enumerate(candidate_names):
            shared = reference_sets[reference_name] & candidate_sets[candidate_name]
            union = reference_sets[reference_name] | candidate_sets[candidate_name]
            intersection[row, column] = len(shared)
            jaccard[row, column] = len(shared) / len(union) if union else 0.0
    row_indices, column_indices = linear_sum_assignment(-jaccard)
    mapping = {candidate_names[column]: reference_names[row] for row, column in zip(row_indices, column_indices, strict=True)}
    rows: list[dict[str, object]] = []
    for row, reference_name in enumerate(reference_names):
        for column, candidate_name in enumerate(candidate_names):
            rows.append(
                {
                    "reference_module": reference_name,
                    "candidate_module": candidate_name,
                    "intersection": int(intersection[row, column]),
                    "union": len(reference_sets[reference_name] | candidate_sets[candidate_name]),
                    "jaccard": float(jaccard[row, column]),
                    "is_optimal_match": mapping.get(candidate_name) == reference_name,
                }
            )
    return mapping, pd.DataFrame(rows)


def _classification(value: float, settings: ValidationConfig) -> str:
    if value >= settings.robust_jaccard:
        return "robust"
    if value >= settings.moderate_jaccard:
        return "moderate"
    return "unstable"


def _bootstrap_indices(samples: pd.DataFrame, columns: pd.Index, rng: np.random.Generator) -> np.ndarray:
    ordered = samples.set_index("sample_id").reindex(columns)
    if "condition" not in ordered:
        return rng.choice(np.arange(len(columns)), size=len(columns), replace=True)
    conditions = ordered["condition"].fillna("unknown").astype(str).to_numpy()
    groups = [np.flatnonzero(conditions == value) for value in pd.unique(conditions)]
    return np.concatenate([rng.choice(group, size=len(group), replace=True) for group in groups])


def _bootstrap_validation(
    matrix: pd.DataFrame,
    baseline: pd.Series,
    samples: pd.DataFrame,
    network: NetworkConfig,
    settings: ValidationConfig,
    progress: Progress,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    baseline = baseline.reindex(matrix.index).fillna("grey").astype(str)
    baseline_sets = _module_sets(baseline)
    membership_hits = pd.Series(0.0, index=matrix.index)
    kme_values: dict[str, list[float]] = {str(key): [] for key in matrix.index}
    hub_hits = pd.Series(0.0, index=matrix.index)
    module_jaccards: dict[str, list[float]] = {module: [] for module in baseline_sets}
    run_rows: list[dict[str, object]] = []
    rng = np.random.default_rng(settings.random_seed)
    for iteration in range(1, settings.bootstrap_iterations + 1):
        indices = _bootstrap_indices(samples, matrix.columns, rng)
        boot_matrix = matrix.iloc[:, indices].copy()
        boot_matrix.columns = [f"boot_{iteration:03d}_{index:03d}" for index in range(len(indices))]
        boot_labels = pd.Series(discover_module_labels(boot_matrix, network), index=matrix.index, dtype=str)
        mapping, matches = _match_modules(baseline, boot_labels)
        mapped = boot_labels.map(mapping).fillna("unmatched")
        membership_hits += mapped.eq(baseline).astype(float)
        fixed_scores = _program_scores(boot_matrix, baseline.to_numpy(dtype=object))
        for module, event_keys in baseline_sets.items():
            if module not in fixed_scores or not event_keys:
                continue
            keys = sorted(event_keys)
            score = fixed_scores[module]
            module_kme = boot_matrix.loc[keys].T.apply(
                lambda values: values.corr(score)
            ).fillna(0.0)
            for event_key, value in module_kme.items():
                kme_values[str(event_key)].append(float(value))
            top_n = max(1, int(np.ceil(len(module_kme) * 0.10)))
            hub_hits.loc[module_kme.abs().nlargest(top_n).index] += 1.0
        optimal = matches.loc[matches["is_optimal_match"]]
        for module in baseline_sets:
            values = optimal.loc[optimal["reference_module"].eq(module), "jaccard"]
            module_jaccards[module].append(float(values.iloc[0]) if not values.empty else 0.0)
        counts = boot_labels.loc[boot_labels.ne("grey")].value_counts()
        run_rows.append(
            {
                "bootstrap_iteration": iteration,
                "adjusted_rand_index": adjusted_rand_score(baseline, boot_labels),
                "n_modules": int(len(counts)),
                "n_grey_events": int(boot_labels.eq("grey").sum()),
                "largest_module_fraction": float(counts.max() / len(boot_labels)) if not counts.empty else 0.0,
            }
        )
        if progress and (iteration == 1 or iteration % 5 == 0 or iteration == settings.bootstrap_iterations):
            progress(f"Module bootstrap {iteration}/{settings.bootstrap_iterations}")
    membership = membership_hits / settings.bootstrap_iterations
    event_stability = pd.DataFrame(
        {
            "event_key": matrix.index,
            "baseline_module": baseline.to_numpy(),
            "membership_stability": membership.to_numpy(),
            "n_bootstraps": settings.bootstrap_iterations,
        }
    )
    event_stability["stability_class"] = event_stability["membership_stability"].map(lambda value: _classification(float(value), settings))
    event_stability["bootstrap_kME_mean"] = event_stability["event_key"].map(
        lambda key: float(np.mean(kme_values[str(key)]))
        if kme_values[str(key)]
        else np.nan
    )
    event_stability["bootstrap_kME_sd"] = event_stability["event_key"].map(
        lambda key: float(np.std(kme_values[str(key)], ddof=1))
        if len(kme_values[str(key)]) > 1
        else np.nan
    )
    event_stability["bootstrap_kME_sign_consistency"] = event_stability["event_key"].map(
        lambda key: (
            float(
                max(
                    np.mean(np.asarray(kme_values[str(key)]) >= 0),
                    np.mean(np.asarray(kme_values[str(key)]) < 0),
                )
            )
            if kme_values[str(key)]
            else np.nan
        )
    )
    event_stability["bootstrap_hub_top_decile_fraction"] = (
        event_stability["event_key"].map(hub_hits).fillna(0.0)
        / settings.bootstrap_iterations
    )
    module_rows: list[dict[str, object]] = []
    for module, event_keys in baseline_sets.items():
        values = np.asarray(module_jaccards[module], dtype=float)
        module_membership = membership.reindex(list(event_keys)).fillna(0.0)
        combined = min(float(values.mean()), float(module_membership.mean()))
        module_rows.append(
            {
                "module": module,
                "n_events": len(event_keys),
                "mean_best_jaccard": float(values.mean()),
                "median_best_jaccard": float(np.median(values)),
                "min_best_jaccard": float(values.min()),
                "mean_event_membership_stability": float(module_membership.mean()),
                "robust_event_fraction": float(module_membership.ge(settings.robust_jaccard).mean()),
                "median_bootstrap_kME_sd": float(
                    event_stability.loc[
                        event_stability["baseline_module"].eq(module),
                        "bootstrap_kME_sd",
                    ].median()
                ),
                "max_hub_top_decile_fraction": float(
                    event_stability.loc[
                        event_stability["baseline_module"].eq(module),
                        "bootstrap_hub_top_decile_fraction",
                    ].max()
                ),
                "validation_score": combined,
                "stability_class": _classification(combined, settings),
            }
        )
    module_stability = pd.DataFrame(
        module_rows,
        columns=[
            "module",
            "n_events",
            "mean_best_jaccard",
            "median_best_jaccard",
            "min_best_jaccard",
            "mean_event_membership_stability",
            "robust_event_fraction",
            "median_bootstrap_kME_sd",
            "max_hub_top_decile_fraction",
            "validation_score",
            "stability_class",
        ],
    )
    return event_stability, module_stability, pd.DataFrame(run_rows)


def _parameter_sensitivity(
    matrix: pd.DataFrame,
    baseline: pd.Series,
    network: NetworkConfig,
    settings: ValidationConfig,
    progress: Progress,
) -> pd.DataFrame:
    baseline = baseline.reindex(matrix.index).fillna("grey").astype(str)
    rows: list[dict[str, object]] = []
    if network.engine == "tom_hierarchical":
        heights = sorted(
            set([0.70, 0.75, network.module_cut_height, 0.85, 0.90])
        )
        for index, height in enumerate(heights, start=1):
            active = np.isclose(height, network.module_cut_height)
            labels = (
                baseline.copy()
                if active
                else pd.Series(
                    discover_module_labels(
                        matrix,
                        replace(network, module_cut_height=float(height)),
                    ),
                    index=matrix.index,
                    dtype=str,
                )
            )
            _mapping, matches = _match_modules(baseline, labels)
            optimal = matches.loc[matches["is_optimal_match"], "jaccard"]
            counts = labels.loc[labels.ne("grey")].value_counts()
            rows.append(
                {
                    "parameter_type": "tom_cut_height",
                    "parameter_value": float(height),
                    "module_cut_height": float(height),
                    "knn_neighbors": np.nan,
                    "community_resolution": np.nan,
                    "adjusted_rand_index": adjusted_rand_score(
                        baseline, labels
                    ),
                    "mean_matched_module_jaccard": (
                        float(optimal.mean()) if not optimal.empty else 0.0
                    ),
                    "min_matched_module_jaccard": (
                        float(optimal.min()) if not optimal.empty else 0.0
                    ),
                    "n_modules": int(len(counts)),
                    "n_grey_events": int(labels.eq("grey").sum()),
                    "largest_module_fraction": (
                        float(counts.max() / len(labels))
                        if not counts.empty
                        else 0.0
                    ),
                    "is_active_setting": active,
                }
            )
            if progress:
                progress(
                    f"TOM cut-height sensitivity {index}/{len(heights)}: "
                    f"height={height:g}"
                )
        return pd.DataFrame(rows)

    combinations = [(int(k), float(resolution)) for k in settings.knn_neighbors for resolution in settings.community_resolutions]
    for index, (neighbors, resolution) in enumerate(combinations, start=1):
        active = neighbors == network.knn_neighbors and np.isclose(resolution, network.community_resolution)
        if active:
            labels = baseline.copy()
        else:
            # Isolate the declared k/resolution perturbation instead of silently
            # reusing the active multi-k and bootstrap-consensus graph.
            variant = replace(
                network,
                knn_neighbors=neighbors,
                community_resolution=resolution,
                multi_k_consensus=False,
                consensus_bootstrap_iterations=0,
            )
            labels = pd.Series(discover_module_labels(matrix, variant, resolution), index=matrix.index, dtype=str)
        _mapping, matches = _match_modules(baseline, labels)
        optimal = matches.loc[matches["is_optimal_match"], "jaccard"]
        counts = labels.loc[labels.ne("grey")].value_counts()
        rows.append(
            {
                "knn_neighbors": neighbors,
                "community_resolution": resolution,
                "adjusted_rand_index": adjusted_rand_score(baseline, labels),
                "mean_matched_module_jaccard": float(optimal.mean()) if not optimal.empty else 0.0,
                "min_matched_module_jaccard": float(optimal.min()) if not optimal.empty else 0.0,
                "n_modules": int(len(counts)),
                "n_grey_events": int(labels.eq("grey").sum()),
                "largest_module_fraction": float(counts.max() / len(labels)) if not counts.empty else 0.0,
                "is_active_setting": active,
            }
        )
        if progress:
            progress(f"Parameter sensitivity {index}/{len(combinations)}: k={neighbors}, resolution={resolution:g}")
    return pd.DataFrame(rows)


def _program_score_confidence(
    matrix: pd.DataFrame,
    labels: pd.Series,
    iterations: int,
    random_seed: int,
) -> pd.DataFrame:
    """Estimate score uncertainty by resampling events within fixed modules."""

    labels = labels.reindex(matrix.index).fillna("grey").astype(str)
    baseline = _program_scores(matrix, labels.to_numpy(dtype=object))
    if baseline.empty:
        return pd.DataFrame()
    rng = np.random.default_rng(random_seed + 50_000)
    draws: dict[tuple[str, str], list[float]] = {
        (str(sample), str(module)): []
        for sample in matrix.columns
        for module in baseline.columns
    }
    for _ in range(iterations):
        for module in baseline.columns:
            members = labels.index[labels.eq(str(module))]
            if len(members) < 2:
                continue
            sampled = rng.choice(members.to_numpy(), size=len(members), replace=True)
            boot = matrix.loc[sampled].copy()
            boot.index = [f"event_{index}" for index in range(len(boot))]
            score = _program_scores(
                boot,
                np.asarray([str(module)] * len(boot), dtype=object),
            )[str(module)]
            if score.corr(baseline[str(module)]) < 0:
                score = -score
            for sample, value in score.items():
                draws[(str(sample), str(module))].append(float(value))
    rows: list[dict[str, object]] = []
    for (sample, module), values in draws.items():
        array = np.asarray(values, dtype=float)
        rows.append(
            {
                "sample_id": sample,
                "module": module,
                "program_score": float(baseline.loc[sample, module]),
                "composition_bootstrap_mean": float(np.mean(array)) if len(array) else np.nan,
                "composition_bootstrap_sd": float(np.std(array, ddof=1)) if len(array) > 1 else np.nan,
                "composition_bootstrap_CI_low": float(np.quantile(array, 0.025)) if len(array) else np.nan,
                "composition_bootstrap_CI_high": float(np.quantile(array, 0.975)) if len(array) else np.nan,
                "bootstrap_iterations": iterations,
                "uncertainty_scope": "event-composition uncertainty within fixed module",
            }
        )
    return pd.DataFrame(rows)


def _preservation(observed: pd.Series, residual: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame]:
    shared = observed.index.intersection(residual.index)
    observed = observed.reindex(shared).astype(str)
    residual = residual.reindex(shared).astype(str)
    observed_sets = _module_sets(observed)
    residual_sets = _module_sets(residual)
    rows: list[dict[str, object]] = []
    total = len(shared)
    for observed_module, observed_events in observed_sets.items():
        for residual_module, residual_events in residual_sets.items():
            overlap = len(observed_events & residual_events)
            union = len(observed_events | residual_events)
            pvalue = hypergeom.sf(overlap - 1, total, len(residual_events), len(observed_events))
            rows.append(
                {
                    "observed_module": observed_module,
                    "residual_module": residual_module,
                    "observed_events": len(observed_events),
                    "residual_events": len(residual_events),
                    "overlap_events": overlap,
                    "jaccard": overlap / union if union else 0.0,
                    "observed_recall": overlap / len(observed_events) if observed_events else 0.0,
                    "residual_precision": overlap / len(residual_events) if residual_events else 0.0,
                    "hypergeometric_pvalue": pvalue,
                }
            )
    pairs = pd.DataFrame(rows)
    if pairs.empty:
        return pairs, pd.DataFrame()
    pairs["FDR"] = multipletests(pairs["hypergeometric_pvalue"].fillna(1.0), method="fdr_bh")[1]
    pairs["is_best_jaccard_for_observed"] = pairs["jaccard"].eq(pairs.groupby("observed_module", observed=True)["jaccard"].transform("max"))
    best = pairs.loc[pairs["is_best_jaccard_for_observed"]].sort_values(
        ["observed_module", "jaccard", "FDR"], ascending=[True, False, True]
    ).drop_duplicates("observed_module")
    summary = best.rename(columns={"residual_module": "best_residual_module", "jaccard": "best_jaccard"})[
        [
            "observed_module",
            "best_residual_module",
            "observed_events",
            "residual_events",
            "overlap_events",
            "best_jaccard",
            "observed_recall",
            "residual_precision",
            "FDR",
        ]
    ].reset_index(drop=True)
    summary["observed_residual_program_class"] = np.select(
        [
            summary["FDR"].le(0.05) & summary["best_jaccard"].ge(0.75),
            summary["FDR"].le(0.05) & summary["best_jaccard"].ge(0.50),
            summary["FDR"].le(0.05),
        ],
        [
            "residual_preserved_intrinsic_program",
            "partially_residual_preserved_program",
            "weak_residual_overlap",
        ],
        default="condition_dominated_or_unstable_program",
    )
    summary["classification_note"] = (
        "Descriptive observed-vs-residual preservation; use continuous Jaccard, recall, "
        "precision and FDR for inference."
    )
    return pairs, summary


def validate_network_modules(
    matrix: pd.DataFrame,
    observed_assignments: pd.DataFrame,
    samples: pd.DataFrame,
    network: NetworkConfig,
    settings: ValidationConfig,
    residual_assignments: pd.DataFrame | None = None,
    progress: Progress = None,
) -> ValidationResult:
    observed = observed_assignments.set_index("event_key")["module"].reindex(matrix.index).fillna("grey").astype(str)
    event_stability, module_stability, bootstrap_runs = _bootstrap_validation(
        matrix, observed, samples, network, settings, progress
    )
    parameter_sensitivity = _parameter_sensitivity(matrix, observed, network, settings, progress)
    if residual_assignments is not None and not residual_assignments.empty:
        residual = residual_assignments.set_index("event_key")["module"]
        preservation_pairs, preservation_summary = _preservation(observed, residual)
    else:
        preservation_pairs = pd.DataFrame()
        preservation_summary = pd.DataFrame()
    program_score_confidence = _program_score_confidence(
        matrix,
        observed,
        settings.bootstrap_iterations,
        settings.random_seed,
    )
    return ValidationResult(
        event_stability,
        module_stability,
        bootstrap_runs,
        parameter_sensitivity,
        preservation_pairs,
        preservation_summary,
        program_score_confidence,
    )


def render_validation_plots(result: ValidationResult, output_dir: Path, settings: ValidationConfig) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    generated: list[Path] = []

    table = result.module_stability.sort_values("validation_score", ascending=False)
    fig, ax = plt.subplots(figsize=(max(9, len(table) * 0.8), 6))
    if table.empty:
        ax.text(
            0.5,
            0.5,
            "No non-grey module was available for bootstrap validation",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        ax.set_xticks([])
    else:
        colors = table["stability_class"].map({"robust": "#087E8B", "moderate": "#F3A712", "unstable": "#C44536"})
        ax.bar(table["module"], table["validation_score"], color=colors, edgecolor="#17324D")
    ax.axhline(settings.robust_jaccard, color="#087E8B", linestyle="--", linewidth=1.2, label="robust threshold")
    ax.axhline(settings.moderate_jaccard, color="#F3A712", linestyle="--", linewidth=1.2, label="moderate threshold")
    ax.set_ylim(0, 1.02)
    ax.set_ylabel("Validation score")
    ax.set_title("Co-splicing module bootstrap stability", fontweight="bold")
    ax.legend(loc="upper right")
    fig.tight_layout()
    path = output_dir / "module_bootstrap_stability.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    generated.append(path)

    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    sns.histplot(result.bootstrap_runs["adjusted_rand_index"], bins=min(12, len(result.bootstrap_runs)), color="#2D6CDF", ax=ax)
    ax.axvline(result.bootstrap_runs["adjusted_rand_index"].median(), color="#C44536", linewidth=2, label="median")
    ax.set_xlim(-0.05, 1.02)
    ax.set_xlabel("Adjusted Rand index vs. baseline modules")
    ax.set_title("Bootstrap recovery of the complete module partition", fontweight="bold")
    ax.legend()
    fig.tight_layout()
    path = output_dir / "bootstrap_adjusted_rand_index.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    generated.append(path)

    if (
        "parameter_type" in result.parameter_sensitivity
        and result.parameter_sensitivity["parameter_type"]
        .eq("tom_cut_height")
        .all()
    ):
        sensitivity = (
            result.parameter_sensitivity.set_index("module_cut_height")[
                "adjusted_rand_index"
            ]
            .to_frame()
            .T
        )
        sensitivity.index = ["TOM"]
        sensitivity_title = "Module sensitivity to TOM cut height"
        sensitivity_xlabel = "TOM dissimilarity cut height"
    else:
        sensitivity = result.parameter_sensitivity.pivot(
            index="community_resolution",
            columns="knn_neighbors",
            values="adjusted_rand_index",
        )
        sensitivity_title = "Module sensitivity to graph parameters"
        sensitivity_xlabel = "kNN neighbors"
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    sns.heatmap(sensitivity, annot=True, fmt=".2f", cmap="YlGnBu", vmin=0, vmax=1, ax=ax, cbar_kws={"label": "ARI vs. active network"})
    ax.set_title(sensitivity_title, fontweight="bold")
    ax.set_xlabel(sensitivity_xlabel)
    fig.tight_layout()
    path = output_dir / "parameter_sensitivity_heatmap.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    generated.append(path)

    if not result.preservation_pairs.empty:
        preservation = result.preservation_pairs.pivot(index="observed_module", columns="residual_module", values="jaccard")
        fig, ax = plt.subplots(figsize=(max(8, preservation.shape[1] * 0.75), max(6, preservation.shape[0] * 0.65)))
        sns.heatmap(preservation, annot=True, fmt=".2f", cmap="mako", vmin=0, vmax=1, ax=ax, cbar_kws={"label": "Event-set Jaccard"})
        ax.set_title("Observed-to-residual module preservation", fontweight="bold")
        fig.tight_layout()
        path = output_dir / "observed_residual_module_preservation.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        generated.append(path)
    return generated


def render_validation_report(result: ValidationResult, output_path: Path) -> Path:
    counts = result.module_stability["stability_class"].value_counts()
    rows = "".join(
        "<tr>" + "".join(f"<td>{escape(str(row[column]))}</td>" for column in ["module", "n_events", "validation_score", "stability_class"]) + "</tr>"
        for _, row in result.module_stability.sort_values("validation_score", ascending=False).iterrows()
    )
    preservation = "Not available because a cached residual network was not supplied."
    if not result.preservation_summary.empty:
        preservation = f"{len(result.preservation_summary)} observed modules were compared with residual-network modules."
    html = f"""<!doctype html><html><head><meta charset="utf-8"><title>CoSplice module validation</title>
<style>body{{font-family:Arial,sans-serif;margin:36px;max-width:1100px;color:#17212b}}table{{border-collapse:collapse;width:100%}}td,th{{border:1px solid #cbd5df;padding:8px}}th{{background:#e9f1f2}}h1{{color:#087e8b}}.note{{background:#f5f8f8;border-left:4px solid #087e8b;padding:12px}}</style></head><body>
<h1>CoSplice module validation</h1>
<p class="note">Bootstrap and parameter sensitivity evaluate recovery of the discovery partition in the same cohort. They measure robustness, not external biological validation.</p>
<p>Robust modules: <strong>{int(counts.get('robust', 0))}</strong>; moderate: <strong>{int(counts.get('moderate', 0))}</strong>; unstable: <strong>{int(counts.get('unstable', 0))}</strong>.</p>
<p>Median bootstrap ARI: <strong>{result.bootstrap_runs['adjusted_rand_index'].median():.3f}</strong>. {escape(preservation)}</p>
<table><thead><tr><th>Module</th><th>Events</th><th>Validation score</th><th>Class</th></tr></thead><tbody>{rows}</tbody></table>
</body></html>"""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html, encoding="utf-8")
    return output_path
