from __future__ import annotations

from dataclasses import dataclass, replace
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.covariance import LedoitWolf, graphical_lasso
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LinearRegression
from sklearn.metrics import adjusted_rand_score, silhouette_score
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

from katana_cosplice.config import BenchmarkConfig, NetworkConfig
from katana_cosplice.network import correlation_matrix, discover_module_labels


METHOD_LABELS = {
    "katana_sparse": "Katana-CoSplice sparse signed network",
    "pearson_hierarchical": "Pearson signed hierarchical baseline",
    "graphical_lasso": "Graphical lasso partial-correlation network",
    "wgcna_exact": "WGCNA blockwiseModules (exact R package)",
}


@dataclass
class BenchmarkResult:
    universe: pd.DataFrame
    assignments: pd.DataFrame
    method_summary: pd.DataFrame
    bootstrap_stability: pd.DataFrame
    module_matching: pd.DataFrame
    method_audit: pd.DataFrame
    metric_ranking: pd.DataFrame
    method_interpretation: pd.DataFrame


def _module_labels(raw_labels: np.ndarray, min_size: int) -> np.ndarray:
    counts = pd.Series(raw_labels).value_counts()
    retained = [label for label, count in counts.items() if int(count) >= min_size and str(label) != "grey"]
    retained.sort(key=lambda label: (-int(counts[label]), str(label)))
    names = {label: f"M{index:02d}" for index, label in enumerate(retained, start=1)}
    return np.asarray([names.get(label, "grey") for label in raw_labels], dtype=object)


def select_benchmark_universe(
    stable_psi: pd.DataFrame,
    qc_metrics: pd.DataFrame,
    max_events: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if not stable_psi.index.is_unique or not stable_psi.columns.is_unique:
        raise ValueError("Benchmark matrix requires unique event and sample identifiers.")
    if stable_psi.isna().any().any() or not np.isfinite(stable_psi.to_numpy(dtype=float)).all():
        raise ValueError("Benchmark matrix must be finite.")
    metrics = qc_metrics.copy()
    if "event_key" in metrics.columns:
        metrics = metrics.set_index("event_key", drop=True)
    if "stability_score" not in metrics.columns:
        raise ValueError("Stable-AS QC metrics must contain stability_score.")
    ranking = metrics.reindex(stable_psi.index)["stability_score"].fillna(0.0).sort_values(ascending=False)
    selected_keys = ranking.index[: min(max_events, len(ranking))]
    selected = stable_psi.loc[selected_keys].copy()
    universe = metrics.reindex(selected_keys).reset_index().rename(columns={metrics.index.name or "index": "event_key"})
    if "event_key" not in universe.columns:
        universe = universe.rename(columns={universe.columns[0]: "event_key"})
    universe.insert(1, "benchmark_rank", np.arange(1, len(universe) + 1))
    universe["selection_rule"] = f"top_{max_events}_stable_events_by_precomputed_stability_score"
    return selected, universe


def _fit_katana(matrix: pd.DataFrame, network: NetworkConfig) -> np.ndarray:
    settings = replace(
        network,
        engine="sparse_knn",
        max_events=None,
        knn_neighbors=min(network.knn_neighbors, max(2, len(matrix) - 1)),
        bootstrap_iterations=0,
    )
    return discover_module_labels(matrix, settings)


def _fit_hierarchical(matrix: pd.DataFrame, network: NetworkConfig) -> np.ndarray:
    correlation = correlation_matrix(matrix, "pearson")
    distance = 1.0 - correlation if network.signed else 1.0 - np.abs(correlation)
    distance = np.clip((distance + distance.T) / 2.0, 0.0, 2.0)
    np.fill_diagonal(distance, 0.0)
    hierarchy = linkage(squareform(distance, checks=False), method="average")
    raw = fcluster(hierarchy, t=network.module_cut_height, criterion="distance")
    return _module_labels(raw, network.min_module_size)


def _fit_graphical_lasso(
    matrix: pd.DataFrame,
    network: NetworkConfig,
    settings: BenchmarkConfig,
) -> tuple[np.ndarray, str]:
    values = StandardScaler().fit_transform(matrix.T.to_numpy(dtype=float))
    covariance_model = LedoitWolf().fit(values)
    empirical_covariance = covariance_model.covariance_
    base_alpha = settings.graphical_lasso_alpha
    alpha_candidates = sorted(set([base_alpha, base_alpha * 2.0, base_alpha * 3.0, base_alpha * 5.0, 1.0]))
    precision = None
    caught: list[warnings.WarningMessage] = []
    failures: list[str] = []
    effective_alpha = base_alpha
    for effective_alpha in alpha_candidates:
        try:
            with warnings.catch_warnings(record=True) as attempt_warnings:
                warnings.simplefilter("always", ConvergenceWarning)
                _covariance, precision = graphical_lasso(
                    empirical_covariance,
                    alpha=effective_alpha,
                    max_iter=settings.graphical_lasso_max_iter,
                    tol=1e-4,
                )
            caught = list(attempt_warnings)
            break
        except (FloatingPointError, np.linalg.LinAlgError) as exc:
            failures.append(f"alpha={effective_alpha:g}: {exc}")
    if precision is None:
        raise RuntimeError("Graphical lasso failed across the declared alpha sequence: " + " | ".join(failures))
    scale = np.sqrt(np.outer(np.diag(precision), np.diag(precision)))
    partial = -precision / np.where(scale <= 0, 1.0, scale)
    np.fill_diagonal(partial, 0.0)
    weights = np.maximum(partial, 0.0) if network.signed else np.abs(partial)
    graph = nx.Graph()
    graph.add_nodes_from(range(len(matrix)))
    rows, cols = np.where(np.triu(weights >= network.edge_threshold, k=1))
    graph.add_weighted_edges_from((int(row), int(col), float(weights[row, col])) for row, col in zip(rows, cols, strict=True))
    if graph.number_of_edges() == 0:
        raw = np.arange(len(matrix), dtype=int)
    else:
        communities = nx.community.louvain_communities(
            graph,
            weight="weight",
            resolution=network.community_resolution,
            seed=settings.random_seed,
        )
        raw = np.zeros(len(matrix), dtype=int)
        for label, community in enumerate(communities, start=1):
            raw[np.fromiter(community, dtype=int)] = label
    warning_parts = []
    warning_parts.append(f"Ledoit-Wolf covariance shrinkage={covariance_model.shrinkage_:.4g}")
    if not np.isclose(effective_alpha, base_alpha):
        warning_parts.append(f"adaptive regularization used alpha={effective_alpha:g} after alpha={base_alpha:g} was ill-conditioned")
    convergence = [str(item.message) for item in caught if issubclass(item.category, ConvergenceWarning)]
    if convergence:
        warning_parts.append(f"{len(convergence)} convergence warnings; first: {convergence[0]}")
    warning_text = "; ".join(warning_parts)
    return _module_labels(raw, network.min_module_size), warning_text


def _rscript_candidates() -> list[Path]:
    candidates: list[Path] = []
    configured = os.environ.get("KATANA_RSCRIPT", "").strip()
    if configured:
        candidates.append(Path(configured))
    on_path = shutil.which("Rscript")
    if on_path:
        candidates.append(Path(on_path))
    for prefix in filter(None, [os.environ.get("CONDA_PREFIX"), sys.prefix]):
        root = Path(prefix)
        candidates.extend([root / "Scripts" / "Rscript.exe", root / "Lib" / "R" / "bin" / "Rscript.exe"])
    conda_exe = os.environ.get("CONDA_EXE", "").strip()
    if conda_exe:
        conda_root = Path(conda_exe).resolve().parent.parent
        candidates.extend(sorted((conda_root / "envs").glob("*/Scripts/Rscript.exe")))
    unique: list[Path] = []
    seen: set[str] = set()
    for candidate in candidates:
        key = str(candidate).casefold()
        if key not in seen and candidate.exists():
            unique.append(candidate)
            seen.add(key)
    return unique


def _r_subprocess_env(executable: str | Path) -> dict[str, str]:
    env = os.environ.copy()
    candidate = Path(executable).resolve()
    conda_prefix = next((parent for parent in candidate.parents if (parent / "conda-meta").is_dir()), None)
    if conda_prefix is not None:
        additions = [
            conda_prefix,
            conda_prefix / "Library" / "bin",
            conda_prefix / "Scripts",
            conda_prefix / "Lib" / "R" / "bin",
        ]
        env["PATH"] = os.pathsep.join([*(str(path) for path in additions), env.get("PATH", "")])
    return env


def _rscript_with_wgcna() -> tuple[str | None, str]:
    diagnostics: list[str] = []
    for candidate in _rscript_candidates():
        check = subprocess.run(
            [
                str(candidate),
                "-e",
                "cat(ifelse(requireNamespace('WGCNA', quietly=TRUE), as.character(packageVersion('WGCNA')), 'MISSING'))",
            ],
            capture_output=True,
            text=True,
            timeout=60,
            check=False,
            env=_r_subprocess_env(candidate),
        )
        version = check.stdout.strip()
        if check.returncode == 0 and version and version != "MISSING":
            return str(candidate), f"R WGCNA {version}; Rscript={candidate}"
        diagnostics.append(f"{candidate}: {version or check.stderr.strip() or 'WGCNA unavailable'}")
    if not diagnostics:
        return None, "Rscript was not found on PATH, in the active environment, or in sibling Conda environments."
    return None, "Rscript candidates were found, but exact WGCNA was unavailable: " + " | ".join(diagnostics)


def _fit_wgcna(
    matrix: pd.DataFrame,
    network: NetworkConfig,
    settings: BenchmarkConfig,
    executable: str,
) -> np.ndarray:
    with tempfile.TemporaryDirectory(prefix="katana_cosplice_wgcna_") as temp_dir:
        temp = Path(temp_dir)
        matrix_path = temp / "matrix.csv"
        output_path = temp / "modules.csv"
        script_path = temp / "run_wgcna.R"
        matrix.T.to_csv(matrix_path, index=True, index_label="sample_id")
        script_path.write_text(
            """
args <- commandArgs(trailingOnly=TRUE)
suppressPackageStartupMessages(library(WGCNA))
dat <- read.csv(args[1], row.names=1, check.names=FALSE)
net <- blockwiseModules(
  dat,
  power=as.numeric(args[3]),
  corType='bicor',
  networkType='signed',
  TOMType='signed',
  minModuleSize=as.integer(args[4]),
  mergeCutHeight=as.numeric(args[5]),
  numericLabels=TRUE,
  pamRespectsDendro=FALSE,
  maxBlockSize=ncol(dat),
  verbose=0
)
write.csv(data.frame(event_key=colnames(dat), raw_module=net$colors), args[2], row.names=FALSE, quote=TRUE)
""".strip(),
            encoding="utf-8",
        )
        completed = subprocess.run(
            [
                executable,
                str(script_path),
                str(matrix_path),
                str(output_path),
                str(network.soft_power),
                str(network.min_module_size),
                str(settings.wgcna_merge_cut_height),
            ],
            capture_output=True,
            text=True,
            timeout=1800,
            check=False,
            env=_r_subprocess_env(executable),
        )
        if completed.returncode != 0 or not output_path.exists():
            raise RuntimeError(f"WGCNA failed: {completed.stderr.strip() or completed.stdout.strip()}")
        result = pd.read_csv(output_path).set_index("event_key").reindex(matrix.index)
        if result["raw_module"].isna().any():
            raise RuntimeError("WGCNA output did not preserve the common event universe.")
        raw = result["raw_module"].astype(str).replace("0", "grey").to_numpy(dtype=object)
        return _module_labels(raw, network.min_module_size)


def _fit_method(
    method: str,
    matrix: pd.DataFrame,
    network: NetworkConfig,
    settings: BenchmarkConfig,
    wgcna_executable: str | None,
) -> tuple[np.ndarray, str]:
    if method == "katana_sparse":
        return _fit_katana(matrix, network), ""
    if method == "pearson_hierarchical":
        return _fit_hierarchical(matrix, network), ""
    if method == "graphical_lasso":
        return _fit_graphical_lasso(matrix, network, settings)
    if method == "wgcna_exact":
        if wgcna_executable is None:
            raise FileNotFoundError("Exact WGCNA is unavailable; see method_availability.tsv.")
        return _fit_wgcna(matrix, network, settings, wgcna_executable), ""
    raise ValueError(f"Unsupported benchmark method: {method}")


def _neutral_graph(matrix: pd.DataFrame, neighbors: int) -> nx.Graph:
    values = matrix.to_numpy(dtype=float)
    values -= values.mean(axis=1, keepdims=True)
    norms = np.sqrt(np.sum(values**2, axis=1, keepdims=True))
    profiles = values / np.where(norms < 1e-12, 1.0, norms)
    model = NearestNeighbors(
        n_neighbors=min(len(matrix), neighbors + 1),
        metric="cosine",
        algorithm="brute",
        n_jobs=-1,
    ).fit(profiles)
    distances, indices = model.kneighbors(profiles)
    graph = nx.Graph()
    graph.add_nodes_from(range(len(matrix)))
    for source in range(len(matrix)):
        for distance, target in zip(distances[source], indices[source], strict=True):
            target = int(target)
            if source == target:
                continue
            weight = max(0.0, 1.0 - float(distance))
            if weight > 0:
                graph.add_edge(source, target, weight=weight)
    return graph


def _program_score(values: np.ndarray) -> np.ndarray:
    centered = values - values.mean(axis=1, keepdims=True)
    scale = centered.std(axis=1, ddof=1, keepdims=True)
    standardized = centered / np.where(scale < 1e-12, 1.0, scale)
    _, _, vt = np.linalg.svd(standardized, full_matrices=False)
    score = vt[0]
    if np.corrcoef(score, standardized.mean(axis=0))[0, 1] < 0:
        score = -score
    return score


def _method_metrics(
    matrix: pd.DataFrame,
    labels: np.ndarray,
    samples: pd.DataFrame,
    neutral_graph: nx.Graph,
) -> dict[str, float | int]:
    label_series = pd.Series(labels, index=matrix.index)
    retained = label_series.ne("grey")
    retained_labels = label_series.loc[retained]
    sizes = retained_labels.value_counts()
    n_modules = int(len(sizes))
    largest_fraction = float(sizes.max() / len(matrix)) if not sizes.empty else 0.0
    silhouette = float("nan")
    if n_modules >= 2 and retained.sum() > n_modules:
        usable = retained_labels[retained_labels.map(sizes).ge(2)]
        if usable.nunique() >= 2 and len(usable) > usable.nunique():
            silhouette = float(silhouette_score(matrix.loc[usable.index], usable.to_numpy(), metric="cosine"))

    correlation = np.corrcoef(matrix.to_numpy(dtype=float))
    upper = np.triu_indices_from(correlation, k=1)
    same = labels[upper[0]] == labels[upper[1]]
    non_grey_pair = (labels[upper[0]] != "grey") & (labels[upper[1]] != "grey")
    within_values = correlation[upper][same & non_grey_pair]
    between_values = correlation[upper][~same & non_grey_pair]
    within = float(np.nanmean(within_values)) if len(within_values) else float("nan")
    between = float(np.nanmean(between_values)) if len(between_values) else float("nan")

    communities = [set(np.flatnonzero(labels == module)) for module in pd.unique(labels)]
    modularity = float(nx.community.modularity(neutral_graph, communities, weight="weight")) if neutral_graph.number_of_edges() else float("nan")

    condition_r2: list[tuple[float, int]] = []
    if "condition" in samples.columns and "sample_id" in samples.columns:
        ordered = samples.set_index("sample_id").reindex(matrix.columns)
        design = pd.get_dummies(ordered["condition"].astype(str), drop_first=False, dtype=float).to_numpy()
        for module, size in sizes.items():
            score = _program_score(matrix.loc[label_series.eq(module)].to_numpy(dtype=float))
            condition_r2.append((float(LinearRegression().fit(design, score).score(design, score)), int(size)))
    weighted_r2 = (
        float(np.average([value for value, _ in condition_r2], weights=[size for _, size in condition_r2]))
        if condition_r2
        else float("nan")
    )
    median_r2 = float(np.median([value for value, _ in condition_r2])) if condition_r2 else float("nan")
    return {
        "n_modules": n_modules,
        "n_grey_events": int((~retained).sum()),
        "grey_fraction": float((~retained).mean()),
        "largest_module_fraction": largest_fraction,
        "silhouette_cosine": silhouette,
        "neutral_graph_modularity": modularity,
        "within_module_mean_pearson": within,
        "between_module_mean_pearson": between,
        "within_between_separation": within - between,
        "weighted_condition_pattern_R2": weighted_r2,
        "median_condition_pattern_R2": median_r2,
    }


def _module_matching(assignments: pd.DataFrame) -> pd.DataFrame:
    methods = [column.removeprefix("module_") for column in assignments if column.startswith("module_")]
    rows: list[dict[str, object]] = []
    for reference_method in methods:
        reference = assignments.set_index("event_key")[f"module_{reference_method}"]
        for candidate_method in methods:
            if reference_method == candidate_method:
                continue
            candidate = assignments.set_index("event_key")[f"module_{candidate_method}"]
            for reference_module in sorted(set(reference).difference({"grey"})):
                reference_events = set(reference.index[reference.eq(reference_module)])
                best_module = ""
                best_overlap = 0
                best_jaccard = 0.0
                best_size = 0
                for candidate_module in sorted(set(candidate).difference({"grey"})):
                    candidate_events = set(candidate.index[candidate.eq(candidate_module)])
                    overlap = len(reference_events.intersection(candidate_events))
                    union = len(reference_events.union(candidate_events))
                    jaccard = overlap / union if union else 0.0
                    if jaccard > best_jaccard:
                        best_module = candidate_module
                        best_overlap = overlap
                        best_jaccard = jaccard
                        best_size = len(candidate_events)
                rows.append(
                    {
                        "reference_method": reference_method,
                        "candidate_method": candidate_method,
                        "reference_module": reference_module,
                        "candidate_best_module": best_module,
                        "reference_module_size": len(reference_events),
                        "candidate_module_size": best_size,
                        "overlap_events": best_overlap,
                        "best_jaccard": best_jaccard,
                    }
                )
    return pd.DataFrame(rows)


BENCHMARK_SCORE_METRICS = {
    "silhouette_cosine": True,
    "neutral_graph_modularity": True,
    "within_between_separation": True,
    "mean_bootstrap_ARI": True,
    "runtime_seconds": False,
}


def _metric_ranking(summary: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    valid = summary.loc[summary.get("valid_partition", False).fillna(False).astype(bool)].copy()
    for metric, higher_is_better in BENCHMARK_SCORE_METRICS.items():
        if metric not in valid:
            continue
        values = pd.to_numeric(valid[metric], errors="coerce")
        usable = values.dropna()
        if usable.empty:
            continue
        minimum = float(usable.min())
        maximum = float(usable.max())
        if np.isclose(minimum, maximum):
            normalized = pd.Series(1.0, index=usable.index)
        elif higher_is_better:
            normalized = (usable - minimum) / (maximum - minimum)
        else:
            normalized = (maximum - usable) / (maximum - minimum)
        ranks = usable.rank(method="min", ascending=not higher_is_better)
        for index in usable.index:
            rows.append(
                {
                    "method": valid.loc[index, "method"],
                    "metric": metric,
                    "raw_value": float(usable.loc[index]),
                    "higher_is_better": higher_is_better,
                    "metric_rank": int(ranks.loc[index]),
                    "normalized_descriptive_score": float(normalized.loc[index]),
                }
            )
    ranking = pd.DataFrame(rows)
    score_lookup = (
        ranking.groupby("method", observed=True)["normalized_descriptive_score"].mean()
        if not ranking.empty
        else pd.Series(dtype=float)
    )
    rank_lookup = score_lookup.rank(method="min", ascending=False)
    metric_winners = (
        ranking.loc[ranking["metric_rank"].eq(1)].groupby("method", observed=True)["metric"].apply(lambda values: ", ".join(values))
        if not ranking.empty
        else pd.Series(dtype=object)
    )
    metric_non_winners = (
        ranking.loc[ranking["metric_rank"].gt(1)].groupby("method", observed=True)["metric"].apply(lambda values: ", ".join(values))
        if not ranking.empty
        else pd.Series(dtype=object)
    )
    interpretation_rows: list[dict[str, object]] = []
    for row in summary.itertuples(index=False):
        method = str(row.method)
        status = str(row.status)
        strengths = metric_winners.get(method, "")
        limitations: list[str] = []
        if status == "degenerate":
            limitations.append("Partition retained fewer than two non-grey modules; quality ranking is not interpretable.")
        elif status in {"failed", "unavailable"}:
            limitations.append("Method did not produce a comparable partition in this run.")
        else:
            grey_fraction = float(getattr(row, "grey_fraction", np.nan))
            largest = float(getattr(row, "largest_module_fraction", np.nan))
            if np.isfinite(grey_fraction) and grey_fraction > 0.25:
                limitations.append(f"High grey-event fraction ({grey_fraction:.1%}).")
            if np.isfinite(largest) and largest > 0.75:
                limitations.append(f"Largest module contains {largest:.1%} of the shared event universe.")
            non_winners = metric_non_winners.get(method, "")
            if non_winners:
                limitations.append(f"Ranked below the best valid method on: {non_winners}.")
            elif not strengths:
                limitations.append("Did not rank first on the declared neutral scorecard metrics.")
        interpretation_rows.append(
            {
                "method": method,
                "status": status,
                "valid_partition": bool(getattr(row, "valid_partition", False)),
                "descriptive_composite_score": score_lookup.get(method, np.nan),
                "descriptive_overall_rank": int(rank_lookup[method]) if method in rank_lookup else np.nan,
                "strengths_in_this_dataset": strengths or "No valid metric win recorded.",
                "limitations_in_this_dataset": " ".join(limitations) or "No automatic structural warning; inspect raw metrics and biology.",
                "claim_scope": "Internal same-cohort computational comparison; not external biological validation.",
            }
        )
    return ranking, pd.DataFrame(interpretation_rows)


def run_method_benchmark(
    stable_psi: pd.DataFrame,
    qc_metrics: pd.DataFrame,
    samples: pd.DataFrame,
    network: NetworkConfig,
    settings: BenchmarkConfig,
    progress=None,
) -> BenchmarkResult:
    matrix, universe = select_benchmark_universe(stable_psi, qc_metrics, settings.max_events)
    if len(matrix) < network.min_module_size * 2:
        raise ValueError("Too few stable events for a shared benchmark universe.")
    if set(matrix.columns) != set(samples["sample_id"]):
        raise ValueError("Benchmark sample metadata does not match the stable PSI matrix.")
    samples = samples.set_index("sample_id").reindex(matrix.columns)
    samples.index.name = "sample_id"
    samples = samples.reset_index()
    neutral_graph = _neutral_graph(matrix, settings.neutral_knn_neighbors)
    wgcna_executable, wgcna_note = _rscript_with_wgcna()

    method_labels: dict[str, np.ndarray] = {}
    audit_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for method in settings.methods:
        if progress:
            progress(f"Benchmark fitting {METHOD_LABELS[method]} on {len(matrix)} shared events.")
        if method == "wgcna_exact" and wgcna_executable is None:
            audit_rows.append(
                {
                    "method": method,
                    "display_name": METHOD_LABELS[method],
                    "status": "unavailable",
                    "implementation": "exact external R WGCNA",
                    "message": wgcna_note,
                    "runtime_seconds": np.nan,
                }
            )
            summary_rows.append({"method": method, "status": "unavailable"})
            continue
        started = time.perf_counter()
        try:
            labels, note = _fit_method(method, matrix, network, settings, wgcna_executable)
            runtime = time.perf_counter() - started
            method_labels[method] = labels
            metrics = _method_metrics(matrix, labels, samples, neutral_graph)
            method_status = "finished" if int(metrics["n_modules"]) >= 2 else "degenerate"
            summary_rows.append(
                {
                    "method": method,
                    "status": method_status,
                    "valid_partition": method_status == "finished",
                    "n_events": len(matrix),
                    "n_samples": matrix.shape[1],
                    "runtime_seconds": runtime,
                    **metrics,
                }
            )
            audit_rows.append(
                {
                    "method": method,
                    "display_name": METHOD_LABELS[method],
                    "status": "degenerate" if method_status == "degenerate" else ("finished_with_warning" if note else "finished"),
                    "implementation": "exact R WGCNA" if method == "wgcna_exact" else "Katana Python implementation",
                    "message": note or (wgcna_note if method == "wgcna_exact" else ""),
                    "runtime_seconds": runtime,
                }
            )
        except Exception as exc:
            runtime = time.perf_counter() - started
            audit_rows.append(
                {
                    "method": method,
                    "display_name": METHOD_LABELS[method],
                    "status": "failed",
                    "implementation": "exact R WGCNA" if method == "wgcna_exact" else "Katana Python implementation",
                    "message": str(exc),
                    "runtime_seconds": runtime,
                }
            )
            summary_rows.append({"method": method, "status": "failed", "runtime_seconds": runtime})

    if "katana_sparse" not in method_labels:
        raise RuntimeError("Katana-CoSplice benchmark method failed; comparison cannot continue.")
    assignments = pd.DataFrame({"event_key": matrix.index})
    for method, labels in method_labels.items():
        assignments[f"module_{method}"] = labels

    bootstrap_rows: list[dict[str, object]] = []
    rng = np.random.default_rng(settings.random_seed)
    condition_values = samples["condition"].astype(str).to_numpy() if "condition" in samples else np.repeat("all", len(samples))
    groups = [np.flatnonzero(condition_values == condition) for condition in pd.unique(condition_values)]
    for iteration in range(1, settings.bootstrap_iterations + 1):
        sampled = np.concatenate([rng.choice(group, size=len(group), replace=True) for group in groups])
        bootstrap_matrix = pd.DataFrame(
            matrix.to_numpy(dtype=float)[:, sampled],
            index=matrix.index,
            columns=[f"bootstrap_{iteration:03d}_{index:03d}" for index in range(len(sampled))],
        )
        for method, reference_labels in method_labels.items():
            started = time.perf_counter()
            reference_modules = set(reference_labels).difference({"grey"})
            if len(reference_modules) < 2:
                bootstrap_rows.append(
                    {
                        "iteration": iteration,
                        "method": method,
                        "adjusted_rand_index": np.nan,
                        "runtime_seconds": 0.0,
                        "status": "not_applicable_degenerate_reference",
                        "error": "Reference partition has fewer than two retained modules.",
                    }
                )
                continue
            try:
                labels, _ = _fit_method(method, bootstrap_matrix, network, settings, wgcna_executable)
                bootstrap_rows.append(
                    {
                        "iteration": iteration,
                        "method": method,
                        "adjusted_rand_index": adjusted_rand_score(reference_labels, labels),
                        "runtime_seconds": time.perf_counter() - started,
                        "status": "finished",
                        "error": "",
                    }
                )
            except Exception as exc:
                bootstrap_rows.append(
                    {
                        "iteration": iteration,
                        "method": method,
                        "adjusted_rand_index": np.nan,
                        "runtime_seconds": time.perf_counter() - started,
                        "status": "failed",
                        "error": str(exc),
                    }
                )
        if progress:
            progress(f"Benchmark bootstrap {iteration}/{settings.bootstrap_iterations} finished.")

    bootstrap = pd.DataFrame(bootstrap_rows)
    summary = pd.DataFrame(summary_rows)
    if not bootstrap.empty:
        bootstrap_summary = bootstrap.groupby("method", observed=True)["adjusted_rand_index"].agg(
            mean_bootstrap_ARI="mean",
            median_bootstrap_ARI="median",
            min_bootstrap_ARI="min",
        )
        summary = summary.merge(bootstrap_summary, left_on="method", right_index=True, how="left")
    metric_ranking, method_interpretation = _metric_ranking(summary)
    return BenchmarkResult(
        universe=universe,
        assignments=assignments,
        method_summary=summary,
        bootstrap_stability=bootstrap,
        module_matching=_module_matching(assignments),
        method_audit=pd.DataFrame(audit_rows),
        metric_ranking=metric_ranking,
        method_interpretation=method_interpretation,
    )


def render_benchmark_plots(result: BenchmarkResult, output_dir: Path) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    generated: list[Path] = []
    finished = result.method_summary.loc[result.method_summary["status"].isin(["finished", "degenerate"])].copy()
    if not finished.empty:
        metrics = [
            "silhouette_cosine",
            "neutral_graph_modularity",
            "within_between_separation",
            "weighted_condition_pattern_R2",
            "mean_bootstrap_ARI",
        ]
        available = [column for column in metrics if column in finished and finished[column].notna().any()]
        if available:
            heatmap = finished.set_index("method")[available].apply(pd.to_numeric, errors="coerce")
            fig, ax = plt.subplots(figsize=(max(9, len(available) * 1.7), max(4.5, len(heatmap) * 0.8)))
            sns.heatmap(heatmap, annot=True, fmt=".3f", cmap="crest", ax=ax, cbar_kws={"label": "Raw metric value"})
            ax.set_title("Co-splicing method benchmark: neutral quality and stability metrics", fontweight="bold")
            ax.set_xlabel("Metric (higher is better within this benchmark)")
            ax.set_ylabel("Method")
            fig.tight_layout()
            path = output_dir / "method_quality_heatmap.png"
            fig.savefig(path, dpi=220, bbox_inches="tight")
            plt.close(fig)
            generated.append(path)

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        axes[0].bar(finished["method"], finished["n_modules"], color="#087e8b")
        axes[0].set_ylabel("Retained modules")
        axes[0].set_title("Module count")
        axes[1].bar(finished["method"], finished["largest_module_fraction"], color="#ef8354")
        axes[1].set_ylabel("Fraction of shared event universe")
        axes[1].set_title("Largest-module fraction")
        for ax in axes:
            ax.tick_params(axis="x", rotation=25)
            for label in ax.get_xticklabels():
                label.set_horizontalalignment("right")
        fig.suptitle("Head-to-head module structure on the same stable-AS events", fontweight="bold")
        fig.tight_layout()
        path = output_dir / "method_module_structure.png"
        fig.savefig(path, dpi=220, bbox_inches="tight")
        plt.close(fig)
        generated.append(path)

    if not result.module_matching.empty:
        matching = result.module_matching.groupby(["reference_method", "candidate_method"], observed=True)["best_jaccard"].mean().unstack()
        fig, ax = plt.subplots(figsize=(8, 6))
        sns.heatmap(matching, annot=True, fmt=".2f", cmap="mako", vmin=0, vmax=1, ax=ax, cbar_kws={"label": "Mean best-module Jaccard"})
        ax.set_title("Cross-method module membership agreement", fontweight="bold")
        fig.tight_layout()
        path = output_dir / "cross_method_module_agreement.png"
        fig.savefig(path, dpi=220, bbox_inches="tight")
        plt.close(fig)
        generated.append(path)

    valid = result.method_summary.loc[result.method_summary["status"].eq("finished")].copy()
    if not valid.empty and not result.metric_ranking.empty:
        fig, axes = plt.subplots(2, 2, figsize=(15, 10))
        normalized = result.metric_ranking.pivot(index="method", columns="metric", values="normalized_descriptive_score")
        sns.heatmap(
            normalized,
            annot=True,
            fmt=".2f",
            cmap="YlGnBu",
            vmin=0,
            vmax=1,
            ax=axes[0, 0],
            cbar_kws={"label": "Within-benchmark normalized score"},
        )
        axes[0, 0].set_title("Neutral metric scorecard", fontweight="bold")
        axes[0, 0].set_xlabel("")
        axes[0, 0].set_ylabel("")

        bootstrap = result.bootstrap_stability.loc[
            result.bootstrap_stability["status"].eq("finished")
            & result.bootstrap_stability["method"].isin(valid["method"])
        ]
        if not bootstrap.empty:
            sns.boxplot(data=bootstrap, x="method", y="adjusted_rand_index", color="#64b6ac", ax=axes[0, 1])
            sns.stripplot(data=bootstrap, x="method", y="adjusted_rand_index", color="#17212b", size=4, ax=axes[0, 1])
        axes[0, 1].set_ylim(-0.05, 1.05)
        axes[0, 1].set_title("Condition-stratified bootstrap stability", fontweight="bold")
        axes[0, 1].set_xlabel("")
        axes[0, 1].set_ylabel("Adjusted Rand index")

        axes[1, 0].scatter(
            valid["largest_module_fraction"],
            valid["neutral_graph_modularity"],
            s=160,
            c=valid["silhouette_cosine"],
            cmap="viridis",
            edgecolor="white",
            linewidth=1.2,
        )
        for row in valid.itertuples(index=False):
            axes[1, 0].annotate(row.method, (row.largest_module_fraction, row.neutral_graph_modularity), xytext=(6, 5), textcoords="offset points")
        axes[1, 0].set_xlabel("Largest-module fraction (lower avoids collapse)")
        axes[1, 0].set_ylabel("Neutral KNN graph modularity")
        axes[1, 0].set_title("Partition structure", fontweight="bold")

        runtime = valid.sort_values("runtime_seconds")
        axes[1, 1].barh(runtime["method"], runtime["runtime_seconds"], color="#ef8354")
        axes[1, 1].set_xlabel("Runtime (seconds; lower is better)")
        axes[1, 1].set_title("Observed runtime on the same input", fontweight="bold")
        for ax in axes.flat:
            ax.tick_params(axis="x", rotation=20)
        fig.suptitle(
            "CoSplice head-to-head dashboard\nSame stable-AS events and samples; score is descriptive, not external validation",
            fontweight="bold",
            fontsize=15,
        )
        fig.tight_layout(rect=[0, 0, 1, 0.94])
        path = output_dir / "head_to_head_scorecard.png"
        fig.savefig(path, dpi=220, bbox_inches="tight")
        plt.close(fig)
        generated.append(path)
    return generated


def render_benchmark_report(result: BenchmarkResult, output_path: Path) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    summary = result.method_summary.to_html(index=False, border=0, float_format=lambda value: f"{value:.4g}")
    audit = result.method_audit.to_html(index=False, border=0, escape=True)
    interpretation = result.method_interpretation.to_html(index=False, border=0, escape=True, float_format=lambda value: f"{value:.4g}")
    ranking = result.metric_ranking.to_html(index=False, border=0, escape=True, float_format=lambda value: f"{value:.4g}")
    output_path.write_text(
        f"""<!doctype html>
<html><head><meta charset="utf-8"><title>CoSplice method benchmark</title>
<style>body{{font-family:Arial,sans-serif;margin:36px;color:#17212b;max-width:1400px}}table{{border-collapse:collapse;width:100%;margin:12px 0 24px}}td,th{{border:1px solid #cbd5df;padding:7px;text-align:left}}th{{background:#e9f1f2}}section{{background:#f5f8f8;border-left:4px solid #087e8b;padding:10px 16px;margin:14px 0}}code{{background:#edf2f4;padding:2px 4px}}</style></head>
<body><h1>Katana-CoSplice head-to-head benchmark</h1>
<section><p>All methods used the same <strong>{len(result.universe):,}</strong> stable-AS events and the same samples. The event universe was selected before method fitting using the cached Stable-AS <code>stability_score</code>.</p>
<p>This is an internal same-cohort benchmark. It can compare computational stability, module separation and descriptive condition structure, but it cannot prove external biological superiority without an independent cohort.</p></section>
<h2>Method summary</h2>{summary}
<h2>Strengths, limitations and descriptive ranking</h2>{interpretation}
<p>The composite score is the equal-weight mean of normalized silhouette, neutral graph modularity, within/between separation, bootstrap ARI and inverse runtime. It is a navigation aid, not a statistical test.</p>
<h2>Per-metric ranking</h2>{ranking}
<h2>Method availability and implementation audit</h2>{audit}
<h2>Interpretation</h2><ul>
<li>Bootstrap ARI measures partition recovery after condition-stratified sample resampling.</li>
<li>Neutral graph modularity and cosine silhouette are calculated independently of each method's own optimization score.</li>
<li>Condition-pattern R2 is descriptive in-sample separation, not held-out validation or causality.</li>
<li>WGCNA is reported only when the exact R package runs; no Python approximation is labeled WGCNA.</li>
</ul></body></html>""",
        encoding="utf-8",
    )
    return output_path
