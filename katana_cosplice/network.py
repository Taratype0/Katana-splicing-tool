from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import networkx as nx
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.sparse import csr_matrix
from scipy.spatial.distance import squareform
from sklearn.neighbors import NearestNeighbors

from katana_cosplice.config import NetworkConfig
from katana_cosplice.dynamic_tree import dynamic_tree_cut_labels

if TYPE_CHECKING:
    from katana_cosplice.gene_prior import GeneMantelPrior


@dataclass
class NetworkResult:
    assignments: pd.DataFrame
    edges: pd.DataFrame
    program_scores: pd.DataFrame
    diagnostics: pd.DataFrame
    cut_height_sensitivity: pd.DataFrame
    selected_matrix: pd.DataFrame
    linkage_matrix: np.ndarray | None


def residualize_matrix(matrix: pd.DataFrame, samples: pd.DataFrame, covariates: list[str]) -> pd.DataFrame:
    missing = [column for column in covariates if column not in samples.columns]
    if missing:
        raise ValueError(f"Residual covariates are missing from sample metadata: {missing}")
    if not covariates:
        raise ValueError("Residual network requested but residual_covariates is empty")
    if "sample_id" not in samples.columns or not samples["sample_id"].is_unique:
        raise ValueError("Residualization metadata requires a unique sample_id column.")
    missing_samples = sorted(set(matrix.columns).difference(samples["sample_id"]))
    if missing_samples:
        raise ValueError(f"Residualization metadata is missing matrix samples: {missing_samples}")
    ordered = samples.set_index("sample_id").reindex(matrix.columns)
    design = pd.get_dummies(ordered[covariates].astype(str), drop_first=True, dtype=float)
    design.insert(0, "intercept", 1.0)
    x = design.to_numpy(dtype=float)
    y = matrix.T.to_numpy(dtype=float)
    coefficients = np.linalg.pinv(x) @ y
    residuals = y - x @ coefficients
    return pd.DataFrame(residuals.T, index=matrix.index, columns=matrix.columns)


def _robust_standardize(matrix: np.ndarray) -> np.ndarray:
    median = np.median(matrix, axis=1, keepdims=True)
    mad = np.median(np.abs(matrix - median), axis=1, keepdims=True)
    mad = np.where(mad < 1e-9, np.std(matrix, axis=1, keepdims=True), mad)
    mad = np.where(mad < 1e-9, 1.0, mad)
    u = (matrix - median) / (9.0 * mad)
    weights = np.where(np.abs(u) < 1.0, (1.0 - u**2) ** 2, 0.0)
    centered = (matrix - median) * weights
    norm = np.sqrt(np.sum(centered**2, axis=1, keepdims=True))
    return centered / np.where(norm < 1e-12, 1.0, norm)


def correlation_matrix(matrix: pd.DataFrame, method: str) -> np.ndarray:
    values = matrix.to_numpy(dtype=float)
    if method == "bicor":
        standardized = _robust_standardize(values)
        correlation = standardized @ standardized.T
    elif method == "spearman":
        ranks = pd.DataFrame(values.T).rank(axis=0).to_numpy().T
        correlation = np.corrcoef(ranks)
    elif method == "pearson":
        correlation = np.corrcoef(values)
    else:
        raise ValueError(f"Unsupported network correlation method: {method}")
    correlation = np.nan_to_num(correlation, nan=0.0, posinf=0.0, neginf=0.0)
    np.fill_diagonal(correlation, 1.0)
    return np.clip(correlation, -1.0, 1.0)


def signed_soft_adjacency(
    correlation: np.ndarray,
    settings: NetworkConfig,
) -> np.ndarray:
    """Build the WGCNA-style soft adjacency used by the production network."""

    base = (
        (1.0 + correlation) / 2.0
        if settings.signed
        else np.abs(correlation)
    )
    adjacency = np.power(
        np.clip(base, 0.0, 1.0),
        settings.soft_power,
    ).astype(np.float32, copy=False)
    np.fill_diagonal(adjacency, 0.0)
    return adjacency


def topological_overlap_matrix(adjacency: np.ndarray) -> np.ndarray:
    """Return the weighted TOM similarity for a symmetric adjacency matrix."""

    if adjacency.ndim != 2 or adjacency.shape[0] != adjacency.shape[1]:
        raise ValueError("TOM adjacency must be a square matrix.")
    adjacency = np.asarray(adjacency, dtype=np.float32)
    adjacency = np.clip((adjacency + adjacency.T) / 2.0, 0.0, 1.0)
    np.fill_diagonal(adjacency, 0.0)
    connectivity = adjacency.sum(axis=1, dtype=np.float64)
    shared = adjacency @ adjacency
    denominator = (
        np.minimum.outer(connectivity, connectivity).astype(np.float32)
        + 1.0
        - adjacency
    )
    tom = np.divide(
        shared + adjacency,
        denominator,
        out=np.zeros_like(shared, dtype=np.float32),
        where=denominator > 1e-12,
    )
    tom = np.clip((tom + tom.T) / 2.0, 0.0, 1.0)
    np.fill_diagonal(tom, 1.0)
    return tom


def _tom_hierarchy(
    matrix: pd.DataFrame,
    settings: NetworkConfig,
    *,
    event_reliability: np.ndarray | None = None,
    reliability_power: float = 0.50,
    gene_prior: GeneMantelPrior | None = None,
    gene_prior_weight: float = 0.0,
    precomputed_correlation: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Fit the signed soft-adjacency/TOM backbone without truncating events."""

    if len(matrix) > settings.tom_dense_max_events:
        raise ValueError(
            f"TOM network requires {len(matrix)} events, exceeding "
            f"network.tom_dense_max_events={settings.tom_dense_max_events}. "
            "Raise the explicit audited memory cap or use exact R blockwise WGCNA; "
            "CoSplice will not silently select a top-N subset."
        )
    correlation = (
        correlation_matrix(matrix, settings.method)
        if precomputed_correlation is None
        else np.asarray(precomputed_correlation, dtype=float)
    )
    if correlation.shape != (len(matrix), len(matrix)):
        raise ValueError("precomputed correlation must match the selected event universe")
    correlation = np.clip((correlation + correlation.T) / 2.0, -1.0, 1.0)
    np.fill_diagonal(correlation, 1.0)
    adjacency = signed_soft_adjacency(correlation, settings)
    if event_reliability is not None:
        reliability_factor = np.power(
            np.outer(event_reliability, event_reliability),
            reliability_power / 2.0,
        ).astype(np.float32, copy=False)
        adjacency *= reliability_factor
    base_as_adjacency = adjacency.copy()
    upper_rows, upper_cols = np.triu_indices(len(matrix), k=1)
    reinforced, gene_prior_values = _reinforce_with_gene_prior(
        adjacency[upper_rows, upper_cols],
        matrix.index,
        upper_rows,
        upper_cols,
        gene_prior,
        gene_prior_weight,
    )
    if len(reinforced):
        adjacency[upper_rows, upper_cols] = reinforced
        adjacency[upper_cols, upper_rows] = reinforced
    np.fill_diagonal(adjacency, 0.0)
    tom = topological_overlap_matrix(adjacency)
    distance = np.clip(1.0 - tom, 0.0, 1.0)
    np.fill_diagonal(distance, 0.0)
    hierarchy = linkage(squareform(distance, checks=False), method="average")
    return (
        correlation,
        adjacency,
        tom,
        hierarchy,
        base_as_adjacency,
        gene_prior_values,
    )


def discover_tom_label_grid(
    matrix: pd.DataFrame,
    settings: NetworkConfig,
    cut_heights: list[float],
    *,
    gene_prior: GeneMantelPrior | None = None,
    gene_prior_weight: float = 0.0,
) -> dict[float, np.ndarray]:
    """Cut one TOM hierarchy at several predeclared heights."""

    (
        _correlation,
        _adjacency,
        _tom,
        hierarchy,
        _base_as_adjacency,
        _gene_prior_values,
    ) = (
        _tom_hierarchy(
            matrix,
            settings,
            gene_prior=gene_prior,
            gene_prior_weight=gene_prior_weight,
        )
    )
    output: dict[float, np.ndarray] = {}
    for height in sorted(set(float(value) for value in cut_heights)):
        raw_labels = fcluster(hierarchy, t=height, criterion="distance")
        labels, _ = _module_labels(raw_labels, settings.min_module_size)
        if settings.module_refinement:
            labels, _ = _refine_module_labels(matrix, labels, settings)
        output[height] = labels
    return output


def discover_dynamic_label_grid(
    matrix: pd.DataFrame,
    settings: NetworkConfig,
    deep_splits: list[int],
    *,
    event_reliability: np.ndarray | None = None,
    reliability_power: float = 0.50,
    gene_prior: GeneMantelPrior | None = None,
    gene_prior_weight: float = 0.0,
) -> dict[int, np.ndarray]:
    """Apply published Dynamic Tree Cut variants to one exact TOM hierarchy."""

    (
        _correlation,
        _adjacency,
        tom,
        hierarchy,
        _base_as_adjacency,
        _gene_prior_values,
    ) = _tom_hierarchy(
        matrix,
        settings,
        event_reliability=event_reliability,
        reliability_power=reliability_power,
        gene_prior=gene_prior,
        gene_prior_weight=gene_prior_weight,
    )
    distance = np.clip(1.0 - tom, 0.0, 1.0)
    np.fill_diagonal(distance, 0.0)
    output: dict[int, np.ndarray] = {}
    for deep_split in sorted(set(int(value) for value in deep_splits)):
        raw_labels = dynamic_tree_cut_labels(
            hierarchy,
            distance,
            min_module_size=settings.min_module_size,
            deep_split=deep_split,
            cut_height=settings.dynamic_cut_height,
            pam_respects_dendro=settings.dynamic_pam_respects_dendro,
        )
        labels, _ = _module_labels(raw_labels, settings.min_module_size)
        if settings.module_refinement:
            labels, _ = _refine_module_labels(matrix, labels, settings)
        output[deep_split] = labels
    return output


def _standardized_profiles(matrix: pd.DataFrame, method: str) -> np.ndarray:
    values = matrix.to_numpy(dtype=float)
    if method == "bicor":
        return _robust_standardize(values)
    if method == "spearman":
        values = pd.DataFrame(values.T).rank(axis=0).to_numpy().T
    elif method != "pearson":
        raise ValueError(f"Unsupported network correlation method: {method}")
    values = values - values.mean(axis=1, keepdims=True)
    norm = np.sqrt(np.sum(values**2, axis=1, keepdims=True))
    return values / np.where(norm < 1e-12, 1.0, norm)


def _sparse_knn_edges(
    matrix: pd.DataFrame,
    settings: NetworkConfig,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    profiles = _standardized_profiles(matrix, settings.method)
    n_events = len(matrix)
    n_neighbors = min(n_events, settings.knn_neighbors + 1)
    model = NearestNeighbors(n_neighbors=n_neighbors, metric="cosine", algorithm="brute", n_jobs=-1)
    model.fit(profiles)
    distances, indices = model.kneighbors(profiles, return_distance=True)
    negative_queries: tuple[np.ndarray, np.ndarray] | None = None
    if not settings.signed:
        negative_queries = model.kneighbors(-profiles, return_distance=True)
    rows, cols, correlations, adjacencies = _raw_knn_edges_from_queries(
        distances,
        indices,
        settings,
        negative_queries,
    )
    overlaps = np.zeros(len(rows), dtype=float)
    if settings.edge_smoothing == "shared_neighbor" and len(rows):
        rows, cols, correlations, adjacencies, overlaps = _shared_neighbor_smoothing(
            profiles,
            rows,
            cols,
            correlations,
            adjacencies,
            settings,
        )
    return rows, cols, correlations, adjacencies, overlaps


def _raw_knn_edges_from_queries(
    distances: np.ndarray,
    indices: np.ndarray,
    settings: NetworkConfig,
    negative_queries: tuple[np.ndarray, np.ndarray] | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Build raw sparse edges from a reusable nearest-neighbor query."""

    edge_map: dict[tuple[int, int], tuple[float, float]] = {}

    def collect(query_distances: np.ndarray, query_indices: np.ndarray, negate: bool = False) -> None:
        for source in range(len(query_distances)):
            for distance, target in zip(query_distances[source], query_indices[source], strict=True):
                target = int(target)
                if source == target:
                    continue
                correlation = 1.0 - float(distance)
                if negate:
                    correlation = -correlation
                base = (1.0 + correlation) / 2.0 if settings.signed else abs(correlation)
                adjacency = float(np.clip(base, 0.0, 1.0) ** settings.soft_power)
                if adjacency < settings.edge_threshold:
                    continue
                pair = (source, target) if source < target else (target, source)
                previous = edge_map.get(pair)
                if previous is None or adjacency > previous[1]:
                    edge_map[pair] = (correlation, adjacency)

    collect(distances, indices)
    if negative_queries is not None:
        negative_distances, negative_indices = negative_queries
        collect(negative_distances, negative_indices, negate=True)

    if not edge_map:
        empty_int = np.asarray([], dtype=int)
        empty_float = np.asarray([], dtype=float)
        return empty_int, empty_int, empty_float, empty_float
    pairs = np.asarray(list(edge_map), dtype=int)
    values = np.asarray(list(edge_map.values()), dtype=float)
    rows, cols = pairs[:, 0], pairs[:, 1]
    correlations, adjacencies = values[:, 0], values[:, 1]
    return rows, cols, correlations, adjacencies


def _multi_k_consensus_edges(
    matrix: pd.DataFrame,
    settings: NetworkConfig,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Retain locally reproducible edges across multiple k-neighborhood scales."""
    neighbor_values = sorted(
        {
            min(int(value), max(2, len(matrix) - 1))
            for value in settings.consensus_knn_neighbors
        }
    )
    profiles = _standardized_profiles(matrix, settings.method)
    max_neighbors = max(neighbor_values)
    model = NearestNeighbors(
        n_neighbors=min(len(matrix), max_neighbors + 1),
        metric="cosine",
        algorithm="brute",
        n_jobs=-1,
    ).fit(profiles)
    all_distances, all_indices = model.kneighbors(profiles, return_distance=True)
    all_negative: tuple[np.ndarray, np.ndarray] | None = None
    if not settings.signed:
        all_negative = model.kneighbors(-profiles, return_distance=True)
    edge_values: dict[tuple[int, int], list[tuple[float, float]]] = {}
    for neighbors in neighbor_values:
        query_width = min(all_distances.shape[1], neighbors + 1)
        negative = (
            (all_negative[0][:, :query_width], all_negative[1][:, :query_width])
            if all_negative is not None
            else None
        )
        variant = replace(
            settings,
            knn_neighbors=neighbors,
            multi_k_consensus=False,
            edge_smoothing="none",
        )
        rows, cols, correlations, adjacencies = _raw_knn_edges_from_queries(
            all_distances[:, :query_width],
            all_indices[:, :query_width],
            variant,
            negative,
        )
        for row, col, correlation, adjacency in zip(
            rows,
            cols,
            correlations,
            adjacencies,
            strict=True,
        ):
            edge_values.setdefault((int(row), int(col)), []).append(
                (float(correlation), float(adjacency))
            )
    retained: list[tuple[int, int, float, float, float, float]] = []
    n_scales = len(neighbor_values)
    for (row, col), values in edge_values.items():
        support = len(values) / n_scales
        if support + 1e-12 < settings.consensus_edge_min_fraction:
            continue
        array = np.asarray(values, dtype=float)
        adjacency = float(array[:, 1].mean() * (0.5 + 0.5 * support))
        if adjacency < settings.edge_threshold:
            continue
        retained.append(
            (
                row,
                col,
                float(array[:, 0].mean()),
                adjacency,
                0.0,
                support,
            )
        )
    if not retained:
        empty_int = np.asarray([], dtype=int)
        empty_float = np.asarray([], dtype=float)
        return empty_int, empty_int, empty_float, empty_float, empty_float, empty_float
    values = np.asarray(retained, dtype=float)
    rows = values[:, 0].astype(int)
    cols = values[:, 1].astype(int)
    correlations = values[:, 2]
    adjacencies = values[:, 3]
    overlaps = values[:, 4]
    support_lookup = {
        (int(row), int(col)): float(support)
        for row, col, support in zip(rows, cols, values[:, 5], strict=True)
    }
    if settings.edge_smoothing == "shared_neighbor" and len(rows):
        rows, cols, correlations, adjacencies, overlaps = _shared_neighbor_smoothing(
            profiles,
            rows,
            cols,
            correlations,
            adjacencies,
            settings,
            allow_new_edges=False,
        )
    support = np.asarray(
        [
            support_lookup[(min(int(row), int(col)), max(int(row), int(col)))]
            for row, col in zip(rows, cols, strict=True)
        ],
        dtype=float,
    )
    return (
        rows,
        cols,
        correlations,
        adjacencies,
        overlaps,
        support,
    )


def _network_edges(
    matrix: pd.DataFrame,
    settings: NetworkConfig,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if settings.multi_k_consensus:
        return _multi_k_consensus_edges(matrix, settings)
    rows, cols, correlations, adjacencies, overlaps = _sparse_knn_edges(matrix, settings)
    return rows, cols, correlations, adjacencies, overlaps, np.ones(len(rows), dtype=float)


def _shared_neighbor_smoothing(
    profiles: np.ndarray,
    rows: np.ndarray,
    cols: np.ndarray,
    correlations: np.ndarray,
    adjacencies: np.ndarray,
    settings: NetworkConfig,
    *,
    allow_new_edges: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Add sparse weighted topological-overlap edges without materializing dense TOM."""
    n_events = len(profiles)
    symmetric_rows = np.concatenate([rows, cols])
    symmetric_cols = np.concatenate([cols, rows])
    symmetric_values = np.concatenate([adjacencies, adjacencies])
    adjacency = csr_matrix(
        (symmetric_values, (symmetric_rows, symmetric_cols)),
        shape=(n_events, n_events),
    )
    connectivity = np.asarray(adjacency.sum(axis=1)).ravel()
    shared_product = (adjacency @ adjacency).tocsr()
    direct_dense = np.zeros(n_events, dtype=float)
    edge_map: dict[tuple[int, int], tuple[float, float]] = {
        (int(row), int(col)): (float(weight), 0.0)
        for row, col, weight in zip(rows, cols, adjacencies, strict=True)
    }

    for source in range(n_events):
        direct_start, direct_end = adjacency.indptr[source : source + 2]
        direct_columns = adjacency.indices[direct_start:direct_end]
        direct_dense[direct_columns] = adjacency.data[direct_start:direct_end]
        shared_start, shared_end = shared_product.indptr[source : source + 2]
        candidates = shared_product.indices[shared_start:shared_end]
        shared = shared_product.data[shared_start:shared_end]
        valid = candidates != source
        candidates = candidates[valid]
        shared = shared[valid]
        if len(candidates):
            direct = direct_dense[candidates]
            denominator = np.minimum(connectivity[source], connectivity[candidates]) + 1.0 - direct
            overlap = np.divide(
                shared + direct,
                denominator,
                out=np.zeros_like(shared, dtype=float),
                where=denominator > 1e-12,
            )
            smoothed = np.maximum(direct, settings.shared_neighbor_weight * overlap)
            retained = np.flatnonzero(smoothed >= settings.edge_threshold)
            if len(retained) > settings.shared_neighbor_top_k:
                top = np.argpartition(smoothed[retained], -settings.shared_neighbor_top_k)[
                    -settings.shared_neighbor_top_k :
                ]
                retained = retained[top]
            for index in retained:
                target = int(candidates[index])
                pair = (source, target) if source < target else (target, source)
                previous = edge_map.get(pair)
                if previous is None and not allow_new_edges:
                    continue
                candidate_value = (float(smoothed[index]), float(overlap[index]))
                if previous is None or candidate_value[0] > previous[0]:
                    edge_map[pair] = candidate_value
                elif candidate_value[1] > previous[1]:
                    edge_map[pair] = (previous[0], candidate_value[1])
        direct_dense[direct_columns] = 0.0

    pairs = np.asarray(list(edge_map), dtype=int)
    weights = np.asarray([value[0] for value in edge_map.values()], dtype=float)
    overlaps = np.asarray([value[1] for value in edge_map.values()], dtype=float)
    smoothed_correlations = np.sum(profiles[pairs[:, 0]] * profiles[pairs[:, 1]], axis=1)
    if not settings.signed:
        signs = np.sign(smoothed_correlations)
        smoothed_correlations = np.abs(smoothed_correlations) * signs
    return pairs[:, 0], pairs[:, 1], smoothed_correlations, weights, overlaps


def _sparse_edge_table(
    rows: np.ndarray,
    cols: np.ndarray,
    correlations: np.ndarray,
    adjacencies: np.ndarray,
    overlaps: np.ndarray,
    consensus_support: np.ndarray,
    graph_bootstrap_support: np.ndarray,
    event_keys: list[str],
    max_edges: int,
    base_as_adjacencies: np.ndarray | None = None,
    gene_prior_values: np.ndarray | None = None,
) -> pd.DataFrame:
    if len(rows) == 0:
        return pd.DataFrame(
            columns=[
                "event_a",
                "event_b",
                "correlation",
                "adjacency",
                "shared_neighbor_overlap",
                "knn_consensus_support",
                "graph_bootstrap_support",
            ]
        )
    selected = np.arange(len(rows))
    if len(selected) > max_edges:
        selected = np.argpartition(adjacencies, -max_edges)[-max_edges:]
    keys = np.asarray(event_keys, dtype=object)
    result = pd.DataFrame(
        {
            "event_a": keys[rows[selected]],
            "event_b": keys[cols[selected]],
            "correlation": correlations[selected],
            "adjacency": adjacencies[selected],
            "shared_neighbor_overlap": overlaps[selected],
            "knn_consensus_support": consensus_support[selected],
            "graph_bootstrap_support": graph_bootstrap_support[selected],
        }
    )
    if base_as_adjacencies is not None:
        result["base_AS_adjacency"] = base_as_adjacencies[selected]
    if gene_prior_values is not None:
        result["gene_mantel_prior"] = gene_prior_values[selected]
        result["gene_prior_adjacency_boost"] = (
            result["adjacency"] - result.get("base_AS_adjacency", result["adjacency"])
        )
    return result.sort_values("adjacency", ascending=False).reset_index(drop=True)


def _reinforce_with_gene_prior(
    adjacencies: np.ndarray,
    matrix_index: pd.Index,
    rows: np.ndarray,
    cols: np.ndarray,
    gene_prior: GeneMantelPrior | None,
    weight: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Apply a bounded prior boost without creating a non-AS edge.

    The base AS adjacency remains unchanged when no gene-level support exists.
    A fully supported prior can add at most ``weight * a * (1-a)`` to an
    existing AS adjacency ``a``.
    """
    base = np.asarray(adjacencies, dtype=float).copy()
    values = np.zeros(len(base), dtype=float)
    if gene_prior is None or weight <= 0 or not len(base):
        return base, values
    from katana_cosplice.gene_prior import edge_gene_prior

    values = edge_gene_prior(gene_prior, matrix_index, rows, cols)
    reinforced = base + float(weight) * base * (1.0 - base) * values
    return np.clip(reinforced, 0.0, 1.0), values


def _module_labels(raw_labels: np.ndarray, min_size: int) -> tuple[np.ndarray, dict[int, str]]:
    counts = pd.Series(raw_labels).value_counts()
    retained = [
        int(label)
        for label, count in counts.items()
        if int(label) != 0 and count >= min_size
    ]
    retained.sort(key=lambda label: (-int(counts[label]), label))
    names = {label: f"M{index:02d}" for index, label in enumerate(retained, start=1)}
    labels = np.asarray([names.get(int(label), "grey") for label in raw_labels], dtype=object)
    return labels, names


def _program_scores(matrix: pd.DataFrame, module_labels: np.ndarray) -> pd.DataFrame:
    output: dict[str, np.ndarray] = {}
    for module in sorted(set(module_labels).difference({"grey"})):
        subset = matrix.loc[module_labels == module]
        values = subset.to_numpy(dtype=float)
        values = values - values.mean(axis=1, keepdims=True)
        scale = values.std(axis=1, ddof=1, keepdims=True)
        values = values / np.where(scale < 1e-12, 1.0, scale)
        _, _, vt = np.linalg.svd(values, full_matrices=False)
        score = vt[0]
        mean_pattern = values.mean(axis=0)
        if np.corrcoef(score, mean_pattern)[0, 1] < 0:
            score = -score
        output[module] = score
    frame = pd.DataFrame(output, index=matrix.columns)
    frame.index.name = "sample_id"
    return frame


def _canonicalize_labels(labels: np.ndarray, min_size: int) -> np.ndarray:
    """Rename retained modules by size while preserving explicit grey assignments."""
    values = pd.Series(np.asarray(labels, dtype=object), dtype=object)
    counts = values.loc[values.ne("grey")].value_counts()
    retained = [str(label) for label, count in counts.items() if int(count) >= min_size]
    retained.sort(key=lambda label: (-int(counts[label]), label))
    names = {label: f"M{index:02d}" for index, label in enumerate(retained, start=1)}
    return np.asarray([names.get(str(label), "grey") for label in values], dtype=object)


def _merge_correlated_modules(
    matrix: pd.DataFrame,
    labels: np.ndarray,
    settings: NetworkConfig,
) -> np.ndarray:
    modules = sorted(set(labels).difference({"grey"}))
    if len(modules) < 2:
        return labels.copy()
    scores = _program_scores(matrix, labels)
    correlation = scores.corr(method="pearson").fillna(0.0)
    graph = nx.Graph()
    graph.add_nodes_from(modules)
    for index, module_a in enumerate(modules):
        for module_b in modules[index + 1 :]:
            value = float(correlation.loc[module_a, module_b])
            similarity = value if settings.signed else abs(value)
            if similarity >= settings.module_merge_correlation:
                graph.add_edge(module_a, module_b)
    mapping: dict[str, str] = {}
    for component_index, component in enumerate(nx.connected_components(graph), start=1):
        merged_name = f"merge_{component_index:03d}"
        mapping.update({str(module): merged_name for module in component})
    merged = np.asarray([mapping.get(str(label), "grey") for label in labels], dtype=object)
    return _canonicalize_labels(merged, settings.min_module_size)


def _membership_statistics(
    matrix: pd.DataFrame,
    labels: np.ndarray,
    settings: NetworkConfig,
) -> pd.DataFrame:
    modules = sorted(set(labels).difference({"grey"}))
    output = pd.DataFrame(index=matrix.index)
    output.attrs["n_samples"] = matrix.shape[1]
    output["assigned_module"] = labels
    output["module_membership"] = np.nan
    output["best_module"] = "grey"
    output["best_module_membership"] = np.nan
    output["second_best_module_membership"] = np.nan
    output["module_membership_margin"] = np.nan
    if not modules:
        return output
    scores = _program_scores(matrix, labels)
    event_profiles = _standardized_profiles(matrix, settings.method)
    score_profiles = _standardized_profiles(scores.T, settings.method)
    correlations = event_profiles @ score_profiles.T
    if not settings.signed:
        correlations = np.abs(correlations)
    module_index = {module: index for index, module in enumerate(modules)}
    best_indices = np.argmax(correlations, axis=1)
    best_values = correlations[np.arange(len(matrix)), best_indices]
    if len(modules) > 1:
        second_values = np.partition(correlations, -2, axis=1)[:, -2]
    else:
        second_values = np.full(len(matrix), -1.0)
    assigned_values = np.asarray(
        [correlations[row, module_index[label]] if label in module_index else np.nan for row, label in enumerate(labels)],
        dtype=float,
    )
    output["module_membership"] = assigned_values
    output["best_module"] = np.asarray(modules, dtype=object)[best_indices]
    output["best_module_membership"] = best_values
    output["second_best_module_membership"] = second_values
    output["module_membership_margin"] = best_values - second_values
    return output


def _refine_module_labels(
    matrix: pd.DataFrame,
    initial_labels: np.ndarray,
    settings: NetworkConfig,
) -> tuple[np.ndarray, pd.DataFrame]:
    """Merge redundant programs and reject weak/ambiguous event memberships."""
    initial = np.asarray(initial_labels, dtype=object)
    labels = initial.copy()
    for _ in range(settings.refinement_iterations):
        labels = _merge_correlated_modules(matrix, labels, settings)
        statistics = _membership_statistics(matrix, labels, settings)
        best = statistics["best_module"].astype(str).to_numpy()
        membership = statistics["best_module_membership"].to_numpy(dtype=float)
        margin = statistics["module_membership_margin"].to_numpy(dtype=float)
        cutoffs = _adaptive_membership_cutoffs(statistics, settings)
        if settings.adaptive_membership:
            current = statistics["assigned_module"].astype(str).to_numpy()
            assigned_membership = statistics["module_membership"].to_numpy(dtype=float)
            retain_current = (
                (current != "grey")
                & (assigned_membership >= cutoffs)
                & (margin >= settings.min_module_membership_margin)
            )
            fixed_reassignment = (
                (membership >= settings.min_module_membership)
                & (margin >= settings.min_module_membership_margin)
            )
            labels = np.where(
                retain_current,
                current,
                np.where(fixed_reassignment, best, "grey"),
            ).astype(object)
        else:
            accepted = (membership >= settings.min_module_membership) & (
                margin >= settings.min_module_membership_margin
            )
            labels = np.where(accepted, best, "grey").astype(object)
        labels = _canonicalize_labels(labels, settings.min_module_size)

    final_statistics = _membership_statistics(matrix, labels, settings)
    final_statistics.insert(0, "module_before_refinement", initial)
    final_statistics["module_after_refinement"] = labels
    membership = final_statistics["best_module_membership"].to_numpy(dtype=float)
    margin = final_statistics["module_membership_margin"].to_numpy(dtype=float)
    cutoffs = _adaptive_membership_cutoffs(final_statistics, settings)
    final_statistics["membership_cutoff_used"] = cutoffs
    action = np.where(labels == "grey", "grey_small_initial_community", "retained").astype(object)
    action[(labels == "grey") & (membership < cutoffs)] = "grey_low_membership"
    action[(labels == "grey") & (membership >= cutoffs) & (margin < settings.min_module_membership_margin)] = "grey_ambiguous_membership"
    action[(labels == "grey") & ~np.isfinite(membership)] = "grey_no_retained_module"
    initial_to_final: dict[str, str] = {}
    for initial_module in sorted(set(initial).difference({"grey"})):
        destinations = pd.Series(labels[initial == initial_module])
        destinations = destinations.loc[destinations.ne("grey")]
        if not destinations.empty:
            initial_to_final[str(initial_module)] = str(destinations.value_counts().index[0])
    final_sources: dict[str, set[str]] = {}
    for initial_module, final_module in initial_to_final.items():
        final_sources.setdefault(final_module, set()).add(initial_module)
    for index, (before, after) in enumerate(zip(initial, labels, strict=True)):
        before = str(before)
        after = str(after)
        if after == "grey":
            continue
        if before == "grey":
            action[index] = "assigned_from_grey"
        elif initial_to_final.get(before) != after:
            action[index] = "reassigned_module"
        elif len(final_sources.get(after, set())) > 1:
            action[index] = "merged_module"
    final_statistics["module_refinement_action"] = action
    return labels, final_statistics


def _adaptive_membership_cutoffs(
    statistics: pd.DataFrame,
    settings: NetworkConfig,
) -> np.ndarray:
    """Estimate robust module-specific kME cutoffs within declared safety bounds."""
    if not settings.adaptive_membership:
        return np.full(len(statistics), settings.min_module_membership, dtype=float)
    assigned_modules = statistics["assigned_module"].astype(str)
    memberships = pd.to_numeric(statistics["module_membership"], errors="coerce")
    sampling_noise_floor = 1.0 / np.sqrt(max(1, int(statistics.attrs.get("n_samples", 1))))
    effective_floor = min(
        settings.min_module_membership,
        max(settings.adaptive_membership_floor, sampling_noise_floor),
    )
    module_cutoffs: dict[str, float] = {}
    for module in sorted(set(assigned_modules).difference({"grey"})):
        # Calibrate from the current module members only; previously grey events
        # must not lower the module's acceptance threshold.
        values = memberships.loc[assigned_modules.eq(module)].dropna().to_numpy(dtype=float)
        if len(values) < settings.min_module_size:
            module_cutoffs[module] = settings.min_module_membership
            continue
        median = float(np.median(values))
        mad = float(np.median(np.abs(values - median)))
        robust_lower_bound = median - settings.adaptive_membership_mad_scale * 1.4826 * mad
        module_cutoffs[module] = float(
            np.clip(
                robust_lower_bound,
                effective_floor,
                settings.min_module_membership,
            )
        )
    # Adaptive relaxation applies only to current module members. Events that
    # are already grey must clear the original fixed cutoff before assignment,
    # preventing weak bridge events from merging otherwise distinct programs.
    return np.asarray(
        [
            module_cutoffs.get(assigned, settings.min_module_membership)
            if assigned != "grey"
            else settings.min_module_membership
            for assigned in assigned_modules
        ],
        dtype=float,
    )


def _edge_table(
    correlation: np.ndarray,
    adjacency: np.ndarray,
    event_keys: list[str],
    threshold: float,
    max_edges: int = 100_000,
    *,
    base_as_adjacency: np.ndarray | None = None,
    gene_prior_matrix: np.ndarray | None = None,
) -> pd.DataFrame:
    upper = np.triu_indices_from(adjacency, k=1)
    scores = adjacency[upper]
    selected = np.flatnonzero(scores >= threshold)
    if len(selected) > max_edges:
        top = np.argpartition(scores[selected], -max_edges)[-max_edges:]
        selected = selected[top]
    rows = upper[0][selected]
    cols = upper[1][selected]
    result = pd.DataFrame(
        {
            "event_a": np.asarray(event_keys, dtype=object)[rows],
            "event_b": np.asarray(event_keys, dtype=object)[cols],
            "correlation": correlation[rows, cols],
            "adjacency": adjacency[rows, cols],
        }
    )
    if base_as_adjacency is not None:
        result["base_AS_adjacency"] = base_as_adjacency[rows, cols]
    if gene_prior_matrix is not None:
        result["gene_mantel_prior"] = gene_prior_matrix[rows, cols]
        result["gene_prior_adjacency_boost"] = (
            result["adjacency"] - result.get("base_AS_adjacency", result["adjacency"])
        )
    return result.sort_values("adjacency", ascending=False).reset_index(drop=True)


def _edge_bootstrap_support(
    matrix: pd.DataFrame,
    edges: pd.DataFrame,
    settings: NetworkConfig,
    samples: pd.DataFrame | None,
) -> np.ndarray:
    if edges.empty or settings.bootstrap_iterations <= 0:
        return np.full(len(edges), np.nan)
    event_index = {event_key: index for index, event_key in enumerate(matrix.index)}
    rows = edges["event_a"].map(event_index).to_numpy(dtype=int)
    cols = edges["event_b"].map(event_index).to_numpy(dtype=int)
    return _edge_bootstrap_support_arrays(
        matrix,
        rows,
        cols,
        settings,
        samples,
        settings.bootstrap_iterations,
    )


def _edge_bootstrap_support_arrays(
    matrix: pd.DataFrame,
    rows: np.ndarray,
    cols: np.ndarray,
    settings: NetworkConfig,
    samples: pd.DataFrame | None,
    iterations: int,
) -> np.ndarray:
    if len(rows) == 0 or iterations <= 0:
        return np.full(len(rows), np.nan)
    rng = np.random.default_rng(settings.random_seed)
    from katana_cosplice.resampling import bootstrap_plan
    plan = bootstrap_plan(samples, matrix.columns)
    plan.require_resampleable()
    support = np.zeros(len(rows), dtype=float)
    values = matrix.to_numpy(dtype=float)
    for _ in range(iterations):
        sampled = plan.draw(rng)
        boot = values[:, sampled]
        if settings.method == "bicor":
            standardized = _robust_standardize(boot)
        elif settings.method == "spearman":
            ranks = pd.DataFrame(boot.T).rank(axis=0).to_numpy().T
            ranks -= ranks.mean(axis=1, keepdims=True)
            norm = np.sqrt(np.sum(ranks**2, axis=1, keepdims=True))
            standardized = ranks / np.where(norm < 1e-12, 1.0, norm)
        else:
            centered = boot - boot.mean(axis=1, keepdims=True)
            norm = np.sqrt(np.sum(centered**2, axis=1, keepdims=True))
            standardized = centered / np.where(norm < 1e-12, 1.0, norm)
        correlations = np.sum(standardized[rows] * standardized[cols], axis=1)
        base = (1.0 + correlations) / 2.0 if settings.signed else np.abs(correlations)
        boot_adjacency = np.power(np.clip(base, 0.0, 1.0), settings.soft_power)
        support += boot_adjacency >= settings.edge_threshold
    return support / float(iterations)


def _bootstrap_weighted_edges(
    matrix: pd.DataFrame,
    rows: np.ndarray,
    cols: np.ndarray,
    adjacencies: np.ndarray,
    settings: NetworkConfig,
    samples: pd.DataFrame | None,
) -> tuple[np.ndarray, np.ndarray]:
    iterations = _effective_consensus_bootstrap_iterations(len(matrix), settings)
    if not len(rows) or iterations <= 0 or settings.consensus_bootstrap_weight <= 0:
        return adjacencies, np.full(len(rows), np.nan)
    support = _edge_bootstrap_support_arrays(
        matrix,
        rows,
        cols,
        settings,
        samples,
        iterations,
    )
    weight = settings.consensus_bootstrap_weight
    adjusted = adjacencies * ((1.0 - weight) + weight * support)
    return adjusted, support


def _effective_consensus_bootstrap_iterations(
    n_events: int,
    settings: NetworkConfig,
) -> int:
    """Return the declared graph-resampling effort without a size-based shortcut."""
    del n_events
    return max(0, settings.consensus_bootstrap_iterations)


def _cut_height_sensitivity(linkage_matrix: np.ndarray, settings: NetworkConfig) -> pd.DataFrame:
    if settings.engine == "tom_hierarchical":
        heights = sorted(
            set([0.70, 0.75, settings.module_cut_height, 0.85, 0.90])
        )
    else:
        heights = sorted(
            set(
                [
                    0.35,
                    0.45,
                    0.55,
                    0.65,
                    settings.module_cut_height,
                    0.75,
                    0.85,
                ]
            )
        )
    rows: list[dict[str, object]] = []
    n_events = linkage_matrix.shape[0] + 1
    for height in heights:
        raw = fcluster(linkage_matrix, t=height, criterion="distance")
        labels, _ = _module_labels(raw, settings.min_module_size)
        counts = pd.Series(labels).value_counts()
        retained = counts.drop(labels="grey", errors="ignore")
        largest = int(retained.max()) if not retained.empty else 0
        rows.append(
            {
                "module_cut_height": height,
                "n_modules": int(len(retained)),
                "n_grey_events": int(counts.get("grey", 0)),
                "largest_module_events": largest,
                "largest_module_fraction": largest / n_events,
                "is_active_setting": bool(np.isclose(height, settings.module_cut_height)),
            }
        )
    return pd.DataFrame(rows)


def _dynamic_tree_sensitivity(
    linkage_matrix: np.ndarray,
    tom: np.ndarray,
    settings: NetworkConfig,
) -> pd.DataFrame:
    distance = np.clip(1.0 - tom, 0.0, 1.0)
    np.fill_diagonal(distance, 0.0)
    rows: list[dict[str, object]] = []
    n_events = len(distance)
    for deep_split in sorted(set(settings.dynamic_deep_splits)):
        raw = dynamic_tree_cut_labels(
            linkage_matrix,
            distance,
            min_module_size=settings.min_module_size,
            deep_split=deep_split,
            cut_height=settings.dynamic_cut_height,
            pam_respects_dendro=settings.dynamic_pam_respects_dendro,
        )
        labels, _ = _module_labels(raw, settings.min_module_size)
        counts = pd.Series(labels).value_counts()
        retained = counts.drop(labels="grey", errors="ignore")
        largest = int(retained.max()) if not retained.empty else 0
        rows.append(
            {
                "dynamic_deep_split": deep_split,
                "dynamic_cut_height": settings.dynamic_cut_height,
                "n_modules": int(len(retained)),
                "n_grey_events": int(counts.get("grey", 0)),
                "largest_module_events": largest,
                "largest_module_fraction": largest / n_events,
                "is_active_setting": deep_split == settings.dynamic_deep_split,
            }
        )
    return pd.DataFrame(rows)


def _graph_communities(
    n_events: int,
    rows: np.ndarray,
    cols: np.ndarray,
    adjacencies: np.ndarray,
    settings: NetworkConfig,
    resolution: float | None = None,
) -> tuple[np.ndarray, nx.Graph]:
    graph = nx.Graph()
    graph.add_nodes_from(range(n_events))
    graph.add_weighted_edges_from(
        (int(row), int(col), float(weight))
        for row, col, weight in zip(rows, cols, adjacencies, strict=True)
    )
    if graph.number_of_edges() == 0:
        return np.arange(1, n_events + 1, dtype=int), graph
    active_resolution = float(resolution if resolution is not None else settings.community_resolution)

    def run_louvain(active_graph: nx.Graph, seed: int) -> list[set[int]]:
        return nx.community.louvain_communities(
            active_graph,
            weight="weight",
            resolution=active_resolution,
            seed=seed,
        )

    if settings.community_consensus_runs > 1 and graph.number_of_edges():
        edge_pairs = list(graph.edges())
        agreement = np.zeros(len(edge_pairs), dtype=float)
        for run in range(settings.community_consensus_runs):
            partition = run_louvain(graph, settings.random_seed + run)
            labels_for_run = np.zeros(n_events, dtype=int)
            for label, community in enumerate(partition, start=1):
                labels_for_run[np.fromiter(community, dtype=int)] = label
            agreement += np.asarray(
                [labels_for_run[source] == labels_for_run[target] for source, target in edge_pairs],
                dtype=float,
            )
        agreement /= settings.community_consensus_runs
        consensus_graph = nx.Graph()
        consensus_graph.add_nodes_from(graph.nodes())
        consensus_weight = settings.community_consensus_weight
        for (source, target), coassignment in zip(edge_pairs, agreement, strict=True):
            base_weight = float(graph[source][target]["weight"])
            adjusted_weight = base_weight * (
                (1.0 - consensus_weight) + consensus_weight * float(coassignment)
            )
            consensus_graph.add_edge(source, target, weight=adjusted_weight)
        graph = consensus_graph
    communities = run_louvain(graph, settings.random_seed + settings.community_consensus_runs)
    raw_labels = np.zeros(n_events, dtype=int)
    for label, community in enumerate(communities, start=1):
        raw_labels[np.fromiter(community, dtype=int)] = label
    return raw_labels, graph


def _community_resolution_sensitivity(
    n_events: int,
    rows: np.ndarray,
    cols: np.ndarray,
    adjacencies: np.ndarray,
    settings: NetworkConfig,
) -> pd.DataFrame:
    resolutions = sorted(set([0.50, 0.75, settings.community_resolution, 1.25, 1.50, 2.00]))
    output: list[dict[str, object]] = []
    # Sensitivity points test the resolution perturbation itself. Repeating the
    # full community-consensus procedure at every diagnostic point multiplies
    # runtime without changing that question; the fitted network above still
    # uses the declared consensus runs.
    sensitivity_settings = replace(settings, community_consensus_runs=1)
    for resolution in resolutions:
        raw_labels, _ = _graph_communities(
            n_events,
            rows,
            cols,
            adjacencies,
            sensitivity_settings,
            resolution,
        )
        labels, _ = _module_labels(raw_labels, settings.min_module_size)
        counts = pd.Series(labels).value_counts()
        retained = counts.drop(labels="grey", errors="ignore")
        largest = int(retained.max()) if not retained.empty else 0
        output.append(
            {
                "community_resolution": resolution,
                "n_modules": int(len(retained)),
                "n_grey_events": int(counts.get("grey", 0)),
                "largest_module_events": largest,
                "largest_module_fraction": largest / n_events,
                "is_active_setting": bool(np.isclose(resolution, settings.community_resolution)),
            }
        )
    return pd.DataFrame(output)


def discover_module_labels(
    matrix: pd.DataFrame,
    settings: NetworkConfig,
    community_resolution: float | None = None,
    event_reliability: pd.Series | None = None,
    reliability_floor: float = 0.10,
    reliability_power: float = 0.50,
    gene_prior: GeneMantelPrior | None = None,
    gene_prior_weight: float = 0.0,
) -> np.ndarray:
    """Discover module labels without program scores, edge bootstrap or sensitivity sweeps."""
    if settings.engine == "tom_hierarchical":
        reliability = (
            pd.to_numeric(event_reliability.reindex(matrix.index), errors="coerce")
            .fillna(0.0)
            .clip(reliability_floor, 1.0)
            .to_numpy(dtype=float)
            if event_reliability is not None
            else None
        )
        if settings.module_detection == "wgcna_blockwise":
            from katana_cosplice.r_wgcna import exact_wgcna_labels

            return exact_wgcna_labels(
                matrix,
                soft_power=settings.soft_power,
                min_module_size=settings.min_module_size,
                merge_cut_height=settings.wgcna_merge_cut_height,
                deep_split=settings.dynamic_deep_split,
                pam_respects_dendro=settings.dynamic_pam_respects_dendro,
                cor_type=settings.method,
                signed=settings.signed,
                rscript=settings.wgcna_rscript,
                event_reliability=event_reliability,
                reliability_power=reliability_power,
            )
        if settings.module_detection == "dynamic_tree_cut":
            return discover_dynamic_label_grid(
                matrix,
                settings,
                [settings.dynamic_deep_split],
                event_reliability=reliability,
                reliability_power=reliability_power,
                gene_prior=gene_prior,
                gene_prior_weight=gene_prior_weight,
            )[int(settings.dynamic_deep_split)]
        return discover_tom_label_grid(
            matrix,
            settings,
            [settings.module_cut_height],
            gene_prior=gene_prior,
            gene_prior_weight=gene_prior_weight,
        )[float(settings.module_cut_height)]
    elif settings.engine == "sparse_knn":
        rows, cols, _correlations, adjacencies, _overlaps, _consensus = _network_edges(
            matrix,
            settings,
        )
        adjacencies, _bootstrap_support = _bootstrap_weighted_edges(
            matrix,
            rows,
            cols,
            adjacencies,
            settings,
            None,
        )
        if len(_bootstrap_support) and np.isfinite(_bootstrap_support).any():
            retained = _bootstrap_support >= settings.consensus_bootstrap_min_support
            rows = rows[retained]
            cols = cols[retained]
            adjacencies = adjacencies[retained]
        adjacencies, _ = _reinforce_with_gene_prior(
            adjacencies,
            matrix.index,
            rows,
            cols,
            gene_prior,
            gene_prior_weight,
        )
        raw_labels, _ = _graph_communities(
            len(matrix),
            rows,
            cols,
            adjacencies,
            settings,
            community_resolution,
        )
    else:
        correlation = correlation_matrix(matrix, settings.method)
        if gene_prior is not None and gene_prior_weight > 0:
            base = (1.0 + correlation) / 2.0 if settings.signed else np.abs(correlation)
            upper_rows, upper_cols = np.triu_indices(len(matrix), k=1)
            reinforced, _ = _reinforce_with_gene_prior(
                base[upper_rows, upper_cols],
                matrix.index,
                upper_rows,
                upper_cols,
                gene_prior,
                gene_prior_weight,
            )
            fused = base.copy()
            fused[upper_rows, upper_cols] = reinforced
            fused[upper_cols, upper_rows] = reinforced
            distance = 1.0 - fused
        else:
            distance = 1.0 - correlation if settings.signed else 1.0 - np.abs(correlation)
        distance = np.clip((distance + distance.T) / 2.0, 0.0, 2.0)
        np.fill_diagonal(distance, 0.0)
        hierarchy = linkage(squareform(distance, checks=False), method="average")
        raw_labels = fcluster(hierarchy, t=settings.module_cut_height, criterion="distance")
    labels, _ = _module_labels(raw_labels, settings.min_module_size)
    if settings.module_refinement:
        labels, _ = _refine_module_labels(matrix, labels, settings)
    return labels


def fit_network(
    matrix: pd.DataFrame,
    qc_metrics: pd.DataFrame,
    settings: NetworkConfig,
    samples: pd.DataFrame | None = None,
    event_reliability: pd.Series | None = None,
    reliability_floor: float = 0.10,
    reliability_power: float = 0.50,
    gene_prior: GeneMantelPrior | None = None,
    gene_prior_weight: float = 0.0,
    precomputed_correlation: np.ndarray | None = None,
    correlation_sd: np.ndarray | None = None,
    edge_sign_probability: np.ndarray | None = None,
    min_edge_sign_probability: float = 0.80,
    precomputed_event_keys: pd.Index | np.ndarray | None = None,
) -> NetworkResult:
    if not matrix.index.is_unique or not matrix.columns.is_unique:
        raise ValueError("Network matrix event and sample labels must be unique.")
    if matrix.isna().any().any() or not np.isfinite(matrix.to_numpy(dtype=float)).all():
        raise ValueError("Network matrix must be finite after Stable-AS imputation.")
    if "stability_score" not in qc_metrics.columns:
        raise ValueError("QC metrics must contain stability_score.")
    missing_metrics = sorted(set(matrix.index).difference(qc_metrics.index))
    if missing_metrics:
        raise ValueError(f"QC metrics are missing network events: {missing_metrics[:10]}")
    if samples is not None:
        if "sample_id" not in samples.columns or not samples["sample_id"].is_unique:
            raise ValueError("Network sample metadata requires a unique sample_id column.")
        if set(matrix.columns) != set(samples["sample_id"]):
            raise ValueError("Network matrix columns must match sample metadata exactly.")
    minimum_samples = settings.minimum_network_samples
    if isinstance(minimum_samples, bool) or not isinstance(minimum_samples, int) or minimum_samples < 3:
        raise ValueError("minimum_network_samples must be an integer of at least 3")
    if matrix.shape[1] < minimum_samples:
        raise ValueError(f"At least {minimum_samples} samples are required by the configured network gate; observed {matrix.shape[1]}")
    if matrix.shape[1] < 15:
        import warnings
        warnings.warn(
            f"Exploratory network with {matrix.shape[1]} samples: the default 15-sample gate was explicitly overridden. "
            "Bootstrap stability does not replace independent validation.",
            UserWarning,
            stacklevel=2,
        )
    if matrix.shape[0] < settings.min_module_size * 2:
        raise ValueError(f"Too few stable events for network fitting: {matrix.shape[0]}")

    stability = qc_metrics.reindex(matrix.index)["stability_score"].fillna(0.0)
    reliability: pd.Series | None = None
    if event_reliability is not None:
        reliability = (
            pd.to_numeric(event_reliability.reindex(matrix.index), errors="coerce")
            .fillna(0.0)
            .clip(reliability_floor, 1.0)
        )
        ranking = (stability * reliability.pow(reliability_power)).sort_values(
            ascending=False
        )
    else:
        ranking = stability.sort_values(ascending=False)
    selected_keys = ranking.index if settings.max_events is None else ranking.index[: settings.max_events]
    selected = matrix.loc[selected_keys].copy()
    if (
        precomputed_correlation is not None
        or correlation_sd is not None
        or edge_sign_probability is not None
    ):
        source_keys = (
            pd.Index(np.asarray(precomputed_event_keys, dtype=str))
            if precomputed_event_keys is not None
            else matrix.index.astype(str)
        )
        full_index = source_keys.get_indexer(selected.index.astype(str))
        if np.any(full_index < 0):
            raise ValueError("precomputed correlation is missing selected network events")
        if precomputed_correlation is not None:
            precomputed_correlation = np.asarray(precomputed_correlation)[
                np.ix_(full_index, full_index)
            ]
        if correlation_sd is not None:
            correlation_sd = np.asarray(correlation_sd)[np.ix_(full_index, full_index)]
        if edge_sign_probability is not None:
            edge_sign_probability = np.asarray(edge_sign_probability)[
                np.ix_(full_index, full_index)
            ]
    linkage_matrix: np.ndarray | None = None
    if settings.engine == "tom_hierarchical":
        selected_reliability = (
            reliability.reindex(selected.index).to_numpy(dtype=float)
            if reliability is not None
            else None
        )
        (
            correlation,
            adjacency,
            tom,
            linkage_matrix,
            base_as_adjacency,
            gene_prior_values,
        ) = _tom_hierarchy(
            selected,
            settings,
            event_reliability=selected_reliability,
            reliability_power=reliability_power,
            gene_prior=gene_prior,
            gene_prior_weight=gene_prior_weight,
            precomputed_correlation=precomputed_correlation,
        )
        if settings.module_detection == "wgcna_blockwise":
            from katana_cosplice.r_wgcna import exact_wgcna_labels

            module_labels = exact_wgcna_labels(
                selected,
                soft_power=settings.soft_power,
                min_module_size=settings.min_module_size,
                merge_cut_height=settings.wgcna_merge_cut_height,
                deep_split=settings.dynamic_deep_split,
                pam_respects_dendro=settings.dynamic_pam_respects_dendro,
                cor_type=settings.method,
                signed=settings.signed,
                rscript=settings.wgcna_rscript,
                event_reliability=(
                    reliability.reindex(selected.index)
                    if reliability is not None
                    else None
                ),
                reliability_power=reliability_power,
                precomputed_correlation=precomputed_correlation,
            )
            raw_labels = None
        elif settings.module_detection == "dynamic_tree_cut":
            distance = np.clip(1.0 - tom, 0.0, 1.0)
            np.fill_diagonal(distance, 0.0)
            raw_labels = dynamic_tree_cut_labels(
                linkage_matrix,
                distance,
                min_module_size=settings.min_module_size,
                deep_split=settings.dynamic_deep_split,
                cut_height=settings.dynamic_cut_height,
                pam_respects_dendro=settings.dynamic_pam_respects_dendro,
            )
        else:
            raw_labels = fcluster(
                linkage_matrix,
                t=settings.module_cut_height,
                criterion="distance",
            )
        if raw_labels is not None:
            module_labels, _ = _module_labels(raw_labels, settings.min_module_size)
        connectivity = adjacency.sum(axis=1)
        gene_prior_matrix = np.zeros_like(adjacency)
        upper_rows, upper_cols = np.triu_indices(len(selected), k=1)
        if len(gene_prior_values):
            gene_prior_matrix[upper_rows, upper_cols] = gene_prior_values
            gene_prior_matrix[upper_cols, upper_rows] = gene_prior_values
        edges = _edge_table(
            correlation,
            adjacency,
            selected.index.tolist(),
            settings.edge_threshold,
            settings.max_output_edges,
            base_as_adjacency=base_as_adjacency,
            gene_prior_matrix=gene_prior_matrix if gene_prior is not None else None,
        )
        if not edges.empty and correlation_sd is not None:
            event_index = {key: index for index, key in enumerate(selected.index)}
            rows = edges["event_a"].map(event_index).to_numpy(dtype=int)
            cols = edges["event_b"].map(event_index).to_numpy(dtype=int)
            edges["posterior_correlation_sd"] = correlation_sd[rows, cols]
        if not edges.empty and edge_sign_probability is not None:
            event_index = {key: index for index, key in enumerate(selected.index)}
            rows = edges["event_a"].map(event_index).to_numpy(dtype=int)
            cols = edges["event_b"].map(event_index).to_numpy(dtype=int)
            edges["posterior_sign_probability"] = edge_sign_probability[rows, cols]
            edges["posterior_direction_stable"] = edges[
                "posterior_sign_probability"
            ].ge(min_edge_sign_probability)
        sensitivity = (
            _dynamic_tree_sensitivity(linkage_matrix, tom, settings)
            if settings.module_detection in {"wgcna_blockwise", "dynamic_tree_cut"}
            else _cut_height_sensitivity(linkage_matrix, settings)
        )
        n_graph_edges = int(
            np.count_nonzero(
                np.triu(adjacency >= settings.edge_threshold, k=1)
            )
        )
    elif settings.engine == "exact_hierarchical":
        correlation = correlation_matrix(selected, settings.method)
        base = (1.0 + correlation) / 2.0 if settings.signed else np.abs(correlation)
        adjacency = np.power(np.clip(base, 0.0, 1.0), settings.soft_power)
        if reliability is not None:
            selected_reliability = reliability.reindex(selected.index).to_numpy(dtype=float)
            reliability_factor = np.power(
                np.outer(selected_reliability, selected_reliability),
                reliability_power / 2.0,
            )
            adjacency *= reliability_factor
        base_as_adjacency = adjacency.copy()
        upper_rows, upper_cols = np.triu_indices(len(selected), k=1)
        reinforced, gene_prior_values = _reinforce_with_gene_prior(
            adjacency[upper_rows, upper_cols],
            selected.index,
            upper_rows,
            upper_cols,
            gene_prior,
            gene_prior_weight,
        )
        if len(reinforced):
            adjacency[upper_rows, upper_cols] = reinforced
            adjacency[upper_cols, upper_rows] = reinforced
        np.fill_diagonal(adjacency, 0.0)
        distance = (
            1.0 - adjacency
            if gene_prior is not None and gene_prior_weight > 0
            else (1.0 - correlation if settings.signed else 1.0 - np.abs(correlation))
        )
        distance = np.clip((distance + distance.T) / 2.0, 0.0, 2.0)
        np.fill_diagonal(distance, 0.0)
        linkage_matrix = linkage(squareform(distance, checks=False), method="average")
        raw_labels = fcluster(linkage_matrix, t=settings.module_cut_height, criterion="distance")
        module_labels, _ = _module_labels(raw_labels, settings.min_module_size)
        connectivity = adjacency.sum(axis=1)
        gene_prior_matrix = np.zeros_like(adjacency)
        if len(gene_prior_values):
            gene_prior_matrix[upper_rows, upper_cols] = gene_prior_values
            gene_prior_matrix[upper_cols, upper_rows] = gene_prior_values
        edges = _edge_table(
            correlation,
            adjacency,
            selected.index.tolist(),
            settings.edge_threshold,
            settings.max_output_edges,
            base_as_adjacency=base_as_adjacency,
            gene_prior_matrix=gene_prior_matrix if gene_prior is not None else None,
        )
        sensitivity = _cut_height_sensitivity(linkage_matrix, settings)
        n_graph_edges = int(np.count_nonzero(np.triu(adjacency >= settings.edge_threshold, k=1)))
    else:
        rows, cols, correlations, adjacencies, overlaps, consensus_support = _network_edges(
            selected,
            settings,
        )
        if reliability is not None and len(rows):
            selected_reliability = reliability.reindex(selected.index).to_numpy(dtype=float)
            reliability_factor = np.power(
                selected_reliability[rows] * selected_reliability[cols],
                reliability_power / 2.0,
            )
            adjacencies *= reliability_factor
            retained_reliability = adjacencies >= settings.edge_threshold
            rows = rows[retained_reliability]
            cols = cols[retained_reliability]
            correlations = correlations[retained_reliability]
            adjacencies = adjacencies[retained_reliability]
            overlaps = overlaps[retained_reliability]
            consensus_support = consensus_support[retained_reliability]
        adjacencies, graph_bootstrap_support = _bootstrap_weighted_edges(
            selected,
            rows,
            cols,
            adjacencies,
            settings,
            samples,
        )
        base_as_adjacencies = adjacencies.copy()
        adjacencies, gene_prior_values = _reinforce_with_gene_prior(
            adjacencies,
            selected.index,
            rows,
            cols,
            gene_prior,
            gene_prior_weight,
        )
        n_unstable_graph_edges_removed = 0
        if len(graph_bootstrap_support) and np.isfinite(graph_bootstrap_support).any():
            retained = graph_bootstrap_support >= settings.consensus_bootstrap_min_support
            n_unstable_graph_edges_removed = int((~retained).sum())
            rows = rows[retained]
            cols = cols[retained]
            correlations = correlations[retained]
            adjacencies = adjacencies[retained]
            overlaps = overlaps[retained]
            consensus_support = consensus_support[retained]
            graph_bootstrap_support = graph_bootstrap_support[retained]
            base_as_adjacencies = base_as_adjacencies[retained]
            gene_prior_values = gene_prior_values[retained]
        raw_labels, graph = _graph_communities(len(selected), rows, cols, adjacencies, settings)
        module_labels, _ = _module_labels(raw_labels, settings.min_module_size)
        connectivity = np.asarray([graph.degree(index, weight="weight") for index in range(len(selected))], dtype=float)
        edges = _sparse_edge_table(
            rows,
            cols,
            correlations,
            adjacencies,
            overlaps,
            consensus_support,
            graph_bootstrap_support,
            selected.index.tolist(),
            settings.max_output_edges,
            base_as_adjacencies,
            gene_prior_values,
        )
        sensitivity = _community_resolution_sensitivity(len(selected), rows, cols, adjacencies, settings)
        n_graph_edges = len(rows)
    initial_module_labels = module_labels.copy()
    if settings.module_refinement and settings.module_detection != "wgcna_blockwise":
        module_labels, membership = _refine_module_labels(selected, module_labels, settings)
    else:
        membership = _membership_statistics(selected, module_labels, settings)
        membership.insert(0, "module_before_refinement", initial_module_labels)
        membership["module_after_refinement"] = module_labels
        membership["membership_cutoff_used"] = settings.min_module_membership
        membership["module_refinement_action"] = "disabled"
    scores = _program_scores(selected, module_labels)
    assignments = pd.DataFrame(
        {
            "event_key": selected.index,
            "module": module_labels,
            "network_connectivity": connectivity,
            "stability_score": stability.reindex(selected.index).to_numpy(),
        }
    )
    assignments["network_priority_score"] = ranking.reindex(selected.index).to_numpy()
    assignments["network_reliability"] = (
        reliability.reindex(selected.index).to_numpy()
        if reliability is not None
        else np.nan
    )
    if gene_prior is not None:
        selected_genes = gene_prior.event_genes.reindex(selected.index.astype(str))
        assignments["gene_prior_gene_id"] = selected_genes.to_numpy()
        assignments["gene_prior_available"] = selected_genes.isin(
            gene_prior.gene_adjacency.index
        ).to_numpy()
    else:
        assignments["gene_prior_gene_id"] = ""
        assignments["gene_prior_available"] = False
    assignments = assignments.join(
        membership[
            [
                "module_before_refinement",
                "module_membership",
                "best_module_membership",
                "second_best_module_membership",
                "module_membership_margin",
                "membership_cutoff_used",
                "module_refinement_action",
            ]
        ],
        on="event_key",
    )
    assignments["module_membership_rank"] = assignments.groupby("module", observed=True)["network_connectivity"].rank(method="dense", ascending=False)
    edges["bootstrap_support"] = _edge_bootstrap_support(selected, edges, settings, samples)
    if "graph_bootstrap_support" not in edges:
        # Dense TOM has no separate sparse-graph bootstrap. Edge survival under
        # the same condition-stratified resampling is the equivalent provenance.
        edges["graph_bootstrap_support"] = edges["bootstrap_support"]
    edges["bootstrap_stable"] = edges["bootstrap_support"].ge(settings.bootstrap_edge_threshold)
    if not edges.empty:
        incident = pd.concat(
            [
                edges[["event_a", "bootstrap_support"]].rename(columns={"event_a": "event_key"}),
                edges[["event_b", "bootstrap_support"]].rename(columns={"event_b": "event_key"}),
            ],
            ignore_index=True,
        )
        mean_support = incident.groupby("event_key", observed=True)["bootstrap_support"].mean()
        assignments["mean_incident_bootstrap_support"] = assignments["event_key"].map(mean_support)
    else:
        assignments["mean_incident_bootstrap_support"] = np.nan
    module_sizes = assignments.groupby("module", observed=True).size().rename("n_events")
    diagnostics = module_sizes.reset_index()
    diagnostics["n_samples"] = selected.shape[1]
    diagnostics["n_input_stable_events"] = matrix.shape[0]
    diagnostics["n_selected_events"] = selected.shape[0]
    diagnostics["selection_rule"] = (
        "all_qc_passed"
        if settings.max_events is None
        else (
            f"top_{settings.max_events}_by_stability_x_reliability"
            if reliability is not None
            else f"top_{settings.max_events}_by_stability_score"
        )
    )
    diagnostics["uncertainty_weighting"] = reliability is not None
    diagnostics["posterior_ensemble_correlation"] = precomputed_correlation is not None
    diagnostics["posterior_ensemble_edge_diagnostics"] = (
        correlation_sd is not None and edge_sign_probability is not None
    )
    diagnostics["min_edge_sign_probability"] = (
        min_edge_sign_probability if edge_sign_probability is not None else np.nan
    )
    diagnostics["posterior_direction_stable_edge_fraction"] = (
        float(edges["posterior_direction_stable"].mean())
        if "posterior_direction_stable" in edges and not edges.empty
        else np.nan
    )
    diagnostics["reliability_floor"] = reliability_floor if reliability is not None else np.nan
    diagnostics["reliability_power"] = reliability_power if reliability is not None else np.nan
    diagnostics["median_network_reliability"] = (
        float(reliability.reindex(selected.index).median())
        if reliability is not None
        else np.nan
    )
    diagnostics["gene_mantel_prior_enabled"] = gene_prior is not None
    diagnostics["gene_mantel_prior_weight"] = (
        gene_prior_weight if gene_prior is not None else 0.0
    )
    diagnostics["gene_mantel_prior_event_coverage"] = (
        float(assignments["gene_prior_available"].mean())
        if gene_prior is not None
        else 0.0
    )
    diagnostics["gene_mantel_prior_edge_coverage"] = (
        float(np.mean(gene_prior_values > 0))
        if gene_prior is not None and len(gene_prior_values)
        else 0.0
    )
    diagnostics["engine"] = settings.engine
    diagnostics["method"] = settings.method
    diagnostics["soft_power"] = settings.soft_power
    diagnostics["signed_network"] = settings.signed
    diagnostics["bootstrap_iterations"] = settings.bootstrap_iterations
    diagnostics["network_backbone"] = (
        "signed_soft_adjacency_exact_TOM_average_linkage"
        if settings.engine == "tom_hierarchical"
        else settings.engine
    )
    diagnostics["module_detection"] = settings.module_detection
    diagnostics["dynamic_tree_backend"] = (
        "official_R_WGCNA"
        if settings.module_detection == "wgcna_blockwise"
        else (
            "dynamicTreeCut_python_0.1.1"
            if settings.module_detection == "dynamic_tree_cut"
            else "not_used"
        )
    )
    diagnostics["dynamic_cut_height"] = (
        settings.dynamic_cut_height
        if settings.module_detection in {"wgcna_blockwise", "dynamic_tree_cut"}
        else np.nan
    )
    diagnostics["dynamic_deep_split"] = (
        settings.dynamic_deep_split
        if settings.module_detection in {"wgcna_blockwise", "dynamic_tree_cut"}
        else np.nan
    )
    diagnostics["dynamic_pam_respects_dendro"] = (
        settings.dynamic_pam_respects_dendro
        if settings.module_detection in {"wgcna_blockwise", "dynamic_tree_cut"}
        else False
    )
    diagnostics["tom_dense_max_events"] = settings.tom_dense_max_events
    diagnostics["knn_neighbors"] = settings.knn_neighbors if settings.engine == "sparse_knn" else np.nan
    diagnostics["edge_smoothing"] = settings.edge_smoothing if settings.engine == "sparse_knn" else "none"
    diagnostics["shared_neighbor_weight"] = settings.shared_neighbor_weight if settings.engine == "sparse_knn" else np.nan
    diagnostics["shared_neighbor_top_k"] = settings.shared_neighbor_top_k if settings.engine == "sparse_knn" else np.nan
    diagnostics["multi_k_consensus"] = settings.multi_k_consensus if settings.engine == "sparse_knn" else False
    diagnostics["consensus_knn_neighbors"] = (
        ",".join(str(value) for value in settings.consensus_knn_neighbors)
        if settings.engine == "sparse_knn"
        else ""
    )
    diagnostics["consensus_edge_min_fraction"] = settings.consensus_edge_min_fraction
    diagnostics["consensus_bootstrap_iterations"] = settings.consensus_bootstrap_iterations
    diagnostics["effective_consensus_bootstrap_iterations"] = (
        _effective_consensus_bootstrap_iterations(len(selected), settings)
        if settings.engine == "sparse_knn"
        else 0
    )
    diagnostics["consensus_bootstrap_weight"] = settings.consensus_bootstrap_weight
    diagnostics["consensus_bootstrap_min_support"] = settings.consensus_bootstrap_min_support
    diagnostics["n_unstable_graph_edges_removed"] = (
        n_unstable_graph_edges_removed if settings.engine == "sparse_knn" else 0
    )
    diagnostics["community_consensus_runs"] = settings.community_consensus_runs
    diagnostics["community_consensus_weight"] = settings.community_consensus_weight
    diagnostics["mean_graph_bootstrap_support"] = (
        float(np.nanmean(graph_bootstrap_support))
        if settings.engine == "sparse_knn" and len(graph_bootstrap_support) and np.isfinite(graph_bootstrap_support).any()
        else np.nan
    )
    diagnostics["community_resolution"] = settings.community_resolution if settings.engine == "sparse_knn" else np.nan
    diagnostics["module_refinement"] = settings.module_refinement
    diagnostics["min_module_membership"] = settings.min_module_membership
    diagnostics["min_module_membership_margin"] = settings.min_module_membership_margin
    diagnostics["adaptive_membership"] = settings.adaptive_membership
    diagnostics["adaptive_membership_floor"] = settings.adaptive_membership_floor
    diagnostics["effective_adaptive_membership_floor"] = min(
        settings.min_module_membership,
        max(settings.adaptive_membership_floor, 1.0 / np.sqrt(selected.shape[1])),
    )
    diagnostics["adaptive_membership_rule"] = (
        "module_robust_lower_bound_bounded_by_sample_size_floor_and_fixed_ceiling"
        if settings.adaptive_membership
        else "fixed_min_module_membership"
    )
    diagnostics["mean_membership_cutoff_used"] = float(
        pd.to_numeric(assignments["membership_cutoff_used"], errors="coerce").mean()
    )
    diagnostics["module_merge_correlation"] = settings.module_merge_correlation
    diagnostics["refinement_iterations"] = settings.refinement_iterations
    refinement_actions = assignments["module_refinement_action"].value_counts()
    diagnostics["n_events_greyed_by_refinement"] = int(
        np.count_nonzero((initial_module_labels != "grey") & (module_labels == "grey"))
    )
    diagnostics["n_events_assigned_from_grey"] = int(
        np.count_nonzero((initial_module_labels == "grey") & (module_labels != "grey"))
    )
    diagnostics["n_events_reassigned_by_refinement"] = int(refinement_actions.get("reassigned_module", 0))
    diagnostics["n_events_in_merged_modules"] = int(refinement_actions.get("merged_module", 0))
    diagnostics["n_initial_modules"] = len(set(initial_module_labels).difference({"grey"}))
    diagnostics["n_final_modules"] = len(set(module_labels).difference({"grey"}))
    diagnostics["n_graph_edges"] = n_graph_edges
    diagnostics["n_output_edges"] = len(edges)
    diagnostics["n_bootstrap_stable_edges"] = int(edges["bootstrap_stable"].sum())
    retained_sizes = diagnostics.loc[diagnostics["module"].ne("grey"), "n_events"]
    largest_retained = int(retained_sizes.max()) if not retained_sizes.empty else 0
    diagnostics["largest_module_fraction"] = largest_retained / selected.shape[0]
    diagnostics["dominant_module_warning"] = diagnostics["largest_module_fraction"].ge(0.60)
    return NetworkResult(assignments, edges, scores, diagnostics, sensitivity, selected, linkage_matrix)


def render_network_plots(result: NetworkResult, samples: pd.DataFrame, output_dir: Path, network_name: str) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    generated: list[Path] = []

    fig, ax = plt.subplots(figsize=(11, 6.5))
    sizes = result.assignments.loc[result.assignments["module"].ne("grey")].groupby("module", observed=True).size().sort_values(ascending=False)
    if sizes.empty:
        ax.text(
            0.5,
            0.5,
            "No non-grey module passed the production confidence rules",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        ax.set_xticks([])
    else:
        sizes.plot.bar(ax=ax, color="#087E8B", edgecolor="#17324D")
    ax.set_title(f"{network_name.title()} co-splicing module sizes", fontweight="bold")
    ax.set_xlabel("Module")
    ax.set_ylabel("Stable AS events")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    path = output_dir / f"{network_name}_module_sizes.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    generated.append(path)

    if not result.program_scores.empty:
        ordered_samples = samples.set_index("sample_id").reindex(result.program_scores.index)
        conditions = ordered_samples["condition"].astype(str)
        order = np.argsort(conditions.to_numpy(), kind="stable")
        plot_frame = result.program_scores.iloc[order].T
        width = max(12.0, plot_frame.shape[1] * 0.32)
        height = max(5.5, plot_frame.shape[0] * 0.42)
        fig, ax = plt.subplots(figsize=(width, height))
        sns.heatmap(plot_frame, cmap="vlag", center=0.0, ax=ax, cbar_kws={"label": "Program score"})
        ax.set_title(f"{network_name.title()} module program scores", fontweight="bold", pad=14)
        ax.set_xlabel("Samples ordered by condition")
        ax.set_ylabel("Co-splicing module")
        ax.tick_params(axis="x", rotation=90)
        fig.tight_layout()
        path = output_dir / f"{network_name}_program_score_heatmap.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        generated.append(path)

    if result.linkage_matrix is not None:
        fig, ax = plt.subplots(figsize=(13, 6.5))
        dendrogram(result.linkage_matrix, no_labels=True, color_threshold=0.7, ax=ax)
        ax.set_title(f"{network_name.title()} event dendrogram", fontweight="bold")
        ax.set_xlabel("Stable AS events")
        ax.set_ylabel("Correlation distance")
        fig.tight_layout()
        path = output_dir / f"{network_name}_event_dendrogram.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        generated.append(path)
    return generated
