from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import TYPE_CHECKING, Callable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.metrics import adjusted_rand_score

try:
    import igraph as ig
    import leidenalg
except ImportError:  # pragma: no cover - exercised only in minimal installations
    ig = None
    leidenalg = None

from katana_cosplice.config import ConsensusConfig, NetworkConfig
from katana_cosplice.network import _program_scores, discover_module_labels

if TYPE_CHECKING:
    from katana_cosplice.gene_prior import GeneMantelPrior


Progress = Callable[[str], None] | None
CandidateScorer = Callable[[np.ndarray], dict[str, float | int]]

@dataclass
class ConsensusResult:
    assignments: pd.DataFrame
    module_summary: pd.DataFrame
    program_scores: pd.DataFrame
    bootstrap_runs: pd.DataFrame
    coassignment_edges: pd.DataFrame
    coassignment_matrix: np.ndarray
    event_keys: np.ndarray
    diagnostics: pd.DataFrame
    parameter_search: pd.DataFrame
    selection_audit: pd.DataFrame


def _bootstrap_indices(samples: pd.DataFrame, columns: pd.Index, rng: np.random.Generator) -> np.ndarray:
    ordered = samples.set_index("sample_id").reindex(columns)
    if "condition" not in ordered:
        return rng.choice(np.arange(len(columns)), size=len(columns), replace=True)
    conditions = ordered["condition"].fillna("unknown").astype(str).to_numpy()
    groups = [np.flatnonzero(conditions == value) for value in pd.unique(conditions)]
    return np.concatenate([rng.choice(group, size=len(group), replace=True) for group in groups])


def _canonical_modules(raw_labels: np.ndarray, min_module_size: int) -> np.ndarray:
    counts = pd.Series(raw_labels).value_counts()
    retained = [int(label) for label, count in counts.items() if int(count) >= min_module_size]
    retained.sort(key=lambda label: (-int(counts[label]), label))
    names = {label: f"C{index:02d}" for index, label in enumerate(retained, start=1)}
    return np.asarray([names.get(int(label), "grey") for label in raw_labels], dtype=object)


def _final_modules(coassignment: np.ndarray, threshold: float, min_module_size: int) -> np.ndarray:
    n_events = len(coassignment)
    if n_events == 1:
        return np.asarray(["grey"], dtype=object)
    distance = np.clip(1.0 - coassignment, 0.0, 1.0)
    distance = (distance + distance.T) / 2.0
    np.fill_diagonal(distance, 0.0)
    hierarchy = linkage(squareform(distance, checks=False), method="average")
    raw = fcluster(hierarchy, t=1.0 - threshold, criterion="distance")
    return _canonical_modules(raw, min_module_size)


def _confidence_tables(
    coassignment: np.ndarray,
    event_keys: np.ndarray,
    labels: np.ndarray,
    baseline: pd.Series,
    membership_role: np.ndarray | None = None,
    candidate_labels: np.ndarray | None = None,
    confidence: np.ndarray | None = None,
    second_best: np.ndarray | None = None,
    margin: np.ndarray | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    modules = sorted(set(labels).difference({"grey"}))
    indices = {module: np.flatnonzero(labels == module) for module in modules}
    rows: list[dict[str, object]] = []
    for event_index, event_key in enumerate(event_keys):
        module = str(labels[event_index])
        own_values = np.asarray([], dtype=float)
        if module != "grey":
            own = indices[module]
            own = own[own != event_index]
            own_values = coassignment[event_index, own]
        own_mean = float(own_values.mean()) if len(own_values) else 0.0
        other_means = [
            float(coassignment[event_index, module_indices].mean())
            for other_module, module_indices in indices.items()
            if other_module != module and len(module_indices)
        ]
        second = max(other_means, default=0.0)
        if confidence is not None:
            own_mean = float(confidence[event_index])
        if second_best is not None:
            second = float(second_best[event_index])
        rows.append(
            {
                "event_key": str(event_key),
                "observed_module": str(baseline.get(str(event_key), "grey")),
                "candidate_graph_module": (
                    str(candidate_labels[event_index]) if candidate_labels is not None else module
                ),
                "consensus_module": module,
                "membership_role": (
                    str(membership_role[event_index]) if membership_role is not None else "core"
                ),
                "consensus_membership_confidence": own_mean,
                "next_best_module_coassignment": second,
                "consensus_confidence_margin": (
                    float(margin[event_index]) if margin is not None else own_mean - second
                ),
                "assignment_uncertainty": 1.0 - own_mean,
            }
        )
    assignments = pd.DataFrame(rows)
    assignments["analysis_module"] = assignments["consensus_module"]
    uncertain = assignments["consensus_module"].eq("grey") & assignments[
        "candidate_graph_module"
    ].ne("grey")
    assignments.loc[uncertain, "analysis_module"] = assignments.loc[
        uncertain, "candidate_graph_module"
    ]
    assignments["analysis_assignment_scope"] = np.where(
        assignments["consensus_module"].ne("grey"),
        "conservative_consensus",
        np.where(uncertain, "exploratory_uncertain", "unassigned"),
    )
    module_rows: list[dict[str, object]] = []
    for module, module_indices in indices.items():
        if len(module_indices) > 1:
            within = coassignment[np.ix_(module_indices, module_indices)]
            within = within[np.triu_indices(len(module_indices), k=1)]
        else:
            within = np.asarray([], dtype=float)
        event_confidence = assignments.loc[
            assignments["consensus_module"].eq(module), "consensus_membership_confidence"
        ]
        module_rows.append(
            {
                "consensus_module": module,
                "n_events": int(len(module_indices)),
                "mean_within_module_coassignment": float(within.mean()) if len(within) else 0.0,
                "median_within_module_coassignment": float(np.median(within)) if len(within) else 0.0,
                "min_within_module_coassignment": float(within.min()) if len(within) else 0.0,
                "mean_event_confidence": float(event_confidence.mean()),
                "min_event_confidence": float(event_confidence.min()),
                "n_core_events": int(
                    assignments.loc[
                        assignments["consensus_module"].eq(module), "membership_role"
                    ].eq("core").sum()
                ),
                "n_peripheral_events": int(
                    assignments.loc[
                        assignments["consensus_module"].eq(module), "membership_role"
                    ].eq("peripheral").sum()
                ),
                "n_uncertain_candidate_events": int(
                    (
                        assignments["candidate_graph_module"].eq(module)
                        & assignments["membership_role"].eq("uncertain")
                    ).sum()
                ),
            }
        )
    return assignments, pd.DataFrame(module_rows)


def _top_k_edges(coassignment: np.ndarray, top_k: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n_events = len(coassignment)
    if n_events < 2:
        return np.asarray([], dtype=int), np.asarray([], dtype=int), np.asarray([], dtype=float)
    k = min(int(top_k), n_events - 1)
    edge_weights: dict[tuple[int, int], float] = {}
    for left in range(n_events):
        row = coassignment[left].copy()
        row[left] = -np.inf
        neighbors = np.argpartition(row, -k)[-k:]
        for right in neighbors:
            source, target = sorted((left, int(right)))
            if source == target:
                continue
            edge_weights[(source, target)] = max(
                edge_weights.get((source, target), 0.0), float(coassignment[source, target])
            )
    pairs = sorted(edge_weights)
    return (
        np.asarray([pair[0] for pair in pairs], dtype=int),
        np.asarray([pair[1] for pair in pairs], dtype=int),
        np.asarray([edge_weights[pair] for pair in pairs], dtype=float),
    )


def _graph_partition(
    n_events: int,
    source: np.ndarray,
    target: np.ndarray,
    weights: np.ndarray,
    threshold: float,
    resolution: float,
    seed: int,
) -> tuple[np.ndarray, float, int, str]:
    retained = weights >= threshold
    edges = list(zip(source[retained].tolist(), target[retained].tolist(), strict=False))
    selected_weights = weights[retained].astype(float).tolist()
    if ig is not None and leidenalg is not None:
        graph = ig.Graph(n=n_events, edges=edges, directed=False)
        graph.es["weight"] = selected_weights
        partition = leidenalg.find_partition(
            graph,
            leidenalg.RBConfigurationVertexPartition,
            weights="weight" if edges else None,
            resolution_parameter=float(resolution),
            seed=int(seed),
        )
        labels = np.asarray(partition.membership, dtype=int)
        modularity = float(graph.modularity(labels.tolist(), weights="weight")) if edges else 0.0
        return labels, modularity, len(edges), "leiden"

    graph = nx.Graph()
    graph.add_nodes_from(range(n_events))
    graph.add_weighted_edges_from(
        (left, right, weight)
        for (left, right), weight in zip(edges, selected_weights, strict=False)
    )
    communities = nx.community.louvain_communities(
        graph, weight="weight", resolution=float(resolution), seed=int(seed)
    )
    labels = np.full(n_events, -1, dtype=int)
    for label, community in enumerate(communities):
        labels[list(community)] = label
    modularity = (
        float(nx.community.modularity(graph, communities, weight="weight", resolution=float(resolution)))
        if edges and communities
        else 0.0
    )
    return labels, modularity, len(edges), "networkx_louvain_fallback"


def _soft_membership(
    coassignment: np.ndarray,
    candidate_labels: np.ndarray,
    min_module_size: int,
    settings: ConsensusConfig,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    modules = sorted(set(candidate_labels).difference({"grey"}))
    if not modules:
        empty = np.zeros(len(candidate_labels), dtype=float)
        return (
            np.full(len(candidate_labels), "grey", dtype=object),
            np.full(len(candidate_labels), "grey", dtype=object),
            empty,
            empty.copy(),
            empty.copy(),
            candidate_labels,
        )
    profiles = np.zeros((len(candidate_labels), len(modules)), dtype=float)
    for module_index, module in enumerate(modules):
        members = np.flatnonzero(candidate_labels == module)
        sums = coassignment[:, members].sum(axis=1)
        profiles[:, module_index] = sums / max(len(members), 1)
        if len(members) > 1:
            profiles[members, module_index] = (sums[members] - 1.0) / (len(members) - 1)
        else:
            profiles[members, module_index] = 0.0
    order = np.argsort(profiles, axis=1)
    best_index = order[:, -1]
    best = profiles[np.arange(len(profiles)), best_index]
    second = profiles[np.arange(len(profiles)), order[:, -2]] if len(modules) > 1 else np.zeros(len(profiles))
    margin = best - second
    best_module = np.asarray([modules[index] for index in best_index], dtype=object)
    labels = np.full(len(candidate_labels), "grey", dtype=object)
    # A graph assignment and a conservative consensus assignment answer
    # different questions. Keep low-confidence graph membership auditable
    # instead of presenting every rejected event as if it had no candidate
    # program at all.
    roles = np.full(len(candidate_labels), "uncertain", dtype=object)
    roles[candidate_labels == "grey"] = "grey"
    is_core = (candidate_labels == best_module) & (best >= settings.core_membership_min)
    is_peripheral = (~is_core) & (best >= settings.peripheral_membership_min) & (
        margin >= settings.peripheral_membership_margin
    )
    labels[is_core] = best_module[is_core]
    roles[is_core] = "core"
    labels[is_peripheral] = best_module[is_peripheral]
    roles[is_peripheral] = "peripheral"

    counts = pd.Series(labels[labels != "grey"]).value_counts()
    removed = set(counts[counts < min_module_size].index.astype(str))
    if removed:
        drop = np.isin(labels, list(removed))
        labels[drop] = "grey"
        roles[drop] = "uncertain"
    counts = pd.Series(labels[labels != "grey"]).value_counts()
    ordered = sorted(counts.index.astype(str), key=lambda value: (-int(counts[value]), value))
    rename = {module: f"C{index:02d}" for index, module in enumerate(ordered, start=1)}
    labels = np.asarray([rename.get(str(module), "grey") for module in labels], dtype=object)
    candidate_display = np.asarray(
        [rename.get(str(module), str(module)) if str(module) != "grey" else "grey" for module in candidate_labels],
        dtype=object,
    )
    return labels, roles, best, second, margin, candidate_display


def _candidate_metrics(
    labels: np.ndarray,
    confidence: np.ndarray,
    margin: np.ndarray,
    modularity: float,
    settings: ConsensusConfig,
    bootstrap_labels: list[np.ndarray] | None = None,
) -> dict[str, float | int | bool]:
    non_grey = labels != "grey"
    sizes = pd.Series(labels[non_grey]).value_counts()
    n_events = len(labels)
    n_modules = int(len(sizes))
    coverage = float(non_grey.mean())
    grey_fraction = 1.0 - coverage
    largest = float(sizes.max() / n_events) if not sizes.empty else 0.0
    if n_modules > 1:
        proportions = sizes.to_numpy(dtype=float) / sizes.sum()
        balance = float(-(proportions * np.log(proportions)).sum() / np.log(n_modules))
    else:
        balance = 0.0
    purity = float(confidence[non_grey].mean()) if non_grey.any() else 0.0
    separation = float(np.clip(margin[non_grey], 0.0, 1.0).mean()) if non_grey.any() else 0.0
    modularity_scaled = float(np.clip((modularity + 1.0) / 2.0, 0.0, 1.0))
    bootstrap_ari = (
        np.asarray([adjusted_rand_score(labels, run) for run in bootstrap_labels], dtype=float)
        if bootstrap_labels
        else np.asarray([], dtype=float)
    )
    mean_bootstrap_ari = float(bootstrap_ari.mean()) if len(bootstrap_ari) else float("nan")
    median_bootstrap_ari = (
        float(np.median(bootstrap_ari)) if len(bootstrap_ari) else float("nan")
    )
    min_bootstrap_ari = float(bootstrap_ari.min()) if len(bootstrap_ari) else float("nan")
    if np.isfinite(mean_bootstrap_ari):
        # Mean ARI alone can select a partition with a poor worst-case
        # bootstrap. Weight the lower tail explicitly so the selected module
        # backbone is reproducible rather than merely good on average.
        stability_score = (
            0.50 * mean_bootstrap_ari
            + 0.30 * median_bootstrap_ari
            + 0.20 * min_bootstrap_ari
        )
    else:
        stability_score = purity
    feasible = (
        settings.min_candidate_modules <= n_modules <= settings.max_candidate_modules
        and grey_fraction <= settings.max_grey_fraction
        and largest <= settings.max_largest_module_fraction
    )
    violation = (
        max(0.0, grey_fraction - settings.max_grey_fraction)
        + max(0.0, largest - settings.max_largest_module_fraction)
        + max(0, settings.min_candidate_modules - n_modules) * 0.10
        + max(0, n_modules - settings.max_candidate_modules) * 0.02
    )
    score = (
        0.22 * purity
        + 0.24 * stability_score
        + 0.19 * coverage
        + 0.15 * balance
        + 0.12 * separation
        + 0.08 * modularity_scaled
        - violation
    )
    return {
        "n_modules": n_modules,
        "grey_fraction": grey_fraction,
        "largest_module_fraction": largest,
        "mean_membership_confidence": purity,
        "mean_confidence_margin": separation,
        "module_balance_entropy": balance,
        "graph_modularity": modularity,
        "mean_bootstrap_ari": mean_bootstrap_ari,
        "median_bootstrap_ari": median_bootstrap_ari,
        "min_bootstrap_ari": min_bootstrap_ari,
        "robust_bootstrap_stability": float(stability_score),
        "selection_score": float(score),
        "constraints_satisfied": bool(feasible),
    }


def _select_final_modules(
    coassignment: np.ndarray,
    settings: ConsensusConfig,
    min_module_size: int,
    bootstrap_labels: list[np.ndarray] | None = None,
    candidate_scorer: CandidateScorer | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, pd.DataFrame]:
    if settings.clustering_method == "hierarchical":
        labels = _final_modules(coassignment, settings.min_coassignment, min_module_size)
        roles = np.where(labels == "grey", "grey", "core").astype(object)
        assignments, _ = _confidence_tables(
            coassignment,
            np.asarray([str(index) for index in range(len(labels))]),
            labels,
            pd.Series(dtype=object),
        )
        confidence = assignments["consensus_membership_confidence"].to_numpy(float)
        second = assignments["next_best_module_coassignment"].to_numpy(float)
        margin = confidence - second
        metrics = _candidate_metrics(labels, confidence, margin, 0.0, settings, bootstrap_labels)
        if candidate_scorer is not None:
            metrics.update(candidate_scorer(labels))
        search = pd.DataFrame([{"method": "hierarchical", "selected": True, **metrics}])
        search = _finalize_candidate_objectives(search, candidate_scorer is not None)
        return labels, roles, confidence, second, margin, labels.copy(), search

    source, target, weights = _top_k_edges(coassignment, settings.graph_top_k)
    candidates: list[dict[str, object]] = []
    payloads: list[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = []
    for threshold in sorted(set(float(value) for value in settings.graph_edge_thresholds)):
        for resolution in sorted(set(float(value) for value in settings.graph_resolutions)):
            raw, modularity, n_edges, method = _graph_partition(
                len(coassignment), source, target, weights, threshold, resolution, settings.random_seed
            )
            candidate = _canonical_modules(raw, min_module_size)
            labels, roles, confidence, second, margin, candidate_display = _soft_membership(
                coassignment, candidate, min_module_size, settings
            )
            metrics = _candidate_metrics(
                labels, confidence, margin, modularity, settings, bootstrap_labels
            )
            if candidate_scorer is not None:
                metrics.update(candidate_scorer(labels))
            candidates.append(
                {
                    "method": method,
                    "edge_threshold": threshold,
                    "resolution": resolution,
                    "graph_edges": n_edges,
                    "max_grey_fraction_allowed": settings.max_grey_fraction,
                    "max_largest_module_fraction_allowed": settings.max_largest_module_fraction,
                    **metrics,
                }
            )
            payloads.append((labels, roles, confidence, second, margin, candidate_display))
    search = pd.DataFrame(candidates)
    search = _finalize_candidate_objectives(search, candidate_scorer is not None)
    feasible = search.index[search["constraints_satisfied"]]
    pool = feasible if len(feasible) else search.index
    if candidate_scorer is not None:
        pareto_pool = search.index[
            search["pareto_optimal"] & search.index.isin(pool)
        ]
        selection_pool = pareto_pool if len(pareto_pool) else pool
        selected = int(
            search.loc[selection_pool, "multiobjective_selection_score"].idxmax()
        )
    else:
        selected = int(search.loc[pool, "selection_score"].idxmax())
    search["selected"] = False
    search.loc[selected, "selected"] = True
    sort_score = (
        "multiobjective_selection_score"
        if candidate_scorer is not None
        else "selection_score"
    )
    search = search.sort_values(
        ["selected", sort_score, "edge_threshold", "resolution"],
        ascending=[False, False, True, True],
    ).reset_index(drop=True)
    return (*payloads[selected], search)


def _finalize_candidate_objectives(
    search: pd.DataFrame,
    multiobjective_enabled: bool,
) -> pd.DataFrame:
    search = search.copy()
    if not multiobjective_enabled:
        search["selection_method"] = "legacy_internal_consensus"
        return search
    if (
        "n_preservation_module_cohort_tests" not in search
        or pd.to_numeric(
            search["n_preservation_module_cohort_tests"], errors="coerce"
        ).fillna(0).max()
        <= 0
    ):
        raise ValueError(
            "Multi-objective module selection has no eligible module/cohort "
            "preservation tests. Check exact event overlap and "
            "min_shared_events_per_module."
        )
    functional_mode = (
        str(search["functional_weight_mode"].iloc[0])
        if "functional_weight_mode" in search
        else "always"
    )
    functional_modules = (
        pd.to_numeric(search["n_functional_modules_scored"], errors="coerce").fillna(0)
        if "n_functional_modules_scored" in search
        else pd.Series(0.0, index=search.index)
    )
    significant_modules = (
        pd.to_numeric(
            search["n_significant_functional_modules"], errors="coerce"
        ).fillna(0)
        if "n_significant_functional_modules" in search
        else functional_modules
    )
    minimum_significant = (
        int(search["min_significant_functional_modules"].iloc[0])
        if "min_significant_functional_modules" in search
        else 1
    )
    functional_values = pd.to_numeric(
        search.get(
            "objective_functional_coherence",
            pd.Series(0.0, index=search.index),
        ),
        errors="coerce",
    ).fillna(0.0)
    has_functional_signal = bool(
        functional_modules.max() > 0
        and functional_values.max() > 0
        and significant_modules.max() >= minimum_significant
    )
    if functional_mode == "always" and not has_functional_signal:
        raise ValueError(
            "Multi-objective module selection requires functional coherence, "
            "but no candidate reached the configured functional evidence gate."
        )
    functional_active = functional_mode == "always" or (
        functional_mode == "auto_significance" and has_functional_signal
    )

    robust_stability = pd.to_numeric(
        search["robust_bootstrap_stability"], errors="coerce"
    ).clip(0.0, 1.0).fillna(0.0)
    coverage = 1.0 - pd.to_numeric(
        search.get("grey_fraction", pd.Series(0.0, index=search.index)),
        errors="coerce",
    ).clip(0.0, 1.0).fillna(1.0)
    balance = pd.to_numeric(
        search.get(
            "module_balance_entropy", pd.Series(1.0, index=search.index)
        ),
        errors="coerce",
    ).clip(0.0, 1.0).fillna(0.0)
    purity = pd.to_numeric(
        search.get(
            "mean_membership_confidence", pd.Series(1.0, index=search.index)
        ),
        errors="coerce",
    ).clip(0.0, 1.0).fillna(0.0)
    separation = pd.to_numeric(
        search.get(
            "mean_confidence_margin", pd.Series(1.0, index=search.index)
        ),
        errors="coerce",
    ).clip(0.0, 1.0).fillna(0.0)
    n_modules = pd.to_numeric(
        search.get("n_modules", pd.Series(2, index=search.index)),
        errors="coerce",
    ).fillna(0)
    partition_validity = (
        0.35 * coverage + 0.25 * balance + 0.25 * purity + 0.15 * separation
    )
    partition_validity = partition_validity.where(
        n_modules >= 2, partition_validity * 0.25
    )
    search["objective_partition_validity"] = partition_validity.clip(0.0, 1.0)
    search["objective_bootstrap_stability"] = robust_stability
    search["objective_internal_stability"] = np.sqrt(
        np.clip(
            robust_stability * search["objective_partition_validity"],
            0.0,
            1.0,
        )
    )
    objective_columns = [
        "objective_internal_stability",
        "objective_external_preservation",
    ]
    weights = [
        float(search["objective_stability_weight"].iloc[0]),
        float(search["objective_preservation_weight"].iloc[0]),
    ]
    if functional_active:
        objective_columns.append("objective_functional_coherence")
        weights.append(float(search["objective_functional_weight"].iloc[0]))
    weights = np.asarray(weights, dtype=float)
    weights = weights / weights.sum()
    search["functional_objective_active"] = functional_active
    search["effective_stability_weight"] = weights[0]
    search["effective_preservation_weight"] = weights[1]
    search["effective_functional_weight"] = (
        weights[2] if functional_active else 0.0
    )
    floor = float(search["objective_score_floor"].iloc[0])
    objective_values = (
        search[objective_columns]
        .apply(pd.to_numeric, errors="coerce")
        .fillna(0.0)
        .clip(lower=floor, upper=1.0)
        .to_numpy(dtype=float)
    )
    search["multiobjective_selection_score"] = np.exp(
        np.log(objective_values) @ weights
    )
    pareto = np.ones(len(search), dtype=bool)
    raw_values = (
        search[objective_columns]
        .apply(pd.to_numeric, errors="coerce")
        .fillna(0.0)
        .to_numpy(dtype=float)
    )
    for index in range(len(search)):
        dominates = np.all(raw_values >= raw_values[index], axis=1) & np.any(
            raw_values > raw_values[index], axis=1
        )
        if np.any(dominates):
            pareto[index] = False
    search["pareto_optimal"] = pareto
    search["selection_method"] = "pareto_weighted_geometric_mean"
    return search


def _edge_table(
    coassignment: np.ndarray,
    event_keys: np.ndarray,
    threshold: float,
    max_edges: int,
) -> pd.DataFrame:
    rows: list[tuple[str, str, float]] = []
    for left in range(len(event_keys) - 1):
        right_indices = np.flatnonzero(coassignment[left, left + 1 :] >= threshold) + left + 1
        for right in right_indices:
            rows.append((str(event_keys[left]), str(event_keys[right]), float(coassignment[left, right])))
        if len(rows) >= max_edges:
            break
    table = pd.DataFrame(rows[:max_edges], columns=["event_a", "event_b", "coassignment_probability"])
    if table.empty:
        return table
    return table.sort_values(
        ["coassignment_probability", "event_a", "event_b"], ascending=[False, True, True]
    ).reset_index(drop=True)


def build_bootstrap_coassignment_consensus(
    matrix: pd.DataFrame,
    observed_assignments: pd.DataFrame,
    samples: pd.DataFrame,
    network: NetworkConfig,
    settings: ConsensusConfig,
    gene_prior: GeneMantelPrior | None = None,
    gene_prior_weight: float = 0.0,
    candidate_scorer: CandidateScorer | None = None,
    progress: Progress = None,
) -> ConsensusResult:
    if matrix.empty:
        raise ValueError("Cannot build consensus modules from an empty event matrix.")
    if len(matrix) > settings.dense_max_events:
        raise ValueError(
            f"Consensus matrix requires {len(matrix)} events, exceeding consensus.dense_max_events="
            f"{settings.dense_max_events}. Increase the explicit memory cap or use a smaller stable-AS universe."
        )
    if not matrix.index.is_unique or not matrix.columns.is_unique:
        raise ValueError("Consensus input event and sample labels must be unique.")
    baseline = observed_assignments.set_index("event_key")["module"].reindex(matrix.index).fillna("grey").astype(str)
    event_keys = matrix.index.astype(str).to_numpy()
    n_events = len(event_keys)
    counts = np.zeros((n_events, n_events), dtype=np.uint16)
    run_rows: list[dict[str, object]] = []
    run_labels: list[np.ndarray] = []
    rng = np.random.default_rng(settings.random_seed)
    # Outer co-assignment bootstrap is the stability estimator. Disable the
    # expensive inner edge/community bootstraps to avoid nested resampling.
    partition_settings = replace(
        network,
        consensus_bootstrap_iterations=0,
        community_consensus_runs=1,
        bootstrap_iterations=0,
    )
    for iteration in range(1, settings.bootstrap_iterations + 1):
        indices = _bootstrap_indices(samples, matrix.columns, rng)
        bootstrap = matrix.iloc[:, indices].copy()
        bootstrap.columns = [
            f"coassign_{iteration:03d}_{position:03d}"
            for position in range(len(indices))
        ]
        labels = discover_module_labels(
            bootstrap,
            partition_settings,
            gene_prior=gene_prior,
            gene_prior_weight=gene_prior_weight,
        )
        run_labels.append(labels.copy())
        for module in set(labels).difference({"grey"}):
            module_indices = np.flatnonzero(labels == module)
            counts[np.ix_(module_indices, module_indices)] += 1
        non_grey_counts = pd.Series(labels[labels != "grey"]).value_counts()
        run_rows.append(
            {
                "bootstrap_iteration": iteration,
                "adjusted_rand_index_vs_observed": adjusted_rand_score(baseline, labels),
                "n_modules": int(len(non_grey_counts)),
                "n_grey_events": int(np.count_nonzero(labels == "grey")),
                "largest_module_fraction": (
                    float(non_grey_counts.max() / n_events) if not non_grey_counts.empty else 0.0
                ),
            }
        )
        if progress and (
            iteration == 1
            or iteration % 5 == 0
            or iteration == settings.bootstrap_iterations
        ):
            progress(
                f"Co-assignment bootstrap {iteration}/"
                f"{settings.bootstrap_iterations}"
            )
    coassignment = counts.astype(np.float32) / float(settings.bootstrap_iterations)
    np.fill_diagonal(coassignment, 1.0)
    min_module_size = settings.min_module_size or network.min_module_size
    (
        labels,
        membership_role,
        confidence,
        second_best,
        confidence_margin,
        candidate_labels,
        parameter_search,
    ) = _select_final_modules(
        coassignment,
        settings,
        min_module_size,
        run_labels,
        candidate_scorer=candidate_scorer,
    )
    assignments, module_summary = _confidence_tables(
        coassignment,
        event_keys,
        labels,
        baseline,
        membership_role,
        candidate_labels,
        confidence,
        second_best,
        confidence_margin,
    )
    scores = _program_scores(matrix, labels)
    bootstrap_runs = pd.DataFrame(run_rows)
    edges = _edge_table(
        coassignment,
        event_keys,
        settings.output_edge_threshold,
        settings.max_output_edges,
    )
    selected_parameters = parameter_search.loc[parameter_search["selected"]].iloc[0]
    selection_audit = (
        candidate_scorer.audit(labels)
        if candidate_scorer is not None and hasattr(candidate_scorer, "audit")
        else pd.DataFrame(
            columns=[
                "audit_type",
                "cohort",
                "module",
                "objective_contribution",
            ]
        )
    )
    diagnostics = pd.DataFrame(
        [
            {
                "n_events": n_events,
                "n_samples": matrix.shape[1],
                "bootstrap_iterations": settings.bootstrap_iterations,
                "min_coassignment": settings.min_coassignment,
                "min_module_size": min_module_size,
                "final_clustering_method": selected_parameters.get("method", settings.clustering_method),
                "selected_graph_edge_threshold": selected_parameters.get("edge_threshold", np.nan),
                "selected_graph_resolution": selected_parameters.get("resolution", np.nan),
                "selected_graph_edges": selected_parameters.get("graph_edges", np.nan),
                "parameter_constraints_satisfied": selected_parameters["constraints_satisfied"],
                "parameter_selection_score": selected_parameters["selection_score"],
                "parameter_selection_method": selected_parameters.get(
                    "selection_method", "legacy_internal_consensus"
                ),
                "multiobjective_selection_score": selected_parameters.get(
                    "multiobjective_selection_score", np.nan
                ),
                "objective_internal_stability": selected_parameters.get(
                    "objective_internal_stability", np.nan
                ),
                "objective_external_preservation": selected_parameters.get(
                    "objective_external_preservation", np.nan
                ),
                "objective_functional_coherence": selected_parameters.get(
                    "objective_functional_coherence", np.nan
                ),
                "core_membership_min": settings.core_membership_min,
                "peripheral_membership_min": settings.peripheral_membership_min,
                "peripheral_membership_margin": settings.peripheral_membership_margin,
                "n_consensus_modules": len(set(labels).difference({"grey"})),
                "n_grey_events": int(np.count_nonzero(labels == "grey")),
                "grey_fraction": float(np.mean(labels == "grey")),
                "n_core_events": int(np.count_nonzero(membership_role == "core")),
                "n_peripheral_events": int(np.count_nonzero(membership_role == "peripheral")),
                "n_uncertain_events": int(np.count_nonzero(membership_role == "uncertain")),
                "adjusted_rand_index_observed_vs_consensus": adjusted_rand_score(baseline, labels),
                "median_bootstrap_ari_vs_observed": float(bootstrap_runs["adjusted_rand_index_vs_observed"].median()),
                "coassignment_storage": "dense_npz",
                "partition_inner_bootstrap_disabled": True,
                "gene_mantel_prior_enabled": gene_prior is not None,
                "gene_mantel_prior_weight": (
                    gene_prior_weight if gene_prior is not None else 0.0
                ),
            }
        ]
    )
    return ConsensusResult(
        assignments,
        module_summary,
        scores,
        bootstrap_runs,
        edges,
        coassignment,
        event_keys,
        diagnostics,
        parameter_search,
        selection_audit,
    )


def save_coassignment_matrix(result: ConsensusResult, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, event_keys=result.event_keys, coassignment=result.coassignment_matrix)
    return path


def render_consensus_plots(result: ConsensusResult, output_dir: Path) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    generated: list[Path] = []
    assignments = result.assignments.sort_values(
        ["consensus_module", "consensus_membership_confidence"], ascending=[True, False]
    )
    visible = assignments.loc[assignments["consensus_module"].ne("grey")].head(250)
    if not visible.empty:
        positions = pd.Series(np.arange(len(result.event_keys)), index=result.event_keys)
        indices = positions.reindex(visible["event_key"]).dropna().astype(int).to_numpy()
        matrix = result.coassignment_matrix[np.ix_(indices, indices)]
        fig, ax = plt.subplots(figsize=(9, 8))
        sns.heatmap(matrix, cmap="mako", vmin=0, vmax=1, ax=ax, xticklabels=False, yticklabels=False)
        ax.set_title("Bootstrap event co-assignment consensus (up to 250 events)", fontweight="bold")
        fig.tight_layout()
        path = output_dir / "bootstrap_coassignment_heatmap.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        generated.append(path)

    fig, ax = plt.subplots(figsize=(9, 5.5))
    sns.histplot(
        result.assignments["consensus_membership_confidence"], bins=20, color="#087E8B", ax=ax
    )
    ax.set_xlim(0, 1)
    ax.set_xlabel("Mean co-assignment to final module")
    ax.set_title("Event consensus confidence", fontweight="bold")
    fig.tight_layout()
    path = output_dir / "event_consensus_confidence.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    generated.append(path)

    if not result.module_summary.empty:
        table = result.module_summary.sort_values("n_events", ascending=False)
        fig, ax = plt.subplots(figsize=(max(8, len(table) * 0.8), 5.5))
        bars = ax.bar(table["consensus_module"], table["n_events"], color="#2D6CDF", edgecolor="#17324D")
        for bar, confidence in zip(bars, table["mean_event_confidence"], strict=False):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f"{confidence:.2f}", ha="center", va="bottom")
        ax.set_ylabel("Events")
        ax.set_title("Final consensus modules (labels show mean confidence)", fontweight="bold")
        fig.tight_layout()
        path = output_dir / "consensus_module_sizes.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        generated.append(path)

    if not result.parameter_search.empty and "edge_threshold" in result.parameter_search:
        table = result.parameter_search.copy()
        color_column = (
            "multiobjective_selection_score"
            if "multiobjective_selection_score" in table
            else "selection_score"
        )
        fig, ax = plt.subplots(figsize=(9, 6))
        scatter = ax.scatter(
            table["grey_fraction"],
            table["largest_module_fraction"],
            c=table[color_column],
            s=55 + 12 * table["n_modules"].clip(upper=20),
            cmap="viridis",
            edgecolor="#17324D",
            linewidth=0.6,
        )
        selected = table.loc[table["selected"]]
        if not selected.empty:
            ax.scatter(
                selected["grey_fraction"],
                selected["largest_module_fraction"],
                marker="*",
                s=320,
                color="#D1495B",
                edgecolor="white",
                linewidth=1.2,
                label="Selected partition",
            )
            ax.legend(loc="best", frameon=False)
        ax.axvline(
            float(table["max_grey_fraction_allowed"].iloc[0]),
            color="#6B7280",
            linestyle="--",
            linewidth=1,
        )
        ax.axhline(
            float(table["max_largest_module_fraction_allowed"].iloc[0]),
            color="#6B7280",
            linestyle="--",
            linewidth=1,
        )
        ax.set_xlabel("Grey event fraction")
        ax.set_ylabel("Largest module fraction")
        ax.set_title("Consensus partition parameter search", fontweight="bold")
        fig.colorbar(
            scatter,
            ax=ax,
            label=(
                "Multi-objective selection score"
                if color_column == "multiobjective_selection_score"
                else "Selection score"
            ),
        )
        fig.tight_layout()
        path = output_dir / "consensus_parameter_search.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        generated.append(path)

    objective_columns = {
        "objective_internal_stability",
        "objective_external_preservation",
        "objective_functional_coherence",
        "multiobjective_selection_score",
    }
    if objective_columns.issubset(result.parameter_search.columns):
        table = result.parameter_search.copy()
        fig, ax = plt.subplots(figsize=(8.5, 6.5))
        sizes = 70 + 260 * pd.to_numeric(
            table["objective_internal_stability"], errors="coerce"
        ).fillna(0.0).clip(0.0, 1.0)
        pareto = table.get(
            "pareto_optimal", pd.Series(False, index=table.index)
        ).astype(bool)
        scatter = ax.scatter(
            table["objective_external_preservation"],
            table["objective_functional_coherence"],
            s=sizes,
            c=table["multiobjective_selection_score"],
            cmap="viridis",
            edgecolor=np.where(pareto, "#111827", "#9CA3AF"),
            linewidth=np.where(pareto, 1.4, 0.5),
            alpha=0.85,
        )
        selected = table.loc[table["selected"]]
        if not selected.empty:
            ax.scatter(
                selected["objective_external_preservation"],
                selected["objective_functional_coherence"],
                marker="*",
                s=420,
                color="#D1495B",
                edgecolor="white",
                linewidth=1.4,
                label="Selected partition",
                zorder=5,
            )
            ax.legend(loc="best", frameon=False)
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.set_xlabel("External preservation objective (tuning cohorts)")
        ax.set_ylabel("Functional coherence objective (tuning GMT)")
        ax.set_title(
            "Multi-objective consensus partition selection\n"
            "Point size = internal bootstrap stability",
            fontweight="bold",
        )
        fig.colorbar(
            scatter, ax=ax, label="Weighted geometric objective"
        )
        fig.tight_layout()
        path = output_dir / "multiobjective_candidate_tradeoff.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        generated.append(path)
    return generated
