from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from igraph import Graph


@dataclass
class MultiViewFusionResult:
    fused_edges: pd.DataFrame
    view_audit: pd.DataFrame
    assignments: pd.DataFrame


def _numeric_similarity(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    valid = np.isfinite(left) & np.isfinite(right)
    distance = np.abs(np.nan_to_num(left - right, nan=0.0))
    scale = np.nanmedian(distance[valid]) if valid.any() else 1.0
    scale = max(float(scale), 1e-6)
    return np.where(valid, np.exp(-distance / scale), np.nan)


def fuse_event_edge_views(
    as_edges: pd.DataFrame,
    event_features: pd.DataFrame,
    *,
    auxiliary_columns: list[str],
    as_weight: float = 0.80,
    auxiliary_weights: dict[str, float] | None = None,
    edge_threshold: float = 0.10,
    min_module_size: int = 20,
    resolution: float = 1.0,
) -> MultiViewFusionResult:
    """Fuse event-level auxiliary views while keeping AS adjacency dominant.

    Auxiliary data reweights existing AS edges only; it cannot create an edge
    unsupported by the splice-usage network.
    """
    required = {"event_a", "event_b", "adjacency"}
    if not required.issubset(as_edges.columns):
        raise ValueError(f"AS edge table is missing columns: {sorted(required - set(as_edges))}")
    if "event_key" not in event_features:
        raise ValueError("event_features requires event_key.")
    if not 0.5 <= as_weight <= 1.0:
        raise ValueError("as_weight must be within [0.5, 1.0].")
    missing = set(auxiliary_columns).difference(event_features.columns)
    if missing:
        raise ValueError(f"Missing auxiliary feature columns: {sorted(missing)}")
    features = event_features.drop_duplicates("event_key").set_index("event_key")
    edges = as_edges.copy()
    raw_weights = auxiliary_weights or {column: 1.0 for column in auxiliary_columns}
    positive = {column: max(0.0, float(raw_weights.get(column, 0.0))) for column in auxiliary_columns}
    total = sum(positive.values())
    normalized = (
        {column: value / total for column, value in positive.items()}
        if total > 0
        else {column: 0.0 for column in auxiliary_columns}
    )
    aggregate = np.zeros(len(edges), dtype=float)
    observed_weight = np.zeros(len(edges), dtype=float)
    audit: list[dict[str, object]] = []
    for column in auxiliary_columns:
        left = pd.to_numeric(edges["event_a"].map(features[column]), errors="coerce").to_numpy()
        right = pd.to_numeric(edges["event_b"].map(features[column]), errors="coerce").to_numpy()
        similarity = _numeric_similarity(left, right)
        weight = normalized[column]
        valid = np.isfinite(similarity) & (weight > 0)
        aggregate[valid] += weight * similarity[valid]
        observed_weight[valid] += weight
        edges[f"{column}_similarity"] = similarity
        audit.append(
            {
                "view": column,
                "configured_weight": raw_weights.get(column, 0.0),
                "normalized_auxiliary_weight": weight,
                "n_edges_with_view": int(valid.sum()),
                "coverage_fraction": float(valid.mean()) if len(valid) else 0.0,
            }
        )
    auxiliary = np.divide(
        aggregate,
        observed_weight,
        out=np.zeros_like(aggregate),
        where=observed_weight > 0,
    )
    base = pd.to_numeric(edges["adjacency"], errors="coerce").fillna(0.0).to_numpy()
    effective_aux_weight = (1.0 - as_weight) * np.clip(observed_weight, 0.0, 1.0)
    edges["auxiliary_consensus_similarity"] = auxiliary
    edges["auxiliary_view_coverage"] = observed_weight
    edges["AS_primary_weight"] = 1.0 - effective_aux_weight
    edges["fused_adjacency"] = (
        (1.0 - effective_aux_weight) * base + effective_aux_weight * auxiliary
    )
    edges["fusion_rule"] = "AS_edge_reweight_only_no_auxiliary_edge_creation"
    edges["retained_for_clustering"] = edges["fused_adjacency"].ge(edge_threshold)
    retained = edges.loc[edges["retained_for_clustering"]].copy()
    event_keys = sorted(
        set(retained["event_a"].astype(str)).union(retained["event_b"].astype(str))
    )
    event_index = {event: index for index, event in enumerate(event_keys)}
    graph = Graph(n=len(event_keys), directed=False)
    if not retained.empty:
        graph.add_edges(
            [
                (event_index[str(left)], event_index[str(right)])
                for left, right in zip(
                    retained["event_a"], retained["event_b"], strict=True
                )
            ]
        )
        graph.es["weight"] = retained["fused_adjacency"].to_list()
        partition = graph.community_leiden(
            objective_function="modularity",
            weights="weight",
            resolution=resolution,
        )
        membership = np.asarray(partition.membership, dtype=int)
    else:
        membership = np.arange(len(event_keys), dtype=int)
    counts = pd.Series(membership).value_counts()
    valid_labels = [
        label for label, count in counts.sort_values(ascending=False).items()
        if count >= min_module_size
    ]
    label_map = {label: f"MV{index + 1:02d}" for index, label in enumerate(valid_labels)}
    assignments = pd.DataFrame(
        {
            "event_key": event_keys,
            "multiview_module": [
                label_map.get(label, "grey") for label in membership
            ],
        }
    )
    if not assignments.empty:
        assignments["multiview_connectivity"] = assignments["event_key"].map(
            pd.concat(
                [
                    retained[["event_a", "fused_adjacency"]].rename(
                        columns={"event_a": "event_key"}
                    ),
                    retained[["event_b", "fused_adjacency"]].rename(
                        columns={"event_b": "event_key"}
                    ),
                ]
            )
            .groupby("event_key", observed=True)["fused_adjacency"]
            .sum()
        ).fillna(0.0)
    return MultiViewFusionResult(
        edges.sort_values("fused_adjacency", ascending=False).reset_index(drop=True),
        pd.DataFrame(audit),
        assignments,
    )
