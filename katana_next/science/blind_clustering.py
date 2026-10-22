"""Cluster pre-UMAP features without annotations; reveal labels separately."""
import re

import numpy as np
import pandas as pd
from sklearn.metrics import (adjusted_mutual_info_score, adjusted_rand_score,
                             normalized_mutual_info_score, silhouette_score)
from sklearn.neighbors import NearestNeighbors


def cluster_latent(latent, *, n_neighbors=15, resolution=1., random_seed=17):
    columns = sorted((c for c in latent if re.fullmatch(r"pc_\d+", c)),
                     key=lambda c: int(c.split("_")[1]))
    if "cell_id" not in latent or not columns or latent.cell_id.isna().any() or latent.cell_id.duplicated().any():
        raise ValueError("Unique cell IDs and pre-UMAP pc_N columns are required")
    if not np.isfinite(resolution) or resolution <= 0 or n_neighbors < 2 or random_seed < 0:
        raise ValueError("Resolution must be positive; neighbors >= 2 and seed >= 0")
    # Only numeric features and identity enter the clustering, never annotations.
    frame = latent[["cell_id", *columns]].copy()
    frame["cell_id"] = frame.cell_id.astype(str)
    if frame.cell_id.str.strip().eq("").any() or frame.cell_id.duplicated().any():
        raise ValueError("Cell identities must be nonempty and unique after conversion")
    frame = frame.sort_values("cell_id").reset_index(drop=True)
    x = frame[columns].to_numpy(dtype=float)
    if len(x) < 5 or not np.isfinite(x).all() or not np.any(np.var(x, axis=0) > 0):
        raise ValueError("At least five cells with finite, variable pre-UMAP features are required")
    k = min(int(n_neighbors), len(x) - 1)
    distances, neighbors = NearestNeighbors(n_neighbors=k + 1).fit(x).kneighbors(x)
    pairs = [([(int(j), float(d)) for j, d in zip(ids, ds) if j != i])[:k]
             for i, (ids, ds) in enumerate(zip(neighbors, distances))]
    scales = np.array([max(p[-1][1], np.finfo(float).eps) for p in pairs])
    edges = {}
    for i, entries in enumerate(pairs):
        for j, distance in entries:
            key = (min(i, j), max(i, j))
            weight = float(np.exp(-distance**2 / (scales[i] * scales[j])))
            edges[key] = max(edges.get(key, 0.), weight)
    import igraph as ig
    import leidenalg
    graph = ig.Graph(n=len(x), edges=list(edges), directed=False)
    graph.es["weight"] = list(edges.values())
    partition = leidenalg.find_partition(graph, leidenalg.RBConfigurationVertexPartition,
        weights="weight", resolution_parameter=float(resolution), seed=int(random_seed), n_iterations=-1)
    memberships = np.array(partition.membership)
    ordered = sorted(np.unique(memberships), key=lambda group: frame.loc[memberships == group, "cell_id"].min())
    names = {group: f"C{i + 1:02d}" for i, group in enumerate(ordered)}
    assignments = pd.DataFrame({"cell_id": frame.cell_id, "cluster_id": [names[g] for g in memberships]})
    assignments["row_id"] = assignments.cell_id
    metrics = {"cells": len(x), "features": len(columns), "neighbors": k,
        "resolution": float(resolution), "seed": int(random_seed), "clusters": len(ordered),
        "graph_components": len(graph.connected_components()), "edges": len(edges),
        "min_cluster_cells": int(assignments.cluster_id.value_counts().min()),
        "max_cluster_fraction": float(assignments.cluster_id.value_counts().max() / len(x)),
        "silhouette_PCA": float(silhouette_score(x, memberships)) if 1 < len(ordered) < len(x) else np.nan,
        "method": "symmetric union kNN, local-scaled RBF weights, Leiden RBConfiguration",
        "labels_used": False, "clustering_space": "pre-UMAP PCA; not the two-dimensional map"}
    return assignments, metrics


def reveal_partition(assignments, annotations, *, label_column="cell_type"):
    if not {"cell_id", label_column} <= set(annotations) or annotations.cell_id.duplicated().any():
        raise ValueError("Unique evaluation annotations are required")
    if assignments.cell_id.duplicated().any() or set(assignments.cell_id) != set(annotations.cell_id):
        raise ValueError("Reveal must use exactly the frozen cells, without silent intersection")
    joined = assignments.merge(annotations[["cell_id", label_column]], on="cell_id", validate="one_to_one")
    if joined[label_column].isna().any() or joined[label_column].astype(str).str.strip().eq("").any():
        raise ValueError("Evaluation labels cannot be missing")
    counts = joined.groupby(["cluster_id", label_column]).size().rename("cells").reset_index()
    counts["fraction_of_cluster"] = counts.cells / counts.groupby("cluster_id").cells.transform("sum")
    counts["row_id"] = counts.cluster_id + "::" + counts[label_column].astype(str)
    metrics = {"cells": len(joined), "clusters": joined.cluster_id.nunique(),
        "label_classes": joined[label_column].nunique(),
        "adjusted_rand_index": float(adjusted_rand_score(joined[label_column], joined.cluster_id)),
        "adjusted_mutual_information": float(adjusted_mutual_info_score(joined[label_column], joined.cluster_id)),
        "normalized_mutual_information": float(normalized_mutual_info_score(joined[label_column], joined.cluster_id)),
        "cluster_purity": float(counts.groupby("cluster_id").cells.max().sum() / len(joined)),
        "label_column": label_column, "interpretation": "Post-freeze annotation agreement, not biological accuracy or donor-level inference"}
    return counts, metrics
