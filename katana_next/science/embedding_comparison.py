"""Cell-matched neighborhood agreement, not a test of biological equivalence."""
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import spearmanr


def _neighbors(values, k):
    distances = squareform(pdist(values))
    np.fill_diagonal(distances, np.inf)
    return np.argsort(distances, axis=1, kind="stable")[:, :k]


def _overlap(a, b):
    return np.array([len(set(x) & set(y)) for x, y in zip(a, b)])


def compare_embeddings(expression, splicing, *, k=15, random_seed=17, permutations=200):
    """Align by explicit cell IDs; compare rotation-invariant distances and KNNs.

    The permutation reference relabels cells on one fixed neighbor graph. It is
    a descriptive chance reference, not a donor-level biological hypothesis test.
    """
    for frame in (expression, splicing):
        if "cell_id" not in frame or frame.cell_id.isna().any() or frame.cell_id.duplicated().any():
            raise ValueError("Unique nonmissing cell IDs are required")
    cells = sorted(set(expression.cell_id) & set(splicing.cell_id))
    if k < 1 or len(cells) <= k or permutations < 1:
        raise ValueError("Need more matched cells than neighbors and at least one reference permutation")
    a = expression.set_index("cell_id").loc[cells]
    b = splicing.set_index("cell_id").loc[cells]
    for label in ("cell_type", "donor_id", "condition"):
        if label in a and label in b and not a[label].equals(b[label]):
            raise ValueError(f"Paired cell metadata disagree: {label}")
    def coordinates(frame):
        columns = [c for c in frame if c.startswith("pc_")]
        if not columns:
            columns = ["umap_1", "umap_2"]
        values = frame[columns].to_numpy(dtype=float)
        if not np.isfinite(values).all():
            raise ValueError("Coordinates must be finite")
        return values
    x, y = coordinates(a), coordinates(b)
    nx, ny = _neighbors(x, k), _neighbors(y, k)
    intersections = _overlap(nx, ny)
    per_cell = pd.DataFrame({"cell_id": cells, "shared_neighbors": intersections,
        "neighbor_overlap_fraction": intersections / k,
        "neighbor_jaccard": intersections / (2 * k - intersections)})
    for label in ("cell_type", "donor_id"):
        if label in a:
            per_cell[label] = a[label].to_numpy()
    rng = np.random.default_rng(random_seed)
    chance = []
    for _ in range(permutations):
        permutation = rng.permutation(len(cells))
        relabeled = np.empty_like(ny)
        relabeled[permutation] = permutation[ny]
        chance.append(float(_overlap(nx, relabeled).mean() / k))
    summary = {"matched_cells": len(cells), "expression_only_cells": len(expression) - len(cells),
        "splicing_only_cells": len(splicing) - len(cells), "neighbors": k,
        "mean_neighbor_overlap_fraction": float(per_cell.neighbor_overlap_fraction.mean()),
        "mean_neighbor_jaccard": float(per_cell.neighbor_jaccard.mean()),
        "chance_overlap_expected": k / (len(cells) - 1),
        "chance_overlap_permutation_mean": float(np.mean(chance)),
        "chance_overlap_permutation_95_interval": np.quantile(chance, [.025, .975]).tolist(),
        "pairwise_distance_spearman": float(spearmanr(pdist(x), pdist(y)).statistic),
        "random_seed": random_seed, "reference_permutations": permutations,
        "interpretation": "Descriptive geometry agreement; no significance or causal inference from cell pairs"}
    for label in ("cell_type", "donor_id"):
        if label in a:
            labels = a[label].to_numpy()
            summary[f"expression_same_{label}_neighbor_fraction"] = float((labels[nx] == labels[:, None]).mean())
            summary[f"splicing_same_{label}_neighbor_fraction"] = float((labels[ny] == labels[:, None]).mean())
    return summary, per_cell
