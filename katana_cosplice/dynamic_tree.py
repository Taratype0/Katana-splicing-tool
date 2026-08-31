from __future__ import annotations

import numpy as np
import warnings
import dynamicTreeCut.df_apply as _df_apply
import dynamicTreeCut.dynamicTreeCut as _dynamic_tree_module
from dynamicTreeCut import cutreeHybrid
from scipy.spatial.distance import squareform


# dynamicTreeCut 0.1.1 imports this helper module but references it through an
# unqualified name in the PAM leftover-assignment branch. Keep the compatibility
# fix local to Katana rather than mutating the installed package on disk.
if not hasattr(_dynamic_tree_module, "df_apply"):
    _dynamic_tree_module.df_apply = _df_apply


def dynamic_tree_cut_labels(
    linkage_matrix: np.ndarray,
    distance_matrix: np.ndarray,
    *,
    min_module_size: int,
    deep_split: int = 2,
    cut_height: float = 0.995,
    pam_respects_dendro: bool = False,
) -> np.ndarray:
    """Run the published hybrid Dynamic Tree Cut on one hierarchy."""

    distance = np.asarray(distance_matrix, dtype=float)
    if distance.ndim != 2 or distance.shape[0] != distance.shape[1]:
        raise ValueError("Dynamic Tree Cut requires a square distance matrix.")
    if linkage_matrix.shape != (len(distance) - 1, 4):
        raise ValueError("Linkage matrix does not match the distance universe.")
    if min_module_size < 2:
        raise ValueError("Dynamic Tree Cut min_module_size must be at least 2.")
    if deep_split not in {0, 1, 2, 3, 4}:
        raise ValueError("Dynamic Tree Cut deep_split must be within 0..4.")
    if not 0.0 < cut_height <= 1.0:
        raise ValueError("Dynamic Tree Cut cut_height must be within (0, 1].")

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Bitwise inversion '~' on bool is deprecated",
            category=DeprecationWarning,
            module="dynamicTreeCut",
        )
        result = cutreeHybrid(
            linkage_matrix,
            squareform(distance, checks=False),
            cutHeight=float(cut_height),
            minClusterSize=int(min_module_size),
            deepSplit=int(deep_split),
            pamStage=True,
            pamRespectsDendro=bool(pam_respects_dendro),
            respectSmallClusters=False,
            verbose=0,
        )
    labels = np.asarray(result["labels"], dtype=int)
    if len(labels) != len(distance):
        raise RuntimeError("Dynamic Tree Cut did not preserve the event universe.")
    return labels
