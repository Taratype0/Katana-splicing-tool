from __future__ import annotations

from datetime import datetime, timezone
import sys
import time
import numpy as np
import warnings
import dynamicTreeCut.df_apply as _df_apply
import dynamicTreeCut.dynamicTreeCut as _dynamic_tree_module
from dynamicTreeCut import cutreeHybrid
from scipy.spatial.distance import squareform
from scipy.stats import rankdata as _scipy_rankdata


# dynamicTreeCut 0.1.1 imports this helper module but references it through an
# unqualified name in the PAM leftover-assignment branch. Keep the compatibility
# fix local to Katana rather than mutating the installed package on disk.
if not hasattr(_dynamic_tree_module, "df_apply"):
    _dynamic_tree_module.df_apply = _df_apply


def _rankdata_for_dynamic_tree_cut(values, method="average", *args, **kwargs):
    """Restore integer ordinal indices expected by dynamicTreeCut 0.1.1.

    Recent SciPy returns floating-point ordinal ranks. Their numeric ordering
    is unchanged, but this older package uses those ranks as NumPy indices.
    Other ranking methods must retain their original values and dtype.
    """
    ranks = _scipy_rankdata(values, method, *args, **kwargs)
    if method != "ordinal":
        return ranks
    array = np.asarray(ranks)
    if (not np.all(np.isfinite(array)) or not np.all(array == np.floor(array))
            or np.any(array < 1) or np.any(array > np.size(values))):
        raise ValueError("Dynamic Tree Cut received invalid ordinal rank indices")
    return array.astype(np.intp)


# Confine the compatibility binding to this dependency's namespace; do not
# modify scipy.stats.rankdata or any installed source files.
_dynamic_tree_module.rankdata = _rankdata_for_dynamic_tree_cut


if not hasattr(np, "in1d"):
    class _DynamicTreeNumpyCompatibility:
        """Keep the removed NumPy helper local to dynamicTreeCut 0.1.1."""

        def __getattr__(self, name):
            return getattr(np, name)

        @staticmethod
        def in1d(ar1, ar2, assume_unique=False, invert=False):
            return np.isin(
                np.asarray(ar1).ravel(), ar2,
                assume_unique=assume_unique, invert=invert,
            )

    # The legacy helper was np.isin(...).ravel() with these same defaults.
    # Do not add attributes to the process-wide NumPy module.
    _dynamic_tree_module.np = _DynamicTreeNumpyCompatibility()


# The dependency accepts square distances, but its default two-stage indexing
# materialises len(indices) * n_events values before selecting the columns.
# Its condensed-distance alternative is much worse: it reconstructs every
# pair through nested Python loops and scipy.special.binom.  Keep its exact
# ordering, duplicates and numeric values while indexing only the requested
# submatrix.  The compatibility binding remains local to this dependency.
if not hasattr(_dynamic_tree_module, "_katana_original_dist_multi_index"):
    _dynamic_tree_module._katana_original_dist_multi_index = _dynamic_tree_module.dist_multi_index
_original_dist_multi_index = _dynamic_tree_module._katana_original_dist_multi_index


def _dist_multi_index_for_dynamic_tree_cut(indices, matrix):
    if np.ndim(matrix) == 2:
        return matrix[np.ix_(indices, indices)]
    return _original_dist_multi_index(indices, matrix)


_dynamic_tree_module.dist_multi_index = _dist_multi_index_for_dynamic_tree_cut


def dynamic_tree_cut_label_grid(
    linkage_matrix: np.ndarray,
    distance_matrix: np.ndarray,
    *,
    min_module_size: int,
    deep_splits: list[int],
    cut_height: float = 0.995,
    pam_respects_dendro: bool = False,
) -> dict[int, np.ndarray]:
    """Run declared cuts on one unchanged hierarchy and one distance matrix.

    Preparing the distance once is an execution optimisation only: every cut
    still calls the same published cutreeHybrid implementation with the same
    numeric matrix, hierarchy, PAM settings and deepSplit value.
    """

    distance = np.asarray(distance_matrix, dtype=float)
    if distance.ndim != 2 or distance.shape[0] != distance.shape[1]:
        raise ValueError("Dynamic Tree Cut requires a square distance matrix.")
    if linkage_matrix.shape != (len(distance) - 1, 4):
        raise ValueError("Linkage matrix does not match the distance universe.")
    if min_module_size < 2:
        raise ValueError("Dynamic Tree Cut min_module_size must be at least 2.")
    splits = sorted(set(deep_splits))
    if not splits or any(split not in {0, 1, 2, 3, 4} for split in splits):
        raise ValueError("Dynamic Tree Cut deep_split must be within 0..4.")
    if not 0.0 < cut_height <= 1.0:
        raise ValueError("Dynamic Tree Cut cut_height must be within (0, 1].")

    # The old condensed route implicitly discarded the diagonal/lower triangle.
    # Normal network TOM distances are symmetric with a zero diagonal. Preserve
    # the old public behaviour for noncanonical caller inputs rather than
    # silently interpreting their lower triangle or diagonal differently.
    prepared = distance
    if np.any(np.diag(distance) != 0) or not np.array_equal(distance, distance.T):
        prepared = squareform(distance, checks=False)

    output: dict[int, np.ndarray] = {}
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Bitwise inversion '~' on bool is deprecated",
            category=DeprecationWarning,
            module="dynamicTreeCut",
        )
        for split in splits:
            started = time.perf_counter()
            prefix = (f"[Katana network] {datetime.now(timezone.utc).isoformat()} "
                      f"dynamic_tree_cut deep_split={split}")
            print(prefix + " started", file=sys.stderr, flush=True)
            result = cutreeHybrid(
                linkage_matrix,
                prepared,
                cutHeight=float(cut_height),
                minClusterSize=int(min_module_size),
                deepSplit=int(split),
                pamStage=True,
                pamRespectsDendro=bool(pam_respects_dendro),
                respectSmallClusters=False,
                verbose=0,
            )
            labels = np.asarray(result["labels"], dtype=int)
            if len(labels) != len(distance):
                raise RuntimeError("Dynamic Tree Cut did not preserve the event universe.")
            output[int(split)] = labels
            print(prefix + f" completed elapsed_seconds={time.perf_counter() - started:.3f}",
                  file=sys.stderr, flush=True)
    return output


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
    return dynamic_tree_cut_label_grid(
        linkage_matrix,
        distance_matrix,
        min_module_size=min_module_size,
        deep_splits=[deep_split],
        cut_height=cut_height,
        pam_respects_dendro=pam_respects_dendro,
    )[deep_split]
