"""Exact scientific regressions for square-distance Dynamic Tree Cut execution."""
from dataclasses import replace
import importlib

import numpy as np
import pandas as pd
import pytest
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist, squareform

import katana_cosplice.dynamic_tree as dt
import katana_cosplice.network as network
from katana_cosplice.config import NetworkConfig


def _fixture(kind):
    rng = np.random.default_rng(61024)
    values = np.vstack([
        rng.normal(0, .04, (24, 6)),
        rng.normal(1, .05, (24, 6)),
        rng.normal(-1, .06, (24, 6)),
        rng.normal(0, .8, (12, 6)),
    ])
    if kind == "ties":
        values[1:8] = values[0]
        values = np.round(values, 1)
    elif kind == "grey":
        values[-20:] = rng.uniform(-2, 2, (20, 6))
    condensed = pdist(values)
    condensed /= condensed.max()
    return linkage(condensed, method="average"), squareform(condensed)


def _legacy_labels(tree, distance, *, min_module_size, deep_split=2,
                   cut_height=.995, pam_respects_dendro=False):
    # Frozen pre-optimisation call, including the exact upper-triangle
    # conversion and all published algorithm/PAM arguments.
    return np.asarray(dt.cutreeHybrid(
        tree, squareform(np.asarray(distance, dtype=float), checks=False),
        cutHeight=float(cut_height), minClusterSize=int(min_module_size),
        deepSplit=int(deep_split), pamStage=True,
        pamRespectsDendro=bool(pam_respects_dendro),
        respectSmallClusters=False, verbose=0,
    )["labels"], dtype=int)


def _legacy_grid(tree, distance, *, deep_splits, **kwargs):
    return {int(split): _legacy_labels(tree, distance, deep_split=split, **kwargs)
            for split in sorted(set(deep_splits))}


@pytest.mark.parametrize("kind", ["clustered", "ties", "grey"])
@pytest.mark.parametrize("pam", [False, True])
@pytest.mark.parametrize("deep_split", range(5))
def test_labels_exactly_equal_frozen_condensed_algorithm(kind, pam, deep_split):
    tree, distance = _fixture(kind)
    before = distance.copy()
    expected = _legacy_labels(tree, distance, min_module_size=8,
                              deep_split=deep_split, pam_respects_dendro=pam)
    actual = dt.dynamic_tree_cut_labels(tree, distance, min_module_size=8,
                                       deep_split=deep_split, pam_respects_dendro=pam)
    np.testing.assert_array_equal(actual, expected)
    np.testing.assert_array_equal(distance, before)


@pytest.mark.parametrize("kind", ["lower", "diagonal"])
def test_noncanonical_distances_keep_legacy_upper_triangle_semantics(kind):
    tree, distance = _fixture("grey")
    if kind == "lower":
        distance[np.tril_indices(len(distance), -1)] = 999
    else:
        np.fill_diagonal(distance, 999)
    expected = _legacy_labels(tree, distance, min_module_size=8)
    actual = dt.dynamic_tree_cut_labels(tree, distance, min_module_size=8)
    np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize("indices", [
    np.array([7, 1, 7, 4, 0]),
    np.array([False, True, False, True, False, False, True, False]),
    np.array([], dtype=int),
])
def test_square_helper_preserves_duplicate_boolean_and_empty_indices(indices):
    _, distance = _fixture("clustered")
    distance = distance[:8, :8]
    expected = distance[indices, :][:, indices]
    actual = dt._dist_multi_index_for_dynamic_tree_cut(indices, distance)
    np.testing.assert_array_equal(actual, expected)


def test_grid_keeps_every_declared_depth_and_prepares_distance_only_once(monkeypatch):
    tree, distance = _fixture("grey")
    distance = distance.astype(np.float32)
    seen = []
    original = dt.cutreeHybrid

    def observe(tree, prepared, **kwargs):
        seen.append((kwargs["deepSplit"], id(prepared), prepared.ndim, prepared.dtype))
        return original(tree, prepared, **kwargs)

    monkeypatch.setattr(dt, "cutreeHybrid", observe)
    result = dt.dynamic_tree_cut_label_grid(tree, distance, min_module_size=8,
                                          deep_splits=[4, 2, 1, 3, 2])
    assert list(result) == [1, 2, 3, 4]
    assert len(set(item[1] for item in seen)) == 1
    assert all(item[2:] == (2, np.dtype(float)) for item in seen)
    for split, labels in result.items():
        np.testing.assert_array_equal(labels, _legacy_labels(
            tree, distance, min_module_size=8, deep_split=split))


def test_import_reload_does_not_recurse_on_legacy_condensed_fallback():
    module = importlib.reload(dt)
    _, distance = _fixture("clustered")
    indices = np.array([3, 1, 5])
    result = module._dist_multi_index_for_dynamic_tree_cut(indices, squareform(distance))
    np.testing.assert_array_equal(result, distance[np.ix_(indices, indices)])


@pytest.mark.parametrize("pam", [False, True])
def test_full_dynamic_fit_equal_frozen_condensed_execution(monkeypatch, pam):
    rng = np.random.default_rng(449)
    patterns = rng.uniform(.15, .85, (3, 6))
    values = np.vstack([np.clip(pattern + rng.normal(0, .025, (24, 6)), 0, 1)
                        for pattern in patterns])
    matrix = pd.DataFrame(values, index=[f"event{i}" for i in range(len(values))],
                          columns=[f"sample{i}" for i in range(6)])
    qc = pd.DataFrame({"stability_score": np.linspace(.9, .2, len(matrix))},
                      index=matrix.index)
    settings = replace(NetworkConfig(), minimum_network_samples=6,
                       module_detection="dynamic_tree_cut", min_module_size=8,
                       dynamic_pam_respects_dendro=pam, bootstrap_iterations=0,
                       consensus_bootstrap_iterations=0)
    actual = network.fit_network(matrix, qc, settings)
    monkeypatch.setattr(network, "dynamic_tree_cut_labels", _legacy_labels)
    monkeypatch.setattr(network, "dynamic_tree_cut_label_grid", _legacy_grid)
    expected = network.fit_network(matrix, qc, settings)
    for name in ["assignments", "edges", "program_scores", "diagnostics",
                 "cut_height_sensitivity", "selected_matrix"]:
        pd.testing.assert_frame_equal(getattr(actual, name), getattr(expected, name),
                                      check_exact=True)
    np.testing.assert_array_equal(actual.linkage_matrix, expected.linkage_matrix)
