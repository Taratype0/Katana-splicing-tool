"""Numerical and export regressions for full-event network execution."""
from dataclasses import replace
import numpy as np
import pandas as pd
import pytest
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

from katana_cosplice.config import NetworkConfig
from katana_cosplice.network import (
    _bounded_edge_indices, _tom_hierarchy, correlation_matrix,
    signed_soft_adjacency, topological_overlap_matrix,
)


@pytest.mark.parametrize('method', ['bicor', 'pearson', 'spearman'])
@pytest.mark.parametrize('weighted', [False, True])
def test_full_tom_and_linkage_equal_prechange_math(method, weighted):
    rng = np.random.default_rng(291)
    matrix = pd.DataFrame(rng.uniform(size=(96, 6)), index=[f'e{i}' for i in range(96)])
    settings = replace(NetworkConfig(), method=method)
    reliability = rng.uniform(.1, 1, 96) if weighted else None
    cor = correlation_matrix(matrix, method)
    cor = np.clip((cor + cor.T) / 2, -1, 1)
    np.fill_diagonal(cor, 1)
    adj = signed_soft_adjacency(cor, settings)
    if weighted:
        adj *= np.power(np.outer(reliability, reliability), .25).astype(np.float32)
    base = adj.copy()
    tom = topological_overlap_matrix(adj)
    distance = np.clip(1 - tom, 0, 1)
    np.fill_diagonal(distance, 0)
    tree = linkage(squareform(distance, checks=False), method='average')
    actual = _tom_hierarchy(matrix, settings, event_reliability=reliability)
    for observed, expected in zip(actual[:5], [cor, adj, tom, tree, base]):
        np.testing.assert_array_equal(observed, expected)


@pytest.mark.parametrize('cap', [1, 11, 1000])
@pytest.mark.parametrize('threshold', [0., .7, 2.])
def test_streamed_export_preserves_exact_largest_edge_set(cap, threshold):
    rng = np.random.default_rng(1)
    matrix = rng.uniform(size=(41, 41))
    matrix = (matrix + matrix.T) / 2
    rows, cols = np.triu_indices(41, 1)
    chosen = [(int(r), int(c)) for r, c in zip(rows, cols) if matrix[r,c] >= threshold]
    expected = sorted(chosen, key=lambda rc: (-matrix[rc], rc))[:cap]
    for block in [1, 7, 64]:
        r, c = _bounded_edge_indices(matrix, threshold, cap, block_rows=block)
        assert set(zip(r,c)) == set(expected)


def test_export_ties_are_deterministic_without_dropping_network_events():
    matrix = np.ones((55, 55), dtype=np.float32)
    expected = {(0, i) for i in range(1, 24)}
    for block in [1, 9, 64]:
        r,c = _bounded_edge_indices(matrix, .1, 23, block_rows=block)
        assert set(zip(r,c)) == expected
    np.testing.assert_array_equal(matrix, np.ones((55,55)))


def test_legacy_event_cap_does_not_truncate_or_reject_the_universe():
    matrix = pd.DataFrame(np.random.default_rng(17).uniform(size=(41, 6)))
    actual = _tom_hierarchy(matrix, replace(NetworkConfig(), tom_dense_max_events=1))
    reference = _tom_hierarchy(matrix, NetworkConfig())
    for observed, expected in zip(actual[:5], reference[:5]):
        np.testing.assert_array_equal(observed, expected)
    assert actual[0].shape == (41, 41)
