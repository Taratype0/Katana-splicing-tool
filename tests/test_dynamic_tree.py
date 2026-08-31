from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from pathlib import Path
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from sklearn.metrics import adjusted_rand_score

from katana_cosplice.benchmark import _rscript_with_wgcna
from katana_cosplice.dynamic_tree import dynamic_tree_cut_labels
from katana_cosplice.mantel_cosplicex import (
    MantelCoSplicExConfig,
    cluster_published_cosplicex,
)


def _synthetic_adjacency() -> pd.DataFrame:
    rng = np.random.default_rng(4)
    samples = 30
    latent = [rng.normal(size=samples) for _ in range(5)]
    values = np.vstack(
        [
            signal + rng.normal(0.0, 0.12, samples)
            for signal, size in zip(latent, [18, 16, 15, 14, 12], strict=True)
            for _ in range(size)
        ]
    )
    correlation = np.corrcoef(values)
    adjacency = np.power(np.abs(correlation), 2.0)
    np.fill_diagonal(adjacency, 1.0)
    names = [f"event_{index:03d}" for index in range(len(adjacency))]
    return pd.DataFrame(adjacency, index=names, columns=names)


def _python_labels(adjacency: pd.DataFrame) -> np.ndarray:
    distance = 1.0 - adjacency.to_numpy(dtype=float)
    np.fill_diagonal(distance, 0.0)
    tree = linkage(squareform(distance, checks=False), method="average")
    return dynamic_tree_cut_labels(
        tree,
        distance,
        min_module_size=10,
        deep_split=2,
        cut_height=0.995,
    )


def test_dynamic_tree_cut_recovers_nontrivial_modules() -> None:
    labels = _python_labels(_synthetic_adjacency())
    non_grey = labels[labels > 0]

    assert len(labels) == 75
    assert len(non_grey) == 75
    assert len(np.unique(non_grey)) == 5
    assert sorted(np.unique(non_grey, return_counts=True)[1].tolist()) == [
        12,
        14,
        15,
        16,
        18,
    ]


def test_python_dynamic_tree_cut_matches_exact_r_partition() -> None:
    rscript, diagnostic = _rscript_with_wgcna()
    local_rscript = Path(
        "D:/anaconda3/envs/katana-splicing-tool/Scripts/Rscript.exe"
    )
    if rscript is None and local_rscript.exists():
        rscript = str(local_rscript)
    if rscript is None:
        pytest.skip(diagnostic)

    adjacency = _synthetic_adjacency()
    python_labels = _python_labels(adjacency)
    r_labels = cluster_published_cosplicex(
        adjacency,
        MantelCoSplicExConfig(
            min_module_size=10,
            deep_split=2,
            cut_height=0.995,
            max_core_scatter=0.82,
        ),
        rscript=rscript,
    )

    assert adjusted_rand_score(python_labels, r_labels.to_numpy()) == pytest.approx(
        1.0
    )
