from __future__ import annotations

import numpy as np
import pandas as pd

from katana_next.science.single_cell import donor_celltype_pseudobulk
from katana_next.science.statistics import (
    beta_binomial_group_test,
    dirichlet_multinomial_two_group,
)


def test_pseudobulk_collapses_samples_to_the_independent_donor() -> None:
    metadata = pd.DataFrame(
        [
            {
                "cell_id": "C1",
                "sample_id": "library-a",
                "donor_id": "D1",
                "cell_type": "Neuron",
                "condition": "WT",
                "batch": "B1",
                "assay": "full_length",
            },
            {
                "cell_id": "C2",
                "sample_id": "library-b",
                "donor_id": "D1",
                "cell_type": "Neuron",
                "condition": "WT",
                "batch": "B1",
                "assay": "full_length",
            },
        ]
    )
    counts = pd.DataFrame(
        [
            {
                "cell_id": cell,
                "event_id": "J1",
                "cluster_id": "C1",
                "gene_id": "G1",
                "gene_symbol": "Gene1",
                "event_type": "junction_choice",
                "annotated": True,
                "inclusion": inclusion,
                "exclusion": exclusion,
                "total": inclusion + exclusion,
            }
            for cell, inclusion, exclusion in (("C1", 3, 2), ("C2", 4, 1))
        ]
    )
    pseudobulk, _coverage = donor_celltype_pseudobulk(metadata, counts)
    assert len(pseudobulk) == 1
    assert pseudobulk.loc[0, "donor_id"] == "D1"
    assert pseudobulk.loc[0, "sample_id"] == "library-a;library-b"
    assert pseudobulk.loc[0, "inclusion"] == 7
    assert pseudobulk.loc[0, "exclusion"] == 3


def test_batch_adjusted_beta_binomial_detects_condition_effect() -> None:
    successes = np.array([8, 9, 7, 8, 42, 44, 41, 43], dtype=float)
    totals = np.repeat(50.0, len(successes))
    groups = np.array([0, 0, 0, 0, 1, 1, 1, 1])
    batches = np.array([[0], [1], [0], [1], [0], [1], [0], [1]], dtype=float)
    result = beta_binomial_group_test(successes, totals, groups, covariates=batches)
    assert result["status"] == "fit"
    assert float(result["dpsi"]) > 0.5
    assert float(result["pvalue"]) < 0.01


def test_dirichlet_multinomial_detects_multi_junction_usage_shift() -> None:
    counts = np.array(
        [
            [85, 10, 5],
            [82, 12, 6],
            [88, 8, 4],
            [80, 15, 5],
            [12, 80, 8],
            [10, 84, 6],
            [15, 78, 7],
            [8, 86, 6],
        ],
        dtype=float,
    )
    groups = np.array([0, 0, 0, 0, 1, 1, 1, 1])
    batches = np.array([[0], [1], [0], [1], [0], [1], [0], [1]], dtype=float)
    result = dirichlet_multinomial_two_group(counts, groups, covariates=batches)
    assert result["status"] == "fit"
    assert float(result["pvalue"]) < 0.01
    delta = np.asarray(result["delta_proportions"], dtype=float)
    assert delta[0] < -0.5
    assert delta[1] > 0.5
