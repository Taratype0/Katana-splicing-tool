from __future__ import annotations

import numpy as np
import pandas as pd

from katana_cosplice.shared_control import (
    blocked_shared_control_permutation,
    fit_shared_control_beta_binomial,
    welch_equivalence_test,
)


def test_blocked_shared_control_permutation_enumerates_all_assignments() -> None:
    psi = pd.DataFrame(
        [
            [0.2, 0.5, 0.48, 0.21, 0.52, 0.49, 0.19, 0.51, 0.50],
            [0.7, 0.4, 0.42, 0.72, 0.39, 0.43, 0.69, 0.41, 0.40],
            [0.4, 0.42, 0.39, 0.41, 0.43, 0.40, 0.39, 0.44, 0.38],
        ]
    )
    result = blocked_shared_control_permutation(psi, n_mut=2, n_drug=2)

    assert len(result.null_distribution) == 216
    assert 0 <= result.observed_rank_overlap <= 2
    assert 0 < result.empirical_overlap_pvalue <= 1


def test_welch_equivalence_requires_confidence_interval_inside_margin() -> None:
    equivalent = welch_equivalence_test(
        np.asarray([0.50, 0.51, 0.49]),
        np.asarray([0.52, 0.50, 0.51]),
        margin=0.10,
    )
    different = welch_equivalence_test(
        np.asarray([0.75, 0.76, 0.74]),
        np.asarray([0.30, 0.31, 0.29]),
        margin=0.10,
    )

    assert equivalent["equivalence_pvalue"] < 0.05
    assert equivalent["practical_similarity"] is True
    assert different["equivalence_pvalue"] > 0.05
    assert different["practical_similarity"] is False


def test_joint_beta_binomial_distinguishes_direct_group_difference() -> None:
    ijc = np.asarray(
        [
            20,
            22,
            21,
            80,
            82,
            78,
            76,
            79,
            77,
        ],
        dtype=float,
    )
    sjc = 100.0 - ijc
    psi = ijc / 100.0

    result = fit_shared_control_beta_binomial(
        ijc,
        sjc,
        psi,
        equivalence_margin=0.10,
    )

    assert result["joint_model_status"] in {"converged", "optimizer_warning"}
    assert result["joint_same_direction"] is True
    assert result["direct_MUT_vs_inhibitor_pvalue"] > 0.05
    assert result["practical_similarity"] is True
