import numpy as np
import pandas as pd
import pytest

from katana_next.science.motiftrace_continuous import continuous_inputs, gene_contributions, score_from_contributions, calibrated_tail


def test_shared_gene_flips_are_not_absolute_weight_randomization():
    u = np.array([[1., -1.], [-2., 2.], [3., 3.]])
    signs = np.array([[1., -1., 1.], [-1., 1., -1.]])
    actual = score_from_contributions(u, signs)
    np.testing.assert_allclose(actual[0], -actual[1])
    assert actual[0, 0] != actual[0, 1]


def test_gene_blocks_keep_multiple_events_together():
    x = np.arange(60).reshape(30, 2)
    y = np.arange(30.)
    u, genes, valid = gene_contributions(x, y, np.ones(30), [str(j//2) for j in range(30)])
    assert np.isnan(u).all()  # 15 genes, not 30 independent units.
    u, genes, valid = gene_contributions(x, y, np.ones(30), [str(j) for j in range(30)])
    assert len(genes) == 30 and np.isfinite(u).all()


def test_calibration_uses_whole_family_and_conservative_ties():
    p = calibrated_tail([0, 2, 5], [1, 2, 2, 4])
    np.testing.assert_allclose(p, [1, .8, .2])
    assert np.isnan(calibrated_tail([np.nan], [1, 2, 3])[0])


def test_no_group_identity_guess_or_missing_measurement_fill():
    e = pd.DataFrame({"event_id": ["e"], "gene_id": ["g"], "role": ["background"], "dpsi": [0], "fdr": [1],
        "IncLevel1": ["NA,0.8"], "IncLevel2": ["0.2,0.3"]})
    t = pd.DataFrame([[1, 2, 3, 4, 10, 10, 10, 10]], index=["e"])
    with pytest.raises(ValueError):
        continuous_inputs(e, t, source_group1_is_numerator=False)
    result = continuous_inputs(e, t, source_group1_is_numerator=True)
    assert not result.PSI_eligible.iloc[0]
    assert np.isnan(result.delta_PSI_KD_minus_Ctrl.iloc[0])
    assert result.expression_weight.iloc[0] == .5


def test_primitive_arm_swap_equals_shared_signed_gene_null_and_keeps_weights():
    n = 24
    e = pd.DataFrame({"event_id": [f"e{i}" for i in range(n)], "gene_id": [f"g{i}" for i in range(n)],
        "role": ["background"]*n, "dpsi": [0.]*n, "fdr": [1.]*n,
        "IncLevel1": ["0.7,0.8"]*n, "IncLevel2": ["0.2,0.3"]*n})
    totals = np.arange(4*n, dtype=float).reshape(n, 4)+10
    fractions = np.tile([.2, .3, .7, .8], (n, 1))
    table = pd.DataFrame(np.c_[totals*fractions, totals], index=e.event_id)
    original = continuous_inputs(e, table, source_group1_is_numerator=True)
    swapped, altered = table.copy(), e.copy()
    mask = np.arange(n) % 2 == 0
    swapped.loc[mask] = table.loc[mask].to_numpy()[:, [2, 3, 0, 1, 6, 7, 4, 5]]
    altered.loc[mask, "IncLevel1"] = "0.2,0.3"
    altered.loc[mask, "IncLevel2"] = "0.7,0.8"
    null = continuous_inputs(altered, swapped, source_group1_is_numerator=True)
    signs = np.where(mask, -1., 1.)
    np.testing.assert_allclose(original.expression_weight, null.expression_weight)
    for column in ("delta_PSI_KD_minus_Ctrl", "delta_transcript_fraction_KD_minus_Ctrl"):
        np.testing.assert_allclose(null[column], original[column]*signs, atol=1e-15)
