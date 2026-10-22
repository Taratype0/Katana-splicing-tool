import numpy as np
import pandas as pd
import pytest

from katana_next.science.full_frozen_replication import prepare_full_frozen_families, strict_frozen_scores


def inputs():
    keys = [f"e{i}" for i in range(25)]
    a_samples = pd.DataFrame({"sample_id": [f"a{i}" for i in range(18)],
        "patient_id": [f"A{i}" for i in range(18)], "condition": "tumor"})
    b_samples = pd.DataFrame({"sample_id": [f"b{i}" for i in range(18)],
        "patient_id": [f"B{i}" for i in range(18)], "condition": "tumor"})
    values = np.tile(np.linspace(.1, .9, 18), (25, 1))
    a = pd.DataFrame(values, index=keys, columns=a_samples.sample_id)
    b = pd.DataFrame(values, index=keys, columns=b_samples.sample_id)
    labels = pd.DataFrame({"Katana": ["C01"]*20 + ["grey"]*5,
                           "SUPPA2_kernel": ["S01"]*20 + ["grey"]*5}, index=keys)
    return a, b, labels, a_samples, b_samples, dict.fromkeys(keys, True)


def test_same_complete_members_and_background_for_both_methods():
    attempts, families, losses = prepare_full_frozen_families(*inputs())
    assert attempts.full_module_input_eligible.all()
    assert attempts.background_events.eq(25).all()
    assert len(families) == 2
    for family in families.values():
        assert len(family.members) == 20
        assert family.labels.wgcna_module.ne("grey").sum() == 20
        assert np.isfinite(family.reference.to_numpy()).all()
        assert np.isfinite(family.target.to_numpy()).all()
    assert losses.retained_background.all()
    assert attempts.preservation_pvalue.isna().all()


def test_missing_member_never_becomes_core_replication():
    args = list(inputs())
    args[1] = args[1].iloc[1:]
    attempts, families, _ = prepare_full_frozen_families(*args)
    assert not attempts.full_module_input_eligible.any()
    assert not families
    assert attempts.original_members.eq(20).all()
    assert attempts.target_complete_patients.eq(0).all()


def test_incomplete_patients_can_be_lost_but_not_original_nodes():
    args = list(inputs())
    args[1].iloc[0, 0] = np.nan
    attempts, families, losses = prepare_full_frozen_families(*args)
    assert attempts.target_complete_patients.eq(17).all()
    assert attempts.original_members.eq(20).all()
    assert attempts.imputed_cells.eq(0).all()
    assert all(f.target.shape == (25, 17) for f in families.values())


def test_unusable_background_is_reported_not_silently_filled():
    args = list(inputs())
    args[1].iloc[-1, 0] = np.nan
    attempts, families, losses = prepare_full_frozen_families(*args)
    assert attempts.background_events.eq(24).all()
    assert len(losses.loc[~losses.retained_background]) == 2
    assert all(f.target.shape == (24, 18) for f in families.values())


def test_no_random_background_is_not_estimable():
    args = list(inputs())
    args[2] = args[2].iloc[:20]
    args[0], args[1] = args[0].iloc[:20], args[1].iloc[:20]
    attempts, families, _ = prepare_full_frozen_families(*args)
    assert not families
    assert attempts.reasons.eq("no_external_permutation_background").all()


def test_technical_duplicates_rejected_before_testing():
    args = list(inputs())
    args[4].loc[1, "patient_id"] = args[4].loc[0, "patient_id"]
    with pytest.raises(ValueError, match="Technical replicates"):
        prepare_full_frozen_families(*args)


def test_shared_biological_ids_cannot_be_independent_validation():
    args = list(inputs())
    args[4].loc[0, "patient_id"] = args[3].loc[0, "patient_id"]
    with pytest.raises(ValueError, match="Biological IDs overlap"):
        prepare_full_frozen_families(*args)


def test_full_score_preserves_weights_and_rejects_partial_sample():
    a, *_ = inputs()
    loads = pd.DataFrame({"event_key": a.index[:20], "center": .4, "scale": .2, "loading": .1})
    expected = .1 * ((a.iloc[:20] - .4) / .2).sum(axis=0)
    score = strict_frozen_scores(a, loads)
    np.testing.assert_allclose(score, expected)
    a.iloc[0, 0] = np.nan
    score = strict_frozen_scores(a, loads)
    assert np.isnan(score.iloc[0])
    np.testing.assert_allclose(score.iloc[1:], expected.iloc[1:])


def test_missing_score_node_is_not_imputed_or_renormalized():
    a, *_ = inputs()
    loads = pd.DataFrame({"event_key": a.index[:20], "center": .4, "scale": .2, "loading": .1})
    assert strict_frozen_scores(a.iloc[1:], loads).isna().all()


def test_false_identity_rejects_both_peers_without_false_pvalues():
    args = list(inputs())
    args[-1]["e0"] = False
    attempts, families, _ = prepare_full_frozen_families(*args)
    assert not families
    assert attempts.target_exact_members.eq(19).all()
    assert attempts.preservation_pvalue.isna().all()
