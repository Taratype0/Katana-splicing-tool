import numpy as np
import pandas as pd
import pytest

from katana_cosplice.bulk_stability import (
    align_and_classify, boundary_stages, core_kme_rescue, module_topology, retain_ranked,
)


def test_label_renaming_is_not_switching():
    _, status = align_and_classify(["a", "a", "b", "grey"], ["z", "z", "y", "grey"])
    assert status.tolist() == ["retained_same"] * 3 + ["shared_grey"]


def test_grey_is_not_a_matching_module():
    _, status = align_and_classify(["a", "a", "grey"], ["grey", "grey", "z"])
    assert status.tolist() == ["dropped", "dropped", "newly_assigned"]


def test_pair_losses_decompose_exactly():
    result = module_topology(["a"] * 4, ["x", "x", "y", "grey"]).iloc[0]
    assert result.pair_preserved_fraction == pytest.approx(1 / 6)
    assert result.pair_split_fraction == pytest.approx(2 / 6)
    assert result.pair_loss_to_grey_fraction == pytest.approx(3 / 6)


def test_whole_module_loss():
    result = module_topology(["a", "a"], ["grey", "grey"]).iloc[0]
    assert result.module_completely_lost
    assert result.pair_loss_to_grey_fraction == 1


def test_boundary_stages_include_module_size_loss():
    m = pd.DataFrame({"observed_module": ["a"] * 3 + ["b"] * 2,
                      "consensus_module": ["C01", "C01", "grey", "grey", "grey"],
                      "consensus_membership_confidence": [.8, .6, .2, .8, .1],
                      "next_best_module_coassignment": [.1] * 5,
                      "consensus_confidence_margin": [.7, .5, .1, .7, 0]})
    settings = {"core_membership_min": .7, "peripheral_membership_min": .45,
                "peripheral_membership_margin": .05, "min_module_size": 2}
    pre, final, reason = boundary_stages(m, settings)
    assert pre[3] == "b" and final[3] == "grey"
    assert reason[3] == "small_module_after_filter"
    m.loc[4, "consensus_module"] = "C01"
    with pytest.raises(ValueError, match="reconstructed"):
        boundary_stages(m, settings)


def test_equal_coverage_ranked_control():
    labels = ["a", "a", "b", "b", "grey"]
    result = retain_ranked(labels, [1, 2, 4, 3, 100], target=["x", "grey", "grey", "y", "grey"])
    assert result.tolist() == ["grey", "a", "b", "grey", "grey"]


def test_invalid_or_singleton_inputs():
    with pytest.raises(ValueError):
        align_and_classify(["a"], ["a", "b"])
    with pytest.raises(ValueError):
        retain_ranked(["a"], [1], fraction=2)
    assert np.isnan(module_topology(["a"], ["a"]).iloc[0].pair_split_fraction)


def test_merges_are_detected_without_label_mapping():
    result = module_topology(["a", "a", "b", "b"], ["x"] * 4)
    assert result.merging_cross_reference_pairs.tolist() == [4, 4]
    assert result.pair_preserved_fraction.tolist() == [1., 1.]


def test_public_biology_output_contract(tmp_path):
    from scripts.audit_bulk_boundaries import write_known_signal
    matrix = pd.DataFrame([[.1, .2], [.2, .8]],
                          index=["g1;SE:chr1:10-20:30-50:+", "g2;SE:chr1:10-100:200-300:+"],
                          columns=["a", "b"])
    write_known_signal(matrix, {"arm": np.array(["M1", "M1"])}, pd.Series(["a", "b"]), tmp_path)
    assert pd.read_csv(tmp_path / "known_microexon_signal.tsv", sep="\t").iloc[0].microexons == 1
    assert (tmp_path / "microexon_members.tsv").is_file()


def test_core_rescue_never_promotes_grey_or_ambiguous_events():
    m = pd.DataFrame({"observed_module": ["a", "a", "a", "grey"],
                      "consensus_module": ["grey"] * 4,
                      "consensus_membership_confidence": [.2] * 4,
                      "next_best_module_coassignment": [.1, .1, .3, .1],
                      "consensus_confidence_margin": [.1, .1, -.1, .1],
                      "bootstrap_mean_abs_kME": [.8] * 4,
                      "bootstrap_q10_abs_kME": [.6] * 4,
                      "bootstrap_kME_sign_consistency": [.9] * 4})
    settings = {"core_membership_min": .7, "peripheral_membership_min": .45,
                "peripheral_membership_margin": .05, "min_module_size": 2,
                "core_kme_mean_min": .7, "core_kme_lower_quantile_min": .5,
                "core_kme_sign_consistency_min": .8}
    assert core_kme_rescue(m, settings).tolist() == ["a", "a", "grey", "grey"]
    m.loc[1, "bootstrap_q10_abs_kME"] = np.nan
    assert core_kme_rescue(m, settings).tolist() == ["grey"] * 4
