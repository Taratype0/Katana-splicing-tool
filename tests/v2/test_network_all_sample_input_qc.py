"""Approved cohort-wide input eligibility, independent of phenotype thresholds."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from katana_cosplice.config import StableASConfig
from katana_cosplice.qc import (
    NETWORK_INPUT_POLICY_VERSION,
    NETWORK_MIN_READS,
    NETWORK_NONCONSTANT_RULE,
    network_input_policy,
    stable_as_qc,
)


def _inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    columns = ["A1", "A2", "B1", "B2"]
    psi = pd.DataFrame(
        [[0.40, 0.41, 0.40, 0.41], [0.1, 0.2, 0.6, 0.8]],
        index=pd.Index(["low_variance_no_condition_effect", "dynamic"], name="event_key"),
        columns=columns,
    )
    ijc = pd.DataFrame(5.0, index=psi.index, columns=columns)
    sjc = pd.DataFrame(5.0, index=psi.index, columns=columns)
    samples = pd.DataFrame({"sample_id": columns, "condition": ["A", "A", "B", "B"]})
    return psi, ijc, sjc, samples


def _qc(inputs, settings=None):
    return stable_as_qc(*inputs, settings or StableASConfig())


def test_every_sample_requires_ten_without_pooling_or_group_average():
    psi, ijc, sjc, samples = _inputs()
    result = _qc((psi, ijc, sjc, samples))
    assert result.metrics["network_input_eligible"].all()
    ijc.loc["dynamic", "B2"] = 4.0
    ijc.loc["dynamic", "B1"] = 1000.0
    result = _qc((psi, ijc, sjc, samples))
    assert not result.metrics.loc["dynamic", "network_input_eligible"]
    assert result.metrics.loc["dynamic", "failed_rules"] == "sample_below_10_reads"
    assert result.metrics.loc["dynamic", "n_read_supported_samples"] == 3
    assert "dynamic" not in result.stable_psi.index


@pytest.mark.parametrize("bad_value", [np.nan, pd.NA, np.inf, -np.inf, -0.01, 1.01, "invalid", 0.5 + 0.1j])
def test_any_invalid_psi_excludes_the_event_with_an_explicit_reason(bad_value):
    psi, ijc, sjc, samples = _inputs()
    psi = psi.astype(object)
    psi.loc["dynamic", "B2"] = bad_value
    result = _qc((psi, ijc, sjc, samples))
    assert not result.metrics.loc["dynamic", "network_input_eligible"]
    assert "invalid_or_missing_psi" in result.metrics.loc["dynamic", "failed_rules"]
    assert result.metrics.loc["low_variance_no_condition_effect", "network_input_eligible"]


@pytest.mark.parametrize("matrix_name", ["ijc", "sjc"])
@pytest.mark.parametrize("bad_value", [np.nan, pd.NA, np.inf, -np.inf, -1.0, "invalid", 10 + 1j])
def test_both_count_components_must_be_valid_even_if_the_other_is_ten(matrix_name, bad_value):
    psi, ijc, sjc, samples = _inputs()
    ijc = ijc.astype(object)
    sjc = sjc.astype(object)
    target, other = (ijc, sjc) if matrix_name == "ijc" else (sjc, ijc)
    target.loc["dynamic", "A1"] = bad_value
    other.loc["dynamic", "A1"] = 10.0
    result = _qc((psi, ijc, sjc, samples))
    assert not result.metrics.loc["dynamic", "network_input_eligible"]
    assert "invalid_or_missing_counts" in result.metrics.loc["dynamic", "failed_rules"]


def test_overflowing_total_is_invalid_not_infinite_read_support():
    psi, ijc, sjc, samples = _inputs()
    ijc.loc["dynamic", "A1"] = np.finfo(float).max
    sjc.loc["dynamic", "A1"] = np.finfo(float).max
    result = _qc((psi, ijc, sjc, samples))
    assert not result.metrics.loc["dynamic", "network_input_eligible"]
    assert "invalid_or_missing_counts" in result.metrics.loc["dynamic", "failed_rules"]


def test_low_variance_repeatability_and_condition_range_are_not_hard_gates():
    result = _qc(_inputs())
    row = result.metrics.loc["low_variance_no_condition_effect"]
    assert row.psi_variance < 0.0025
    assert row.replicate_repeatability < 0.30
    assert row.condition_mean_psi_range < 0.25
    assert row.sample_psi_range < 0.50
    assert row.network_input_eligible
    assert row.stable_as
    assert row.failed_rules == ""
    assert "replicate_repeatability" in row.diagnostic_only_metrics
    assert not {"passes_repeatability", "passes_condition_support", "passes_dynamic_range"} & set(result.summary.rule)


def test_only_exactly_constant_events_are_removed_without_a_variance_floor():
    psi, ijc, sjc, samples = _inputs()
    psi.loc["dynamic"] = 0.5
    psi.loc["low_variance_no_condition_effect"] = [0.5, 0.5, 0.5, np.nextafter(0.5, 1.0)]
    result = _qc((psi, ijc, sjc, samples))
    assert result.metrics.loc["dynamic", "failed_rules"] == "constant_psi"
    assert not result.metrics.loc["dynamic", "network_input_eligible"]
    assert result.metrics.loc["low_variance_no_condition_effect", "network_input_eligible"]
    assert result.metrics.loc["low_variance_no_condition_effect", "psi_variance"] < 1e-12


def test_one_sample_per_condition_does_not_restore_a_minimum_two_filter():
    psi, ijc, sjc, samples = _inputs()
    columns = ["A1", "B1"]
    result = _qc((psi[columns], ijc[columns], sjc[columns], samples[samples.sample_id.isin(columns)]))
    assert result.metrics.loc["dynamic", "network_input_eligible"]
    assert result.condition_metrics.n_replicates.eq(1).all()
    assert result.condition_metrics.passes_all_samples.all()


def test_condition_relabeling_reassignment_and_metadata_order_do_not_change_eligibility_or_cap_ranking():
    psi, ijc, sjc, samples = _inputs()
    expected = _qc((psi, ijc, sjc, samples))
    renamed = samples.iloc[::-1].copy()
    renamed["condition"] = ["x", "y", "x", "y"]
    actual = _qc((psi, ijc, sjc, renamed))
    pd.testing.assert_frame_equal(expected.stable_psi, actual.stable_psi)
    for field in ("network_input_eligible", "stable_as", "failed_rules", "stability_score"):
        pd.testing.assert_series_equal(expected.metrics[field], actual.metrics[field])


def test_complete_cohort_cannot_be_silently_shrunk_by_metadata_or_count_axes():
    psi, ijc, sjc, samples = _inputs()
    with pytest.raises(ValueError, match="metadata and matrix columns differ"):
        _qc((psi, ijc, sjc, samples.iloc[:-1]))
    with pytest.raises(ValueError, match="same ordered sample columns"):
        _qc((psi, ijc.iloc[:, :-1], sjc, samples))
    with pytest.raises(ValueError, match="same ordered sample columns"):
        _qc((psi, ijc.iloc[:, ::-1], sjc, samples))
    with pytest.raises(ValueError, match="same ordered event_key rows"):
        _qc((psi, ijc.iloc[::-1], sjc, samples))


@pytest.mark.parametrize("source", ["psi_columns", "psi_index", "metadata", "metadata_normalized"])
def test_duplicate_sample_or_event_identities_are_rejected(source):
    psi, ijc, sjc, samples = _inputs()
    if source == "psi_columns":
        psi.columns = ["A1", "A1", "B1", "B2"]
    elif source == "psi_index":
        psi.index = ["duplicate", "duplicate"]
    elif source == "metadata":
        samples.loc[1, "sample_id"] = "A1"
    else:
        samples.loc[1, "sample_id"] = " A1 "
    with pytest.raises(ValueError, match="must be unique"):
        _qc((psi, ijc, sjc, samples))


def test_legacy_presets_and_numeric_overrides_are_recorded_but_never_applied():
    inputs = _inputs()
    ordinary = _qc(inputs)
    settings = StableASConfig(
        preset="strict",
        min_total_junction_reads=1000.0,
        min_read_supported_fraction=0.0,
        min_supported_replicates_per_condition=20.0,
        min_psi_variance=1.0,
        min_repeatability=1.0,
        min_condition_mean_psi_range=1.0,
    )
    strict = _qc(inputs, settings)
    pd.testing.assert_frame_equal(ordinary.stable_psi, strict.stable_psi)
    assert strict.metrics.network_input_eligible.all()
    assert strict.metrics.min_total_junction_reads.eq(10.0).all()
    ignored = json.loads(strict.metrics.iloc[0].ignored_legacy_qc_settings)
    assert ignored["preset"] == "strict"
    assert ignored["min_total_junction_reads"] == 1000.0


def test_native_input_matrices_and_metadata_remain_unchanged_and_no_imputation_occurs():
    inputs = _inputs()
    originals = [frame.copy(deep=True) for frame in inputs]
    result = _qc(inputs)
    for frame, original in zip(inputs, originals):
        pd.testing.assert_frame_equal(frame, original)
    pd.testing.assert_frame_equal(result.stable_psi, originals[0])
    assert result.stable_psi.columns.tolist() == originals[0].columns.tolist()


def test_policy_describes_fixed_counts_complete_cohort_and_phenotype_blind_cap_ranking():
    policy = network_input_policy()
    assert policy["policy_id"] == NETWORK_INPUT_POLICY_VERSION
    assert policy["min_count"] == NETWORK_MIN_READS == 10
    assert policy["constant_criterion"] == NETWORK_NONCONSTANT_RULE == "sample_max_gt_min"
    assert policy["differential_preselection"] is False
    assert policy["optional_cap_ranking"] == "cohort_PSI_variance; condition independent"
    policy["min_count"] = 0
    assert network_input_policy()["min_count"] == 10


def test_empty_event_matrix_is_handled_but_an_empty_cohort_is_not_valid():
    psi, ijc, sjc, samples = _inputs()
    result = _qc((psi.iloc[:0], ijc.iloc[:0], sjc.iloc[:0], samples))
    assert result.metrics.empty
    assert result.stable_psi.empty
    assert result.summary.n_passed.eq(0).all()
    with pytest.raises(ValueError, match="cohort must contain declared samples"):
        _qc((psi.iloc[:, :0], ijc.iloc[:, :0], sjc.iloc[:, :0], samples.iloc[:0]))
