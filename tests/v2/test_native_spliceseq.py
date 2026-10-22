import numpy as np
import pandas as pd
import pytest

from katana_next.science.native_spliceseq import read_native_spliceseq_counts


def data(tmp_path, rows="1\tT1\t101\t1\t0.23\t8\t2\n"):
    path = tmp_path / "native.tsv"
    path.write_text("sample_id\tname\tas_id\tsample_id\tpsi\tincluded_counts\texcluded_counts\n" + rows)
    events = pd.DataFrame({"event_key": ["a", "b"], "native_as_id": [101, 102],
                           "exact_full_path_identity": [True, False]})
    samples = pd.DataFrame({"sample_id": ["T1", "T2"], "patient_id": ["p1", "p2"],
                            "condition": ["tumor", "tumor"]})
    return path, events, samples


def test_native_psi_is_not_counts_fraction_and_all_nodes_remain(tmp_path):
    path, events, samples = data(tmp_path)
    result = read_native_spliceseq_counts([path], events, samples)
    assert result.supported_psi.loc["a", "T1"] == .23
    assert result.supported_psi.shape == (2, 2)
    assert np.isnan(result.supported_psi.loc["b", "T1"])
    assert not result.sample_audit.set_index("sample_id").loc["T2", "library_observed_in_export"]


@pytest.mark.parametrize("row", ["1\tT1\t101\t1\t0.23\t7\t2\n", "1\tT1\t101\t1\t0.23\t\\N\t2\n"])
def test_low_or_unknown_support_is_missing_not_zero(tmp_path, row):
    path, events, samples = data(tmp_path, row)
    result = read_native_spliceseq_counts([path], events, samples)
    assert result.raw_psi.loc["a", "T1"] == .23
    assert np.isnan(result.supported_psi.loc["a", "T1"])


@pytest.mark.parametrize("row,match", [
    ("1\tT1\t101\t2\t0.23\t8\t2\n", "disagree"),
    ("1\tunknown\t101\t1\t0.23\t8\t2\n", "Undeclared"),
    ("1\tT1\t101\t1\t1.2\t8\t2\n", "outside"),
    ("1\tT1\t101\t1\t0.2\t-1\t2\n", "Invalid"),
    ("1\tT1\t101\t1\t0.2\t8.5\t2\n", "Invalid"),
])
def test_invalid_or_unmapped_measurements_fail(tmp_path, row, match):
    path, events, samples = data(tmp_path, row)
    with pytest.raises(ValueError, match=match):
        read_native_spliceseq_counts([path], events, samples)


def test_duplicate_measurement_cannot_be_averaged(tmp_path):
    path, events, samples = data(tmp_path)
    with pytest.raises(ValueError, match="Duplicate frozen"):
        read_native_spliceseq_counts([path, path], events, samples)


def test_paired_tissues_are_not_duplicate_patients(tmp_path):
    path, events, samples = data(tmp_path)
    samples.loc[1, ["patient_id", "condition"]] = ["p1", "normal"]
    assert len(read_native_spliceseq_counts([path], events, samples).sample_audit) == 2
    samples.loc[1, "condition"] = "tumor"
    with pytest.raises(ValueError, match="One declared"):
        read_native_spliceseq_counts([path], events, samples)
