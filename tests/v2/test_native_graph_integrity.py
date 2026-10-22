import copy
import csv

import pytest

from katana_next.science.native_graph_integrity import (
    TABLE_COLUMNS, build_native_definition_snapshot, compare_native_definition_snapshots)


def fixture(tmp_path):
    records = {
        "graph": [["1", "H_sapiens", "TEST", "1", "+", "100", "500"]],
        "as_ref": [["11", "1", "ES", "1", "4", "0", "2:3", "hash", "0.2"]],
        "as_ref_exon": [["11", "2"], ["11", "3"]],
        "exon": [["1", str(i), str(100 + i * 10), str(109 + i * 10), str(i),
                  str(i * 10), str(i * 10 + 9), "1", "ACGT"] for i in range(1, 5)],
        "version": [["20140108"]],
    }
    rows = [dict(event_key="fixture|ES|11", native_as_id="11", original_as_id="11",
                 composite_event="True", exact_full_path_identity="True", reasons="")]
    write(tmp_path, records)
    return records, rows


def write(path, records):
    for table, values in records.items():
        with (path / (table + ".tsv")).open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(TABLE_COLUMNS[table])
            writer.writerows(values)


def test_full_composite_and_reference_version_are_retained(tmp_path):
    records, rows = fixture(tmp_path)
    a = build_native_definition_snapshot(tmp_path, rows)
    b = build_native_definition_snapshot(tmp_path, rows)
    result = compare_native_definition_snapshots(a, b, rows)
    assert result["all_original_definitions_unchanged"]
    assert result["original_event_count"] == 1
    assert not result["external_preservation_tested"]


@pytest.mark.parametrize("table,row,column,value", [
    ("as_ref", 0, 6, "2"),
    ("as_ref", 0, 5, "1"),
    ("as_ref", 0, 7, "different_hash"),
    ("as_ref_exon", 1, 1, "9"),
    ("exon", 2, 2, "999"),
    ("exon", 2, 8, "TTTT"),
    ("graph", 0, 4, "-"),
    ("graph", 0, 2, "OTHER"),
    ("version", 0, 0, "20160101"),
])
def test_original_path_coordinate_sequence_flags_and_version_changes_detected(tmp_path, table, row, column, value):
    records, rows = fixture(tmp_path)
    baseline = build_native_definition_snapshot(tmp_path, rows)
    records[table][row][column] = value
    write(tmp_path, records)
    current = build_native_definition_snapshot(tmp_path, rows)
    result = compare_native_definition_snapshots(baseline, current, rows)
    assert result["changed_original_event_count"] == 1
    assert not result["all_original_definitions_unchanged"]


def test_sample_dependent_expression_does_not_redefine_event(tmp_path):
    records, rows = fixture(tmp_path)
    baseline = build_native_definition_snapshot(tmp_path, rows)
    records["as_ref"][0][-1] = "0.9"
    write(tmp_path, records)
    assert compare_native_definition_snapshots(baseline,
        build_native_definition_snapshot(tmp_path, rows), rows)["all_original_definitions_unchanged"]


def test_existing_mismatch_is_not_silently_repaired(tmp_path):
    records, rows = fixture(tmp_path)
    rows[0].update(exact_full_path_identity="False", reasons="original_mismatch")
    snapshot = build_native_definition_snapshot(tmp_path, rows)
    result = compare_native_definition_snapshots(snapshot, snapshot, rows)
    assert result["inherited_identity_mismatches"] == 1
    assert not result["missing_values_or_mismatches_repaired"]


def test_duplicate_original_inclusive_link_rejected(tmp_path):
    records, rows = fixture(tmp_path)
    records["as_ref_exon"].append(records["as_ref_exon"][0])
    write(tmp_path, records)
    with pytest.raises(ValueError, match="Duplicate original inclusive"):
        build_native_definition_snapshot(tmp_path, rows)


def test_deleted_original_event_detected_without_deleting_membership(tmp_path):
    records, rows = fixture(tmp_path)
    baseline = build_native_definition_snapshot(tmp_path, rows)
    records["as_ref"] = []
    write(tmp_path, records)
    current = build_native_definition_snapshot(tmp_path, rows)
    assert current["missing_native_as_ids"] == ["11"]
    assert len(current["event_definition_sha256"]) == 1
    assert not compare_native_definition_snapshots(baseline, current, rows)["all_original_definitions_unchanged"]


def test_partial_event_universe_refused(tmp_path):
    records, rows = fixture(tmp_path)
    baseline = build_native_definition_snapshot(tmp_path, rows)
    partial = copy.deepcopy(baseline)
    partial["event_definition_sha256"] = {}
    with pytest.raises(ValueError, match="every original event"):
        compare_native_definition_snapshots(baseline, partial, rows)
