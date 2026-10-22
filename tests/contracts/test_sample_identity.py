"""Sample-sheet identity checks at the nf-core / Katana metadata boundary.

All reads are synthetic.  These checks exercise record identity, lane grouping,
and explicit experimental metadata; they never merge reads, launch a workflow,
infer pairings, or generate comparisons.
"""
from __future__ import annotations

from copy import deepcopy
import csv
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from katana_next.contracts.sample_identity import (
    check_project_sample_mapping,
    validate_sample_rows,
    validate_samplesheets,
)


def reads(sample="S1", lane="L1", *, paired=True, strand="reverse"):
    return {
        "sample": sample,
        "fastq_1": f"reads/{sample}_{lane}_R1.fastq.gz",
        "fastq_2": f"reads/{sample}_{lane}_R2.fastq.gz" if paired else "",
        "strandedness": strand,
    }


def metadata(sample="S1", condition="treated", unit=None, donor="", **extra):
    return {
        "sample": sample,
        "condition": condition,
        "biological_unit_id": unit if unit is not None else f"bio_{sample}",
        "donor_id": donor,
        **extra,
    }


def groups(report):
    return {item["sample_id"]: item for item in report["upstream_groups"]}


def samples(report):
    return {item["sample_id"]: item for item in report["samples"]}


def accepted(upstream=None, meta=None, **kwargs):
    report = validate_sample_rows(
        upstream if upstream is not None else [reads()],
        meta if meta is not None else [metadata()],
        **kwargs,
    )
    assert report["valid"] is True, report
    assert not [item for item in report["issues"] if item["severity"] == "error"], report
    assert report["execution_authorized"] is False
    assert report["scientific_tools_executed"] == 0
    assert report["merges_executed"] == 0
    assert report["contrasts_generated"] == 0
    return report


def rejected(upstream=None, meta=None, **kwargs):
    report = validate_sample_rows(
        upstream if upstream is not None else [reads()],
        meta if meta is not None else [metadata()],
        **kwargs,
    )
    assert report["valid"] is False, report
    assert any(item["severity"] == "error" for item in report["issues"]), report
    assert report["execution_authorized"] is False
    assert report["scientific_tools_executed"] == 0
    assert report["merges_executed"] == 0
    assert report["contrasts_generated"] == 0
    for item in report["issues"]:
        assert {"code", "path", "message", "severity"} <= item.keys()
    return report


def has_warning(report):
    assert any(item["severity"] == "warning" for item in report["issues"]), report


def write_sheet(path, rows, columns=None, *, encoding="utf-8"):
    columns = columns or list(rows[0])
    with path.open("w", encoding=encoding, newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    return path


def project_from(report):
    return {
        "kind": "project",
        "schema_version": "1.0.0",
        "project_id": "sample_identity_fixture",
        "display_name": "Synthetic sample identity mapping",
        "assay_mode": report["assay_mode"],
        "analysis_level": "bulk",
        "input_kind": "processed_results",
        "intent": "import",
        "samples": [
            {
                **{
                    key: row[key]
                    for key in (
                        "sample_id", "group", "biological_unit_id", "donor_id", "cell_type"
                    )
                },
                "display_name": row["display_name"],
                "source_assay": report["assay_mode"],
                "covariates": {},
            }
            for row in report["samples"]
        ],
        "contrasts": [],
        "references": [],
        "assets": [],
        "methods": [],
        "runs": [],
        "results": [],
    }


def test_one_canonical_sample_preserves_explicit_condition_and_unknown_donor():
    report = accepted([reads("control_named_file")], [metadata("control_named_file", "drug")])
    row = samples(report)["control_named_file"]
    assert row["group"] == row["condition"] == "drug"
    assert row["donor_id"] is None
    assert row["biological_unit_id"] == "bio_control_named_file"
    assert report["counts"]["sample_count"] == 1
    assert report["counts"]["donor_count"] == 0


def test_multiple_lanes_group_by_stable_id_without_sorting_pairs_or_merging():
    upstream = [reads("S2", "L2"), reads("S1", "L9"), reads("S2", "L1")]
    meta = [metadata("S1"), metadata("S2")]
    original = deepcopy((upstream, meta))
    report = accepted(upstream, meta)
    grouped = groups(report)
    assert grouped["S2"]["rows"] == [upstream[0], upstream[2]]
    assert grouped["S1"]["rows"] == [upstream[1]]
    assert grouped["S2"]["row_count"] == 2
    assert grouped["S2"]["layout"] == "paired"
    assert grouped["S2"]["merge_owner"] == "nf-core"
    assert grouped["S2"]["merge_performed"] is False
    assert report["counts"]["sample_count"] == 2
    assert report["counts"]["upstream_row_count"] == 3
    assert (upstream, meta) == original
    grouped["S2"]["rows"][0]["fastq_1"] = "result_is_a_copy"
    assert (upstream, meta) == original


def test_single_end_lane_group_is_explicit_and_supported():
    report = accepted([reads(lane="L3", paired=False), reads(lane="L1", paired=False)])
    assert groups(report)["S1"]["layout"] == "single"
    assert groups(report)["S1"]["row_count"] == 2


@pytest.mark.parametrize("strand", ["forward", "reverse", "unstranded", "auto"])
def test_supported_nfcore_strandedness_is_preserved(strand):
    report = accepted([reads(strand=strand)])
    assert groups(report)["S1"]["rows"][0]["strandedness"] == strand


def test_identical_upstream_row_is_duplicate_evidence_not_an_additional_lane():
    row = reads()
    rejected([row, deepcopy(row)])


@pytest.mark.parametrize("same_sample", [True, False])
def test_same_read_path_cannot_be_counted_twice(same_sample):
    row1 = reads("S1", "L1")
    sample2 = "S1" if same_sample else "S2"
    row2 = reads(sample2, "L2")
    row2["fastq_1"] = row1["fastq_1"]
    meta = [metadata("S1")] + ([] if same_sample else [metadata("S2")])
    rejected([row1, row2], meta)


def test_read_one_and_read_two_cannot_reference_the_same_file():
    row = reads()
    row["fastq_2"] = row["fastq_1"]
    rejected([row])


def test_reusing_a_mate_as_another_samples_read_one_is_rejected():
    row1, row2 = reads("S1"), reads("S2")
    row2["fastq_1"] = row1["fastq_2"]
    rejected([row1, row2], [metadata("S1"), metadata("S2")])


@pytest.mark.parametrize(
    "first,second",
    [
        ("/reads/one/../lane_R1.fastq.gz", "/reads/lane_R1.fastq.gz"),
        (r"C:\Reads\one\..\Lane_R1.fastq.gz", "c:/reads/lane_R1.fastq.gz"),
    ],
)
def test_duplicate_local_paths_use_platform_appropriate_lexical_normalization(first, second):
    row1, row2 = reads("S1"), reads("S2")
    row1["fastq_1"], row2["fastq_1"] = first, second
    rejected([row1, row2], [metadata("S1"), metadata("S2")])


def test_malformed_remote_uri_returns_an_invalid_report_instead_of_crashing():
    row = reads()
    row["fastq_1"] = "s3://[bad/file.fastq.gz"
    rejected([row])


@pytest.mark.parametrize("missing", ["", "   "])
def test_fastq_entry_requires_read_one(missing):
    row = reads()
    row["fastq_1"] = missing
    rejected([row])


def test_same_sample_cannot_mix_single_and_paired_end_lanes():
    rejected([reads(lane="L1"), reads(lane="L2", paired=False)])


def test_same_sample_cannot_mix_known_strandedness_across_lanes():
    rejected([reads(lane="L1", strand="reverse"), reads(lane="L2", strand="forward")])


def test_invalid_strandedness_is_not_silently_guessed():
    rejected([reads(strand="guess-from-name")])


def test_genome_bam_adapter_uses_explicit_ids_and_does_not_need_a_local_bai():
    report = accepted(
        [{"sample": "S1", "bam": "s3://synthetic/S1.bam", "bai": ""}],
        input_format="bam",
    )
    assert groups(report)["S1"]["layout"] == "bam"
    assert groups(report)["S1"]["rows"][0]["bam"] == "s3://synthetic/S1.bam"


def test_official_nfcore_genome_bam_column_is_supported_without_renaming():
    row = {"sample": "S1", "genome_bam": "/synthetic/S1.bam"}
    report = accepted([row], input_format="bam")
    assert groups(report)["S1"]["rows"] == [row]


def test_conflicting_genome_bam_and_legacy_bam_alias_are_ambiguous():
    rejected(
        [{"sample": "S1", "genome_bam": "/synthetic/a.bam", "bam": "/synthetic/b.bam"}],
        input_format="bam",
    )


def test_transcriptome_bam_does_not_replace_the_required_genome_bam():
    rejected([{"sample": "S1", "transcriptome_bam": "/synthetic/S1.transcriptome.bam"}], input_format="bam")


def test_multiple_bam_rows_do_not_silently_reuse_fastq_lane_merging_semantics():
    rejected(
        [{"sample": "S1", "genome_bam": "/synthetic/a.bam"}, {"sample": "S1", "genome_bam": "/synthetic/b.bam"}],
        input_format="bam",
    )


@pytest.mark.parametrize("assay", ["single_cell_rna", "single_nucleus_rna"])
def test_existing_single_cell_bam_sheet_check_confers_no_processing_eligibility(assay):
    report = accepted(
        [{"sample": "S1", "genome_bam": "/synthetic/aligned.bam"}],
        [metadata(donor="D1", cell_type="T_cell")],
        input_format="bam", assay_mode=assay,
    )
    assert report["method_eligibility_checked"] is False
    assert report["execution_authorized"] is False
    assert report["scientific_tools_executed"] == 0


def test_sample_ids_adapter_works_for_processed_inputs_without_fabricated_raw_files():
    report = accepted([{"sample": "S1"}], input_format="sample_ids")
    assert groups(report)["S1"]["layout"] == "sample_ids"
    assert groups(report)["S1"]["rows"] == [{"sample": "S1"}]


@pytest.mark.parametrize("assay", ["single_cell_rna", "single_nucleus_rna"])
def test_bulk_upstream_does_not_expand_to_raw_single_cell_fastq(assay):
    rejected(assay_mode=assay)


@pytest.mark.parametrize("assay", ["single_cell_rna", "single_nucleus_rna"])
def test_processed_single_cell_identities_keep_donor_and_cell_type(assay):
    report = accepted(
        [{"sample": "S1"}],
        [metadata(donor="D1", cell_type="T_cell")],
        input_format="sample_ids", assay_mode=assay,
    )
    assert samples(report)["S1"]["donor_id"] == "D1"
    assert samples(report)["S1"]["cell_type"] == "T_cell"
    assert samples(report)["S1"]["source_assay"] == assay


@pytest.mark.parametrize("duplicate", [metadata(), metadata(condition="other")])
def test_even_identical_duplicate_metadata_rows_need_resolution(duplicate):
    rejected(meta=[metadata(), duplicate])


def test_missing_metadata_sample_is_an_error():
    rejected([reads("S1"), reads("S2")], [metadata("S1")])


def test_metadata_cannot_silently_add_a_sample_absent_from_the_inputs():
    rejected(meta=[metadata("S1"), metadata("S2")])


def test_empty_condition_is_not_inferred_from_sample_or_file_name():
    rejected([reads("control_rep1")], [metadata("control_rep1", condition="")])


def test_empty_biological_unit_is_retained_as_unknown_with_a_warning():
    report = accepted(meta=[metadata(unit="")])
    assert samples(report)["S1"]["biological_unit_id"] is None
    assert report["counts"]["unknown_biological_unit_count"] == 1
    has_warning(report)


def test_two_unknown_units_are_not_inferred_independent_from_sample_names():
    report = accepted(
        [reads("mouse1_rep1"), reads("mouse2_rep1")],
        [metadata("mouse1_rep1", unit=""), metadata("mouse2_rep1", unit="")],
    )
    assert report["counts"]["unknown_biological_unit_count"] == 2
    assert report["counts"]["biological_unit_count"] == 0
    assert all(item["biological_unit_id"] is None for item in report["samples"])
    has_warning(report)


def test_distinct_biological_units_can_share_a_donor_without_being_merged():
    report = accepted(
        [reads("S1"), reads("S2")],
        [metadata("S1", donor="D1"), metadata("S2", donor="D1")],
    )
    assert report["counts"]["sample_count"] == 2
    assert report["counts"]["biological_unit_count"] == 2
    assert report["counts"]["donor_count"] == 1
    assert report["related_donors"]
    has_warning(report)


def test_one_biological_unit_cannot_map_to_two_conflicting_known_donors():
    rejected(
        [reads("S1"), reads("S2")],
        [metadata("S1", unit="B1", donor="D1"), metadata("S2", unit="B1", donor="D2")],
    )


def test_shared_unit_across_conditions_is_visible_without_automatic_pairing():
    report = accepted(
        [reads("S1"), reads("S2")],
        [metadata("S1", "before", unit="B1", donor="D1"), metadata("S2", "after", unit="B1", donor="D1")],
    )
    assert report["counts"]["biological_unit_count"] == 1
    has_warning(report)
    assert report["contrasts_generated"] == 0
    assert not report.get("contrasts")


def test_duplicate_display_names_warn_but_do_not_merge_stable_ids():
    report = accepted(
        [reads("S1"), reads("S2")],
        [metadata("S1", display_name="A mouse"), metadata("S2", display_name="A mouse")],
    )
    assert set(samples(report)) == {"S1", "S2"}
    has_warning(report)


def test_valid_technical_replicates_remain_separate_records_with_explicit_parent():
    report = accepted(
        [reads("S1"), reads("S1_repeat")],
        [
            metadata("S1", unit="B1", donor="D1"),
            metadata("S1_repeat", unit="B1", donor="D1", technical_replicate_of="S1"),
        ],
    )
    assert report["counts"]["sample_count"] == 2
    assert report["counts"]["biological_unit_count"] == 1
    assert report["counts"]["technical_replicate_count"] == 1
    assert samples(report)["S1_repeat"]["technical_replicate_of"] == "S1"


@pytest.mark.parametrize("parent", ["not_in_project", "S1"])
def test_technical_parent_must_exist_and_cannot_be_self(parent):
    rejected(meta=[metadata(technical_replicate_of=parent)])


def test_technical_replicate_cycles_are_not_accepted():
    rejected(
        [reads("S1"), reads("S2")],
        [
            metadata("S1", unit="B1", technical_replicate_of="S2"),
            metadata("S2", unit="B1", technical_replicate_of="S1"),
        ],
    )


@pytest.mark.parametrize(
    "field,value",
    [
        ("condition", "different_condition"),
        ("biological_unit_id", "different_unit"),
        ("donor_id", "different_donor"),
        ("cell_type", "different_cell_type"),
    ],
)
def test_technical_replicate_metadata_must_be_compatible(field, value):
    parent = metadata("S1", unit="B1", donor="D1", cell_type="T_cell")
    child = metadata("S2", unit="B1", donor="D1", cell_type="T_cell", technical_replicate_of="S1")
    child[field] = value
    rejected([reads("S1"), reads("S2")], [parent, child])


def test_explicit_covariates_and_batch_are_preserved_as_metadata():
    report = accepted(meta=[metadata(batch="batch_2", covariate_sex="F", covariate_age="12")])
    row = samples(report)["S1"]
    assert row["batch"] == "batch_2"
    assert row["covariates"]["sex"] == "F"
    assert row["covariates"]["age"] == "12"


def test_unrecognized_metadata_column_warns_instead_of_discarding_a_valid_project():
    report = accepted(meta=[metadata(operator_note="keep this explicit note")])
    has_warning(report)


@pytest.mark.parametrize(
    "field,value",
    [
        ("condition", "wrong"),
        ("group", "wrong"),
        ("biological_unit_id", "wrong"),
        ("donor_id", "wrong"),
        ("cell_type", "wrong"),
    ],
)
def test_metadata_present_in_both_sheets_must_agree(field, value):
    row = reads()
    row[field] = value
    rejected([row])


def test_matching_upstream_metadata_does_not_require_removing_nfcore_columns():
    row = reads()
    row.update(condition="treated", biological_unit_id="bio_S1", donor_id="D1", cell_type="T_cell")
    accepted([row], [metadata(donor="D1", cell_type="T_cell")])


def test_file_reader_accepts_utf8_bom_and_quoted_commas_without_rewriting(tmp_path):
    upstream = tmp_path / "reads.csv"
    meta = tmp_path / "metadata.csv"
    write_sheet(upstream, [reads()], encoding="utf-8-sig")
    write_sheet(meta, [metadata(condition="drug, high dose", display_name="小鼠, one")], encoding="utf-8-sig")
    before = {path.name: path.read_bytes() for path in (upstream, meta)}
    report = validate_samplesheets(upstream, meta)
    assert report["valid"] is True, report
    assert samples(report)["S1"]["group"] == "drug, high dose"
    assert samples(report)["S1"]["display_name"] == "小鼠, one"
    assert {path.name: path.read_bytes() for path in (upstream, meta)} == before


@pytest.mark.parametrize("target", ["upstream", "metadata"])
def test_duplicate_csv_headers_are_rejected_before_dict_reader_loses_evidence(tmp_path, target):
    upstream = write_sheet(tmp_path / "reads.csv", [reads()])
    meta = write_sheet(tmp_path / "metadata.csv", [metadata()])
    path = upstream if target == "upstream" else meta
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle))
    rows[0].append(rows[0][0])
    rows[1].append(rows[1][0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        csv.writer(handle).writerows(rows)
    report = validate_samplesheets(upstream, meta)
    assert report["valid"] is False, report


@pytest.mark.parametrize("target", ["upstream", "metadata"])
@pytest.mark.parametrize("extra_cell", [True, False])
def test_csv_rows_with_extra_or_missing_cells_are_rejected(tmp_path, target, extra_cell):
    upstream = write_sheet(tmp_path / "reads.csv", [reads()])
    meta = write_sheet(tmp_path / "metadata.csv", [metadata()])
    path = upstream if target == "upstream" else meta
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle))
    if extra_cell:
        rows[1].append("unheaded-value")
    else:
        rows[1].pop()
    with path.open("w", encoding="utf-8", newline="") as handle:
        csv.writer(handle).writerows(rows)
    report = validate_samplesheets(upstream, meta)
    assert report["valid"] is False, report


@pytest.mark.parametrize("missing_header", ["condition", "biological_unit_id", "donor_id"])
def test_metadata_headers_are_required_even_when_values_can_be_unknown(tmp_path, missing_header):
    upstream = write_sheet(tmp_path / "reads.csv", [reads()])
    row = metadata()
    del row[missing_header]
    meta = write_sheet(tmp_path / "metadata.csv", [row])
    report = validate_samplesheets(upstream, meta)
    assert report["valid"] is False, report


def test_project_mapping_accepts_same_identities_without_mutating_either_document():
    sample_report = accepted()
    project = project_from(sample_report)
    before = deepcopy((sample_report, project))
    report = check_project_sample_mapping(project, sample_report)
    assert report["valid"] is True, report
    assert (sample_report, project) == before


@pytest.mark.parametrize("field,value", [("group", "wrong"), ("biological_unit_id", "wrong"), ("donor_id", "wrong"), ("cell_type", "wrong")])
def test_project_and_sheet_identity_conflicts_are_rejected(field, value):
    sample_report = accepted()
    project = project_from(sample_report)
    project["samples"][0][field] = value
    report = check_project_sample_mapping(project, sample_report)
    assert report["valid"] is False, report


def test_project_mapping_uses_stable_id_not_display_name():
    sample_report = accepted(meta=[metadata(display_name="misleading_name")])
    project = project_from(sample_report)
    project["samples"][0]["sample_id"] = "misleading_name"
    report = check_project_sample_mapping(project, sample_report)
    assert report["valid"] is False, report


def test_project_mapping_rejects_a_duplicate_canonical_sample():
    sample_report = accepted()
    project = project_from(sample_report)
    project["samples"].append(deepcopy(project["samples"][0]))
    report = check_project_sample_mapping(project, sample_report)
    assert report["valid"] is False, report


def test_project_mapping_rejects_samples_omitted_from_the_project():
    sample_report = accepted([reads("S1"), reads("S2")], [metadata("S1"), metadata("S2")])
    project = project_from(sample_report)
    project["samples"].pop()
    report = check_project_sample_mapping(project, sample_report)
    assert report["valid"] is False, report


def test_samplesheet_project_check_is_read_only_and_does_not_add_comparisons(tmp_path):
    upstream = write_sheet(tmp_path / "reads.csv", [reads("S1"), reads("S2")])
    meta = write_sheet(tmp_path / "metadata.csv", [metadata("S1", "control"), metadata("S2", "drug")])
    base = accepted([reads("S1"), reads("S2")], [metadata("S1", "control"), metadata("S2", "drug")])
    project = project_from(base)
    project["contrasts"] = [{
        "contrast_id": "chosen_only",
        "display_name": "Drug versus control",
        "kind": "group_difference",
        "sample_ids": ["S1", "S2"],
        "effect_direction": "numerator_minus_denominator",
        "model_formula": None,
        "paired": False,
        "numerator_sample_ids": ["S2"],
        "denominator_sample_ids": ["S1"],
    }]
    project_path = tmp_path / "project.json"
    project_path.write_text(json.dumps(project), encoding="utf-8")
    before = {path.name: path.read_bytes() for path in (upstream, meta, project_path)}
    report = validate_samplesheets(upstream, meta, project_path=project_path)
    assert report["valid"] is True, report
    assert report["contrasts_generated"] == 0
    assert {path.name: path.read_bytes() for path in (upstream, meta, project_path)} == before


def test_invalid_samplesheet_format_returns_validation_failure():
    rejected(input_format="automatically-guess-anything")


def test_empty_input_and_metadata_cannot_form_a_valid_project():
    rejected([], [])


@pytest.mark.parametrize("bad", [False, True])
def test_samplesheet_cli_exit_code_and_json_report_are_consistent(tmp_path, bad):
    upstream = write_sheet(tmp_path / "reads.csv", [reads()])
    meta = write_sheet(tmp_path / "metadata.csv", [metadata(condition="" if bad else "drug")])
    env = dict(os.environ, PYTHONPATH=str(Path(__file__).resolve().parents[2]), PYTHONDONTWRITEBYTECODE="1")
    result = subprocess.run(
        [sys.executable, "-B", "-m", "katana_next.contracts", "samplesheets", "--input", str(upstream), "--metadata", str(meta)],
        cwd=tmp_path, env=env, capture_output=True, text=True, encoding="utf-8", timeout=20,
    )
    assert result.returncode == (2 if bad else 0), result.stderr
    report = json.loads(result.stdout)
    assert report["valid"] is not bad
    assert report["contrasts_generated"] == report["merges_executed"] == report["scientific_tools_executed"] == 0
    assert report["execution_authorized"] is False


def test_diagnostic_order_is_stable_across_python_hash_seeds(tmp_path):
    upstream = write_sheet(tmp_path / "reads.csv", [reads("S1")])
    meta = write_sheet(tmp_path / "metadata.csv", [metadata("S1"), metadata("S2"), metadata("S3"), metadata("S4")])
    reports = []
    for seed in ("1", "2", "3"):
        env = dict(
            os.environ, PYTHONPATH=str(Path(__file__).resolve().parents[2]),
            PYTHONDONTWRITEBYTECODE="1", PYTHONHASHSEED=seed,
        )
        result = subprocess.run(
            [sys.executable, "-B", "-m", "katana_next.contracts", "samplesheets", "--input", str(upstream), "--metadata", str(meta)],
            cwd=tmp_path, env=env, capture_output=True, text=True, encoding="utf-8", timeout=20,
        )
        assert result.returncode == 2, result.stderr
        reports.append(json.loads(result.stdout))
    assert reports[0] == reports[1] == reports[2]


@pytest.mark.parametrize(
    "upstream,meta",
    [
        (None, [metadata()]),
        ("sample,fastq_1", [metadata()]),
        (["not-a-row"], [metadata()]),
        ([{"sample": "S1", "fastq_1": 42, "fastq_2": "", "strandedness": "reverse"}], [metadata()]),
        ([reads()], [{"sample": "S1", "condition": ["drug"], "biological_unit_id": "B1", "donor_id": ""}]),
    ],
)
def test_malformed_direct_inputs_report_errors_instead_of_crashing(upstream, meta):
    report = validate_sample_rows(upstream, meta)
    assert report["valid"] is False, report
    assert any(item["severity"] == "error" for item in report["issues"])
