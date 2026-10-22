"""User-selected comparison preflight, independent of scientific computation.

These synthetic experiments test experimental-unit and model decisions.  No test
fits a statistical model, generates all pairs, imports a GUI, or launches a tool.
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

from katana_next.contracts.selected_contrasts import (
    validate_contrast_rows,
    validate_contrastsheets,
)


def sample(sample_id, condition, *, unit=None, donor=None, **extra):
    return {
        "sample": sample_id,
        "condition": condition,
        "biological_unit_id": unit if unit is not None else f"bio_{sample_id}",
        "donor_id": donor if donor is not None else f"donor_{sample_id}",
        **extra,
    }


def independent(*, third_condition=False):
    rows = [
        sample("C2", "control"),
        sample("T1", "treatment"),
        sample("C1", "control"),
        sample("T2", "treatment"),
    ]
    if third_condition:
        rows += [sample("R1", "rescue"), sample("R2", "rescue")]
    return rows


def comparison(contrast="selected", treatment="treatment", control="control", **extra):
    return {"contrast": contrast, "treatment": treatment, "control": control, **extra}


def paired(n=3, *, pairing="donor_id", unmatched=False):
    rows = []
    for index in range(1, n + 1):
        for condition, prefix in (("control", "C"), ("treatment", "T")):
            kwargs = {"donor": f"D{index}"}
            if pairing == "biological_unit_id":
                kwargs["unit"] = f"B{index}"
            rows.append(sample(f"{prefix}{index}", condition, **kwargs))
    if unmatched:
        rows.append(sample("T_extra", "treatment", donor="D_extra"))
    return rows


def invariant(report):
    assert report["execution_authorized"] is False
    assert report["scientific_tools_executed"] == 0
    assert report["statistical_models_fitted"] == 0
    assert report["all_pairs_generated"] is False
    assert isinstance(report["issues"], list)
    for issue in report["issues"]:
        assert {"code", "path", "message", "severity"} <= issue.keys()
        assert issue["severity"] in {"warning", "error", "info"}
        assert isinstance(issue["code"], str) and issue["code"]
        assert isinstance(issue["message"], str) and issue["message"]
    return report


def accepted(metadata=None, comparisons=None, **kwargs):
    report = invariant(validate_contrast_rows(
        independent() if metadata is None else metadata,
        [comparison()] if comparisons is None else comparisons,
        **kwargs,
    ))
    assert report["valid"] is True, report
    assert not [item for item in report["issues"] if item["severity"] == "error"], report
    return report


def rejected(metadata=None, comparisons=None, **kwargs):
    report = invariant(validate_contrast_rows(
        independent() if metadata is None else metadata,
        [comparison()] if comparisons is None else comparisons,
        **kwargs,
    ))
    assert report["valid"] is False, report
    assert any(item["severity"] == "error" for item in report["issues"]), report
    return report


def first(report):
    return report["selected_contrasts"][0]


def write_sheet(path, rows, columns=None):
    columns = columns or list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    return path


def test_only_requested_pair_is_selected_from_three_conditions():
    report = accepted(independent(third_condition=True))
    assert report["contrast_count"] == 1
    assert len(report["selected_contrasts"]) == 1
    selected = first(report)
    assert selected["contrast_id"] == "selected"
    assert selected["treatment"] == "treatment"
    assert selected["control"] == "control"
    assert selected["numerator_sample_ids"] == ["T1", "T2"]
    assert selected["denominator_sample_ids"] == ["C2", "C1"]
    assert selected["effect_direction"] == "treatment_minus_control"
    assert selected["design"]["row_count"] == 4


def test_selected_row_order_and_metadata_sample_order_are_preserved():
    report = accepted(independent(third_condition=True), [
        comparison("rescue_vs_control", "rescue", "control"),
        comparison("treatment_vs_control"),
    ])
    assert [item["contrast_id"] for item in report["selected_contrasts"]] == [
        "rescue_vs_control", "treatment_vs_control"
    ]
    assert report["selected_contrasts"][0]["denominator_sample_ids"] == ["C2", "C1"]
    assert report["contrast_count"] == 2


def test_unselected_condition_is_not_forced_into_selected_replication_or_model_checks():
    rows = independent() + [sample("unselected_singleton", "not_selected")]
    selected = first(accepted(rows))
    assert selected["design"]["row_count"] == 4
    assert "unselected_singleton" not in selected["numerator_sample_ids"] + selected["denominator_sample_ids"]


def test_imported_numeric_condition_labels_remain_literal_stable_labels():
    rows = [sample("S1", "0"), sample("S2", "0"), sample("S3", "1"), sample("S4", "1")]
    selected = first(accepted(rows, [comparison(treatment="1", control="0")],
                              route="import_only", intent="import"))
    assert selected["treatment"] == "1" and selected["control"] == "0"
    assert selected["design"]["contrast_estimable"] is None


def test_native_route_rejects_unsupported_condition_ids_without_silent_renaming():
    rows = [sample("S1", "0"), sample("S2", "0"), sample("S3", "1"), sample("S4", "1")]
    contrasts = [comparison(treatment="1", control="0")]
    original = deepcopy((rows, contrasts))
    rejected(rows, contrasts)
    assert (rows, contrasts) == original


def test_independent_design_has_actual_rank_and_positive_residual_degrees_of_freedom():
    design = first(accepted())["design"]
    assert design["requested_formula"] == "~ condition"
    assert design["terms"] == ["condition"]
    assert design["row_count"] == 4
    assert design["column_count"] == len(design["matrix_columns"]) == 2
    assert design["rank"] == 2
    assert design["residual_df"] == 2
    assert design["contrast_estimable"] is True
    assert design["paired_by"] is None
    assert design["pair_map"] == {}


def test_labels_are_not_inferred_from_ids_or_filename_like_strings():
    rows = [
        sample("KO_r1.fastq.gz", "drug"), sample("KO_r2.fastq.gz", "drug"),
        sample("WT_r1.fastq.gz", "vehicle"), sample("WT_r2.fastq.gz", "vehicle"),
    ]
    selected = first(accepted(rows, [comparison(treatment="vehicle", control="drug")]))
    assert selected["numerator_sample_ids"] == ["WT_r1.fastq.gz", "WT_r2.fastq.gz"]
    assert selected["denominator_sample_ids"] == ["KO_r1.fastq.gz", "KO_r2.fastq.gz"]


def test_explicit_reciprocal_pair_is_kept_with_a_warning():
    report = accepted(comparisons=[comparison("forward"), comparison("reverse", "control", "treatment")])
    assert report["contrast_count"] == 2
    assert first(report)["numerator_sample_ids"] == report["selected_contrasts"][1]["denominator_sample_ids"]
    assert any(item["severity"] == "warning" for item in report["issues"]), report


def test_same_directed_pair_cannot_be_disguised_with_different_ids():
    rejected(comparisons=[comparison("one"), comparison("two")])


def test_same_contrast_id_cannot_name_two_different_comparisons():
    rejected(independent(third_condition=True), [comparison("same"), comparison("same", "rescue", "control")])


@pytest.mark.parametrize("field,value", [("contrast", ""), ("treatment", ""), ("control", ""),
                                          ("treatment", "unknown"), ("control", "unknown")])
def test_missing_or_unknown_comparison_identity_is_rejected(field, value):
    row = comparison()
    row[field] = value
    rejected(comparisons=[row])


def test_self_comparison_is_rejected():
    rejected(comparisons=[comparison(treatment="control", control="control")])


def test_empty_selection_never_generates_all_pairs():
    report = rejected(independent(third_condition=True), [])
    assert report["contrast_count"] == 0
    assert report["selected_contrasts"] == []


def test_one_independent_unit_in_an_arm_is_insufficient_for_compute():
    rejected(independent()[:-1])


@pytest.mark.parametrize("field", ["sample", "condition", "biological_unit_id"])
def test_required_sample_identity_cannot_be_blank(field):
    rows = independent()
    rows[0][field] = ""
    rejected(rows)


def test_duplicate_sample_id_is_rejected_even_if_conditions_differ():
    rows = independent()
    rows[1]["sample"] = rows[0]["sample"]
    rejected(rows)


def test_duplicate_biological_unit_within_arm_is_not_an_independent_replicate():
    rows = independent()
    rows[2]["biological_unit_id"] = rows[0]["biological_unit_id"]
    rejected(rows)


def test_two_samples_from_one_donor_within_arm_need_a_formal_nested_model():
    rows = independent()
    rows[2]["donor_id"] = rows[0]["donor_id"]
    rejected(rows)


def test_explicit_technical_repeat_is_not_accepted_as_an_independent_replicate():
    rows = independent()
    rows[2]["technical_replicate_of"] = rows[0]["sample"]
    rejected(rows)


def test_technical_repeat_link_must_reference_an_existing_sample():
    rows = independent()
    rows[2]["technical_replicate_of"] = "not_registered"
    rejected(rows)


def test_unknown_donors_do_not_invent_pairing_for_independent_samples():
    rows = independent()
    for row in rows:
        row["donor_id"] = ""
    selected = first(accepted(rows))
    assert selected["design"]["paired_by"] is None
    assert selected["design"]["pair_map"] == {}


def test_shared_donors_cannot_be_silently_treated_as_independent():
    rejected(paired())


def test_shared_biological_units_cannot_be_silently_treated_as_independent():
    rejected(paired(pairing="biological_unit_id"))


@pytest.mark.parametrize("pairing", ["donor_id", "biological_unit_id"])
def test_explicit_pairing_builds_complete_subject_blocked_model(pairing):
    selected = first(accepted(paired(pairing=pairing), route="deseq2_adapter",
                              design="~ subject + condition", paired_by=pairing))
    design = selected["design"]
    assert design["paired_by"] == pairing
    assert len(design["pair_map"]) == 3
    assert design["row_count"] == 6
    assert design["rank"] == design["column_count"] == 4
    assert design["residual_df"] == 2
    assert design["contrast_estimable"] is True
    serialized = json.dumps(design["pair_map"])
    for index in range(1, 4):
        assert f"C{index}" in serialized and f"T{index}" in serialized


def test_two_complete_pairs_have_one_residual_degree_of_freedom():
    design = first(accepted(paired(n=2), route="deseq2_adapter",
                            design="~ subject + condition", paired_by="donor_id"))["design"]
    assert design["row_count"] == 4
    assert design["column_count"] == design["rank"] == 3
    assert design["residual_df"] == 1


def test_one_pair_is_not_sufficient_even_when_subject_id_is_known():
    rejected(paired(n=1), route="deseq2_adapter", design="~ subject + condition", paired_by="donor_id")


def test_unmatched_subject_is_rejected_without_dropping_its_sample():
    rows = paired(unmatched=True)
    original = deepcopy(rows)
    rejected(rows, route="deseq2_adapter", design="~ subject + condition", paired_by="donor_id")
    assert rows == original


def test_missing_subject_identifier_is_not_imputed():
    rows = paired()
    rows[0]["donor_id"] = ""
    rejected(rows, route="deseq2_adapter", design="~ subject + condition", paired_by="donor_id")


def test_repeated_subject_within_one_arm_cannot_be_fixed_by_pairing_flag():
    rows = paired()
    rows.append(sample("C_extra", "control", donor="D1"))
    rejected(rows, route="deseq2_adapter", design="~ subject + condition", paired_by="donor_id")


def test_subject_design_requires_explicit_pairing_key():
    rejected(paired(), route="deseq2_adapter", design="~ subject + condition")


def test_pairing_key_requires_subject_in_the_design():
    rejected(paired(), route="deseq2_adapter", design="~ condition", paired_by="donor_id")


def test_bulk_complete_route_does_not_claim_a_paired_engine():
    rejected(paired(), design="~ subject + condition", paired_by="donor_id")


def test_pairing_field_cannot_be_a_convenient_arbitrary_metadata_column():
    rejected(paired(), route="deseq2_adapter", design="~ subject + condition", paired_by="batch")


@pytest.mark.parametrize("paired_value", ["false", "False", "FALSE"])
def test_row_paired_false_conflicts_with_explicit_subject_model(paired_value):
    rejected(paired(), [comparison(paired=paired_value)], route="deseq2_adapter",
             design="~ subject + condition", paired_by="donor_id")


def test_row_paired_true_requires_an_explicit_pairing_key_and_subject_model():
    rejected(comparisons=[comparison(paired="true")])


def test_row_level_design_and_pairing_overrides_are_explicit():
    selected = first(accepted(paired(), [comparison(design="~ subject + condition", paired_by="donor_id", paired="true")],
                              route="deseq2_adapter"))
    assert selected["design"]["paired_by"] == "donor_id"
    assert selected["design"]["requested_formula"] == "~ subject + condition"


@pytest.mark.parametrize("paired_value", ["false", "0", "no", "NO"])
def test_native_false_boolean_aliases_do_not_request_pairing(paired_value):
    selected = first(accepted(comparisons=[comparison(paired=paired_value)]))
    assert selected["design"]["declared_paired"] is False
    assert selected["design"]["paired_by"] is None


@pytest.mark.parametrize("paired_value", ["true", "1", "yes", "YES"])
def test_native_true_boolean_aliases_require_and_preserve_explicit_pairing(paired_value):
    selected = first(accepted(paired(), [comparison(paired=paired_value)], route="deseq2_adapter",
                              design="~ subject + condition", paired_by="donor_id"))
    assert selected["design"]["declared_paired"] is True
    assert selected["design"]["paired_by"] == "donor_id"


def test_batch_adjusted_adapter_checks_a_crossed_design_without_fitting():
    rows = independent()
    for row, batch in zip(rows, ["B1", "B1", "B2", "B2"]):
        row["batch"] = batch
    selected = first(accepted(rows, route="deseq2_adapter", design="~ batch + condition"))
    design = selected["design"]
    assert design["row_count"] == 4
    assert design["rank"] == design["column_count"] == 3
    assert design["residual_df"] == 1
    assert design["contrast_estimable"] is True
    assert selected["methods"]
    assert all(method["profile_status"] == "adapter_present_route_unwired" for method in selected["methods"])
    assert all(method["runtime_ready"] is False for method in selected["methods"])


def test_bulk_complete_route_does_not_silently_drop_requested_batch_adjustment():
    rows = independent()
    for row, batch in zip(rows, ["B1", "B1", "B2", "B2"]):
        row["batch"] = batch
    rejected(rows, design="~ batch + condition")


def test_batch_per_condition_is_rank_deficient_and_cannot_be_dropped():
    rows = independent()
    for row in rows:
        row["batch"] = "B1" if row["condition"] == "control" else "B2"
    rejected(rows, route="deseq2_adapter", design="~ batch + condition")


@pytest.mark.parametrize("batch_values", [["B1", "B1", "B1", "B1"], ["B1", "B1", "B2", ""]])
def test_batch_single_level_or_missing_value_is_rejected(batch_values):
    rows = independent()
    for row, value in zip(rows, batch_values):
        row["batch"] = value
    rejected(rows, route="deseq2_adapter", design="~ batch + condition")


def test_design_with_no_residual_degrees_of_freedom_is_rejected():
    rows = independent()
    for index, row in enumerate(rows):
        row["batch"] = f"B{index}"
    rejected(rows, route="deseq2_adapter", design="~ batch + condition")


@pytest.mark.parametrize("formula", ["~ age + condition", "~ covariate_age + condition", "~ batch * condition",
                                     "~ batch:condition", "~ condition + I(age**2)", "~ 0 + condition",
                                     "~ condition + (1|subject)", "~ unknown", "condition", "", "~ condition + condition"])
def test_unsupported_compute_models_never_fall_back_to_a_simple_comparison(formula):
    rows = independent()
    for index, row in enumerate(rows):
        row["batch"] = "B1" if index < 2 else "B2"
        row["covariate_age"] = str(20 + index)
    rejected(rows, route="deseq2_adapter", design=formula)


def test_unknown_model_column_in_contrast_sheet_is_not_silently_ignored():
    rejected(comparisons=[comparison(model_formula="~ batch + condition")])


def test_import_only_retains_opaque_native_formula_without_compute_eligibility():
    formula = "~ batch * condition + splines::ns(age, df=3)"
    selected = first(accepted(route="import_only", intent="import", design=formula))
    design = selected["design"]
    assert design["requested_formula"] == formula
    assert design["contrast_estimable"] is None
    assert design["rank"] is None
    assert design["residual_df"] is None
    assert selected["numerator_sample_ids"] == ["T1", "T2"]
    assert selected["denominator_sample_ids"] == ["C2", "C1"]
    assert all(method["runtime_ready"] is False for method in selected["methods"])


def test_import_only_does_not_require_a_new_native_job_or_compute_replication():
    rows = [sample("C", "control", unit=""), sample("T", "treatment", unit="")]
    report = accepted(rows, route="import_only", intent="import", design="native model recorded externally")
    assert report["contrast_count"] == 1
    assert first(report)["design"]["contrast_estimable"] is None


def test_native_formula_is_never_evaluated_as_code(tmp_path):
    marker = tmp_path / "must_not_be_created.txt"
    formula = f"__import__('pathlib').Path({str(marker)!r}).write_text('unexpected')"
    accepted(route="import_only", intent="import", design=formula)
    rejected(route="deseq2_adapter", design=formula)
    assert not marker.exists()


def test_import_only_route_cannot_authorize_new_compute():
    rejected(route="import_only", intent="compute")


@pytest.mark.parametrize("paired_value", ["maybe", "2", "NaN", "?"])
def test_import_cannot_record_a_malformed_typed_pairing_flag(paired_value):
    rejected(comparisons=[comparison(paired=paired_value)], route="import_only", intent="import")


def test_import_row_pairing_key_is_validated_like_the_global_option():
    rejected(comparisons=[comparison(paired_by="banana")], route="import_only", intent="import")


def test_import_preserves_valid_native_pairing_fields_without_certifying_the_native_model():
    raw = comparison(paired="yes", paired_by="donor_id", design="~ subject * condition + native_term")
    selected = first(accepted(comparisons=[raw], route="import_only", intent="import"))
    assert selected["source_contrast_fields"] == raw
    assert selected["design"]["declared_paired"] is True
    assert selected["design"]["paired_by"] == "donor_id"
    assert selected["design"]["requested_formula"] == raw["design"]
    assert selected["design"]["contrast_estimable"] is None
    assert selected["design"]["pair_map"] == {}
    assert selected["fit_sample_ids"] is None
    assert selected["fit_scope"] == "native_fit_unverified"
    assert selected["selected_sample_ids"] == ["C2", "T1", "C1", "T2"]
    selected["source_contrast_fields"]["paired"] = "modified_return_value"
    assert raw["paired"] == "yes"


def test_import_selection_does_not_invent_a_native_fit_universe():
    selected = first(accepted(independent(third_condition=True), route="import_only", intent="import"))
    assert selected["selected_sample_ids"] == ["C2", "T1", "C1", "T2"]
    assert selected["fit_sample_ids"] is None
    assert selected["fit_scope"] == "native_fit_unverified"
    assert selected["design"]["declared_paired"] is None


@pytest.mark.parametrize("assay", ["single_cell_rna", "single_nucleus_rna"])
def test_bulk_complete_does_not_advertise_raw_single_cell_compute(assay):
    rejected(assay_mode=assay)


@pytest.mark.parametrize("assay", ["single_cell_rna", "single_nucleus_rna"])
def test_single_cell_import_mapping_does_not_relabel_units_as_valid_bulk_compute(assay):
    report = accepted(assay_mode=assay, intent="import", route="import_only")
    assert first(report)["design"]["contrast_estimable"] is None
    assert all(method["runtime_ready"] is False for method in first(report)["methods"])


def test_multiple_cell_types_require_explicit_selection_even_for_import():
    rows = [dict(row, cell_type=cell_type) for row, cell_type in zip(independent(), ["Tcell", "Tcell", "Bcell", "Bcell"])]
    rejected(rows, route="import_only", intent="import", assay_mode="single_cell_rna")


def test_explicit_cell_type_does_not_pool_other_cell_types():
    rows = []
    for cell_type in ("Tcell", "Bcell"):
        for row in independent():
            rows.append(dict(row, sample=f"{row['sample']}_{cell_type}", cell_type=cell_type))
    selected = first(accepted(rows, route="import_only", intent="import", assay_mode="single_cell_rna", cell_type="Tcell"))
    assert selected["numerator_sample_ids"] == ["T1_Tcell", "T2_Tcell"]
    assert selected["denominator_sample_ids"] == ["C2_Tcell", "C1_Tcell"]


def test_selected_cell_type_must_exist():
    rows = [dict(row, cell_type="Tcell") for row in independent()]
    rejected(rows, route="import_only", intent="import", assay_mode="single_cell_rna", cell_type="unknown")


def test_native_effect_conventions_remain_method_specific_and_not_executed():
    methods = first(accepted())["methods"]
    assert methods
    normalized = {row["method_id"].lower(): row for row in methods}
    for fragment in ("deseq", "rmats", "suppa", "dexseq"):
        found = [row for name, row in normalized.items() if fragment in name]
        assert found, methods
        assert all(isinstance(row["native_effect"], str) and row["native_effect"] for row in found)
    assert all(row["runtime_ready"] is False for row in methods)
    assert len({row["native_effect"] for row in methods}) > 1
    assert normalized["deseq2"]["native_effect"] == "log2(treatment/control)"
    assert normalized["rmats"]["native_effect"] == "group1_minus_group2"
    assert normalized["suppa"]["native_effect"] == "second_minus_first; adapter_control_then_treatment"
    assert normalized["dexseq_exon"]["native_effect"] == "native_header_authoritative_control_over_treatment"
    assert "qvalue_unsigned" in normalized["dexseq_dtu_stager"]["native_effect"]


@pytest.mark.parametrize("metadata_rows,contrast_rows", [
    (None, [comparison()]), ("not rows", [comparison()]), ({"sample": "S1"}, [comparison()]),
    ([None], [comparison()]), (["bad"], [comparison()]), (independent(), None),
    (independent(), "not rows"), (independent(), [None]), (independent(), ["bad"]),
    (independent(), [{"contrast": [], "treatment": "treatment", "control": "control"}]),
])
def test_malformed_rows_return_structured_diagnostics(metadata_rows, contrast_rows):
    report = invariant(validate_contrast_rows(metadata_rows, contrast_rows))
    assert report["valid"] is False
    assert any(item["severity"] == "error" for item in report["issues"])


def test_inputs_are_not_mutated_and_reports_are_independent_copies():
    rows, contrasts = independent(), [comparison()]
    before = deepcopy((rows, contrasts))
    report = accepted(rows, contrasts)
    assert (rows, contrasts) == before
    first(report)["numerator_sample_ids"].append("not_in_input")
    assert (rows, contrasts) == before


def test_csv_entry_point_matches_in_memory_selection(tmp_path):
    metadata_path = write_sheet(tmp_path / "metadata.csv", independent())
    contrast_path = write_sheet(tmp_path / "contrasts.csv", [comparison()])
    before = {path: path.read_bytes() for path in (metadata_path, contrast_path)}
    report = invariant(validate_contrastsheets(metadata_path, contrast_path))
    assert report["valid"] is True, report
    assert first(report)["numerator_sample_ids"] == ["T1", "T2"]
    assert all(path.read_bytes() == data for path, data in before.items())


@pytest.mark.parametrize("which", ["metadata", "contrasts"])
def test_missing_sheet_has_structured_error(tmp_path, which):
    metadata_path = write_sheet(tmp_path / "metadata.csv", independent())
    contrast_path = write_sheet(tmp_path / "contrasts.csv", [comparison()])
    missing = tmp_path / "does_not_exist.csv"
    report = invariant(validate_contrastsheets(missing if which == "metadata" else metadata_path,
                                              missing if which == "contrasts" else contrast_path))
    assert report["valid"] is False, report


@pytest.mark.parametrize("which", ["metadata", "contrasts"])
def test_csv_duplicate_headers_are_rejected(tmp_path, which):
    metadata_path = write_sheet(tmp_path / "metadata.csv", independent())
    contrast_path = write_sheet(tmp_path / "contrasts.csv", [comparison()])
    if which == "metadata":
        metadata_path.write_text("sample,sample,condition,biological_unit_id\nS1,S2,control,B1\n", encoding="utf-8")
    else:
        contrast_path.write_text("contrast,treatment,control,control\nC,treatment,control,control\n", encoding="utf-8")
    report = invariant(validate_contrastsheets(metadata_path, contrast_path))
    assert report["valid"] is False, report


@pytest.mark.parametrize("line", ["selected,treatment\n", "selected,treatment,control,extra\n"])
def test_csv_row_width_errors_are_not_silently_truncated_or_filled(tmp_path, line):
    metadata_path = write_sheet(tmp_path / "metadata.csv", independent())
    contrast_path = tmp_path / "contrasts.csv"
    contrast_path.write_text("contrast,treatment,control\n" + line, encoding="utf-8")
    report = invariant(validate_contrastsheets(metadata_path, contrast_path))
    assert report["valid"] is False, report


@pytest.mark.parametrize("is_valid,expected_exit", [(True, 0), (False, 2)])
def test_cli_returns_strict_json_and_correct_exit_status(tmp_path, is_valid, expected_exit):
    metadata_path = write_sheet(tmp_path / "metadata.csv", independent())
    contrast_path = write_sheet(tmp_path / "contrasts.csv", [comparison(control="control" if is_valid else "missing")])
    env = dict(os.environ, PYTHONPATH=str(Path(__file__).resolve().parents[2]), PYTHONDONTWRITEBYTECODE="1")
    completed = subprocess.run([
        sys.executable, "-B", "-m", "katana_next.contracts", "contrasts",
        "--metadata", str(metadata_path), "--contrasts", str(contrast_path),
        "--route", "bulk_complete", "--design", "~ condition", "--intent", "compute", "--assay", "bulk_rna",
    ], capture_output=True, text=True, encoding="utf-8", env=env, timeout=30)
    assert completed.returncode == expected_exit, (completed.stdout, completed.stderr)
    report = invariant(json.loads(completed.stdout))
    assert report["valid"] is is_valid
    assert not completed.stderr.strip(), completed.stderr


def test_contract_module_does_not_eagerly_import_ui_or_scientific_fitting_modules():
    env = dict(os.environ, PYTHONPATH=str(Path(__file__).resolve().parents[2]), PYTHONDONTWRITEBYTECODE="1")
    code = (
        "import json, sys; from katana_next.contracts.selected_contrasts import validate_contrast_rows; "
        "blocked=('PySide6', 'PyQt5', 'PyQt6', 'matplotlib', 'statsmodels', 'scanpy', 'rpy2'); "
        "print(json.dumps([name for name in blocked if name in sys.modules]))"
    )
    completed = subprocess.run([sys.executable, "-B", "-c", code], capture_output=True,
                               text=True, encoding="utf-8", env=env, timeout=30)
    assert completed.returncode == 0, completed.stderr
    assert json.loads(completed.stdout) == []
