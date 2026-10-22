"""Independent metadata/provenance contract tests; no scientific tools are run.

These tests assert externally useful rejection and preservation behavior.  They
do not confer scientific acceptance on the synthetic examples or a new package.
"""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from katana_next.contracts import validate_document, validate_file


SOURCE = Path(__file__).resolve().parents[2]
EXAMPLES = SOURCE / "examples" / "contracts" / "v1"
EXAMPLE_NAMES = (
    "bulk-native-import.project.json",
    "missing-result.project.json",
    "single-cell-pseudobulk.project.json",
)


def example(name=EXAMPLE_NAMES[0]):
    return json.loads((EXAMPLES / name).read_text(encoding="utf-8"))


def issue_codes(report):
    return {issue["code"] for issue in report["issues"]}


def rejected(document, code=None, **kwargs):
    report = validate_document(document, base_dir=EXAMPLES, **kwargs)
    assert report["valid"] is False, report
    assert report["execution_authorized"] is False
    assert report["scientific_tools_executed"] == 0
    if code:
        assert code in issue_codes(report), report
    return report


def accepted(document, **kwargs):
    report = validate_document(document, base_dir=EXAMPLES, **kwargs)
    assert report["valid"] is True, report
    assert report["execution_authorized"] is False
    assert report["scientific_tools_executed"] == 0
    return report


def computed_project():
    """A complete *metadata claim*, not evidence of actual scientific execution."""
    doc = example()
    doc["intent"] = "analyze"
    method = doc["methods"][0]
    method.update(operation="compute", statistics_origin="computed", source_sha256="a" * 64)
    run = doc["runs"][0]
    run["execution_kind"] = "new_computation"
    run["software"]["source_commit"] = "b" * 40
    run["provenance"].update(
        source_uri=None, source_schema=None, source_version=None,
        source_artifact_id=None, source_completion="not_applicable",
    )
    doc["references"][0].update(identity_status="known", taxon_id=10090, genome_build="GRCm39")
    doc["results"][0]["statistics_origin"] = "computed"
    return doc


def tree_identity(root):
    return {
        path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in root.rglob("*") if path.is_file()
    }


def cli(*arguments):
    env = dict(os.environ, PYTHONPATH=str(SOURCE), PYTHONDONTWRITEBYTECODE="1")
    return subprocess.run(
        [sys.executable, "-B", "-m", "katana_next.contracts", *map(str, arguments)],
        cwd=SOURCE, env=env, capture_output=True, text=True,
        encoding="utf-8", timeout=30, check=False,
    )


@pytest.mark.parametrize("name", EXAMPLE_NAMES)
def test_shipped_examples_validate_and_bind_local_bytes(name):
    report = validate_file(EXAMPLES / name, verify_hashes=True)
    assert report["valid"], report
    assert report["execution_authorized"] is False
    assert report["scientific_tools_executed"] == 0
    assert not [i for i in report["issues"] if i["severity"] == "error"]
    if name == "missing-result.project.json":
        assert report["file_verification"]["local_files_checked"] == 0
    else:
        assert report["file_verification"]["hashes_checked"] >= 1


def test_document_validation_is_nonmutating():
    doc = example()
    original = deepcopy(doc)
    accepted(doc, verify_hashes=True)
    assert doc == original


def test_duplicate_display_names_do_not_redefine_sample_identity():
    doc = example()
    for sample in doc["samples"]:
        sample["display_name"] = "display name may repeat"
    accepted(doc)


def test_duplicate_sample_identity_is_rejected():
    doc = example()
    doc["samples"].append(deepcopy(doc["samples"][0]))
    rejected(doc, "DUPLICATE_ID")


def test_unknown_comparison_sample_is_rejected():
    doc = example()
    doc["contrasts"][0]["sample_ids"].append("sample_unknown")
    doc["contrasts"][0]["numerator_sample_ids"].append("sample_unknown")
    rejected(doc, "UNKNOWN_SAMPLE")


def test_comparison_arms_cannot_share_a_sample():
    doc = example()
    contrast = doc["contrasts"][0]
    contrast["numerator_sample_ids"].append(contrast["denominator_sample_ids"][0])
    rejected(doc, "CONTRAST_ARM_OVERLAP")


def test_comparison_scope_cannot_hide_an_arm_member():
    doc = example()
    doc["contrasts"][0]["sample_ids"].pop()
    rejected(doc, "CONTRAST_SCOPE_MISMATCH")


@pytest.mark.parametrize("direction", [None, "alphabetic", "denominator_minus_numerator"])
def test_group_effect_direction_must_be_explicit(direction):
    doc = example()
    doc["contrasts"][0]["effect_direction"] = direction
    rejected(doc, "SCHEMA")


def test_paired_design_requires_explicit_complete_one_to_one_mapping():
    doc = example()
    doc["contrasts"][0]["paired"] = True
    rejected(doc, "PAIR_MAPPING_REQUIRED")
    doc["contrasts"][0]["pair_ids"] = {
        "sample_A1": "pair1", "sample_B1": "pair1",
        "sample_A2": "pair2", "sample_B2": "pair2",
    }
    accepted(doc)


@pytest.mark.parametrize("source_completion", ["failed", "cancelled", "unknown"])
def test_completed_import_is_not_completed_source_science(source_completion):
    doc = example()
    doc["runs"][0]["provenance"]["source_completion"] = source_completion
    rejected(doc, "SOURCE_COMPLETION_UNKNOWN")
    result = doc["results"][0]
    result.update(
        status="unverified", row_count=None, tested_count=None,
        significant_count=None, findings="not_assessed",
        reason={"code": "source_unverified", "detail": "Original source science is not completed."},
    )
    accepted(doc)


def test_import_does_not_require_a_new_native_job_or_source_commit():
    doc = example()
    run = doc["runs"][0]
    run["execution_evidence"] = {"status": "not_applicable", "receipt_asset_ids": []}
    assert run["software"]["source_commit"] is None
    assert "native_job_id" not in run
    accepted(doc)


def test_foreign_native_artifact_identity_is_not_an_unbound_project_reference():
    doc = example()
    doc["runs"][0]["provenance"]["source_artifact_id"] = "asset_not_registered"
    rejected(doc, "UNKNOWN_ASSET")


@pytest.mark.parametrize("field,value", [
    ("role", "input"), ("availability", "missing"), ("sha256", None), ("bytes", None),
])
def test_completed_new_compute_receipt_is_not_a_placeholder(field, value):
    doc = computed_project()
    doc["assets"][1][field] = value
    rejected(doc, "COMPLETION_RECEIPT_IDENTITY")


@pytest.mark.parametrize("new_compute", [False, True])
def test_worker_operation_matches_import_or_new_compute(new_compute):
    doc = computed_project() if new_compute else example()
    doc["methods"][0]["operation"] = "import" if new_compute else "compute"
    rejected(doc, "METHOD_EXECUTION_KIND")


@pytest.mark.parametrize("field", ["row_count", "tested_count", "significant_count"])
def test_missing_result_does_not_become_zero(field):
    doc = example("missing-result.project.json")
    doc["results"][0][field] = 0
    rejected(doc, "MISSING_IS_NOT_ZERO")


def test_completed_zero_findings_need_a_nonempty_tested_family():
    doc = example()
    result = doc["results"][0]
    result.update(significant_count=0, findings="zero_findings")
    accepted(doc)
    result["tested_count"] = 0
    rejected(doc, "FINDINGS_REQUIRE_TESTED_UNIVERSE")


@pytest.mark.parametrize("field", ["family_id", "tested_universe_asset_id"])
def test_zero_findings_cannot_omit_tested_background_identity(field):
    doc = example()
    doc["results"][0].update(significant_count=0, findings="zero_findings")
    doc["results"][0][field] = None
    rejected(doc, "FINDINGS_REQUIRE_TESTED_UNIVERSE")


def test_a_missing_tested_universe_cannot_establish_an_assessed_family():
    doc = example()
    doc["assets"][0]["availability"] = "missing"
    rejected(doc, "TESTED_UNIVERSE_IDENTITY")


def test_empty_completed_table_is_not_automatically_a_zero_finding_claim():
    doc = example()
    doc["results"][0].update(
        row_count=0, tested_count=None, significant_count=None,
        findings="not_assessed", family_id=None, tested_universe_asset_id=None,
    )
    accepted(doc)


def test_significant_count_cannot_exceed_tested_count():
    doc = example()
    doc["results"][0]["significant_count"] = 3
    rejected(doc, "INVALID_FINDING_COUNTS")


@pytest.mark.parametrize("field,value", [
    ("measurement_unit", "TPM"),
    ("inference_unit", "donor"),
    ("statistics_origin", "computed"),
])
def test_result_cannot_relabel_method_unit_or_statistical_origin(field, value):
    doc = example()
    doc["results"][0][field] = value
    rejected(doc, "RESULT_METHOD_UNIT_OR_ORIGIN")


def test_output_asset_native_scale_must_match_result_scale():
    doc = example()
    doc["assets"][0]["measurement_unit"] = "PSI_percent"
    rejected(doc, "ASSET_RESULT_UNIT")


def test_an_input_asset_has_an_explicit_measurement_unit_field():
    doc = example("single-cell-pseudobulk.project.json")
    assert doc["assets"][0]["measurement_unit"] == "UMI_count"
    del doc["assets"][0]["measurement_unit"]
    rejected(doc, "SCHEMA")


@pytest.mark.parametrize("assay", ["single_cell_rna", "single_nucleus_rna"])
def test_raw_single_cell_fastq_is_not_a_bulk_upstream_entry(assay):
    doc = example("single-cell-pseudobulk.project.json")
    doc["assay_mode"] = assay
    for sample in doc["samples"]:
        sample["source_assay"] = assay
        sample["aggregation"]["source_assay"] = assay
    doc["input_kind"] = "bulk_fastq"
    rejected(doc, "RAW_SINGLE_CELL_UNSUPPORTED")


def test_pseudobulk_preserves_assay_instead_of_becoming_bulk():
    doc = example("single-cell-pseudobulk.project.json")
    doc["analysis_level"] = "bulk"
    rejected(doc, "ASSAY_LEVEL_MISMATCH")


@pytest.mark.parametrize("field", ["donor_id", "cell_type", "aggregation"])
def test_pseudobulk_requires_donor_celltype_and_source_lineage(field):
    doc = example("single-cell-pseudobulk.project.json")
    doc["samples"][0][field] = None
    rejected(doc, "PSEUDOBULK_LINEAGE_REQUIRED")


def test_pseudobulk_grouping_must_include_donor_and_cell_type():
    doc = example("single-cell-pseudobulk.project.json")
    doc["samples"][0]["aggregation"]["group_by"] = ["sample_id"]
    rejected(doc, "PSEUDOBULK_GROUPING")


def test_active_sc_statistical_model_cannot_use_cells_as_biological_replicates():
    doc = computed_project()
    doc.update(assay_mode="single_cell_rna", analysis_level="cell", input_kind="processed_counts")
    for sample in doc["samples"]:
        sample["source_assay"] = "single_cell_rna"
        sample["cell_type"] = "neuronal_fixture"
        sample["donor_id"] = None
    rejected(doc, "DONOR_INFERENCE_REQUIRED")


def test_imported_cell_group_statistics_remain_cell_group_statistics():
    doc = example()
    doc.update(assay_mode="single_cell_rna", analysis_level="cell")
    for sample in doc["samples"]:
        sample["source_assay"] = "single_cell_rna"
    doc["methods"][0]["inference_unit"] = "cell"
    doc["results"][0]["inference_unit"] = "cell"
    original = deepcopy(doc)
    accepted(doc)
    assert doc == original
    assert doc["results"][0]["inference_unit"] == "cell"


def test_new_cell_unit_contract_is_not_silently_rewritten_to_donor():
    doc = computed_project()
    doc.update(assay_mode="single_cell_rna", analysis_level="cell", input_kind="processed_counts")
    for sample in doc["samples"]:
        sample["source_assay"] = "single_cell_rna"
    doc["methods"][0]["inference_unit"] = "cell"
    doc["results"][0]["inference_unit"] = "cell"
    report = accepted(doc)
    assert report["execution_authorized"] is False
    assert doc["results"][0]["inference_unit"] == "cell"


def test_declared_donor_inference_requires_donor_mapping():
    doc = computed_project()
    doc.update(assay_mode="single_cell_rna", analysis_level="cell", input_kind="processed_counts")
    for sample in doc["samples"]:
        sample["source_assay"] = "single_cell_rna"
    doc["methods"][0]["inference_unit"] = "donor"
    doc["results"][0]["inference_unit"] = "donor"
    rejected(doc, "DONOR_INFERENCE_REQUIRED")
    for sample in doc["samples"]:
        sample["donor_id"] = "donor_" + sample["sample_id"]
    accepted(doc)


def test_new_completed_compute_requires_receipt_and_source_identity():
    doc = computed_project()
    accepted(doc)
    doc["runs"][0]["execution_evidence"] = {"status": "reported", "receipt_asset_ids": []}
    rejected(doc, "COMPLETION_EVIDENCE_REQUIRED")
    doc = computed_project()
    doc["runs"][0]["software"]["source_commit"] = None
    rejected(doc, "SOURCE_IDENTITY_REQUIRED")


@pytest.mark.parametrize("intent", ["import", "report"])
def test_import_and_report_intents_do_not_authorize_new_scientific_compute(intent):
    doc = computed_project()
    doc["intent"] = intent
    rejected(doc, "INTENT_COMPUTE_CONFLICT")


def test_missing_new_compute_input_is_not_ready():
    doc = computed_project()
    doc["assets"][0]["availability"] = "missing"
    rejected(doc, "UNAVAILABLE_COMPUTE_INPUT")


def test_unresolved_reference_may_be_imported_but_not_newly_computed():
    accepted(example())
    doc = computed_project()
    doc["references"][0].update(identity_status="unresolved", taxon_id=None, genome_build=None)
    rejected(doc, "UNRESOLVED_COMPUTE_REFERENCE")


def test_run_ancestry_cannot_form_a_cycle():
    doc = example()
    second = deepcopy(doc["runs"][0])
    second["run_id"] = "run_second"
    second["parent_run_ids"] = [doc["runs"][0]["run_id"]]
    doc["runs"][0]["parent_run_ids"] = [second["run_id"]]
    doc["runs"].append(second)
    rejected(doc, "RUN_CYCLE")


@pytest.mark.parametrize("path", ["../outside.csv", "data/../../outside.csv", "/absolute.csv", "C:/absolute.csv", "data\\native-events.csv", "~user/file.csv", "."])
def test_relative_asset_locations_cannot_escape_the_project(path):
    doc = example()
    doc["assets"][0]["location"]["path"] = path
    rejected(doc, "INVALID_RELATIVE_PATH")


def test_declared_file_size_and_hash_are_actually_checked():
    doc = example()
    doc["assets"][0]["bytes"] += 1
    rejected(doc, "SIZE_MISMATCH", verify_hashes=True)
    doc = example()
    doc["assets"][0]["sha256"] = "0" * 64
    rejected(doc, "HASH_MISMATCH", verify_hashes=True)


def test_metadata_only_does_not_claim_a_file_exists_or_has_matching_hash():
    doc = example()
    doc["assets"][0]["location"]["path"] = "data/not-really-present.csv"
    report = accepted(doc)
    assert report["scope"] == "metadata_only"
    assert report["file_verification"]["hashes_checked"] == 0
    rejected(doc, "FILE_MISSING", check_files=True)


def test_remote_assets_are_never_fetched(monkeypatch):
    def prohibited(*args, **kwargs):
        raise AssertionError("Contract validation must not fetch remote content")
    import urllib.request
    monkeypatch.setattr(urllib.request, "urlopen", prohibited)
    asset = example()["assets"][0]
    asset["location"] = {"kind": "remote", "path": "s3://synthetic-fixture/native-events.csv"}
    asset["availability"] = "remote"
    report = accepted(asset, verify_hashes=True)
    assert "REMOTE_NOT_FETCHED" in issue_codes(report)
    assert report["file_verification"]["remote_assets_not_fetched"] == 1
    assert report["file_verification"]["hashes_checked"] == 0


@pytest.mark.parametrize("uri", ["https://user:secret@example.org/file", "http://example.org/file", "s3:///missing-bucket"])
def test_remote_locator_rejects_embedded_credentials_or_unsupported_forms(uri):
    asset = example()["assets"][0]
    asset["location"] = {"kind": "remote", "path": uri}
    asset["availability"] = "remote"
    rejected(asset, "INVALID_REMOTE_PATH")


@pytest.mark.parametrize("kind", ["asset", "run", "result"])
def test_standalone_records_explicitly_warn_that_project_links_are_unchecked(kind):
    field = {"asset": "assets", "run": "runs", "result": "results"}[kind]
    report = accepted(deepcopy(example()[field][0]))
    assert "PROJECT_CONTEXT_NOT_CHECKED" in issue_codes(report)


@pytest.mark.parametrize("schema_version", ["2.0", "9.0.0"])
def test_legacy_or_unknown_contract_requires_an_explicit_adapter(schema_version):
    doc = example()
    doc["schema_version"] = schema_version
    rejected(doc, "UNSUPPORTED_CONTRACT")


@pytest.mark.parametrize("kind", [[], {"bad": "kind"}])
def test_invalid_kind_type_is_a_diagnostic_not_an_uncaught_exception(kind):
    rejected({"kind": kind, "schema_version": "1.0.0"}, "UNSUPPORTED_CONTRACT")


@pytest.mark.parametrize("collection,field", [
    ("runs", "status"), ("results", "status"), ("assets", "availability"),
])
def test_malformed_record_state_does_not_crash_state_summary(collection, field):
    doc = example()
    doc[collection][0][field] = ["not-a-scalar"]
    rejected(doc, "SCHEMA")


@pytest.mark.parametrize("collection", ["runs", "results", "assets"])
def test_nonarray_collection_is_a_schema_diagnostic_not_summary_crash(collection):
    doc = example()
    doc[collection] = None
    rejected(doc, "SCHEMA")


def test_recorded_state_summary_preserves_missing_instead_of_inventing_success():
    doc = example("missing-result.project.json")
    report = accepted(doc)
    assert report["recorded_states"] == {
        "runs": {"blocked": 1},
        "results": {"not_computed": 1},
        "assets": {"missing": 1},
    }


def test_nonfinite_in_memory_values_are_rejected():
    doc = example()
    doc["extensions"]["fixture"]["not_a_number"] = float("nan")
    rejected(doc, "NONFINITE_NUMBER")


@pytest.mark.parametrize("text", [
    '{"kind":"project","kind":"asset","schema_version":"1.0.0"}',
    '{"kind":"project","schema_version":"1.0.0","bad":NaN}',
    "{ malformed json",
])
def test_invalid_json_duplicate_keys_and_nonjson_numbers_are_rejected(tmp_path, text):
    path = tmp_path / "invalid.json"
    path.write_text(text, encoding="utf-8")
    report = validate_file(path)
    assert report["valid"] is False
    assert "DOCUMENT_LOAD" in issue_codes(report)


def test_cli_json_exit_codes_and_no_project_side_effects(tmp_path):
    project = tmp_path / "portable-project"
    shutil.copytree(EXAMPLES, project)
    before = tree_identity(project)
    valid = cli("validate", project / EXAMPLE_NAMES[0], "--verify-hashes")
    assert valid.returncode == 0, valid.stderr
    report = json.loads(valid.stdout)
    assert report["valid"] is True
    assert report["execution_authorized"] is False
    assert report["scientific_tools_executed"] == 0
    invalid_doc = example()
    invalid_doc["schema_version"] = "2.0"
    invalid = tmp_path / "legacy-project.json"
    invalid.write_text(json.dumps(invalid_doc), encoding="utf-8")
    failed = cli("validate", invalid)
    assert failed.returncode == 2
    assert "UNSUPPORTED_CONTRACT" in issue_codes(json.loads(failed.stdout))
    assert tree_identity(project) == before
    assert not list(project.rglob("*.duckdb"))
    assert not (project / "artifacts").exists()
    assert not list(project.rglob("__pycache__"))


def test_fresh_process_avoids_gui_legacy_services_catalog_and_scientific_imports():
    code = r'''
import builtins, json, sys
blocked = ("PySide6", "duckdb", "pandas", "numpy", "scipy", "statsmodels",
           "katana_next.services", "katana_next.infrastructure", "katana_next.api")
original = builtins.__import__
def guard(name, *args, **kwargs):
    if any(name == p or name.startswith(p + ".") for p in blocked):
        raise AssertionError("Forbidden import for metadata validation: " + name)
    return original(name, *args, **kwargs)
builtins.__import__ = guard
from katana_next.contracts import validate_file
report = validate_file(sys.argv[1], verify_hashes=True)
report["forbidden_modules_loaded"] = [m for m in sys.modules if any(m == p or m.startswith(p + ".") for p in blocked)]
print(json.dumps(report))
'''
    env = dict(os.environ, PYTHONPATH=str(SOURCE), PYTHONDONTWRITEBYTECODE="1")
    result = subprocess.run(
        [sys.executable, "-B", "-c", code, str(EXAMPLES / EXAMPLE_NAMES[0])],
        env=env, cwd=SOURCE, capture_output=True, text=True,
        encoding="utf-8", timeout=30, check=False,
    )
    assert result.returncode == 0, result.stderr
    report = json.loads(result.stdout)
    assert report["valid"] is True
    assert report["forbidden_modules_loaded"] == []
