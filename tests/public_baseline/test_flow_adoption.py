"""Historical import is verified transport, never a fabricated scientific run."""
import copy
import json
from pathlib import Path
import shutil

import pytest

from katana_next.flow.adoption import (seal_existing_results, validate_existing_source,
    adopt_existing_results, validate_adoption, materialize_adopted_descriptor, create_identity_projection)
from katana_next.flow.runtime import file_record, read_json, write_json


def fixture_source(tmp_path):
    root = tmp_path / "original results"
    root.mkdir()
    project = write_json(root / "project.json", {"schema_version": "2.0", "project_id": "old-study",
        "name": "Existing result fixture", "input_root": str(root), "output_root": str(root / "results"),
        "annotation": {"genome_build": "test"}, "samples": [
            {"sample_id": "a", "biological_sample_id": "mouse-a", "condition": "A"},
            {"sample_id": "b", "biological_sample_id": "mouse-b", "condition": "B"}],
        "contrasts": [{"contrast_id": "B_vs_A", "label": "B versus A",
            "numerator": {"label": "B", "sample_ids": ["b"]},
            "denominator": {"label": "A", "sample_ids": ["a"]}}]})
    family = root / "old-native-family"
    family.mkdir()
    result = family / "full.tsv"
    result.write_text("event_id\tgene_id\tFDR\tdelta_psi\nE1\tG1\t0.01\t-0.2\nE2\tG2\t0.7\t0.03\n", encoding="utf-8")
    member = file_record(result)
    write_json(family / "COMPLETE.json", {"complete": True, "members": {
        result.name: {"bytes": member["size_bytes"], "sha256": member["sha256"]}}})
    completion = write_json(root / "original-native.receipt.json", {"native_job_id": "old-job-12",
        "status": "COMPLETED", "exit_code": "0:0", "parameters": {"event_set": "RS-DAS", "seed": 7}})
    record = {"artifact_id": "old-artifact", "project_id": "old-study", "analysis_id": "differential.as",
        "name": "All events", "kind": "table", "entity_type": "splicing_event", "status": "complete",
        "fingerprint": "original-scientific-fingerprint", "run_fingerprint": "original-run",
        "path": str(result), "contrast_ids": ["B_vs_A"], "parameters": {"event_set": "RS-DAS", "fdr": 0.05},
        "metadata": {"workflow_version": "old-version", "source_sample_order": ["a", "b"]}}
    descriptor = write_json(root / "source.json", {"project_manifest": str(project), "artifacts": [record],
        "report_contrast_ids": ["B_vs_A"], "extra_views": [{"project_id": "old-study",
            "source_project_manifest_sha256": file_record(project)["sha256"], "status": "completed",
            "panels": [], "tables": [{"chapter": 3, "path": str(result), "title": "Full original table"}],
            "unavailable": [{"chapter": 6, "component": "module_bootstrap", "reason": "not computed"}]}]})
    return root, project, descriptor, family, completion, record


def seal(tmp_path):
    root, project, descriptor, family, completion, record = fixture_source(tmp_path)
    source = seal_existing_results(descriptor, tmp_path / "sealed/source.json", allowed_roots=[root],
        completion_evidence=[completion], family_roots=[family])
    return root, project, descriptor, family, completion, record, source


def test_adoption_preserves_original_products_and_reports_import_truthfully(tmp_path):
    root, project, descriptor, family, completion, record, source = seal(tmp_path)
    before = {str(path): (path.read_bytes(), path.stat().st_mtime_ns) for path in root.rglob("*") if path.is_file()}
    receipt_path = adopt_existing_results(source, tmp_path / "adopted")
    receipt = validate_adoption(receipt_path)
    assert receipt["status"] == "adopted"
    assert receipt["scientific_computation_performed"] is False
    assert receipt["native_jobs_submitted"] == 0
    view = materialize_adopted_descriptor(receipt_path)
    result = view["artifacts"][0]
    for key in ("artifact_id", "fingerprint", "run_fingerprint", "contrast_ids", "parameters", "metadata"):
        # Reader completion bookkeeping is separate from original native evidence.
        if key == "metadata":
            assert all(result[key][name] == value for name, value in record[key].items())
        else:
            assert result[key] == record[key]
    assert Path(result["path"]).read_bytes() == Path(record["path"]).read_bytes()
    assert Path(view["project_manifest"]).read_bytes() == project.read_bytes()
    assert view["extra_views"][0]["unavailable"][0]["component"] == "module_bootstrap"
    assert before == {str(path): (path.read_bytes(), path.stat().st_mtime_ns) for path in root.rglob("*") if path.is_file()}
    receipt_bytes = receipt_path.read_bytes()
    assert adopt_existing_results(source, receipt_path.parent) == receipt_path
    assert receipt_path.read_bytes() == receipt_bytes


def test_portable_snapshot_reloads_without_original_sources_and_preserves_native_json(tmp_path):
    root, project, descriptor, family, completion, record, source = seal(tmp_path)
    original_project = project.read_bytes()
    original_marker = (family / "COMPLETE.json").read_bytes()
    adopted = adopt_existing_results(source, tmp_path / "adopted")
    relocated = tmp_path / "another-machine" / "adopted"
    shutil.copytree(adopted.parent, relocated)
    shutil.move(str(root), str(tmp_path / "disconnected-original"))
    validate_adoption(relocated)
    result = materialize_adopted_descriptor(relocated)
    assert Path(result["artifacts"][0]["path"]).is_file()
    assert Path(result["project_manifest"]).read_bytes() == original_project
    assert next(relocated.rglob("COMPLETE.json")).read_bytes() == original_marker
    # The raw manifest stays unchanged; portability here covers products and
    # views, not remapping external raw-read/reference computation resources.
    assert result["historical_adoption"]["scientific_computation_performed"] is False


@pytest.mark.parametrize("change", ["table", "proof", "project", "descriptor", "family"])
def test_changed_historical_input_is_rejected_without_any_new_success(tmp_path, change):
    root, project, descriptor, family, completion, record, source = seal(tmp_path)
    targets = {"table": family / "full.tsv", "proof": completion, "project": project,
               "descriptor": descriptor, "family": family / "COMPLETE.json"}
    targets[change].write_text("modified", encoding="utf-8")
    with pytest.raises((ValueError, FileNotFoundError)):
        validate_existing_source(source)
    with pytest.raises((ValueError, FileNotFoundError)):
        adopt_existing_results(source, tmp_path / "never-created")
    assert not (tmp_path / "never-created").exists()


@pytest.mark.parametrize("change", ["other_project", "unknown_contrast", "running", "outside"])
def test_unbound_incomplete_and_out_of_bounds_imports_are_rejected(tmp_path, change):
    root, project, descriptor, family, completion, record = fixture_source(tmp_path)
    value = read_json(descriptor)
    if change == "other_project":
        value["artifacts"][0]["project_id"] = "different-study"
    elif change == "unknown_contrast":
        value["artifacts"][0]["contrast_ids"] = ["fake_contrast"]
    elif change == "running":
        value["artifacts"][0]["status"] = "running"
    else:
        alien = tmp_path / "outside.tsv"
        alien.write_text("outside", encoding="utf-8")
        value["artifacts"][0]["path"] = str(alien)
    write_json(descriptor, value)
    with pytest.raises(ValueError):
        seal_existing_results(descriptor, tmp_path / "must-not-exist.json", allowed_roots=[root])
    assert not (tmp_path / "must-not-exist.json").exists()


def test_copied_products_and_receipt_are_authenticated_and_unknown_attempts_preserved(tmp_path):
    *_, source = seal(tmp_path)
    unknown = tmp_path / "unknown"
    unknown.mkdir()
    (unknown / "do-not-delete").write_text("previous attempt", encoding="utf-8")
    with pytest.raises(FileExistsError):
        adopt_existing_results(source, unknown)
    assert (unknown / "do-not-delete").read_text() == "previous attempt"
    adopted = adopt_existing_results(source, tmp_path / "adopted")
    next(adopted.parent.rglob("full.tsv")).write_text("tampered", encoding="utf-8")
    with pytest.raises(ValueError):
        validate_adoption(adopted)


def test_family_marker_cannot_smuggle_external_members(tmp_path):
    root, project, descriptor, family, completion, record = fixture_source(tmp_path)
    write_json(family / "COMPLETE.json", {"complete": True, "members": {"../../outside.tsv": "abc"}})
    with pytest.raises(ValueError, match="escapes"):
        seal_existing_results(descriptor, tmp_path / "bad.json", allowed_roots=[root], family_roots=[family])


def test_transport_task_uses_explicit_overlay_identity_and_retains_old_science(tmp_path):
    from katana_next.flow.runtime import execute_task, validate_receipt
    from katana_next.flow.model import seal_plan
    root, project, descriptor, family, completion, record, source = seal(tmp_path)
    adopted = adopt_existing_results(source, tmp_path / "adopted")
    overlay = read_json(project)
    overlay["output_root"] = str(tmp_path / "new-compute-output")
    compute = write_json(tmp_path / "compute-project.json", overlay)
    projection = create_identity_projection(adopted, compute, tmp_path / "projection.json", contrast_ids=["B_vs_A"])
    task = {"id": "adopt-inputs", "workflow_id": "flow.historical.adopt", "adapter": "historical.adopt",
        "dependencies": [], "scope": {"contrast_ids": ["B_vs_A"]}, "parameters": {},
        "inputs": {"adoption": str(adopted), "identity_projection": str(projection)}}
    plan = seal_plan({"project_manifest": str(compute), "project_manifest_sha256": file_record(compute)["sha256"],
                      "tasks": [task]})
    plan_path = write_json(tmp_path / "plan.json", plan)
    receipt = execute_task(plan_path, task["id"], tmp_path / "import-task")
    actual = validate_receipt(receipt, plan=plan, task=plan["tasks"][0])
    assert actual["status"] == "completed"
    assert actual["workflow_id"] == "flow.historical.adopt"
    assert actual["execution_evidence"]["scientific_computation_performed"] is False
    assert actual["execution_evidence"]["native_jobs_submitted"] == 0
    assert "job" not in actual
    artifact = actual["generated_artifacts"][0]
    for key in ("analysis_id", "artifact_id", "fingerprint", "run_fingerprint", "parameters", "contrast_ids"):
        assert artifact[key] == record[key]
    assert Path(artifact["path"]).read_bytes() == Path(record["path"]).read_bytes()
    before = receipt.read_bytes()
    assert execute_task(plan_path, task["id"], tmp_path / "import-task") == receipt
    assert receipt.read_bytes() == before


@pytest.mark.parametrize("change", ["sample_order", "biological_unit", "arm", "build", "condition"])
def test_identity_projection_rejects_scientific_design_changes(tmp_path, change):
    root, project, descriptor, family, completion, record, source = seal(tmp_path)
    adopted = adopt_existing_results(source, tmp_path / "adopted")
    overlay = read_json(project)
    if change == "sample_order":
        overlay["samples"].reverse()
    elif change == "biological_unit":
        overlay["samples"][0]["biological_sample_id"] = "other-mouse"
    elif change == "arm":
        contrast = overlay["contrasts"][0]
        contrast["numerator"], contrast["denominator"] = contrast["denominator"], contrast["numerator"]
    elif change == "build":
        overlay["annotation"]["genome_build"] = "different-build"
    else:
        overlay["samples"][0]["condition"] = "new-condition"
    compute = write_json(tmp_path / "changed-project.json", overlay)
    with pytest.raises(ValueError, match="changes original"):
        create_identity_projection(adopted, compute, tmp_path / "projection.json", contrast_ids=["B_vs_A"])


def test_actual_old_as_is_imported_for_new_overlap_without_repeating_quality_or_as(tmp_path, monkeypatch):
    from katana_next.flow.runtime import execute_task, validate_receipt
    from katana_next.flow.model import seal_plan
    from katana_next.full_fixture import create_full_fixture
    native_root = tmp_path / "native"
    original_project = create_full_fixture(native_root / "source")
    old_tasks = [{"id": "quality", "workflow_id": "quality.stable_as", "parameters": {}, "scope": {}, "dependencies": []},
                 {"id": "as", "workflow_id": "differential.as", "parameters": {}, "scope": {}, "dependencies": ["quality"]}]
    old_plan = write_json(native_root / "plan.json", seal_plan({"project_manifest": str(original_project),
        "project_manifest_sha256": file_record(original_project)["sha256"], "tasks": old_tasks}))
    quality = execute_task(old_plan, "quality", native_root / "quality")
    old_as = execute_task(old_plan, "as", native_root / "as", [quality])
    old_receipt_bytes = old_as.read_bytes()
    old_receipt = validate_receipt(old_as)
    descriptor = write_json(native_root / "descriptor.json", {"project_manifest": str(original_project),
        "artifacts": old_receipt["artifacts"]})
    source = seal_existing_results(descriptor, tmp_path / "seal/source.json", allowed_roots=[native_root],
        completion_evidence=[quality, old_as])
    adopted = adopt_existing_results(source, tmp_path / "adopted")
    overlay = read_json(original_project)
    overlay["output_root"] = str(tmp_path / "unmodified-future-results")
    compute = write_json(tmp_path / "compute.json", overlay)
    selected = [row["contrast_id"] for row in overlay["contrasts"]]
    projection = create_identity_projection(adopted, compute, tmp_path / "projection.json", contrast_ids=selected)
    tasks = [{"id": "historical", "workflow_id": "flow.historical.adopt", "adapter": "historical.adopt",
              "scope": {"contrast_ids": selected}, "parameters": {}, "dependencies": [],
              "inputs": {"adoption": str(adopted), "identity_projection": str(projection)}},
             {"id": "new-overlap", "workflow_id": "differential.cross_comparison", "scope": {},
              "parameters": {}, "dependencies": ["historical"]}]
    plan = write_json(tmp_path / "new-plan.json", seal_plan({"project_manifest": str(compute),
        "project_manifest_sha256": file_record(compute)["sha256"], "tasks": tasks}))
    # Fail loudly if historical analysis kernels are called by the new path.
    from katana_next.client import KatanaProject
    original_run = KatanaProject.run
    def guarded_run(self, workflow_id, *args, **kwargs):
        if workflow_id in {"quality.stable_as", "differential.as"}:
            raise AssertionError("Completed historical science was rerun")
        return original_run(self, workflow_id, *args, **kwargs)
    monkeypatch.setattr(KatanaProject, "run", guarded_run)
    historical = execute_task(plan, "historical", tmp_path / "historical-task")
    overlap = execute_task(plan, "new-overlap", tmp_path / "overlap-task", [historical])
    new = validate_receipt(overlap)
    assert new["job"]["status"] == "complete"
    assert new["job"]["analysis_id"] == "differential.cross_comparison"
    assert any(row["analysis_id"] == "differential.as" for row in new["artifacts"])
    assert old_as.read_bytes() == old_receipt_bytes


def test_analysis_only_transport_keeps_legacy_policy_and_explicitly_excludes_exports(tmp_path):
    from katana_next.flow.adoption import import_adopted_artifacts
    root, project, descriptor, family, completion, record = fixture_source(tmp_path)
    value = read_json(descriptor)
    legacy = value["artifacts"][0]
    legacy["metadata"].update(workflow_version="2.0.0", effective_evidence_parameters={
        "min_reads": 10, "min_supported_replicates": 2, "min_direction_consistency": 0.67})
    export = root / "human-readable.tsv"
    export.write_text("Unchanged old human-readable duplicate", encoding="utf-8")
    legacy["metadata"]["exports"] = {"tsv": str(export)}
    write_json(descriptor, value)
    sealed = seal_existing_results(descriptor, tmp_path / "source.json", allowed_roots=[root],
                                   purpose="analysis-inputs", include_exports=False)
    snapshot = read_json(sealed)
    assert snapshot["purpose"] == "analysis-inputs"
    assert str(export) not in {row["path"] for row in snapshot["files"]}
    adopted = adopt_existing_results(sealed, tmp_path / "adopted")
    projection = create_identity_projection(adopted, project, tmp_path / "projection.json", contrast_ids=["B_vs_A"])
    result = import_adopted_artifacts(adopted, projection, project, tmp_path / "transport")
    actual = result["artifacts"][0]
    assert actual["metadata"]["workflow_version"] == "2.0.0"
    assert actual["metadata"]["effective_evidence_parameters"] == legacy["metadata"]["effective_evidence_parameters"]
    assert "read_support_policy" not in actual["metadata"]
    assert actual["metadata"]["exports"] == legacy["metadata"]["exports"]
    assert not any(row["path"].endswith(export.name) for row in result["outputs"])
    with pytest.raises(ValueError, match="analysis-inputs"):
        seal_existing_results(descriptor, tmp_path / "bad.json", allowed_roots=[root], include_exports=False)


def test_generic_cli_seals_adopts_and_binds_without_hardcoded_study(tmp_path, capsys):
    from katana_next.flow.cli import main
    root, project, descriptor, family, completion, record = fixture_source(tmp_path)
    sealed = tmp_path / "source.json"
    assert main(["seal-existing", "--descriptor", str(descriptor), "--out", str(sealed),
                 "--allowed-root", str(root), "--completion-evidence", str(completion),
                 "--family-root", str(family)]) == 0
    assert json.loads(capsys.readouterr().out)["receipt"] == str(sealed)
    adopted = tmp_path / "adopted"
    assert main(["adopt-existing", "--source", str(sealed), "--outdir", str(adopted)]) == 0
    assert Path(json.loads(capsys.readouterr().out)["receipt"]).is_file()
    projection = tmp_path / "projection.json"
    assert main(["bind-existing", "--adoption", str(adopted), "--project", str(project),
                 "--contrast-id", "B_vs_A", "--out", str(projection)]) == 0
    assert Path(json.loads(capsys.readouterr().out)["receipt"]).is_file()
    assert main(["seal-existing", "--descriptor", str(descriptor), "--out", str(sealed),
                 "--allowed-root", str(root)]) == 1
    assert json.loads(capsys.readouterr().err)["status"] == "failed"


def test_preflight_authenticates_historical_transport_instead_of_unknown_adapter(tmp_path):
    from katana_next.flow.model import seal_plan
    from katana_next.flow.planner import preflight_plan
    root, project, descriptor, family, completion, record, source = seal(tmp_path)
    adopted = adopt_existing_results(source, tmp_path / "adopted")
    projection = create_identity_projection(adopted, project, tmp_path / "projection.json", contrast_ids=["B_vs_A"])
    task = {"id": "history", "workflow_id": "flow.historical.adopt", "adapter": "historical.adopt",
            "scope": {}, "parameters": {}, "dependencies": [],
            "inputs": {"adoption": str(adopted), "identity_projection": str(projection)}}
    plan = seal_plan({"project_manifest": str(project), "project_manifest_sha256": file_record(project)["sha256"],
                      "configuration": {"design_factors": ["condition"]}, "tasks": [task]})
    # Missing new-project resources remain honest unrelated readiness errors;
    # the explicitly proved historical adapter itself is recognized and checked.
    checked = preflight_plan(plan)
    assert not [row for row in checked["errors"] if row["task_id"] == "history"]
    projection.write_text("{}", encoding="utf-8")
    checked = preflight_plan(plan)
    assert any(row["task_id"] == "history" for row in checked["errors"])
