"""Actual completed-result import and collection, without new native fitting."""
import importlib.util
from pathlib import Path

import pytest

from katana_next.flow.collect import collect_results, materialize_report_descriptor, validate_collection
from katana_next.flow.model import file_sha256, seal_plan, write_json
from katana_next.flow.runtime import execute_task, read_json, validate_receipt


def _source(tmp_path):
    # Reuse the frozen reader's complete structural fixture. Its synthetic
    # calibration numbers are importer test data, not scientific calibration.
    spec = importlib.util.spec_from_file_location("parallel_reader_fixture", Path(__file__).parents[1] / "v2/test_motiftrace_parallel.py")
    fixture_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fixture_module)
    bundle_root = tmp_path / "原始-parallel"
    bundle_root.mkdir()
    _, _, _, seal = fixture_module.bundle.__wrapped__(bundle_root)
    asset = seal()
    from katana_next.domain.models import AnnotationReference, Contrast, ContrastArm, InputAsset, ProjectManifest, Sample
    from katana_next.infrastructure.manifest_store import ManifestStore
    source_root = tmp_path / "source"
    source_root.mkdir()
    manifest = ProjectManifest(project_id="importer-fixture", name="Generated completed-import fixture",
        input_root=source_root, output_root=source_root / "existing-results",
        annotation=AnnotationReference(genome_build="hg19"), metadata={"species": "human", "synthetic_fixture": True},
        samples=tuple(Sample(sample_id=s, biological_sample_id=s, condition="KD" if s.startswith("k") else "CTRL")
                      for s in ("k1", "k2", "c1", "c2")),
        contrasts=(Contrast(contrast_id="KD", label="KD versus CTRL",
            numerator=ContrastArm(label="KD", sample_ids=("k1", "k2")),
            denominator=ContrastArm(label="CTRL", sample_ids=("c1", "c2"))),),
        assets=(InputAsset(asset_id="existing-parallel", kind="motiftrace_parallel_results", path=asset.path, metadata=asset.metadata),))
    project = ManifestStore().save(manifest, source_root / "project.json")
    task = {"id": "parallel", "workflow_id": "prioritize.motiftrace_parallel", "workflow_version": "1.0.0",
            "dependencies": [], "scope": {"contrast_ids": ["KD"]},
            "parameters": {"contrast_ids": ["KD"], "parallel_bundle_asset_id": "existing-parallel"}}
    plan = seal_plan({"project_manifest": str(project), "project_manifest_sha256": file_sha256(project),
                      "contrast_ids": ["KD"], "tasks": [task]})
    plan_path = write_json(tmp_path / "plan.json", plan)
    return project, plan_path, bundle_root


def test_actual_completed_parallel_import_is_task_owned_and_collects_without_source_edits(tmp_path):
    project, plan, bundle_root = _source(tmp_path)
    original = {str(path): (path.read_bytes(), path.stat().st_mtime_ns) for path in [project, *bundle_root.iterdir()]}
    result = execute_task(plan, "parallel", tmp_path / "task")
    receipt = validate_receipt(result)
    assert receipt["job"]["status"] == "complete"
    assert receipt["execution_evidence"]["native_fitting_performed"] is False
    assert receipt["execution_evidence"]["new_inference"] is False
    assert all(Path(row["path"]).is_relative_to(tmp_path / "task") for row in receipt["generated_artifacts"])
    original_manifest = bundle_root / "manifest.json"
    assert any(row["path"] == str(original_manifest) for row in receipt["inputs"])
    copied_manifest = Path(receipt["execution_evidence"]["owned_completed_parallel_snapshot"]["private_manifest"])
    assert copied_manifest.read_bytes() == original_manifest.read_bytes()
    assert any(row["path"] == str(copied_manifest) for row in receipt["outputs"])
    collection = collect_results(plan, [result], tmp_path / "collection")
    assert validate_collection(collection)["status"] == "completed"
    source = materialize_report_descriptor(collection)
    browser = next(row for row in source["artifacts"] if row["entity_type"] == "motiftrace_parallel_report")
    assert Path(browser["path"]).is_relative_to(tmp_path / "collection")
    from katana_next.reports import build_report
    report = build_report(source, tmp_path / "report")
    assert not any(row.get("metadata", {}).get("component") == "mapp_parallel" for row in report["unavailable"])
    assert (tmp_path / "report/REPORT.zip").stat().st_size > 0
    assert original == {str(path): (path.read_bytes(), path.stat().st_mtime_ns) for path in [project, *bundle_root.iterdir()]}
    before = result.read_bytes()
    assert execute_task(plan, "parallel", tmp_path / "task") == result
    assert result.read_bytes() == before


@pytest.mark.parametrize("change", ["tamper", "missing"])
def test_import_reuse_remains_bound_to_original_completed_source(tmp_path, change):
    project, plan, bundle_root = _source(tmp_path)
    result = execute_task(plan, "parallel", tmp_path / "task")
    receipt = result.read_bytes()
    member = bundle_root / "models.tsv"
    if change == "tamper":
        member.write_text("modified original source", encoding="utf-8")
    else:
        member.unlink()
    with pytest.raises((ValueError, FileNotFoundError)):
        execute_task(plan, "parallel", tmp_path / "task")
    assert result.read_bytes() == receipt
    with pytest.raises((ValueError, FileNotFoundError)):
        execute_task(plan, "parallel", tmp_path / "new-task")
    assert not (tmp_path / "new-task/task.receipt.json").exists()
