"""Real frozen bulk workflows through private flow task catalogues."""
import json
from pathlib import Path

import pandas as pd
import pytest

from katana_next.flow.collect import collect_results, materialize_report_descriptor, validate_collection
from katana_next.flow.model import file_sha256, seal_plan, write_json
from katana_next.flow.runtime import execute_task, read_json, validate_receipt
from katana_next.full_fixture import create_full_fixture


def _plan(tmp_path, tasks, project=None):
    project = project or create_full_fixture(tmp_path / "source")
    payload = seal_plan({"project_manifest": str(project),
                         "project_manifest_sha256": file_sha256(project),
                         "tasks": tasks})
    return write_json(tmp_path / "plan.json", payload), project


def _task(identity, workflow, dependencies=(), parameters=None, **extra):
    return {"id": identity, "workflow_id": workflow, "parameters": parameters or {},
            "dependencies": list(dependencies), "scope": {}, **extra}


def test_real_bulk_chain_private_catalogues_and_portable_collection(tmp_path):
    tasks = [_task("quality", "quality.stable_as"),
             _task("as", "differential.as", ["quality"]),
             _task("overlap", "differential.cross_comparison", ["as"])]
    plan, project = _plan(tmp_path, tasks)
    original = project.read_bytes()
    source = read_json(project)
    paths = [Path(a["path"]) for a in source["assets"]]
    hashes = {str(p): file_sha256(p) for p in paths}
    quality = execute_task(plan, "quality", tmp_path / "quality")
    differential = execute_task(plan, "as", tmp_path / "as", [quality])
    # Nextflow passes earlier levels too; only declared parents are consumed.
    overlap = execute_task(plan, "overlap", tmp_path / "overlap", [quality, differential])
    quality_data = validate_receipt(quality)
    as_data = validate_receipt(differential)
    assert quality_data["job"]["status"] == "complete"
    assert as_data["job"]["status"] == "complete"
    assert set(as_data["job"]["artifact_ids"]) == {a["artifact_id"] for a in as_data["generated_artifacts"]}
    actual = next(a for a in as_data["generated_artifacts"] if a["kind"] == "table" and a["name"] == "Differential AS events")
    frame = pd.read_parquet(actual["path"])
    assert {"FDR", "fdr", "dpsi", "significant", "stable_significant"} <= set(frame)
    assert frame["FDR"].equals(frame["fdr"])
    assert len(frame) > 0 and frame["comparison_id"].nunique() == 3
    assert len({Path(a).parent / "results/.katana/catalog.duckdb" for a in (quality, differential, overlap)}) == 3
    assert not (Path(source["output_root"]) / ".katana").exists()
    collection = collect_results(plan, [quality, differential, overlap], tmp_path / "collection")
    receipt = validate_collection(collection)
    assert receipt["status"] == "completed" and len(receipt["tasks"]) == 3
    from katana_next.reports.sources import load_source
    report_source = load_source(materialize_report_descriptor(collection))
    assert len(report_source["artifacts"]) >= 3
    assert all(str(tmp_path / "collection") in a["path"] for a in report_source["artifacts"])
    assert all(Path(a["path"]).is_file() for a in report_source["artifacts"])
    # Copying the collection to another mount/directory does not require any
    # source-catalog edits or rewrites to its recorded scientific products.
    import shutil
    moved = tmp_path / "relocated"
    shutil.copytree(collection.parent, moved)
    relocated = load_source(materialize_report_descriptor(moved))
    assert all(str(moved) in a["path"] for a in relocated["artifacts"])
    assert all(Path(a["path"]).is_file() for a in relocated["artifacts"])
    assert project.read_bytes() == original
    assert hashes == {str(p): file_sha256(p) for p in paths}


def test_complete_receipt_reuse_keeps_original_job_and_products(tmp_path):
    plan, _ = _plan(tmp_path, [_task("quality", "quality.stable_as")])
    receipt = execute_task(plan, "quality", tmp_path / "task")
    before, modified = receipt.read_bytes(), receipt.stat().st_mtime_ns
    assert execute_task(plan, "quality", tmp_path / "task") == receipt
    assert receipt.read_bytes() == before and receipt.stat().st_mtime_ns == modified


def test_completed_product_corruption_does_not_trigger_a_new_fit(tmp_path):
    plan, _ = _plan(tmp_path, [_task("quality", "quality.stable_as")])
    receipt = execute_task(plan, "quality", tmp_path / "task")
    before = receipt.read_bytes()
    data = read_json(receipt)
    table = next(a for a in data["generated_artifacts"] if a["kind"] == "table")
    Path(table["path"]).write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="Verified content changed"):
        execute_task(plan, "quality", tmp_path / "task")
    assert receipt.read_bytes() == before


def test_failed_attempt_retains_claim_and_cannot_be_overwritten(tmp_path):
    plan, _ = _plan(tmp_path, [_task("bad", "not.a.workflow")])
    with pytest.raises(KeyError, match="Unknown workflow"):
        execute_task(plan, "bad", tmp_path / "task")
    assert (tmp_path / "task/execution.claim.json").is_file()
    failure = read_json(tmp_path / "task/failure.receipt.json")
    assert failure["status"] == "failed" and not (tmp_path / "task/task.receipt.json").exists()
    with pytest.raises(FileExistsError, match="unknown or failed attempt"):
        execute_task(plan, "bad", tmp_path / "task")


def test_source_revision_and_plan_tampering_fail_before_execution(tmp_path):
    plan, project = _plan(tmp_path, [_task("quality", "quality.stable_as")])
    project.write_text(project.read_text() + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="changed after planning"):
        execute_task(plan, "quality", tmp_path / "task")
    assert not (tmp_path / "task").exists()
    raw = read_json(plan)
    raw["tasks"][0]["parameters"]["fdr"] = 0.9
    write_json(plan, raw)
    with pytest.raises(ValueError, match="content changed"):
        execute_task(plan, "quality", tmp_path / "other")


def test_missing_planned_producer_cannot_be_collected(tmp_path):
    plan, _ = _plan(tmp_path, [_task("quality", "quality.stable_as"),
                               _task("as", "differential.as", ["quality"])])
    receipt = execute_task(plan, "quality", tmp_path / "quality")
    with pytest.raises(ValueError, match="no completed products"):
        collect_results(plan, [receipt], tmp_path / "collection")
    assert not (tmp_path / "collection").exists()


def test_real_disjoint_contrast_join_keeps_native_statistics(tmp_path):
    project = create_full_fixture(tmp_path / "source")
    ids = [c["contrast_id"] for c in read_json(project)["contrasts"]][:2]
    tasks = [_task("as-a", "differential.as", parameters={"contrast_ids": [ids[0]]}),
             _task("as-b", "differential.as", parameters={"contrast_ids": [ids[1]]}),
             _task("overlap", "differential.cross_comparison", ["as-a", "as-b"], {"contrast_ids": ids})]
    plan, _ = _plan(tmp_path, tasks, project)
    a = execute_task(plan, "as-a", tmp_path / "a")
    b = execute_task(plan, "as-b", tmp_path / "b")
    result = execute_task(plan, "overlap", tmp_path / "overlap", [a, b])
    receipt = validate_receipt(result)
    joined = [artifact for artifact in receipt["artifacts"] if artifact["metadata"].get("flow_lossless_parent_join")]
    assert joined
    for artifact in joined:
        parents = [parent for parent in receipt["artifacts"] if parent["artifact_id"] in artifact["parent_artifact_ids"]]
        expected = pd.concat([pd.read_parquet(parent["path"]) for parent in parents], ignore_index=True)
        actual = pd.read_parquet(artifact["path"])
        pd.testing.assert_frame_equal(actual, expected)
        assert actual["row_id"].is_unique


def test_empty_rs_foreground_is_audited_not_a_fabricated_motif_run(tmp_path):
    tasks = [_task("as", "differential.as", parameters={"fdr": 0}),
             _task("motif", "prioritize.motiftrace_consequences", ["as"],
                   selection={"analysis_id": "differential.as", "name": "Differential AS events",
                              "event_set": "stable_significant"})]
    plan, _ = _plan(tmp_path, tasks)
    first = execute_task(plan, "as", tmp_path / "as")
    result = execute_task(plan, "motif", tmp_path / "motif", [first])
    receipt = validate_receipt(result)
    assert receipt["status"] == "not_applicable"
    assert receipt["selection"]["source_rows"] > 0 and receipt["selection"]["selected_rows"] == 0
    assert "job" not in receipt
    assert receipt["execution_evidence"]["reason"] == "No events satisfy the declared foreground selection"
