"""Execution planning must retain selected comparisons and fail closed on gaps."""
import copy
import json
from pathlib import Path

import pytest

from katana_next.flow.model import seal_plan, validate_plan, write_json
from katana_next.flow.planner import create_plan, preflight_plan


@pytest.fixture
def study(tmp_path):
    from katana_next.full_fixture import create_full_fixture
    return create_full_fixture(tmp_path / "study")


def test_only_requested_comparisons_get_independent_expensive_tasks(study):
    from katana_next import __version__
    plan = create_plan(study, {"contrast_ids": ["perturbation-vs-baseline"]})
    validate_plan(plan)
    assert plan["package_version"] == __version__
    motifs = [t for t in plan["tasks"] if t["workflow_id"] == "prioritize.motiftrace_consequences"]
    assert len(motifs) == 1
    assert motifs[0]["scope"]["contrast_ids"] == ["perturbation-vs-baseline"]
    assert motifs[0]["selection"]["event_set"] == "stable_significant"
    enrichment = [t for t in plan["tasks"] if t["workflow_id"] == "prioritize.enrichment"]
    assert {t["parameters"]["foreground_source"] for t in enrichment} == {"das", "deg", "dtu", "dexseq"}
    assert all(t["parameters"]["background_source"] == "match_foreground" for t in enrichment)
    assert all(t["parameters"]["contrast_ids"] == ["perturbation-vs-baseline"] for t in enrichment)
    assert all(t["level"] < 12 for t in plan["tasks"])
    assert not any(t["id"] == "overlap" for t in plan["tasks"])
    assert plan["applicability"]["cross_comparison"]["status"] == "not_applicable"
    multiple = create_plan(study, {})
    assert any(t["id"] == "overlap" for t in multiple["tasks"])
    assert multiple["applicability"]["cross_comparison"]["status"] == "applicable"


def test_complete_report_does_not_silently_disable_missing_computations(study):
    plan = create_plan(study, {"design_factors": ["background", "treatment"]})
    checked = preflight_plan(plan)
    assert checked["ready"] is False
    descriptions = json.dumps(checked["errors"])
    assert "rbp_binding_peaks" in descriptions
    assert "interaction" in descriptions
    assert "mapp" in descriptions
    assert "genebody" in descriptions
    assert next(t for t in plan["tasks"] if t["id"] == "reference-bed")["adapter"] == "reference.bed12"
    assert {t.get("adapter") for t in plan["tasks"]} >= {"network.fit", "network.bootstrap", "qc.genebody", "mapp.maei"}


def test_changed_plan_and_cycles_are_rejected(study):
    plan = create_plan(study, {})
    broken = copy.deepcopy(plan)
    broken["tasks"][0]["parameters"]["fdr"] = 0.99
    with pytest.raises(ValueError, match="changed"):
        validate_plan(broken)
    with pytest.raises(ValueError, match="cycle"):
        seal_plan({"tasks": [{"id": "a", "workflow_id": "project.audit", "dependencies": ["b"]},
                             {"id": "b", "workflow_id": "project.audit", "dependencies": ["a"]}]})


def test_ordinary_das_switch_does_not_weaken_required_rs_motif_inputs(study):
    plan = create_plan(study, {"event_set": "significant", "contrast_ids": ["perturbation-vs-baseline"]})
    motif = next(t for t in plan["tasks"] if t["workflow_id"] == "prioritize.motiftrace_consequences")
    assert motif["selection"]["event_set"] == "stable_significant"
    multiple = create_plan(study, {"event_set": "significant"})
    overlap = next(t for t in multiple["tasks"] if t["id"] == "overlap")
    assert overlap["parameters"]["event_set"] == "significant"
    network = next(t for t in plan["tasks"] if t.get("adapter") == "network.fit")
    assert network["scope"]["unit"] == "measurable_contrast_PSI"
    assert network["scope"]["contrast_ids"] == ["perturbation-vs-baseline"]


def test_real_audit_quality_and_as_chain_isolated_and_reusable(study, tmp_path):
    from katana_next.flow.task import execute_task
    from katana_next.flow.collect import collect_results, validate_collection
    from katana_next.flow.model import file_sha256
    original = file_sha256(study)
    original_manifest = json.loads(study.read_text("utf-8"))
    original_output = Path(original_manifest["output_root"])
    before_outputs = {str(x.relative_to(original_output)) for x in original_output.rglob("*")}
    plan = create_plan(study, {"contrast_ids": ["perturbation-vs-baseline"]})
    plan["tasks"] = [t for t in plan["tasks"] if t["id"] in {"audit", "quality", "as"}]
    plan["report_profile"] = "integration-test-subgraph"
    plan = seal_plan(plan)
    plan_path = write_json(tmp_path / "plan.json", plan)
    roots, receipts = [], []
    for identity in ("audit", "quality", "as"):
        root = tmp_path / ("task-" + identity)
        execute_task(plan_path, identity, root, dependency_receipts=receipts)
        roots.append(root)
        receipts.append(root / "task.receipt.json")
    saved = file_sha256(receipts[-1])
    execute_task(plan_path, "as", roots[-1], dependency_receipts=receipts[:-1])
    assert file_sha256(receipts[-1]) == saved
    collection = tmp_path / "collected"
    collect_results(plan_path, roots, collection)
    assert validate_collection(collection)["status"] == "completed"
    assert file_sha256(study) == original
    assert {str(x.relative_to(original_output)) for x in original_output.rglob("*")} == before_outputs
    assert not (original_output / ".katana" / "catalog.duckdb").exists()
    with pytest.raises(ValueError, match="no completed"):
        collect_results(plan_path, roots[:-1], tmp_path / "missing-results")


def test_changed_success_output_cannot_be_reused(study, tmp_path):
    from katana_next.flow.task import execute_task
    plan = create_plan(study, {})
    plan["tasks"] = [t for t in plan["tasks"] if t["id"] == "audit"]
    plan = seal_plan(plan)
    plan_path = write_json(tmp_path / "plan.json", plan)
    root = tmp_path / "audit"
    execute_task(plan_path, "audit", root)
    receipt = json.loads((root / "task.receipt.json").read_text("utf-8"))
    Path(receipt["outputs"][0]["path"]).write_text("corrupt", encoding="utf-8")
    with pytest.raises(ValueError, match="changed"):
        execute_task(plan_path, "audit", root)


def test_optional_parallel_absence_does_not_add_import_or_remove_native_fits(study):
    plan = create_plan(study, {"contrast_ids": ["perturbation-vs-baseline"]})
    assert not any(t["workflow_id"] == "prioritize.motiftrace_parallel" for t in plan["tasks"])
    assert len([t for t in plan["tasks"] if t.get("adapter") == "mapp.maei"]) == 14


def test_prepared_maei_inputs_need_one_exact_comparison(study):
    config = {"mapp": {"inputs": {"inclusion_table": "already-prepared.tsv"}}}
    with pytest.raises(ValueError, match="one exact selected comparison"):
        create_plan(study, config)
    plan = create_plan(study, {**config, "contrast_ids": ["perturbation-vs-baseline"]})
    assert next(t for t in plan["tasks"] if t["id"] == "mapp")["scope"]["contrast_ids"] == ["perturbation-vs-baseline"]


def test_explicit_missing_parallel_bundle_fails_readiness(study):
    plan = create_plan(study, {"contrast_ids": ["perturbation-vs-baseline"], "parameters": {
        "prioritize.motiftrace_parallel": {"parallel_bundle_asset_id": "missing-completed-fit"}}})
    task = next(t for t in plan["tasks"] if t["workflow_id"] == "prioritize.motiftrace_parallel")
    assert task["execution_group"] == "motiftrace"
    assert task["parameters"]["parallel_bundle_asset_id"] == "missing-completed-fit"
    assert any(e["task_id"] == task["id"] and "parallel-model bundle" in e["error"]
               for e in preflight_plan(plan)["errors"])


def test_registered_parallel_import_preflight_binds_hashes_and_exact_arms(study, tmp_path):
    from katana_next.flow.model import file_sha256
    from katana_next.infrastructure.manifest_store import ManifestStore
    store = ManifestStore()
    _, manifest = store.load(study)
    contrast = manifest.contrasts[0]
    root = tmp_path / "prior-completed-evidence"
    root.mkdir()
    member = root / "original.tsv"
    member.write_text("native_P\n0.000001\n", encoding="utf-8")
    payload = {
        "schema": "motiftrace.completed_parallel.v1", "state": "complete",
        "species": "test-species", "genome_build": manifest.annotation.genome_build,
        "new_inference": False, "methods": {"MotifTrace": {}, "MAPP": {}},
        "family": {"units_per_case": 1},
        "cases": {contrast.contrast_id: {
            "numerator": list(contrast.numerator.sample_ids),
            "denominator": list(contrast.denominator.sample_ids)}},
        "files": {"original": {"path": member.name, "bytes": member.stat().st_size,
                                "sha256": file_sha256(member)}}}
    bundle_path = write_json(root / "manifest.json", payload)
    raw = json.loads(study.read_text("utf-8"))
    raw["metadata"]["species"] = "test-species"
    raw["assets"].append({"asset_id": "prior-parallel", "kind": "motiftrace_parallel_results",
                          "path": str(bundle_path), "metadata": {"sha256": file_sha256(bundle_path)}})
    write_json(study, raw)
    plan = create_plan(study, {"contrast_ids": [contrast.contrast_id]})
    task = next(t for t in plan["tasks"] if t["workflow_id"] == "prioritize.motiftrace_parallel")
    assert not [e for e in preflight_plan(plan)["errors"] if e["task_id"] == task["id"]]
    member.write_text("changed", encoding="utf-8")
    assert any(e["task_id"] == task["id"] and "Changed completed result member" in e["error"]
               for e in preflight_plan(plan)["errors"])
    member.write_text("native_P\n0.000001\n", encoding="utf-8")
    payload["cases"][contrast.contrast_id]["numerator"] = ["different-biological-unit"]
    write_json(bundle_path, payload)
    raw["assets"][-1]["metadata"]["sha256"] = file_sha256(bundle_path)
    write_json(study, raw)
    plan = create_plan(study, {"contrast_ids": [contrast.contrast_id]})
    assert any(e["task_id"] == task["id"] and "sample/contrast identity mismatch" in e["error"]
               for e in preflight_plan(plan)["errors"])
