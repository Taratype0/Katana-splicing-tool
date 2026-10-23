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
    plan = create_plan(study, {"contrast_ids": ["perturbation-vs-baseline"]})
    validate_plan(plan)
    motifs = [t for t in plan["tasks"] if t["workflow_id"] == "prioritize.motiftrace_consequences"]
    assert len(motifs) == 1
    assert motifs[0]["scope"]["contrast_ids"] == ["perturbation-vs-baseline"]
    assert motifs[0]["selection"]["event_set"] == "stable_significant"
    enrichment = [t for t in plan["tasks"] if t["workflow_id"] == "prioritize.enrichment"]
    assert {t["parameters"]["foreground_source"] for t in enrichment} == {"das", "deg", "dtu", "dexseq"}
    assert all(t["parameters"]["background_source"] == "match_foreground" for t in enrichment)
    assert all(t["parameters"]["contrast_ids"] == ["perturbation-vs-baseline"] for t in enrichment)
    assert all(t["level"] < 12 for t in plan["tasks"])


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
    overlap = next(t for t in plan["tasks"] if t["id"] == "overlap")
    assert overlap["parameters"]["event_set"] == "significant"
    assert next(t for t in plan["tasks"] if t["id"] == "network-fit")["scope"]["unit"] == "measurable_cohort_PSI"


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
