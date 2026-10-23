"""Sealed-plan reuse authenticates old identities without regenerating tasks."""
from pathlib import Path

import pytest

from katana_next.flow.model import file_sha256, read_json, seal_plan, write_json
from katana_next.flow.plan_reuse import validate_existing_plan


@pytest.fixture
def original(tmp_path):
    project = tmp_path / "original" / "project.json"
    write_json(project, {"project_id": "generated-plan-transport-fixture"})
    plan = seal_plan({"project_manifest": str(project),
        "project_manifest_sha256": file_sha256(project), "configuration": {},
        "tasks": [{"id": "audit", "workflow_id": "quality.project_audit", "dependencies": [],
                   "parameters": {}}]})
    target = write_json(tmp_path / "original" / "plan.json", plan)
    return target, project, plan


def _ready(plan):
    return {"plan_id": plan["plan_id"], "ready": True, "errors": [], "warnings": []}


def test_original_plan_path_and_task_parameters_are_preserved_for_staged_project(original, tmp_path, monkeypatch):
    from katana_next.flow import planner
    path, project, sealed = original
    staged = tmp_path / "staged-project.json"
    staged.write_bytes(project.read_bytes())
    before = {p: file_sha256(p) for p in (path, project, staged)}
    observed = []
    def checked(plan):
        observed.append(plan)
        return _ready(plan)
    def forbidden(*args, **kwargs):
        raise AssertionError("Reuse may not regenerate a plan")
    monkeypatch.setattr(planner, "create_plan", forbidden)
    monkeypatch.setattr(planner, "preflight_plan", checked)
    result = validate_existing_plan(path, staged)
    assert result["ready"] and result["replanned"] is False
    assert result["plan_id"] == sealed["plan_id"]
    assert result["plan_path"] == str(path.resolve())
    assert result["project_manifest"] == str(staged.resolve())
    assert observed == [sealed]
    assert before == {p: file_sha256(p) for p in before}


@pytest.mark.parametrize("change", ["task", "current-project", "original-project", "missing-original"])
def test_modified_identity_and_missing_frozen_project_fail_before_preflight(original, tmp_path, monkeypatch, change):
    from katana_next.flow import planner
    path, project, plan = original
    staged = tmp_path / "current.json"
    staged.write_bytes(project.read_bytes())
    if change == "task":
        plan["tasks"][0]["parameters"]["silently_added_default"] = "JCEC"
        write_json(path, plan)
    elif change == "current-project":
        write_json(staged, {"project_id": "different-project"})
    elif change == "original-project":
        write_json(project, {"project_id": "different-original"})
    else:
        project.unlink()
    def forbidden(*args, **kwargs):
        raise AssertionError("Invalid identity must fail before runtime readiness")
    monkeypatch.setattr(planner, "preflight_plan", forbidden)
    with pytest.raises((ValueError, FileNotFoundError)):
        validate_existing_plan(path, staged)


def test_failed_preflight_is_retained_and_cli_returns_nonzero(original, tmp_path, monkeypatch, capsys):
    from katana_next.flow import cli, planner
    path, project, plan = original
    failed = {"plan_id": plan["plan_id"], "ready": False,
              "errors": [{"task_id": "network", "error": "Frozen input path is unavailable"}], "warnings": []}
    monkeypatch.setattr(planner, "preflight_plan", lambda _: failed)
    output = tmp_path / "reuse.json"
    assert cli.main(["validate-plan", "--plan", str(path), "--project", str(project), "--out", str(output)]) == 2
    assert read_json(output)["preflight"] == failed
    assert read_json(output)["ready"] is False
    before = output.read_bytes()
    assert cli.main(["validate-plan", "--plan", str(path), "--project", str(project), "--out", str(output)]) == 1
    assert output.read_bytes() == before


def test_plan_mutation_during_preflight_is_not_admitted(original, monkeypatch):
    from katana_next.flow import planner
    path, project, plan = original
    def mutate(current):
        path.write_bytes(path.read_bytes() + b"\n")
        return _ready(current)
    monkeypatch.setattr(planner, "preflight_plan", mutate)
    with pytest.raises(ValueError, match="changed during reuse"):
        validate_existing_plan(path, project)


def test_success_cli_writes_only_reuse_receipt_and_keeps_plan_bytes(original, tmp_path, monkeypatch):
    from katana_next.flow import cli, planner
    path, project, plan = original
    before = path.read_bytes()
    monkeypatch.setattr(planner, "preflight_plan", _ready)
    target = tmp_path / "checked.json"
    assert cli.main(["validate-plan", "--plan", str(path), "--project", str(project), "--out", str(target)]) == 0
    assert read_json(target)["plan_id"] == plan["plan_id"]
    assert read_json(target)["plan_sha256"] == file_sha256(path)
    assert path.read_bytes() == before
