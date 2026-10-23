"""Validate a sealed plan for reuse without replanning or moving its inputs."""
from pathlib import Path

from .model import file_sha256, load_plan


def validate_existing_plan(plan_path, project_path):
    """Check the original plan bytes, current project and frozen input paths.

    The current manifest may be a byte-identical staged copy. Preflight uses
    the original plan and its bound paths; it must never call the planner or
    silently replace historical task parameters with newer defaults.
    """
    from .planner import preflight_plan

    plan_path, project_path = Path(plan_path).resolve(), Path(project_path).resolve()
    plan_sha = file_sha256(plan_path)
    plan = load_plan(plan_path)
    expected = plan["project_manifest_sha256"]
    current_sha = file_sha256(project_path)
    if current_sha != expected:
        raise ValueError("Current project manifest differs from the original sealed plan")
    original_project = Path(plan["project_manifest"])
    if file_sha256(original_project) != expected:
        raise ValueError("The plan-bound project manifest changed; original frozen input paths must remain valid")
    checked = preflight_plan(plan)
    if (file_sha256(plan_path) != plan_sha or file_sha256(project_path) != current_sha
            or file_sha256(original_project) != expected):
        raise ValueError("The original plan or project changed during reuse validation")
    return {"schema": "katana.flow.plan_reuse/v1", "ready": checked["ready"],
            "plan_id": plan["plan_id"], "plan_path": str(plan_path), "plan_sha256": plan_sha,
            "project_manifest": str(project_path), "project_manifest_sha256": current_sha,
            "bound_project_manifest": str(original_project), "replanned": False,
            "preflight": checked, "errors": checked["errors"], "warnings": checked.get("warnings", []),
            "scope": "Original sealed plan and current manifest verified; original task parameters and paths retained"}
