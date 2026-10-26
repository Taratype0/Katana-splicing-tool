"""Portable, content-addressed execution plans; no project/database mutations."""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

PLAN_SCHEMA = "katana.flow.plan/v1"
MAX_LEVELS = 12


def canonical_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     ensure_ascii=True, allow_nan=False).encode()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=True, allow_nan=False) + "\n", encoding="utf-8")
    return path


def file_sha256(path):
    from .file_io import digest
    return digest(path)


def assign_levels(tasks):
    by_id = {t["id"]: t for t in tasks}
    if len(by_id) != len(tasks):
        raise ValueError("Execution task IDs must be unique")
    active, done = set(), {}

    def visit(identity):
        if identity not in by_id:
            raise ValueError(f"Unknown dependency: {identity}")
        if identity in active:
            raise ValueError("Execution dependencies contain a cycle")
        if identity in done:
            return done[identity]
        active.add(identity)
        level = max((visit(d) + 1 for d in by_id[identity].get("dependencies", [])), default=0)
        active.remove(identity)
        if level >= MAX_LEVELS:
            raise ValueError(f"Execution graph exceeds {MAX_LEVELS} dependency levels")
        done[identity] = level
        return level

    return [{**t, "level": visit(t["id"])} for t in tasks]


def seal_plan(payload):
    payload = {**payload, "schema": PLAN_SCHEMA, "tasks": assign_levels(payload["tasks"])}
    payload.pop("plan_id", None)
    return {**payload, "plan_id": canonical_hash(payload)}


def validate_plan(plan):
    if plan.get("schema") != PLAN_SCHEMA or not isinstance(plan.get("tasks"), list) or not plan["tasks"]:
        raise ValueError("Expected a nonempty Katana flow plan")
    for task in plan["tasks"]:
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,149}", task.get("id", "")):
            raise ValueError("Task IDs must be safe, bounded directory names")
        if not isinstance(task.get("workflow_id"), str) or not task["workflow_id"]:
            raise ValueError("Every task requires a workflow or typed adapter identity")
        if len(task.get("dependencies", [])) != len(set(task.get("dependencies", []))):
            raise ValueError("Duplicate dependency")
    expected = seal_plan(plan)
    if expected != plan:
        raise ValueError("Plan identity, dependency levels or content changed; create a new plan")
    if not plan.get("project_manifest") or not plan.get("project_manifest_sha256"):
        raise ValueError("The exact project manifest must be bound to the plan")
    return plan


def load_plan(path):
    return validate_plan(read_json(path))

