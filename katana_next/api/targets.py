"""Portable target requests and explicitly started local Python executions."""
from __future__ import annotations

import json
import hashlib
from pathlib import Path, PurePosixPath
import re
import shutil
import threading
from typing import Any, Literal
from urllib.parse import quote
from uuid import uuid4

from fastapi import APIRouter, BackgroundTasks, HTTPException, Query
from fastapi.responses import FileResponse
from pydantic import BaseModel, ConfigDict, Field


class TargetSelectionRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")
    source_artifact_id: str
    row_ids: list[str] = Field(min_length=1)
    target_kind: Literal["gene", "event"]
    filters: dict[str, Any] = Field(default_factory=dict)


class TargetBatchRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")
    selections: list[TargetSelectionRequest] = Field(default_factory=list)
    modules: list[str] | None = None
    contrast_ids: list[str] | None = None
    request_id: str | None = None
    resources: dict[str, str] | None = None
    workers: int | None = Field(default=None, ge=1)
    target_ids: list[str] | None = None


class TargetImportRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")
    manifest_path: str


def export_targets(*args, **kwargs):
    from katana_next.targets.bundle import export_targets as export
    return export(*args, **kwargs)


def plan_targets(*args, **kwargs):
    from katana_next.targets.runner import plan_targets as plan
    return plan(*args, **kwargs)


def run_targets(*args, **kwargs):
    from katana_next.targets.runner import run_targets as run
    return run(*args, **kwargs)


def inspect_results(*args, **kwargs):
    from katana_next.targets.runner import inspect_results as inspect
    return inspect(*args, **kwargs)


def _identity(value: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}", value) or value in {".", ".."}:
        raise ValueError("Invalid target request or run identifier")
    return value


def _bounded(root: Path, path: str | Path, *, file=False) -> Path:
    root = root.resolve()
    candidate = Path(path)
    if not candidate.is_absolute():
        candidate = root / candidate
    resolved = candidate.resolve()
    if not resolved.is_relative_to(root):
        raise ValueError("Target path leaves its owned export/result directory")
    for component in (candidate, *candidate.parents):
        if component == root.parent:
            break
        if component.is_symlink():
            raise ValueError("Symlinks are not served as target outputs")
    if file and not resolved.is_file():
        raise FileNotFoundError("Target output is unavailable")
    return resolved


def _write(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".{uuid4().hex}.tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, default=str), encoding="utf-8")
    temporary.replace(path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def create_targets_router(runtime) -> APIRouter:
    """Use the same ProjectService as the parent app; do not start a service/job."""
    router = APIRouter(prefix="/api/projects/{project_id}/targets", tags=["targets"])
    submission_lock = threading.RLock()

    def context(project_id):
        try:
            return runtime.projects.get(project_id)
        except (KeyError, ValueError) as exc:
            raise HTTPException(404, "Open this project before using target reports") from exc

    def workspace(project_id):
        root = context(project_id).manifest.output_root.resolve()
        return _bounded(root, root / "target-workspace")

    def record_path(project_id, kind, identifier):
        root = workspace(project_id)
        return _bounded(root, root / kind / ".registry" / f"{_identity(identifier)}.json")

    def read_record(project_id, kind, identifier):
        path = record_path(project_id, kind, identifier)
        if not path.is_file():
            raise HTTPException(404, "Target request/run not found")
        return json.loads(path.read_text(encoding="utf-8"))

    def export(project_id, payload):
        if not payload.selections:
            raise ValueError("Add exact gene/event rows before exporting")
        root = workspace(project_id) / "exports"
        scope = root / uuid4().hex
        scope.mkdir(parents=True, exist_ok=False)
        result = export_targets(context(project_id), [s.model_dump() for s in payload.selections],
                                scope / "bundle", modules=payload.modules, contrast_ids=payload.contrast_ids)
        request_id = _identity(str(result["request_id"]))
        for key in ("request_path", "archive_path", "bundle_directory"):
            _bounded(scope, result[key], file=key != "bundle_directory")
        result = {**result, "archive_sha256": _sha256(Path(result["archive_path"])),
                  "download_url": f"/api/projects/{quote(project_id, safe='')}/targets/exports/{quote(request_id, safe='')}/download"}
        _write(record_path(project_id, "exports", request_id), result)
        return result

    def request_record(project_id, payload):
        if payload.request_id:
            if payload.selections:
                raise ValueError("Use an exported request_id or fresh selections, not both")
            record = read_record(project_id, "exports", payload.request_id)
            _bounded(workspace(project_id) / "exports", record["request_path"], file=True)
            return record
        return export(project_id, payload)

    def public_run(project_id, record):
        value = dict(record)
        result = value.get("result") or {}
        base = f"/api/projects/{quote(project_id, safe='')}/targets/runs/{quote(value['run_id'], safe='')}/files/"
        root = _bounded(workspace(project_id), value["output_directory"])
        if result:
            result = inspect_results(root)
            if (result.get("verified") is not True or result.get("project_id") != project_id
                    or result.get("request_id") != value["request_id"]):
                raise ValueError("Stored target results no longer pass project/member verification")
            value["result"] = result
        links = {}
        for name, key in (("report_url", "report"), ("archive_url", "archive"), ("manifest_url", "manifest")):
            if result.get(key):
                path = _bounded(root, result[key], file=True)
                links[name] = base + quote(path.relative_to(root).as_posix(), safe="/")
        value.update(links)
        value["execution"] = "ordinary Python on this service host; no scheduler submission"
        return value

    def execute(project_id, run_id, payload, exported, outdir):
        path = record_path(project_id, "runs", run_id)
        record = json.loads(path.read_text(encoding="utf-8"))
        record.update(status="running")
        _write(path, record)
        progress_lock = threading.Lock()
        def progress(event):
            with progress_lock:
                record["progress"] = event
                if event.get("target_id"):
                    item = dict(event)
                    if item.get("type") == "target_started":
                        item["status"] = "running"
                    record.setdefault("target_statuses", {})[str(item["target_id"])] = item
                _write(path, record)
        try:
            result = run_targets(exported["request_path"], outdir, resources=payload.resources,
                                 workers=payload.workers, target_ids=payload.target_ids, progress=progress)
            verified = inspect_results(outdir)
            if (verified.get("verified") is not True or verified.get("project_id") != project_id
                    or verified.get("request_id") != exported["request_id"]):
                raise ValueError("Generated target results failed project/member verification")
            record.update(status=result["status"], result=verified)
        except Exception as exc:
            record.update(status="failed", error=str(exc))
        _write(path, record)

    @router.post("/export")
    def export_route(project_id: str, payload: TargetBatchRequest):
        try:
            if payload.request_id:
                raise ValueError("Export needs explicit selections, not an existing request_id")
            return export(project_id, payload)
        except (ValueError, OSError, KeyError) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.post("/plan")
    def plan_route(project_id: str, payload: TargetBatchRequest):
        try:
            exported = request_record(project_id, payload)
            # Planning with a reused export reads its bundle and writes no result directory.
            plan = plan_targets(exported["request_path"], resources=payload.resources,
                                outdir=workspace(project_id) / "runs" / "products" / exported["request_id"])
            return {**plan, "request_id": exported["request_id"]}
        except (ValueError, OSError, KeyError) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.post("/run", status_code=202)
    def run_route(project_id: str, payload: TargetBatchRequest, background: BackgroundTasks):
        try:
            with submission_lock:
                exported = request_record(project_id, payload)
                options = {"resources": payload.resources, "target_ids": payload.target_ids}
                registry = workspace(project_id) / "runs" / ".registry"
                for saved in registry.glob("*.json") if registry.is_dir() else []:
                    previous = json.loads(saved.read_text(encoding="utf-8"))
                    if previous.get("request_id") == exported["request_id"] and previous.get("status") in {"queued", "running"}:
                        if previous.get("options") != options:
                            raise HTTPException(409, "This request is active with different bindings/targets; wait before changing its execution inputs")
                        return public_run(project_id, previous)
                run_id = uuid4().hex
                # Stable request-owned products allow the runner to reuse verified
                # targets when an explicit repeated run requests missing work.
                outdir = workspace(project_id) / "runs" / "products" / exported["request_id"]
                record = {"run_id": run_id, "request_id": exported["request_id"], "status": "queued",
                          "output_directory": str(outdir), "result": None, "options": options}
                _write(record_path(project_id, "runs", run_id), record)
                background.add_task(execute, project_id, run_id, payload, exported, outdir)
                return public_run(project_id, record)
        except (ValueError, OSError, KeyError) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.get("/exports/{request_id}/download")
    def download_route(project_id: str, request_id: str):
        try:
            record = read_record(project_id, "exports", request_id)
            path = _bounded(workspace(project_id) / "exports", record["archive_path"], file=True)
            if _sha256(path) != record["archive_sha256"]:
                raise ValueError("Exported target archive was modified after verification")
            return FileResponse(path, filename=f"katana-targets-{request_id}.zip", media_type="application/zip")
        except (ValueError, OSError, KeyError) as exc:
            raise HTTPException(404, str(exc)) from exc

    @router.get("/runs/{run_id}")
    def status_route(project_id: str, run_id: str):
        try:
            return public_run(project_id, read_record(project_id, "runs", run_id))
        except (ValueError, OSError, KeyError) as exc:
            raise HTTPException(404, str(exc)) from exc

    @router.get("/runs/{run_id}/files/{path:path}")
    def file_route(project_id: str, run_id: str, path: str, download: bool = Query(False)):
        try:
            relative = PurePosixPath(path.replace("\\", "/"))
            if relative.is_absolute() or ".." in relative.parts:
                raise ValueError("Invalid result-relative file path")
            record = read_record(project_id, "runs", run_id)
            root = _bounded(workspace(project_id), record["output_directory"])
            verified = inspect_results(root)
            if verified.get("verified") is not True or verified.get("project_id") != project_id:
                raise ValueError("Target result verification failed")
            allowed = set(verified.get("members", {})) | {"result-manifest.json"}
            if verified.get("archive"):
                allowed.add(Path(verified["archive"]).name)
            if relative.as_posix() not in allowed:
                raise ValueError("File is not part of the verified target result")
            resolved = _bounded(root, relative.as_posix(), file=True)
            return FileResponse(resolved, filename=resolved.name if download else None)
        except (ValueError, OSError, KeyError) as exc:
            raise HTTPException(404, str(exc)) from exc

    @router.post("/results/import")
    def import_route(project_id: str, payload: TargetImportRequest):
        try:
            source = Path(payload.manifest_path).expanduser().resolve()
            verified = inspect_results(source)
            if verified.get("verified") is not True or verified.get("project_id") != project_id:
                raise ValueError("Import must contain verified results for this project")
            source_root = Path(verified["manifest"]).resolve().parent
            run_id = uuid4().hex
            destination = workspace(project_id) / "imports" / run_id / "result"
            destination.mkdir(parents=True, exist_ok=False)
            names = set(verified["members"]) | {"result-manifest.json"}
            if verified.get("archive"):
                names.add(Path(verified["archive"]).name)
            for name in sorted(names):
                original = _bounded(source_root, name, file=True)
                target = _bounded(destination, name)
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(original, target)
            copied = inspect_results(destination)
            if copied.get("verified") is not True or copied.get("project_id") != project_id:
                raise ValueError("Copied import failed verification")
            record = {"run_id": run_id, "request_id": copied["request_id"], "status": copied["status"],
                      "output_directory": str(destination), "result": copied, "imported": True}
            _write(record_path(project_id, "runs", run_id), record)
            return public_run(project_id, record)
        except (ValueError, OSError, KeyError) as exc:
            raise HTTPException(422, str(exc)) from exc

    return router
