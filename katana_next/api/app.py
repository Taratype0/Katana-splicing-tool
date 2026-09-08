from __future__ import annotations

import asyncio
import json
from pathlib import Path
from typing import Any

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles

from katana_next import __version__
from katana_next.api.schemas import (
    CreateSelectionRequest,
    OpenProjectRequest,
    RunWorkflowRequest,
)
from katana_next.domain.models import ArtifactKind, JobStatus, SelectionSet
from katana_next.infrastructure.jobs import JobManager
from katana_next.services.project_service import ProjectService
from katana_next.services.workflow_service import WorkflowService
from katana_next.workflows.catalog import WORKSPACES


TERMINAL_JOB_STATES = {
    JobStatus.COMPLETE,
    JobStatus.FAILED,
    JobStatus.CANCELLED,
    JobStatus.CACHE_HIT,
}


class PlatformRuntime:
    def __init__(self, application_root: str | Path | None = None) -> None:
        self.projects = ProjectService(application_root=application_root)
        self.jobs = JobManager(max_workers=2)
        self.workflows = WorkflowService(self.projects, self.jobs)


def create_app(runtime: PlatformRuntime | None = None) -> FastAPI:
    active = runtime or PlatformRuntime()
    app = FastAPI(
        title="Katana Next Local API",
        version=__version__,
        description="Typed local API for Katana analysis artifacts and explicit jobs.",
    )
    app.state.runtime = active
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["http://127.0.0.1:4173", "http://127.0.0.1:5173"],
        allow_credentials=False,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    @app.get("/api/health")
    def health() -> dict[str, str]:
        return {"status": "ok", "version": __version__}

    @app.get("/api/workspaces")
    def workspaces() -> list[dict[str, object]]:
        return [workspace.to_dict() for workspace in WORKSPACES]

    @app.post("/api/projects/open")
    def open_project(request: OpenProjectRequest) -> dict[str, Any]:
        try:
            context = active.projects.open_project(request.manifest_path)
        except (FileNotFoundError, ValueError, json.JSONDecodeError) as exc:
            raise HTTPException(status_code=422, detail=str(exc)) from exc
        return _project_payload(context)

    @app.get("/api/projects/{project_id}")
    def project(project_id: str) -> dict[str, Any]:
        try:
            context = active.projects.get(project_id)
        except KeyError as exc:
            raise HTTPException(status_code=404, detail=str(exc)) from exc
        return _project_payload(context)

    @app.get("/api/projects/{project_id}/artifacts")
    def artifacts(
        project_id: str,
        analysis_id: str | None = None,
        kind: str | None = None,
    ) -> list[dict[str, Any]]:
        try:
            records = active.projects.list_artifacts(
                project_id, analysis_id=analysis_id, kind=kind
            )
        except KeyError as exc:
            raise HTTPException(status_code=404, detail=str(exc)) from exc
        return [record.model_dump(mode="json") for record in records]

    @app.get("/api/projects/{project_id}/artifacts/{artifact_id}")
    def artifact(project_id: str, artifact_id: str) -> dict[str, Any]:
        context = _context_or_404(active, project_id)
        record = context.catalog.get_artifact(artifact_id)
        if record is None or record.project_id != project_id:
            raise HTTPException(status_code=404, detail="Artifact not found")
        return record.model_dump(mode="json")

    @app.get("/api/projects/{project_id}/artifacts/{artifact_id}/rows")
    def artifact_rows(
        project_id: str,
        artifact_id: str,
        offset: int = Query(default=0, ge=0),
        limit: int = Query(default=100, ge=1, le=500),
        search: str = "",
        sort_by: str | None = None,
        descending: bool = False,
    ) -> dict[str, Any]:
        context = _context_or_404(active, project_id)
        record = context.catalog.get_artifact(artifact_id)
        if record is None or record.project_id != project_id:
            raise HTTPException(status_code=404, detail="Artifact not found")
        if record.kind not in {ArtifactKind.TABLE, ArtifactKind.SELECTION}:
            raise HTTPException(status_code=409, detail="Artifact is not a table")
        try:
            return context.results.table_page(
                record.path,
                offset=offset,
                limit=limit,
                search=search,
                sort_by=sort_by,
                descending=descending,
            )
        except (FileNotFoundError, ValueError) as exc:
            raise HTTPException(status_code=422, detail=str(exc)) from exc

    @app.get("/api/projects/{project_id}/artifacts/{artifact_id}/content")
    def artifact_content(project_id: str, artifact_id: str) -> dict[str, Any]:
        context = _context_or_404(active, project_id)
        record = context.catalog.get_artifact(artifact_id)
        if record is None or record.project_id != project_id:
            raise HTTPException(status_code=404, detail="Artifact not found")
        if record.kind != ArtifactKind.PLOT_SPEC:
            raise HTTPException(
                status_code=409, detail="Artifact is not a plot specification"
            )
        return context.results.read_plot_spec(record.path)

    @app.get("/api/projects/{project_id}/selections")
    def selections(project_id: str) -> list[dict[str, Any]]:
        context = _context_or_404(active, project_id)
        return [
            item.model_dump(mode="json")
            for item in context.catalog.list_selections(project_id)
        ]

    @app.post("/api/projects/{project_id}/selections")
    def create_selection(
        project_id: str, request: CreateSelectionRequest
    ) -> dict[str, Any]:
        _context_or_404(active, project_id)
        if not request.entity_ids:
            raise HTTPException(status_code=422, detail="Select at least one entity")
        selection = SelectionSet(
            project_id=project_id,
            name=request.name.strip() or "Untitled selection",
            entity_type=request.entity_type,
            source_artifact_id=request.source_artifact_id,
            entity_ids=tuple(dict.fromkeys(request.entity_ids)),
            filters=request.filters,
        )
        try:
            saved = active.projects.save_selection(selection)
        except ValueError as exc:
            raise HTTPException(status_code=422, detail=str(exc)) from exc
        return saved.model_dump(mode="json")

    @app.post("/api/projects/{project_id}/workflows/phase1-profile/run")
    def run_phase1_profile(
        project_id: str, request: RunWorkflowRequest
    ) -> dict[str, Any]:
        try:
            job = active.workflows.run_phase1_profile(
                project_id, parameters=request.parameters, force=request.force
            )
        except (KeyError, ValueError) as exc:
            raise HTTPException(status_code=422, detail=str(exc)) from exc
        return job.model_dump(mode="json")

    @app.get("/api/jobs")
    def jobs(project_id: str | None = None) -> list[dict[str, Any]]:
        return [job.model_dump(mode="json") for job in active.jobs.list(project_id)]

    @app.get("/api/jobs/{job_id}")
    def job(job_id: str) -> dict[str, Any]:
        record = active.jobs.get(job_id)
        if record is None:
            raise HTTPException(status_code=404, detail="Job not found")
        return record.model_dump(mode="json")

    @app.post("/api/jobs/{job_id}/cancel")
    def cancel_job(job_id: str) -> dict[str, Any]:
        try:
            record = active.jobs.cancel(job_id)
        except KeyError as exc:
            raise HTTPException(status_code=404, detail="Job not found") from exc
        return record.model_dump(mode="json")

    @app.get("/api/jobs/{job_id}/events")
    async def job_events(
        job_id: str, after: int = Query(default=0, ge=0)
    ) -> StreamingResponse:
        if active.jobs.get(job_id) is None:
            raise HTTPException(status_code=404, detail="Job not found")

        async def stream():
            cursor = after
            while True:
                record = active.jobs.get(job_id)
                if record is None:
                    return
                for event in record.events:
                    if event.sequence <= cursor:
                        continue
                    cursor = event.sequence
                    yield "event: progress\n"
                    yield f"data: {event.model_dump_json()}\n\n"
                if record.status in TERMINAL_JOB_STATES:
                    yield "event: complete\n"
                    yield f"data: {record.model_dump_json()}\n\n"
                    return
                await asyncio.sleep(0.2)

        return StreamingResponse(stream(), media_type="text/event-stream")

    frontend = Path(__file__).resolve().parents[2] / "apps" / "web" / "dist"
    if frontend.is_dir():
        assets = frontend / "assets"
        if assets.is_dir():
            app.mount("/assets", StaticFiles(directory=assets), name="assets")

        @app.get("/{path:path}", include_in_schema=False)
        def frontend_route(path: str) -> FileResponse:
            requested = (frontend / path).resolve()
            if path and requested.is_file() and frontend in requested.parents:
                return FileResponse(requested)
            return FileResponse(frontend / "index.html")

    return app


def _context_or_404(active: PlatformRuntime, project_id: str):
    try:
        return active.projects.get(project_id)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc


def _project_payload(context) -> dict[str, Any]:
    return {
        "snapshot": context.snapshot.model_dump(mode="json"),
        "manifest": context.manifest.model_dump(mode="json"),
        "artifact_count": len(
            context.catalog.list_artifacts(context.manifest.project_id)
        ),
        "selection_count": len(
            context.catalog.list_selections(context.manifest.project_id)
        ),
    }
