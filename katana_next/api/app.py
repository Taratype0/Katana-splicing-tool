from __future__ import annotations

import asyncio
import json
import os
from pathlib import Path
from typing import Any

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, StreamingResponse, JSONResponse
from fastapi.staticfiles import StaticFiles

from katana_next import __version__
from katana_next.api.schemas import (
    CreateSelectionRequest,
    ConfigureResourcesRequest,
    ConfigureGeneFilterRequest,
    MigrateProjectRequest,
    OpenProjectRequest,
    RunWorkflowRequest,
)
from katana_next.domain.models import ArtifactKind, JobStatus, SelectionSet
from katana_next.infrastructure.jobs import JobManager
from katana_next.infrastructure.catalog import CatalogBusyError
from katana_next.migration import migrate_legacy_project
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

    @app.exception_handler(CatalogBusyError)
    async def catalog_busy(_request, exc: CatalogBusyError):
        return JSONResponse(status_code=409, content={"detail": str(exc)})
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

    @app.get("/api/workflows")
    def workflows() -> list[dict[str, object]]:
        return active.workflows.definitions()

    @app.get("/api/motiftrace/modes")
    def motiftrace_modes() -> dict[str, Any]:
        """One read-only public catalog; execution uses the workflow routes."""
        from katana_motiftrace import mode_catalog
        return mode_catalog()

    @app.post("/api/projects/open")
    def open_project(request: OpenProjectRequest) -> dict[str, Any]:
        try:
            context = active.projects.open_project(request.manifest_path)
            active.jobs.bind_project(context.manifest.project_id, context.manifest.output_root)
        except (FileNotFoundError, ValueError, json.JSONDecodeError) as exc:
            raise HTTPException(status_code=422, detail=str(exc)) from exc
        return _project_payload(context)

    @app.post("/api/projects/migrate")
    def migrate_project(request: MigrateProjectRequest) -> dict[str, Any]:
        try:
            manifest_path = migrate_legacy_project(
                request.legacy_config_path,
                request.output_root,
                request.manifest_path,
            )
            context = active.projects.open_project(manifest_path)
            active.jobs.bind_project(context.manifest.project_id, context.manifest.output_root)
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

    @app.post("/api/projects/configure-resources")
    def configure_resources(request: ConfigureResourcesRequest) -> dict[str, Any]:
        from katana_next.resources import attach_resources
        try:
            root = Path(request.output_root).expanduser().resolve()
            active.projects._validate_output_root(root)
            if root.exists() and any(root.iterdir()):
                raise ValueError("Choose a new empty output directory; existing projects are preserved.")
            manifest_path = attach_resources(Path(request.manifest_path), Path(request.resource_bundle_path),
                                             root / "katana.project.json", root / "results")
            context = active.projects.open_project(manifest_path)
            active.jobs.bind_project(context.manifest.project_id, context.manifest.output_root)
        except (OSError, ValueError, KeyError) as exc:
            raise HTTPException(status_code=422, detail=str(exc)) from exc
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

    @app.post("/api/projects/{project_id}/gene-filter")
    def configure_gene_filter(project_id: str, request: ConfigureGeneFilterRequest):
        from katana_next.services.gene_filter_service import configure_gene_filter
        context = _context_or_404(active, project_id)
        try:
            filtered = configure_gene_filter(active.projects, context, request)
            active.jobs.bind_project(filtered.manifest.project_id, filtered.manifest.output_root)
        except (OSError, ValueError) as exc:
            raise HTTPException(status_code=422, detail=str(exc)) from exc
        return _project_payload(filtered)

    @app.get("/api/projects/{project_id}/design")
    def project_design(project_id: str):
        from katana_next.services.readiness import experimental_design
        return experimental_design(_context_or_404(active, project_id).manifest)

    @app.get("/api/projects/{project_id}/design/edit")
    def design_editor(project_id: str):
        from katana_next.services.design_editor import editable_design
        return editable_design(_context_or_404(active, project_id))

    @app.post("/api/projects/{project_id}/design/contrasts")
    def contrast_plan(project_id: str, request: dict[str, Any]):
        from katana_next.science.contrast_planner import plan_contrasts
        try:
            if set(request) - {"factors", "include_batch", "pairing", "reference_levels"}:
                raise ValueError("Unknown contrast planner option")
            if "include_batch" in request and not isinstance(request["include_batch"], bool):
                raise ValueError("include_batch must be boolean")
            if "factors" in request and (not isinstance(request["factors"], list) or not all(isinstance(f,str) for f in request["factors"])):
                raise ValueError("factors must be a list of declared factor names")
            return plan_contrasts(_context_or_404(active, project_id).manifest, **request)
        except ValueError as exc:
            raise HTTPException(status_code=422, detail=str(exc)) from exc

    @app.post("/api/projects/{project_id}/design/edit")
    def update_design(project_id: str, request: dict[str, Any], apply: bool = False):
        from katana_next.services.design_editor import revise_design
        _context_or_404(active, project_id)
        try:
            with active.workflows._submission_lock:
                if any(j.status not in TERMINAL_JOB_STATES for j in active.jobs.list(project_id)):
                    raise HTTPException(status_code=409, detail="Finish or cancel active project jobs before editing design")
                result = revise_design(active.projects, project_id, request, apply=apply)
                if apply:
                    result["project"] = _project_payload(active.projects.get(project_id))
                return result
        except (ValueError, OSError) as exc:
            raise HTTPException(status_code=422, detail=str(exc)) from exc

    @app.get("/api/projects/{project_id}/artifacts/{artifact_id}/module-evidence")
    def inspect_module(project_id: str, artifact_id: str, module: str, limit: int = Query(default=30, ge=1, le=100)):
        from katana_next.services.module_browser import module_evidence
        try:
            return module_evidence(_context_or_404(active, project_id), artifact_id, module, limit=limit)
        except (ValueError, OSError) as exc:
            raise HTTPException(status_code=422, detail=str(exc)) from exc

    @app.get("/api/projects/{project_id}/artifacts/{artifact_id}/module-dashboard")
    def inspect_module_dashboard(project_id: str, artifact_id: str, module: str, condition: str = "all",
                                 selected_event: str | None = None, hide_same_gene: bool = False,
                                 node_limit: int = Query(default=80, ge=10, le=200)):
        from katana_next.services.module_dashboard import module_dashboard
        try:
            return module_dashboard(_context_or_404(active, project_id), artifact_id, module,
                                    condition=condition, selected_event=selected_event,
                                    hide_same_gene=hide_same_gene, node_limit=node_limit)
        except (ValueError, OSError, KeyError) as exc:
            raise HTTPException(status_code=422, detail=str(exc)) from exc

    @app.get("/api/projects/{project_id}/artifacts/{artifact_id}/lineage")
    def lineage(project_id: str, artifact_id: str):
        from katana_next.services.lineage import artifact_lineage
        try:
            return artifact_lineage(_context_or_404(active,project_id),artifact_id)
        except ValueError as exc:
            raise HTTPException(status_code=422,detail=str(exc)) from exc

    @app.post("/api/projects/{project_id}/workflows/{workflow_id}/plan")
    def plan_workflow(project_id: str, workflow_id: str, request: RunWorkflowRequest):
        from katana_next.services.readiness import preflight
        _context_or_404(active, project_id)
        try:
            return preflight(active.workflows, project_id, workflow_id, request.parameters)
        except KeyError as exc:
            raise HTTPException(status_code=404, detail=str(exc)) from exc

    @app.get("/api/projects/{project_id}/single-cell/conditions")
    def single_cell_conditions(project_id: str) -> list[str]:
        from katana_next.science.single_cell import load_single_cell_metadata
        context = _context_or_404(active, project_id)
        try:
            metadata = load_single_cell_metadata(context)
        except (ValueError, FileNotFoundError) as exc:
            raise HTTPException(status_code=422, detail=str(exc)) from exc
        return sorted(metadata["condition"].dropna().astype(str).unique().tolist())

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

    @app.get("/api/projects/{project_id}/artifacts/{artifact_id}/file")
    def artifact_file(
        project_id: str, artifact_id: str, download: bool = False
    ) -> FileResponse:
        context = _context_or_404(active, project_id)
        record = context.catalog.get_artifact(artifact_id)
        if record is None or record.project_id != project_id:
            raise HTTPException(status_code=404, detail="Artifact not found")
        if not record.path.is_file():
            raise HTTPException(status_code=404, detail="Artifact file is missing")
        media_type = str(record.metadata.get("mime_type", "")) or None
        return FileResponse(
            record.path,
            media_type=media_type,
            filename=record.path.name if download else None,
        )

    @app.get("/api/projects/{project_id}/artifacts/{artifact_id}/export/{file_format}")
    def artifact_export(project_id: str, artifact_id: str, file_format: str) -> FileResponse:
        context = _context_or_404(active, project_id)
        record = context.catalog.get_artifact(artifact_id)
        if record is None or record.project_id != project_id or file_format not in {"png", "svg", "tsv"}:
            raise HTTPException(status_code=404, detail="Plot export not found")
        name = record.metadata.get("exports", {}).get(file_format)
        if not name:
            raise HTTPException(status_code=404, detail="No cached static export; rebuild this plot explicitly")
        path = Path(name).resolve()
        if not path.is_relative_to(context.manifest.output_root.resolve()) or not path.is_file():
            raise HTTPException(status_code=404, detail="Cached export is unavailable")
        return FileResponse(path, filename=path.name)

    @app.get("/api/projects/{project_id}/selections")
    def selections(project_id: str) -> list[dict[str, Any]]:
        context = _context_or_404(active, project_id)
        current = []
        from katana_next.workflows.sources import source_complete
        for item in context.catalog.list_selections(project_id):
            source = context.catalog.get_artifact(item.source_artifact_id)
            if source is None or not source.path.is_file() or not source_complete(context, source):
                continue
            if (
                source.metadata.get("project_snapshot_fingerprint")
                != context.snapshot.manifest_fingerprint
            ):
                continue
            from katana_next.services.selection_routes import compatible_workflows
            current.append({**item.model_dump(mode="json"), "compatible_workflows": compatible_workflows(source)})
        return current

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

    @app.post("/api/projects/{project_id}/workflows/{workflow_id}/run")
    def run_workflow(
        project_id: str, workflow_id: str, request: RunWorkflowRequest
    ) -> dict[str, Any]:
        try:
            job = active.workflows.run(
                workflow_id,
                project_id,
                parameters=request.parameters,
                force=request.force,
            )
        except (KeyError, ValueError, FileNotFoundError) as exc:
            raise HTTPException(status_code=422, detail=str(exc)) from exc
        return job.model_dump(mode="json")

    @app.post("/api/projects/{project_id}/workflows/{workflow_id}/cache")
    def workflow_cache(
        project_id: str, workflow_id: str, request: RunWorkflowRequest
    ) -> dict[str, Any]:
        try:
            run_fingerprint, cached = active.workflows.compatible_cache(
                workflow_id,
                project_id,
                parameters=request.parameters,
            )
        except (KeyError, ValueError, FileNotFoundError) as exc:
            raise HTTPException(status_code=422, detail=str(exc)) from exc
        return {
            "run_fingerprint": run_fingerprint,
            "cache_hit": bool(cached),
            "artifacts": [artifact.model_dump(mode="json") for artifact in cached],
        }

    @app.post("/api/projects/{project_id}/open-output")
    def open_output(project_id: str) -> dict[str, Any]:
        context = _context_or_404(active, project_id)
        output_root = context.manifest.output_root.resolve()
        output_root.mkdir(parents=True, exist_ok=True)
        try:
            os.startfile(output_root)  # type: ignore[attr-defined]
        except (AttributeError, OSError) as exc:
            raise HTTPException(
                status_code=501,
                detail=f"Could not open the local output folder: {exc}",
            ) from exc
        return {"opened": True, "path": str(output_root)}

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
        except ValueError as exc:
            raise HTTPException(status_code=409, detail=str(exc)) from exc
        return record.model_dump(mode="json")

    @app.post("/api/jobs/{job_id}/retry")
    def retry_job(job_id: str):
        try:
            return active.workflows.retry(job_id).model_dump(mode="json")
        except (ValueError, KeyError) as exc:
            raise HTTPException(status_code=409, detail=str(exc)) from exc

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
