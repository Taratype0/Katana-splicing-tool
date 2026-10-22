from __future__ import annotations

import json
import threading
from pathlib import Path
from typing import Any

from katana_next.domain.models import AnalysisArtifact, FileFingerprint, JobRecord, JobStatus
from katana_next.infrastructure.fingerprints import (
    fingerprint_file,
    fingerprint_value,
    workflow_fingerprint,
)
from katana_next.infrastructure.jobs import JobManager
from katana_next.services.project_service import ProjectService
from katana_next.workflows.phase1_profile import Phase1ProfileWorkflow
from katana_next.workflows.registry import WorkflowRegistry


class WorkflowService:
    def __init__(
        self,
        projects: ProjectService,
        jobs: JobManager,
        registry: WorkflowRegistry | None = None,
    ) -> None:
        self.projects = projects
        self.jobs = jobs
        self.registry = registry or WorkflowRegistry()
        self.phase1 = Phase1ProfileWorkflow()
        self._fingerprint_cache: dict[tuple[str, int, int], FileFingerprint] = {}
        self._fingerprint_lock = threading.RLock()
        self._submission_lock = threading.RLock()

    def definitions(self) -> list[dict[str, object]]:
        return self.registry.definitions()

    def run(
        self,
        workflow_id: str,
        project_id: str,
        *,
        parameters: dict[str, Any] | None = None,
        force: bool = False,
    ) -> JobRecord:
        with self._submission_lock:
            return self._run(workflow_id, project_id, parameters=parameters, force=force)

    def _run(self, workflow_id: str, project_id: str, *, parameters: dict[str, Any] | None = None, force: bool = False) -> JobRecord:
        self.jobs.bind_project(project_id, self.projects.get(project_id).manifest.output_root, execute=True)
        context, workflow, normalized, run_fingerprint, cached = self._resolve(
            workflow_id, project_id, parameters
        )
        active = next((job for job in self.jobs.list(project_id)
                       if job.analysis_id == workflow_id and job.run_fingerprint == run_fingerprint
                       and job.status in {JobStatus.QUEUED, JobStatus.RUNNING}), None)
        if active is not None:
            return active
        if not force and len(cached) >= workflow.definition.expected_artifacts:
            return self.jobs.cache_hit(
                project_id=project_id,
                analysis_id=workflow.definition.workflow_id,
                run_fingerprint=run_fingerprint,
                artifact_ids=[artifact.artifact_id for artifact in cached],
                parameters=normalized,
            )
        context.catalog.clear_run_artifacts(project_id, workflow.definition.workflow_id, run_fingerprint)
        context.results.write_plot_spec({"completed": False, "workflow_id": workflow_id, "run_fingerprint": run_fingerprint},
                                        context.results.artifact_directory(workflow_id, run_fingerprint) / "parameters.json")
        return self.jobs.submit(
            project_id=project_id,
            analysis_id=workflow.definition.workflow_id,
            run_fingerprint=run_fingerprint,
            runner=lambda job: workflow.run(
                context, job, run_fingerprint, normalized
            ),
            parameters=normalized,
        )

    def retry(self, job_id):
        with self._submission_lock:
            job = self.jobs.get(job_id)
            if job is None or job.status not in {JobStatus.FAILED, JobStatus.CANCELLED}:
                raise ValueError("Only failed or cancelled recorded jobs can be retried")
            if job.analysis_id == self.phase1.workflow_id:
                raise ValueError("The legacy Phase 1 prototype has no resumable contract. Use its explicit Run action instead.")
            _, _, _, fingerprint, _ = self._resolve(job.analysis_id, job.project_id, job.parameters)
            if fingerprint != job.run_fingerprint:
                raise ValueError("Inputs, source results, parameters or method version changed. Start a new explicit analysis instead of resuming stale work.")
            return self.run(job.analysis_id, job.project_id, parameters=job.parameters)

    def compatible_cache(
        self,
        workflow_id: str,
        project_id: str,
        *,
        parameters: dict[str, Any] | None = None,
    ) -> tuple[str, list[AnalysisArtifact]]:
        """Resolve the exact current run without submitting a job."""
        _context, _workflow, _normalized, run_fingerprint, cached = self._resolve(
            workflow_id, project_id, parameters
        )
        return run_fingerprint, cached

    def _resolve(
        self,
        workflow_id: str,
        project_id: str,
        parameters: dict[str, Any] | None,
    ):
        context = self.projects.get(project_id)
        if context.manifest.metadata.get("gene_filter_identity"):
            from katana_next.science.gene_filter import validate_gene_filter_identity
            validate_gene_filter_identity(context.manifest)
        workflow = self.registry.get(workflow_id)
        workflow.validate_assets(context)
        normalized = workflow.normalize_parameters(parameters)
        from katana_next.services.selection_routes import validate_route
        validate_route(context, workflow_id, normalized)
        normalized = self._normalize_contrast_scope(context, workflow, normalized)
        input_paths = list(
            dict.fromkeys(
                Path(path).resolve()
                for path in workflow.input_paths(context, normalized)
                if Path(path).exists()
            )
        )
        parents = workflow.parent_artifacts(context, normalized)
        normalized["source_artifact_ids"] = {
            f"{artifact.analysis_id}::{artifact.name}": artifact.artifact_id
            for artifact in parents
        }
        input_paths.extend(
            artifact.path.resolve()
            for artifact in parents
            if artifact.path.is_file() and artifact.path.resolve() not in input_paths
        )
        run_fingerprint = workflow_fingerprint(
            workflow_id=workflow.definition.workflow_id,
            workflow_version=workflow.version,
            project_fingerprint=fingerprint_value(
                {
                    "schema": context.manifest.schema_version,
                    "project_id": context.manifest.project_id,
                    "manifest": context.manifest.model_dump(mode="json"),
                }
            ),
            inputs=[self._fingerprint(path) for path in input_paths],
            parameters={
                "artifact_contract_version": "3.1-na-preserving",
                **normalized,
                "parent_artifact_fingerprints": [
                    artifact.fingerprint for artifact in parents
                ],
            },
        )
        cached = context.catalog.cached_artifacts(
            project_id, workflow.definition.workflow_id, run_fingerprint
        )
        if cached:
            completion = context.results.artifact_directory(workflow_id, run_fingerprint) / "parameters.json"
            try:
                marker = json.loads(completion.read_text(encoding="utf-8"))
                expected_ids = {row["artifact_id"] for row in marker["generated_artifacts"]}
                exports = marker.get("generated_plot_files", []) + marker.get("generated_table_files", [])
                if (marker.get("completed") is not True or marker.get("run_fingerprint") != run_fingerprint
                        or expected_ids != {item.artifact_id for item in cached}
                        or not all(Path(path).is_file() for path in exports)):
                    cached = []
            except (OSError, ValueError, KeyError, TypeError):
                cached = []
        return context, workflow, normalized, run_fingerprint, cached

    @staticmethod
    def _normalize_contrast_scope(
        context, workflow, parameters: dict[str, Any]
    ) -> dict[str, Any]:
        normalized = dict(parameters)
        if "contrast_ids" not in normalized:
            return normalized
        if not workflow.definition.supports_contrast_selection:
            raise ValueError(
                f"{workflow.definition.label} does not support contrast subset selection."
            )
        raw = normalized["contrast_ids"]
        if isinstance(raw, str):
            requested = [raw]
        elif isinstance(raw, (list, tuple, set)):
            requested = [str(item) for item in raw]
        else:
            raise ValueError("contrast_ids must be a list of declared contrast IDs.")
        known = [contrast.contrast_id for contrast in context.manifest.contrasts]
        requested_set = set(requested)
        unknown = sorted(requested_set.difference(known))
        if unknown:
            raise ValueError(f"Unknown contrast IDs: {', '.join(unknown)}")
        selected = [contrast_id for contrast_id in known if contrast_id in requested_set]
        if not selected:
            raise ValueError("Select at least one contrast; an empty selection never means all.")
        if selected == known:
            normalized.pop("contrast_ids")
        else:
            normalized["contrast_ids"] = selected
        return normalized

    def _fingerprint(self, path: Path) -> FileFingerprint:
        resolved = path.resolve()
        stat = resolved.stat()
        key = (str(resolved), int(stat.st_size), int(stat.st_mtime_ns))
        with self._fingerprint_lock:
            cached = self._fingerprint_cache.get(key)
        if cached is not None:
            return cached
        observed = fingerprint_file(resolved)
        with self._fingerprint_lock:
            if len(self._fingerprint_cache) >= 4096:
                self._fingerprint_cache.clear()
            self._fingerprint_cache[key] = observed
        return observed

    def run_phase1_profile(
        self,
        project_id: str,
        *,
        parameters: dict[str, Any] | None = None,
        force: bool = False,
    ) -> JobRecord:
        """Compatibility endpoint retained for the Phase 1 acceptance contract."""
        with self._submission_lock:
            return self._run_phase1_profile(project_id, parameters=parameters, force=force)

    def _run_phase1_profile(self, project_id, *, parameters=None, force=False):
        context = self.projects.get(project_id)
        self.jobs.bind_project(project_id, context.manifest.output_root, execute=True)
        normalized = {"minimum_read_support": 20, **(parameters or {})}
        source = self.phase1.source_asset(context)
        run_fingerprint = workflow_fingerprint(
            workflow_id=self.phase1.workflow_id,
            workflow_version=self.phase1.version,
            project_fingerprint=context.snapshot.manifest_fingerprint,
            inputs=[fingerprint_file(source, force_full=True)],
            parameters=normalized,
        )
        cached = context.catalog.cached_artifacts(
            project_id, self.phase1.workflow_id, run_fingerprint
        )
        active = next((job for job in self.jobs.list(project_id)
                       if job.analysis_id == self.phase1.workflow_id and job.run_fingerprint == run_fingerprint
                       and job.status in {JobStatus.QUEUED, JobStatus.RUNNING}), None)
        if active is not None:
            return active
        if not force and len(cached) >= 2:
            return self.jobs.cache_hit(
                project_id=project_id,
                analysis_id=self.phase1.workflow_id,
                run_fingerprint=run_fingerprint,
                artifact_ids=[artifact.artifact_id for artifact in cached],
                parameters=normalized,
            )
        if force:
            context.catalog.clear_run_artifacts(
                project_id, self.phase1.workflow_id, run_fingerprint
            )
        return self.jobs.submit(
            project_id=project_id,
            analysis_id=self.phase1.workflow_id,
            run_fingerprint=run_fingerprint,
            runner=lambda job: self.phase1.run(
                context, job, run_fingerprint, normalized
            ),
            parameters=normalized,
        )
