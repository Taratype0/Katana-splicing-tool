"""In-process Python interface sharing the web workflow and cache services."""
from pathlib import Path

from katana_next.infrastructure.jobs import JobManager
from katana_next.services.project_service import ProjectService
from katana_next.services.workflow_service import WorkflowService
from katana_next.services.readiness import experimental_design, preflight


class KatanaProject:
    def __init__(self, manifest: str | Path):
        self.projects = ProjectService()
        self.context = self.projects.open_project(manifest)
        self.project_id = self.context.manifest.project_id
        self.jobs = JobManager(max_workers=1)
        self.jobs.bind_project(self.project_id, self.context.manifest.output_root)
        self.workflows = WorkflowService(self.projects, self.jobs)

    def inventory(self):
        return self.workflows.definitions()

    def design(self):
        return experimental_design(self.context.manifest)

    def editable_design(self):
        from katana_next.services.design_editor import editable_design
        return editable_design(self.context)

    def revise_design(self, edit, *, apply=False):
        from katana_next.services.design_editor import revise_design
        with self.workflows._submission_lock:
            if any(j.status in {"queued", "running"} for j in self.jobs.list(self.project_id)):
                raise ValueError("Finish active jobs before editing design")
            result = revise_design(self.projects, self.project_id, edit, apply=apply)
            self.context = self.projects.get(self.project_id)
        return result

    def plan(self, workflow: str, parameters=None):
        return preflight(self.workflows, self.project_id, workflow, parameters)

    def contrasts(self, factors=None, *, include_batch=False, pairing="auto", reference_levels=None):
        from katana_next.science.contrast_planner import plan_contrasts
        return plan_contrasts(self.context.manifest, factors, include_batch=include_batch, pairing=pairing, reference_levels=reference_levels)

    def retry(self, job_id, *, timeout=None):
        job = self.workflows.retry(job_id)
        result = self.jobs.wait(job.job_id, timeout=timeout)
        if result.status not in {"complete", "cache_hit"}:
            raise RuntimeError(result.error or result.message)
        return result

    def run(self, workflow: str, parameters=None, *, force=False, timeout=None):
        job = self.workflows.run(workflow, self.project_id, parameters=parameters, force=force)
        try:
            result = self.jobs.wait(job.job_id, timeout=timeout)
        except (TimeoutError, KeyboardInterrupt):
            self.jobs.cancel(job.job_id)
            raise
        if result.status not in {"complete", "cache_hit"}:
            raise RuntimeError(result.error or result.message)
        return result

    def artifacts(self, workflow=None):
        return self.projects.list_artifacts(self.project_id, analysis_id=workflow)

    def module(self, source_artifact_id: str, module: str, *, limit=30):
        from katana_next.services.module_browser import module_evidence
        return module_evidence(self.context, source_artifact_id, module, limit=limit)

    def close(self):
        self.jobs.shutdown()

    def module_dashboard(self, source_artifact_id: str, module: str, **options):
        from katana_next.services.module_dashboard import module_dashboard
        return module_dashboard(self.context, source_artifact_id, module, **options)

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()
