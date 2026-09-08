from __future__ import annotations

from typing import Any

from katana_next.domain.models import JobRecord
from katana_next.infrastructure.fingerprints import (
    fingerprint_file,
    workflow_fingerprint,
)
from katana_next.infrastructure.jobs import JobManager
from katana_next.services.project_service import ProjectService
from katana_next.workflows.phase1_profile import Phase1ProfileWorkflow


class WorkflowService:
    def __init__(self, projects: ProjectService, jobs: JobManager) -> None:
        self.projects = projects
        self.jobs = jobs
        self.phase1 = Phase1ProfileWorkflow()

    def run_phase1_profile(
        self,
        project_id: str,
        *,
        parameters: dict[str, Any] | None = None,
        force: bool = False,
    ) -> JobRecord:
        context = self.projects.get(project_id)
        normalized_parameters = {"minimum_read_support": 20, **(parameters or {})}
        source = self.phase1.source_asset(context)
        source_fingerprint = fingerprint_file(source, force_full=True)
        run_fingerprint = workflow_fingerprint(
            workflow_id=self.phase1.workflow_id,
            workflow_version=self.phase1.version,
            project_fingerprint=context.snapshot.manifest_fingerprint,
            inputs=[source_fingerprint],
            parameters=normalized_parameters,
        )
        if not force:
            cached = context.catalog.cached_artifacts(
                project_id, self.phase1.workflow_id, run_fingerprint
            )
            if len(cached) >= 2:
                return self.jobs.cache_hit(
                    project_id=project_id,
                    analysis_id=self.phase1.workflow_id,
                    run_fingerprint=run_fingerprint,
                    artifact_ids=[artifact.artifact_id for artifact in cached],
                )

        return self.jobs.submit(
            project_id=project_id,
            analysis_id=self.phase1.workflow_id,
            run_fingerprint=run_fingerprint,
            runner=lambda job: self.phase1.run(
                context, job, run_fingerprint, normalized_parameters
            ),
        )
