"""Register the approved seven-section report without executing analyses."""
from __future__ import annotations

import json
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

from katana_next.domain.models import ArtifactKind
from katana_next.workflows.base import ArtifactWriter, ParameterDefinition, WorkflowBase, WorkflowDefinition


class StudyReportWorkflow(WorkflowBase):
    version = "1.0.0"
    definition = WorkflowDefinition(
        workflow_id="reports.study", workspace_id="reports",
        label="Seven-section analysis report",
        question="What do the completed design, quality, differential, network and MotifTrace results show?",
        method="Read-only assembly of completed project results, preserving native effects, P/FDR, comparison identities and missing states. No analysis or upstream execution.",
        output_summary="Portable seven-section HTML report, complete-result browsers, figures, provenance, machine-readable manifest and ZIP package.",
        analysis_unit="completed analysis artifact and declared comparison", expected_artifacts=4,
        parameters=(
            ParameterDefinition("event_set", "DAS event set", "select", "stable_significant",
                                "RS-DAS is the default. Switchable completed downstream results must match; RS-required and measured-PSI modules retain their own scope.",
                                ("stable_significant", "significant")),
        ),
        limitations="Unavailable analyses remain not computed. The report does not fit new models, rerun upstream tasks, pool comparisons or assign biological causality.",
    )

    def parent_artifacts(self, context, parameters):
        del parameters
        selected, seen = [], set()
        for artifact in context.catalog.list_artifacts(context.manifest.project_id):
            if (artifact.analysis_id.startswith("reports.") or str(artifact.status) != "complete"
                    or artifact.project_id != context.manifest.project_id or not artifact.path.is_file()
                    or artifact.metadata.get("project_snapshot_fingerprint") != context.snapshot.manifest_fingerprint):
                continue
            if artifact.metadata.get("artifact_contract_version") == "3.0-correctness":
                marker = (Path(context.manifest.output_root) / "artifacts" / artifact.analysis_id
                          / artifact.run_fingerprint[:16] / "parameters.json")
                try:
                    complete = json.loads(marker.read_text(encoding="utf-8-sig"))
                except (ValueError, OSError):
                    continue
                if (complete.get("completed") is not True or complete.get("workflow_id") != artifact.analysis_id
                        or complete.get("run_fingerprint") != artifact.run_fingerprint
                        or artifact.artifact_id not in {row.get("artifact_id") for row in complete.get("generated_artifacts", [])}):
                    continue
            scope = artifact.metadata.get("das_selection", {}).get("event_set", artifact.parameters.get("event_set"))
            key = (artifact.analysis_id, artifact.name, str(artifact.kind), tuple(sorted(artifact.contrast_ids)),
                   artifact.parameters.get("sequence_scope"), scope)
            if key not in seen:
                seen.add(key)
                selected.append(artifact)
        return selected

    def input_paths(self, context, parameters):
        from katana_next.reports import engine, sources, rich_results
        paths = [Path(context.manifest_path), Path(__file__)]
        paths.extend(path for path in Path(engine.__file__).parent.rglob("*.py")
                     if "__pycache__" not in path.parts)
        paths.extend(path for path in engine.ASSET_ROOT.rglob("*") if path.is_file())
        for artifact in self.parent_artifacts(context, parameters):
            paths.append(artifact.path)
            paths.extend(Path(path) for path in artifact.metadata.get("exports", {}).values() if Path(path).is_file())
        return list(dict.fromkeys(path.resolve() for path in paths if path.is_file()))

    def run(self, context, job, run_fingerprint, parameters):
        from katana_next.reports import build_report
        parents = self.parent_artifacts(context, parameters)
        if not parents:
            raise ValueError("No completed results from this project snapshot are available. Import or finish the required analysis explicitly.")
        job.check_cancelled()
        job.report(0.05, "Reading completed analysis results; no upstream tasks are started")
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint,
                                parameters, self.version, common_metadata={
                                    "report_contract": "seven-section-v1", "no_upstream_rerun": True,
                                    "no_source_statistics_refit": True,
                                })
        destination = writer.directory / "study-report"
        if destination.exists():
            raise FileExistsError("Existing study report is preserved; create a new report run.")
        descriptor = {"project": context.manifest.model_dump(mode="json"),
                      "artifacts": [artifact.model_dump(mode="json") for artifact in parents]}
        # The public builder always requires output outside original results.
        # Render there, then register the finished portable package as a new
        # workflow output. Existing parent result directories are never edited.
        with TemporaryDirectory(prefix="katana-seven-section-report-") as temporary:
            rendered = Path(temporary) / "report"
            build_report(descriptor, rendered, event_set=parameters.get("event_set", "stable_significant"))
            job.check_cancelled()
            job.report(0.8, "Registering portable report, full results and provenance")
            shutil.copytree(rendered, destination)
        parent_ids = tuple(artifact.artifact_id for artifact in parents)
        exports = {path.relative_to(destination).as_posix(): str(path)
                   for path in destination.rglob("*") if path.is_file()}
        writer.existing_file("Seven-section Katana report", destination / "REPORT.html",
                             kind=ArtifactKind.REPORT, entity_type="project_report", parents=parent_ids,
                             metadata={"mime_type": "text/html", "preserve_relative_report_links": True,
                                       "exports": exports})
        writer.existing_file("Portable report and complete results", destination / "REPORT.zip",
                             kind=ArtifactKind.REPORT, entity_type="report_package", parents=parent_ids,
                             metadata={"mime_type": "application/zip"})
        writer.existing_file("Seven-section report provenance", destination / "provenance.json",
                             entity_type="report_provenance", parents=parent_ids,
                             metadata={"mime_type": "application/json"})
        writer.existing_file("Seven-section report manifest", destination / "report-manifest.json",
                             entity_type="report_manifest", parents=parent_ids,
                             metadata={"mime_type": "application/json"})
        job.report(1.0, "Seven-section report and complete portable outputs registered")
        return writer.finish()
