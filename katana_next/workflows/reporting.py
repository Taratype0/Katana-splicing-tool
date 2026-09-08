from __future__ import annotations

import html
import json
import platform
import sys
import zipfile
from pathlib import Path
from typing import Any

import pandas as pd

from katana_next.domain.models import AnalysisArtifact, ArtifactKind
from katana_next.infrastructure.fingerprints import fingerprint_file
from katana_next.infrastructure.jobs import JobContext
from katana_next.services.project_service import ProjectContext
from katana_next.workflows.base import (
    ArtifactWriter,
    ParameterDefinition,
    WorkflowBase,
    WorkflowDefinition,
)


class ReproducibilityReportWorkflow(WorkflowBase):
    version = "2.0.1"
    definition = WorkflowDefinition(
        workflow_id="reports.bundle",
        workspace_id="reports",
        label="Reproducibility evidence package",
        question="What results, methods, parameters and provenance can be shared together?",
        method="Cache-only assembly of immutable artifacts; no upstream workflow is executed.",
        output_summary="HTML and PDF summaries, artifact inventory, machine-readable provenance and ZIP package.",
        analysis_unit="analysis artifact",
        expected_artifacts=5,
        parameters=(
            ParameterDefinition(
                "latest_only",
                "Latest result per output",
                "boolean",
                True,
                "Keep the newest compatible artifact for each workflow output name.",
            ),
            ParameterDefinition(
                "include_result_files",
                "Include result files in ZIP",
                "boolean",
                True,
                "Package tables, figures and reports; raw input files are never copied.",
            ),
        ),
        limitations="The package records observed results and claim boundaries; it does not convert exploratory evidence into causal proof.",
    )

    def parent_artifacts(
        self, context: ProjectContext, parameters: dict[str, Any]
    ) -> list[AnalysisArtifact]:
        artifacts = [
            item
            for item in context.catalog.list_artifacts(context.manifest.project_id)
            if item.analysis_id != self.definition.workflow_id
            and item.path.is_file()
            and item.metadata.get("project_snapshot_fingerprint")
            == context.snapshot.manifest_fingerprint
        ]
        if not bool(parameters.get("latest_only", True)):
            return artifacts
        selected: list[AnalysisArtifact] = []
        seen: set[tuple[str, str, str]] = set()
        for artifact in artifacts:
            key = (artifact.analysis_id, artifact.name, artifact.kind.value)
            if key in seen:
                continue
            seen.add(key)
            selected.append(artifact)
        return selected

    def input_paths(
        self, context: ProjectContext, parameters: dict[str, Any]
    ) -> list[Path]:
        return [artifact.path for artifact in self.parent_artifacts(context, parameters)]

    def run(
        self,
        context: ProjectContext,
        job: JobContext,
        run_fingerprint: str,
        parameters: dict[str, Any],
    ) -> list[str]:
        parents = self.parent_artifacts(context, parameters)
        if not parents:
            raise ValueError("No cached analysis artifacts are available to report.")
        writer = ArtifactWriter(
            context,
            self.definition.workflow_id,
            run_fingerprint,
            parameters,
            self.version,
            common_metadata={
                "question": self.definition.question,
                "method": self.definition.method,
                "analysis_unit": self.definition.analysis_unit,
                "interpretation": "Every result remains traceable to its own source artifact and run fingerprint.",
                "limitations": self.definition.limitations,
            },
        )
        job.report(0.12, "Indexing cached artifacts")
        inventory = _inventory(parents)
        inventory_artifact = writer.table(
            "Evidence package inventory",
            "evidence_inventory.parquet",
            inventory,
            entity_type="analysis_artifact",
            parents=tuple(item.artifact_id for item in parents),
        )
        provenance = {
            "schema_version": "2.0",
            "project": context.snapshot.model_dump(mode="json"),
            "manifest": context.manifest.model_dump(mode="json"),
            "python": sys.version,
            "platform": platform.platform(),
            "run_fingerprint": run_fingerprint,
            "parameters": parameters,
            "artifacts": [item.model_dump(mode="json") for item in parents],
        }
        provenance_path = writer.directory / "provenance.json"
        provenance_path.write_text(
            json.dumps(provenance, ensure_ascii=False, indent=2), encoding="utf-8"
        )
        provenance_artifact = writer.existing_file(
            "Machine-readable provenance",
            provenance_path,
            entity_type="analysis_provenance",
            metadata={"mime_type": "application/json"},
            parents=(inventory_artifact.artifact_id,),
        )
        job.report(0.35, "Writing HTML and PDF summaries")
        html_path = writer.directory / "KATANA_EVIDENCE_REPORT.html"
        html_path.write_text(_html_report(context, parents), encoding="utf-8")
        html_artifact = writer.existing_file(
            "Katana evidence report (HTML)",
            html_path,
            kind=ArtifactKind.REPORT,
            entity_type="project_report",
            metadata={"mime_type": "text/html"},
            parents=(inventory_artifact.artifact_id,),
        )
        pdf_path = writer.directory / "KATANA_EVIDENCE_REPORT.pdf"
        _pdf_report(context, parents, pdf_path)
        pdf_artifact = writer.existing_file(
            "Katana evidence report (PDF)",
            pdf_path,
            kind=ArtifactKind.REPORT,
            entity_type="project_report",
            metadata={"mime_type": "application/pdf"},
            parents=(inventory_artifact.artifact_id,),
        )
        job.report(0.68, "Packaging result files and provenance")
        archive_path = writer.directory / "KATANA_EVIDENCE_PACKAGE.zip"
        _write_archive(
            archive_path,
            context,
            parents,
            html_path,
            pdf_path,
            provenance_path,
            include_files=bool(parameters.get("include_result_files", True)),
        )
        writer.existing_file(
            "Katana reproducibility package",
            archive_path,
            entity_type="reproducibility_bundle",
            metadata={"mime_type": "application/zip"},
            parents=(
                inventory_artifact.artifact_id,
                provenance_artifact.artifact_id,
                html_artifact.artifact_id,
                pdf_artifact.artifact_id,
            ),
        )
        job.report(0.96, "Evidence package complete")
        return writer.finish()


def _inventory(artifacts: list[AnalysisArtifact]) -> pd.DataFrame:
    rows = []
    for item in artifacts:
        file_fingerprint = fingerprint_file(item.path, force_full=True)
        rows.append(
            {
                "row_id": item.artifact_id,
                "workflow_id": item.analysis_id,
                "artifact_name": item.name,
                "kind": item.kind.value,
                "entity_type": item.entity_type,
                "rows": item.row_count,
                "path": str(item.path),
                "sha256": file_fingerprint.digest,
                "run_fingerprint": item.run_fingerprint,
                "created_at": item.created_at.isoformat(),
            }
        )
    return pd.DataFrame(rows)


def _html_report(context: ProjectContext, artifacts: list[AnalysisArtifact]) -> str:
    groups: dict[str, list[AnalysisArtifact]] = {}
    for artifact in artifacts:
        groups.setdefault(artifact.analysis_id, []).append(artifact)
    sections = []
    for analysis_id, values in groups.items():
        metadata = values[0].metadata
        rows = "".join(
            "<tr>"
            f"<td>{html.escape(value.name)}</td>"
            f"<td>{html.escape(value.kind.value)}</td>"
            f"<td>{value.row_count if value.row_count is not None else ''}</td>"
            f"<td><code>{html.escape(value.run_fingerprint[:12])}</code></td>"
            "</tr>"
            for value in values
        )
        sections.append(
            f"<section><p class='kicker'>{html.escape(analysis_id)}</p>"
            f"<h2>{html.escape(str(metadata.get('question', analysis_id)))}</h2>"
            f"<p><strong>Method.</strong> {html.escape(str(metadata.get('method', 'Recorded in artifact metadata.')))}</p>"
            f"<p><strong>Interpretation.</strong> {html.escape(str(metadata.get('interpretation', 'Review source artifacts.')))}</p>"
            f"<p class='boundary'><strong>Boundary.</strong> {html.escape(str(metadata.get('limitations', 'No additional causal claim.')))}</p>"
            "<table><thead><tr><th>Artifact</th><th>Kind</th><th>Rows</th><th>Run</th></tr></thead>"
            f"<tbody>{rows}</tbody></table></section>"
        )
    return f"""<!doctype html><html><head><meta charset='utf-8'><title>Katana evidence report</title>
<style>body{{font-family:Arial,sans-serif;margin:0;color:#17252a;background:#f4f4ef}}main{{max-width:1080px;margin:auto;padding:52px}}header{{padding:34px;background:#12383b;color:white}}h1{{font-size:40px;margin:0 0 8px}}h2{{font-size:23px}}section{{margin:22px 0;padding:26px;background:white;border:1px solid #d9ded8}}.kicker{{color:#087d70;text-transform:uppercase;letter-spacing:.1em;font-size:12px;font-weight:700}}.boundary{{padding:12px;border-left:3px solid #da6f32;background:#fff4e8}}table{{width:100%;border-collapse:collapse}}th,td{{padding:9px;text-align:left;border-bottom:1px solid #e3e6e1;font-size:13px}}code{{font-size:12px}}</style></head>
<body><header><h1>{html.escape(context.manifest.name)}</h1><p>Katana Next reproducibility evidence report</p></header><main>
<p><strong>Genome:</strong> {html.escape(context.manifest.annotation.genome_build)} &nbsp; <strong>Samples:</strong> {len(context.manifest.samples)} &nbsp; <strong>Contrasts:</strong> {len(context.manifest.contrasts)}</p>
{''.join(sections)}</main></body></html>"""


def _pdf_report(
    context: ProjectContext, artifacts: list[AnalysisArtifact], target: Path
) -> None:
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import getSampleStyleSheet
    from reportlab.lib.units import mm
    from reportlab.platypus import Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle
    from reportlab.lib import colors

    styles = getSampleStyleSheet()
    story = [
        Paragraph(context.manifest.name, styles["Title"]),
        Paragraph("Katana Next reproducibility evidence report", styles["Heading2"]),
        Spacer(1, 5 * mm),
        Paragraph(
            f"Genome {context.manifest.annotation.genome_build}; "
            f"{len(context.manifest.samples)} samples; {len(context.manifest.contrasts)} contrasts.",
            styles["BodyText"],
        ),
    ]
    groups: dict[str, list[AnalysisArtifact]] = {}
    for artifact in artifacts:
        groups.setdefault(artifact.analysis_id, []).append(artifact)
    for analysis_id, values in groups.items():
        story.extend([Spacer(1, 6 * mm), Paragraph(analysis_id, styles["Heading2"])])
        metadata = values[0].metadata
        for label, key in (("Question", "question"), ("Method", "method"), ("Interpretation", "interpretation"), ("Boundary", "limitations")):
            if metadata.get(key):
                story.append(Paragraph(f"<b>{label}.</b> {html.escape(str(metadata[key]))}", styles["BodyText"]))
        data = [["Artifact", "Kind", "Rows"]] + [
            [item.name[:72], item.kind.value, str(item.row_count or "")]
            for item in values
        ]
        table = Table(data, colWidths=[120 * mm, 30 * mm, 20 * mm], repeatRows=1)
        table.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#12383b")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#cbd5d1")),
            ("FONTSIZE", (0, 0), (-1, -1), 8),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ]))
        story.append(table)
    document = SimpleDocTemplate(str(target), pagesize=A4, rightMargin=18 * mm, leftMargin=18 * mm, topMargin=18 * mm, bottomMargin=18 * mm)
    document.build(story)


def _write_archive(
    target: Path,
    context: ProjectContext,
    artifacts: list[AnalysisArtifact],
    html_path: Path,
    pdf_path: Path,
    provenance_path: Path,
    *,
    include_files: bool,
) -> None:
    with zipfile.ZipFile(target, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.write(html_path, "report/KATANA_EVIDENCE_REPORT.html")
        archive.write(pdf_path, "report/KATANA_EVIDENCE_REPORT.pdf")
        archive.write(provenance_path, "provenance/provenance.json")
        archive.write(context.manifest_path, "provenance/katana.project.json")
        if include_files:
            seen: set[Path] = set()
            for artifact in artifacts:
                source = artifact.path.resolve()
                if source in seen or not source.is_file():
                    continue
                seen.add(source)
                safe_analysis = artifact.analysis_id.replace("/", "_").replace("\\", "_")
                archive.write(
                    source,
                    f"artifacts/{safe_analysis}/{artifact.artifact_id[:8]}_{source.name}",
                )
