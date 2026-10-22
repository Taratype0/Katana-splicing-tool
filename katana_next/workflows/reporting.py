from __future__ import annotations

import html
import base64
import json
import platform
import sys
import zipfile
from pathlib import Path
from typing import Any

import pandas as pd
from katana_next.science.report_language import is_retired_report
from katana_next.science.network_label_report import NETWORK_LABEL_METHODS

from katana_cosplice.branding import display_module_name

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
    version = "2.3.0+das_selection.20261022"
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
            ParameterDefinition("event_set", "DAS event set", "select", "stable_significant",
                "RS-DAS is the default. Include only matching completed switchable downstream results; absent results are not recomputed or relabelled. RS-required and measured-PSI modules keep their own declared scope.",
                ("stable_significant", "significant")),
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
        from katana_next.workflows.sources import source_complete
        artifacts = [
            item
            for item in context.catalog.list_artifacts(context.manifest.project_id)
            if item.analysis_id != self.definition.workflow_id
            and item.path.is_file()
            and not is_retired_report(item.path)
            and source_complete(context, item)
            and item.metadata.get("project_snapshot_fingerprint")
            == context.snapshot.manifest_fingerprint
        ]
        artifacts = [item for item in artifacts if _report_artifact_matches(item, parameters.get("event_set"))]
        if not bool(parameters.get("latest_only", True)):
            return artifacts
        selected: list[AnalysisArtifact] = []
        seen: set[tuple] = set()
        for artifact in artifacts:
            from katana_next.science.rbp_evidence_layers import artifact_layer
            scope = artifact.parameters.get("sequence_scope") if artifact_layer(artifact) else None
            key = (artifact.analysis_id, artifact.name, artifact.kind.value, tuple(sorted(artifact.contrast_ids)), scope)
            if key in seen:
                continue
            seen.add(key)
            selected.append(artifact)
        return selected

    def input_paths(
        self, context: ProjectContext, parameters: dict[str, Any]
    ) -> list[Path]:
        from katana_next.science.rbp_evidence_layers import CONTRACT_PATH
        return [artifact.path for artifact in self.parent_artifacts(context, parameters)] + [
            Path(__file__), CONTRACT_PATH, CONTRACT_PATH.with_suffix('.py')]

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
        if context.manifest.metadata.get("gene_filter"):
            from katana_next.science.gene_filter import filter_adapter
            from src.analysis.gene_filter import project_gene_filter
            annotation_audit = project_gene_filter(filter_adapter(context.manifest)).annotation_audit()
            annotation_audit["row_id"] = annotation_audit.gene_id
            writer.table("Gene exclusion annotation decisions", "gene_filter_decisions.parquet",
                         annotation_audit, entity_type="gene_filter_audit")
            counts_path = context.manifest.output_root / "gene_filter" / "input_filter_counts.tsv"
            if counts_path.is_file():
                counts = pd.read_csv(counts_path, sep="\t")
                counts["row_id"] = counts.asset_id
                writer.table("Gene exclusion input counts", "gene_filter_input_counts.parquet",
                             counts, entity_type="gene_filter_audit")
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
        from katana_next.science.das_selection import selection_metadata
        das_policy = selection_metadata(parameters.get("event_set"))
        das_policy["scope"] = "Matching saved downstream artifacts only; no analysis or caller refitting. RS-required and measured-PSI modules retain their separate declared populations."
        provenance["das_report_selection"] = das_policy
        from katana_next.science.rbp_evidence_layers import CONTRACT, layer_inventory
        provenance["rbp_evidence_layers"] = {
            "contract_version": CONTRACT["version"],
            "inventory": layer_inventory(parents),
            "display_only_no_combined_P": True,
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
        from katana_next.services.module_dashboard import module_dashboard
        from katana_next.science.module_dashboard_report import render_dashboard_png, dashboard_section
        import re
        payloads = []
        for artifact in parents:
            if artifact.analysis_id not in {"programs.cosplice", "programs.psi_network"} or Path(artifact.path).stem != "consensus_modules__assignments":
                continue
            membership = pd.read_parquet(artifact.path)
            label = "consensus_module" if "consensus_module" in membership else "module"
            for module in sorted(set(membership[label].dropna().astype(str)) - {"grey"}):
                for condition in ["all", *sorted({s.condition for s in context.manifest.samples})]:
                    d = module_dashboard(context, artifact.artifact_id, module, condition=condition, include_measurements=True)
                    payloads.append(d)
                    if condition != "all":
                        continue
                    image = writer.directory / f"module_{artifact.run_fingerprint[:12]}_{re.sub(r'[^a-zA-Z0-9_-]', '_', module)}.png"
                    render_dashboard_png(d, image)
                    plot = writer.existing_file(f"Module {module}: event network, scaled PSI and condition activity",image,
                        entity_type="module_dashboard", metadata={"mime_type":"image/png","interpretation":d["interpretation"]},parents=(artifact.artifact_id,))
                    parents.append(plot)
        html_path = writer.directory / "KATANA_EVIDENCE_REPORT.html"
        content = _html_report(context, parents, das_policy=das_policy)
        if payloads:
            content = content.replace("</main>", dashboard_section(payloads)+"</main>")
            (writer.directory / "module_dashboards.json").write_text(json.dumps(payloads,allow_nan=False),encoding="utf-8")
            payload_artifact = writer.existing_file("Event module dashboard payloads", writer.directory / "module_dashboards.json", entity_type="module_dashboard_payload",
                                 metadata={"mime_type":"application/json"},parents=tuple(a.artifact_id for a in parents if a.analysis_id in {"programs.cosplice","programs.psi_network"}))
            parents.append(payload_artifact)
        html_path.write_text(content, encoding="utf-8")
        html_artifact = writer.existing_file(
            "Katana evidence report (HTML)",
            html_path,
            kind=ArtifactKind.REPORT,
            entity_type="project_report",
            metadata={"mime_type": "text/html"},
            parents=(inventory_artifact.artifact_id,),
        )
        pdf_path = writer.directory / "KATANA_EVIDENCE_REPORT.pdf"
        _pdf_report(context, parents, pdf_path, das_policy=das_policy)
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


def _report_artifact_matches(artifact: AnalysisArtifact, event_set: Any = None) -> bool:
    """Reporting never gives a saved plot/table another selection identity.

    Even a full DAS table may contain a recorded selected_das column. Requiring
    the stored context here prevents an ordinary-DAS report from quietly showing
    that RS-selected column under a new heading. Source files remain untouched.
    """
    from katana_next.science.das_selection import (
        DAS_SELECTION_POLICY_VERSION, normalize_event_set, workflow_das_capability,
    )
    if workflow_das_capability(artifact.analysis_id)["mode"] != "switchable":
        return True
    recorded = artifact.metadata.get("das_selection")
    return (isinstance(recorded, dict)
            and recorded.get("policy_version") == DAS_SELECTION_POLICY_VERSION
            and recorded.get("event_set") == normalize_event_set(event_set))


def _artifact_scope_label(artifact: AnalysisArtifact) -> str:
    from katana_next.science.das_selection import workflow_das_capability
    capability = workflow_das_capability(artifact.analysis_id)
    recorded = artifact.metadata.get("das_selection", {})
    saved = recorded.get("label") if isinstance(recorded, dict) else None
    if capability["mode"] == "rs_required":
        return "RS-DAS required for new analysis; saved scope: " + (str(saved) if saved else "not recorded")
    if capability["mode"] == "measurement_population":
        return "Declared measured PSI / read-support population; DAS switch does not apply"
    return str(saved) if saved else "Own recorded result scope; no DAS relabelling"


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
                "das_scope": _artifact_scope_label(item),
            }
        )
    return pd.DataFrame(rows)


def _html_report(context: ProjectContext, artifacts: list[AnalysisArtifact], *, das_policy=None) -> str:
    from katana_next.science.gene_filter import gene_filter_footnote
    from katana_next.science.rbp_evidence_layers import report_groups, overview_html, layer_note_html
    filter_note = html.escape(gene_filter_footnote(context.manifest))
    sections = []
    for analysis_id, layer, values in report_groups(artifacts):
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
            + layer_note_html(layer)
            +
            f"<h2>{html.escape(str(metadata.get('question', analysis_id)))}</h2>"
            + f"<p class='das-scope'><strong>DAS scope.</strong> {html.escape(_artifact_scope_label(values[0]))}</p>"
            + (_das_top10_html(values, das_policy["event_set"]) if analysis_id == "differential.as" and das_policy else "")
            + "".join(_embedded_figure(value) for value in values)
            + "".join(_summary_table(value) for value in values)
            + _network_module_parameters_html(analysis_id, values)
            + _network_label_policy_html(analysis_id, values)
            + "<details><summary>Methods, interpretation and source inventory</summary>"
            f"<p><strong>Method.</strong> {html.escape(str(metadata.get('method', 'Recorded in artifact metadata.')))}</p>"
            f"<p><strong>Interpretation.</strong> {html.escape(str(metadata.get('interpretation', 'Review source artifacts.')))}</p>"
            f"<p class='boundary'><strong>Boundary.</strong> {html.escape(str(metadata.get('limitations', 'No additional causal claim.')))}</p>"
            "<table><thead><tr><th>Artifact</th><th>Kind</th><th>Rows</th><th>Run</th></tr></thead>"
            f"<tbody>{rows}</tbody></table></details></section>"
        )
    contrasts = "".join(f"<tr><td>{html.escape(c.contrast_id)}</td><td>{html.escape(c.label)}</td><td>{html.escape(c.numerator.label)}</td><td>{html.escape(c.denominator.label)}</td></tr>" for c in context.manifest.contrasts)
    selected_scope = "" if das_policy is None else (
        "<section id='das-report-selection'><h2>DAS report selection</h2>"
        f"<p><strong>{html.escape(str(das_policy['label']))}</strong>. {html.escape(str(das_policy['scope']))}</p>"
        "<p>A switchable module without a matching completed result is not computed for this view. "
        "No result from another DAS selection is substituted. Full original results remain in their source output directories.</p></section>"
    )
    return f"""<!doctype html><html><head><meta charset='utf-8'><title>Katana evidence report</title>
<style>body{{font-family:Arial,sans-serif;margin:0;color:#17252a;background:#f4f4ef}}main{{max-width:1080px;margin:auto;padding:52px}}header{{padding:34px;background:#12383b;color:white}}h1{{font-size:40px;margin:0 0 8px}}h2{{font-size:23px}}section{{margin:22px 0;padding:26px;background:white;border:1px solid #d9ded8}}.kicker{{color:#087d70;text-transform:uppercase;letter-spacing:.1em;font-size:12px;font-weight:700}}.boundary{{padding:12px;border-left:3px solid #da6f32;background:#fff4e8}}table{{width:100%;border-collapse:collapse}}th,td{{padding:9px;text-align:left;border-bottom:1px solid #e3e6e1;font-size:13px}}code{{font-size:12px}}</style></head>
<body><header><h1>{html.escape(context.manifest.name)}</h1><p>Katana Next reproducibility evidence report</p></header><main>
<p><strong>Genome:</strong> {html.escape(context.manifest.annotation.genome_build)} &nbsp; <strong>Samples:</strong> {len(context.manifest.samples)} &nbsp; <strong>Contrasts:</strong> {len(context.manifest.contrasts)}</p>
{selected_scope}
<section><h2>Comparisons and effect direction</h2><p>All standardized effects use numerator minus denominator. These are distinct experimental contrasts, not independent replications of one effect.</p><table><tr><th>ID</th><th>Comparison</th><th>Numerator</th><th>Denominator</th></tr>{contrasts}</table></section>
{overview_html(artifacts)}
{''.join(sections)}<footer id="gene-exclusion-footnote"><h2>Gene exclusion</h2><p>{filter_note}</p></footer></main></body></html>"""


def _png_path(artifact: AnalysisArtifact) -> Path | None:
    name = artifact.metadata.get("exports", {}).get("png")
    path = Path(name) if name else artifact.path if artifact.path.suffix.lower() == ".png" else None
    return path if path and path.is_file() else None


def _das_top10_frames(artifacts: list[AnalysisArtifact], event_set: str):
    """Read-only per-comparison ranking from the saved full DAS source."""
    from katana_next.science.das_selection import selected_event_mask
    ranked = []
    for artifact in artifacts:
        if artifact.name != "Differential AS events" or artifact.kind != ArtifactKind.TABLE:
            continue
        frame = pd.read_parquet(artifact.path) if artifact.path.suffix == ".parquet" else pd.read_csv(artifact.path, sep="\t")
        required = {"comparison_id", "event_id", "dpsi", "fdr"}
        if not required.issubset(frame.columns):
            continue
        chosen = frame.loc[selected_event_mask(frame, event_set)].copy()
        chosen["_rank_fdr"] = pd.to_numeric(chosen["fdr"], errors="coerce")
        chosen["_rank_effect"] = pd.to_numeric(chosen["dpsi"], errors="coerce").abs()
        chosen = chosen.sort_values(["_rank_fdr", "_rank_effect", "event_id"], ascending=[True, False, True], na_position="last", kind="stable")
        for comparison in frame["comparison_id"].dropna().astype(str).drop_duplicates():
            top = chosen.loc[chosen["comparison_id"].astype(str).eq(comparison)].head(10).drop(columns=["_rank_fdr", "_rank_effect"])
            ranked.append((comparison, top, artifact))
    return ranked


def _das_top10_html(artifacts: list[AnalysisArtifact], event_set: str) -> str:
    sections = []
    for comparison, top, artifact in _das_top10_frames(artifacts, event_set):
        label = "RS-DAS" if event_set == "stable_significant" else "DAS"
        columns = [key for key in ("gene_symbol", "gene_id", "event_id", "event_type", "dpsi", "fdr") if key in top]
        detail = "".join(
            f"<details><summary>{html.escape(str(row.get('gene_symbol', row['event_id'])))}: event fields</summary>"
            + pd.DataFrame({"Field": list(row.index), "Value": [str(value) for value in row]}).to_html(index=False, escape=True, border=0)
            + "</details>" for _, row in top.iterrows()
        )
        sections.append(
            f"<div class='das-top10'><h3>Top {len(top)} {label} events | {html.escape(comparison)}</h3>"
            "<p>Ranked within this comparison by native FDR, absolute PSI effect and event ID. "
            f"The saved full table remains unchanged (source run {html.escape(artifact.run_fingerprint[:12])}).</p>"
            + (top[columns].to_html(index=False, escape=True, border=0) if len(top) else "<p>No selected events; this is a computed empty foreground.</p>")
            + detail + "</div>"
        )
    return "".join(sections)


def _network_module_parameter_records(
    analysis_id: str, artifacts: list[AnalysisArtifact]
) -> list[dict[str, str]]:
    """Read actual saved values only; requested controls do not prove execution."""
    network_ids = {"programs.cosplice", "programs.psi_network", "single_cell.programs"}
    if analysis_id not in network_ids and not any(
        "network_module_parameters" in item.metadata for item in artifacts
    ):
        return []
    rows = []
    seen = set()
    for item in artifacts:
        recorded = item.metadata.get("network_module_parameters")
        recorded = recorded if isinstance(recorded, dict) else {}
        minimum = recorded.get("min_module_size")
        merge = recorded.get("merge_cut_height")
        row = {
            "run": item.run_fingerprint,
            "min_module_size": str(minimum) if minimum is not None else "not recorded",
            "merge_cut_height": str(merge) if merge is not None else "not recorded",
        }
        key = tuple(row.values())
        if key not in seen:
            seen.add(key)
            rows.append(row)
    return rows


def _network_module_parameters_html(
    analysis_id: str, artifacts: list[AnalysisArtifact]
) -> str:
    rows = _network_module_parameter_records(analysis_id, artifacts)
    if not rows:
        return ""
    body = "".join(
        f"<tr><td><code>{html.escape(row['run'])}</code></td>"
        f"<td>{html.escape(row['min_module_size'])}</td>"
        f"<td>{html.escape(row['merge_cut_height'])}</td></tr>"
        for row in rows
    )
    return (
        "<div class='network-module-parameters'><h3>Recorded actual module parameters</h3>"
        "<p>Values are taken from the saved network run. Missing historical values are not inferred from current defaults or requested controls.</p>"
        "<table><thead><tr><th>Source run</th><th>Minimum module size (events)</th>"
        "<th>Module merge cut height</th></tr></thead>"
        f"<tbody>{body}</tbody></table></div>"
    )


def _network_label_policy_records(analysis_id, artifacts):
    network_ids = {"programs.cosplice", "programs.psi_network", "programs.hub_events", "single_cell.programs"}
    if analysis_id not in network_ids and not any("network_label_policy" in item.metadata for item in artifacts):
        return []
    rows, seen = [], set()
    for item in artifacts:
        saved = item.metadata.get("network_label_policy")
        policy = saved if isinstance(saved, dict) else {}
        enabled = policy.get("apply_threshold_labels")
        filtering = policy.get("hard_filter")
        row = {"run": item.run_fingerprint,
               "policy_id": str(policy.get("policy_id") or "not recorded"),
               "apply_threshold_labels": "on" if enabled is True else "off" if enabled is False else "not recorded",
               "hard_filter": "yes" if filtering is True else "no" if filtering is False else "not recorded",
               "label_thresholds": json.dumps(policy["label_thresholds"], sort_keys=True)
               if isinstance(policy.get("label_thresholds"), dict) else "not recorded"}
        key = tuple(row.values())
        if key not in seen:
            rows.append(row)
            seen.add(key)
    return rows


def _network_label_policy_html(analysis_id, artifacts):
    rows = _network_label_policy_records(analysis_id, artifacts)
    if not rows:
        return ""
    body = "".join(
        "<tr>" + "".join(f"<td>{html.escape(row[key])}</td>" for key in
                         ("run", "policy_id", "apply_threshold_labels", "hard_filter", "label_thresholds")) + "</tr>"
        for row in rows
    )
    return ("<div class='network-label-policy'><h3>Recorded optional network labels</h3>"
            "<p>Off means labels were not applied; it does not mean that events failed a label rule. "
            "Missing historical policy records are not replaced by current defaults.</p>"
            "<table><thead><tr><th>Source run</th><th>Policy</th><th>Threshold labels</th>"
            f"<th>Hard filtering</th><th>Recorded label rules</th></tr></thead><tbody>{body}</tbody></table>"
            f"<p>{html.escape(NETWORK_LABEL_METHODS)}</p></div>")


def _embedded_figure(artifact: AnalysisArtifact) -> str:
    path = _png_path(artifact)
    if path is None:
        return ""
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return (f"<figure><img style='width:100%;height:auto' alt='{html.escape(display_module_name(artifact.name), quote=True)}' src='data:image/png;base64,{encoded}'/>"
            f"<figcaption><strong>{html.escape(display_module_name(artifact.name))}</strong>. Source run: {artifact.run_fingerprint[:12]}. "
            f"{html.escape(str(artifact.metadata.get('interpretation', '')))}</figcaption></figure>")


def _summary_table(artifact: AnalysisArtifact) -> str:
    if artifact.kind != ArtifactKind.TABLE or not artifact.row_count or artifact.row_count > 100:
        return ""
    frame = pd.read_parquet(artifact.path) if artifact.path.suffix == ".parquet" else pd.read_csv(artifact.path, sep="\t")
    excluded = {"row_id", "path", "source_path"}
    columns = [c for c in frame if c not in excluded][:10]
    return (f"<h3>{html.escape(display_module_name(artifact.name))}</h3><p>First {min(12, len(frame))} of {len(frame)} rows; complete TSV is included.</p>"
            + frame[columns].head(12).to_html(index=False, escape=True, float_format=lambda x: f"{x:.4g}", border=0))


def _pdf_report(
    context: ProjectContext, artifacts: list[AnalysisArtifact], target: Path, *, das_policy=None
) -> None:
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import getSampleStyleSheet
    from reportlab.lib.units import mm
    from reportlab.platypus import Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle, Image, PageBreak
    from PIL import Image as PillowImage
    from reportlab.lib import colors

    styles = getSampleStyleSheet()
    from katana_next.science.gene_filter import gene_filter_footnote
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
    from katana_next.science.rbp_evidence_layers import CONTRACT, layer_inventory, report_groups
    if das_policy is not None:
        story.append(Paragraph("DAS report selection: " + html.escape(str(das_policy["label"])), styles["Heading2"]))
        story.append(Paragraph(html.escape(str(das_policy["scope"])), styles["BodyText"]))
        story.append(Paragraph("Missing matching downstream results are not computed for this view; no other selection is substituted.", styles["BodyText"]))
    story.append(Paragraph(html.escape(gene_filter_footnote(context.manifest)), styles["BodyText"]))
    story.append(Paragraph(CONTRACT['title_en'], styles['Heading2']))
    for layer in layer_inventory(artifacts):
        story.append(Paragraph(f"{layer['number']}. {html.escape(layer['title_en'])}", styles['Heading3']))
        for key in ('value_en', 'unit_en', 'boundary_en'):
            story.append(Paragraph(html.escape(layer[key]), styles['BodyText']))
        status = (f"{len(layer['artifact_ids'])} saved outputs; not gene/event counts." if layer['artifact_ids'] else
                  "No saved output of this layer is included; not a negative result and no analysis was started.")
        story.append(Paragraph(status, styles['BodyText']))
    story.append(Paragraph(CONTRACT['shared_boundary_en'], styles['BodyText']))
    for analysis_id, layer, values in report_groups(artifacts):
        story.extend([Spacer(1, 6 * mm), Paragraph(analysis_id, styles["Heading2"])])
        story.append(Paragraph("DAS scope: " + html.escape(_artifact_scope_label(values[0])), styles["BodyText"]))
        if analysis_id == "differential.as" and das_policy:
            label = "RS-DAS" if das_policy["event_set"] == "stable_significant" else "DAS"
            for comparison, top, source in _das_top10_frames(values, das_policy["event_set"]):
                story.append(Paragraph(f"Top {len(top)} {label} events | {html.escape(comparison)}", styles["Heading3"]))
                story.append(Paragraph("Ranked within the comparison by native FDR, absolute PSI effect and event ID. Complete fields remain in the HTML and original table.", styles["BodyText"]))
                for rank, (_, row) in enumerate(top.iterrows(), 1):
                    text = (f"{rank}. {row.get('gene_symbol', '')} | {row['event_id']} | "
                            f"dPSI {row['dpsi']} | FDR {row['fdr']}")
                    story.append(Paragraph(html.escape(text), styles["BodyText"]))
        if layer:
            story.append(Paragraph(f"{layer['number']}. {html.escape(layer['title_en'])}", styles['Heading3']))
            story.append(Paragraph(html.escape(layer['boundary_en']), styles['BodyText']))
        label_rows = _network_label_policy_records(analysis_id, values)
        if label_rows:
            story.append(Paragraph("Recorded optional network labels", styles["Heading3"]))
            for row in label_rows:
                story.append(Paragraph(
                    f"Source run: {html.escape(row['run'])}; policy: {html.escape(row['policy_id'])}; "
                    f"threshold labels: {row['apply_threshold_labels']}; hard filtering: {row['hard_filter']}. "
                    f"Recorded label rules: {html.escape(row['label_thresholds'])}.",
                    styles["BodyText"],
                ))
            story.append(Paragraph(html.escape(NETWORK_LABEL_METHODS), styles["BodyText"]))
        parameter_rows = _network_module_parameter_records(analysis_id, values)
        if parameter_rows:
            story.append(Paragraph("Recorded actual module parameters", styles["Heading3"]))
            story.append(Paragraph(
                "Values are taken from the saved network run. Missing historical values are not inferred from current defaults or requested controls.",
                styles["BodyText"],
            ))
            for row in parameter_rows:
                story.append(Paragraph(
                    f"<b>Source run:</b> {html.escape(row['run'])}<br/>"
                    f"<b>Minimum module size (events):</b> {html.escape(row['min_module_size'])}<br/>"
                    f"<b>Module merge cut height:</b> {html.escape(row['merge_cut_height'])}",
                    styles["BodyText"],
                ))
        for item in values:
            image_path = _png_path(item)
            if image_path:
                with PillowImage.open(image_path) as image:
                    width, height = image.size
                scale = min(174 * mm / width, 210 * mm / height)
                story.extend([Paragraph(html.escape(item.name), styles["Heading3"]),
                              Image(str(image_path), width=width * scale, height=height * scale),
                              Paragraph(f"Source run: {item.run_fingerprint[:12]}", styles["BodyText"]), PageBreak()])
        metadata = values[0].metadata
        for label, key in (("Question", "question"), ("Method", "method"), ("Interpretation", "interpretation"), ("Boundary", "limitations")):
            if metadata.get(key):
                story.append(Paragraph(f"<b>{label}.</b> {html.escape(str(metadata[key]))}", styles["BodyText"]))
        data = [["Artifact", "Kind", "Rows"]] + [
            [Paragraph(html.escape(item.name), styles["BodyText"]), item.kind.value, str(item.row_count or "")]
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
        if context.manifest.metadata.get("gene_filter"):
            for name in ("gene_filter.json", "FOOTNOTE.txt", "input_filter_counts.tsv", "annotation_decisions.tsv"):
                path = context.manifest.output_root / "gene_filter" / name
                if path.is_file():
                    archive.write(path, f"provenance/gene_filter/{name}")
        if include_files:
            seen: set[Path] = set()
            for artifact in artifacts:
                source = artifact.path.resolve()
                if source in seen or not source.is_file() or is_retired_report(source):
                    continue
                seen.add(source)
                safe_analysis = artifact.analysis_id.replace("/", "_").replace("\\", "_")
                run_directory = context.results.artifact_directory(artifact.analysis_id, artifact.run_fingerprint).resolve()
                preserve_links = artifact.metadata.get("preserve_relative_report_links") is True
                if preserve_links:
                    if not source.is_relative_to(run_directory):
                        raise ValueError("Portable report HTML must be inside its immutable run directory")
                    archive_name = f"artifacts/{safe_analysis}/{artifact.run_fingerprint[:12]}/{source.relative_to(run_directory).as_posix()}"
                else:
                    archive_name = f"artifacts/{safe_analysis}/{artifact.artifact_id[:8]}_{source.name}"
                archive.write(
                    source,
                    archive_name,
                )
                extras = [Path(path) for path in artifact.metadata.get("exports", {}).values()]
                extras.append(run_directory / "parameters.json")
                for path in extras:
                    path = path.resolve()
                    if path in seen or not path.is_file() or is_retired_report(path) or not path.is_relative_to(context.manifest.output_root.resolve()):
                        continue
                    seen.add(path)
                    relative = path.relative_to(run_directory).as_posix() if path.is_relative_to(run_directory) else f"{artifact.artifact_id[:8]}_{path.name}"
                    archive.write(path, f"artifacts/{safe_analysis}/{artifact.run_fingerprint[:12]}/{relative}")
