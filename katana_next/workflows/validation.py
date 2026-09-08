from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import pandas as pd

from katana_next.domain.models import AnalysisArtifact
from katana_next.infrastructure.jobs import JobContext
from katana_next.services.project_service import ProjectContext
from katana_next.workflows.base import (
    ArtifactWriter,
    ParameterDefinition,
    WorkflowBase,
    WorkflowDefinition,
)
from katana_next.workflows.bulk import _latest_table
from src.adapters.native_sashimi_adapter import NativeSashimiAdapter


def _safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")[:100] or "event"


class NativeSashimiWorkflow(WorkflowBase):
    version = "2.0.0"
    definition = WorkflowDefinition(
        workflow_id="validate.native_sashimi",
        workspace_id="validate",
        label="Native Sashimi read validation",
        question="Do BAM coverage and junction reads support selected differential AS events?",
        method="Per-replicate BAM coverage and junction counting with rMATS event geometry and optional GTF exon labels.",
        output_summary="Sashimi SVG/PNG previews, replicate tracks, junction tables and failure audit.",
        analysis_unit="selected exact splice event and biological replicate",
        expected_artifacts=1,
        parameters=(
            ParameterDefinition("selection_id", "Saved event selection", "selection", "", "A Selection Basket set created from the Differential AS table."),
            ParameterDefinition("max_events", "Maximum events", "number", 10, "Batch safety limit; each event remains an independent plot.", minimum=1, maximum=100, step=1),
            ParameterDefinition("view", "Primary view", "select", "readable", "Readable compresses long introns while retaining aligned exon anchors.", ("readable", "genomic")),
        ),
        required_asset_kinds=("rmats_result",),
        limitations="Coverage is read support, not an independent statistical replication of the rMATS test.",
    )

    def parent_artifacts(self, context: ProjectContext, parameters: dict[str, Any]) -> list[AnalysisArtifact]:
        del parameters
        artifact = _latest_table(context, "differential.as", "Differential AS events")
        return [artifact] if artifact else []

    def input_paths(self, context: ProjectContext, parameters: dict[str, Any]) -> list[Path]:
        paths = [artifact.path for artifact in self.parent_artifacts(context, parameters)]
        paths.extend(
            path
            for sample in context.manifest.samples
            for kind, path in sample.assets.items()
            if kind in {"bam", "bai"} and path.exists()
        )
        if context.manifest.annotation.gtf and context.manifest.annotation.gtf.exists():
            paths.append(context.manifest.annotation.gtf)
        return paths

    def run(self, context: ProjectContext, job: JobContext, run_fingerprint: str, parameters: dict[str, Any]) -> list[str]:
        differential = _latest_table(context, "differential.as", "Differential AS events")
        if differential is None:
            raise ValueError("Run Differential alternative splicing first. Native Sashimi never launches it implicitly.")
        selection_id = str(parameters.get("selection_id", "")).strip()
        selection = next(
            (
                item
                for item in context.catalog.list_selections(context.manifest.project_id)
                if item.selection_id == selection_id
            ),
            None,
        )
        if selection is None:
            raise ValueError("Choose a saved splice-event selection from the Selection Basket.")
        if selection.source_artifact_id != differential.artifact_id:
            source = context.catalog.get_artifact(selection.source_artifact_id)
            if source is None or source.entity_type != "splice_event":
                raise ValueError("The selected set is not derived from a splice-event table.")
            event_source = source
        else:
            event_source = differential
        events = pd.read_parquet(event_source.path)
        selected_ids = set(selection.entity_ids)
        id_column = "row_id" if selected_ids & set(events.get("row_id", pd.Series(dtype=str)).astype(str)) else "event_id"
        events = events[events[id_column].astype(str).isin(selected_ids)].head(int(parameters.get("max_events", 10)))
        if events.empty:
            raise ValueError("The saved selection has no rows in its source event artifact.")

        adapter = NativeSashimiAdapter(samtools_path=str(context.manifest.metadata.get("samtools_path", "samtools")))
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
                "interpretation": "Arc width is proportional to junction support. Coverage is normalized within the Native Sashimi implementation and shown by replicate.",
                "limitations": self.definition.limitations,
            },
        )
        manifest_rows = []
        event_roots: list[tuple[pd.Series, Path]] = []
        for index, (_, event) in enumerate(events.iterrows(), start=1):
            job.report(0.05 + 0.15 * index / len(events), f"Preparing event {index} of {len(events)}")
            contrast = next(
                item
                for item in context.manifest.contrasts
                if item.contrast_id == str(event["comparison_id"])
            )
            numerator_bams = _bam_paths(context, contrast.numerator.sample_ids)
            denominator_bams = _bam_paths(context, contrast.denominator.sample_ids)
            if not numerator_bams or not denominator_bams:
                raise ValueError(f"BAM assets are missing for contrast {contrast.contrast_id}.")
            slug = _safe_name(f"{event.get('gene_symbol', '')}_{event['event_type']}_{index}")
            input_dir = writer.directory / "inputs" / slug
            output_dir = writer.directory / "plots" / slug
            input_dir.mkdir(parents=True, exist_ok=True)
            event_file = input_dir / "event.tsv"
            pd.DataFrame([event]).to_csv(event_file, sep="\t", index=False)
            b1 = input_dir / "numerator_bams.txt"
            b2 = input_dir / "denominator_bams.txt"
            b1.write_text("\n".join(str(path) for path in numerator_bams), encoding="utf-8")
            b2.write_text("\n".join(str(path) for path in denominator_bams), encoding="utf-8")
            manifest_rows.append(
                {
                    "geneSymbol": event.get("gene_symbol", ""),
                    "GeneID": event.get("gene_id", ""),
                    "event_type": event["event_type"],
                    "event_uid": event["event_id"],
                    "comparison_id": contrast.contrast_id,
                    "comparison_name": contrast.label,
                    "dPSI": event.get("dpsi"),
                    "FDR": event.get("fdr"),
                    "event_file": str(event_file),
                    "b1_txt": str(b1),
                    "b2_txt": str(b2),
                    "label1": contrast.numerator.label,
                    "label2": contrast.denominator.label,
                    "outdir": str(output_dir),
                }
            )
            event_roots.append((event, output_dir / "native_sashimi"))

        job.report(0.22, f"Rendering {len(events)} selected event(s)")
        _, failures = adapter.run_manifest(
            pd.DataFrame(manifest_rows),
            writer.directory,
            context.manifest.annotation.gtf,
        )
        if not failures.empty:
            failure_artifact = writer.table(
                "Native Sashimi failure audit",
                "native_sashimi_failures.parquet",
                failures.assign(row_id=lambda frame: frame.index.astype(str)),
                entity_type="validation_failure",
                parents=(event_source.artifact_id,),
            )
            if len(failures) == len(events):
                raise RuntimeError(
                    f"All selected Native Sashimi events failed. Inspect {failure_artifact.path}."
                )
        view = str(parameters.get("view", "readable"))
        primary_name = "native_summary_readable.png" if view == "readable" else "native_summary.png"
        for event, root in event_roots:
            primary = root / primary_name
            if not primary.exists():
                continue
            writer.existing_file(
                f"{event.get('gene_symbol', event.get('gene_id', 'Event'))} {event['event_type']} Sashimi",
                primary,
                entity_type="splice_event_validation",
                metadata={
                    "mime_type": "image/png",
                    "event_id": str(event["event_id"]),
                    "comparison_id": str(event["comparison_id"]),
                    "source_event_artifact": event_source.artifact_id,
                },
                parents=(event_source.artifact_id,),
            )
            evidence = root / "native_rmats_evidence.png"
            if evidence.exists():
                writer.existing_file(
                    f"{event.get('gene_symbol', event.get('gene_id', 'Event'))} rMATS evidence",
                    evidence,
                    entity_type="splice_event_validation",
                    metadata={"mime_type": "image/png", "event_id": str(event["event_id"])},
                    parents=(event_source.artifact_id,),
                )
        if not writer.artifacts:
            raise RuntimeError("Native Sashimi produced no preview artifact.")
        job.report(0.95, "Native Sashimi validation complete")
        return writer.finish()


def _bam_paths(context: ProjectContext, sample_ids: tuple[str, ...]) -> list[Path]:
    lookup = {sample.sample_id: sample for sample in context.manifest.samples}
    paths = []
    for sample_id in sample_ids:
        sample = lookup.get(sample_id)
        if sample is None:
            continue
        bam = sample.assets.get("bam")
        if bam and bam.is_file():
            paths.append(bam)
    return paths
