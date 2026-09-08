from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from katana_next.domain.models import AnalysisArtifact
from katana_next.infrastructure.jobs import JobContext
from katana_next.science.bulk import (
    STABLE_PRESETS,
    apply_stable_event_qc,
    build_candidate_evidence,
    cross_comparison_summary,
    landscape_tables,
    load_bulk_events,
    load_deg_results,
    load_dtu_results,
    selected_contrast_ids,
    significance_columns,
)
from katana_next.services.project_service import ProjectContext
from katana_next.workflows.base import (
    ArtifactWriter,
    ParameterDefinition,
    WorkflowBase,
    WorkflowDefinition,
)


INK = "#17222d"
TEAL = "#0f766e"
ORANGE = "#d97706"
BLUE = "#2563eb"
RED = "#dc2626"
GRID = "#dfe5e3"


def _base_config() -> dict[str, Any]:
    return {
        "background": "transparent",
        "view": {"stroke": None},
        "axis": {
            "domainColor": "#64748b",
            "gridColor": GRID,
            "labelColor": INK,
            "titleColor": INK,
            "labelFont": "IBM Plex Sans",
            "titleFont": "IBM Plex Sans",
        },
        "legend": {
            "labelFont": "IBM Plex Sans",
            "titleFont": "IBM Plex Sans",
            "labelColor": INK,
            "titleColor": INK,
        },
        "title": {"font": "Newsreader", "color": INK, "fontSize": 18},
    }


def _plot(
    values: list[dict[str, Any]],
    mark: dict[str, Any] | str,
    encoding: dict[str, Any],
    *,
    title: str,
    description: str,
    height: int = 320,
) -> dict[str, Any]:
    return {
        "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
        "description": description,
        "autosize": {"type": "fit", "contains": "padding"},
        "width": "container",
        "height": height,
        "title": title,
        "data": {"values": values},
        "mark": mark,
        "encoding": encoding,
        "config": _base_config(),
    }


COMMON_THRESHOLDS = (
    ParameterDefinition(
        "stability_profile",
        "Stable-AS profile",
        "select",
        "balanced",
        "Three declared measurement-quality presets plus Custom.",
        ("exploratory", "balanced", "stringent", "custom"),
    ),
    ParameterDefinition(
        "fdr", "rMATS FDR", "number", 0.05, "BH-adjusted caller FDR cutoff.", minimum=0, maximum=1, step=0.01
    ),
    ParameterDefinition(
        "min_abs_dpsi", "Minimum |dPSI|", "number", 0.10, "Minimum absolute standardized effect.", minimum=0, maximum=1, step=0.01
    ),
)


class ProjectAuditWorkflow(WorkflowBase):
    version = "2.0.0"
    definition = WorkflowDefinition(
        workflow_id="project.audit",
        workspace_id="project",
        label="Audit project design",
        question="Are inputs, samples, references and contrast directions explicit and usable?",
        method="Schema, path and experimental-unit validation without running an analysis.",
        output_summary="Input audit, sample sheet, contrast design and sample-count plot.",
        analysis_unit="input asset and biological sample",
        expected_artifacts=4,
    )

    def input_paths(self, context: ProjectContext, parameters: dict[str, Any]) -> list[Path]:
        del parameters
        return [item.path for item in context.snapshot.input_fingerprints if item.exists]

    def run(self, context: ProjectContext, job: JobContext, run_fingerprint: str, parameters: dict[str, Any]) -> list[str]:
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
                "interpretation": "Passing confirms input contracts only; it is not evidence of biological signal.",
            },
        )
        job.report(0.08, "Checking input paths and schemas")
        audit_rows = []
        for asset in context.manifest.assets:
            status, message = _audit_asset(asset.kind, asset.path)
            audit_rows.append(
                {
                    "row_id": asset.asset_id,
                    "asset_id": asset.asset_id,
                    "kind": asset.kind,
                    "path": str(asset.path),
                    "required": asset.required,
                    "status": status,
                    "message": message,
                    "contrast_id": str(asset.metadata.get("contrast_id", "")),
                }
            )
        for sample in context.manifest.samples:
            for kind, path in sample.assets.items():
                status = "pass" if path.is_file() else "fail"
                audit_rows.append(
                    {
                        "row_id": f"{sample.sample_id}::{kind}",
                        "asset_id": sample.sample_id,
                        "kind": kind,
                        "path": str(path),
                        "required": True,
                        "status": status,
                        "message": "Readable sample asset" if status == "pass" else "Missing sample asset",
                        "contrast_id": "",
                    }
                )
        audit = pd.DataFrame(audit_rows)
        blocking = audit["required"] & audit["status"].eq("fail") if not audit.empty else pd.Series(dtype=bool)
        if not audit.empty and bool(blocking.any()):
            raise ValueError("Project audit found required missing or incompatible inputs.")
        job.check_cancelled()
        samples = pd.DataFrame(
            [
                {
                    "row_id": sample.sample_id,
                    "sample_id": sample.sample_id,
                    "biological_sample_id": sample.biological_sample_id,
                    "condition": sample.condition,
                    "batch": sample.batch or "",
                    **{f"attribute_{key}": value for key, value in sample.attributes.items()},
                    "available_assets": ", ".join(sorted(sample.assets)),
                }
                for sample in context.manifest.samples
            ]
        )
        contrasts = pd.DataFrame(
            [
                {
                    "row_id": contrast.contrast_id,
                    "contrast_id": contrast.contrast_id,
                    "label": contrast.label,
                    "numerator": contrast.numerator.label,
                    "denominator": contrast.denominator.label,
                    "numerator_samples": ", ".join(contrast.numerator.sample_ids),
                    "denominator_samples": ", ".join(contrast.denominator.sample_ids),
                    "effect_definition": f"{contrast.numerator.label} minus {contrast.denominator.label}",
                }
                for contrast in context.manifest.contrasts
            ]
        )
        job.report(0.45, "Writing validated design artifacts")
        audit_artifact = writer.table("Input contract audit", "input_audit.parquet", audit, entity_type="input_asset")
        sample_artifact = writer.table("Biological samples", "samples.parquet", samples, entity_type="sample")
        writer.table("Declared contrasts", "contrasts.parquet", contrasts, entity_type="contrast")
        counts = (
            samples.groupby("condition", dropna=False).size().rename("biological_samples").reset_index()
            if not samples.empty
            else pd.DataFrame(columns=["condition", "biological_samples"])
        )
        writer.plot(
            "Samples per condition",
            "sample_design.vl.json",
            _plot(
                counts.to_dict("records"),
                {"type": "bar", "cornerRadiusEnd": 4, "color": TEAL},
                {
                    "x": {"field": "condition", "type": "nominal", "axis": {"title": "Condition", "labelAngle": -25}},
                    "y": {"field": "biological_samples", "type": "quantitative", "axis": {"title": "Biological samples"}},
                    "tooltip": [{"field": "condition"}, {"field": "biological_samples"}],
                },
                title="Declared biological design",
                description="Biological sample counts from the immutable project manifest.",
            ),
            entity_type="sample_summary",
            parents=(sample_artifact.artifact_id, audit_artifact.artifact_id),
        )
        job.report(0.92, "Project design audit complete")
        return writer.finish()


def _audit_asset(kind: str, path: Path) -> tuple[str, str]:
    if not path.is_file():
        return "fail", "File does not exist"
    try:
        if kind == "rmats_result":
            columns = set(pd.read_csv(path, sep="\t", nrows=2).columns)
            missing = {"FDR", "IncLevelDifference", "IncLevel1", "IncLevel2"} - columns
            return ("fail", f"Missing rMATS columns: {sorted(missing)}") if missing else ("pass", "Compatible rMATS event table")
        if kind == "deg_result":
            columns = set(pd.read_csv(path, sep="\t", nrows=2).columns)
            if not ({"log2FoldChange", "log2FC"} & columns) or not ({"padj", "FDR", "qvalue"} & columns):
                return "fail", "No recognized DEG effect/FDR columns"
            return "pass", "Compatible DEG table"
        if kind in {"dtu_result", "dexseq_result"}:
            columns = set(pd.read_csv(path, sep="\t", nrows=2).columns)
            recognized = bool({"padj", "qvalue", "adj_pval", "FDR"} & columns)
            stager = {"geneID", "txID", "gene", "transcript"}.issubset(columns)
            return (
                ("pass", "Compatible stageR transcript-usage table")
                if stager
                else ("pass", "Compatible transcript/exon usage table")
                if recognized
                else ("fail", "No recognized adjusted-p column")
            )
    except Exception as exc:
        return "fail", f"Could not parse header: {exc}"
    return "pass", "Readable input"


class StableAsWorkflow(WorkflowBase):
    version = "2.0.0"
    definition = WorkflowDefinition(
        workflow_id="quality.stable_as",
        workspace_id="quality",
        label="Stable-AS measurement QC",
        question="Which events have adequate read support and reproducible replicate behavior?",
        method="Replicate support, pairwise and leave-one-out direction consistency, and count-informed PSI uncertainty.",
        output_summary="Event-level QC, attrition table and interactive attrition plot.",
        analysis_unit="rMATS event within one declared contrast",
        expected_artifacts=3,
        parameters=COMMON_THRESHOLDS + (
            ParameterDefinition("min_reads", "Custom reads per replicate", "number", 15.0, "Used only by Custom profile.", minimum=0, step=1),
            ParameterDefinition("min_supported_replicates", "Custom supported replicates", "number", 2, "Used only by Custom profile.", minimum=1, step=1),
            ParameterDefinition("min_direction_consistency", "Custom direction consistency", "number", 0.75, "Used only by Custom profile.", minimum=0, maximum=1, step=0.05),
            ParameterDefinition("max_measurement_se", "Custom maximum PSI SE", "number", 0.12, "Used only by Custom profile.", minimum=0, maximum=0.5, step=0.01),
        ),
        required_asset_kinds=("rmats_result",),
        limitations="rMATS IJC/SJC are caller-derived evidence counts; uncertainty is a measurement diagnostic, not a replacement for rMATS inference.",
    )

    def run(self, context: ProjectContext, job: JobContext, run_fingerprint: str, parameters: dict[str, Any]) -> list[str]:
        job.report(0.05, "Reading selected rMATS comparisons")
        data = load_bulk_events(context, parameters)
        job.report(0.25, f"Normalized rMATS source cache: {data.source_cache_status}")
        job.report(0.35, "Evaluating replicate stability and read support")
        evaluated = significance_columns(apply_stable_event_qc(data.events, parameters), parameters)
        stages = [
            ("All tested", pd.Series(True, index=evaluated.index)),
            ("Replicate count", evaluated["passes_replicate_count"]),
            ("Read support", evaluated["passes_replicate_count"] & evaluated["passes_read_support"]),
            ("Direction consistency", evaluated["passes_replicate_count"] & evaluated["passes_read_support"] & evaluated["passes_direction_consistency"] & evaluated["passes_leave_one_out"]),
            ("Measurement uncertainty", evaluated["stable_event"]),
            ("Stable + significant", evaluated["stable_significant"]),
        ]
        attrition = pd.concat(
            [
                evaluated.loc[mask].groupby("comparison_id").size().rename("events").reset_index().assign(stage=stage, stage_order=index)
                for index, (stage, mask) in enumerate(stages)
            ],
            ignore_index=True,
        )
        attrition["row_id"] = attrition["comparison_id"] + "::" + attrition["stage"]
        writer = ArtifactWriter(
            context, self.definition.workflow_id, run_fingerprint, parameters, self.version,
            common_metadata={
                "question": self.definition.question,
                "method": self.definition.method,
                "analysis_unit": self.definition.analysis_unit,
                "interpretation": "Stable means measurable and replicate-supported. Stable does not itself mean differential.",
                "limitations": self.definition.limitations,
                "preset_values": STABLE_PRESETS,
                "source_cache_status": data.source_cache_status,
            },
        )
        table = writer.table("Stable-AS event QC", "stable_as_events.parquet", evaluated, entity_type="splice_event")
        writer.table("Stable-AS attrition", "stable_as_attrition.parquet", attrition, entity_type="qc_summary", parents=(table.artifact_id,))
        writer.plot(
            "Stable-AS attrition",
            "stable_as_attrition.vl.json",
            _plot(
                attrition.to_dict("records"),
                {"type": "line", "point": {"filled": True, "size": 70}, "strokeWidth": 2.5},
                {
                    "x": {"field": "stage", "type": "ordinal", "sort": {"field": "stage_order"}, "axis": {"title": None, "labelAngle": -20}},
                    "y": {"field": "events", "type": "quantitative", "axis": {"title": "Events retained"}},
                    "color": {"field": "comparison_id", "type": "nominal", "scale": {"scheme": "tableau10"}, "legend": {"title": "Contrast"}},
                    "tooltip": [{"field": "comparison_id"}, {"field": "stage"}, {"field": "events"}],
                },
                title="Stable-AS QC attrition",
                description="Independent measurement-quality gates followed by differential significance.",
                height=360,
            ),
            entity_type="qc_summary",
            parents=(table.artifact_id,),
        )
        job.report(0.92, "Stable-AS QC complete")
        return writer.finish()


class LandscapeWorkflow(WorkflowBase):
    version = "2.0.0"
    definition = WorkflowDefinition(
        workflow_id="landscape.bulk",
        workspace_id="landscape",
        label="Bulk splicing landscape",
        question="Do samples and event classes show coherent broad splicing structure?",
        method="Deduplicated sample-level PSI matrix, Spearman similarity and PCA on variable measured events.",
        output_summary="PSI matrix, duplicate-source audit, PCA, similarity and event composition.",
        analysis_unit="biological sample and exact splice event",
        expected_artifacts=7,
        parameters=(ParameterDefinition("max_pca_events", "Variable events for PCA", "number", 2000, "Highest-variance measurable events only.", minimum=50, maximum=20000, step=50),),
        required_asset_kinds=("rmats_result",),
        limitations="Pairwise rMATS outputs can repeat a sample/event estimate. Katana reports disagreement and uses the median rather than silently duplicating observations.",
    )

    def run(self, context: ProjectContext, job: JobContext, run_fingerprint: str, parameters: dict[str, Any]) -> list[str]:
        job.report(0.05, "Building deduplicated sample-level PSI matrix")
        data = load_bulk_events(context, parameters)
        job.report(0.22, f"Normalized rMATS source cache: {data.source_cache_status}")
        summary, pca, similarity = landscape_tables(
            data, max_pca_events=int(parameters.get("max_pca_events", 2000))
        )
        sample_condition = {sample.sample_id: sample.condition for sample in context.manifest.samples}
        if not pca.empty:
            pca["condition"] = pca["sample_id"].map(sample_condition).fillna("unmapped")
        writer = ArtifactWriter(
            context, self.definition.workflow_id, run_fingerprint, parameters, self.version,
            common_metadata={
                "question": self.definition.question,
                "method": self.definition.method,
                "analysis_unit": self.definition.analysis_unit,
                "interpretation": "PCA and similarity are descriptive QC. Separation can reflect biology or technical structure and requires metadata review.",
                "limitations": self.definition.limitations,
                "source_cache_status": data.source_cache_status,
            },
        )
        psi = data.psi_matrix.reset_index()
        psi["row_id"] = psi["event_id"]
        psi_artifact = writer.table("Sample-level PSI matrix", "sample_level_psi.parquet", psi, entity_type="splice_event")
        writer.table("Duplicate PSI source audit", "psi_source_audit.parquet", data.matrix_audit, entity_type="measurement_audit", parents=(psi_artifact.artifact_id,))
        writer.table("Event composition", "event_composition.parquet", summary.assign(row_id=lambda x: x["comparison_id"] + "::" + x["event_type"]), entity_type="event_summary")
        pca_artifact = writer.table("Sample PSI PCA", "sample_pca.parquet", pca, entity_type="sample")
        similarity_artifact = writer.table("Sample PSI similarity", "sample_similarity.parquet", similarity, entity_type="sample_pair")
        writer.plot(
            "Sample PSI PCA",
            "sample_pca.vl.json",
            _plot(
                pca.to_dict("records"),
                {"type": "point", "filled": True, "size": 120, "stroke": "white", "strokeWidth": 1},
                {
                    "x": {"field": "PC1", "type": "quantitative", "axis": {"title": "PC1"}},
                    "y": {"field": "PC2", "type": "quantitative", "axis": {"title": "PC2"}},
                    "color": {"field": "condition", "type": "nominal", "scale": {"scheme": "tableau10"}},
                    "tooltip": [{"field": "sample_id"}, {"field": "condition"}, {"field": "PC1", "format": ".3f"}, {"field": "PC2", "format": ".3f"}, {"field": "events_used"}],
                },
                title="Sample-level splicing structure",
                description="PCA on the most variable deduplicated event PSI measurements.",
                height=390,
            ),
            entity_type="sample_summary",
            parents=(pca_artifact.artifact_id, psi_artifact.artifact_id),
        )
        writer.plot(
            "Sample PSI similarity",
            "sample_similarity.vl.json",
            _plot(
                similarity.to_dict("records"),
                "rect",
                {
                    "x": {"field": "sample_a", "type": "nominal", "axis": {"title": None, "labelAngle": -45}},
                    "y": {"field": "sample_b", "type": "nominal", "axis": {"title": None}},
                    "color": {"field": "spearman_r", "type": "quantitative", "scale": {"domain": [-1, 1], "scheme": "redblue", "reverse": True}, "legend": {"title": "Spearman r"}},
                    "tooltip": [{"field": "sample_a"}, {"field": "sample_b"}, {"field": "spearman_r", "format": ".3f"}, {"field": "shared_events"}],
                },
                title="Pairwise sample splicing similarity",
                description="Spearman correlation on shared measured exact events.",
                height=470,
            ),
            entity_type="sample_pair_summary",
            parents=(similarity_artifact.artifact_id, psi_artifact.artifact_id),
        )
        job.report(0.92, "Bulk splicing landscape complete")
        return writer.finish()


class DifferentialAsWorkflow(WorkflowBase):
    version = "2.0.0"
    definition = WorkflowDefinition(
        workflow_id="differential.as",
        workspace_id="differential",
        label="Differential alternative splicing",
        question="Which exact AS events change in each declared contrast?",
        method="Imported rMATS effect/FDR with Katana replicate-stability annotation and declared numerator-minus-denominator direction.",
        output_summary="One event table, contrast summaries, volcano and event-count plots.",
        analysis_unit="exact splice event within a contrast",
        expected_artifacts=4,
        parameters=COMMON_THRESHOLDS,
        required_asset_kinds=("rmats_result",),
        limitations="The inferential p-value remains the imported rMATS result; Katana stability is an independent measurement-support annotation.",
    )

    def run(self, context: ProjectContext, job: JobContext, run_fingerprint: str, parameters: dict[str, Any]) -> list[str]:
        job.report(0.05, "Normalizing rMATS exact events")
        data = load_bulk_events(context, parameters)
        job.report(0.24, f"Normalized rMATS source cache: {data.source_cache_status}")
        events = significance_columns(apply_stable_event_qc(data.events, parameters), parameters)
        events["minus_log10_fdr"] = -np.log10(events["fdr"].clip(lower=np.finfo(float).tiny))
        summary = (
            events.groupby(["comparison_id", "event_type"], dropna=False)
            .agg(tested_events=("event_id", "size"), significant_events=("significant", "sum"), stable_significant_events=("stable_significant", "sum"))
            .reset_index()
        )
        summary["row_id"] = summary["comparison_id"] + "::" + summary["event_type"]
        writer = ArtifactWriter(
            context, self.definition.workflow_id, run_fingerprint, parameters, self.version,
            common_metadata={
                "question": self.definition.question,
                "method": self.definition.method,
                "analysis_unit": self.definition.analysis_unit,
                "direction_definition": "dPSI = numerator minus denominator from the manifest contrast",
                "interpretation": "Use stable_significant for the strict event set. significant retains caller significance even when replicate-support QC fails.",
                "limitations": self.definition.limitations,
                "source_cache_status": data.source_cache_status,
            },
        )
        table = writer.table("Differential AS events", "differential_as_events.parquet", events, entity_type="splice_event")
        summary_artifact = writer.table("Differential AS summary", "differential_as_summary.parquet", summary, entity_type="event_summary", parents=(table.artifact_id,))
        plot_values = events.loc[events["fdr"].notna() & events["dpsi"].notna(), ["row_id", "comparison_id", "event_type", "gene_symbol", "dpsi", "fdr", "minus_log10_fdr", "stable_significant", "significant"]]
        if len(plot_values) > 50000:
            plot_values = pd.concat([plot_values[plot_values["significant"]], plot_values[~plot_values["significant"]].sample(n=min(30000, int((~plot_values["significant"]).sum())), random_state=17)]).drop_duplicates("row_id")
        writer.plot(
            "Differential AS volcano",
            "differential_as_volcano.vl.json",
            _plot(
                plot_values.to_dict("records"),
                {"type": "point", "filled": True, "opacity": 0.58, "size": 35},
                {
                    "x": {"field": "dpsi", "type": "quantitative", "axis": {"title": "dPSI (numerator minus denominator)"}},
                    "y": {"field": "minus_log10_fdr", "type": "quantitative", "axis": {"title": "-log10 rMATS FDR"}},
                    "color": {"field": "stable_significant", "type": "nominal", "scale": {"domain": [True, False], "range": [RED, "#aab3ba"]}, "legend": {"title": "Stable + significant"}},
                    "facet": {"field": "comparison_id", "type": "nominal", "columns": 2, "header": {"title": None}},
                    "tooltip": [{"field": "gene_symbol"}, {"field": "event_type"}, {"field": "comparison_id"}, {"field": "dpsi", "format": ".3f"}, {"field": "fdr", "format": ".2e"}, {"field": "stable_significant"}],
                },
                title="Differential alternative splicing",
                description="Imported rMATS significance with independent Katana stability annotation.",
                height=260,
            ),
            entity_type="event_summary",
            parents=(table.artifact_id,),
        )
        writer.plot(
            "Stable significant events by type",
            "differential_as_counts.vl.json",
            _plot(
                summary.to_dict("records"),
                {"type": "bar", "cornerRadiusEnd": 3},
                {
                    "x": {"field": "comparison_id", "type": "nominal", "axis": {"title": None, "labelAngle": -25}},
                    "y": {"field": "stable_significant_events", "type": "quantitative", "axis": {"title": "Stable significant events"}},
                    "color": {"field": "event_type", "type": "nominal", "scale": {"scheme": "tableau10"}},
                    "xOffset": {"field": "event_type"},
                    "tooltip": [{"field": "comparison_id"}, {"field": "event_type"}, {"field": "tested_events"}, {"field": "significant_events"}, {"field": "stable_significant_events"}],
                },
                title="Strict AS event set by contrast",
                description="Events pass rMATS FDR, effect-size and Katana stable-AS gates.",
                height=360,
            ),
            entity_type="event_summary",
            parents=(summary_artifact.artifact_id,),
        )
        job.report(0.92, "Differential AS complete")
        return writer.finish()


class DegWorkflow(WorkflowBase):
    version = "2.0.0"
    definition = WorkflowDefinition(
        workflow_id="differential.deg",
        workspace_id="differential",
        label="Differential gene expression",
        question="Which genes change expression in the selected contrasts?",
        method="Transparent import of DESeq2-like effect sizes and adjusted p-values.",
        output_summary="DEG table, summary and volcano plot.",
        analysis_unit="gene within a contrast",
        expected_artifacts=3,
        parameters=(
            ParameterDefinition("fdr", "DEG FDR", "number", 0.05, "Adjusted p-value cutoff.", minimum=0, maximum=1, step=0.01),
            ParameterDefinition("min_abs_log2fc", "Minimum |log2FC|", "number", 1.0, "Minimum expression effect size.", minimum=0, step=0.1),
        ),
        required_asset_kinds=("deg_result",),
    )

    def run(self, context: ProjectContext, job: JobContext, run_fingerprint: str, parameters: dict[str, Any]) -> list[str]:
        contrast_ids = selected_contrast_ids(context, parameters)
        job.report(0.10, "Reading DEG results")
        frame = load_deg_results(context, contrast_ids)
        if frame.empty:
            raise ValueError("No compatible DEG result is available for selected contrasts.")
        frame["significant"] = frame["deg_fdr"].le(float(parameters.get("fdr", 0.05))) & frame["log2fc"].abs().ge(float(parameters.get("min_abs_log2fc", 1.0)))
        frame["direction"] = np.where(frame["log2fc"].ge(0), "higher in numerator", "lower in numerator")
        frame["minus_log10_fdr"] = -np.log10(frame["deg_fdr"].clip(lower=np.finfo(float).tiny))
        summary = frame.groupby(["comparison_id", "direction"], dropna=False).agg(genes=("gene_id", "nunique"), significant_genes=("significant", "sum")).reset_index()
        summary["row_id"] = summary["comparison_id"] + "::" + summary["direction"]
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version, common_metadata={"question": self.definition.question, "method": self.definition.method, "analysis_unit": self.definition.analysis_unit, "direction_definition": "log2FC = numerator minus denominator on the log2 expression scale"})
        table = writer.table("Differentially expressed genes", "deg_results.parquet", frame, entity_type="gene")
        writer.table("DEG summary", "deg_summary.parquet", summary, entity_type="gene_summary", parents=(table.artifact_id,))
        values = frame.dropna(subset=["log2fc", "minus_log10_fdr"])
        if len(values) > 50000:
            nonsignificant = values[~values["significant"]]
            nonsignificant = nonsignificant.sample(
                min(30_000, len(nonsignificant)), random_state=17
            )
            values = pd.concat(
                [values[values["significant"]], nonsignificant]
            ).drop_duplicates("row_id")
        writer.plot(
            "DEG volcano",
            "deg_volcano.vl.json",
            _plot(values.to_dict("records"), {"type": "point", "filled": True, "opacity": 0.55, "size": 32}, {"x": {"field": "log2fc", "type": "quantitative", "axis": {"title": "log2 fold change (numerator / denominator)"}}, "y": {"field": "minus_log10_fdr", "type": "quantitative", "axis": {"title": "-log10 adjusted p-value"}}, "color": {"field": "significant", "type": "nominal", "scale": {"domain": [True, False], "range": [ORANGE, "#aab3ba"]}}, "facet": {"field": "comparison_id", "columns": 2}, "tooltip": [{"field": "gene_symbol"}, {"field": "gene_id"}, {"field": "comparison_id"}, {"field": "log2fc", "format": ".3f"}, {"field": "deg_fdr", "format": ".2e"}]}, title="Differential gene expression", description="Imported DESeq2-like gene-level statistics.", height=260),
            entity_type="gene_summary",
            parents=(table.artifact_id,),
        )
        job.report(0.92, "DEG import complete")
        return writer.finish()


class DtuWorkflow(WorkflowBase):
    version = "2.0.0"
    definition = WorkflowDefinition(
        workflow_id="differential.dtu",
        workspace_id="differential",
        label="Differential transcript usage",
        question="Which genes and transcripts show DTU or exon-usage support?",
        method="Schema-normalized import of DTU, stageR or DEXSeq adjusted p-values.",
        output_summary="Transcript/exon evidence table and gene-level support summary.",
        analysis_unit="transcript or exon feature within a contrast",
        expected_artifacts=2,
        parameters=(ParameterDefinition("fdr", "DTU/DEXSeq FDR", "number", 0.05, "Adjusted p-value cutoff.", minimum=0, maximum=1, step=0.01),),
        optional_asset_kinds=("dtu_result", "dexseq_result"),
    )

    def run(self, context: ProjectContext, job: JobContext, run_fingerprint: str, parameters: dict[str, Any]) -> list[str]:
        contrast_ids = selected_contrast_ids(context, parameters)
        frame = load_dtu_results(context, contrast_ids)
        if frame.empty:
            raise ValueError("No compatible DTU or DEXSeq result is available for selected contrasts.")
        frame["significant"] = frame["dtu_fdr"].le(float(parameters.get("fdr", 0.05)))
        summary = frame.groupby(["comparison_id", "gene_id"], dropna=False).agg(best_fdr=("dtu_fdr", "min"), tested_features=("row_id", "size"), significant_features=("significant", "sum")).reset_index()
        summary["row_id"] = summary["comparison_id"] + "::" + summary["gene_id"]
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version, common_metadata={"question": self.definition.question, "method": self.definition.method, "analysis_unit": self.definition.analysis_unit, "interpretation": "DTU/DEXSeq supports transcript or exon usage. It does not replace event-level DAS."})
        detail = writer.table("DTU and exon-usage evidence", "dtu_features.parquet", frame, entity_type="transcript_feature")
        writer.table("Gene-level DTU support", "dtu_gene_summary.parquet", summary, entity_type="gene", parents=(detail.artifact_id,))
        return writer.finish()


def _latest_table(context: ProjectContext, analysis_id: str, name: str | None = None) -> AnalysisArtifact | None:
    artifacts = context.catalog.list_artifacts(context.manifest.project_id, analysis_id=analysis_id, kind="table")
    artifacts = [
        artifact
        for artifact in artifacts
        if artifact.metadata.get("project_snapshot_fingerprint")
        == context.snapshot.manifest_fingerprint
    ]
    if name:
        artifacts = [artifact for artifact in artifacts if artifact.name == name]
    return artifacts[0] if artifacts else None


class CrossComparisonWorkflow(WorkflowBase):
    version = "2.0.0"
    definition = WorkflowDefinition(
        workflow_id="differential.cross_comparison",
        workspace_id="differential",
        label="Exact cross-comparison AS",
        question="How much of each strict AS signature is shared, directionally concordant or contrast-specific?",
        method="Exact normalized event identity, declared tested-event universe, hypergeometric overlap, direction concordance and effect correlation.",
        output_summary="Pair summary, event membership, overlap profile and similarity heatmap.",
        analysis_unit="pair of contrasts and exact splice event",
        expected_artifacts=4,
        limitations="Same-gene and nearby-coordinate events are not silently treated as the same event.",
    )

    def parent_artifacts(self, context: ProjectContext, parameters: dict[str, Any]) -> list[AnalysisArtifact]:
        del parameters
        artifact = _latest_table(context, "differential.as", "Differential AS events")
        return [artifact] if artifact else []

    def input_paths(self, context: ProjectContext, parameters: dict[str, Any]) -> list[Path]:
        return [artifact.path for artifact in self.parent_artifacts(context, parameters)]

    def run(self, context: ProjectContext, job: JobContext, run_fingerprint: str, parameters: dict[str, Any]) -> list[str]:
        parent = _latest_table(context, "differential.as", "Differential AS events")
        if parent is None:
            raise ValueError("Run Differential alternative splicing first. Cross-comparison never starts it implicitly.")
        events = pd.read_parquet(parent.path)
        selected = parameters.get("contrast_ids")
        if selected:
            events = events[events["comparison_id"].isin(selected)]
        if events["comparison_id"].nunique() < 2:
            raise ValueError("Select at least two contrasts for cross-comparison.")
        job.report(0.25, "Matching exact event structures")
        summary, membership = cross_comparison_summary(events)
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version, common_metadata={"question": self.definition.question, "method": self.definition.method, "analysis_unit": self.definition.analysis_unit, "interpretation": "Overlap, direction agreement and effect correlation answer different questions and are kept separate.", "limitations": self.definition.limitations})
        summary_artifact = writer.table("Cross-comparison exact-event summary", "cross_comparison_summary.parquet", summary, entity_type="contrast_pair", parents=(parent.artifact_id,))
        writer.table("Cross-comparison event membership", "cross_comparison_events.parquet", membership, entity_type="splice_event", parents=(summary_artifact.artifact_id, parent.artifact_id))
        long_values = []
        for row in summary.to_dict("records"):
            for category, field in (("Shared exact", "shared_exact_events"), ("A-specific", "a_only_events"), ("B-specific", "b_only_events")):
                long_values.append({"pair": f"{row['comparison_a']}\nvs\n{row['comparison_b']}", "category": category, "events": row[field]})
        writer.plot(
            "Exact-event overlap profile",
            "cross_comparison_overlap.vl.json",
            _plot(long_values, {"type": "bar", "cornerRadiusEnd": 3}, {"x": {"field": "events", "type": "quantitative", "axis": {"title": "Stable significant exact events"}}, "y": {"field": "pair", "type": "nominal", "axis": {"title": None, "labelLimit": 320}}, "color": {"field": "category", "type": "nominal", "scale": {"domain": ["Shared exact", "A-specific", "B-specific"], "range": [TEAL, BLUE, ORANGE]}}, "tooltip": [{"field": "pair"}, {"field": "category"}, {"field": "events"}]}, title="Exact AS overlap across selected contrasts", description="A-specific and B-specific sets are shown separately; direction is not part of event identity.", height=max(300, 65 * len(summary))),
            entity_type="contrast_pair_summary",
            parents=(summary_artifact.artifact_id,),
        )
        matrix_values = []
        for row in summary.to_dict("records"):
            label_a, label_b = row["comparison_a"], row["comparison_b"]
            matrix_values.extend([
                {"comparison_a": label_a, "comparison_b": label_b, "jaccard": row["jaccard"], "direction_concordance": row["direction_concordance"]},
                {"comparison_a": label_b, "comparison_b": label_a, "jaccard": row["jaccard"], "direction_concordance": row["direction_concordance"]},
            ])
        writer.plot(
            "AS similarity matrix",
            "cross_comparison_similarity.vl.json",
            _plot(matrix_values, "rect", {"x": {"field": "comparison_a", "type": "nominal", "axis": {"title": None, "labelAngle": -45, "labelLimit": 220}}, "y": {"field": "comparison_b", "type": "nominal", "axis": {"title": None, "labelLimit": 220}}, "color": {"field": "jaccard", "type": "quantitative", "scale": {"domain": [0, 1], "scheme": "tealblues"}, "legend": {"title": "Jaccard"}}, "tooltip": [{"field": "comparison_a"}, {"field": "comparison_b"}, {"field": "jaccard", "format": ".3f"}, {"field": "direction_concordance", "format": ".1%"}]}, title="Strict AS signature similarity", description="Jaccard quantifies set similarity; tooltip separately reports same-sign dPSI among shared exact events.", height=500),
            entity_type="contrast_pair_summary",
            parents=(summary_artifact.artifact_id,),
        )
        job.report(0.92, "Exact cross-comparison complete")
        return writer.finish()


class CandidateEvidenceWorkflow(WorkflowBase):
    version = "2.0.0"
    definition = WorkflowDefinition(
        workflow_id="prioritize.evidence",
        workspace_id="prioritize",
        label="Transparent candidate evidence",
        question="Which genes have convergent AS, DEG and transcript-usage evidence?",
        method="Gene-by-contrast evidence matrix with no hidden weighted ranking score.",
        output_summary="Candidate matrix and evidence heatmap linked to source artifacts.",
        analysis_unit="gene within a contrast",
        expected_artifacts=2,
        parameters=(ParameterDefinition("fdr", "Evidence FDR", "number", 0.05, "Common display threshold for imported evidence.", minimum=0, maximum=1, step=0.01), ParameterDefinition("deg_min_abs_log2fc", "DEG minimum |log2FC|", "number", 1.0, "DEG effect-size gate.", minimum=0, step=0.1)),
    )

    def parent_artifacts(self, context: ProjectContext, parameters: dict[str, Any]) -> list[AnalysisArtifact]:
        del parameters
        artifacts = []
        for analysis_id, name in (("differential.as", "Differential AS events"), ("differential.deg", "Differentially expressed genes"), ("differential.dtu", "DTU and exon-usage evidence")):
            artifact = _latest_table(context, analysis_id, name)
            if artifact:
                artifacts.append(artifact)
        return artifacts

    def input_paths(self, context: ProjectContext, parameters: dict[str, Any]) -> list[Path]:
        return [artifact.path for artifact in self.parent_artifacts(context, parameters)]

    def run(self, context: ProjectContext, job: JobContext, run_fingerprint: str, parameters: dict[str, Any]) -> list[str]:
        as_artifact = _latest_table(context, "differential.as", "Differential AS events")
        if as_artifact is None:
            raise ValueError("Run Differential alternative splicing first. Candidate evidence does not trigger it implicitly.")
        deg_artifact = _latest_table(context, "differential.deg", "Differentially expressed genes")
        dtu_artifact = _latest_table(context, "differential.dtu", "DTU and exon-usage evidence")
        events = pd.read_parquet(as_artifact.path)
        deg = pd.read_parquet(deg_artifact.path) if deg_artifact else pd.DataFrame()
        dtu = pd.read_parquet(dtu_artifact.path) if dtu_artifact else pd.DataFrame()
        selected = parameters.get("contrast_ids")
        if selected:
            events = events[events["comparison_id"].isin(selected)]
            if not deg.empty:
                deg = deg[deg["comparison_id"].isin(selected)]
            if not dtu.empty:
                dtu = dtu[dtu["comparison_id"].isin(selected)]
        job.report(0.35, "Joining transparent evidence dimensions")
        candidates = build_candidate_evidence(events, deg, dtu, parameters)
        long = candidates.melt(id_vars=["row_id", "comparison_id", "gene_symbol", "gene_id"], value_vars=["as_supported", "deg_supported", "dtu_supported"], var_name="evidence", value_name="supported")
        long = long[long["supported"]].copy()
        top_ids = set(candidates.head(60)["row_id"])
        long = long[long["row_id"].isin(top_ids)]
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version, common_metadata={"question": self.definition.question, "method": self.definition.method, "analysis_unit": self.definition.analysis_unit, "interpretation": "Convergent evidence increases prioritization confidence but does not merge distinct statistical tests or prove mechanism.", "missing_evidence": "Missing DEG or DTU remains missing and is never treated as a negative result."})
        table = writer.table("Candidate evidence matrix", "candidate_evidence.parquet", candidates, entity_type="gene", parents=tuple(artifact.artifact_id for artifact in self.parent_artifacts(context, parameters)))
        writer.plot(
            "Candidate evidence heatmap",
            "candidate_evidence.vl.json",
            _plot(long.to_dict("records"), {"type": "rect", "cornerRadius": 2}, {"x": {"field": "evidence", "type": "nominal", "axis": {"title": None, "labelExpr": "replace(datum.label, '_supported', '')"}}, "y": {"field": "gene_symbol", "type": "nominal", "sort": {"field": "comparison_id"}, "axis": {"title": None, "labelLimit": 130}}, "color": {"field": "comparison_id", "type": "nominal", "scale": {"scheme": "tableau10"}, "legend": {"title": "Contrast"}}, "tooltip": [{"field": "gene_symbol"}, {"field": "gene_id"}, {"field": "comparison_id"}, {"field": "evidence"}]}, title="Top candidate evidence dimensions", description="Up to 60 highest-evidence gene-by-contrast rows; the full table is never truncated.", height=max(360, min(1100, 18 * max(1, long["gene_symbol"].nunique())))),
            entity_type="gene_summary",
            parents=(table.artifact_id,),
        )
        return writer.finish()
