from __future__ import annotations

from pathlib import Path
from dataclasses import replace
from typing import Any

import numpy as np
import pandas as pd

from katana_next.domain.models import AnalysisArtifact
from katana_next.infrastructure.jobs import JobContext
from katana_next.science.bulk import (
    STABLE_PRESETS,
    apply_stable_event_qc,
    build_candidate_evidence,
    candidate_display_labels,
    candidate_pair_plot_counts,
    pairwise_candidate_membership,
    cross_comparison_summary,
    landscape_tables,
    load_bulk_events,
    load_deg_results,
    load_dtu_results,
    selected_contrast_ids,
    significance_columns,
)
from katana_next.services.project_service import ProjectContext
from katana_next.workflows.sources import resolve_table, validate_source
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
    result = {
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
    if "facet" in encoding:
        facet = dict(encoding["facet"])
        columns = facet.pop("columns", 2)
        result.pop("mark")
        result.pop("encoding")
        result.pop("width")
        result.pop("height")
        result["facet"] = facet
        result["columns"] = columns
        result["spec"] = {"width": 430, "height": height, "mark": mark,
                          "encoding": {key: value for key, value in encoding.items() if key != "facet"}}
        result["autosize"] = {"type": "pad"}
    return result


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

QUALITY_SOURCE = ParameterDefinition(
    "quality_artifact_id", "Measurement policy from Quality", "artifact", "",
    "Use a completed Quality result. Empty selects the newest scope-compatible result; no upstream job is started.",
    ("quality.stable_as", "Stable-AS event QC"),
)


def _quality_source(context: ProjectContext, parameters: dict[str, Any]) -> AnalysisArtifact | None:
    identity = str(parameters.get("quality_artifact_id", "")).strip()
    if identity:
        source = context.catalog.get_artifact(identity)
        validate_source(context, source)
        if source.analysis_id != "quality.stable_as" or source.name != "Stable-AS event QC":
            raise ValueError("Select a Stable-AS event QC artifact.")
        bound = {**parameters, "source_artifact_ids": {"quality.stable_as::Stable-AS event QC": identity}}
        return resolve_table(context, "quality.stable_as", "Stable-AS event QC", bound)
    return resolve_table(context, "quality.stable_as", "Stable-AS event QC", parameters)


def _quality_events(context: ProjectContext, parameters: dict[str, Any], events: pd.DataFrame) -> tuple[pd.DataFrame, AnalysisArtifact | None]:
    source = _quality_source(context, parameters)
    if source is None:
        return significance_columns(apply_stable_event_qc(events, parameters), parameters), None
    evaluated = pd.read_parquet(source.path)
    evaluated = evaluated[evaluated["comparison_id"].isin(selected_contrast_ids(context, parameters))].copy()
    if "measurement_eligible" not in evaluated:
        raise ValueError("Quality cache predates the separate measurement policy. Rebuild Quality explicitly.")
    if set(events["row_id"]) - set(evaluated["row_id"]):
        raise ValueError("Quality cache does not cover the selected source events. Rebuild Quality explicitly.")
    return evaluated, source


class ProjectAuditWorkflow(WorkflowBase):
    version = "2.2.0"
    definition = WorkflowDefinition(
        workflow_id="project.audit",
        workspace_id="project",
        label="Audit project design",
        question="Are inputs, samples, references and contrast directions explicit and usable?",
        method="Schema, path and experimental-unit validation without running an analysis.",
        output_summary="Input audit, sample sheet, contrast design and sample-count plot.",
        analysis_unit="input asset and biological sample",
        expected_artifacts=4,
        parameters=(ParameterDefinition("design_factors", "Declared design factors", "text", "condition", "Comma-separated declared categorical factors, at most three. No guessing from names."),
                    ParameterDefinition("design_references", "Factor references (JSON)", "text", "{}", "Explicit baseline per factor, e.g. {\"genotype\":\"WT\",\"DP\":\"noDP\"}. Otherwise alphabetical ordering is labeled explicitly."),
                    ParameterDefinition("include_batch", "Include declared batch", "boolean", False, "Rejects incomplete batch; confounded contrasts remain non-estimable."),
                    ParameterDefinition("pairing", "Biological-unit blocking", "select", "auto", "Auto blocks repeated biological units; technical samples are not independent replicates.", ("auto", "paired", "unpaired"))),
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
        from katana_next.services.readiness import experimental_design
        design = experimental_design(context.manifest)
        if context.manifest.samples:
            from katana_next.science.contrast_planner import plan_contrasts
            import json
            planned = plan_contrasts(context.manifest, [f.strip() for f in parameters.get("design_factors", "condition").split(",")],
                include_batch=bool(parameters.get("include_batch", False)), pairing=parameters.get("pairing", "auto"),
                reference_levels=json.loads(parameters.get("design_references", "{}")))
            rows = [{k: json.dumps(v, sort_keys=True) if isinstance(v,(dict,list)) else v for k,v in row.items()} for row in planned["contrasts"]]
            writer.table("Estimable design contrasts (not fitted results)", "planned_contrasts.parquet", pd.DataFrame(rows), entity_type="design_contrast",
                         metadata={"design": {k:v for k,v in planned.items() if k != "contrasts"}}, parents=(sample_artifact.artifact_id,))
        pair_review = pd.DataFrame(design["comparison_pairs"])
        if len(pair_review):
            writer.table("Shared controls and comparison dependence", "comparison_dependence.parquet", pair_review, entity_type="contrast_pair")
        counts = (
            samples.groupby("condition", dropna=False).biological_sample_id.nunique().rename("biological_samples").reset_index()
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
    version = "2.2.0"
    definition = WorkflowDefinition(
        workflow_id="quality.stable_as",
        workspace_id="quality",
        label="DAS read support and replicate consistency",
        question="Which events have adequate read support and reproducible replicate behavior?",
        method="Replicate read support and direction consistency; approximate PSI uncertainty is diagnostic only by default.",
        output_summary="Which events pass the read-count and replicate-consistency checks, and how many are significant; full event-level audit and filtering plot.",
        analysis_unit="rMATS event within one declared contrast",
        expected_artifacts=3,
        supports_contrast_selection=True,
        parameters=COMMON_THRESHOLDS + (
            ParameterDefinition("min_reads", "Custom reads per replicate", "number", 10.0, "Used only by Custom profile.", minimum=0, step=1),
            ParameterDefinition("min_supported_replicates", "Custom supported replicates", "number", 2, "Used only by Custom profile.", minimum=1, step=1),
            ParameterDefinition("min_direction_consistency", "Custom direction consistency", "number", 0.75, "Used only by Custom profile.", minimum=0, maximum=1, step=0.05),
            ParameterDefinition("apply_measurement_se_gate", "Experimental SE gate (off by default)", "boolean", False, "Uncalibrated approximate SE; do not treat as a literature-validated stability criterion."),
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
            ("Measurement eligible", evaluated["measurement_eligible"]),
            ("Repeat-supported difference", evaluated["stable_event"]),
            ("Stable + significant", evaluated["stable_significant"]),
        ]
        from katana_next.science.attrition import cumulative_attrition
        attrition, first_failure = cumulative_attrition(evaluated, stages)
        evaluated["first_exclusion_stage"] = first_failure
        writer = ArtifactWriter(
            context, self.definition.workflow_id, run_fingerprint, parameters, self.version,
            common_metadata={
                "question": self.definition.question,
                "method": self.definition.method,
                "analysis_unit": self.definition.analysis_unit,
                "interpretation": "measurement_eligible uses read/replicate/uncertainty gates without requiring a group effect. differential_repeat_supported records direction and leave-one-out consistency. stable_significant requires both plus caller FDR/effect thresholds.",
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
        supports_contrast_selection=True,
        parameters=(QUALITY_SOURCE, ParameterDefinition("max_pca_events", "Variable events for PCA", "number", 2000, "Highest-variance measurable events only.", minimum=50, maximum=20000, step=50)),
        required_asset_kinds=("rmats_result",),
        limitations="Pairwise rMATS outputs can repeat a sample/event estimate. Katana reports disagreement and uses the median rather than silently duplicating observations.",
    )

    def parent_artifacts(self, context: ProjectContext, parameters: dict[str, Any]) -> list[AnalysisArtifact]:
        source = _quality_source(context, parameters)
        return [source] if source else []

    def run(self, context: ProjectContext, job: JobContext, run_fingerprint: str, parameters: dict[str, Any]) -> list[str]:
        job.report(0.05, "Building deduplicated sample-level PSI matrix")
        data = load_bulk_events(context, parameters)
        evaluated, quality = _quality_events(context, parameters, data.events)
        eligible = set(evaluated.loc[evaluated["measurement_eligible"], "event_id"])
        data = replace(data, events=evaluated, psi_matrix=data.psi_matrix.loc[data.psi_matrix.index.isin(eligible)],
                       read_matrix=data.read_matrix.loc[data.read_matrix.index.isin(eligible)])
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
                "measurement_policy_source": quality.artifact_id if quality else "declared defaults; no cached Quality result",
            },
        )
        psi = data.psi_matrix.reset_index()
        psi["row_id"] = psi["event_id"]
        psi_artifact = writer.table("Sample-level PSI matrix", "sample_level_psi.parquet", psi, entity_type="splice_event", parents=(quality.artifact_id,) if quality else ())
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
        supports_contrast_selection=True,
        parameters=(QUALITY_SOURCE,),
        required_asset_kinds=("rmats_result",),
        limitations="The inferential p-value remains the imported rMATS result; Katana read/replicate flags are supporting diagnostics, not independent validation.",
    )

    def parent_artifacts(self, context: ProjectContext, parameters: dict[str, Any]) -> list[AnalysisArtifact]:
        source = _quality_source(context, parameters)
        return [source] if source else []

    def run(self, context: ProjectContext, job: JobContext, run_fingerprint: str, parameters: dict[str, Any]) -> list[str]:
        job.report(0.05, "Normalizing rMATS exact events")
        data = load_bulk_events(context, parameters)
        job.report(0.24, f"Normalized rMATS source cache: {data.source_cache_status}")
        events, quality = _quality_events(context, parameters, data.events)
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
                "measurement_policy_source": quality.artifact_id if quality else "declared defaults; no cached Quality result",
                "effective_evidence_parameters": quality.parameters if quality else parameters,
                "interpretation": "Use stable_significant for the strict event set. significant retains caller significance even when replicate-support QC fails.",
                "limitations": self.definition.limitations,
                "source_cache_status": data.source_cache_status,
            },
        )
        table = writer.table("Differential AS events", "differential_as_events.parquet", events, entity_type="splice_event", parents=(quality.artifact_id,) if quality else ())
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
        supports_contrast_selection=True,
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
    version = "2.2.0"
    definition = WorkflowDefinition(
        workflow_id="differential.dtu",
        workspace_id="differential",
        label="Transcript / exon usage evidence",
        question="Which genes and transcripts show DTU or exon-usage support?",
        method="Schema-normalized import of DTU, stageR or DEXSeq adjusted p-values.",
        output_summary="Transcript/exon evidence table and gene-level support summary.",
        analysis_unit="transcript or exon feature within a contrast",
        expected_artifacts=2,
        supports_contrast_selection=True,
        parameters=(ParameterDefinition("fdr", "DTU/DEXSeq FDR", "number", 0.05, "Adjusted p-value cutoff.", minimum=0, maximum=1, step=0.01),),
        optional_asset_kinds=("dtu_result", "dexseq_result"),
    )

    def run(self, context: ProjectContext, job: JobContext, run_fingerprint: str, parameters: dict[str, Any]) -> list[str]:
        contrast_ids = selected_contrast_ids(context, parameters)
        frame = load_dtu_results(context, contrast_ids)
        if frame.empty:
            raise ValueError("No compatible DTU or DEXSeq result is available for selected contrasts.")
        from katana_next.science.isoform import transcript_gene_map
        from katana_next.science.enrichment import fill_gene_symbols
        if context.manifest.annotation.gtf:
            annotation = transcript_gene_map(context.manifest.annotation.gtf,
                context.manifest.output_root / ".katana" / "report_gene_annotation.parquet")
            frame = fill_gene_symbols(frame, annotation)
        frame["feature_type"] = np.where(frame.source_kind.eq("dexseq_result"), "exon_bin", "transcript")
        frame["significant"] = frame["dtu_fdr"].le(float(parameters.get("fdr", 0.05)))
        summary = frame.groupby(["comparison_id", "gene_id", "source_kind"], dropna=False).agg(best_fdr=("dtu_fdr", "min"), tested_features=("row_id", "size"), significant_features=("significant", "sum")).reset_index()
        summary["row_id"] = summary["comparison_id"] + "::" + summary["gene_id"] + "::" + summary["source_kind"]
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version, common_metadata={"question": self.definition.question, "method": self.definition.method, "analysis_unit": self.definition.analysis_unit, "interpretation": "DTU/DEXSeq supports transcript or exon usage. It does not replace event-level DAS."})
        detail = writer.table("DTU and exon-usage evidence", "dtu_features.parquet", frame, entity_type="transcript_feature")
        writer.table("Gene-level usage support by source", "dtu_gene_summary.parquet", summary, entity_type="gene", parents=(detail.artifact_id,))
        return writer.finish()


def _latest_table(context: ProjectContext, analysis_id: str, name: str | None = None,
                  parameters: dict[str, Any] | None = None) -> AnalysisArtifact | None:
    return resolve_table(context, analysis_id, name, parameters)


class CrossComparisonWorkflow(WorkflowBase):
    version = "2.1.1"
    definition = WorkflowDefinition(
        workflow_id="differential.cross_comparison",
        workspace_id="differential",
        label="Exact cross-comparison AS",
        question="How much of each strict AS signature is shared, directionally concordant or contrast-specific?",
        method="Exact normalized event identity, declared tested-event universe, hypergeometric overlap, direction concordance and effect correlation.",
        output_summary="Pair summary, event membership, overlap profile and similarity heatmap.",
        analysis_unit="pair of contrasts and exact splice event",
        expected_artifacts=6,
        supports_contrast_selection=True,
        parameters=(
            ParameterDefinition("membership_mode", "Intersection membership", "select", "exact_event", "Exact event ignores sign for membership. Event and direction requires both exact identity and equal nonzero standardized dPSI sign; zero-effect events are excluded in that mode.", ("exact_event", "event_and_direction")),
            ParameterDefinition("event_set", "Significant event set", "select", "stable_significant", "Both require caller significance; strict also uses the saved replicate QC policy.", ("stable_significant", "significant")),
            ParameterDefinition("min_reads", "Minimum IJC + SJC per replicate", "number", 0, "Optional extra read gate in each arm; 0 keeps the saved event set unchanged.", minimum=0, step=1),
            ParameterDefinition("min_supported_replicates", "Supported replicates in each arm", "number", 2, "Required in numerator and denominator when the extra read gate is enabled.", minimum=1, step=1),
        ),
        limitations="Same-gene and nearby-coordinate events are not silently treated as the same event.",
    )

    def parent_artifacts(self, context: ProjectContext, parameters: dict[str, Any]) -> list[AnalysisArtifact]:
        artifact = _latest_table(context, "differential.as", "Differential AS events", parameters)
        return [artifact] if artifact else []

    def input_paths(self, context: ProjectContext, parameters: dict[str, Any]) -> list[Path]:
        return [artifact.path for artifact in self.parent_artifacts(context, parameters)]

    def run(self, context: ProjectContext, job: JobContext, run_fingerprint: str, parameters: dict[str, Any]) -> list[str]:
        parent = _latest_table(context, "differential.as", "Differential AS events", parameters)
        if parent is None:
            raise ValueError("Run Differential alternative splicing first. Cross-comparison never starts it implicitly.")
        events = pd.read_parquet(parent.path)
        selected = parameters.get("contrast_ids")
        if selected:
            events = events[events["comparison_id"].isin(selected)]
        expected = set(selected or [c.contrast_id for c in context.manifest.contrasts])
        missing = expected - set(events["comparison_id"])
        if missing:
            raise ValueError(f"Source DAS result is missing comparisons: {', '.join(sorted(missing))}. Run them explicitly first.")
        if events["comparison_id"].nunique() < 2:
            raise ValueError("Select at least two contrasts for cross-comparison.")
        event_set = str(parameters.get("event_set", "stable_significant"))
        if event_set not in {"stable_significant", "significant"}:
            raise ValueError("Choose significant or stable_significant events.")
        events = events.copy()
        events["stable_significant"] = events[event_set].fillna(False).astype(bool)
        if float(parameters.get("min_reads", 0)) > 0:
            from katana_next.science.bulk import _arm_read_totals
            for arm in ("numerator", "denominator"):
                supported = events.apply(lambda row: sum(value >= float(parameters["min_reads"])
                                                          for value in _arm_read_totals(row, arm)), axis=1)
                events["stable_significant"] &= supported.ge(int(parameters.get("min_supported_replicates", 2)))
        job.report(0.25, "Matching exact event structures")
        summary, membership = cross_comparison_summary(events)
        mode = parameters.get("membership_mode", "exact_event")
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version, common_metadata={"question": self.definition.question, "method": self.definition.method, "analysis_unit": self.definition.analysis_unit, "interpretation": "Overlap, direction agreement and effect correlation answer different questions and are kept separate.", "limitations": self.definition.limitations})
        summary_artifact = writer.table("Cross-comparison exact-event summary", "cross_comparison_summary.parquet", summary, entity_type="contrast_pair", parents=(parent.artifact_id,))
        writer.table("Cross-comparison event membership", "cross_comparison_events.parquet", membership, entity_type="splice_event", parents=(summary_artifact.artifact_id, parent.artifact_id))
        from katana_next.workflows.intersection_plots import exact_intersections, render_intersections, event_membership_sets
        labels = {contrast.contrast_id: contrast.label for contrast in context.manifest.contrasts}
        sets, set_members = event_membership_sets(events, mode)
        writer.table("Intersection member identities and directions", "intersection_members.parquet", set_members, entity_type="splice_event", parents=(parent.artifact_id,))
        intersections = writer.table("All exact set intersections", "all_exact_intersections.parquet", exact_intersections(sets), entity_type="intersection_summary", parents=(summary_artifact.artifact_id,))
        exports = render_intersections(sets, labels, writer.directory / "plots", mode=mode)
        writer.existing_file("Selected-contrast Venn" if len(sets) <= 3 else "Selected-contrast UpSet", Path(exports["png"]), entity_type="intersection_plot",
                             metadata={"mime_type": "image/png", "exports": exports, "identity_rule": mode, "zero_direction": "excluded in event_and_direction; retained in exact_event"}, parents=(intersections.artifact_id,))
        long_values = []
        shared_label = "Shared event + direction" if mode == "event_and_direction" else "Shared exact"
        for row in summary.to_dict("records"):
            a, b = sets[row['comparison_a']], sets[row['comparison_b']]
            for category, count in ((shared_label, len(a & b)), ("A-specific", len(a - b)), ("B-specific", len(b - a))):
                long_values.append({"pair": f"{labels.get(row['comparison_a'], row['comparison_a'])}\ncompared with\n{labels.get(row['comparison_b'], row['comparison_b'])}", "category": category, "events": count, "comparison_a": row['comparison_a'], "comparison_b": row['comparison_b']})
        writer.plot(
            "Exact-event overlap profile",
            "cross_comparison_overlap.vl.json",
            _plot(long_values, {"type": "bar", "cornerRadiusEnd": 3}, {"x": {"field": "events", "type": "quantitative", "axis": {"title": "Selected AS set members"}}, "y": {"field": "pair", "type": "nominal", "axis": {"title": None, "labelLimit": 320}}, "color": {"field": "category", "type": "nominal", "scale": {"domain": [shared_label, "A-specific", "B-specific"], "range": [TEAL, BLUE, ORANGE]}}, "tooltip": [{"field": "pair"}, {"field": "category"}, {"field": "events"}]}, title=f"AS overlap: {mode.replace('_', ' ')}", description="All plot counts use the selected membership rule. A/B identify the first/second named contrast; specificity is relative to this rule, not biological absence.", height=max(300, 65 * len(summary))),
            entity_type="contrast_pair_summary",
            parents=(summary_artifact.artifact_id,),
        )
        matrix_values = []
        for row in summary.to_dict("records"):
            label_a, label_b = labels.get(row["comparison_a"], row["comparison_a"]), labels.get(row["comparison_b"], row["comparison_b"])
            sa, sb = sets[row["comparison_a"]], sets[row["comparison_b"]]
            display_jaccard = len(sa & sb) / len(sa | sb) if sa | sb else None
            matrix_values.extend([
                {"comparison_a": label_a, "comparison_b": label_b, "jaccard": display_jaccard, "direction_concordance": row["direction_concordance"], "membership_mode": mode},
                {"comparison_a": label_b, "comparison_b": label_a, "jaccard": display_jaccard, "direction_concordance": row["direction_concordance"], "membership_mode": mode},
            ])
        for name, event_ids in sets.items():
            matrix_values.append({"comparison_a": labels.get(name, name), "comparison_b": labels.get(name, name), "jaccard": 1.0 if event_ids else None, "direction_concordance": 1.0 if event_ids else None})
        writer.plot(
            "AS similarity matrix",
            "cross_comparison_similarity.vl.json",
            _plot(matrix_values, "rect", {"x": {"field": "comparison_a", "type": "nominal", "axis": {"title": None, "labelAngle": -45, "labelLimit": 220}}, "y": {"field": "comparison_b", "type": "nominal", "axis": {"title": None, "labelLimit": 220}}, "color": {"field": "jaccard", "type": "quantitative", "scale": {"domain": [0, 1], "scheme": "tealblues"}, "legend": {"title": "Jaccard"}}, "tooltip": [{"field": "comparison_a"}, {"field": "comparison_b"}, {"field": "jaccard", "format": ".3f"}, {"field": "direction_concordance", "format": ".1%"}]}, title="Strict AS signature similarity", description="Jaccard quantifies set similarity; tooltip separately reports same-sign dPSI among shared exact events.", height=500),
            entity_type="contrast_pair_summary",
            parents=(summary_artifact.artifact_id,),
        )
        writer.common_metadata["membership_mode"] = mode
        job.report(0.92, "Exact cross-comparison complete")
        return writer.finish()


class CandidateEvidenceWorkflow(WorkflowBase):
    version = "2.1.2"
    definition = WorkflowDefinition(
        workflow_id="prioritize.evidence",
        workspace_id="prioritize",
        label="Transparent candidate evidence",
        question="Which genes have AS, DEG, DTU or DEXSeq support, and which are shared between comparisons?",
        method="One gene per contrast; ranking_score counts AS/DEG/DTU/DEXSeq support (0-4), not independent evidence or a p-value. Pair membership uses at least one supported analysis. Direction changes compare significant DEG or dominant stable AS signs at gene level, not necessarily the same event.",
        output_summary="Candidate matrix, contrast-resolved heatmap and all-pair gene membership/direction table and summary.",
        analysis_unit="gene within a contrast",
        expected_artifacts=4,
        supports_contrast_selection=True,
        parameters=(ParameterDefinition("fdr", "Evidence FDR", "number", 0.05, "Common display threshold for imported evidence.", minimum=0, maximum=1, step=0.01), ParameterDefinition("deg_min_abs_log2fc", "DEG minimum |log2FC|", "number", 1.0, "DEG effect-size gate.", minimum=0, step=0.1)),
    )

    def parent_artifacts(self, context: ProjectContext, parameters: dict[str, Any]) -> list[AnalysisArtifact]:
        artifacts = []
        for analysis_id, name in (("differential.as", "Differential AS events"), ("differential.deg", "Differentially expressed genes"), ("differential.dtu", "DTU and exon-usage evidence")):
            artifact = _latest_table(context, analysis_id, name, parameters)
            if artifact:
                artifacts.append(artifact)
        return artifacts

    def input_paths(self, context: ProjectContext, parameters: dict[str, Any]) -> list[Path]:
        return [artifact.path for artifact in self.parent_artifacts(context, parameters)]

    def run(self, context: ProjectContext, job: JobContext, run_fingerprint: str, parameters: dict[str, Any]) -> list[str]:
        as_artifact = _latest_table(context, "differential.as", "Differential AS events", parameters)
        deg_artifact = _latest_table(context, "differential.deg", "Differentially expressed genes", parameters)
        dtu_artifact = _latest_table(context, "differential.dtu", "DTU and exon-usage evidence", parameters)
        if not any((as_artifact, deg_artifact, dtu_artifact)):
            raise ValueError("Run at least one AS, DEG or DTU evidence workflow explicitly first.")
        events = pd.read_parquet(as_artifact.path) if as_artifact else pd.DataFrame()
        deg = pd.read_parquet(deg_artifact.path) if deg_artifact else pd.DataFrame()
        dtu = pd.read_parquet(dtu_artifact.path) if dtu_artifact else pd.DataFrame()
        selected = parameters.get("contrast_ids")
        if selected:
            if not events.empty:
                events = events[events["comparison_id"].isin(selected)]
            if not deg.empty:
                deg = deg[deg["comparison_id"].isin(selected)]
            if not dtu.empty:
                dtu = dtu[dtu["comparison_id"].isin(selected)]
        job.report(0.35, "Joining transparent evidence dimensions")
        candidates = build_candidate_evidence(events, deg, dtu, parameters)
        labels = {c.contrast_id: c.label for c in context.manifest.contrasts}
        candidates["display_row"] = candidate_display_labels(candidates, labels)
        candidates = candidates[[c for c in candidates if c != "ranking_score"] + ["ranking_score"]]
        long = candidates.melt(id_vars=["row_id", "comparison_id", "gene_symbol", "gene_id", "display_row"], value_vars=["as_supported", "deg_supported", "dtu_supported", "dexseq_supported"], var_name="evidence", value_name="supported")
        top_ids = set(candidates.head(60)["row_id"])
        long = long[long["row_id"].isin(top_ids)]
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version, common_metadata={"question": self.definition.question, "method": self.definition.method, "analysis_unit": self.definition.analysis_unit, "interpretation": "Supporting analyses may be correlated; their count is prioritization, not independent statistical confirmation. Pair count 1 = only one side, 2 = shared. Dominant-AS signs can concern different events; gene-level direction changed is not exact-event reversal.", "missing_evidence": "Missing DEG, DTU or DEXSeq stays unavailable, distinct from a measured nonsignificant result."})
        table = writer.table("Candidate evidence matrix", "candidate_evidence.parquet", candidates, entity_type="gene", parents=tuple(artifact.artifact_id for artifact in self.parent_artifacts(context, parameters)))
        writer.plot(
            "Candidate evidence heatmap",
            "candidate_evidence.vl.json",
            _plot(long.to_dict("records"), {"type": "rect", "cornerRadius": 2}, {"x": {"field": "evidence", "type": "nominal", "axis": {"title": None, "labelExpr": "replace(datum.label, '_supported', '')"}}, "y": {"field": "display_row", "type": "nominal", "sort": None, "axis": {"title": None, "labelLimit": 380}}, "color": {"field": "supported", "type": "nominal", "scale": {"domain": [False, True], "range": ["#e2e8f0", "#0f766e"]}, "legend": {"title": "Threshold support"}}, "tooltip": [{"field": "gene_symbol"}, {"field": "gene_id"}, {"field": "comparison_id"}, {"field": "evidence"}, {"field": "supported"}]}, title="Top gene-by-comparison evidence", description="Up to 60 distinct gene/contrast rows; full evidence and availability flags in the table. Grey means no threshold-passing support, not proof of no effect.", height=max(360, 18 * max(1, long["display_row"].nunique()))),
            entity_type="gene_summary",
            parents=(table.artifact_id,),
        )
        comparison_ids = selected_contrast_ids(context, parameters)
        pairs = pairwise_candidate_membership(candidates, comparison_ids)
        pairs["comparison_A_label"] = pairs.comparison_A.map(labels)
        pairs["comparison_B_label"] = pairs.comparison_B.map(labels)
        pair_table = writer.table("Pairwise candidate gene membership and patterns", "pairwise_candidate_genes.parquet", pairs,
                                  entity_type="gene_pair", parents=(table.artifact_id,))
        if not pairs.empty:
            counts = candidate_pair_plot_counts(pairs, labels)
            counts["membership_label"] = np.select([counts.membership.eq("A_specific"), counts.membership.eq("B_specific")],
                [counts.comparison_A_label + " only", counts.comparison_B_label + " only"], default="Shared candidates")
            writer.plot("Candidate gene membership by comparison pair", "pairwise_candidate_genes.vl.json",
                _plot(counts.to_dict("records"), "bar", {"y": {"field": "pair", "type": "nominal", "axis": {"labelLimit": 420, "title": None}}, "x": {"field": "genes", "type": "quantitative", "title": "Unique candidate genes"}, "color": {"field": "membership_label", "type": "nominal", "legend": {"title": "Membership"}}, "tooltip": [{"field": "comparison_A_label"}, {"field": "comparison_B_label"}, {"field": "membership_label"}, {"field": "genes"}]}, title="Gene-level candidate overlap", description="At least one supported AS/DEG/DTU/DEXSeq analysis. This is gene membership, not event-level interaction.", height=max(300, counts.pair.nunique()*50)),
                entity_type="gene_pair_summary", parents=(pair_table.artifact_id,))
        return writer.finish()
