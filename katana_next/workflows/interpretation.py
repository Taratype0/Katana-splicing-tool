from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from katana_next.domain.models import AnalysisArtifact
from katana_next.infrastructure.jobs import JobContext
from katana_next.science.enrichment import overrepresentation, read_gene_sets
from katana_next.science.isoform import (
    dominant_isoform_switches,
    isoform_contrast_summary,
    load_quant_support,
)
from katana_next.services.project_service import ProjectContext
from katana_next.workflows.base import (
    ArtifactWriter,
    ParameterDefinition,
    WorkflowBase,
    WorkflowDefinition,
)
from katana_next.workflows.bulk import TEAL, _latest_table, _plot


def _selection_genes(context: ProjectContext, parameters: dict[str, Any]) -> tuple[set[str], AnalysisArtifact]:
    source = _latest_table(context, "prioritize.evidence", "Candidate evidence matrix")
    if source is None:
        raise ValueError("Run Transparent candidate evidence first.")
    frame = pd.read_parquet(source.path)
    selection_id = str(parameters.get("selection_id", "")).strip()
    if selection_id:
        selection = next((item for item in context.catalog.list_selections(context.manifest.project_id) if item.selection_id == selection_id), None)
        if selection is None:
            raise ValueError("The selected gene set no longer exists.")
        frame = frame[frame["row_id"].astype(str).isin(selection.entity_ids)]
    else:
        frame = frame[frame["evidence_dimensions"].fillna(0).gt(0)]
    symbols = set(frame["gene_symbol"].dropna().astype(str).str.strip()) - {""}
    identifiers = set(frame["gene_id"].dropna().astype(str).str.strip()) - {""}
    return symbols | identifiers, source


class GeneSetEnrichmentWorkflow(WorkflowBase):
    version = "2.0.0"
    definition = WorkflowDefinition(
        workflow_id="prioritize.enrichment",
        workspace_id="prioritize",
        label="GO / KEGG / custom enrichment",
        question="Which declared pathways are overrepresented among selected candidate genes?",
        method="Hypergeometric overrepresentation with BH FDR against the tested candidate-gene background.",
        output_summary="Enrichment table and dot plot with overlapping genes.",
        analysis_unit="gene set",
        expected_artifacts=2,
        parameters=(ParameterDefinition("selection_id", "Saved gene selection", "selection", "", "Optional saved candidate-gene set; otherwise all evidence-supported candidates are used."), ParameterDefinition("minimum_overlap", "Minimum overlapping genes", "number", 2, "Gene sets below this overlap are not tested.", minimum=1, step=1)),
        required_asset_kinds=("gene_sets_gmt",),
        limitations="Enrichment depends on the declared tested background and supports interpretation, not causality.",
    )

    def parent_artifacts(self, context: ProjectContext, parameters: dict[str, Any]) -> list[AnalysisArtifact]:
        del parameters
        artifact = _latest_table(context, "prioritize.evidence", "Candidate evidence matrix")
        return [artifact] if artifact else []

    def input_paths(self, context: ProjectContext, parameters: dict[str, Any]) -> list[Path]:
        return [*super().input_paths(context, parameters), *(artifact.path for artifact in self.parent_artifacts(context, parameters))]

    def run(self, context: ProjectContext, job: JobContext, run_fingerprint: str, parameters: dict[str, Any]) -> list[str]:
        genes, source = _selection_genes(context, parameters)
        candidates = pd.read_parquet(source.path)
        background = set(candidates["gene_symbol"].dropna().astype(str)) | set(candidates["gene_id"].dropna().astype(str))
        asset = next(item for item in context.manifest.assets if item.kind == "gene_sets_gmt" and item.path.is_file())
        gene_sets = read_gene_sets(asset.path)
        result = overrepresentation(genes, background, gene_sets, minimum_overlap=int(parameters.get("minimum_overlap", 2)))
        if result.empty:
            result = pd.DataFrame(columns=["row_id", "term", "selected_genes", "background_genes", "term_genes_in_background", "overlap_genes", "expected_overlap", "fold_enrichment", "pvalue", "fdr", "genes"])
        result["minus_log10_fdr"] = -np.log10(pd.to_numeric(result.get("fdr"), errors="coerce").clip(lower=np.finfo(float).tiny)) if "fdr" in result else np.nan
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version, common_metadata={"question": self.definition.question, "method": self.definition.method, "analysis_unit": self.definition.analysis_unit, "interpretation": "Terms with low FDR and coherent overlapping genes prioritize biological themes; they do not establish mechanism.", "limitations": self.definition.limitations, "background_definition": "all genes in the candidate evidence artifact"})
        table = writer.table("Gene-set enrichment", "gene_set_enrichment.parquet", result, entity_type="gene_set", parents=(source.artifact_id,))
        top = result.head(30)
        writer.plot("Gene-set enrichment", "gene_set_enrichment.vl.json", _plot(top.to_dict("records"), {"type": "circle", "opacity": 0.85, "stroke": "white"}, {"x": {"field": "fold_enrichment", "type": "quantitative", "axis": {"title": "Fold enrichment"}}, "y": {"field": "term", "type": "nominal", "sort": {"field": "fdr"}, "axis": {"title": None, "labelLimit": 340}}, "size": {"field": "overlap_genes", "type": "quantitative"}, "color": {"field": "minus_log10_fdr", "type": "quantitative", "scale": {"scheme": "tealblues"}, "legend": {"title": "-log10 FDR"}}, "tooltip": [{"field": "term"}, {"field": "overlap_genes"}, {"field": "fold_enrichment", "format": ".2f"}, {"field": "fdr", "format": ".2e"}, {"field": "genes"}]}, title="Functional themes in selected candidates", description="Hypergeometric overrepresentation against the tested candidate background.", height=max(360, min(900, 24 * len(top)))), entity_type="gene_set_summary", parents=(table.artifact_id,))
        return writer.finish()


class RbpSupportWorkflow(WorkflowBase):
    version = "2.0.0"
    definition = WorkflowDefinition(
        workflow_id="prioritize.rbp_support",
        workspace_id="prioritize",
        label="RBP and splicing-factor support",
        question="Which candidate genes overlap declared RBP or splicing-factor target sets?",
        method="Matched-background overlap and hypergeometric enrichment from user-versioned or public gene sets.",
        output_summary="RBP enrichment, matched candidate genes and cross-contrast support matrix.",
        analysis_unit="RBP target set and candidate gene",
        expected_artifacts=3,
        parameters=(ParameterDefinition("selection_id", "Saved gene selection", "selection", "", "Optional candidate-gene subset."), ParameterDefinition("minimum_overlap", "Minimum overlapping genes", "number", 1, "Minimum RBP-set overlap.", minimum=1, step=1)),
        required_asset_kinds=("rbp_gene_sets",),
        limitations="RBP overlap prioritizes regulators or targets and is never reported as causal control.",
    )

    def parent_artifacts(self, context: ProjectContext, parameters: dict[str, Any]) -> list[AnalysisArtifact]:
        del parameters
        artifact = _latest_table(context, "prioritize.evidence", "Candidate evidence matrix")
        return [artifact] if artifact else []

    def input_paths(self, context: ProjectContext, parameters: dict[str, Any]) -> list[Path]:
        return [*super().input_paths(context, parameters), *(artifact.path for artifact in self.parent_artifacts(context, parameters))]

    def run(self, context: ProjectContext, job: JobContext, run_fingerprint: str, parameters: dict[str, Any]) -> list[str]:
        genes, source = _selection_genes(context, parameters)
        candidates = pd.read_parquet(source.path)
        background = set(candidates["gene_symbol"].dropna().astype(str)) | set(candidates["gene_id"].dropna().astype(str))
        asset = next(item for item in context.manifest.assets if item.kind == "rbp_gene_sets" and item.path.is_file())
        gene_sets = read_gene_sets(asset.path)
        enrichment = overrepresentation(genes, background, gene_sets, minimum_overlap=int(parameters.get("minimum_overlap", 1)))
        if enrichment.empty:
            enrichment = pd.DataFrame(
                columns=[
                    "row_id",
                    "term",
                    "selected_genes",
                    "background_genes",
                    "term_genes_in_background",
                    "overlap_genes",
                    "expected_overlap",
                    "fold_enrichment",
                    "pvalue",
                    "fdr",
                    "genes",
                ]
            )
        membership_rows = []
        for rbp, members in gene_sets.items():
            for gene in sorted(genes & members):
                membership_rows.append({"row_id": f"{rbp}::{gene}", "rbp_or_set": rbp, "gene": gene})
        membership = pd.DataFrame(membership_rows, columns=["row_id", "rbp_or_set", "gene"])
        target_genes = set().union(*gene_sets.values()) if gene_sets else set()
        selected_candidates = candidates[
            candidates["gene_symbol"].isin(genes) | candidates["gene_id"].isin(genes)
        ]
        matched = selected_candidates[
            selected_candidates["gene_symbol"].isin(target_genes)
            | selected_candidates["gene_id"].isin(target_genes)
        ].copy()
        matched["supporting_rbp_sets"] = matched.apply(
            lambda row: "; ".join(
                sorted(
                    name
                    for name, members in gene_sets.items()
                    if str(row.get("gene_symbol", "")) in members
                    or str(row.get("gene_id", "")) in members
                )
            ),
            axis=1,
        )
        matched["is_rbp_or_target"] = matched["supporting_rbp_sets"].ne("")
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version, common_metadata={"question": self.definition.question, "method": self.definition.method, "analysis_unit": self.definition.analysis_unit, "interpretation": "Use these results for mechanism prioritization only. Motif/eCLIP and expression association are separate supporting layers.", "limitations": self.definition.limitations})
        enrichment_artifact = writer.table("RBP target-set enrichment", "rbp_enrichment.parquet", enrichment, entity_type="rbp_gene_set", parents=(source.artifact_id,))
        writer.table("RBP matched genes", "rbp_matched_genes.parquet", membership, entity_type="gene", parents=(enrichment_artifact.artifact_id, source.artifact_id))
        writer.table("RBP candidate support matrix", "rbp_candidate_support.parquet", matched, entity_type="gene", parents=(source.artifact_id,))
        return writer.finish()


class IsoformSupportWorkflow(WorkflowBase):
    version = "2.0.0"
    definition = WorkflowDefinition(
        workflow_id="prioritize.isoform_support",
        workspace_id="prioritize",
        label="Isoform expression and switch support",
        question="Do candidate genes have transcript-level expression and dominant-isoform changes?",
        method="Salmon quant.sf aggregation through the declared GTF plus separate DTU evidence.",
        output_summary="Transcript TPM, gene TPM, contrast delta TPM and descriptive dominant-switch candidates.",
        analysis_unit="transcript within a biological sample",
        expected_artifacts=5,
        limitations="A dominant transcript change is descriptive unless supported by an independent DTU/switch test.",
    )

    def input_paths(self, context: ProjectContext, parameters: dict[str, Any]) -> list[Path]:
        del parameters
        paths = [sample.assets["quant_sf"] for sample in context.manifest.samples if "quant_sf" in sample.assets]
        if context.manifest.annotation.gtf:
            paths.append(context.manifest.annotation.gtf)
        return paths

    def run(self, context: ProjectContext, job: JobContext, run_fingerprint: str, parameters: dict[str, Any]) -> list[str]:
        cache = context.manifest.output_root / ".katana" / "transcript_gene_map.parquet"
        job.report(0.08, "Resolving transcript-to-gene annotation")
        transcript, gene = load_quant_support(context, cache)
        contrasts = isoform_contrast_summary(context, transcript)
        switches = dominant_isoform_switches(context, transcript)
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version, common_metadata={"question": self.definition.question, "method": self.definition.method, "analysis_unit": self.definition.analysis_unit, "interpretation": "TPM supports abundance and visualization. DTU significance must remain a separate statistical result.", "limitations": self.definition.limitations})
        transcript_artifact = writer.table("Sample transcript abundance", "transcript_tpm.parquet", transcript, entity_type="transcript")
        writer.table("Sample gene abundance", "gene_tpm.parquet", gene, entity_type="gene", parents=(transcript_artifact.artifact_id,))
        contrast_artifact = writer.table("Isoform delta TPM by contrast", "isoform_delta_tpm.parquet", contrasts, entity_type="transcript", parents=(transcript_artifact.artifact_id,))
        switch_artifact = writer.table("Descriptive dominant isoform switches", "dominant_isoform_switches.parquet", switches, entity_type="gene", parents=(transcript_artifact.artifact_id,))
        plot_values = contrasts.assign(abs_delta_tpm=lambda frame: frame["delta_tpm"].abs()).sort_values("abs_delta_tpm", ascending=False).groupby("comparison_id").head(25)
        writer.plot("Isoform abundance shifts", "isoform_delta_tpm.vl.json", _plot(plot_values.to_dict("records"), {"type": "bar", "cornerRadiusEnd": 2}, {"x": {"field": "delta_tpm", "type": "quantitative", "axis": {"title": "Mean TPM difference (numerator minus denominator)"}}, "y": {"field": "transcript_id", "type": "nominal", "sort": {"field": "delta_tpm", "order": "descending"}, "axis": {"title": None, "labelLimit": 160}}, "color": {"condition": {"test": "datum.delta_tpm >= 0", "value": TEAL}, "value": "#c2410c"}, "facet": {"field": "comparison_id", "columns": 2}, "tooltip": [{"field": "comparison_id"}, {"field": "gene_symbol"}, {"field": "transcript_id"}, {"field": "delta_tpm", "format": ".2f"}]}, title="Largest transcript abundance shifts", description="Descriptive Salmon TPM support; use DTU for significance.", height=500), entity_type="transcript_summary", parents=(contrast_artifact.artifact_id, switch_artifact.artifact_id))
        return writer.finish()
