from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from katana_next.domain.models import AnalysisArtifact
from katana_next.infrastructure.jobs import JobContext
from katana_next.science.enrichment import canonical_gene_mapping, map_gene_sets, overrepresentation, declared_gene_sets, library_overrepresentation, fill_gene_symbols, enrichment_background
from katana_next.science.isoform import (
    dominant_isoform_switches,
    isoform_contrast_summary,
    load_quant_support,
    transcript_gene_map,
)
from katana_next.science.bulk import usage_gene_support
from katana_next.services.project_service import ProjectContext
from katana_next.workflows.base import (
    ArtifactWriter,
    ParameterDefinition,
    WorkflowBase,
    WorkflowDefinition,
)
from katana_next.workflows.bulk import TEAL, _latest_table, _plot
from katana_next.workflows.sources import selected_source, selection_rows


def annotated_gene_source(context, source):
    frame = pd.read_parquet(source.path)
    from katana_next.science.gene_filter import apply_gene_filter
    frame = apply_gene_filter(context, frame)
    gtf = context.manifest.annotation.gtf
    if gtf is not None and gtf.is_file() and "gene_id" in frame:
        annotation = transcript_gene_map(gtf, context.manifest.output_root / ".katana" / "transcript_gene_map.parquet")
        frame = fill_gene_symbols(frame, annotation)
    return frame


def _selection_frame(context: ProjectContext, parameters: dict[str, Any]) -> tuple[pd.DataFrame, AnalysisArtifact]:
    source = selected_source(context, parameters) or _latest_table(context, "prioritize.evidence", "Candidate evidence matrix", parameters)
    if source is None:
        raise ValueError("Run Transparent candidate evidence first.")
    frame = annotated_gene_source(context, source)
    if not {"gene_id", "gene_symbol"}.intersection(frame.columns):
        raise ValueError("This selection has no gene annotation. Select annotated gene/event/module-member rows.")
    canonical, _, _ = canonical_gene_mapping(frame)
    frame = frame.assign(canonical_gene_id=canonical)
    if parameters.get("selection_id"):
        frame = selection_rows(context, parameters, source, frame)
    else:
        frame = frame[frame["evidence_dimensions"].fillna(0).gt(0)]
    return frame, source


def _selection_genes(context: ProjectContext, parameters: dict[str, Any]) -> tuple[set[str], AnalysisArtifact]:
    frame, source = _selection_frame(context, parameters)
    return set(frame["canonical_gene_id"]) - {""}, source


def _comparison_selection_mask(candidates: pd.DataFrame, selected: pd.DataFrame) -> pd.Series:
    keys = ["comparison_id", "canonical_gene_id"] if "comparison_id" in candidates and "comparison_id" in selected else ["canonical_gene_id"]
    chosen = pd.MultiIndex.from_frame(selected[keys].drop_duplicates())
    return pd.Series(pd.MultiIndex.from_frame(candidates[keys]).isin(chosen), index=candidates.index)


def _candidate_background(context, parameters, source, selected):
    candidates = annotated_gene_source(context, source)
    if parameters.get("contrast_ids") and "comparison_id" in candidates:
        candidates = candidates[candidates.comparison_id.isin(parameters["contrast_ids"])].copy()
    if parameters.get("selection_id") and "comparison_id" in selected and "comparison_id" in candidates:
        candidates = candidates[candidates.comparison_id.isin(selected.comparison_id)].copy()
    if candidates.empty:
        raise ValueError("The saved selection does not cover the requested comparison; choose a matching selection.")
    candidates, audit = enrichment_background(candidates)
    canonical, aliases, mapping_audit = canonical_gene_mapping(candidates)
    candidates = candidates.assign(canonical_gene_id=canonical)
    candidates["selected_input"] = _comparison_selection_mask(candidates, selected)
    if not parameters.get("selection_id") and "evidence_dimensions" in candidates:
        candidates["selected_input"] &= candidates.evidence_dimensions.fillna(0).gt(0)
    return candidates, audit, aliases, mapping_audit


class GeneSetEnrichmentWorkflow(WorkflowBase):
    version = "2.3.1"
    definition = WorkflowDefinition(
        workflow_id="prioritize.enrichment",
        workspace_id="prioritize",
        label="GO / KEGG / custom enrichment",
        question="Which declared pathways are overrepresented among selected candidate genes?",
        method="Hypergeometric overrepresentation with BH FDR against genes tested in the contrast and annotated in each declared library.",
        output_summary="All enrichment results and comparison-specific bubble plots: gene ratio, contributor count, FDR and GO category.",
        analysis_unit="gene set",
        supports_contrast_selection=True,
        expected_artifacts=4,
        parameters=(
            ParameterDefinition("selection_id", "Saved gene selection", "selection", "", "Optional saved annotated gene/event/module-member selection; otherwise evidence-supported candidates are used."),
            ParameterDefinition("minimum_overlap", "Minimum overlapping genes (display only)", "number", 2, "Applied after BH correction; zero-overlap eligible gene sets stay in the testing family.", minimum=1, step=1),
            ParameterDefinition("ontology_category", "Show GO / pathway category", "select", "All", "Display filter only; full results and the declared FDR family are unchanged. GO classification requires a versioned OBO resource.", ("All", "GO BP", "GO MF", "GO CC", "KEGG", "Unclassified")),
            ParameterDefinition("bubble_x", "Bubble plot horizontal axis", "select", "gene_ratio", "Gene ratio is contributing genes divided by selected genes mapped to this library.", ("gene_ratio", "fold_enrichment")),
            ParameterDefinition("top_terms", "Terms shown per comparison", "number", 12, "Top terms by FDR for display; all terms remain downloadable.", minimum=1, maximum=40, step=1),
        ),
        required_asset_kinds=("gene_sets_gmt",),
        optional_asset_kinds=("go_ontology_obo",),
        limitations="Enrichment depends on the declared tested background and supports interpretation, not causality.",
    )

    def parent_artifacts(self, context: ProjectContext, parameters: dict[str, Any]) -> list[AnalysisArtifact]:
        artifact = selected_source(context, parameters) or _latest_table(context, "prioritize.evidence", "Candidate evidence matrix", parameters)
        return [artifact] if artifact else []

    def input_paths(self, context: ProjectContext, parameters: dict[str, Any]) -> list[Path]:
        return [*super().input_paths(context, parameters), *(artifact.path for artifact in self.parent_artifacts(context, parameters))]

    def run(self, context: ProjectContext, job: JobContext, run_fingerprint: str, parameters: dict[str, Any]) -> list[str]:
        selected, source = _selection_frame(context, parameters)
        candidates, background_audit, aliases, mapping_audit = _candidate_background(context, parameters, source, selected)
        libraries, resource_audit = declared_gene_sets(context.manifest, "gene_sets_gmt")
        gene_sets, set_audit = map_gene_sets(libraries, aliases)
        mapped_by_library = set_audit.groupby(set_audit.term.str.split(" / ").str[0]).mapped_unique_genes.sum()
        if mapped_by_library.eq(0).any():
            raise ValueError("Annotation libraries have no mapped background genes: " + ", ".join(mapped_by_library[mapped_by_library.eq(0)].index) + ". Check species and identifier mapping, not biological enrichment.")
        results = []
        groups = candidates.groupby("comparison_id", sort=False) if "comparison_id" in candidates else [("selection", candidates)]
        for comparison, group in groups:
            background = set(group.canonical_gene_id) - {""}
            chosen = group.selected_input
            part = library_overrepresentation(set(group.loc[chosen, "canonical_gene_id"]), background, gene_sets,
                                      minimum_overlap=int(parameters.get("minimum_overlap", 2)))
            part["comparison_id"] = comparison
            results.append(part)
        result = pd.concat(results, ignore_index=True) if results else pd.DataFrame()
        if not result.empty:
            from katana_next.science.statistics import bh_adjust
            result["fdr"] = bh_adjust(result.pvalue)
            result["row_id"] = result.comparison_id + "::" + result.term
            result["fdr_family"] = "all eligible comparison-by-library-term tests in this run"
            result = result.sort_values(["fdr", "comparison_id", "term"])
        if result.empty:
            result = pd.DataFrame(columns=["row_id", "term", "selected_genes", "background_genes", "term_genes_in_background", "overlap_genes", "expected_overlap", "fold_enrichment", "pvalue", "fdr", "genes"])
        result["minus_log10_fdr"] = -np.log10(pd.to_numeric(result.get("fdr"), errors="coerce").clip(lower=np.finfo(float).tiny)) if "fdr" in result else np.nan
        from katana_next.science.enrichment_views import annotate_categories, enrichment_bubble_plot
        result = annotate_categories(result, context.manifest)
        category = parameters.get("ontology_category", "All")
        if not result.empty and category not in {"All", "Unclassified"} and not result.ontology_category.eq(category).any():
            raise ValueError(f"No terms classified as {category}; attach a versioned GO ontology or explicitly classified library. This is not a negative enrichment result.")
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version, common_metadata={"question": self.definition.question, "method": self.definition.method, "analysis_unit": self.definition.analysis_unit, "interpretation": "Terms with low FDR and coherent overlapping genes prioritize biological themes; they do not establish mechanism.", "limitations": self.definition.limitations, "background_definition": "genes with available evidence in this contrast intersected with all mapped genes in the same annotation library; pre-annotation counts remain in the table"})
        table = writer.table("Gene-set enrichment", "gene_set_enrichment.parquet", result, entity_type="gene_set", parents=(source.artifact_id,))
        writer.table("Enrichment gene inputs by comparison", "gene_inputs.parquet", candidates, entity_type="gene_input_audit", parents=(source.artifact_id,))
        writer.table("Enrichment background eligibility", "background_eligibility.parquet", background_audit, entity_type="background_audit", parents=(source.artifact_id,))
        writer.table("Gene identifier mapping", "gene_identifier_mapping.parquet", mapping_audit, entity_type="mapping_audit", parents=(source.artifact_id,))
        if "gene_symbol_source" in candidates:
            writer.table("Project GTF gene-name resolution", "gene_name_resolution.parquet", candidates[["gene_id", "gene_symbol", "gene_symbol_source"]].drop_duplicates(), entity_type="mapping_audit", parents=(source.artifact_id,))
        writer.table("Gene-set mapping audit", "gene_set_mapping.parquet", set_audit, entity_type="mapping_audit", parents=(source.artifact_id,))
        writer.table("Annotation resource provenance", "annotation_resources.parquet", resource_audit, entity_type="resource_audit")
        writer.plot("Gene-set enrichment", "gene_set_enrichment.vl.json", enrichment_bubble_plot(
            result, category=parameters.get("ontology_category", "All"),
            x_axis=parameters.get("bubble_x", "gene_ratio"), top_terms=int(parameters.get("top_terms", 12)),
            labels={c.contrast_id: c.label for c in context.manifest.contrasts}),
            entity_type="gene_set_summary", parents=(table.artifact_id,))
        return writer.finish()


class RbpSupportWorkflow(WorkflowBase):
    version = "2.2.0"
    definition = WorkflowDefinition(
        workflow_id="prioritize.rbp_support",
        workspace_id="prioritize",
        label="Auxiliary RBP identity / declared target sets",
        presentation_role="auxiliary",
        question="Which candidates are annotated RBPs/splicing factors, or belong to explicitly declared target sets?",
        method="Comparison-specific hypergeometric enrichment against measured/testable genes, with BH across all comparison-by-set tests. RBP membership and declared target sets retain their resource semantics.",
        output_summary="RBP enrichment, matched candidate genes and cross-contrast support matrix.",
        analysis_unit="declared RBP membership or target set and candidate gene",
        expected_artifacts=5,
        supports_contrast_selection=True,
        parameters=(ParameterDefinition("selection_id", "Saved gene selection", "selection", "", "Optional candidate-gene subset."), ParameterDefinition("minimum_overlap", "Minimum overlapping genes", "number", 1, "Minimum RBP-set overlap.", minimum=1, step=1)),
        required_asset_kinds=("rbp_gene_sets",),
        limitations="RBP overlap prioritizes regulators or targets and is never reported as causal control.",
    )

    def parent_artifacts(self, context: ProjectContext, parameters: dict[str, Any]) -> list[AnalysisArtifact]:
        artifact = selected_source(context, parameters) or _latest_table(context, "prioritize.evidence", "Candidate evidence matrix", parameters)
        return [artifact] if artifact else []

    def input_paths(self, context: ProjectContext, parameters: dict[str, Any]) -> list[Path]:
        return [*super().input_paths(context, parameters), *(artifact.path for artifact in self.parent_artifacts(context, parameters))]

    def run(self, context: ProjectContext, job: JobContext, run_fingerprint: str, parameters: dict[str, Any]) -> list[str]:
        selected, source = _selection_frame(context, parameters)
        candidates, background_audit, aliases, mapping_audit = _candidate_background(context, parameters, source, selected)
        genes = set(candidates.loc[candidates.selected_input, "canonical_gene_id"]) - {""}
        libraries, resource_audit = declared_gene_sets(context.manifest, "rbp_gene_sets")
        gene_sets, set_audit = map_gene_sets(libraries, aliases)
        if not any(gene_sets.values()):
            raise ValueError("RBP resources have no mapped measured-background genes; check species and identifiers, not biological enrichment.")
        groups = candidates.groupby("comparison_id", sort=False) if "comparison_id" in candidates else [("selection", candidates)]
        results, membership_rows = [], []
        for comparison, group in groups:
            selected_genes = set(group.loc[group.selected_input, "canonical_gene_id"]) - {""}
            background = set(group.canonical_gene_id) - {""}
            part = overrepresentation(selected_genes, background, gene_sets, minimum_overlap=int(parameters.get("minimum_overlap", 1)))
            part["comparison_id"] = comparison
            results.append(part)
            for rbp, members in gene_sets.items():
                for gene in sorted(selected_genes & members):
                    membership_rows.append({"row_id": f"{comparison}::{rbp}::{gene}", "comparison_id": comparison, "rbp_or_set": rbp, "gene_id": gene})
        enrichment = pd.concat(results, ignore_index=True) if results else pd.DataFrame()
        if not enrichment.empty:
            from katana_next.science.statistics import bh_adjust
            enrichment["fdr"] = bh_adjust(enrichment.pvalue)
            enrichment["fdr_family"] = "all eligible comparison-by-RBP-set tests in this run"
            enrichment["row_id"] = enrichment.comparison_id + "::" + enrichment.term
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
        membership = pd.DataFrame(membership_rows, columns=["row_id", "comparison_id", "rbp_or_set", "gene_id"])
        target_genes = set().union(*gene_sets.values()) if gene_sets else set()
        selected_candidates = candidates[candidates.selected_input]
        matched = selected_candidates[
            selected_candidates["canonical_gene_id"].isin(target_genes)
        ].copy()
        matched["supporting_rbp_sets"] = matched.apply(
            lambda row: "; ".join(
                sorted(
                    name
                    for name, members in gene_sets.items()
                    if str(row.get("canonical_gene_id", "")) in members
                )
            ),
            axis=1,
        )
        matched["is_rbp_or_target"] = matched["supporting_rbp_sets"].ne("")
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version, common_metadata={"question": self.definition.question, "method": self.definition.method, "analysis_unit": self.definition.analysis_unit, "interpretation": "Use these results for mechanism prioritization only. Motif/eCLIP and expression association are separate supporting layers.", "limitations": self.definition.limitations})
        writer.common_metadata["resource_semantics"] = "RBP membership is not RBP-target regulation. Consult Annotation resource provenance for each declared set."
        writer.table("RBP background eligibility", "background_eligibility.parquet", background_audit, entity_type="background_audit", parents=(source.artifact_id,))
        writer.table("RBP gene inputs by comparison", "gene_inputs.parquet", candidates, entity_type="gene_input_audit", parents=(source.artifact_id,))
        writer.table("Annotation resource provenance", "annotation_resources.parquet", resource_audit, entity_type="resource_audit")
        enrichment_artifact = writer.table("RBP declared-set enrichment", "rbp_enrichment.parquet", enrichment, entity_type="rbp_gene_set", parents=(source.artifact_id,))
        writer.table("RBP identifier mapping", "rbp_identifier_mapping.parquet", mapping_audit, entity_type="mapping_audit", parents=(source.artifact_id,))
        if "gene_symbol_source" in candidates:
            writer.table("Project GTF gene-name resolution", "gene_name_resolution.parquet", candidates[["gene_id", "gene_symbol", "gene_symbol_source"]].drop_duplicates(), entity_type="mapping_audit", parents=(source.artifact_id,))
        writer.table("RBP database mapping audit", "rbp_database_mapping.parquet", set_audit, entity_type="mapping_audit", parents=(source.artifact_id,))
        writer.table("RBP matched genes", "rbp_matched_genes.parquet", membership, entity_type="gene", parents=(enrichment_artifact.artifact_id, source.artifact_id))
        writer.table("RBP candidate support matrix", "rbp_candidate_support.parquet", matched, entity_type="gene", parents=(source.artifact_id,))
        from katana_next.workflows.support_plots import rbp_support_plots
        rbp_support_plots(writer, candidates, matched, genes, gene_sets,
                          {c.contrast_id: c.label for c in context.manifest.contrasts}, (source.artifact_id,), enrichment=enrichment)
        return writer.finish()


class IsoformSupportWorkflow(WorkflowBase):
    version = "2.1.1"
    definition = WorkflowDefinition(
        workflow_id="prioritize.isoform_support",
        workspace_id="prioritize",
        label="Isoform expression and switch support",
        question="Do candidate genes have transcript-level expression and dominant-isoform changes?",
        method="Salmon quant.sf aggregation through the declared GTF plus separate DTU evidence.",
        output_summary="Transcript TPM, gene TPM, contrast delta TPM and descriptive dominant-switch candidates.",
        analysis_unit="transcript within a biological sample",
        expected_artifacts=5,
        supports_contrast_selection=True,
        limitations="A dominant transcript change is descriptive unless supported by an independent DTU/switch test.",
    )

    def input_paths(self, context: ProjectContext, parameters: dict[str, Any]) -> list[Path]:
        paths = [sample.assets["quant_sf"] for sample in context.manifest.samples if "quant_sf" in sample.assets]
        if context.manifest.annotation.gtf:
            paths.append(context.manifest.annotation.gtf)
        return [*paths, *(a.path for a in self.parent_artifacts(context, parameters))]

    def parent_artifacts(self, context: ProjectContext, parameters: dict[str, Any]) -> list[AnalysisArtifact]:
        source = _latest_table(context, "differential.dtu", "DTU and exon-usage evidence", parameters)
        return [source] if source else []

    def run(self, context: ProjectContext, job: JobContext, run_fingerprint: str, parameters: dict[str, Any]) -> list[str]:
        cache = context.manifest.output_root / ".katana" / "transcript_gene_map.parquet"
        job.report(0.08, "Resolving transcript-to-gene annotation")
        included_samples = {s for c in context.manifest.contrasts if not parameters.get("contrast_ids") or c.contrast_id in parameters["contrast_ids"] for s in (*c.numerator.sample_ids, *c.denominator.sample_ids)}
        transcript, gene = load_quant_support(context, cache, sample_ids=included_samples or None)
        contrasts = isoform_contrast_summary(context, transcript)
        switches = dominant_isoform_switches(context, transcript)
        if parameters.get("contrast_ids"):
            contrasts = contrasts[contrasts.comparison_id.isin(parameters["contrast_ids"])].copy()
            switches = switches[switches.comparison_id.isin(parameters["contrast_ids"])].copy()
            included_samples = {sample for c in context.manifest.contrasts if c.contrast_id in parameters["contrast_ids"] for sample in (*c.numerator.sample_ids, *c.denominator.sample_ids)}
            transcript = transcript[transcript.sample_id.isin(included_samples)].copy()
            gene = gene[gene.sample_id.isin(included_samples)].copy()
        sources = self.parent_artifacts(context, parameters)
        if sources:
            dtu = pd.read_parquet(sources[0].path)
            gene_dtu = usage_gene_support(dtu)
            contrasts = contrasts.merge(gene_dtu, on=["comparison_id", "gene_id"], how="left", validate="many_to_one")
            switches = switches.merge(gene_dtu, on=["comparison_id", "gene_id"], how="left", validate="many_to_one")
        else:
            contrasts["best_dtu_fdr"] = np.nan
            switches["best_dtu_fdr"] = np.nan
            contrasts["best_dexseq_fdr"] = np.nan
            switches["best_dexseq_fdr"] = np.nan
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version, common_metadata={"question": self.definition.question, "method": self.definition.method, "analysis_unit": self.definition.analysis_unit, "interpretation": "TPM supports abundance and visualization. DTU significance must remain a separate statistical result.", "limitations": self.definition.limitations})
        writer.common_metadata["dtu_source"] = sources[0].artifact_id if sources else "not available; DTU is not run implicitly"
        writer.common_metadata["dtu_gene_summary"] = "best_dtu_fdr and best_dexseq_fdr are separate minima of imported transcript-usage and exon-usage q-values. Neither is a new gene-level or switch-specific test."
        transcript_artifact = writer.table("Sample transcript abundance", "transcript_tpm.parquet", transcript, entity_type="transcript", parents=tuple(a.artifact_id for a in sources))
        writer.table("Sample gene abundance", "gene_tpm.parquet", gene, entity_type="gene", parents=(transcript_artifact.artifact_id,))
        contrast_artifact = writer.table("Isoform delta TPM by contrast", "isoform_delta_tpm.parquet", contrasts, entity_type="transcript", parents=(transcript_artifact.artifact_id,))
        switch_artifact = writer.table("Descriptive dominant isoform switches", "dominant_isoform_switches.parquet", switches, entity_type="gene", parents=(transcript_artifact.artifact_id,))
        plot_values = contrasts.assign(abs_delta_tpm=lambda frame: frame["delta_tpm"].abs()).sort_values("abs_delta_tpm", ascending=False).groupby("comparison_id").head(25)
        writer.plot("Isoform abundance shifts", "isoform_delta_tpm.vl.json", _plot(plot_values.to_dict("records"), {"type": "bar", "cornerRadiusEnd": 2}, {"x": {"field": "delta_tpm", "type": "quantitative", "axis": {"title": "Mean TPM difference (numerator minus denominator)"}}, "y": {"field": "transcript_id", "type": "nominal", "sort": {"field": "delta_tpm", "order": "descending"}, "axis": {"title": None, "labelLimit": 160}}, "color": {"condition": {"test": "datum.delta_tpm >= 0", "value": TEAL}, "value": "#c2410c"}, "facet": {"field": "comparison_id", "columns": 2}, "tooltip": [{"field": "comparison_id"}, {"field": "gene_symbol"}, {"field": "transcript_id"}, {"field": "delta_tpm", "format": ".2f"}]}, title="Largest transcript abundance shifts", description="Descriptive Salmon TPM support; use DTU for significance.", height=500), entity_type="transcript_summary", parents=(contrast_artifact.artifact_id, switch_artifact.artifact_id))
        from katana_next.workflows.support_plots import isoform_support_plots
        isoform_support_plots(writer, contrasts, switches, {c.contrast_id: c.label for c in context.manifest.contrasts}, (contrast_artifact.artifact_id, switch_artifact.artifact_id))
        return writer.finish()
