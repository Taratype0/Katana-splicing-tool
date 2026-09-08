from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from katana_next.domain.models import AnalysisArtifact
from katana_next.infrastructure.jobs import JobContext
from katana_next.science.single_cell import (
    differential_pseudobulk,
    donor_celltype_pseudobulk,
    extract_junction_choices_from_bams,
    load_single_cell_counts,
    load_single_cell_metadata,
    single_cell_audit,
)
from katana_next.services.project_service import ProjectContext
from katana_next.workflows.base import (
    ArtifactWriter,
    ParameterDefinition,
    WorkflowBase,
    WorkflowDefinition,
)
from katana_next.workflows.bulk import RED, _latest_table, _plot


class SingleCellAuditWorkflow(WorkflowBase):
    version = "2.0.0"
    definition = WorkflowDefinition(
        workflow_id="single_cell.audit",
        workspace_id="single-cell",
        label="Single-cell assay and design audit",
        question="Can this assay support exploratory splicing and donor-aware differential claims?",
        method="Assay declaration, cell identity, donor structure, cell-count and junction-input gates.",
        output_summary="Claim gate table, donor-by-cell-type design table and cell-count plot.",
        analysis_unit="cell for measurement; donor/sample for inference",
        expected_artifacts=3,
        parameters=(ParameterDefinition("min_cells_per_donor_cell_type", "Minimum cells per donor and cell type", "number", 20, "Coverage gate for donor-cell-type pseudobulk.", minimum=1, step=5),),
        required_asset_kinds=("single_cell_metadata",),
        optional_asset_kinds=("single_cell_junction_counts",),
        limitations="Without donor/sample IDs Katana permits exploration only and blocks population-level significance.",
    )

    def run(self, context: ProjectContext, job: JobContext, run_fingerprint: str, parameters: dict[str, Any]) -> list[str]:
        metadata = load_single_cell_metadata(context)
        try:
            counts = load_single_cell_counts(context)
        except ValueError:
            counts = None
        gates, design = single_cell_audit(metadata, counts, min_cells_per_donor_cell_type=int(parameters.get("min_cells_per_donor_cell_type", 20)))
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version, common_metadata={"question": self.definition.question, "method": self.definition.method, "analysis_unit": self.definition.analysis_unit, "interpretation": "A pass establishes an analysis design, not biological significance.", "limitations": self.definition.limitations})
        gates_artifact = writer.table("Single-cell claim gates", "single_cell_claim_gates.parquet", gates, entity_type="quality_gate")
        design_artifact = writer.table("Donor by cell-type design", "single_cell_design.parquet", design, entity_type="donor_cell_type", parents=(gates_artifact.artifact_id,))
        writer.plot(
            "Single-cell donor design",
            "single_cell_design.vl.json",
            _plot(design.to_dict("records"), {"type": "bar", "cornerRadiusEnd": 3}, {"x": {"field": "donor_id", "type": "nominal", "axis": {"title": "Donor", "labelAngle": -45}}, "y": {"field": "cells", "type": "quantitative", "axis": {"title": "Cells"}}, "color": {"field": "condition", "type": "nominal", "scale": {"scheme": "tableau10"}}, "column": {"field": "cell_type", "type": "nominal", "header": {"title": "Cell type"}}, "tooltip": [{"field": "donor_id"}, {"field": "condition"}, {"field": "cell_type"}, {"field": "cells"}, {"field": "passes_cell_count"}]}, title="Independent-unit coverage by cell type", description="Cells are displayed for coverage; donor remains the inferential unit.", height=300),
            entity_type="donor_cell_type_summary",
            parents=(design_artifact.artifact_id,),
        )
        return writer.finish()


class SingleCellQuantifyWorkflow(WorkflowBase):
    version = "2.0.1"
    definition = WorkflowDefinition(
        workflow_id="single_cell.quantify",
        workspace_id="single-cell",
        label="Junction usage and donor pseudobulk",
        question="Which junction choices are measurably used in each donor and cell type?",
        method="UMI-aware input counts or BAM junction extraction, GTF splice-site annotation, donor-by-cell-type count aggregation and posterior PSI.",
        output_summary="Cell-event counts, donor pseudobulk, coverage table and uncertainty heatmap.",
        analysis_unit="junction choice within donor and cell type",
        expected_artifacts=4,
        parameters=(
            ParameterDefinition("source", "Count source", "select", "precomputed", "Use validated long-form counts or extract CB/UB-tagged BAM junctions.", ("precomputed", "bam")),
            ParameterDefinition("min_cell_event_reads", "Cell-event read gate", "number", 1, "Minimum reads to call a cell-event detected.", minimum=0, step=1),
            ParameterDefinition("cell_barcode_tag", "Cell barcode tag", "text", "CB", "SAM tag carrying corrected cell barcode."),
            ParameterDefinition("umi_tag", "UMI tag", "text", "UB", "SAM tag carrying corrected UMI."),
            ParameterDefinition("min_mapping_quality", "Minimum mapping quality", "number", 20, "BAM extraction mapping-quality gate.", minimum=0, maximum=255, step=1),
            ParameterDefinition("max_reads_per_bam", "BAM read cap", "number", 0, "Zero means no cap; positive values are for explicit smoke tests only.", minimum=0, step=10000),
        ),
        required_asset_kinds=("single_cell_metadata",),
        optional_asset_kinds=("single_cell_junction_counts",),
        limitations="Droplet data are sparse. Donor pseudobulk is the confirmatory path; cell-level PSI is descriptive.",
    )

    def input_paths(self, context: ProjectContext, parameters: dict[str, Any]) -> list[Path]:
        paths = super().input_paths(context, parameters)
        if str(parameters.get("source", "precomputed")) == "bam":
            paths.extend(
                path
                for sample in context.manifest.samples
                for kind, path in sample.assets.items()
                if kind in {"single_cell_bam", "bam"} and path.exists()
            )
        return list(dict.fromkeys(path.resolve() for path in paths))

    def run(self, context: ProjectContext, job: JobContext, run_fingerprint: str, parameters: dict[str, Any]) -> list[str]:
        metadata = load_single_cell_metadata(context)
        source = str(parameters.get("source", "precomputed"))
        if source == "bam":
            job.report(0.08, "Extracting UMI-deduplicated junction choices from BAM")
            raw_path = context.results.artifact_directory(self.definition.workflow_id, run_fingerprint) / "extracted_single_cell_junction_counts.tsv"
            extract_junction_choices_from_bams(context, raw_path, cell_barcode_tag=str(parameters.get("cell_barcode_tag", "CB")), umi_tag=str(parameters.get("umi_tag", "UB")), min_mapping_quality=int(parameters.get("min_mapping_quality", 20)), max_reads_per_bam=int(parameters.get("max_reads_per_bam", 0)))
            counts = pd.read_csv(raw_path, sep="\t")
            counts["total"] = counts["inclusion"] + counts["exclusion"]
            counts["psi_posterior"] = (counts["inclusion"] + 0.5) / (counts["total"] + 1.0)
            counts["row_id"] = counts["cell_id"].astype(str) + "::" + counts["event_id"].astype(str)
        else:
            counts = load_single_cell_counts(context)
        unknown = sorted(set(counts["cell_id"]) - set(metadata["cell_id"]))
        if unknown:
            raise ValueError(f"Count input contains cells absent from metadata: {unknown[:5]}")
        job.report(0.38, "Aggregating donor-by-cell-type junction counts")
        pseudobulk, coverage = donor_celltype_pseudobulk(metadata, counts, min_cell_event_reads=float(parameters.get("min_cell_event_reads", 1)))
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version, common_metadata={"question": self.definition.question, "method": self.definition.method, "analysis_unit": self.definition.analysis_unit, "interpretation": "Posterior PSI summarizes count support. It does not create independent replicates from cells.", "limitations": self.definition.limitations})
        raw_artifact = writer.table("Single-cell junction choices", "single_cell_junction_counts.parquet", counts, entity_type="cell_splice_event")
        pseudo_artifact = writer.table("Donor-cell-type pseudobulk", "donor_celltype_pseudobulk.parquet", pseudobulk, entity_type="donor_cell_type_event", parents=(raw_artifact.artifact_id,))
        coverage_artifact = writer.table("Pseudobulk coverage", "pseudobulk_coverage.parquet", coverage, entity_type="donor_cell_type", parents=(pseudo_artifact.artifact_id,))
        writer.plot(
            "Single-cell splicing coverage",
            "single_cell_coverage.vl.json",
            _plot(coverage.to_dict("records"), {"type": "circle", "opacity": 0.85, "stroke": "white"}, {"x": {"field": "donor_id", "type": "nominal", "axis": {"title": "Donor", "labelAngle": -45}}, "y": {"field": "cell_type", "type": "nominal", "axis": {"title": None}}, "size": {"field": "measured_events", "type": "quantitative", "legend": {"title": "Measured events"}}, "color": {"field": "median_measurement_se", "type": "quantitative", "scale": {"scheme": "orangered", "reverse": True}, "legend": {"title": "Median PSI SE"}}, "tooltip": [{"field": "donor_id"}, {"field": "condition"}, {"field": "cell_type"}, {"field": "cells"}, {"field": "measured_events"}, {"field": "median_event_reads", "format": ".1f"}, {"field": "median_measurement_se", "format": ".3f"}]}, title="Donor-aware junction coverage", description="Point size shows measured event count; color shows count-informed uncertainty.", height=max(330, 38 * coverage["cell_type"].nunique())),
            entity_type="donor_cell_type_summary",
            parents=(coverage_artifact.artifact_id,),
        )
        return writer.finish()


class SingleCellDifferentialWorkflow(WorkflowBase):
    version = "2.0.0"
    definition = WorkflowDefinition(
        workflow_id="single_cell.differential",
        workspace_id="single-cell",
        label="Donor-aware cell-type differential splicing",
        question="Which splice choices differ between conditions within each cell type?",
        method="Beta-binomial regression on donor-by-cell-type inclusion and total junction counts with BH FDR per cell type.",
        output_summary="Cell-type event effects, model diagnostics, forest-like effect plot and significance summary.",
        analysis_unit="biological donor within cell type",
        expected_artifacts=3,
        parameters=(
            ParameterDefinition("numerator", "Numerator condition", "condition", "", "Effect direction is dPSI = numerator minus denominator."),
            ParameterDefinition("denominator", "Denominator condition", "condition", "", "Reference condition."),
            ParameterDefinition("min_total_reads", "Minimum donor-event reads", "number", 10, "Pseudobulk event read gate.", minimum=1, step=1),
            ParameterDefinition("min_donors_per_condition", "Minimum donors per condition", "number", 2, "Hard inferential gate; three or more is strongly preferred.", minimum=2, step=1),
            ParameterDefinition("fdr", "FDR", "number", 0.05, "BH FDR within cell type.", minimum=0, maximum=1, step=0.01),
            ParameterDefinition("min_abs_dpsi", "Minimum |dPSI|", "number", 0.10, "Effect-size gate.", minimum=0, maximum=1, step=0.01),
        ),
        limitations="With only two donors per arm, estimates are allowed but labeled low-power; cells never inflate degrees of freedom.",
    )

    def parent_artifacts(self, context: ProjectContext, parameters: dict[str, Any]) -> list[AnalysisArtifact]:
        del parameters
        artifact = _latest_table(context, "single_cell.quantify", "Donor-cell-type pseudobulk")
        return [artifact] if artifact else []

    def input_paths(self, context: ProjectContext, parameters: dict[str, Any]) -> list[Path]:
        return [artifact.path for artifact in self.parent_artifacts(context, parameters)]

    def run(self, context: ProjectContext, job: JobContext, run_fingerprint: str, parameters: dict[str, Any]) -> list[str]:
        parent = _latest_table(context, "single_cell.quantify", "Donor-cell-type pseudobulk")
        if parent is None:
            raise ValueError("Run Junction usage and donor pseudobulk first. Differential analysis never starts it implicitly.")
        pseudobulk = pd.read_parquet(parent.path)
        conditions = sorted(pseudobulk["condition"].dropna().astype(str).unique())
        numerator = str(parameters.get("numerator", "")).strip() or (conditions[1] if len(conditions) > 1 else "")
        denominator = str(parameters.get("denominator", "")).strip() or (conditions[0] if conditions else "")
        if not numerator or not denominator or numerator == denominator:
            raise ValueError("Choose two different single-cell conditions.")
        job.report(0.18, f"Fitting donor-aware beta-binomial models: {numerator} vs {denominator}")
        results = differential_pseudobulk(pseudobulk, numerator=numerator, denominator=denominator, min_total_reads=float(parameters.get("min_total_reads", 10)), min_donors_per_condition=int(parameters.get("min_donors_per_condition", 2)), fdr_cutoff=float(parameters.get("fdr", 0.05)), min_abs_dpsi=float(parameters.get("min_abs_dpsi", 0.10)))
        if results.empty:
            raise ValueError("No event passed the donor/read feasibility gates for this cell-type contrast.")
        summary = results.groupby(["cell_type", "direction"], dropna=False).agg(tested_events=("event_id", "size"), significant_events=("significant", "sum")).reset_index()
        summary["row_id"] = summary["cell_type"] + "::" + summary["direction"]
        plot_values = results.sort_values(["significant", "fdr"], ascending=[False, True]).groupby("cell_type").head(40)
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, {**parameters, "numerator": numerator, "denominator": denominator}, self.version, common_metadata={"question": self.definition.question, "method": self.definition.method, "analysis_unit": self.definition.analysis_unit, "direction_definition": f"dPSI = {numerator} minus {denominator}", "interpretation": "Significant events pass both model FDR and effect-size gates. Rho reports extra-binomial variation.", "limitations": self.definition.limitations})
        table = writer.table("Cell-type differential splicing", "single_cell_differential_as.parquet", results, entity_type="cell_type_splice_event", parents=(parent.artifact_id,))
        writer.table("Cell-type DAS summary", "single_cell_das_summary.parquet", summary, entity_type="cell_type_summary", parents=(table.artifact_id,))
        writer.plot(
            "Cell-type splicing effects",
            "single_cell_das_effects.vl.json",
            _plot(plot_values.to_dict("records"), {"type": "point", "filled": True, "size": 60}, {"x": {"field": "dpsi", "type": "quantitative", "axis": {"title": f"dPSI ({numerator} minus {denominator})"}}, "y": {"field": "gene_symbol", "type": "nominal", "sort": {"field": "dpsi", "order": "descending"}, "axis": {"title": None, "labelLimit": 150}}, "color": {"field": "significant", "type": "nominal", "scale": {"domain": [True, False], "range": [RED, "#9aa5ad"]}}, "facet": {"field": "cell_type", "type": "nominal", "columns": 2}, "tooltip": [{"field": "cell_type"}, {"field": "gene_symbol"}, {"field": "event_id"}, {"field": "dpsi", "format": ".3f"}, {"field": "fdr", "format": ".2e"}, {"field": "numerator_donors"}, {"field": "denominator_donors"}, {"field": "overdispersion_rho", "format": ".3f"}]}, title="Donor-aware cell-type splicing effects", description="Top 40 events per cell type; full inferential table is available separately.", height=500),
            entity_type="cell_type_summary",
            parents=(table.artifact_id,),
        )
        job.report(0.95, "Single-cell differential splicing complete")
        return writer.finish()
