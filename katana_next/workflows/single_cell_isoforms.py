"""Independent single-cell transcript-usage workflow; no CoSplice dependency."""
from pathlib import Path

import pandas as pd

from katana_next.science.single_cell import load_single_cell_metadata
from katana_next.science.single_cell_isoforms import analyze_isoform_usage
from katana_next.workflows.base import ArtifactWriter, ParameterDefinition, WorkflowBase, WorkflowDefinition
from katana_next.workflows.bulk import _plot


class SingleCellIsoformWorkflow(WorkflowBase):
    version = "1.0.2"
    definition = WorkflowDefinition(
        workflow_id="single_cell.isoform_usage", workspace_id="single-cell",
        label="Transcript usage and diversity",
        question="Which genes change transcript composition within the same cell type, and does diversity persist at matched molecule depth?",
        method="Donor-cell-type transcript counts, donor-aware Dirichlet-multinomial LRT, global BH, Renyi-2/Shannon diversity and without-replacement count rarefaction.",
        output_summary="Usage fractions, donor diversity, gene DIU, unadjusted transcript effects, measurement exclusions and two separate plots.",
        analysis_unit="donor-condition-cell-type; cells are not independent experimental replicates",
        expected_artifacts=7,
        required_asset_kinds=("single_cell_metadata",),
        optional_asset_kinds=("single_cell_isoform_counts",),
        parameters=(
            ParameterDefinition("counts_path", "Transcript count TSV", "text", "", "Columns: cell_id, gene_id, transcript_id, count. Unique molecule or unambiguously assigned fragment counts, not TPM. Omitted transcript rows mean zero only in a quantified cell-gene."),
            ParameterDefinition("numerator", "Experimental condition", "condition", "", "Choose an imported metadata condition; effect = experimental minus reference."),
            ParameterDefinition("denominator", "Reference condition", "condition", "", "Choose the reference from imported metadata."),
            ParameterDefinition("min_donors", "Minimum donors per condition", "number", 3, "Exploratory minimum, not a guarantee of calibrated inference. Two-donor designs retain descriptive outputs by default.", minimum=2),
            ParameterDefinition("min_gene_count", "Minimum donor-gene count", "number", 20, "Lower-depth genes stay unmeasured, not zero usage.", minimum=1),
            ParameterDefinition("min_expressing_cells", "Minimum expressing cells", "number", 5, "Gate per donor/cell type/gene.", minimum=1),
            ParameterDefinition("rarefaction_depth", "Matched molecule count", "number", 20, "Without-replacement subsampling; units below this depth are excluded from rarefaction only.", minimum=1),
            ParameterDefinition("rarefaction_draws", "Depth sensitivity draws", "number", 30, "Monte Carlo variation is not a biological confidence interval.", minimum=1, maximum=1000),
            ParameterDefinition("adjust_batch", "Adjust declared batch", "boolean", True, "Optional batch and paired-donor effects; confounded designs remain non-estimable."),
            ParameterDefinition("seed", "Random seed", "number", 17, "Reproducible rarefaction.", minimum=0),
            ParameterDefinition("plot_transcripts", "Transcripts in effect preview", "number", 30, "Display-only largest absolute changes; full tables retain every transcript.", minimum=1, maximum=200),
            ParameterDefinition("plot_genes", "Genes in diversity preview", "number", 6, "Display-only genes with greatest total molecule support; all genes remain in exported diversity tables.", minimum=1, maximum=12),
        ),
        limitations="Requires quantified mutually exclusive transcripts. Does not reconstruct isoforms from short reads or run Iso-Seq/SQANTI3. This donor-level DM model is not the paper's cell-level separate-concentration LRT. Both pooled-molecule and mean-cell-fraction diversity are exported; only pooled counts are rarefied. Above 30 choices, the DM test pools the abundance-ranked tail into one category; all original fractions remain exported. No causal or protein-function claims. Missing genes stay unmeasured.",
    )

    def input_paths(self, context, parameters):
        raw = str(parameters.get("counts_path", "")).strip()
        if not raw:
            candidates = [a.path for a in context.manifest.assets if a.kind == "single_cell_isoform_counts"]
            if len(candidates) == 1:
                raw = str(candidates[0])
        if not raw or not Path(raw).is_file():
            raise ValueError("Choose an existing transcript count TSV; junction counts and quant.sf are not interchangeable with it.")
        return super().input_paths(context, parameters) + [Path(raw).resolve()]

    def run(self, context, job, run_fingerprint, parameters):
        paths = self.input_paths(context, parameters)
        job.report(.05, "Reading transcript counts; no bulk or upstream single-cell workflow is started")
        counts = pd.read_csv(paths[-1], sep="\t", dtype={"cell_id": str, "gene_id": str, "transcript_id": str})
        settings = {k: parameters[k] for k in ("numerator", "denominator", "min_donors", "min_gene_count", "min_expressing_cells", "rarefaction_depth", "rarefaction_draws", "seed", "adjust_batch")}
        result = analyze_isoform_usage(counts, load_single_cell_metadata(context), **settings, checkpoint=job.check_cancelled, progress=job.report)
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version,
            common_metadata={"method": self.definition.method, "limitations": self.definition.limitations,
                "reference": "10.64898/2026.09.11.748700", "fdr_family": "All tested genes and cell types in this contrast",
                "effect_definition": "Equal-donor mean fraction in numerator minus denominator; unadjusted descriptive effect",
                "inference": "Donor-level asymptotic DM LRT; rarefaction does not modify DIU test counts"})
        for name, field in (("Transcript usage", "usage"), ("Donor isoform diversity", "diversity"),
                            ("Gene differential transcript usage", "differential"), ("Transcript effect sizes", "effects"),
                            ("Transcript measurement exclusions", "diagnostics")):
            writer.table(name, f"{field}.parquet", getattr(result, field), entity_type="transcript_usage")
        job.report(.85, "Exporting donor-level plots and complete tables")
        displayed_genes = result.diversity.groupby("gene_id").gene_count.sum().nlargest(int(parameters["plot_genes"])).index
        diversity = result.diversity[result.diversity.gene_id.isin(displayed_genes)].melt(id_vars=["gene_id", "donor_id", "condition", "cell_type"],
            value_vars=["renyi2_bits", "rarefied_renyi2_mean"], var_name="measurement", value_name="bits").dropna(subset=["bits"])
        diversity["measurement"] = diversity.measurement.map({"renyi2_bits": "Observed molecules",
            "rarefied_renyi2_mean": f"Matched depth ({parameters['rarefaction_depth']} molecules)"})
        diversity_plot = _plot(diversity.to_dict("records"), "point", {
            "x": {"field": "condition", "type": "nominal", "axis": {"labelAngle": -45}},
            "y": {"field": "bits", "type": "quantitative", "title": "Renyi-2 diversity (bits)"},
            "color": {"field": "measurement", "type": "nominal"},
            "shape": {"field": "gene_id", "type": "nominal", "title": "Gene"},
            "xOffset": {"field": "gene_id", "type": "nominal"},
            "column": {"field": "cell_type", "type": "nominal", "title": "Cell type"},
            "tooltip": [{"field": v} for v in ["gene_id", "donor_id", "condition", "measurement", "bits"]],
        }, title="Transcript diversity before and after depth matching", description="Each point is one measured donor-gene. Equal-depth subsampling is a sensitivity check, not a biological replicate.", height=320)
        if diversity.empty:
            diversity_plot = empty_plot("No donor-gene passes the measurement gates; inspect exclusions.")
        else:
            diversity_plot["title"] = {"text": "Transcript diversity before and after depth matching",
                "subtitle": f"Each point: one donor-gene. Preview: up to {parameters['plot_genes']} genes with greatest molecule support; full tables retain all genes."}
        writer.plot("Observed and equal-depth diversity", "diversity.vl.json", diversity_plot, entity_type="transcript_usage")
        # A documented display limit only; exports retain every original row.
        effect = result.effects.copy()
        if not effect.empty:
            effect = effect.loc[effect.delta_fraction.abs().nlargest(int(parameters["plot_transcripts"])).index].copy()
        effect["label"] = effect.gene_id + " / " + effect.transcript_id
        effect_plot = _plot(effect.to_dict("records"), "bar", {
            "x": {"field": "delta_fraction", "type": "quantitative", "title": f"Transcript fraction: {parameters['numerator']} minus {parameters['denominator']} (unadjusted)"},
            "y": {"field": "label", "type": "nominal", "title": "Gene / transcript", "axis": {"labelLimit": 300}},
            "color": {"field": "delta_fraction", "type": "quantitative", "scale": {"scheme": "blueorange", "domainMid": 0}, "legend": None},
            "column": {"field": "cell_type", "type": "nominal"},
            "tooltip": [{"field": v} for v in ["gene_id", "transcript_id", "numerator", "denominator", "delta_fraction", "fdr", "model_status"]],
        }, title=f"Transcript usage: {parameters['numerator']} versus {parameters['denominator']}", description=f"Display: at most {parameters['plot_transcripts']} largest transcript changes. Full tables retain all genes. Gene-level FDR is not transcript-specific significance.", height=max(220, 22 * effect.label.nunique()))
        if effect.empty:
            effect_plot = empty_plot("No gene is measured in both conditions; this is not evidence of no change.")
        else:
            effect_plot["title"] = {"text": f"Transcript usage: {parameters['numerator']} versus {parameters['denominator']}",
                "subtitle": f"Up to {parameters['plot_transcripts']} largest absolute changes; all transcripts exported. Bars are descriptive, not significance calls."}
        writer.plot("Transcript composition changes", "effects.vl.json", effect_plot, entity_type="transcript_usage")
        return writer.finish()


def empty_plot(message):
    return _plot([{"message": message}], {"type": "text", "fontSize": 14},
        {"text": {"field": "message", "type": "nominal"}}, title="Transcript usage: not estimable",
        description=message, height=200)
