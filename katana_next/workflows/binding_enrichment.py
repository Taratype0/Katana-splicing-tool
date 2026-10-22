from pathlib import Path

import pandas as pd

from katana_next.science.binding_enrichment import paired_binding_enrichment
from katana_next.science.motif_significance import extract_windows, match_background, measurable_event_background
from katana_next.science.versioned_evidence import binding_support, verified_resource, read_binding_peaks
from katana_next.workflows.base import ArtifactWriter, ParameterDefinition, WorkflowBase, WorkflowDefinition
from katana_next.workflows.figure_style import readable_plot as _plot
from katana_next.workflows.sources import selected_source, selection_rows


class BindingEnrichmentWorkflow(WorkflowBase):
    version = "1.0.2"
    definition = WorkflowDefinition(
        "prioritize.binding_enrichment", "prioritize", "RBP CLIP-site enrichment by comparison",
        "Do selected splice boundaries overlap known experimental RBP binding sites more often than matched measurable events?",
        "One event per gene; exact event type, GC within 0.05 and sequence length within 20%. Exact paired test with BH across declared CLIP resources in one comparison.",
        "Selected events, tested-event universe, matched controls, complete CLIP overlaps, P/FDR and resource provenance.",
        "independent matched event-gene pair and reference CLIP experiment",
        parameters=(ParameterDefinition("selection_id", "Saved event selection", "selection", "", "Choose significant, read-supported events from one exact comparison; no pooling."),
                    ParameterDefinition("window_bp", "Boundary window per side (bp)", "number", 100, "Same boundary window for foreground and controls.", minimum=10, maximum=1000, step=10)),
        required_asset_kinds=("rbp_binding_peaks",),
        presentation_role="advanced",
        limitations="Reference CLIP is experimental binding in its own biosample, not binding in the user's cells. Enrichment is association, not causal regulation. Motif catalogs are not CLIP peaks. Unmatched or unevaluable pairs are retained and excluded from the statistical denominator.",
    )

    def parent_artifacts(self, context, parameters):
        source = selected_source(context, parameters)
        if source is None:
            raise ValueError("Save exact events from an unfiltered, measurable-event result table")
        return [source]

    def input_paths(self, context, parameters):
        paths = [*super().input_paths(context, parameters), *(a.path for a in self.parent_artifacts(context, parameters))]
        fasta = context.manifest.annotation.genome_fasta
        if not fasta or not fasta.is_file() or not Path(str(fasta) + ".fai").is_file():
            raise ValueError("Declared genome FASTA and index are required for matched sequence background")
        paths.extend([fasta, Path(str(fasta) + ".fai")])
        for asset in context.manifest.assets:
            if asset.kind == "rbp_binding_peaks":
                verified_resource(asset, species=context.manifest.metadata.get("species", ""), genome_build=context.manifest.annotation.genome_build)
        return paths

    def run(self, context, job, run_fingerprint, parameters):
        source = self.parent_artifacts(context, parameters)[0]
        universe = pd.read_parquet(source.path)
        foreground = selection_rows(context, parameters, source, universe)
        if "comparison_id" not in foreground or foreground.comparison_id.nunique() != 1:
            raise ValueError("Select exact events from one comparison; no pooled CLIP enrichment")
        comparison = foreground.comparison_id.iloc[0]
        universe = universe.loc[universe.comparison_id.eq(comparison)].copy()
        required = {"event_id", "row_id", "gene_id", "event_type", "chr", "strand"}
        if not required <= set(universe):
            raise ValueError(f"Exact measurable event catalog needs {sorted(required)}")
        universe, eligibility = measurable_event_background(universe)
        universe = universe.drop_duplicates("event_id")
        if not set(foreground.event_id) <= set(universe.event_id):
            raise ValueError("Foreground contains events not eligible in the measurable-event universe")
        assets = [a for a in context.manifest.assets if a.kind == "rbp_binding_peaks"]
        for asset in assets:
            chroms = set(read_binding_peaks(asset.path)[0].astype(str))
            for chrom in universe.chr.astype(str).unique():
                alternate = chrom.removeprefix("chr") if chrom.startswith("chr") else "chr" + chrom
                if chrom not in chroms and alternate in chroms:
                    raise ValueError("CLIP/event chromosome naming differs; harmonize reference names explicitly before analysis")
        job.report(.1, "Matching independent, measurable control genes without examining CLIP hits")
        windows, features = extract_windows(universe, context.manifest.annotation.genome_fasta, int(parameters["window_bp"]))
        pairs = match_background(features, set(foreground.event_id))
        if pairs.empty:
            raise ValueError("No foreground events")
        controls = set(pairs.loc[pairs.status.eq("matched"), "background_event"])
        query = universe.loc[universe.event_id.isin(set(foreground.event_id) | controls)]
        evidence, resources = binding_support(query, assets, species=context.manifest.metadata.get("species", ""),
            genome_build=context.manifest.annotation.genome_build, window_bp=int(parameters["window_bp"]))
        result = paired_binding_enrichment(evidence, pairs)
        pairs["comparison_id"] = comparison
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version,
            common_metadata={"method": self.definition.method, "limitations": self.definition.limitations, "background_source": source.artifact_id})
        for label, name, table in (("CLIP enrichment input events", "input_events", foreground),
                                  ("Measurable-event background eligibility", "eligibility", eligibility),
                                  ("Sequence-matched CLIP controls", "matching", pairs),
                                  ("Measurable-event sequence features", "background_features", features),
                                  ("Queried boundary sequence windows", "windows", windows.drop(columns="sequence")),
                                  ("All reference CLIP overlaps", "overlaps", evidence),
                                  ("Matched-event CLIP enrichment", "enrichment", result),
                                  ("Experimental CLIP provenance", "resources", resources)):
            table = table.copy()
            if "row_id" not in table:
                table["row_id"] = table.index.astype(str)
            writer.table(label, name + ".parquet", table, entity_type="binding_enrichment", parents=(source.artifact_id,))
        writer.plot("CLIP enrichment in matched measurable events", "enrichment.vl.json", _plot(result.to_dict("records"), "circle", {
            "x": {"field": "hit_fraction_difference", "type": "quantitative", "title": "Selected minus matched-control binding fraction"},
            "y": {"field": "RBP", "type": "nominal"}, "color": {"field": "enrichment_fdr", "type": "quantitative", "title": "FDR"},
            "size": {"field": "evaluable_gene_pairs", "type": "quantitative"}, "row": {"field": "biosample", "type": "nominal"},
            "tooltip": [{"field": c} for c in ("comparison_id", "RBP", "accession", "evaluable_gene_pairs", "foreground_hit_genes", "background_hit_genes", "enrichment_pvalue", "enrichment_fdr")],
        }, title="Known RBP binding-site enrichment", description=self.definition.limitations), entity_type="binding_enrichment_plot", parents=(source.artifact_id,))
        return writer.finish()
