import pandas as pd

from katana_next.science.versioned_evidence import binding_support
from katana_next.workflows.base import WorkflowBase, WorkflowDefinition, ParameterDefinition, ArtifactWriter
from katana_next.workflows.sources import selected_source, selection_rows
from katana_next.workflows.bulk import _plot


class BindingEvidenceWorkflow(WorkflowBase):
    version = "1.0.3"
    definition = WorkflowDefinition(
        "prioritize.binding_evidence", "prioritize", "Versioned RBP binding evidence",
        "Do selected event boundaries overlap experimental CLIP peaks in a declared reference biosample?",
        "Strand-specific BED overlap in an identical genome build; accession, assay, biosample, processing version and SHA256 are mandatory.",
        "Every selected event/resource pair, reference-context binding support plot and provenance table.",
        "selected event and reference CLIP experiment", parameters=(
            ParameterDefinition("selection_id", "Saved event selection", "selection", "", "Exact selected event rows; no automatic gene-to-event substitution."),
            ParameterDefinition("window_bp", "Boundary window (bp each side)", "number", 100, "Explicit descriptive overlap window, not a significance threshold.", minimum=1, maximum=1000)),
        required_asset_kinds=("rbp_binding_peaks",),
        presentation_role="advanced",
        limitations="Reference-context binding is not binding measured in the user's cells and does not prove upstream regulation. No overlap is not proof of no binding. No enrichment p-value is invented without an appropriate measurable-event background.")

    def parent_artifacts(self, context, parameters):
        source = selected_source(context, parameters)
        if source is None:
            raise ValueError("Save an event selection before checking binding support")
        return [source]

    def input_paths(self, context, parameters):
        from katana_next.science.versioned_evidence import verified_resource
        for asset in context.manifest.assets:
            if asset.kind == "rbp_binding_peaks":
                verified_resource(asset, species=context.manifest.metadata.get("species", ""),
                                  genome_build=context.manifest.annotation.genome_build)
        return [*super().input_paths(context, parameters), *(a.path for a in self.parent_artifacts(context, parameters))]

    def run(self, context, job, run_fingerprint, parameters):
        source = selected_source(context, parameters)
        if source is None:
            raise ValueError("Save an event selection before checking binding support")
        frame = selection_rows(context, parameters, source, pd.read_parquet(source.path))
        if "comparison_id" not in frame or frame.comparison_id.nunique() != 1:
            raise ValueError("Choose exact events from one comparison; CLIP evidence is not pooled across conditions")
        evidence, resources = binding_support(frame, [a for a in context.manifest.assets if a.kind == "rbp_binding_peaks"],
            species=context.manifest.metadata.get("species", ""), genome_build=context.manifest.annotation.genome_build,
            window_bp=int(parameters["window_bp"]))
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version,
            common_metadata={"method": self.definition.method, "limitations": self.definition.limitations})
        writer.table("Reference-context RBP binding support", "binding_support.parquet", evidence, entity_type="binding_evidence", parents=(source.artifact_id,))
        writer.table("Binding resource provenance", "binding_resources.parquet", resources, entity_type="resource_audit")
        counts = evidence.groupby(["RBP", "biosample", "status"]).size().reset_index(name="event_resource_pairs")
        writer.plot("Experimental binding support by reference biosample", "binding_support.vl.json", _plot(counts.to_dict("records"), "bar",
            {"x": {"field": "event_resource_pairs", "type": "quantitative"}, "y": {"field": "RBP", "type": "nominal"},
             "color": {"field": "status", "type": "nominal"}, "row": {"field": "biosample", "type": "nominal"}, "tooltip": [{"field": c} for c in counts]},
            title="CLIP overlap in reference cells, not causal regulation", description=self.definition.limitations), entity_type="binding_summary", parents=(source.artifact_id,))
        return writer.finish()
