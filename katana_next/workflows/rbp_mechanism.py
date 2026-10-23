"""One event-level RBP entry point; evidence layers remain separate."""
from pathlib import Path

import pandas as pd

from katana_next.science.versioned_evidence import verified_resource
from katana_next.workflows.base import ArtifactWriter, ParameterDefinition, WorkflowBase, WorkflowDefinition
from katana_next.workflows.binding_enrichment import BindingEnrichmentWorkflow
from katana_next.workflows.motif import EventMotifWorkflow
from katana_next.workflows.sources import selected_source, selection_rows


class _StageProgress:
    def __init__(self, job, start, span):
        self.job, self.start, self.span = job, start, span

    def check_cancelled(self):
        return self.job.check_cancelled()

    def report(self, progress, message):
        return self.job.report(self.start + self.span * progress, message)


class RbpMechanismWorkflow(WorkflowBase):
    version = "1.6.1"
    definition = WorkflowDefinition(
        "prioritize.rbp_mechanism", "prioritize", "Katana MotifTrace: event motifs and known binding sites",
        "Do the selected DAS events share RBP sequence motifs or overlap experimentally measured RBP binding sites?",
        "Boundary mode uses independent gene-matched controls and separate motif/CLIP families. Completed SE maps reuse the frozen public all-model/window/contrast family, complete short RNA and genuine joins. Optional user models default to position preview; explicit custom enrichment uses a separate BH family. No quantification or CLIP validation is started by completed maps.",
        "RNA sequence/position browser, native model logos and source classes, complete unchanged P/FDR, background audit and real RNA/split-coordinate exports; legacy motif/CLIP outputs remain available.",
        "selected splice event, independent matched gene pair and reference RBP experiment",
        required_asset_kinds=(),
        optional_asset_kinds=("rbp_motif_database", "fimo_runtime", "rbp_binding_peaks", "rbp_sequence_results", "motiftrace_matched_windows", "motiftrace_parallel_results"),
        supports_contrast_selection=True,
        parameters=(
            ParameterDefinition("selection_id", "Saved DAS events from one comparison", "selection", "", "Select read-supported DAS from the complete measurable event table. The selected rows are never used as their own background."),
            ParameterDefinition("sequence_scope", "RNA sequence analysis", "select", "boundary_association", "completed_parallel_models reuses original MAPP sample fits and MotifTrace calibration with separate P values and real RNA/usage evidence. No fits or quantification. Other modes remain unchanged.", ("boundary_association", "se_sequence_impact", "completed_sequence_maps", "event_consequences", "completed_parallel_models")),
            ParameterDefinition("parallel_bundle_asset_id", "Completed MAPP / MotifTrace bundle", "text", "", "Parallel mode only: registered SHA-bound result bundle. Choose one original contrast, not a new selected-event background."),
            ParameterDefinition("model_ids", "Consequence mode: explicit motif IDs", "text", "", "Optional native model IDs or hashes, comma-separated; selected before scanning, never by significance."),
            ParameterDefinition("cluster_permutations", "Consequence mode: paired permutations", "number", 999, "New cluster FWER and window BH families; original P/FDR unchanged.", minimum=99, maximum=9999),
            ParameterDefinition("max_tests", "Consequence mode: maximum complete-family tests", "number", 1000000, "Execution budget checked before scoring; never truncate models or windows and never split the statistical family.", minimum=1),
            ParameterDefinition("max_matrix_cells", "Consequence mode: maximum paired-matrix cells", "number", 8000000, "Execution memory budget checked before scoring; 8 bytes per float64 cell plus tables and permutation arrays.", minimum=1),
            ParameterDefinition("bundle_asset_id", "Completed result asset ID", "text", "", "Completed maps only: empty selects the unique registered completed RBP bundle."),
            ParameterDefinition("window_bundle_asset_id", "Matched count / coverage result asset ID (optional)", "text", "", "Completed maps only: registered original count and coverage windows, with separate P/FDR and external-support labels; no rescoring."),
            ParameterDefinition("view_scope", "Display comparisons", "select", "selected_comparison", "Completed maps only; original multi-comparison FDR is not recalculated.", ("selected_comparison", "all_frozen_comparisons")),
            ParameterDefinition("custom_motif_path", "User motif library path (optional)", "text", "", "Completed maps: local consensus TSV or Minimal MEME plus metadata. Species and source checked; original data and results are preserved."),
            ParameterDefinition("custom_motif_mode", "User motif analysis", "select", "preview_only", "Positions only by default. Explicit enrichment uses cached RNA and a separate custom-library FDR family, not an improved-sensitivity claim.", ("preview_only", "fixed_enrichment")),
            ParameterDefinition("window_bp", "Boundary window per side (bp)", "number", 100, "The same strand-oriented window is used for selected and matched control events.", minimum=10, maximum=1000, step=10),
            ParameterDefinition("motif_evidence", "Motif evidence layers", "select", "all_evidence_separate", "Direct experiment, homology inference and model prediction remain separate; no pooled stronger-binding claim.", ("direct_experiment", "all_evidence_separate")),
            ParameterDefinition("site_pvalue", "FIMO site P threshold", "number", .0001, "Screening threshold, not enrichment FDR.", minimum=.00000001, maximum=.01, step=.0001),
            ParameterDefinition("binding_mode", "Experimental binding evidence", "select", "if_configured", "Require a verified species/build-matched CLIP resource, or explicitly record missing input. Never replace CLIP with predicted motif sites.", ("if_configured", "required")),
        ),
        limitations="Motif preference is not measured binding. Reference CLIP shows binding in its documented biosample, not in the user's cells. Neither layer establishes regulatory causality. Auxiliary RBP identity remains available under advanced analyses.",
    )

    def parent_artifacts(self, context, parameters):
        if parameters.get("sequence_scope") == "completed_parallel_models":
            return []
        source = selected_source(context, parameters)
        if source is None:
            raise ValueError("Save exact DAS events from a complete measurable-event table first")
        return [source]

    def input_paths(self, context, parameters):
        if parameters.get("sequence_scope") == "completed_parallel_models":
            from katana_next.workflows.motiftrace_parallel import MotifTraceParallelWorkflow
            workflow = MotifTraceParallelWorkflow()
            return workflow.input_paths(context, parameters) + [Path(__file__)]
        if parameters.get("sequence_scope") == "event_consequences":
            from katana_next.workflows.motiftrace_consequences import MotifTraceConsequencesWorkflow
            workflow = MotifTraceConsequencesWorkflow()
            args = workflow.normalize_parameters({p.key: parameters[p.key] for p in workflow.definition.parameters if p.key in parameters})
            return workflow.input_paths(context, args) + [Path(__file__)]
        if parameters.get("sequence_scope") == "completed_sequence_maps":
            from katana_next.workflows.rbp_sequence_results import CompletedRbpSequenceWorkflow
            workflow = CompletedRbpSequenceWorkflow()
            return workflow.input_paths(context, workflow.normalize_parameters({key: parameters[key] for key in
                ("selection_id", "bundle_asset_id", "window_bundle_asset_id", "view_scope", "custom_motif_path", "custom_motif_mode")})) + [Path(__file__)]
        if parameters.get("sequence_scope") == "se_sequence_impact":
            from katana_next.workflows.se_rbp_impact import SeRbpImpactWorkflow
            return SeRbpImpactWorkflow().input_paths(context, parameters)
        EventMotifWorkflow().validate_assets(context)
        paths = EventMotifWorkflow().input_paths(context, parameters)
        for asset in context.manifest.assets:
            if asset.kind in {"rbp_motif_database", "rbp_binding_peaks"}:
                verified_resource(asset, species=context.manifest.metadata.get("species", ""),
                                  genome_build=context.manifest.annotation.genome_build if asset.kind == "rbp_binding_peaks" else None)
        peaks = [a for a in context.manifest.assets if a.kind == "rbp_binding_peaks"]
        if parameters.get("binding_mode") == "required" and not peaks:
            raise ValueError("Known-site enrichment requires experimental CLIP peaks")
        paths.extend(a.path for a in peaks)
        if peaks:
            paths.extend(BindingEnrichmentWorkflow().input_paths(context, parameters))
        paths.extend(Path(__file__).with_name(n) for n in ("motif.py", "binding_enrichment.py"))
        science = Path(__file__).resolve().parents[1] / "science"
        paths.extend(science / n for n in ("motif_significance.py", "binding_enrichment.py", "versioned_evidence.py"))
        paths.append(Path(__file__).resolve().parents[2] / "katana_cosplice" / "motif_scanner.py")
        return list(dict.fromkeys(paths))

    def run(self, context, job, run_fingerprint, parameters):
        if parameters.get("sequence_scope") == "completed_parallel_models":
            from katana_next.workflows.motiftrace_parallel import MotifTraceParallelWorkflow
            workflow = MotifTraceParallelWorkflow()
            child_ids = workflow.run(context, job, run_fingerprint + "-parallel", parameters)
            writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version,
                common_metadata={"new_inference": False, "quantification_rerun": False, "method_P_values_separate": True})
            for child_id in child_ids:
                artifact = context.catalog.get_artifact(child_id)
                writer.existing_file(artifact.name, artifact.path, kind=artifact.kind, entity_type=artifact.entity_type,
                    metadata={**artifact.metadata, "evidence_source_workflow": artifact.analysis_id}, parents=(child_id,))
            return writer.finish()
        if parameters.get("sequence_scope") == "event_consequences":
            from katana_next.workflows.motiftrace_consequences import MotifTraceConsequencesWorkflow
            workflow = MotifTraceConsequencesWorkflow()
            args = workflow.normalize_parameters({p.key: parameters[p.key] for p in workflow.definition.parameters if p.key in parameters})
            child_ids = workflow.run(context, job, run_fingerprint + "-consequences", args)
            source = self.parent_artifacts(context, parameters)[0]
            writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version,
                common_metadata={"new_inference": True, "historical_P_FDR_unchanged": True, "quantification_rerun": False})
            for child_id in child_ids:
                artifact = context.catalog.get_artifact(child_id)
                writer.existing_file(artifact.name, artifact.path, kind=artifact.kind, entity_type=artifact.entity_type,
                    metadata=artifact.metadata, parents=(source.artifact_id, child_id))
            return writer.finish()
        if parameters.get("sequence_scope") == "completed_sequence_maps":
            from katana_next.workflows.rbp_sequence_results import CompletedRbpSequenceWorkflow
            workflow = CompletedRbpSequenceWorkflow()
            args = workflow.normalize_parameters({key: parameters[key] for key in
                ("selection_id", "bundle_asset_id", "window_bundle_asset_id", "view_scope", "custom_motif_path", "custom_motif_mode")})
            child_ids = workflow.run(context, job, run_fingerprint + "-completed-maps", args)
            source = self.parent_artifacts(context, parameters)[0]
            writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version,
                common_metadata={"reuse_completed_results": True, "new_inference": bool(parameters.get("custom_motif_path")) and parameters["custom_motif_mode"] == "fixed_enrichment",
                                 "public_inference_unchanged": True, "tool_name": "Katana MotifTrace"})
            for child_id in child_ids:
                artifact = context.catalog.get_artifact(child_id)
                writer.existing_file(artifact.name, artifact.path, kind=artifact.kind, entity_type=artifact.entity_type,
                    metadata={**artifact.metadata, "evidence_source_workflow": artifact.analysis_id, "workflow_version": self.version},
                    parents=(source.artifact_id, child_id), warnings=artifact.warnings)
            return writer.finish()
        if parameters.get("sequence_scope") == "se_sequence_impact":
            from katana_next.workflows.se_rbp_impact import SeRbpImpactWorkflow
            workflow = SeRbpImpactWorkflow()
            return workflow.run(context, job, run_fingerprint+"-se", workflow.normalize_parameters(parameters))
        source = self.parent_artifacts(context, parameters)[0]
        events = selection_rows(context, parameters, source, pd.read_parquet(source.path))
        if events.empty or "comparison_id" not in events or events.comparison_id.nunique() != 1:
            raise ValueError("RBP analysis needs nonempty exact events from one comparison; no pooling")
        comparison = str(events.comparison_id.iloc[0])
        child_ids = []
        motif = EventMotifWorkflow()
        motif_parameters = motif.normalize_parameters({k: v for k, v in parameters.items() if k in {p.key for p in motif.definition.parameters}})
        child_ids.extend(motif.run(context, _StageProgress(job, 0, .65), run_fingerprint + "-motif", motif_parameters))
        stages = [{"row_id": "motif", "comparison_id": comparison, "evidence_layer": "sequence_motif",
                   "status": "complete", "meaning": "Selected events versus sequence-matched measurable controls; evidence classes are separate."}]
        if any(a.kind == "rbp_binding_peaks" for a in context.manifest.assets):
            binding = BindingEnrichmentWorkflow()
            args = binding.normalize_parameters({"selection_id": parameters["selection_id"], "window_bp": parameters["window_bp"]})
            child_ids.extend(binding.run(context, _StageProgress(job, .65, .30), run_fingerprint + "-binding", args))
            stages.append({"row_id": "binding", "comparison_id": comparison, "evidence_layer": "experimental_CLIP",
                           "status": "complete", "meaning": "Measured reference binding near splice boundaries; each experiment retains its biosample and assay."})
        else:
            if parameters["binding_mode"] == "required":
                raise ValueError("Experimental CLIP peaks are missing")
            stages.append({"row_id": "binding", "comparison_id": comparison, "evidence_layer": "experimental_CLIP",
                           "status": "missing_input", "meaning": "Not computed; absence of a configured CLIP resource is not a negative binding result."})
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version,
                                common_metadata={"comparison_id": comparison, "method": self.definition.method,
                                                 "limitations": self.definition.limitations, "no_combined_evidence_pvalue": True})
        writer.table("RBP analysis stages and how to read them", "stages.parquet", pd.DataFrame(stages),
                     entity_type="rbp_evidence_stage", parents=(source.artifact_id,))
        by_id = {a.artifact_id: a for a in context.catalog.list_artifacts(context.manifest.project_id)}
        # Main-entry proxies keep every child file in the cache verification contract.
        for child_id in child_ids:
            artifact = by_id[child_id]
            writer.existing_file(artifact.name, artifact.path, kind=artifact.kind, entity_type=artifact.entity_type,
                                 metadata={"exports": artifact.metadata.get("exports", {}),
                                           "column_definitions": artifact.metadata.get("column_definitions", {}),
                                           "evidence_source_workflow": artifact.analysis_id},
                                 parents=(source.artifact_id, child_id), warnings=artifact.warnings)
        job.report(.98, "RBP evidence layers complete; all nonsignificant results and matching losses retained")
        return writer.finish()
