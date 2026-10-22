from __future__ import annotations
from pathlib import Path
from dataclasses import asdict
import json
import pandas as pd
from katana_next.workflows.base import WorkflowBase, WorkflowDefinition, ParameterDefinition, ArtifactWriter
from katana_next.science.psi_network import fit_psi_network, psi_filter_audit, PsiFilterConfig
from katana_next.science.frozen_psi import freeze_psi_modules, save_frozen_model, _samples
from katana_cosplice.config import NetworkConfig


class PsiNetworkWorkflow(WorkflowBase):
    version = "1.3.0"
    definition = WorkflowDefinition("programs.psi_network", "programs", "cosplice: imported PSI networks",
        "Which measured splice events vary together across biological samples?",
        "Signed WGCNA event network and 50 biological-block confidence refits. PSI-only quality checks; no invented junction counts.",
        "Event modules, bootstrap confidence, edges, sample scores, PSI matrix and input-loss audit", "splice event",
        expected_artifacts=13, required_asset_kinds=("psi_matrix", "event_annotation"),
        parameters=(ParameterDefinition("rscript", "Rscript executable", "text", "", "Use an installed R runtime with WGCNA."),
            ParameterDefinition('filter_profile','PSI eligibility','select','coverage_variable','Coverage and nonconstant variation are distinct from large-effect PSI gates. Both profiles retain per-event reasons.',options=('coverage_variable','legacy')),
            ParameterDefinition('max_features','Maximum variable events','number',2000,'Training variance ranking; no outcome selection.',minimum=40,maximum=3000),
            ParameterDefinition("method", "Correlation", "select", "bicor", "Passed to official R WGCNA.", options=("bicor", "pearson")),
            ParameterDefinition("signed", "Signed network", "boolean", True, "Keep the existing signed default."),
            ParameterDefinition("soft_power", "Soft power", "number", 6, "Existing fixed power default.", minimum=1, maximum=30),
            *(ParameterDefinition(k, label, "number", value, "Discovery samples only; all failed checks are recorded.", minimum=0, maximum=1)
              for k, label, value in (("min_valid_fraction", "Minimum observed fraction", .8), ("min_variance", "Minimum PSI variance", .0025),
                  ("min_range", "Minimum PSI range", .5), ("min_alternative_fraction", "Minimum alternative fraction", .1))),
            ParameterDefinition("feature_space_id", "PSI feature definition", "text", "", "Build + annotation + event definition + quantification version. Blank uses a project-scoped identity."),
            ParameterDefinition("training_traits", "Optional numeric training traits TSV", "text", "", "Rows indexed by sample_id; coefficients and module selection are learned on discovery only."),
            ParameterDefinition("trait_prediction_method", "Trait prediction", "select", "validated_hierarchy", "Training-only model choice; unavailable modules can use a frozen reference or a validated alternative, with each route recorded.", options=("single_module", "multimodule_ridge",'validated_hierarchy')),
            ParameterDefinition('training_reference_scores','Optional training reference scores TSV','text','','Numeric sample-by-feature scores, such as expression PCs fitted only on training samples. This is a prediction reference, not a splice module.'),
            ParameterDefinition('reference_space_id','Reference transform identity','text','','Required with reference scores; projection must use the same frozen transform.')),
        limitations="Imported PSI cannot establish junction-read support, library coverage QC or repeat-supported DAS. Network associations do not prove regulation.")

    def input_paths(self, context, parameters):
        paths = super().input_paths(context, parameters)
        if parameters.get("training_traits"):
            paths.append(Path(parameters["training_traits"]).resolve())
        if parameters.get('training_reference_scores'):paths.append(Path(parameters['training_reference_scores']).resolve())
        return paths

    def run(self, context, job, run_fingerprint, parameters):
        assets = {a.kind: a.path for a in context.manifest.assets}
        matrix = pd.read_csv(assets["psi_matrix"], sep="\t", index_col=0)
        samples = pd.DataFrame([{"sample_id": s.sample_id, "biological_sample_id": s.biological_sample_id, "condition": s.condition} for s in context.manifest.samples])
        if set(matrix.columns) != set(samples.sample_id):
            raise ValueError("PSI matrix columns must match all declared project samples")
        _samples(matrix, samples)
        events = pd.read_csv(assets["event_annotation"], sep="\t")
        filters = PsiFilterConfig(**{k: parameters.get(k, v) for k, v in asdict(PsiFilterConfig()).items()})
        network = NetworkConfig(wgcna_rscript=parameters.get("rscript") or None,
            method=parameters.get("method", "bicor"), signed=parameters.get("signed", True), soft_power=parameters.get("soft_power", 6))
        qc_all = psi_filter_audit(matrix, filters)
        profile=parameters.get('filter_profile','coverage_variable');cap=int(parameters.get('max_features',2000))
        if profile=='coverage_variable':
            from katana_next.science.psi_robust import coverage_variable_audit
            qc_all=coverage_variable_audit(matrix,filters,max_features=cap)
        qc_all["failed_gene_filter"] = False
        if context.manifest.metadata.get("gene_filter"):
            from katana_next.science.gene_filter import apply_gene_filter
            events = apply_gene_filter(context, events)
            key = next((c for c in ("event_key", "event_id") if c in events), None)
            if key is None:
                raise ValueError("Filtered PSI networks require event_key/event_id gene annotations.")
            matrix = matrix.loc[matrix.index.isin(events[key])]
            qc_all["failed_gene_filter"] = ~qc_all.index.isin(matrix.index)
            events = events.loc[events[key].isin(matrix.index)].copy()
        if profile=='coverage_variable':
            effective=coverage_variable_audit(matrix,filters,max_features=cap)
            qc_all.loc[effective.index,effective.columns]=effective
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version,
            common_metadata={"input_scope": "imported_PSI_only", "bootstrap_iterations": 50,
                             "read_support": "not_assessed", "repeat_supported_DAS": "not_assessed"})
        qc_all["passed"] &= ~qc_all.failed_gene_filter
        qc_all["failure_count"] += qc_all.failed_gene_filter.astype(int)
        qc_all.loc[qc_all.failed_gene_filter, "filter_reasons"] = qc_all.loc[qc_all.failed_gene_filter, "filter_reasons"].replace("passed", "").map(lambda x: (x+";" if x else "")+"gene_filter")
        # Persist all rows even if too few events remain to fit a network.
        qc_export = qc_all.rename_axis("event_key").reset_index()
        qc_export.insert(0, "row_id", [f"psi_qc:{i}" for i in range(len(qc_export))])
        writer.table("stable_qc__metrics", "stable_qc__metrics.parquet", qc_export, entity_type="splice_event")
        job.report(.02, "Validating measured PSI and biological units")
        traits = pd.read_csv(parameters["training_traits"], sep="\t", index_col=0) if parameters.get("training_traits") else None
        reference=pd.read_csv(parameters['training_reference_scores'],sep='\t',index_col=0) if parameters.get('training_reference_scores') else None
        if reference is not None and (traits is None or parameters.get('trait_prediction_method','validated_hierarchy')!='validated_hierarchy'):
            raise ValueError('Reference fallback requires training traits and validated_hierarchy')
        if traits is not None and (not traits.index.is_unique or set(traits.index) != set(matrix.columns)):
            raise ValueError("Training trait rows must match discovery samples exactly")
        measured, qc, means, baseline, confidence = fit_psi_network(matrix, samples, network=network, filters=filters,
            progress=lambda msg: job.report(.4, msg),filter_profile=profile,max_features=cap)
        feature_space = parameters.get("feature_space_id") or context.manifest.metadata.get("psi_feature_space_id") or f"{context.manifest.annotation.genome_build}|project:{context.manifest.project_id}"
        frozen = freeze_psi_modules(measured, means, baseline, samples, feature_space_id=feature_space, network=network, filters=filters, traits=traits,
            trait_prediction_method=parameters.get("trait_prediction_method", "validated_hierarchy") if traits is not None else 'single_module',
            reference_scores=reference,reference_space_id=parameters.get('reference_space_id'),filter_profile=profile)
        model_path = writer.directory/"frozen_psi_modules.json"
        save_frozen_model(frozen, model_path)
        writer.existing_file("frozen_psi__model", model_path, entity_type="module", metadata={"feature_space_id": feature_space})
        resolved_path = writer.directory/"resolved_psi_configuration.json"
        resolved_path.write_text(json.dumps({"network": asdict(network), "filters": asdict(filters), "feature_space_id": feature_space,
            'filter_profile':profile,'max_features':cap,'network_status':frozen['network_status']}, indent=2), encoding="utf-8")
        writer.existing_file("frozen_psi__configuration", resolved_path, entity_type="module")
        gate_summary = pd.DataFrame([{"reason": c.removeprefix("failed_"), "failed_events": int(qc_all[c].sum()), "counting": "overlapping_failures"} for c in qc_all if c.startswith("failed_")])
        gate_summary.insert(0, "row_id", [f"psi_gate:{i}" for i in range(len(gate_summary))])
        writer.table("stable_qc__reason_summary", "stable_qc__reason_summary.parquet", gate_summary, entity_type="splice_event")
        outputs = {"consensus_modules__assignments": confidence.assignments,
            "consensus_modules__modules": confidence.module_summary, "consensus_modules__bootstrap": confidence.bootstrap_runs,
            "network_observed__edges": baseline.edges, "network_observed__scores": baseline.program_scores.reset_index(),
            "consensus_modules__scores": confidence.program_scores.reset_index(),
            "stable_qc__matrix": measured.rename_axis("event_key").reset_index(),
            "matrix__events": events, "matrix__imputation_means": means.rename("discovery_mean_psi").rename_axis("event_key").reset_index()}
        for name, frame in outputs.items():
            if "row_id" not in frame:
                frame.insert(0, "row_id", [f"{name}:{i}" for i in range(len(frame))])
            writer.table(name, name+".parquet", frame, entity_type="splice_event" if "event_key" in frame else "module")
        job.report(.95, "Registering immutable module views")
        return writer.finish()
