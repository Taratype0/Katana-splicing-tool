from pathlib import Path
import pandas as pd

from katana_next.workflows.base import WorkflowBase, WorkflowDefinition, ParameterDefinition, ArtifactWriter
from katana_next.science.frozen_psi import load_frozen_model, project_psi_modules


class PsiProjectionWorkflow(WorkflowBase):
    version = "1.2.0"
    definition = WorkflowDefinition("programs.psi_projection", "programs", "cosplice: validate modules in new samples",
        "Do saved PSI modules have measured support in independent samples?",
        "Freeze discovery event membership, imputation, scaling, PC loadings and optional trait coefficients; project without refitting.",
        "New-sample scores, per-event/sample coverage, measured coherence and optional held-out trait prediction.", "module",
        expected_artifacts=7,
        parameters=(
            ParameterDefinition("model", "Frozen PSI model JSON", "text", "", "Use frozen_psi_modules.json from the PSI network workflow."),
            ParameterDefinition("psi", "New-sample PSI TSV", "text", "", "Rows are exact event IDs; columns are sample IDs. Missing values are allowed."),
            ParameterDefinition("samples", "New-sample metadata TSV", "text", "", "sample_id and biological_sample_id; no discovery samples or technical duplicates."),
            ParameterDefinition("feature_space_id", "Matching PSI feature definition", "text", "", "Must equal the discovery definition, including genome build and quantification."),
            ParameterDefinition("traits", "Optional numeric new-sample traits TSV", "text", "", "Only evaluates saved coefficients; never selects a new module."),
            ParameterDefinition("min_observed_fraction", "Required measured module fraction", "number", .8, "Per sample; inadequate scores are missing, never zero.", minimum=.01, maximum=1),
            ParameterDefinition("projection_method", "Missing-event scoring", "select", "auto", "Auto retains mean scoring for legacy models and guarded scoring for robust models. Guarded scoring preserves mean-score coverage.", options=('auto',"mean_impute", "observed_loading_ls",'guarded_loading')),
            ParameterDefinition('reference_scores','Optional new-sample reference scores TSV','text','','Use exactly the reference transform frozen on training samples; no refit on validation data.'),
            ParameterDefinition('reference_space_id','Reference transform identity','text','','Must match the frozen predictor when reference scores are provided.'),
            ParameterDefinition("min_loading_fraction", "Required measured loading information", "number", .8, "Observed-loading LS additionally requires this fraction of squared PC1 loadings.", minimum=.01, maximum=1)),
        limitations="Projection and coherence do not establish regulation or formal WGCNA/NetRep preservation significance. Fewer than 15 independent samples cannot estimate coherence.")

    def input_paths(self, context, parameters):
        return [Path(parameters[k]).expanduser().resolve() for k in ("model", "psi", "samples", "traits",'reference_scores') if parameters.get(k)]

    def run(self, context, job, run_fingerprint, parameters):
        for key in ("model", "psi", "samples", "feature_space_id"):
            if not parameters.get(key):
                raise ValueError(f"Supply {key}")
        model = load_frozen_model(parameters["model"])
        matrix = pd.read_csv(parameters["psi"], sep="\t", index_col=0)
        samples = pd.read_csv(parameters["samples"], sep="\t")
        traits = pd.read_csv(parameters["traits"], sep="\t", index_col=0) if parameters.get("traits") else None
        reference=pd.read_csv(parameters['reference_scores'],sep='\t',index_col=0) if parameters.get('reference_scores') else None
        projection=parameters.get('projection_method','auto')
        if projection=='auto':projection='guarded_loading' if model.get('schema')=='katana-frozen-psi-1.2' else 'mean_impute'
        job.report(.1, "Checking independence and projecting frozen modules")
        outputs = project_psi_modules(model, matrix, samples, feature_space_id=parameters["feature_space_id"],
            min_observed_fraction=parameters["min_observed_fraction"], traits=traits,
            projection_method=projection, min_loading_fraction=parameters.get("min_loading_fraction", .8),
            reference_scores=reference,reference_space_id=parameters.get('reference_space_id'))
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version,
            common_metadata={"frozen_model_sha256": model["sha256"], "feature_space_id": model["feature_space_id"],
                "module_refit": False, "regulatory_accuracy": "not_tested",
                "projection_method": projection,
                "trait_prediction_method": model.get("trait_prediction_method", "single_module")})
        labels = {"scores": "New-sample scores", "event_audit": "Event coverage", "sample_coverage": "Sample coverage",
                  "module_validation": "Module validation", "predictions": "Predictions", "trait_validation": "Trait prediction checks"}
        for name, frame in outputs.items():
            if name == "scores":
                frame = frame.reset_index()
            frame.insert(0, "row_id", [f"psi_projection_{name}:{i}" for i in range(len(frame))])
            writer.table(labels[name], name+".parquet", frame, entity_type="module", metadata={"output_key": name})
        model_path = writer.directory/"frozen_psi_modules.json"
        from katana_next.science.frozen_psi import save_frozen_model
        save_frozen_model(model, model_path)
        writer.existing_file("Frozen model", model_path, entity_type="module")
        return writer.finish()
