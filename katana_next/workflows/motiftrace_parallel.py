"""Register completed native MAPP fits beside unchanged MotifTrace results."""
from pathlib import Path

from katana_next.science.motiftrace_parallel import ASSET_KIND, ParallelModelBundle, render_parallel
from katana_next.workflows.base import ArtifactWriter, ParameterDefinition, WorkflowBase, WorkflowDefinition


class MotifTraceParallelWorkflow(WorkflowBase):
    version = "1.0.0"
    definition = WorkflowDefinition(
        workflow_id="prioritize.motiftrace_parallel", workspace_id="prioritize",
        label="Katana MotifTrace: parallel MAPP sample model",
        question="How do the existing MAPP per-sample activity model and MotifTrace usage tests differ on the same frozen motif/window family?",
        method="Read verified completed native MAPP fits and calibration beside the original MotifTrace branches. Keep each P, correction family, non-applicability and sample identity; no combined P and no fitting or quantification.",
        output_summary="Independent method P values, MAPP sample activity/SE/Z, continuous PSI and transcript use, actual RNA/motif coordinates and descriptive paired-form locus changes.",
        analysis_unit="population motif/window test; exact RNA event and transcript locus are separate descriptive units",
        required_asset_kinds=(ASSET_KIND,), supports_contrast_selection=True,
        parameters=(ParameterDefinition("parallel_bundle_asset_id", "Completed parallel result asset", "text", "", "Choose a registered SHA-bound completed bundle; existing fits and calibration are reused, never rerun."),),
        limitations="Reuse adapter, not a new arbitrary-input MAPP fitting engine. Cohort P is not an event/site P. Missing DEXSeq and missing paired RNA stay unavailable. Reference motif matches/predicted availability are not binding or causal proof; strict FWER equivalence may fail.",
    )

    def _load(self, context, parameters):
        candidates = [a for a in context.manifest.assets if a.kind == ASSET_KIND]
        if parameters.get("parallel_bundle_asset_id"):
            candidates = [a for a in candidates if a.asset_id == parameters["parallel_bundle_asset_id"]]
        if len(candidates) != 1:
            raise ValueError("Select exactly one completed parallel-model bundle")
        ids = parameters.get("contrast_ids", [])
        if len(ids) != 1:
            raise ValueError("Choose exactly one original contrast; backgrounds cannot be pooled")
        contrast = next((c for c in context.manifest.contrasts if c.contrast_id == ids[0]), None)
        if contrast is None:
            raise ValueError("Unknown project contrast")
        bundle = ParallelModelBundle(candidates[0], species=context.manifest.metadata.get("species", ""),
                                     genome_build=context.manifest.annotation.genome_build)
        return bundle, bundle.bind(contrast)

    def input_paths(self, context, parameters):
        bundle, _ = self._load(context, parameters)
        science = Path(__file__).resolve().parents[1] / "science"
        return [bundle.manifest_path, *bundle.files.values(), Path(__file__),
                science / "motiftrace_parallel.py", science / "motiftrace_parallel.html"]

    def run(self, context, job, run_fingerprint, parameters):
        bundle, case = self._load(context, parameters)
        job.report(.1, "Verifying existing MAPP and MotifTrace results; no new fitting")
        payload = bundle.payload([case], job.check_cancelled)
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version,
            common_metadata={"tool_name": "Katana MotifTrace", "comparison_id": case,
                             "method_P_values_separate": True, "new_inference": False,
                             "quantification_rerun": False, "MAPP_native_fits_reused": True,
                             "family": payload["family"], "limitations": self.definition.limitations})
        links = []
        for name, path in bundle.files.items():
            if ":" in name and not name.startswith(case + ":"):
                continue
            artifact = writer.existing_file("MotifTrace / unchanged " + name, path, entity_type="motiftrace_parallel_source")
            links.append({"label": name, "url": f"/api/projects/{context.manifest.project_id}/artifacts/{artifact.artifact_id}/file?download=true"})
        report = render_parallel(payload, writer.directory / "MOTIFTRACE_PARALLEL_MAPP.html", links)
        writer.existing_file("Katana MotifTrace / MAPP parallel methods and actual RNA", report,
                             entity_type="motiftrace_parallel_report", metadata={"mime_type": "text/html"})
        job.report(1, "Completed result reuse and exact sequence/usage projection accepted")
        return writer.finish()
