"""Discoverable postprocessing entry for MotifTrace features 1, 3 and 4."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from katana_cosplice.motif_scanner import load_motif_database
from katana_next.science.motiftrace_execution import analyze_consequences, EXECUTION_ADAPTER_VERSION
from katana_next.science.motiftrace_custom import load_custom_library
from katana_next.science.rbp_sequence_results import digest
from katana_next.science.versioned_evidence import verified_resource
from katana_next.workflows.base import ArtifactWriter, ParameterDefinition, WorkflowBase, WorkflowDefinition
from katana_next.workflows.sources import selected_source, selection_rows


class MotifTraceConsequencesWorkflow(WorkflowBase):
    version = "1.0.1"
    definition = WorkflowDefinition(
        workflow_id="prioritize.motiftrace_consequences", workspace_id="prioritize",
        label="Katana MotifTrace: paired RNA and region statistics",
        question="Which motif instances change between the two local RNA forms, and which matched-gene regions show sequence association?",
        method="Existing DAS/PSI only; SE/A3SS/A5SS/MXE/RI reference RNA reconstruction; locus-paired instances; fixed motif-blind same-type gene matching; whole-pair swaps and maximum cluster-mass FWER. New window BH is separate from all historical results.",
        output_summary="Actual local RNA FASTA, gained/lost/retained loci, split genomic coordinates, replicate PSI-weighted opportunities, matching losses, complete new window P/FDR and region FWER P, portable event report.",
        analysis_unit="local splice event and independent matched gene pair, not a new biological replicate",
        required_asset_kinds=(), optional_asset_kinds=("rbp_motif_database",),
        parameters=(
            ParameterDefinition("selection_id", "Existing supported DAS events", "selection", "", "One comparison from its full measured catalog; all five rMATS event types supported. Successful quantification is never rerun."),
            ParameterDefinition("custom_motif_path", "User motif library (optional)", "text", "", "Species-matched consensus TSV or Minimal MEME plus metadata. Otherwise use registered verified public motif databases."),
            ParameterDefinition("model_ids", "Explicit model IDs (optional)", "text", "", "Comma-separated native IDs or exact model hashes. Empty means all supplied models, up to 256, subject to explicit execution budgets; never silently truncate or select by significance."),
            ParameterDefinition("cluster_permutations", "Whole-pair permutations", "number", 999, "Fixed before scanning; shared swaps across all models, methods and regions. Cluster P is FWER, not FDR.", minimum=99, maximum=9999),
            ParameterDefinition("max_tests", "Maximum complete-family tests", "number", 1000000, "Execution budget, not a scientific filter. Exact scalar/window count is checked before motif scoring. All models and windows remain in one family; exceeding this value fails without truncation.", minimum=1),
            ParameterDefinition("max_matrix_cells", "Maximum paired-matrix cells", "number", 8000000, "Execution memory budget, checked before scoring. Float64 cells use 8 bytes each; tables and permutation arrays require additional memory. Never remove matched pairs to fit this budget.", minimum=1),
        ),
        presentation_role="advanced",
        limitations="Local reference forms are not phased complete transcripts. Motif availability/PSI weighting is not measured occupancy or causality. PWM candidate positions have no calibrated site P. Matching exchangeability is an assumption; fewer than 20 independent pairs yields NA. No HPC, CLIP, read quantification or source-result mutation.",
    )

    def parent_artifacts(self, context, parameters):
        source = selected_source(context, parameters)
        if source is None:
            raise ValueError("Save supported DAS from its full measured source catalog first")
        return [source]

    def _library(self, context, parameters):
        species = context.manifest.metadata.get("species", "")
        if not species:
            raise ValueError("Explicit project species required")
        if parameters.get("custom_motif_path"):
            library = load_custom_library(parameters["custom_motif_path"], species=species)
        else:
            models, inputs = {}, []
            import hashlib
            for asset in context.manifest.assets:
                if asset.kind != "rbp_motif_database":
                    continue
                verified_resource(asset, species=species)
                inputs.append({"path": str(asset.path), "sha256": digest(asset.path)})
                for motif in load_motif_database(asset.path, organism=species):
                    if motif.pwm is not None:
                        matrix = np.asarray(motif.pwm, float)
                        matrix = matrix/matrix.sum(axis=1, keepdims=True)
                        sha = hashlib.sha256(matrix.astype("<f8").tobytes()).hexdigest()
                        model = {"kind": "PWM", "matrix": matrix.tolist(), "width": len(matrix)}
                    elif motif.consensus:
                        from katana_next.science.motiftrace_custom import fixed_pattern
                        pattern, width = fixed_pattern(motif.consensus)
                        sha = hashlib.sha256(pattern.encode()).hexdigest()
                        model = {"kind": "consensus", "pattern": pattern, "width": width}
                    else:
                        continue
                    token = "p_" + model["kind"] + "_" + sha[:20]
                    entry = models.setdefault(token, {**model, "model_SHA256": sha, "assignments": []})
                    entry["assignments"].append({"native_id": motif.motif_id, "RBP": motif.rbp,
                        "source": motif.source, "evidence": motif.evidence_class, "species": motif.target_species})
            for model in models.values():
                model.update(RBPs=";".join(sorted({a["RBP"] for a in model["assignments"]})),
                    native_ids=";".join(sorted({a["native_id"] for a in model["assignments"]})))
            library = {"models": models, "inputs": inputs}
        requested = {s.strip() for s in str(parameters.get("model_ids", "")).split(",") if s.strip()}
        if requested:
            found, chosen = set(), {}
            for token, model in library["models"].items():
                matches = requested & ({token} | set(model.get("native_ids", "").split(";")))
                if matches:
                    found.update(matches)
                    chosen[token] = model
            if found != requested:
                raise ValueError("Unknown requested motif IDs: " + ", ".join(sorted(requested-found)))
            library = {**library, "models": chosen}
        from katana_next.science.motiftrace_events import validate_models
        validate_models(library["models"])
        return library

    def input_paths(self, context, parameters):
        source = self.parent_artifacts(context, parameters)[0]
        fasta = context.manifest.annotation.genome_fasta
        if not fasta or not context.manifest.annotation.genome_build:
            raise ValueError("Declared genomic FASTA/build required; do not use transcript FASTA")
        library = self._library(context, parameters)
        science = Path(__file__).resolve().parents[1]/"science"
        return [source.path, Path(fasta), Path(str(fasta)+".fai"), Path(__file__),
            Path(__file__).resolve().parents[1]/"reports/motif_native_index.py",
            *(Path(row["path"]) for row in library["inputs"]),
            *(science/name for name in ("motiftrace_events.py", "motiftrace_clusters.py", "motiftrace_consequences.py",
                "motiftrace_consequences.html", "motiftrace_budget.py", "motiftrace_execution.py", "motiftrace_custom.py", "se_rbp.py", "statistics.py")),
            Path(__file__).resolve().parents[2]/"katana_cosplice/motif_scanner.py"]

    def run(self, context, job, run_fingerprint, parameters):
        source = self.parent_artifacts(context, parameters)[0]
        complete = pd.read_parquet(source.path)
        selected = selection_rows(context, parameters, source, complete)
        if selected.empty or selected.comparison_id.nunique() != 1:
            raise ValueError("Select events from exactly one comparison")
        complete = complete[complete.comparison_id.eq(selected.comparison_id.iloc[0])].copy()
        library = self._library(context, parameters)
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version,
            common_metadata={"tool_name": "Katana MotifTrace", "quantification_rerun": False,
                "execution_adapter_version": EXECUTION_ADAPTER_VERSION,
                "historical_P_FDR_unchanged": True, "new_inference": False, "cluster_correction": "maximum cluster mass FWER",
                "genome_build": context.manifest.annotation.genome_build, "limitations": self.definition.limitations})
        job.report(.1, "Freezing five-type RNA consequences and a separate new paired-swap statistical family")
        summary = analyze_consequences(complete, set(selected.event_id), library["models"],
            context.manifest.annotation.genome_fasta, writer.directory/"event_consequences",
            permutations=parameters["cluster_permutations"], check_cancelled=job.check_cancelled,
            max_tests=parameters["max_tests"], max_matrix_cells=parameters["max_matrix_cells"],
            provenance={"source_artifact_id": source.artifact_id, "source_table_SHA256": digest(source.path),
                "genome_build": context.manifest.annotation.genome_build, "species": context.manifest.metadata["species"],
                "genome_FASTA": str(context.manifest.annotation.genome_fasta), "motif_inputs": library["inputs"],
                "implementation_SHA256": {p.name: digest(p) for p in (
                    Path(__file__), Path(__file__).resolve().parents[1]/"reports/motif_native_index.py", *(Path(__file__).resolve().parents[1]/"science"/n for n in
                        ("motiftrace_events.py", "motiftrace_clusters.py", "motiftrace_consequences.py", "motiftrace_consequences.html", "motiftrace_budget.py", "motiftrace_execution.py")))}})
        writer.common_metadata.update(new_inference=summary["eligible_new_window_tests"] > 0,
            inference_status=summary["inference_status"])
        directory = writer.directory/"event_consequences"
        job.report(.9, "Registering actual RNA, locus changes, complete NA and separate region/window statistics")
        urls = {}
        for path in sorted(directory.iterdir()):
            if path.suffix in {".tsv", ".gz"}:
                writer.table("Katana MotifTrace: "+path.stem, path.name.replace(".tsv.gz", "").replace(".tsv", "")+".parquet",
                    pd.read_csv(path, sep="\t"), entity_type="motiftrace_event_evidence", parents=(source.artifact_id,))
            artifact = writer.existing_file("Katana MotifTrace: "+path.name, path, entity_type="motiftrace_consequences",
                parents=(source.artifact_id,), metadata={"mime_type": "text/html"} if path.suffix == ".html" else None)
            urls[path.name] = f"/api/projects/{context.manifest.project_id}/artifacts/{artifact.artifact_id}/file?download=true"
        # Keep the portable, sealed HTML unchanged; the app copy resolves every
        # sibling attachment through its registered, permission-checked route.
        html = (directory/"MOTIFTRACE_EVENT_CONSEQUENCES.html").read_text(encoding="utf-8")
        serialized = json.dumps(urls).replace("<", "\\u003c")
        html = html.replace("const D=", "const ATTACHMENT_URLS="+serialized+";\nconst D=", 1)
        html = html.replace("encodeURIComponent(n)", "(ATTACHMENT_URLS[n]||encodeURIComponent(n))")
        app_report = writer.directory/"MOTIFTRACE_EVENT_CONSEQUENCES_APP.html"
        app_report.write_text(html, encoding="utf-8")
        writer.existing_file("Katana MotifTrace: event consequences and regional statistics", app_report,
            entity_type="motiftrace_consequences", parents=(source.artifact_id,), metadata={"mime_type": "text/html"})
        return writer.finish()
