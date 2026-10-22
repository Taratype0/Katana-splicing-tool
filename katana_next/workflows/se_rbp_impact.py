"""Registered SE-specific sequence consequence layer of the RBP module."""
from pathlib import Path
import json

import pandas as pd

from katana_cosplice.motif_scanner import load_motif_database
from katana_next.science.se_rbp import POLICY, run_se_rbp
from katana_next.science.se_rbp_report import build_report, digest
from katana_next.science.public_rbp_resources import write_coverage
from katana_next.science.versioned_evidence import verified_resource
from katana_next.workflows.base import ArtifactWriter, ParameterDefinition, WorkflowBase, WorkflowDefinition
from katana_next.workflows.sources import selected_source, selection_rows


def _brand_new_report(directory):
    """Brand newly generated display files, not historical evidence or kernels."""
    directory = Path(directory)
    path = directory / "REPORT.html"
    text = path.read_text(encoding="utf-8")
    text = text.replace("<title>Katana SE RBP real evidence</title>",
                        "<title>Katana MotifTrace | legacy fixed-site SE evidence</title>")
    text = text.replace("<h1>Katana：SE 序列影响与 RBP motif</h1>",
                        "<h1>Katana MotifTrace：SE 序列影响与 RBP motif</h1>")
    path.write_text(text, encoding="utf-8")
    conclusion = directory / "CONCLUSION.md"
    lines = conclusion.read_text(encoding="utf-8").splitlines(keepends=True)
    if lines:
        lines[0] = "# Katana MotifTrace / " + lines[0].removeprefix("# ")
        conclusion.write_text("".join(lines), encoding="utf-8")


class SeRbpImpactWorkflow(WorkflowBase):
    version = "1.0.1"
    definition = WorkflowDefinition(
        workflow_id="prioritize.se_rbp_impact", workspace_id="prioritize",
        label="Katana MotifTrace: legacy SE fixed-site enrichment",
        question="Are sequences included/skipped by significant SE events enriched for known RBP motifs?",
        method="Original DAS only; transcript-oriented cassette and crossing-junction sequences; separate regulatory introns; gene-matched non-DAS controls; joint BH over RBP/evidence/region/direction.",
        output_summary="Enrichment P/FDR including negatives, positional maps, PWM logos, real RNA hit FASTA/TSV and split genomic coordinates; no refit or causal claim.",
        analysis_unit="independent matched gene pair, RNA sequence and database motif assignment",
        required_asset_kinds=("rbp_motif_database", "fimo_runtime"),
        parameters=(ParameterDefinition("selection_id", "Read-supported DAS SE selection", "selection", "", "Select only SE events from one complete measured-event catalog. Other event types fail explicitly."),
                    ParameterDefinition("motif_evidence", "Motif evidence layers", "select", "all_evidence_separate", "Separate direct experiment, homology inference and model prediction.", ("direct_experiment", "all_evidence_separate")),
                    ParameterDefinition("site_pvalue", "FIMO site P threshold", "number", .0001, "Screening threshold, not enrichment FDR.", minimum=.00000001, maximum=.01)),
        presentation_role="advanced",
        limitations="Predicted site availability, not measured binding, altered affinity, or biological causality. Pre-mRNA intronic associations are not sequences deleted by SE. Partial motif database coverage is explicit.",
    )

    def parent_artifacts(self, context, parameters):
        source = selected_source(context, parameters)
        if source is None:
            raise ValueError("Save exact read-supported SE DAS from a full measured catalog first")
        return [source]

    def input_paths(self, context, parameters):
        paths = super().input_paths(context, parameters)
        paths.extend(a.path for a in self.parent_artifacts(context, parameters))
        fasta = context.manifest.annotation.genome_fasta
        if not fasta:
            raise ValueError("SE sequence analysis requires the declared genomic FASTA, not transcript FASTA")
        paths.extend([fasta, Path(str(fasta)+".fai")])
        science = Path(__file__).resolve().parents[1]/"science"
        paths.extend(science/n for n in ("se_rbp.py", "se_rbp_report.py", "public_rbp_resources.py", "motif_significance.py", "statistics.py"))
        paths.append(Path(__file__))
        paths.extend([Path(__file__).resolve().parents[2]/"scripts/fimo_worker.py",
                      Path(__file__).resolve().parents[2]/"katana_cosplice/motif_scanner.py"])
        return list(dict.fromkeys(paths))

    def run(self, context, job, run_fingerprint, parameters):
        source = self.parent_artifacts(context, parameters)[0]
        complete = pd.read_parquet(source.path)
        selected = selection_rows(context, parameters, source, complete)
        if selected.empty or selected.comparison_id.nunique() != 1 or not selected.event_type.eq("SE").all():
            raise ValueError("SE RNA impact requires a nonempty SE-only selection from exactly one comparison")
        complete = complete[complete.comparison_id.eq(selected.comparison_id.iloc[0])].copy()
        motifs, resources = [], []
        species = context.manifest.metadata.get("species", "")
        for asset in context.manifest.assets:
            if asset.kind == "rbp_motif_database":
                verified_resource(asset, species=species)
                target = asset.metadata.get("species", species)
                loaded = load_motif_database(asset.path, organism=target)
                motifs.extend(loaded)
                resources.append({"asset_id": asset.asset_id, "path": str(asset.path), **asset.metadata})
        if parameters["motif_evidence"] == "direct_experiment":
            motifs = [m for m in motifs if m.evidence_class == "direct_experiment"]
        if not motifs:
            raise ValueError("No numeric species-matched motifs available")
        runtime = next(a.path for a in context.manifest.assets if a.kind == "fimo_runtime")
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version,
            common_metadata={"comparison_id": str(selected.comparison_id.iloc[0]), "genome_build": context.manifest.annotation.genome_build,
                             "resources": resources, "tool_name": "Katana MotifTrace", "quantification_rerun": False, "limitations": self.definition.limitations})
        # The reusable scanner writes only to this newly fingerprinted artifact
        # directory. It never rewrites the source selection or old RBP results.
        directory = writer.directory/"sequence_evidence"
        directory.mkdir(parents=True, exist_ok=True)
        inputs = self.input_paths(context, parameters)
        frozen = {"contrast": str(selected.comparison_id.iloc[0]), "genome_build": context.manifest.annotation.genome_build, "species": species,
                  "numerator": str(selected.numerator.iloc[0]) if "numerator" in selected else "numerator",
                  "denominator": str(selected.denominator.iloc[0]) if "denominator" in selected else "denominator",
                  "genome_fasta": str(context.manifest.annotation.genome_fasta),
                  "policy": {**POLICY, "site_pvalue": parameters["site_pvalue"]}, "databases": resources,
                  "selected_event_ids": sorted(selected.event_id), "parameters": parameters,
                  "input_files": [{"path": str(p), "sha256": digest(p), "bytes": p.stat().st_size} for p in inputs]}
        (directory/"ANALYSIS_FREEZE.json").write_text(json.dumps(frozen, indent=2), encoding="utf-8")
        selected.to_csv(directory/"frozen_selected_DAS_SE.tsv", sep="\t", index=False)
        write_coverage(directory)
        run_se_rbp(complete, set(selected.event_id), motifs, context.manifest.annotation.genome_fasta,
                   runtime, directory, threshold=parameters["site_pvalue"],
                   progress=lambda message: (job.check_cancelled(), job.report(.35, message)))
        job.report(.8, "Rendering SE RNA maps/logos and reloading every motif genomic block")
        build_report(directory)
        _brand_new_report(directory)
        for name in ("matching", "sequence_audit", "motif_catalog", "enrichment", "sites"):
            writer.table("Katana MotifTrace: "+name, name+".parquet", pd.read_parquet(directory/(name+".parquet")),
                         entity_type="motif_hit" if name == "sites" else "rbp_evidence", parents=(source.artifact_id,))
        for path in sorted(directory.rglob("*")):
            if path.is_file() and path.suffix != ".parquet":
                writer.existing_file("Katana MotifTrace: "+path.name, path, entity_type="rbp_sequence_evidence", parents=(source.artifact_id,))
        return writer.finish()
