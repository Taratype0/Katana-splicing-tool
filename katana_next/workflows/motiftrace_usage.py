"""Explicit integration of pre-existing DTU, DEXSeq, expression and RNA sequence."""
import json
import re
import shutil
from pathlib import Path
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from katana_next.science import motiftrace_usage as science
from katana_next.workflows.base import ArtifactWriter, ParameterDefinition, WorkflowBase, WorkflowDefinition
from katana_next.workflows.sources import resolve_table


def _contrast(context, parameters):
    ids = parameters.get("contrast_ids") or []
    if len(ids) != 1:
        raise ValueError("Choose exactly one existing comparison; conditions are never pooled")
    return next(c for c in context.manifest.contrasts if c.contrast_id == ids[0])


def _asset(context, kind, contrast):
    found = [a for a in context.manifest.assets if a.kind == kind and a.metadata.get("contrast_id") == contrast.contrast_id]
    if len(found) != 1:
        raise ValueError("Exactly one matched original "+kind+" result is required")
    return found[0]


def _file(path, full=True):
    p = Path(path)
    stat = p.stat()
    return {"path": str(p), "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns,
            "SHA256": science.sha256(p) if full else "deferred_until_sequence_validation"}


def prepare_existing_project_case(context, parameters, directory, job):
    """Read existing tables/quant.sf and validate sequence; never launch Salmon/R."""
    directory = Path(directory)
    directory.mkdir()
    contrast = _contrast(context, parameters)
    annotation = context.manifest.annotation
    if not all((annotation.gtf, annotation.genome_fasta, annotation.transcript_fasta, annotation.genome_build)):
        raise ValueError("Matched GTF, genomic FASTA/build and transcript FASTA are required")
    source = resolve_table(context, "prioritize.isoform_support", "Sample transcript abundance", parameters)
    table = Path(parameters.get("usage_table_path") or (source.path if source else ""))
    if not table.is_file():
        raise ValueError("Provide existing Sample transcript abundance/TPM table; no implicit quantification")
    original = {kind: _asset(context, kind, contrast) for kind in ("dtu_result", "dexseq_result", "deg_result")}
    sample_ids = list(contrast.numerator.sample_ids) + list(contrast.denominator.sample_ids)
    sample_assets = {s.sample_id: s.assets.get("quant_sf") for s in context.manifest.samples if s.sample_id in sample_ids}
    if any(not sample_assets.get(s) for s in sample_ids):
        raise ValueError("Original selected quant.sf files are required to verify the cached values")
    se_assets = [a for a in context.manifest.assets if a.kind == "rmats_result" and a.metadata.get("contrast_id") == contrast.contrast_id and a.metadata.get("event_type") == "SE"]
    if len(se_assets) != 1 or se_assets[0].metadata.get("reverse_direction", False):
        raise ValueError("This preparation currently requires one explicitly forward-mapped SE source")
    bundle = Path(parameters["motif_bundle_path"])
    paths = [table, *(a.path for a in original.values()), se_assets[0].path, *sample_assets.values(),
             *(bundle/n for n in ("UNIQUE_TESTED_MODELS.tsv", "NATIVE_CONSENSUS_PATTERNS.json", "NATIVE_NUMERIC_MATRICES.npz"))]
    policy = {**science.POLICY, "UTC_before_input_selection_or_motif_scores": datetime.now(timezone.utc).isoformat(),
              "comparison_id": contrast.contrast_id, "species": context.manifest.metadata.get("species", "unknown"),
              "genome_build": annotation.genome_build,
              "numerator": contrast.numerator.model_dump(mode="json"), "denominator": contrast.denominator.model_dump(mode="json"),
              "additional_case_eligibility": "original SE FDR<=0.05 and abs(dPSI)>=0.1 target exon overlaps a significant same-strand DEXSeq bin; no motif inspection",
              "models": len(pd.read_csv(bundle/"UNIQUE_TESTED_MODELS.tsv", sep="\t")),
              "source_inputs": [_file(p) for p in dict.fromkeys(paths)],
              "source_statistics_preserved": True, "quantification_rerun": False}
    science.write_json(directory/"PLAN_BEFORE_CASE_SELECTION_AND_SCANNING.json", policy)
    usage = science.usage_fractions(science.read_table(table), contrast.numerator.sample_ids, contrast.denominator.sample_ids)
    dtu = science.canonical_dtu(science.read_table(original["dtu_result"].path))
    dexseq = science.canonical_dexseq(science.read_table(original["dexseq_result"].path), contrast.contrast_id)
    se = pd.read_csv(se_assets[0].path, sep="\t", low_memory=False)
    se["gene_id"] = se.GeneID.astype(str).str.strip('"').map(science.normalized_id)
    se["fdr"], se["dpsi"] = pd.to_numeric(se.FDR, errors="coerce"), pd.to_numeric(se.IncLevelDifference, errors="coerce")
    eligible = set()
    for row in se[se.fdr.le(.05) & se.dpsi.abs().ge(.1)].to_dict("records"):
        matched = dexseq[dexseq.gene_id.eq(row["gene_id"]) & dexseq.unambiguous_gene & dexseq.DEXSeq_bin_FDR.le(.05)]
        if any(b["genomicData.seqnames"] == row["chr"] and b["genomicData.strand"] == row["strand"] and min(row["exonEnd"], b["end0"]) > max(row["exonStart_0base"], b["start0"]) for b in matched.to_dict("records")):
            eligible.add(row["gene_id"])
    case, selection = science.select_case(usage, dtu, dexseq, requested_gene=parameters.get("gene_id", ""), eligible_genes=eligible)
    science.write_json(directory/"CASE_SELECTED_BEFORE_MOTIF_SCANNING.json", case)
    selection.to_csv(directory/"ALL_CANDIDATE_DENOMINATORS_AND_SELECTION.tsv.gz", sep="\t", index=False)
    gene_usage = usage[usage.gene_id.eq(case["gene_id"])].copy()
    checks = []
    for sample in sample_ids:
        quant = pd.read_csv(sample_assets[sample], sep="\t")
        quant["transcript_id"] = quant.Name.map(science.normalized_id)
        merged = gene_usage[gene_usage.sample_id.eq(sample)].merge(quant, on="transcript_id", how="left", validate="one_to_one")
        if not np.allclose(merged.tpm, merged.TPM, rtol=0, atol=1e-9):
            raise ValueError("Cached transcript TPM differs from original quant.sf")
        checks.append({"sample_id": sample, "original_TPM_equal": True, "transcripts": len(merged)})
    science.write_json(directory/"ORIGINAL_QUANT_VALUES_RELOADED_WITHOUT_REQUANTIFICATION.json", checks)
    gene_usage.to_parquet(directory/"EXISTING_GENE_TRANSCRIPT_USAGE.parquet", index=False)
    gene_usage.to_csv(directory/"EXISTING_GENE_TRANSCRIPT_USAGE.tsv", sep="\t", index=False)
    dtu[dtu.gene_id.eq(case["gene_id"])].to_csv(directory/"ORIGINAL_STAGER_GENE_AND_TRANSCRIPT_Q.tsv", sep="\t", index=False)
    bins = dexseq[dexseq.gene_id.eq(case["gene_id"]) & dexseq.unambiguous_gene]
    bins.to_csv(directory/"ORIGINAL_DEXSEQ_BINS_P_FDR_AND_COUNTS.tsv", sep="\t", index=False)
    deg = science.read_table(original["deg_result"].path)
    deg[deg.gene_id.map(science.normalized_id).eq(case["gene_id"])].to_csv(directory/"ORIGINAL_DEG_P_FDR_AND_EFFECT.tsv", sep="\t", index=False)
    se[se.gene_id.eq(case["gene_id"])].to_csv(directory/"ORIGINAL_SE_P_FDR_AND_PSI.tsv", sep="\t", index=False)
    job.report(.12, "Validating full annotated RNA against both genome/GTF and transcript FASTA")
    rna = science.reconstruct_transcripts(case["gene_id"], set(gene_usage.transcript_id), annotation.gtf,
                                        annotation.genome_fasta, annotation.transcript_fasta, job.check_cancelled)
    rna.to_parquet(directory/"REFERENCE_TRANSCRIPT_RNA.parquet", index=False)
    rna.drop(columns=["sequence", "sequence_RNA"]).to_csv(directory/"TRANSCRIPT_COORDINATES_AND_SEQUENCE_VALIDATION.tsv", sep="\t", index=False)
    with (directory/"REAL_REFERENCE_TRANSCRIPTS_RNA.fasta").open("w", encoding="ascii") as handle:
        for row in rna.itertuples():
            handle.write(f">{row.transcript_id}\n{row.sequence_RNA}\n")
    sample_bins, mapping = science.dexseq_bin_support(bins, gene_usage, rna)
    sample_bins.to_csv(directory/"DEXSEQ_PER_SAMPLE_COUNTS_AND_WITHIN_GENE_FRACTIONS.tsv", sep="\t", index=False)
    mapping.to_csv(directory/"DEXSEQ_BIN_TO_TRANSCRIPT_COORDINATE_MAP.tsv", sep="\t", index=False)
    science.write_json(directory/"PREPARATION_COMPLETE.json", {"case": case, "samples": sample_ids,
        "transcripts": len(rna), "reference_sources": [_file(p) for p in (annotation.gtf, annotation.genome_fasta, annotation.transcript_fasta)],
        "plan_SHA256": science.sha256(directory/"PLAN_BEFORE_CASE_SELECTION_AND_SCANNING.json")})
    science.seal_prepared_case(directory)
    return directory


def receive_completed_case(source, output, prepared, bundle):
    """Verify and register an already completed scan; never rescore it."""
    source, output, prepared, bundle = map(Path, (source, output, prepared, bundle))
    manifest = json.loads((source/"REPORT_MANIFEST.json").read_text())
    complete = json.loads((source/"USAGE_MOTIF_COMPLETE.json").read_text())
    freeze = json.loads((source/"ANALYSIS_BEFORE_MODEL_SCORING.json").read_text())
    if complete["prepared_manifest_SHA256"] != science.sha256(prepared/"PREPARED_CASE_MANIFEST.json"):
        raise ValueError("Completed scan belongs to different prepared inputs")
    if freeze["execution_implementation_SHA256"] != science.sha256(science.__file__):
        raise ValueError("Completed scan's numerical engine is not the current verified implementation")
    for name, sha in freeze["model_source_hashes"].items():
        if science.sha256(bundle/name) != sha:
            raise ValueError("Completed scan used a different motif library")
    for name, record in manifest["files"].items():
        if Path(name).name != name:
            raise ValueError("Result manifest must contain bounded local filenames")
        p = source/name
        if not p.is_file() or p.stat().st_size != record["bytes"] or science.sha256(p) != record["SHA256"]:
            raise ValueError("Completed result content changed: "+name)
    if output.exists():
        raise ValueError("Existing import destination preserved")
    output.mkdir()
    for name in [*manifest["files"], "REPORT_MANIFEST.json"]:
        shutil.copyfile(source/name, output/name)
    science.write_json(output/"VERIFIED_RESULT_REUSE.json", {"completed_source": str(source),
        "report_manifest_SHA256": science.sha256(source/"REPORT_MANIFEST.json"), "models_rescored": 0,
        "all_result_hashes_verified": True, "quantification_rerun": False})
    return complete


def validate_prepared_context(context, contrast, plan, complete):
    if plan["comparison_id"] != contrast.contrast_id:
        raise ValueError("Prepared comparison mismatch")
    for arm in ("numerator", "denominator"):
        if set(plan[arm]["sample_ids"]) != set(getattr(contrast, arm).sample_ids):
            raise ValueError("Prepared comparison arm/sample identity mismatch: "+arm)
    if set(complete["samples"]) != set((*contrast.numerator.sample_ids, *contrast.denominator.sample_ids)):
        raise ValueError("Prepared sample universe mismatch")
    if plan["genome_build"] != context.manifest.annotation.genome_build or plan["species"] != context.manifest.metadata.get("species"):
        raise ValueError("Prepared reference/species mismatch")
    reference_paths = {Path(row["path"]).resolve() for row in complete["reference_sources"]}
    if any(Path(p).resolve() not in reference_paths for p in (context.manifest.annotation.gtf, context.manifest.annotation.genome_fasta, context.manifest.annotation.transcript_fasta)):
        raise ValueError("Prepared reference annotation identity mismatch")


class MotifTraceUsageWorkflow(WorkflowBase):
    version = "1.0.0+english.20261021"
    definition = WorkflowDefinition(
        workflow_id="prioritize.motiftrace_usage", workspace_id="prioritize",
        label="Katana MotifTrace: DTU, exon usage and expression to RNA motifs",
        question="Which actual motif loci become more or less available as RNA composition changes, separately from gene abundance?",
        method="Reuse existing Salmon values, stageR, coordinate-resolved DEXSeq, DEG and DAS; validate genome/GTF/transcript FASTA equality; scan all supplied motifs; weight shared genomic/junction loci by all gene isoform fractions. No new single-case motif P/FDR.",
        analysis_unit="annotated gene isoforms and independent biological samples within one comparison",
        output_summary="Real RNA FASTA, split motif coordinates, per-sample transcript fractions/expression/exon counts, locus gain/loss, separate original P/FDR and an HTML case report.",
        supports_contrast_selection=True,
        parameters=(
            ParameterDefinition("usage_table_path", "Existing transcript TPM table (optional)", "text", "", "Otherwise reuse the Sample transcript abundance artifact; no implicit upstream execution."),
            ParameterDefinition("motif_bundle_path", "Verified complete motif-model bundle", "text", "", "Directory containing UNIQUE_TESTED_MODELS.tsv, NATIVE_CONSENSUS_PATTERNS.json and NATIVE_NUMERIC_MATRICES.npz; every model is retained."),
            ParameterDefinition("gene_id", "Gene ID (optional)", "text", "", "Otherwise use the prespecified motif-blind DTU/DEXSeq/DAS coordinate-supported ranking."),
            ParameterDefinition("prepared_case_path", "Verified existing prepared case (optional)", "text", "", "Reuse an immutable prior preparation after source/hash/contrast checks, without reconstructing it again."),
            ParameterDefinition("completed_analysis_path", "Receive a completed MotifTrace usage run (optional)", "text", "", "Validate all result/model/source hashes and register existing output; never rescore completed models."),
        ),
        presentation_role="advanced",
        limitations="Sequence/usage association, not measured binding, motif enrichment or causal regulation. The increasing transcript need not itself pass DTU; its own q remains separate. DEXSeq bins are not whole-exon PSI. Reference RNA is not sample-specific allelic sequence. CLIP is unavailable unless independently supplied.",
    )

    def parent_artifacts(self, context, parameters):
        source = resolve_table(context, "prioritize.isoform_support", "Sample transcript abundance", parameters)
        return [source] if source else []

    def input_paths(self, context, parameters):
        contrast = _contrast(context, parameters)
        bundle = Path(parameters.get("motif_bundle_path", ""))
        paths = [bundle/n for n in ("UNIQUE_TESTED_MODELS.tsv", "NATIVE_CONSENSUS_PATTERNS.json", "NATIVE_NUMERIC_MATRICES.npz")]
        if parameters.get("prepared_case_path"):
            paths.append(Path(parameters["prepared_case_path"])/"PREPARED_CASE_MANIFEST.json")
        else:
            paths.extend(a.path for a in self.parent_artifacts(context, parameters))
            if parameters.get("usage_table_path"):
                paths.append(Path(parameters["usage_table_path"]))
            paths.extend(_asset(context, k, contrast).path for k in ("dtu_result", "dexseq_result", "deg_result"))
            paths.extend(p for p in (context.manifest.annotation.gtf, context.manifest.annotation.genome_fasta, context.manifest.annotation.transcript_fasta) if p)
            paths.extend(s.assets["quant_sf"] for s in context.manifest.samples if s.sample_id in set((*contrast.numerator.sample_ids, *contrast.denominator.sample_ids)) and "quant_sf" in s.assets)
        if parameters.get("completed_analysis_path"):
            paths.append(Path(parameters["completed_analysis_path"])/"REPORT_MANIFEST.json")
        return [*paths, Path(__file__), Path(science.__file__)]

    def run(self, context, job, run_fingerprint, parameters):
        contrast = _contrast(context, parameters)
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version,
            common_metadata={"tool_name": "Katana MotifTrace", "quantification_rerun": False,
                             "new_motif_P_FDR": "not_computed_single_gene_descriptive_case", "limitations": self.definition.limitations})
        prepared = Path(parameters["prepared_case_path"]) if parameters.get("prepared_case_path") else prepare_existing_project_case(context, parameters, writer.directory/"prepared", job)
        plan, complete = science.verify_prepared_case(prepared)
        validate_prepared_context(context, contrast, plan, complete)
        output = writer.directory/"usage_sequence"
        if parameters.get("completed_analysis_path"):
            result = receive_completed_case(parameters["completed_analysis_path"], output, prepared, parameters["motif_bundle_path"])
            writer.common_metadata["verified_completed_scan_reuse"] = True
        else:
            result = science.analyze_prepared_case(prepared, output, parameters["motif_bundle_path"],
                check_cancelled=job.check_cancelled, progress=lambda n, total: job.report(.15+.7*n/total, f"MotifTrace actual models complete: {n}/{total}"))
            from katana_next.science.motiftrace_usage_report import build_report
            build_report(prepared, output)
        parents = tuple(a.artifact_id for a in self.parent_artifacts(context, parameters))
        urls = {}
        for path in sorted(output.iterdir()):
            if path.is_file():
                artifact = writer.existing_file("MotifTrace RNA usage: "+path.name, path, entity_type="motiftrace_usage_evidence", parents=parents,
                                                metadata={"mime_type": "text/html"} if path.suffix == ".html" else None)
                urls[path.name] = f"/api/projects/{context.manifest.project_id}/artifacts/{artifact.artifact_id}/file?download=true"
        portable = (output/"INDEX.html").read_text(encoding="utf-8")
        app_html = re.sub(r'(href|src)="([^"]+)"', lambda match: match.group(1)+'="'+urls.get(match.group(2), match.group(2))+'"', portable)
        app_path = writer.directory/"MOTIFTRACE_USAGE_APP.html"
        app_path.write_text(app_html, encoding="utf-8")
        writer.existing_file("Katana MotifTrace: RNA usage, actual sequences and independent evidence", app_path,
                             entity_type="motiftrace_usage_report", parents=parents, metadata={"mime_type": "text/html"})
        writer.common_metadata.update(case=result["case"], models=result["models"], transcripts=result["transcripts"])
        return writer.finish()
