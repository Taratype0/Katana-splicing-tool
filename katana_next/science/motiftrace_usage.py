"""Trace existing RNA usage into motif availability, without re-quantification.

All coordinates are BED-like 0-based half-open. Motif opportunities are
sequence predictions, not binding probabilities, and receive no invented P.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from katana_cosplice.motif_scanner import IndexedFasta
from katana_next.science.isoform_consequences import load_models, load_sequences, normalized_id
from katana_next.science.motiftrace_events import coalesce, score_positions, validate_models
from katana_next.science.se_rbp import _chromosome, mapped_blocks, sequence_from_blocks


POLICY = {
    "version": "1.0.0",
    "contrast": "one declared numerator-minus-denominator comparison",
    "gene_expression_eligibility": "gene TPM >=1 in every selected biological sample",
    "candidate_rule": "at least one stageR transcript confirmation q<=0.05 and one unambiguous DEXSeq bin padj<=0.05; select largest absolute transcript-fraction difference among confirmed transcripts, ties gene then transcript ID; never inspect motif scores",
    "pair_rule": "largest increase versus largest decrease in mean per-sample transcript fraction; ties transcript ID; all quantified annotated gene isoforms remain in weighting",
    "library": "all supplied models, no truncation or model selection by result",
    "consensus": "exact finite IUPAC pattern, all overlapping starts",
    "PWM_display_candidates": "existing MotifTrace relative likelihood product >=0.8 after (20*PWM+0.1*background)/20.1; not a calibrated site P",
    "base_background": [0.25, 0.25, 0.25, 0.25],
    "new_motif_P_FDR": "not_computed_single_gene_descriptive_case; no population-enrichment or occupancy claim",
    "source_statistics": "DAS, stageR gene screen, stageR transcript confirmation, DEXSeq bin and DEG remain separate; no combined P/FDR",
    "decoy_or_enrichment_background": "not_used_for_single_gene_sequence_consequences; existing backgrounds and historical results untouched",
    "missing_values": "preserved; no missing transcript, count, proportion or reference zero filling",
}


def sha256(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def read_table(path):
    path = Path(path)
    return pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path, sep="\t", low_memory=False)


def usage_fractions(frame, numerator, denominator):
    required = {"transcript_id", "gene_id", "sample_id", "tpm"}
    if not required <= set(frame):
        raise ValueError("Usage needs transcript_id/gene_id/sample_id/tpm, not only DTU P-values")
    numerator, denominator = list(numerator), list(denominator)
    if not numerator or not denominator or set(numerator) & set(denominator):
        raise ValueError("Distinct, explicit biological sample arms required")
    if len(set(numerator + denominator)) != len(numerator + denominator):
        raise ValueError("Duplicate sample identity")
    selected = numerator + denominator
    result = frame[frame.sample_id.isin(selected)].copy()
    if set(result.sample_id) != set(selected):
        raise ValueError("Missing selected sample; never replace it by zero")
    for column in ("transcript_id", "gene_id"):
        if result[column].isna().any():
            raise ValueError("Missing transcript/gene mapping")
        result[column] = result[column].map(normalized_id)
    if result.duplicated(["sample_id", "transcript_id"]).any():
        raise ValueError("Duplicate normalized transcript in a biological sample")
    if result.groupby("transcript_id").gene_id.nunique().gt(1).any():
        raise ValueError("Ambiguous transcript to gene mapping")
    catalogs = result.groupby("sample_id").transcript_id.agg(set)
    if any(value != catalogs.iloc[0] for value in catalogs):
        raise ValueError("Unequal transcript catalogs; no missing-transcript zero filling")
    result["tpm"] = pd.to_numeric(result.tpm, errors="raise")
    if not np.isfinite(result.tpm).all() or result.tpm.lt(0).any():
        raise ValueError("Invalid TPM")
    result["gene_TPM"] = result.groupby(["sample_id", "gene_id"]).tpm.transform("sum")
    result["transcript_fraction"] = result.tpm / result.gene_TPM.where(result.gene_TPM.gt(0))
    result["arm"] = np.where(result.sample_id.isin(numerator), "numerator", "denominator")
    return result


def canonical_dtu(frame):
    if not {"geneID", "txID", "gene", "transcript"} <= set(frame):
        raise ValueError("This adapter requires native stageR gene/transcript confirmation columns")
    result = frame.copy()
    result["gene_id"] = result.geneID.map(normalized_id)
    result["transcript_id"] = result.txID.map(normalized_id)
    result["stageR_gene_screen_q"] = pd.to_numeric(result.gene, errors="coerce")
    result["stageR_transcript_confirmation_q"] = pd.to_numeric(result.transcript, errors="coerce")
    if result.duplicated(["gene_id", "transcript_id"]).any():
        raise ValueError("Duplicate DTU transcript result")
    return result


def canonical_dexseq(frame, comparison_id):
    needed = {"groupID", "featureID", "pvalue", "padj", "genomicData.seqnames", "genomicData.start", "genomicData.end", "genomicData.strand"}
    if not needed <= set(frame):
        raise ValueError("DEXSeq needs native bin coordinates, strand and original P/padj")
    result = frame.copy()
    if "contrast_id" in result and not result.contrast_id.eq(comparison_id).all():
        raise ValueError("DEXSeq contrast mismatch")
    # Multi-gene flattened groups are retained as ambiguous, never assigned to one gene.
    result["unambiguous_gene"] = ~result.groupID.astype(str).str.contains("+", regex=False)
    result["gene_id"] = result.groupID.map(normalized_id)
    result["start0"] = pd.to_numeric(result["genomicData.start"], errors="raise") - 1
    result["end0"] = pd.to_numeric(result["genomicData.end"], errors="raise")
    if ((result.start0 < 0) | (result.end0 <= result.start0) | (result.start0 % 1 != 0) | (result.end0 % 1 != 0)).any():
        raise ValueError("Invalid one-based inclusive DEXSeq coordinates")
    result["start0"] = result.start0.astype(int)
    result["end0"] = result.end0.astype(int)
    result["DEXSeq_P"] = pd.to_numeric(result.pvalue, errors="coerce")
    result["DEXSeq_bin_FDR"] = pd.to_numeric(result.padj, errors="coerce")
    return result


def select_case(usage, dtu, dexseq, *, requested_gene="", eligible_genes=None):
    means = usage.groupby(["gene_id", "transcript_id", "arm"]).transcript_fraction.mean().unstack("arm")
    means["delta_fraction"] = means.numerator - means.denominator
    means = means.reset_index().merge(dtu[["gene_id", "transcript_id", "stageR_gene_screen_q", "stageR_transcript_confirmation_q"]], on=["gene_id", "transcript_id"], how="left", validate="one_to_one")
    eligibility = usage.groupby(["gene_id", "sample_id"]).gene_TPM.first().groupby("gene_id").min().ge(1)
    dex_genes = set(dexseq.loc[dexseq.unambiguous_gene & dexseq.DEXSeq_bin_FDR.le(.05), "gene_id"])
    means["eligible_gene_expression"] = means.gene_id.map(eligibility).fillna(False)
    means["DEXSeq_supported_gene"] = means.gene_id.isin(dex_genes)
    means["candidate_eligible"] = means.eligible_gene_expression & means.DEXSeq_supported_gene & means.stageR_transcript_confirmation_q.le(.05) & means.delta_fraction.notna()
    if eligible_genes is not None:
        means["candidate_eligible"] &= means.gene_id.isin(eligible_genes)
    candidates = means[means.candidate_eligible].assign(abs_delta=lambda x: x.delta_fraction.abs()).sort_values(["abs_delta", "gene_id", "transcript_id"], ascending=[False, True, True])
    if requested_gene:
        candidates = candidates[candidates.gene_id.eq(normalized_id(requested_gene))]
    if candidates.empty:
        raise ValueError("No predeclared eligible DTU/DEXSeq/expression case; do not change selection after motif inspection")
    gene = candidates.iloc[0].gene_id
    selected = means[means.gene_id.eq(gene)]
    up = selected.sort_values(["delta_fraction", "transcript_id"], ascending=[False, True]).iloc[0]
    down = selected.sort_values(["delta_fraction", "transcript_id"], ascending=[True, True]).iloc[0]
    if up.delta_fraction <= 0 or down.delta_fraction >= 0 or up.transcript_id == down.transcript_id:
        raise ValueError("No opposing transcript-fraction change")
    return {"gene_id": gene, "increased_transcript": up.transcript_id, "decreased_transcript": down.transcript_id,
            "increased_delta_fraction": float(up.delta_fraction), "decreased_delta_fraction": float(down.delta_fraction),
            "selection_rank": "largest absolute confirmed transcript-fraction change, motif-blind"}, means


def reconstruct_transcripts(gene, transcript_ids, gtf, genome_fasta, transcript_fasta, check_cancelled=lambda: None):
    requested = set(transcript_ids)
    models = load_models(Path(gtf), requested, check_cancelled)
    if set(models) != requested or any(m.gene_id != gene for m in models.values()):
        raise ValueError("Incomplete/discordant GTF transcript to gene mapping")
    reference = load_sequences(Path(transcript_fasta), requested, check_cancelled, expected_models=models)
    if set(reference) != requested:
        raise ValueError("Incomplete declared transcript FASTA; no reference repair")
    fasta = IndexedFasta(Path(genome_fasta))
    rows = []
    with fasta.fasta.open("rb") as handle:
        for tx in sorted(requested):
            check_cancelled()
            model = models[tx]
            blocks = coalesce(model.exons, model.strand)
            if not blocks:
                raise ValueError("Transcript has no exons")
            chrom = _chromosome(fasta, model.chrom)
            sequence = sequence_from_blocks(fasta, handle, chrom, model.strand, blocks)
            if sequence != reference[tx]:
                raise ValueError(f"Genome+GTF reconstructed RNA differs from declared transcript FASTA: {tx}")
            if set(sequence) - set("ACGT"):
                raise ValueError(f"Ambiguous sequence retained as failure, not zero motif: {tx}")
            rows.append({"transcript_id": tx, "annotation_transcript_id": model.transcript_id, "gene_id": gene,
                         "chr": chrom, "strand": model.strand, "blocks_json": json.dumps(blocks),
                         "junctions_json": json.dumps(np.cumsum([b-a for a, b in blocks])[:-1].tolist()),
                         "sequence": sequence, "sequence_RNA": sequence.replace("T", "U"), "length_nt": len(sequence),
                         "sequence_SHA256": hashlib.sha256(sequence.encode()).hexdigest(),
                         "reference_validation": "genome_GTF_and_independent_transcript_FASTA_exact_equal"})
    return pd.DataFrame(rows)


def read_model_bundle(directory):
    directory = Path(directory)
    catalog = pd.read_csv(directory / "UNIQUE_TESTED_MODELS.tsv", sep="\t", keep_default_na=False)
    patterns = json.loads((directory / "NATIVE_CONSENSUS_PATTERNS.json").read_text())
    models = {}
    with np.load(directory / "NATIVE_NUMERIC_MATRICES.npz", allow_pickle=False) as matrices:
        for row in catalog.sort_values("model_id").itertuples():
            token = row.model_id
            if token in models:
                raise ValueError("Duplicate motif model")
            model = {"kind": row.kind, "width": int(row.width), "RBPs": row.RBPs}
            if row.kind == "consensus":
                model["pattern"] = patterns[token]
                content_hash = hashlib.sha256(model["pattern"].encode()).hexdigest()
            elif row.kind == "PWM":
                matrix = matrices[token]
                model["matrix"] = matrix.tolist()
                content_hash = hashlib.sha256(matrix.astype("<f8").tobytes()).hexdigest()
            else:
                raise ValueError("Unknown motif kind")
            if token[2:] != content_hash[:20]:
                raise ValueError("Native model content hash mismatch")
            validate_models({token: model})
            models[token] = model
    return catalog, models


HIT_COLUMNS = ["transcript_id", "gene_id", "model_id", "RBP", "model_kind", "RNA_start0", "RNA_end0", "chr", "strand", "genomic_blocks_json", "hit_RNA", "crosses_splice_join", "score", "site_key", "evidence"]


def candidate_positions(sequence, model):
    """Same existing candidate definition, vectorized for verified A/C/G/T RNA."""
    if set(sequence) - set("ACGT"):
        scores, candidate, _ = score_positions(sequence, model, background=POLICY["base_background"])
        starts = np.flatnonzero(candidate)
        return starts, scores[starts], int(np.isfinite(scores).sum())
    width = model["width"]
    n = max(0, len(sequence)-width+1)
    if model["kind"] == "consensus":
        starts = np.asarray([m.start(1) for m in re.finditer("(?=("+model["pattern"]+"))", sequence)], dtype=int)
        return starts, np.ones(len(starts)), n
    if not n:
        return np.empty(0, dtype=int), np.empty(0), 0
    bg = np.asarray(POLICY["base_background"])
    matrix = (20*np.asarray(model["matrix"], float)+.1*bg)/20.1
    relative, odds = matrix/matrix.max(axis=1, keepdims=True), matrix/bg
    lookup = np.full(256, 4, dtype=np.uint8)
    lookup[np.frombuffer(b"ACGT", dtype=np.uint8)] = np.arange(4)
    encoded = lookup[np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)]
    rel, values = np.ones(n), np.ones(n)
    for i in range(width):
        rel *= relative[i, encoded[i:i+n]]
        values *= odds[i, encoded[i:i+n]]
    starts = np.flatnonzero(rel >= .8)
    return starts, values[starts], n


def scan_transcript_models(rna, models, check_cancelled=lambda: None, progress=lambda *_: None):
    hits, counts = [], []
    forms = rna.to_dict("records")
    for index, (token, model) in enumerate(models.items()):
        check_cancelled()
        validate_models({token: model})
        for form in forms:
            starts, scores, valid = candidate_positions(form["sequence"], model)
            width = model["width"]
            counts.append({"transcript_id": form["transcript_id"], "model_id": token, "RBP": model.get("RBPs", ""),
                           "valid_starts": valid, "candidate_instances": len(starts) if valid else np.nan,
                           "status": "scanned" if valid else "sequence_shorter_than_motif_not_evaluable"})
            for hit_index, start in enumerate(starts):
                start = int(start)
                blocks = coalesce(mapped_blocks(json.loads(form["blocks_json"]), form["strand"], start, start+width), form["strand"])
                encoded = json.dumps(blocks)
                site_key = token + ":" + form["chr"] + ":" + form["strand"] + ":" + encoded
                hits.append({"transcript_id": form["transcript_id"], "gene_id": form["gene_id"], "model_id": token,
                             "RBP": model.get("RBPs", ""), "model_kind": model["kind"], "RNA_start0": start, "RNA_end0": start+width,
                             "chr": form["chr"], "strand": form["strand"], "genomic_blocks_json": encoded,
                             "hit_RNA": form["sequence_RNA"][start:start+width],
                             "crosses_splice_join": any(start < j < start+width for j in json.loads(form["junctions_json"])),
                             "score": float(scores[hit_index]), "site_key": site_key,
                             "evidence": "finite_pattern_match" if model["kind"] == "consensus" else "relative_likelihood_candidate_not_site_P"})
        if (index+1) % 100 == 0 or index+1 == len(models):
            progress(index+1, len(models))
    return pd.DataFrame(hits, columns=HIT_COLUMNS), pd.DataFrame(counts)


def weighted_model_opportunities(usage, counts):
    if usage.gene_id.nunique() != 1:
        raise ValueError("One complete gene transcript universe is required")
    frame = usage.merge(counts, on="transcript_id", how="left", validate="many_to_many")
    if frame.model_id.isna().any():
        raise ValueError("Unscanned transcript; no missing-score zero filling")
    frame["fraction_weighted_instances"] = frame.transcript_fraction * frame.candidate_instances
    frame["TPM_weighted_instances"] = frame.tpm * frame.candidate_instances
    rows = []
    for (sample, arm, token), group in frame.groupby(["sample_id", "arm", "model_id"], sort=False):
        complete = group.candidate_instances.notna().all() and group.transcript_fraction.notna().all()
        rows.append({"sample_id": sample, "arm": arm, "model_id": token, "RBP": group.RBP.iloc[0],
                     "gene_TPM": float(group.gene_TPM.iloc[0]),
                     "instances_per_RNA": float(group.fraction_weighted_instances.sum()) if complete else np.nan,
                     "TPM_weighted_instances": float(group.TPM_weighted_instances.sum()) if complete else np.nan,
                     "status": "descriptive_sequence_opportunity" if complete else "not_evaluable",
                     "motif_P": np.nan, "motif_FDR": np.nan})
    sample_scores = pd.DataFrame(rows)
    summary = []
    for token, group in sample_scores.groupby("model_id", sort=False):
        means = group.groupby("arm")[["gene_TPM", "instances_per_RNA", "TPM_weighted_instances"]].mean()
        if len(means) != 2:
            raise ValueError("Missing comparison arm")
        a1, a0 = means.loc["numerator", "gene_TPM"], means.loc["denominator", "gene_TPM"]
        b1, b0 = means.loc["numerator", "TPM_weighted_instances"], means.loc["denominator", "TPM_weighted_instances"]
        c1, c0 = b1/a1 if a1 > 0 else np.nan, b0/a0 if a0 > 0 else np.nan
        composition, abundance = (a1+a0)/2*(c1-c0), (c1+c0)/2*(a1-a0)
        summary.append({"model_id": token, "RBP": group.RBP.iloc[0],
                        "delta_mean_instances_per_RNA": means.loc["numerator", "instances_per_RNA"]-means.loc["denominator", "instances_per_RNA"],
                        "delta_mean_TPM_weighted_instances": b1-b0,
                        "composition_component_TPM_weighted": composition, "gene_abundance_component_TPM_weighted": abundance,
                        "decomposition_residual": b1-b0-composition-abundance,
                        "decomposition_definition": "symmetric exact product difference on arm-mean TPM, not mean individual fractions",
                        "motif_P": np.nan, "motif_FDR": np.nan})
    return sample_scores, pd.DataFrame(summary)


def weighted_locus_opportunities(usage, hits, pair):
    rows, summaries = [], []
    sample_order = usage.sample_id.drop_duplicates().tolist()
    pivot_fraction = usage.pivot(index="transcript_id", columns="sample_id", values="transcript_fraction")
    pivot_tpm = usage.pivot(index="transcript_id", columns="sample_id", values="tpm")
    arms = usage.drop_duplicates("sample_id").set_index("sample_id").arm
    for key, group in hits.groupby("site_key", sort=False):
        if group.hit_RNA.nunique() != 1 or group.duplicated("transcript_id").any():
            raise ValueError("Ambiguous shared motif locus")
        txs = sorted(set(group.transcript_id))
        fractions = pivot_fraction.loc[txs].sum(axis=0, min_count=len(txs))
        abundance = pivot_tpm.loc[txs].sum(axis=0, min_count=len(txs))
        up, down = pair["increased_transcript"] in txs, pair["decreased_transcript"] in txs
        state = "retained_in_pair" if up and down else "increased_form_only" if up else "decreased_form_only" if down else "other_gene_isoforms_only"
        source = group.iloc[0]
        record = {name: source[name] for name in ("site_key", "model_id", "RBP", "model_kind", "chr", "strand", "genomic_blocks_json", "hit_RNA", "evidence")}
        record.update(present_transcript_ids=";".join(txs), pair_state=state,
                      crosses_splice_join=bool(group.crosses_splice_join.any()),
                      delta_fraction_availability=float(fractions[arms[arms.eq("numerator")].index].mean()-fractions[arms[arms.eq("denominator")].index].mean()),
                      delta_TPM_availability=float(abundance[arms[arms.eq("numerator")].index].mean()-abundance[arms[arms.eq("denominator")].index].mean()))
        summaries.append(record)
        for sample in sample_order:
            rows.append({"site_key": key, "model_id": source.model_id, "sample_id": sample, "arm": arms[sample],
                         "RNA_fraction_with_locus": fractions[sample], "transcript_TPM_with_locus": abundance[sample]})
    return pd.DataFrame(rows), pd.DataFrame(summaries)


def dexseq_bin_support(bins, usage, rna):
    samples = usage.sample_id.drop_duplicates().tolist()
    bins = bins.copy()
    for sample in samples:
        if "countData."+sample not in bins:
            raise ValueError("DEXSeq lacks explicit selected sample counts: "+sample)
    rows, mapping = [], []
    for row in bins.to_dict("records"):
        for form in rna.to_dict("records"):
            if row["genomicData.seqnames"] != form["chr"] or row["genomicData.strand"] != form["strand"]:
                overlap = 0
            else:
                overlap = sum(max(0, min(b, row["end0"])-max(a, row["start0"])) for a, b in json.loads(form["blocks_json"]))
            mapping.append({"featureID": row["featureID"], "transcript_id": form["transcript_id"], "overlap_nt": overlap,
                            "bin_width_nt": row["end0"]-row["start0"], "bin_exonic_fraction": overlap/(row["end0"]-row["start0"])})
    mapping = pd.DataFrame(mapping)
    for sample in samples:
        values = pd.to_numeric(bins["countData."+sample], errors="raise")
        if not np.isfinite(values).all() or values.lt(0).any():
            raise ValueError("Invalid native DEXSeq bin counts")
        total = values.sum()
        for i, row in enumerate(bins.to_dict("records")):
            rows.append({"sample_id": sample, "featureID": row["featureID"], "original_bin_count": float(values.iloc[i]),
                         "within_gene_bin_count_fraction": float(values.iloc[i]/total) if total > 0 else np.nan,
                         "meaning": "descriptive fraction of this gene's flattened-bin counts, not PSI or DEXSeq fitted usage"})
    return pd.DataFrame(rows), mapping


def bin_hit_overlap(hits, bins):
    rows = []
    for item in hits.to_dict("records"):
        blocks = json.loads(item["genomic_blocks_json"])
        for row in bins.to_dict("records"):
            if item["chr"] != row["genomicData.seqnames"] or item["strand"] != row["genomicData.strand"]:
                continue
            n = sum(max(0, min(b, row["end0"])-max(a, row["start0"])) for a, b in blocks)
            if n:
                rows.append({"site_key": item["site_key"], "transcript_id": item["transcript_id"], "model_id": item["model_id"],
                             "featureID": row["featureID"], "motif_overlap_nt": n, "motif_width_nt": sum(b-a for a, b in blocks),
                             "DEXSeq_P": row["DEXSeq_P"], "DEXSeq_bin_FDR": row["DEXSeq_bin_FDR"],
                             "meaning": "coordinate overlap with tested bin, not a combined P-value or causal binding"})
    return pd.DataFrame(rows)


def write_json(path, value):
    path = Path(path)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=True, allow_nan=False), encoding="utf-8")


def seal_prepared_case(directory):
    directory = Path(directory)
    marker = directory/"PREPARED_CASE_MANIFEST.json"
    names = ["PLAN_BEFORE_CASE_SELECTION_AND_SCANNING.json", "PREPARATION_COMPLETE.json",
             "CASE_SELECTED_BEFORE_MOTIF_SCANNING.json", "EXISTING_GENE_TRANSCRIPT_USAGE.parquet",
             "REFERENCE_TRANSCRIPT_RNA.parquet", "REAL_REFERENCE_TRANSCRIPTS_RNA.fasta",
             "ORIGINAL_STAGER_GENE_AND_TRANSCRIPT_Q.tsv", "ORIGINAL_DEXSEQ_BINS_P_FDR_AND_COUNTS.tsv",
             "ORIGINAL_DEG_P_FDR_AND_EFFECT.tsv", "ORIGINAL_SE_P_FDR_AND_PSI.tsv",
             "DEXSEQ_PER_SAMPLE_COUNTS_AND_WITHIN_GENE_FRACTIONS.tsv", "DEXSEQ_BIN_TO_TRANSCRIPT_COORDINATE_MAP.tsv",
             "ORIGINAL_QUANT_VALUES_RELOADED_WITHOUT_REQUANTIFICATION.json",
             "TRANSCRIPT_COORDINATES_AND_SEQUENCE_VALIDATION.tsv"]
    value = {"schema": "motiftrace_usage_prepared_case_v1", "files": {
        name: {"bytes": (directory/name).stat().st_size, "SHA256": sha256(directory/name)} for name in names}}
    if marker.exists():
        if json.loads(marker.read_text()) != value:
            raise ValueError("Existing prepared-case seal differs")
    else:
        write_json(marker, value)
    return marker


def verify_prepared_case(directory, *, rehash_originals=True):
    directory = Path(directory)
    seal = json.loads((directory/"PREPARED_CASE_MANIFEST.json").read_text())
    if seal.get("schema") != "motiftrace_usage_prepared_case_v1":
        raise ValueError("Unknown usage-case schema")
    for name, info in seal["files"].items():
        path = directory/name
        if not path.is_file() or path.stat().st_size != info["bytes"] or sha256(path) != info["SHA256"]:
            raise ValueError("Prepared input changed: "+name)
    plan = json.loads((directory/"PLAN_BEFORE_CASE_SELECTION_AND_SCANNING.json").read_text())
    complete = json.loads((directory/"PREPARATION_COMPLETE.json").read_text())
    if sha256(directory/"PLAN_BEFORE_CASE_SELECTION_AND_SCANNING.json") != complete["plan_SHA256"]:
        raise ValueError("Preparation plan identity mismatch")
    for row in plan["source_inputs"]:
        path = Path(row["path"])
        stat = path.stat()
        if stat.st_size != row["bytes"] or stat.st_mtime_ns != row["mtime_ns"]:
            raise ValueError("Original source metadata changed: "+str(path))
        if rehash_originals and sha256(path) != row["SHA256"]:
            raise ValueError("Original source bytes changed: "+str(path))
    for row in complete["reference_sources"]:
        stat = Path(row["path"]).stat()
        if stat.st_size != row["bytes"] or stat.st_mtime_ns != row["mtime_ns"]:
            raise ValueError("Previously validated reference changed")
    return plan, complete


def analyze_prepared_case(directory, output, model_bundle, *, check_cancelled=lambda: None, progress=lambda *_: None):
    """Resumable product engine; each fixed chunk is validated before reuse."""
    directory, output, model_bundle = Path(directory), Path(output), Path(model_bundle)
    plan, complete = verify_prepared_case(directory)
    catalog, models = read_model_bundle(model_bundle)
    if len(models) != plan["models"]:
        raise ValueError("Motif universe differs from the prespecified preparation plan")
    output.mkdir(parents=True, exist_ok=True)
    binding = {"UTC_before_first_model_scan": datetime.now(timezone.utc).isoformat(),
               "prepared_case_manifest_SHA256": sha256(directory/"PREPARED_CASE_MANIFEST.json"),
               "parameters": {key: plan[key] for key in POLICY}, "models": list(models),
               "model_source_hashes": {name: sha256(model_bundle/name) for name in (
                   "UNIQUE_TESTED_MODELS.tsv", "NATIVE_CONSENSUS_PATTERNS.json", "NATIVE_NUMERIC_MATRICES.npz")},
               "execution_implementation_SHA256": sha256(__file__),
               "candidate_definition_source_SHA256": sha256(Path(__file__).with_name("motiftrace_events.py"))}
    freeze = output/"ANALYSIS_BEFORE_MODEL_SCORING.json"
    if freeze.exists():
        original = json.loads(freeze.read_text())
        for key in binding:
            if key != "UTC_before_first_model_scan" and binding[key] != original[key]:
                raise ValueError("Cannot resume changed execution/scientific inputs: "+key)
    else:
        write_json(freeze, binding)
    if (output/"USAGE_MOTIF_COMPLETE.json").exists():
        raise ValueError("Completed analysis is not re-executed")
    usage = pd.read_parquet(directory/"EXISTING_GENE_TRANSCRIPT_USAGE.parquet")
    rna = pd.read_parquet(directory/"REFERENCE_TRANSCRIPT_RNA.parquet")
    bins = pd.read_csv(directory/"ORIGINAL_DEXSEQ_BINS_P_FDR_AND_COUNTS.tsv", sep="\t")
    case = complete["case"]
    pieces = output/"model_chunks"
    pieces.mkdir(exist_ok=True)
    order = list(models)
    paths = []
    for offset in range(0, len(order), 64):
        check_cancelled()
        tokens = order[offset:offset+64]
        folder = pieces/f"{offset:05d}"
        folder.mkdir(exist_ok=True)
        marker = folder/"COMPLETE.json"
        if marker.exists():
            saved = json.loads(marker.read_text())
            if saved["models"] != tokens:
                raise ValueError("Checkpoint motif identity changed")
            for name, info in saved["files"].items():
                p = folder/name
                if p.stat().st_size != info["bytes"] or sha256(p) != info["SHA256"]:
                    raise ValueError("Corrupt completed motif chunk")
        else:
            if any(folder.iterdir()):
                raise ValueError("Partial chunk retained; inspect before separately resuming that chunk")
            hits, counts = scan_transcript_models(rna, {t: models[t] for t in tokens}, check_cancelled)
            samples, summary = weighted_model_opportunities(usage, counts)
            if len(hits):
                locus_samples, locus_summary = weighted_locus_opportunities(usage, hits, case)
                overlap = bin_hit_overlap(hits, bins)
            else:
                locus_samples, locus_summary, overlap = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
            tables = {"hits": hits, "transcript_scores": counts, "sample_scores": samples, "model_summary": summary,
                      "locus_samples": locus_samples, "locus_summary": locus_summary, "DEXSeq_motif_coordinate_overlap": overlap}
            for name, frame in tables.items():
                frame.to_parquet(folder/(name+".parquet"), index=False)
            write_json(marker, {"models": tokens, "files": {p.name: {"bytes": p.stat().st_size, "SHA256": sha256(p)} for p in folder.glob("*.parquet")}})
        paths.append(folder)
        progress(min(offset+len(tokens), len(order)), len(order))
    exports = {
        "hits": "ALL_REAL_MOTIF_HITS_AND_SPLIT_GENOMIC_POSITIONS.tsv.gz",
        "transcript_scores": "ALL_TRANSCRIPT_MODEL_COUNTS_INCLUDING_ZERO_AND_NA.tsv.gz",
        "sample_scores": "ALL_MODEL_PER_SAMPLE_RNA_AND_TPM_OPPORTUNITIES.tsv.gz",
        "model_summary": "ALL_MODEL_USAGE_VS_EXPRESSION_COMPONENTS.tsv",
        "locus_samples": "ALL_LOCUS_PER_SAMPLE_AVAILABILITY.tsv.gz",
        "locus_summary": "ALL_LOCUS_GAIN_LOSS_RETAINED_AND_COORDINATES.tsv.gz",
        "DEXSeq_motif_coordinate_overlap": "ALL_MOTIF_TO_DEXSEQ_BIN_COORDINATE_OVERLAP.tsv.gz",
    }
    counts_out = {}
    for kind, name in exports.items():
        destination = output/name
        opener = gzip.open if name.endswith(".gz") else open
        n, header = 0, True
        with opener(destination, "wt", encoding="utf-8", newline="") as handle:
            for folder in paths:
                frame = pd.read_parquet(folder/(kind+".parquet"))
                if len(frame.columns):
                    frame.to_csv(handle, sep="\t", index=False, header=header)
                    n += len(frame)
                    header = False
        counts_out[kind] = n
    catalog.to_csv(output/"FULL_MODEL_CATALOG_AND_SOURCE_EVIDENCE.tsv", sep="\t", index=False)
    for name in json.loads((directory/"PREPARED_CASE_MANIFEST.json").read_text())["files"]:
        if name.endswith(".parquet"):
            continue
        shutil.copyfile(directory/name, output/name)
    shutil.copyfile(directory/"PREPARED_CASE_MANIFEST.json", output/"PREPARED_CASE_MANIFEST.json")
    verify_prepared_case(directory)
    result = {"UTC": datetime.now(timezone.utc).isoformat(), "state": "complete_descriptive_usage_sequence_case",
              "case": case, "models": len(models), "transcripts": len(rna), "samples": usage.sample_id.nunique(),
              "rows": counts_out, "original_inputs_unchanged": True,
              "new_motif_P_FDR": "not_computed_single_gene_not_population_enrichment",
              "quantification_rerun": False, "remote_or_other_job_actions": False,
              "source_statistics_preserved_separately": True,
              "PWM_candidates_not_calibrated_binding_sites": True,
              "prepared_manifest_SHA256": binding["prepared_case_manifest_SHA256"],
              "execution_freeze_SHA256": sha256(freeze)}
    write_json(output/"USAGE_MOTIF_COMPLETE.json", result)
    return result
