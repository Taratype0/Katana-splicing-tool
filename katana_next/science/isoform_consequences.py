"""Annotation-constrained transcript consequences, not a de novo ORF predictor."""
from __future__ import annotations

from dataclasses import dataclass, field
import gzip
import hashlib
from pathlib import Path
import re

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

ATTR = re.compile(r'(\S+)\s+"([^"]*)"')


def normalized_id(value: str) -> str:
    return re.sub(r"\.\d+$", "", str(value).split("|")[0])


@dataclass
class TranscriptModel:
    transcript_id: str
    gene_id: str
    chrom: str
    strand: str
    biotype: str = ""
    tags: set[str] = field(default_factory=set)
    features: dict[str, list[tuple[int, int, str]]] = field(default_factory=dict)

    @property
    def exons(self):
        return sorted(set((a, b) for a, b, _ in self.features.get("exon", [])), reverse=self.strand == "-")

    def mapped(self, kind):
        """Map genomic half-open intervals to spliced, 5'-to-3' transcript space."""
        result, offset = [], 0
        for lo, hi in self.exons:
            for a, b, phase in self.features.get(kind, []):
                start, end = max(a, lo), min(b, hi)
                if end > start:
                    left = offset + (start - lo if self.strand == "+" else hi - end)
                    result.append((left, left + end - start, phase))
            offset += hi - lo
        return sorted(set(result))


def load_models(gtf: Path, requested: set[str] | None, check_cancelled=lambda: None):
    models = {}
    opener = gzip.open if gtf.suffix == ".gz" else open
    with opener(gtf, "rt", encoding="utf-8") as handle:
        for i, line in enumerate(handle):
            if i % 100000 == 0:
                check_cancelled()
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] not in {"transcript", "exon", "CDS", "start_codon", "stop_codon"}:
                continue
            pairs = ATTR.findall(fields[8])
            attr = dict(pairs)
            raw_id = attr.get("transcript_id", "")
            tx = normalized_id(raw_id)
            if not tx or (requested is not None and tx not in requested):
                continue
            model = models.setdefault(tx, TranscriptModel(raw_id, normalized_id(attr.get("gene_id", "")),
                                      fields[0], fields[6], attr.get("transcript_type", attr.get("transcript_biotype", ""))))
            if model.transcript_id != raw_id or (model.chrom, model.strand) != (fields[0], fields[6]):
                raise ValueError(f"Ambiguous transcript identity in GTF: {tx}")
            model.tags.update(v for k, v in pairs if k == "tag")
            if fields[2] != "transcript":
                model.features.setdefault(fields[2], []).append((int(fields[3]) - 1, int(fields[4]), fields[7]))
    return models


def load_sequences(fasta: Path, requested: set[str], check_cancelled=lambda: None, expected_models=None):
    sequences = {}
    opener = gzip.open if fasta.suffix == ".gz" else open
    with opener(fasta, "rt", encoding="utf-8") as handle:
        for i, record in enumerate(SeqIO.parse(handle, "fasta")):
            if i % 10000 == 0:
                check_cancelled()
            key = normalized_id(record.id)
            if key in requested:
                if key in sequences:
                    raise ValueError(f"Ambiguous transcript FASTA ID: {key}")
                raw_id = record.id.split("|")[0]
                expected = expected_models.get(key) if expected_models else None
                if expected and re.search(r"\.\d+$", raw_id) and re.search(r"\.\d+$", expected.transcript_id) and raw_id != expected.transcript_id:
                    raise ValueError(f"Transcript version mismatch between FASTA ({raw_id}) and GTF ({expected.transcript_id})")
                sequences[key] = str(record.seq).upper()
    return sequences


def transcript_features(key, model, sequence, *, nmd_distance=55):
    result = dict(transcript_id=key, sequence_status="missing_annotation", protein_status="not_evaluable",
                  nmd_prediction="not_evaluable", protein_sequence="", protein_sha256="", protein_length=None,
                  transcript_length=None, cds_length=None, utr5_length=None, utr3_length=None,
                  stop_to_last_junction_nt=None, exon_count=None, biotype="", gene_id="", exon_intervals="")
    if model is None:
        return result
    exons = model.exons
    length = sum(b - a for a, b in exons)
    result.update(gene_id=model.gene_id, biotype=model.biotype, exon_count=len(exons), transcript_length=length,
                  exon_intervals=";".join(f"{model.chrom}:{a}-{b}:{model.strand}" for a, b in exons))
    if not exons or model.strand not in {"+", "-"} or any(a < 0 or b <= a for a, b in exons):
        result["sequence_status"] = "invalid_exon_annotation"
        return result
    genomic = sorted(exons)
    if any(a[1] > b[0] for a, b in zip(genomic, genomic[1:])):
        result["sequence_status"] = "overlapping_exons"
        return result
    result["sequence_status"] = "missing_sequence" if sequence is None else "length_mismatch" if len(sequence) != length else "length_matched"
    cds, stops, starts = model.mapped("CDS"), model.mapped("stop_codon"), model.mapped("start_codon")
    result["cds_length"] = sum(b - a for a, b, _ in cds)
    if not cds:
        result["protein_status"] = "no_annotated_CDS"
        return result
    # A CDS is contiguous in spliced transcript space. A missing GTF segment is not repaired.
    contiguous = all(a[1] == b[0] for a, b in zip(cds, cds[1:]))
    if not contiguous or result["cds_length"] != sum(b-a for a, b, _ in model.features.get("CDS", [])):
        result["protein_status"] = "invalid_CDS_mapping"
        return result
    cumulative = 0
    for a, b, phase in cds:
        if phase != str((3-cumulative % 3) % 3):
            result["protein_status"] = "inconsistent_or_partial_CDS_phase"
            return result
        cumulative += b-a
    if sequence is None or result["sequence_status"] != "length_matched":
        return result
    partial = bool(model.tags & {"cds_start_NF", "cds_end_NF", "mRNA_start_NF", "mRNA_end_NF"})
    start, cds_end = cds[0][0], cds[-1][1]
    stop_positions = sorted({p for a, b, _ in stops for p in range(a, b)})
    start_positions = sorted({p for a, b, _ in starts for p in range(a, b)})
    if partial or cds[0][2] != "0" or start_positions != list(range(start, start+3)) or len(stop_positions) != 3:
        result["protein_status"] = "partial_or_missing_codon_annotation"
        return result
    stop = stop_positions[0]
    if stop_positions != list(range(stop, stop+3)) or stop not in {cds_end, cds_end-3}:
        result["protein_status"] = "inconsistent_stop_annotation"
        return result
    coding = sequence[start:stop+3]
    result.update(utr5_length=start, utr3_length=length-stop-3)
    if len(coding) % 3 or coding[:3] != "ATG" or coding[-3:] not in {"TAA", "TAG", "TGA"} or set(coding) - set("ACGT"):
        result["protein_status"] = "invalid_or_ambiguous_CDS_sequence"
        return result
    protein = str(Seq(coding).translate())
    if "*" in protein[:-1]:
        result["protein_status"] = "internal_stop_requires_review"
        return result
    protein = protein[:-1]
    result.update(protein_status="validated_annotated_ORF", protein_sequence=protein, protein_length=len(protein),
                  protein_sha256=hashlib.sha256(protein.encode()).hexdigest())
    distance = length - (exons[-1][1]-exons[-1][0]) - (stop+3) if len(exons) > 1 else None
    result["stop_to_last_junction_nt"] = distance
    result["nmd_prediction"] = "rule_positive" if distance is not None and distance > nmd_distance else "rule_negative_not_proof_of_escape"
    return result


def consequences(switches, models, sequences, *, nmd_distance=55):
    ids = set(switches.numerator_dominant_transcript.map(normalized_id)) | set(switches.denominator_dominant_transcript.map(normalized_id))
    features = {tx: transcript_features(tx, models.get(tx), sequences.get(tx), nmd_distance=nmd_distance) for tx in sorted(ids)}
    rows = []
    for row in switches.to_dict("records"):
        a = features[normalized_id(row["numerator_dominant_transcript"])]
        b = features[normalized_id(row["denominator_dominant_transcript"])]
        valid = a["protein_status"] == b["protein_status"] == "validated_annotated_ORF"
        row.update(protein_change=("same_sequence" if a["protein_sha256"] == b["protein_sha256"] else "different_sequence") if valid else "not_evaluable",
                   protein_length_delta=a["protein_length"]-b["protein_length"] if valid else None,
                   domain_consequence="not_assessed_no_protein_domain_mapping",
                   consequence_scope="annotation-based support; dominant-TPM switch is descriptive, not a switch-level statistical test")
        for side, feature in (("numerator", a), ("denominator", b)):
            for key in ("protein_status", "nmd_prediction", "sequence_status", "biotype", "exon_count", "cds_length", "utr5_length", "utr3_length"):
                row[f"{side}_{key}"] = feature[key]
        row["exon_structure_change"] = a["exon_intervals"] != b["exon_intervals"] if a["exon_intervals"] and b["exon_intervals"] else None
        row["annotation_gene_match"] = a["gene_id"] == b["gene_id"] == normalized_id(row["gene_id"])
        if not row["annotation_gene_match"]:
            row.update(protein_change="not_evaluable_gene_mismatch", protein_length_delta=None)
        rows.append(row)
    return pd.DataFrame(rows), pd.DataFrame(features.values())
