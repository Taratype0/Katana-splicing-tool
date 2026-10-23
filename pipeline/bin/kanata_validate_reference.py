#!/usr/bin/env python3
"""Stream FASTA/GTF compatibility checks before extracting transcript sequences.

Coordinates are one-based and inclusive, as documented by Ensembl:
https://static.ensembl.org/info/website/upload/gff.html
Build and annotation-release labels are user declarations, not inferred identities.
"""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from typing import TextIO


class ReferenceError(ValueError):
    pass


ATTRIBUTE = re.compile(r'\s*([A-Za-z_][A-Za-z0-9_]*)\s+("(?:[^"\\]|\\.)*"|[^;"\s]+)\s*(?:;|$)')


def open_text(path: Path) -> TextIO:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8-sig")
    return path.open("r", encoding="utf-8-sig")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fasta_lengths(path: Path) -> dict[str, int]:
    """Store one length per contig, without retaining genome sequence in memory."""
    lengths: dict[str, int] = {}
    current: str | None = None
    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None and lengths[current] == 0:
                    raise ReferenceError(f"FASTA contig {current!r} has no sequence")
                header = line[1:].strip()
                if not header:
                    raise ReferenceError(f"FASTA line {line_number}: empty contig identifier")
                current = header.split()[0]
                if current in lengths:
                    raise ReferenceError(f"FASTA line {line_number}: duplicate contig identifier {current!r}")
                lengths[current] = 0
            else:
                if current is None:
                    raise ReferenceError(f"FASTA line {line_number}: sequence before the first header")
                sequence = "".join(line.split())
                if not sequence.isascii() or not sequence.isalpha():
                    raise ReferenceError(f"FASTA line {line_number}: sequence must contain letters only")
                lengths[current] += len(sequence)
    if not lengths:
        raise ReferenceError("FASTA contains no contigs")
    if current is not None and lengths[current] == 0:
        raise ReferenceError(f"FASTA contig {current!r} has no sequence")
    return lengths


def exon_identifiers(attributes: str, line_number: int) -> None:
    """Accept repeated non-ID attributes such as GENCODE's tag annotations."""
    ids: dict[str, str] = {}
    position = 0
    while position < len(attributes):
        if not attributes[position:].strip():
            break
        match = ATTRIBUTE.match(attributes, position)
        if match is None:
            raise ReferenceError(f"GTF line {line_number}: malformed exon attributes")
        key, value = match.group(1, 2)
        if key in {"gene_id", "transcript_id"}:
            if key in ids:
                raise ReferenceError(f"GTF line {line_number}: repeated exon {key}")
            value = value[1:-1] if value.startswith('"') else value
            if not value.strip():
                raise ReferenceError(f"GTF line {line_number}: empty exon {key}")
            ids[key] = value
        position = match.end()
    missing = sorted({"gene_id", "transcript_id"} - set(ids))
    if missing:
        raise ReferenceError(f"GTF line {line_number}: exon missing {', '.join(missing)}")


def annotation_summary(path: Path, contigs: dict[str, int]) -> dict:
    feature_counts: Counter[str] = Counter()
    exon_counts: Counter[str] = Counter()
    annotated_contigs: set[str] = set()
    rows = 0
    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\r\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9 or any(not value for value in fields):
                raise ReferenceError(f"GTF line {line_number}: expected nine nonempty tab-separated fields")
            contig, _, feature, raw_start, raw_end, score, strand, frame, attributes = fields
            if contig not in contigs:
                raise ReferenceError(f"GTF line {line_number}: contig {contig!r} is absent from FASTA; contig names must match exactly")
            try:
                start, end = int(raw_start), int(raw_end)
            except ValueError:
                raise ReferenceError(f"GTF line {line_number}: coordinates must be integers") from None
            if start < 1 or end < start or end > contigs[contig]:
                raise ReferenceError(f"GTF line {line_number}: coordinates {start}-{end} are outside 1-{contigs[contig]} on {contig}")
            if strand not in {"+", "-", "."} or frame not in {"0", "1", "2", "."}:
                raise ReferenceError(f"GTF line {line_number}: invalid strand or frame")
            if score != ".":
                try:
                    valid_score = math.isfinite(float(score))
                except ValueError:
                    valid_score = False
                if not valid_score:
                    raise ReferenceError(f"GTF line {line_number}: score must be a finite number or '.'")
            if feature == "exon":
                exon_identifiers(attributes, line_number)
                exon_counts[contig] += 1
            feature_counts[feature] += 1
            annotated_contigs.add(contig)
            rows += 1
    if not exon_counts:
        raise ReferenceError("GTF contains no exon features")
    return {
        "feature_rows": rows,
        "exon_rows": sum(exon_counts.values()),
        "annotated_contig_count": len(annotated_contigs),
        "annotation_contigs": sorted(annotated_contigs),
        "exons_per_contig": dict(sorted(exon_counts.items())),
        "feature_counts": dict(sorted(feature_counts.items())),
    }


def declaration(value: str | None) -> dict:
    normalized = value.strip() if value else ""
    return {"value": normalized or None, "status": "user_declared" if normalized else "not_declared"}


def validate_reference(fasta: Path, gtf: Path, *, reference_build: str | None = None,
                       annotation_release: str | None = None) -> dict:
    fasta, gtf = fasta.resolve(), gtf.resolve()
    for kind, path in (("FASTA", fasta), ("GTF", gtf)):
        if not path.is_file() or path.stat().st_size == 0:
            raise ReferenceError(f"{kind} file is missing or empty: {path}")
    contigs = fasta_lengths(fasta)
    annotation = annotation_summary(gtf, contigs)
    return {
        "schema_version": "kanata.reference_audit.v1",
        "status": "passed",
        "scope": "fasta_gtf_structural_and_coordinate_compatibility",
        "contig_name_policy": "exact_matching_without_renaming",
        "coordinate_system": "one_based_inclusive",
        "declarations": {
            "reference_build": declaration(reference_build),
            "annotation_release": declaration(annotation_release),
        },
        "declaration_verification": "not_independently_verified",
        "biological_annotation_validity": "not_established_by_structural_checks",
        "fasta": {
            "path": str(fasta), "bytes": fasta.stat().st_size, "sha256": sha256(fasta),
            "contig_count": len(contigs), "total_bases": sum(contigs.values()),
            "contig_lengths": contigs,
            "contigs_without_annotation": sorted(set(contigs) - set(annotation["annotation_contigs"])),
        },
        "gtf": {"path": str(gtf), "bytes": gtf.stat().st_size, "sha256": sha256(gtf), **annotation},
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reference-build")
    parser.add_argument("--annotation-release")
    args = parser.parse_args(argv)
    try:
        audit = validate_reference(args.fasta, args.gtf, reference_build=args.reference_build,
                                   annotation_release=args.annotation_release)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({"status": audit["status"], "contigs": audit["fasta"]["contig_count"], "exons": audit["gtf"]["exon_rows"]}))
        return 0
    except (ReferenceError, OSError, EOFError, UnicodeError) as exc:
        print(f"Reference validation failed: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
