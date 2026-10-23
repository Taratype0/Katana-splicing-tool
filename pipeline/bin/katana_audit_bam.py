#!/usr/bin/env python3
"""Gate genome-BAM eligibility using exact primary-record metrics, not averages.

The compact input is produced by katana_bam_primary_metrics.awk from a complete
samtools view stream. Direct BAM mode invokes samtools without modifying BAM or
reference files. This establishes technical compatibility, not biological validity.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Iterable


class AuditError(ValueError):
    pass


PAIR_PARTITIONS = tuple(
    f"{mapping}_mate_{mate}_{side}_records"
    for mapping in ("mapped", "unmapped")
    for mate in ("mapped", "unmapped")
    for side in ("read1", "read2")
)
REQUIRED_METRICS = (
    "primary_records", "paired_records", "unpaired_records", "read1_records",
    "read2_records", "invalid_pair_flag_records", "mapped_records",
    "unmapped_records", "missing_seq_records", "invalid_sequence_records",
    "malformed_cigar_records", "cigar_query_length_mismatch_records",
    "mapped_without_cigar_records", "soft_clipped_records", "hard_clipped_records",
    "duplicate_records", "qcfail_records", "spliced_records",
) + PAIR_PARTITIONS
CIGAR = re.compile(r"([1-9][0-9]*)([MIDNSHP=X])")
BASES = re.compile(r"^[ACMGRSVTWYHKDBNacmgrsvtwyhkdbn]+$")


def open_text(path: str | Path):
    path = Path(path)
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def parse_header(lines: Iterable[str]) -> dict:
    sequences, programs, header = {}, [], {}
    for line in lines:
        fields = line.rstrip("\r\n").split("\t")
        if not fields or not fields[0].startswith("@"):
            raise AuditError("BAM header contains a non-header record")
        tags = {}
        for field in fields[1:]:
            if ":" in field:
                key, value = field.split(":", 1)
                tags[key] = value
        if fields[0] == "@SQ":
            name = tags.get("SN", "")
            try:
                length = int(tags.get("LN", ""))
            except ValueError as exc:
                raise AuditError("@SQ is missing a valid LN") from exc
            if not name or name in sequences or length <= 0:
                raise AuditError("@SQ names must be unique and lengths positive")
            sequences[name] = dict(length=length, **{key: tags[key] for key in ("M5", "UR", "AS", "SP") if key in tags})
        elif fields[0] == "@PG":
            programs.append(tags)
        elif fields[0] == "@HD":
            header = tags
    if not sequences:
        raise AuditError("BAM has no @SQ sequence dictionary")
    return {"sequences": sequences, "programs": programs, "header": header}


def parse_fai(lines: Iterable[str]) -> dict[str, int]:
    lengths = {}
    for line in lines:
        fields = line.rstrip("\r\n").split("\t")
        try:
            length = int(fields[1])
        except (IndexError, ValueError) as exc:
            raise AuditError("Reference FAI contains an invalid length record") from exc
        if not fields[0] or fields[0] in lengths or length <= 0:
            raise AuditError("Reference FAI names must be unique and lengths positive")
        lengths[fields[0]] = length
    if not lengths:
        raise AuditError("Reference FAI is empty")
    return lengths


def fasta_lengths(path: str | Path) -> dict[str, int]:
    """Read the supplied FASTA without creating a .fai beside an original file."""
    lengths, name, length = {}, None, 0
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    if length == 0:
                        raise AuditError(f"Empty reference sequence: {name}")
                    lengths[name] = length
                words = line[1:].split()
                if not words or words[0] in lengths:
                    raise AuditError("FASTA sequence names must be nonempty and unique")
                name, length = words[0], 0
            elif line.strip():
                if name is None:
                    raise AuditError("Reference FASTA sequence appears before its header")
                length += len(line.strip())
    if name is None or length == 0:
        raise AuditError("Reference FASTA is empty or has an empty sequence")
    lengths[name] = length
    return lengths


def primary_metrics(lines: Iterable[str]) -> dict:
    metrics = Counter({name: 0 for name in REQUIRED_METRICS})
    lengths, references = Counter(), Counter()
    for number, line in enumerate(lines, 1):
        if not line.strip() or line.startswith("@"):
            continue
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) < 11:
            raise AuditError(f"Malformed SAM record at stream line {number}")
        try:
            flag = int(fields[1])
        except ValueError as exc:
            raise AuditError(f"Invalid SAM FLAG at stream line {number}") from exc
        if flag < 0 or flag > 65535:
            raise AuditError(f"Invalid SAM FLAG at stream line {number}")
        if flag & (256 | 2048):
            continue
        metrics["primary_records"] += 1
        paired, first, second, unmapped = bool(flag & 1), bool(flag & 64), bool(flag & 128), bool(flag & 4)
        metrics["paired_records" if paired else "unpaired_records"] += 1
        metrics["read1_records"] += first
        metrics["read2_records"] += second
        mate_unmapped = bool(flag & 8)
        metrics["invalid_pair_flag_records"] += (
            (paired and first == second)
            or (not paired and bool(flag & (2 | 8 | 32 | 64 | 128)))
            or (bool(flag & 2) and (unmapped or mate_unmapped))
        )
        if paired and first != second:
            mapping, mate, side = (
                "unmapped" if unmapped else "mapped",
                "unmapped" if mate_unmapped else "mapped",
                "read1" if first else "read2",
            )
            metrics[f"{mapping}_mate_{mate}_{side}_records"] += 1
        metrics["unmapped_records" if unmapped else "mapped_records"] += 1
        metrics["duplicate_records"] += bool(flag & 1024)
        metrics["qcfail_records"] += bool(flag & 512)
        if not unmapped:
            references[fields[2]] += 1
        sequence, cigar = fields[9], fields[5]
        if sequence == "*":
            metrics["missing_seq_records"] += 1
        else:
            lengths[len(sequence)] += 1
            metrics["invalid_sequence_records"] += not bool(BASES.fullmatch(sequence))
        if cigar == "*":
            metrics["mapped_without_cigar_records"] += not unmapped
            continue
        tokens = CIGAR.findall(cigar)
        invalid_clip_position = any(
            op == "H" and index not in (0, len(tokens) - 1)
            or op == "S" and not (
                all(before == "H" for _, before in tokens[:index])
                or all(after == "H" for _, after in tokens[index + 1:])
            )
            for index, (_, op) in enumerate(tokens)
        )
        if not tokens or "".join(n + op for n, op in tokens) != cigar or invalid_clip_position:
            metrics["malformed_cigar_records"] += 1
            continue
        operations = {op for _, op in tokens}
        metrics["soft_clipped_records"] += "S" in operations
        metrics["hard_clipped_records"] += "H" in operations
        metrics["spliced_records"] += "N" in operations
        query_length = sum(int(n) for n, op in tokens if op in "MIS=X")
        metrics["cigar_query_length_mismatch_records"] += sequence != "*" and query_length != len(sequence)
    return {"metrics": dict(metrics), "read_lengths": dict(lengths), "mapped_references": dict(references)}


def parse_metrics(lines: Iterable[str]) -> dict:
    result = {"metrics": {}, "read_lengths": {}, "mapped_references": {}}
    reader = csv.DictReader(lines, delimiter="\t")
    if reader.fieldnames != ["kind", "name", "count"]:
        raise AuditError("Primary metrics must have columns kind, name, count")
    kinds = {"metric": "metrics", "read_length": "read_lengths", "mapped_reference": "mapped_references"}
    for row in reader:
        if row["kind"] not in kinds:
            raise AuditError(f"Unknown primary metrics kind: {row['kind']}")
        try:
            count = int(row["count"])
            name = int(row["name"]) if row["kind"] == "read_length" else row["name"]
        except (ValueError, TypeError) as exc:
            raise AuditError("Primary metrics contain a non-integer count or read length") from exc
        target = result[kinds[row["kind"]]]
        if count < 0 or name in target or not str(name) or (row["kind"] == "read_length" and name <= 0):
            raise AuditError("Primary metrics contain negative counts, invalid names, or duplicate rows")
        target[name] = count
    missing = set(REQUIRED_METRICS) - result["metrics"].keys()
    if missing:
        raise AuditError(f"Missing exact primary metrics: {', '.join(sorted(missing))}")
    return result


def audit(header: dict, reference: dict[str, int], observed: dict, *, sample: str,
          single_end: bool, read_length: int, allow_variable_read_length: bool = False,
          allow_reference_subset: bool = False, reference_dictionary: dict | None = None) -> dict:
    if not sample or read_length <= 0:
        raise AuditError("Sample ID and positive nominal read length are required")
    sequences = header["sequences"]
    extra = sorted(set(sequences) - set(reference))
    missing = sorted(set(reference) - set(sequences))
    if extra:
        raise AuditError(f"BAM @SQ contigs absent from reference: {', '.join(extra[:10])}")
    if missing and not allow_reference_subset:
        raise AuditError(f"Reference contigs missing from BAM @SQ: {', '.join(missing[:10])}; subset policy not enabled")
    mismatch = [name for name in sequences if sequences[name]["length"] != reference[name]]
    if mismatch:
        raise AuditError(f"BAM @SQ length differs from reference for: {', '.join(mismatch[:10])}")
    md5_checked = []
    if reference_dictionary is not None:
        ref_sequences = reference_dictionary["sequences"]
        if set(ref_sequences) != set(reference) or any(ref_sequences[name]["length"] != reference[name] for name in reference):
            raise AuditError("Reference dictionary and FAI/FASTA names or lengths disagree")
        for name, sequence in sequences.items():
            if "M5" not in sequence:
                continue
            observed_md5, expected_md5 = sequence["M5"], ref_sequences[name].get("M5", "")
            if not re.fullmatch(r"[0-9a-fA-F]{32}", observed_md5) or not re.fullmatch(r"[0-9a-fA-F]{32}", expected_md5):
                raise AuditError(f"Missing or invalid reference/BAM M5 for {name}")
            if observed_md5.lower() != expected_md5.lower():
                raise AuditError(f"BAM reference sequence M5 differs from supplied reference for {name}")
            md5_checked.append(name)
    metrics = observed["metrics"]
    for name in REQUIRED_METRICS:
        if name not in metrics or not isinstance(metrics[name], int) or metrics[name] < 0:
            raise AuditError(f"Invalid or missing exact primary metric: {name}")
    total = metrics["primary_records"]
    if total <= 0 or metrics["mapped_records"] <= 0:
        raise AuditError("BAM must contain primary mapped records")
    if any(metrics[name] > total for name in REQUIRED_METRICS if name != "primary_records"):
        raise AuditError("Primary metric exceeds total primary records")
    if metrics["paired_records"] + metrics["unpaired_records"] != total or metrics["mapped_records"] + metrics["unmapped_records"] != total:
        raise AuditError("Primary layout/mapping totals are inconsistent")
    if sum(observed["read_lengths"].values()) + metrics["missing_seq_records"] != total:
        raise AuditError("Read-length histogram does not cover every primary record")
    if sum(observed["mapped_references"].values()) != metrics["mapped_records"]:
        raise AuditError("Mapped-reference totals do not cover every mapped primary record")
    unknown = set(observed["mapped_references"]) - set(sequences)
    if unknown:
        raise AuditError(f"Mapped primary records refer to contigs absent from @SQ: {', '.join(sorted(unknown))}")
    fatal = {
        "missing_seq_records": "missing primary read sequence",
        "invalid_sequence_records": "invalid/unknown primary read sequence encoding",
        "invalid_pair_flag_records": "contradictory paired/read1/read2 flags",
        "malformed_cigar_records": "malformed primary CIGAR",
        "cigar_query_length_mismatch_records": "CIGAR query length differs from stored sequence",
        "mapped_without_cigar_records": "mapped primary records lack CIGAR",
        "hard_clipped_records": "hard-clipped primary sequence cannot establish complete read length",
    }
    for name, message in fatal.items():
        if metrics[name]:
            raise AuditError(f"BAM contains {metrics[name]} records with {message}")
    if metrics["read1_records"] + metrics["read2_records"] != metrics["paired_records"]:
        raise AuditError("Exact primary read1/read2 counts are inconsistent with paired records")
    if single_end and metrics["paired_records"]:
        raise AuditError("Declared single-end BAM contains paired primary records")
    if not single_end:
        if metrics["unpaired_records"] or not metrics["read1_records"] or not metrics["read2_records"]:
            raise AuditError("Declared paired-end BAM must contain paired primary records and both read1/read2")
        for side in ("read1", "read2"):
            if sum(metrics[name] for name in PAIR_PARTITIONS if f"_{side}_" in name) != metrics[f"{side}_records"]:
                raise AuditError("Primary mate-mapping partitions do not cover every read1/read2 record")
        # A mapped-only genome BAM can legally omit an unmapped mate.  Counts
        # of all R1/R2 records therefore need not match.  Records claiming that
        # BOTH ends map must still be balanced; this is an aggregate check,
        # not proof that every QNAME has its correctly identified mate.
        if metrics["mapped_mate_mapped_read1_records"] != metrics["mapped_mate_mapped_read2_records"]:
            raise AuditError("Primary records claiming both mates mapped have unequal read1/read2 counts")
    lengths = sorted(length for length, count in observed["read_lengths"].items() if count)
    if not lengths or lengths[-1] > read_length:
        raise AuditError("Observed primary sequence length exceeds declared nominal read length or is unknown")
    if not allow_variable_read_length and lengths != [read_length]:
        raise AuditError("Observed sequence lengths differ from declared read length; explicit variable-read-length policy is required")
    return {
        "schema_version": "katana.bam_audit.v1", "status": "passed", "sample": sample,
        "single_end": single_end, "declared_nominal_read_length": read_length,
        "read_length_policy": "shorter_or_equal" if allow_variable_read_length else "exact",
        "allow_reference_subset": allow_reference_subset, "reference_contigs_missing_from_header": missing,
        "header_contigs": sequences, "header": header["header"], "program_records": header["programs"],
        "primary_metrics": metrics, "observed_read_length_histogram": {str(length): observed["read_lengths"][length] for length in lengths},
        "observed_min_read_length": lengths[0], "observed_max_read_length": lengths[-1],
        "variable_observed_read_lengths": len(lengths) > 1,
        "mapped_reference_counts": observed["mapped_references"],
        "reference_identity": "available_BAM_M5_verified" if md5_checked else "contig_names_and_lengths_verified_sequence_identity_not_established",
        "reference_m5_checked_contigs": md5_checked,
        "reference_m5_unavailable_contigs": sorted(set(sequences) - set(md5_checked)),
        "strandedness_validation": "not_established", "unmapped_read_retention": "not_established",
        "complete_fragment_identity": "not_established_by_aggregate_audit",
        "eligibility_policy": "paired_mapping_partition_v2",
        "primary_read1_minus_read2": metrics["read1_records"] - metrics["read2_records"],
        "mate_flag_observations": {
            "mapped_primary_reads_flagging_mate_unmapped": (
                metrics["mapped_mate_unmapped_read1_records"] + metrics["mapped_mate_unmapped_read2_records"]
            ),
            "interpretation": "Flag observations are not an exact count of absent mates or complete fragments. Input reads are preserved; native downstream tools apply their own eligibility rules.",
        },
        "biological_validity": "not_established",
        "alignment_history": "Header PG records retained; removed reads, alignments, and filtering history cannot be recovered or inferred from their absence.",
        "audit_scope": "All primary records, including unmapped records; secondary and supplementary records excluded from layout/read-length gates.",
    }


def parse_bool(value) -> bool:
    if isinstance(value, bool):
        return value
    if str(value).lower() not in ("true", "false"):
        raise AuditError("single_end must be explicitly true or false")
    return str(value).lower() == "true"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam")
    parser.add_argument("--fasta")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--header")
    parser.add_argument("--primary-records")
    parser.add_argument("--primary-metrics")
    parser.add_argument("--fai")
    parser.add_argument("--reference-dict", help="samtools dict SAM header generated from the supplied reference; checks any available BAM M5")
    parser.add_argument("--meta-json")
    parser.add_argument("--sample")
    parser.add_argument("--single-end", choices=("true", "false"))
    parser.add_argument("--read-length", type=int)
    parser.add_argument("--allow-variable-read-length", action="store_true")
    parser.add_argument("--allow-reference-subset", action="store_true")
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    result = {"schema_version": "katana.bam_audit.v1", "status": "failed"}
    try:
        meta = json.loads(Path(args.meta_json).read_text(encoding="utf-8")) if args.meta_json else {}
        sample = args.sample or meta.get("sample") or meta.get("id")
        single_end = parse_bool(args.single_end if args.single_end is not None else meta.get("single_end"))
        read_length = args.read_length if args.read_length is not None else int(meta.get("read_length", 0))
        result["sample"] = sample
        if args.bam:
            if not args.fasta or any((args.header, args.primary_records, args.primary_metrics, args.fai, args.reference_dict)):
                raise AuditError("Direct BAM mode requires --fasta and cannot combine captured-file options")
            version = subprocess.run([args.samtools, "--version"], check=True, capture_output=True, text=True).stdout.splitlines()[0]
            header_text = subprocess.run([args.samtools, "view", "-H", args.bam], check=True, capture_output=True, text=True).stdout
            header = parse_header(header_text.splitlines())
            reference = fasta_lengths(args.fasta)
            reference_dictionary = parse_header(subprocess.run([args.samtools, "dict", args.fasta], check=True, capture_output=True, text=True).stdout.splitlines())
            # Inherit stderr so large diagnostic logs cannot deadlock an unread pipe.
            with subprocess.Popen([args.samtools, "view", "-F", "2304", args.bam], stdout=subprocess.PIPE, text=True) as process:
                try:
                    observed = primary_metrics(process.stdout)
                except Exception:
                    process.terminate()
                    raise
                if process.wait() != 0:
                    raise AuditError("samtools view failed before completing the primary-record audit")
        else:
            if not args.header or not args.fai or bool(args.primary_records) == bool(args.primary_metrics) or args.fasta:
                raise AuditError("Captured mode requires --header --fai and exactly one of --primary-records/--primary-metrics")
            with open_text(args.header) as handle:
                header = parse_header(handle)
            with open_text(args.fai) as handle:
                reference = parse_fai(handle)
            with open_text(args.primary_records or args.primary_metrics) as handle:
                observed = primary_metrics(handle) if args.primary_records else parse_metrics(handle)
            reference_dictionary = None
            if args.reference_dict:
                with open_text(args.reference_dict) as handle:
                    reference_dictionary = parse_header(handle)
            version = "captured_by_upstream_samtools_process_see_versions_yml"
        result = audit(header, reference, observed, sample=sample, single_end=single_end,
                       read_length=read_length, allow_variable_read_length=args.allow_variable_read_length,
                       allow_reference_subset=args.allow_reference_subset, reference_dictionary=reference_dictionary)
        result["samtools_version"] = version
    except (AuditError, ValueError, OSError, subprocess.SubprocessError) as exc:
        result["error"] = str(exc)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if result["status"] != "passed":
        print(f"BAM eligibility audit failed: {result.get('error', 'unknown error')}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
