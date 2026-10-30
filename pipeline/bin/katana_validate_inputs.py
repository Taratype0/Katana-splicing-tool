#!/usr/bin/env python3
"""Validate the public input contract without importing Katana analysis code."""
from __future__ import annotations

import argparse
import base64
import csv
import hashlib
import json
import re
import sys
from pathlib import Path

SAMPLE_COLUMNS = ["sample", "condition", "source", "assay", "fastq_1", "fastq_2", "genome_bam", "transcriptome_bam", "strandedness", "single_end", "read_length", "subject", "batch", "fastq_lanes"]
CONTRAST_COLUMNS = ["contrast", "treatment", "control"]
SUPPORTED_CONTRAST_COLUMNS = frozenset(CONTRAST_COLUMNS + ["paired", "blocking", "design", "formula", "description"])
SAFE_ID = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]*$")
SAFE_CONDITION = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")
REMOTE_URI = re.compile(r"^[A-Za-z][A-Za-z0-9+.-]*://")


class InputError(ValueError):
    pass


def read_sheet(path: Path, required: set[str]) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t" if path.suffix.lower() == ".tsv" else ",")
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)):
            raise InputError(f"Duplicate column names in {path}")
        missing = sorted(required - set(fields))
        if missing:
            raise InputError(f"{path.name} is missing columns: {', '.join(missing)}")
        rows = []
        for line, raw in enumerate(reader, start=2):
            if None in raw:
                raise InputError(f"{path.name}:{line}: too many columns")
            row = {key: (value or "").strip() for key, value in raw.items()}
            if any(row.values()):
                rows.append(row)
        if not rows:
            raise InputError(f"No records in {path}")
        return fields, rows


def boolean(value: str, *, context: str) -> bool:
    if value.lower() in {"true", "1", "yes"}:
        return True
    if value.lower() in {"false", "0", "no"}:
        return False
    raise InputError(f"{context}: expected true or false, got {value!r}")


def input_path(value: str, base: Path, *, check: bool, context: str) -> str:
    if not value:
        return ""
    if REMOTE_URI.match(value):
        raise InputError(f"{context}: v1 requires a local or mounted file path; download {value} first")
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = base / path
    path = path.resolve()
    if check and (not path.is_file() or path.stat().st_size == 0):
        raise InputError(f"{context}: missing or empty input file: {path}")
    return str(path)


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def validate_samples(rows: list[dict[str, str]], base: Path, *, check: bool, read_length: int | None) -> list[dict[str, str]]:
    seen: set[str] = set()
    observed_paths: dict[str, str] = {}
    result = []
    for source_row in rows:
        row = dict(source_row)
        sample = row["sample"]
        if row.get("transcriptome_bam"):
            raise InputError(f"{sample}: v1 quantifies a genome BAM using Salmon genome projection; transcriptome_bam is not a supported input")
        if not SAFE_ID.fullmatch(sample) or (sample in seen and not row.get("lane")):
            raise InputError(f"Invalid or duplicate sample ID {sample!r}; use unique filesystem-safe IDs")
        if not SAFE_CONDITION.fullmatch(row["condition"]):
            raise InputError(f"{sample}: condition must start with a letter and contain letters, digits or underscores")
        assay = (row.get("assay") or "").strip().lower() or "bulk_rna_seq"
        if assay != "bulk_rna_seq":
            raise InputError(f"{sample}: unsupported assay {assay!r}; the complete route currently accepts only bulk_rna_seq. A declared assay is not inferred from the input files.")
        row["assay"] = assay
        seen.add(sample)
        fq1, fq2, bam = (row.get(field, "") for field in ("fastq_1", "fastq_2", "genome_bam"))
        if bool(fq1) == bool(bam):
            raise InputError(f"{sample}: provide either fastq_1 (and optional fastq_2) or genome_bam")
        if fq2 and not fq1:
            raise InputError(f"{sample}: fastq_2 requires fastq_1")
        detected_source = "fastq" if fq1 else "genome_bam"
        if row.get("source", detected_source) not in {"", "auto", detected_source}:
            raise InputError(f"{sample}: source disagrees with supplied files")
        row["source"] = detected_source
        supplied_single = row.get("single_end", "")
        if fq1:
            single = not bool(fq2)
            if supplied_single and boolean(supplied_single, context=sample) != single:
                raise InputError(f"{sample}: single_end conflicts with FASTQ pairing")
        else:
            if not supplied_single:
                raise InputError(f"{sample}: BAM input requires explicit single_end=true/false")
            single = boolean(supplied_single, context=sample)
            if not bam.lower().endswith(".bam"):
                raise InputError(f"{sample}: genome_bam must be a BAM file")
        row["single_end"] = str(single).lower()
        strand = row.get("strandedness", "").lower()
        if strand not in {"forward", "reverse", "unstranded", "auto"}:
            raise InputError(f"{sample}: strandedness must be forward, reverse, unstranded or auto")
        row["strandedness"] = strand
        nominal_length = row.get("read_length") or read_length
        try:
            length = int(nominal_length)  # type: ignore[arg-type]
        except (ValueError, TypeError):
            raise InputError(f"{sample}: provide a positive read_length or --read-length") from None
        if length < 1:
            raise InputError(f"{sample}: read_length must be positive")
        row["read_length"] = str(length)
        for field in ("fastq_1", "fastq_2", "genome_bam", "transcriptome_bam"):
            row[field] = input_path(row.get(field, ""), base, check=check, context=f"{sample}/{field}")
            if row[field]:
                owner = observed_paths.get(row[field])
                if owner is not None:
                    raise InputError(f"{sample}: input {row[field]} is also assigned to {owner}; do not count one library twice")
                observed_paths[row[field]] = sample
        for field in ("subject", "batch"):
            row.setdefault(field, "")
        # This field is generated here; users supply one CSV row per lane.
        row["fastq_lanes"] = ""
        result.append(row)
    if len({row["read_length"] for row in result}) > 1:
        raise InputError("Complete mode requires a single declared read_length across the study, including samples in different contrasts")
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in result:
        grouped.setdefault(row["sample"], []).append(row)
    collapsed = []
    for sample, lanes in grouped.items():
        if len(lanes) == 1:
            collapsed.append(lanes[0])
            continue
        lane_ids = [row.get("lane", "") for row in lanes]
        if not all(lane_ids) or len(set(lane_ids)) != len(lane_ids):
            raise InputError(f"{sample}: repeated sample rows require distinct nonempty lane IDs")
        if any(row["source"] != "fastq" for row in lanes):
            raise InputError(f"{sample}: lane merging is supported for FASTQ inputs only")
        ignored = {"fastq_1", "fastq_2", "lane", "fastq_lanes"}
        metadata = [{k: v for k, v in row.items() if k not in ignored} for row in lanes]
        if any(row != metadata[0] for row in metadata[1:]):
            raise InputError(f"{sample}: lane rows must have identical biological metadata, pairing and library settings")
        lanes = sorted(lanes, key=lambda row: row["lane"])
        row = dict(lanes[0])
        row["fastq_lanes"] = base64.b64encode(json.dumps([
            [lane["fastq_1"], lane["fastq_2"]] if lane["fastq_2"] else [lane["fastq_1"]]
            for lane in lanes
        ]).encode()).decode('ascii')
        row["lane"] = ";".join(lane["lane"] for lane in lanes)
        collapsed.append(row)
    return collapsed


def validate_contrasts(rows: list[dict[str, str]], samples: list[dict[str, str]], *, design: str, analysis_purpose: str) -> tuple[list[dict[str, str]], list[dict]]:
    unknown_columns = sorted({column for row in rows for column in row} - SUPPORTED_CONTRAST_COLUMNS)
    if unknown_columns:
        raise InputError(f"Unsupported contrastsheet columns: {', '.join(unknown_columns)}. Supported columns are contrast,treatment,control,paired,blocking,design,formula,description; unsupported models, aliases and sample selections are not silently ignored.")
    supported = {"~condition": False, "~subject+condition": True}
    if re.sub(r"\s+", "", design) not in supported:
        raise InputError("Complete entry supports '~ condition' or '~ subject + condition'; other models need method-specific support")
    groups: dict[str, list[dict[str, str]]] = {}
    for row in samples:
        groups.setdefault(row["condition"], []).append(row)
    seen: set[str] = set()
    seen_pairs: dict[tuple[str, str], str] = {}
    audit = []
    validated = []
    for source_row in rows:
        row = dict(source_row)
        name, treatment, control = (row[field] for field in CONTRAST_COLUMNS)
        if not SAFE_ID.fullmatch(name) or name in seen:
            raise InputError(f"Invalid or duplicate contrast ID {name!r}")
        seen.add(name)
        if treatment == control or treatment not in groups or control not in groups:
            raise InputError(f"{name}: treatment/control must name two different existing conditions")
        pair = (treatment, control)
        if pair in seen_pairs:
            raise InputError(f"{name}: repeated treatment/control pair {treatment}/{control} already declared as {seen_pairs[pair]}; use one contrast ID per directional comparison")
        seen_pairs[pair] = name
        declared = [re.sub(r"\s+", "", row[field]) for field in ("design", "formula") if row.get(field)]
        if len(set(declared)) > 1 or any(value not in supported for value in declared):
            raise InputError(f"{name}: conflicting or unsupported design/formula")
        flag = row.get("paired", "").strip().lower()
        if flag not in {"", "true", "1", "yes", "false", "0", "no"}:
            raise InputError(f"{name}: paired must be true or false")
        model = declared[0] if declared else ("~subject+condition" if flag in {"true", "1", "yes"} else re.sub(r"\s+", "", design))
        paired = supported[model]
        if flag and paired != (flag in {"true", "1", "yes"}):
            raise InputError(f"{name}: paired flag conflicts with design")
        if row.get("blocking") and (not paired or row["blocking"] != "subject"):
            raise InputError(f"{name}: the supported paired blocking column is subject")
        row.update(paired="true" if paired else "false", design="~ subject + condition" if paired else "~ condition")
        selection = groups[treatment] + groups[control]
        selected_subjects = [member["subject"] for member in selection]
        if any(selected_subjects) and not all(selected_subjects):
            raise InputError(f"{name}: subject identity is only partially specified in this contrast; provide every selected subject or leave all selected subject values empty")
        for condition in (treatment, control):
            members = groups[condition]
            if len(members) < 2:
                raise InputError(f"{name}: condition {condition} has fewer than two samples; differential testing is unavailable")
            subjects = [member["subject"] for member in members if member["subject"]]
            if analysis_purpose == "analysis" and subjects and len(subjects) != len(set(subjects)):
                raise InputError(f"{name}: repeated subjects within {condition}; technical partitions are not biological replicates")
        for field in ("single_end", "strandedness", "read_length"):
            # auto is an unresolved declaration, not a fourth library type.
            # Known conflicting values still fail here; auto samples are checked
            # again against every selected arm after native strand inference.
            values = {member[field] for member in selection if field != "strandedness" or member[field] != "auto"}
            if len(values) > 1:
                raise InputError(f"{name}: mixed {field} is unsupported in the shared rMATS complete-bundle route")
        left_subjects = {member["subject"] for member in groups[treatment] if member["subject"]}
        right_subjects = {member["subject"] for member in groups[control] if member["subject"]}
        if paired:
            if not all(selected_subjects) or left_subjects != right_subjects:
                raise InputError(f"{name}: paired analysis requires complete treatment/control subject pairs")
            if len(left_subjects) != len(groups[treatment]) or len(right_subjects) != len(groups[control]):
                raise InputError(f"{name}: paired analysis requires exactly one sample per subject per condition")
        elif analysis_purpose == "analysis" and left_subjects & right_subjects:
            raise InputError(f"{name}: shared subjects across conditions require an explicit paired design")
        ordered = {condition: sorted(groups[condition], key=lambda member: member["subject"]) if paired else groups[condition]
                   for condition in (treatment, control)}
        validated.append(row)
        audit.append({"contrast": name, "treatment": treatment, "control": control,
                      "rmats_group1_samples": [member["sample"] for member in ordered[treatment]],
                      "rmats_group2_samples": [member["sample"] for member in ordered[control]],
                      "effect_direction": "treatment_minus_control", "design": row["design"], "paired": paired,
                      "pair_ids": sorted(left_subjects) if paired else [],
                      "covariates_adjusted": ["subject"] if paired else [], "analysis_purpose": analysis_purpose})
    return validated, audit


def write_sheet(path: Path, priority: list[str], original: list[str], rows: list[dict[str, str]]) -> None:
    fields = priority + [field for field in original if field not in priority]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samplesheet", type=Path, required=True)
    parser.add_argument("--samplesheet-base", type=Path, help="Original samplesheet directory for resolving relative inputs when the sheet was copied or staged.")
    parser.add_argument("--contrasts", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--read-length", type=int)
    parser.add_argument("--design", default="~ condition")
    parser.add_argument("--analysis-purpose", choices=["analysis", "engine_smoke"], default="analysis")
    parser.add_argument("--skip-file-check", action="store_true", help="Only validate metadata; downstream file/reference audit is still required.")
    args = parser.parse_args(argv)
    try:
        sample_fields, samples = read_sheet(args.samplesheet, {"sample", "condition", "strandedness"})
        contrast_fields, contrasts = read_sheet(args.contrasts, set(CONTRAST_COLUMNS))
        samplesheet_base = args.samplesheet_base.expanduser().resolve() if args.samplesheet_base else args.samplesheet.resolve().parent
        samples = validate_samples(samples, samplesheet_base, check=not args.skip_file_check, read_length=args.read_length)
        contrasts, comparisons = validate_contrasts(contrasts, samples, design=args.design, analysis_purpose=args.analysis_purpose)
        # One deterministic subject order feeds BAM lists, SUPPA columns and
        # exported PSI identities. Keep conditions contiguous for native ranges.
        paired_groups = {c[key] for c in comparisons if c["paired"] for key in ("treatment", "control")}
        if paired_groups:
            groups = list(dict.fromkeys(row["condition"] for row in samples))
            positions = {row["sample"]: i for i, row in enumerate(samples)}
            samples.sort(key=lambda row: (groups.index(row["condition"]),
                         row["subject"] if row["condition"] in paired_groups else str(positions[row["sample"]]).zfill(12)))
            contrasts, comparisons = validate_contrasts(contrasts, samples, design=args.design, analysis_purpose=args.analysis_purpose)
        references = {}
        for kind, value in (("fasta", args.fasta), ("gtf", args.gtf)):
            resolved = input_path(value, Path.cwd(), check=not args.skip_file_check, context=kind)
            references[kind] = {"path": resolved}
            if Path(resolved).is_file():
                references[kind]["sha256"] = digest(Path(resolved))
        args.outdir.mkdir(parents=True, exist_ok=True)
        write_sheet(args.outdir / "samplesheet.valid.csv", SAMPLE_COLUMNS, sample_fields, samples)
        write_sheet(args.outdir / "contrastsheet.valid.csv", CONTRAST_COLUMNS + ["paired", "design"], contrast_fields, contrasts)
        audit = {"schema_version": "1.0", "status": "metadata_validated", "analysis_purpose": args.analysis_purpose,
                 "assay_contract": {"supported_assay": "bulk_rna_seq", "label_origin": "user_declared_when_supplied_otherwise_default", "identity_verification": "not_inferred_from_input_files"},
                 "biological_validity": "not_established_by_metadata_validation", "experimental_unit_assumption": "one sample per subject and condition; paired contrasts use subject blocks",
                 "smoke_test_only": args.analysis_purpose == "engine_smoke",
                 "file_existence_checked": not args.skip_file_check, "bam_content_audit_pending": any(row["source"] == "genome_bam" for row in samples),
                 "references": references, "samples": samples, "comparisons": comparisons,
                 "input_sheet_sha256": {"samplesheet": digest(args.samplesheet), "contrasts": digest(args.contrasts)}}
        (args.outdir / "input_audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({"status": "metadata_validated", "samples": len(samples), "contrasts": len(contrasts), "outdir": str(args.outdir)}))
        return 0
    except (InputError, OSError, csv.Error) as exc:
        print(f"Input validation failed: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
