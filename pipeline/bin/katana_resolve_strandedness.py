#!/usr/bin/env python3
"""Resolve auto libraries from native RSeQC evidence before strand-aware analysis."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
from pathlib import Path

from katana_validate_inputs import InputError, SAMPLE_COLUMNS, boolean, read_sheet, validate_contrasts, write_sheet

RESOLVED = frozenset({"forward", "reverse", "unstranded"})
RULES = {
    True: ("SingleEnd", "++,--", "+-,-+"),
    False: ("PairEnd", "1++,1--,2+-,2-+", "1+-,1-+,2++,2--"),
}


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def thresholds(stranded: float, unstranded: float) -> None:
    if not math.isfinite(stranded) or not 0.5 < stranded <= 1:
        raise InputError("stranded_threshold must be finite, > 0.5 and <= 1")
    if not math.isfinite(unstranded) or not 0 <= unstranded < 1:
        raise InputError("unstranded_threshold must be finite, >= 0 and < 1")
    if unstranded >= 2 * stranded - 1:
        raise InputError("Stranded and unstranded threshold regions overlap")


def infer(text: str, *, single_end: bool, stranded: float = 0.8, unstranded: float = 0.1) -> dict:
    """Classify normalized informative-read fractions, retaining native fractions.

    RSeQC signs describe read vs annotated gene, not genomic plus/minus totals.
    Paired 1++,1--,2+-,2-+ (or single ++,--) is forward/second-strand.
    Opposite signs are reverse/first-strand. Failure fraction is not evidence
    for an unstranded library. Ambiguous/empty/malformed evidence fails closed.
    """
    thresholds(stranded, unstranded)
    layout, forward_rule, reverse_rule = RULES[single_end]
    layouts = re.findall(r"^This is (SingleEnd|PairEnd) Data\s*$", text, re.M)
    if layouts != [layout]:
        raise InputError(f"RSeQC layout does not match declared {layout}: {layouts}")
    extracted = re.findall(r'^Fraction of reads explained by "([^"]+)":\s*(\S+)\s*$', text, re.M)
    failed = re.findall(r"^Fraction of reads failed to determine:\s*(\S+)\s*$", text, re.M)
    if len(extracted) != 2 or {r[0] for r in extracted} != {forward_rule, reverse_rule} or len(failed) != 1:
        raise InputError("Incomplete or duplicate native RSeQC strand fractions")
    try:
        fractions = {rule: float(value) for rule, value in extracted}
        failed_fraction = float(failed[0])
    except ValueError as exc:
        raise InputError("Invalid RSeQC fraction") from exc
    values = [*fractions.values(), failed_fraction]
    if any(not math.isfinite(v) or not 0 <= v <= 1 for v in values) or abs(sum(values) - 1) > 0.002:
        raise InputError("RSeQC fractions must be finite probabilities summing to one")
    informative = fractions[forward_rule] + fractions[reverse_rule]
    if informative <= 0:
        raise InputError("RSeQC found no informative strand evidence; supply the known library type or inspect annotation/mapping")
    forward = fractions[forward_rule] / informative
    reverse = fractions[reverse_rule] / informative
    if forward >= stranded:
        result = "forward"
    elif reverse >= stranded:
        result = "reverse"
    elif abs(forward - reverse) <= unstranded + 1e-12:
        result = "unstranded"
    else:
        raise InputError(f"Undetermined strandedness: informative forward={forward:.6f}, reverse={reverse:.6f}; inspect native evidence and set an explicit library type")
    return {
        "strandedness": result, "method": "rseqc_infer_experiment",
        "single_end": single_end, "forward_rule": forward_rule, "reverse_rule": reverse_rule,
        "native_forward_fraction": fractions[forward_rule], "native_reverse_fraction": fractions[reverse_rule],
        "failed_to_determine_fraction": failed_fraction, "informative_fraction": informative,
        "informative_forward_fraction": forward, "informative_reverse_fraction": reverse,
        "stranded_threshold": stranded, "unstranded_threshold": unstranded,
        "threshold_basis": "nf-core/rnaseq strand classification; fractions normalized over informative reads",
    }


def resolve(samplesheet: Path, contrastsheet: Path, evidence: dict, outdir: Path, *, stranded: float = 0.8,
            unstranded: float = 0.1, design: str = "~ condition", analysis_purpose: str = "analysis") -> dict:
    thresholds(stranded, unstranded)
    fields, rows = read_sheet(samplesheet, {"sample", "condition", "strandedness", "single_end", "read_length"})
    _, contrasts = read_sheet(contrastsheet, {"contrast", "treatment", "control"})
    ids = [r["sample"] for r in rows]
    if len(ids) != len(set(ids)):
        raise InputError("Duplicate sample in validated samplesheet")
    expected = {r["sample"] for r in rows if r["strandedness"] == "auto"}
    if set(evidence) != expected:
        raise InputError(f"Native inference evidence must match auto samples exactly; missing={sorted(expected-set(evidence))}, unexpected={sorted(set(evidence)-expected)}")
    audits = []
    for row in rows:
        for field in SAMPLE_COLUMNS:
            row.setdefault(field, "")
        sample = row["sample"]
        requested = row["strandedness"]
        if requested == "auto":
            record = evidence[sample]
            if not isinstance(record, dict) or set(record) != {"report", "log"}:
                raise InputError(f"{sample}: native report/log binding is required")
            report, log = (Path(record[k]) for k in ("report", "log"))
            if not report.is_file() or not log.is_file():
                raise InputError(f"{sample}: native RSeQC evidence is missing")
            audit = infer(report.read_text(encoding="utf-8"), single_end=boolean(row["single_end"], context=sample), stranded=stranded, unstranded=unstranded)
            audit["native_report"] = {"filename": report.name, "sha256": file_hash(report)}
            audit["native_log"] = {"filename": log.name, "sha256": file_hash(log)}
            sampled = re.findall(r"Total\s+(\d+)\s+usable reads were sampled", log.read_text(encoding="utf-8"))
            if len(sampled) != 1 or int(sampled[0]) == 0:
                raise InputError(f"{sample}: missing, duplicate or zero usable-read count in native RSeQC log")
            audit["usable_reads_sampled"] = int(sampled[0])
        elif requested in RESOLVED:
            audit = {"strandedness": requested, "method": "user_declared_not_inferred"}
        else:
            raise InputError(f"{sample}: unresolved/invalid strandedness {requested!r}")
        row["strandedness_requested"] = requested
        row["strandedness"] = audit["strandedness"]
        row["strandedness_method"] = audit["method"]
        audits.append({"sample": sample, "requested": requested, **audit})
    # rMATS uses a single libType per selected contrast. Resolve all per-sample
    # auto values first, then reject true mixed libraries rather than pool them.
    validate_contrasts(contrasts, rows, design=design, analysis_purpose=analysis_purpose)
    if outdir.exists():
        raise InputError(f"Resolution output already exists: {outdir}")
    outdir.mkdir(parents=True)
    priority = SAMPLE_COLUMNS + ["strandedness_requested", "strandedness_method"]
    write_sheet(outdir / "samplesheet.resolved.csv", priority, fields, rows)
    receipt = {
        "schema": "katana.strandedness_resolution/v1", "status": "resolved_and_contrasts_validated",
        "samples": audits, "auto_samples": len(expected), "declared_samples": len(rows)-len(expected),
        "ambiguous_policy": "fail_without_unstranded_fallback", "scientific_models_refitted": False,
        "input_samplesheet_sha256": file_hash(samplesheet), "contrastsheet_sha256": file_hash(contrastsheet),
    }
    (outdir / "strandedness.audit.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    with (outdir / "strandedness.summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        columns = ["sample", "requested", "strandedness", "method", "informative_forward_fraction", "informative_reverse_fraction", "failed_to_determine_fraction", "usable_reads_sampled"]
        writer = csv.DictWriter(handle, columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(audits)
    return receipt


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samplesheet", type=Path, required=True)
    parser.add_argument("--contrasts", type=Path, required=True)
    parser.add_argument("--evidence-map", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--stranded-threshold", type=float, default=0.8)
    parser.add_argument("--unstranded-threshold", type=float, default=0.1)
    parser.add_argument("--design", default="~ condition")
    parser.add_argument("--analysis-purpose", choices=["analysis", "engine_smoke"], default="analysis")
    args = parser.parse_args(argv)
    try:
        evidence = json.loads(args.evidence_map.read_text(encoding="utf-8"))
        if not isinstance(evidence, dict):
            raise InputError("Inference map must be a sample-to-native-evidence object")
        receipt = resolve(args.samplesheet, args.contrasts, evidence, args.outdir, stranded=args.stranded_threshold,
                          unstranded=args.unstranded_threshold, design=args.design, analysis_purpose=args.analysis_purpose)
        print(json.dumps({"status": receipt["status"], "auto_samples": receipt["auto_samples"]}))
        return 0
    except (InputError, OSError, json.JSONDecodeError) as exc:
        print(f"Strandedness resolution failed: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
