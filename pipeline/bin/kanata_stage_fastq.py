#!/usr/bin/env python3
"""Stage one lane by symlink, or concatenate compressed lanes per sample.

The first 100 records are a bounded format/pair probe, not a full FASTQ audit.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import shutil
import sys
from pathlib import Path

SAFE_SAMPLE = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]*$")
PROBE_RECORDS = 100


class StagingError(ValueError):
    pass


def probe_fastq(path: Path, limit: int = PROBE_RECORDS) -> dict:
    """Read only a bounded prefix; gzip trailers beyond it are not checked."""
    if not path.is_file() or path.stat().st_size == 0:
        raise StagingError(f"Missing or empty FASTQ: {path}")
    with path.open("rb") as handle:
        if handle.read(2) != b"\x1f\x8b":
            raise StagingError(f"FASTQ must be gzip compressed: {path}")
    identifiers: list[str] = []
    lengths: list[int] = []
    try:
        with gzip.open(path, "rt", encoding="ascii", newline="") as handle:
            for position in range(1, limit + 1):
                header = handle.readline()
                if not header:
                    break
                sequence, plus, quality = (handle.readline() for _ in range(3))
                if not sequence or not plus or not quality:
                    raise StagingError(f"Incomplete FASTQ record {position}: {path}")
                header, sequence, plus, quality = (
                    line.rstrip("\r\n") for line in (header, sequence, plus, quality)
                )
                if not header.startswith("@") or not header[1:].strip():
                    raise StagingError(f"Invalid FASTQ header at record {position}: {path}")
                if not plus.startswith("+"):
                    raise StagingError(f"Invalid FASTQ separator at record {position}: {path}")
                if not sequence or any(character.isspace() for character in sequence):
                    raise StagingError(f"Empty or whitespace-containing sequence at record {position}: {path}")
                if len(sequence) != len(quality) or any(not 33 <= ord(character) <= 126 for character in quality):
                    raise StagingError(f"Invalid FASTQ sequence/quality at record {position}: {path}")
                identifier = header[1:].split()[0]
                if identifier.endswith(("/1", "/2")):
                    identifier = identifier[:-2]
                if not identifier:
                    raise StagingError(f"Empty read identifier at record {position}: {path}")
                identifiers.append(identifier)
                lengths.append(len(sequence))
    except (OSError, EOFError, UnicodeError) as exc:
        raise StagingError(f"Cannot read gzip FASTQ prefix {path}: {exc}") from exc
    if not identifiers:
        raise StagingError(f"FASTQ contains no records: {path}")
    return {"identifiers": identifiers, "records_probed": len(identifiers),
            "sequence_lengths": lengths}


def stage_reads(request: dict, *, sample: str, output_dir: Path) -> dict:
    if not SAFE_SAMPLE.fullmatch(sample):
        raise StagingError("Sample ID must start with a letter and use letters, digits, dots, underscores or hyphens")
    if not isinstance(request, dict) or not isinstance(request.get("single_end"), bool):
        raise StagingError("FASTQ request must contain boolean single_end")
    reads = request.get("reads")
    if not isinstance(reads, list) or not all(isinstance(value, str) and value for value in reads):
        raise StagingError("FASTQ request must contain a reads list of nonempty paths")
    expected = 1 if request["single_end"] else 2
    lane_count = request.get("lane_count", 1)
    if type(lane_count) is not int or lane_count < 1:
        raise StagingError("lane_count must be a positive integer")
    if len(reads) != expected * lane_count:
        raise StagingError(f"single_end={request['single_end']} requires {expected} FASTQ file(s) per lane")
    sources = [Path(os.path.abspath(value)) for value in reads]
    if len({path.resolve() for path in sources}) != len(sources):
        raise StagingError("The same FASTQ cannot be assigned to both mates")
    output_dir = output_dir.resolve()
    outputs = [output_dir / f"{sample}_R{number}.fastq.gz" for number in range(1, expected + 1)]
    audit_path = output_dir / "fastq_staging_audit.json"
    temporary = [path.with_suffix(path.suffix + '.partial') for path in outputs]
    for destination in outputs + temporary + [audit_path]:
        if os.path.lexists(destination):
            raise StagingError(f"Refusing to overwrite existing staging output: {destination}")
    probes = [probe_fastq(path) for path in sources]
    for lane in range(lane_count):
        offset = lane * expected
        if expected == 2 and probes[offset]["identifiers"] != probes[offset + 1]["identifiers"]:
            raise StagingError(f"Mate identifiers or record counts differ within the bounded FASTQ prefix for lane {lane + 1}")
    audit = {
        "schema_version": "kanata.fastq_staging.v1",
        "status": "staged", "sample": sample, "single_end": request["single_end"],
        "operation": "relative_symlinks_without_raw_data_copy" if lane_count == 1 else "concatenate_gzip_members_in_lane_order",
        "lane_count": lane_count,
        "probe_scope": "bounded_not_full", "probe_limit_records": PROBE_RECORDS,
        "full_fastq_integrity": "not_established", "full_mate_matching": "not_established",
        "reads": [
            {"mate_number": number % expected + 1, "lane_number": number // expected + 1, "staged_source_path": str(source),
             "resolved_source_path": str(source.resolve()), "source_size_bytes": source.stat().st_size,
             "canonical_filename": outputs[number % expected].name, "records_probed": probe["records_probed"],
             "observed_prefix_sequence_lengths": probe["sequence_lengths"]}
            for number, (source, probe) in enumerate(zip(sources, probes))
        ],
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    created: list[Path] = []
    try:
        for mate, destination in enumerate(outputs):
            if lane_count == 1:
                destination.symlink_to(os.path.relpath(sources[mate], output_dir))
                created.append(destination)
            else:
                part = temporary[mate]
                with part.open('xb') as output:
                    created.append(part)
                    for source in sources[mate::expected]:
                        with source.open('rb') as handle:
                            shutil.copyfileobj(handle, output, length=8 * 1024 * 1024)
                if part.stat().st_size != sum(source.stat().st_size for source in sources[mate::expected]):
                    raise StagingError('Lane inputs changed size while merging')
                part.rename(destination)
                created.remove(part)
                created.append(destination)
        with audit_path.open("x", encoding="utf-8") as handle:
            created.append(audit_path)
            handle.write(json.dumps(audit, ensure_ascii=False, indent=2) + "\n")
    except (OSError, StagingError):
        for destination in reversed(created):
            destination.unlink()
        raise
    return audit


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reads-json", type=Path, required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    args = parser.parse_args(argv)
    try:
        request = json.loads(args.reads_json.read_text(encoding="utf-8"))
        audit = stage_reads(request, sample=args.sample, output_dir=args.output_dir)
        print(json.dumps({"status": audit["status"], "sample": audit["sample"],
                          "probe_scope": audit["probe_scope"], "reads": len(audit["reads"])}))
        return 0
    except (OSError, ValueError, TypeError) as exc:
        print(f"FASTQ staging failed: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
