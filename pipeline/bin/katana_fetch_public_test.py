#!/usr/bin/env python3
"""Fetch a revision- and SHA256-locked public chrX fixture outside workflow source."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import sys
import urllib.request
from collections import Counter
from pathlib import Path


class FixtureError(ValueError):
    pass


def digest(path: Path) -> str:
    checksum = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def checked_path(root: Path, relative: str) -> Path:
    candidate = (root / relative).resolve()
    if candidate == root or root not in candidate.parents:
        raise FixtureError(f"Unsafe fixture path: {relative}")
    return candidate


def fetch(root: Path, record: dict, verify_only: bool) -> Path:
    path = checked_path(root, record["path"])
    if path.exists():
        if path.stat().st_size != record["bytes"] or digest(path) != record["sha256"]:
            raise FixtureError(f"Existing file differs from frozen fixture; preserved unchanged: {path}")
        return path
    if verify_only:
        raise FixtureError(f"Fixture file missing: {path}")
    if not record["url"].startswith("https://raw.githubusercontent.com/nf-core/test-datasets/"):
        raise FixtureError("Only HTTPS nf-core test-dataset files are accepted by this helper")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".part")
    if temporary.exists():
        raise FixtureError(f"Incomplete earlier download exists; inspect or remove this file: {temporary}")
    request = urllib.request.Request(record["url"], headers={"User-Agent": "katana-public-smoke/0.1"})
    try:
        with urllib.request.urlopen(request, timeout=120) as response, temporary.open("xb") as handle:
            while block := response.read(1024 * 1024):
                handle.write(block)
        if temporary.stat().st_size != record["bytes"] or digest(temporary) != record["sha256"]:
            raise FixtureError(f"Downloaded file fails frozen size/SHA256 verification: {record['path']}")
        temporary.replace(path)
    except Exception:
        if temporary.is_file():
            temporary.unlink()
        raise
    return path


def fastq_records(path: Path):
    with gzip.open(path, "rt", encoding="ascii") as handle:
        while header := handle.readline():
            sequence = handle.readline().rstrip("\r\n")
            plus = handle.readline()
            quality = handle.readline().rstrip("\r\n")
            if not header.startswith("@") or not plus.startswith("+") or not sequence or len(sequence) != len(quality):
                raise FixtureError(f"Malformed FASTQ record in {path}")
            identity = header.split()[0][1:]
            if identity.endswith(("/1", "/2")):
                identity = identity[:-2]
            yield identity, len(sequence)


def inspect_pair(first: Path, second: Path) -> dict:
    from itertools import zip_longest
    lengths: Counter[int] = Counter()
    count = 0
    for left, right in zip_longest(fastq_records(first), fastq_records(second)):
        if left is None or right is None or left[0] != right[0]:
            raise FixtureError(f"Read-pair identity/order mismatch: {first.name}, {second.name}")
        lengths.update([left[1], right[1]])
        count += 1
    if count == 0 or len(lengths) != 1:
        raise FixtureError(f"This frozen rMATS fixture needs nonempty, uniform-length paired reads: {dict(lengths)}")
    return {"read_pairs": count, "read_length": next(iter(lengths)), "mate_order_verified": True}


def write_generated(path: Path, content: str) -> None:
    if path.exists():
        if path.read_text(encoding="utf-8") != content:
            raise FixtureError(f"Existing generated metadata differs; preserved unchanged: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8", newline="\n")


def csv_content(columns: list[str], records: list[dict]) -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, columns, lineterminator="\n")
    writer.writeheader()
    writer.writerows(records)
    return stream.getvalue()


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    source = Path(__file__).resolve().parents[1]
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--lock", type=Path, default=source / "assets" / "public_test.lock.json")
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args(argv)
    try:
        root = args.output_dir.expanduser().resolve()
        if root == source or source in root.parents:
            raise FixtureError("Test reads and references must live outside the versioned workflow source directory")
        lock = json.loads(args.lock.read_text(encoding="utf-8"))
        if lock["purpose"] != "engine_smoke":
            raise FixtureError("This helper is restricted to a software smoke fixture")
        for record in lock["files"]:
            if f"/{lock['revision']}/" not in record["url"]:
                raise FixtureError("Fixture URL is not locked to its declared revision")
            fetch(root, record, args.verify_only)
        samples, read_audit = [], {}
        for row in lock["samples"]:
            first, second = (checked_path(root, row[key]) for key in ("fastq_1", "fastq_2"))
            observed = inspect_pair(first, second)
            read_audit[row["sample"]] = observed
            samples.append({"sample": row["sample"], "condition": row["condition"], "source": "fastq", "fastq_1": str(first), "fastq_2": str(second), "genome_bam": "", "transcriptome_bam": "", "strandedness": row["strandedness"], "single_end": "false", "read_length": observed["read_length"], "subject": row["sample"], "batch": ""})
        lengths = {value["read_length"] for value in read_audit.values()}
        if len(lengths) != 1:
            raise FixtureError("rMATS shared read_length cannot represent this mixed-length fixture")
        columns = ["sample", "condition", "source", "fastq_1", "fastq_2", "genome_bam", "transcriptome_bam", "strandedness", "single_end", "read_length", "subject", "batch"]
        write_generated(root / "samplesheet.csv", csv_content(columns, samples))
        write_generated(root / "contrastsheet.csv", csv_content(["contrast", "treatment", "control"], lock["contrasts"]))
        parameters = {"input": str(root / "samplesheet.csv"), "contrasts": str(root / "contrastsheet.csv"), "fasta": str(checked_path(root, lock["references"]["fasta"])), "gtf": str(checked_path(root, lock["references"]["gtf"])), "public_test_dir": str(root), "read_length": next(iter(lengths)), "analysis_purpose": "engine_smoke", "design": "~ condition", "star_sjdb_overhang": next(iter(lengths)) - 1}
        write_generated(root / "fixture.params.json", json.dumps(parameters, indent=2) + "\n")
        audit = {"schema_version": "1.0", "fixture": lock["fixture"], "revision": lock["revision"], "lock_sha256": digest(args.lock), "status": "download_and_fastq_structure_verified", "purpose": "engine_smoke", "biological_validity": lock["biological_validity"], "reference_scope": "chrX subset, insufficient for genome-wide biological inference", "read_audit": read_audit, "files": lock["files"]}
        write_generated(root / "fixture.audit.json", json.dumps(audit, indent=2) + "\n")
        print(json.dumps({"status": audit["status"], "samples": len(samples), "read_length": next(iter(lengths)), "total_bytes": sum(record["bytes"] for record in lock["files"]), "output_dir": str(root), "purpose": "engine_smoke"}))
        return 0
    except (FixtureError, OSError, ValueError, KeyError) as error:
        print(f"Public fixture preparation failed: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
