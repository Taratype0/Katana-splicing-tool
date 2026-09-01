from __future__ import annotations

import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import time


CIGAR_PATTERN = re.compile(r"(\d+)([MIDNSHP=X])")
REFERENCE_CONSUMING = {"M", "D", "N", "=", "X"}
ANCHOR_OPERATIONS = {"M", "=", "X"}
BREAK_OPERATIONS = {"N", "S", "H"}


def _file_signature(path: Path) -> dict[str, object]:
    stat = path.stat()
    return {
        "path": str(path.resolve()),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def alignment_fingerprint(
    *,
    executable: Path,
    index_prefix: Path,
    read1: Path,
    read2: Path,
    threads: int,
    min_mapq: int,
    min_anchor: int,
    max_intron: int,
) -> str:
    index_files = sorted(index_prefix.parent.glob(f"{index_prefix.name}.*.ht2*"))
    payload = {
        "schema": 1,
        "executable": _file_signature(executable),
        "index": [_file_signature(path) for path in index_files],
        "read1": _file_signature(read1),
        "read2": _file_signature(read2),
        "threads": max(1, threads),
        "min_mapq": min_mapq,
        "min_anchor": min_anchor,
        "max_intron": max_intron,
    }
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True).encode("utf-8")
    ).hexdigest()


def _materialize_fastq(source: Path, scratch: Path, label: str) -> Path:
    if source.suffix.casefold() != ".gz":
        return source
    scratch.mkdir(parents=True, exist_ok=True)
    target = scratch / f"{label}.fastq"
    partial = target.with_suffix(".fastq.part")
    with gzip.open(source, "rb") as reader, partial.open("wb") as writer:
        shutil.copyfileobj(reader, writer, length=16 * 1024 * 1024)
    partial.replace(target)
    return target


def _junctions(
    chrom: str,
    position_1based: int,
    cigar: str,
    *,
    min_anchor: int,
    max_intron: int,
) -> list[tuple[str, int, int]]:
    operations = [(operation, int(length)) for length, operation in CIGAR_PATTERN.findall(cigar)]
    if not operations:
        return []
    position = position_1based - 1
    output: list[tuple[str, int, int]] = []
    for index, (operation, length) in enumerate(operations):
        if operation == "N":
            left_anchor = 0
            for left_operation, left_length in reversed(operations[:index]):
                if left_operation in ANCHOR_OPERATIONS:
                    left_anchor = left_length
                    break
                if left_operation in BREAK_OPERATIONS:
                    break
            right_anchor = 0
            for right_operation, right_length in operations[index + 1 :]:
                if right_operation in ANCHOR_OPERATIONS:
                    right_anchor = right_length
                    break
                if right_operation in BREAK_OPERATIONS:
                    break
            if min(left_anchor, right_anchor) >= min_anchor and 0 < length <= max_intron:
                output.append((chrom, position, position + length))
        if operation in REFERENCE_CONSUMING:
            position += length
    return output


def align_and_count(
    *,
    executable: Path,
    index_prefix: Path,
    read1: Path,
    read2: Path,
    output: Path,
    summary: Path,
    log: Path,
    threads: int,
    min_mapq: int,
    min_anchor: int,
    max_intron: int,
) -> dict[str, object]:
    output.parent.mkdir(parents=True, exist_ok=True)
    summary.parent.mkdir(parents=True, exist_ok=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    scratch = output.parent / ".hisat2_fastq_scratch" / output.name.removesuffix(
        ".junction_counts.tsv.gz"
    )
    aligned_read1 = _materialize_fastq(read1, scratch, "read1")
    aligned_read2 = _materialize_fastq(read2, scratch, "read2")
    source_fingerprint = alignment_fingerprint(
        executable=executable,
        index_prefix=index_prefix,
        read1=read1,
        read2=read2,
        threads=threads,
        min_mapq=min_mapq,
        min_anchor=min_anchor,
        max_intron=max_intron,
    )
    command = [
        str(executable),
        "--wrapper",
        "basic-0",
        "-x",
        str(index_prefix),
        "-1",
        str(aligned_read1),
        "-2",
        str(aligned_read2),
        "-p",
        str(max(1, threads)),
        "--add-chrname",
        "--dta",
        "--no-mixed",
        "--no-discordant",
        "--new-summary",
    ]
    started = time.perf_counter()
    counts: Counter[tuple[str, int, int]] = Counter()
    alignments_seen = 0
    alignments_retained = 0
    spliced_alignments = 0
    try:
        with log.open("w", encoding="utf-8") as stderr:
            process = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=stderr,
                text=True,
                encoding="utf-8",
                errors="replace",
                bufsize=1024 * 1024,
            )
            assert process.stdout is not None
            for line in process.stdout:
                if line.startswith("@"):
                    continue
                fields = line.rstrip("\n").split("\t", 11)
                if len(fields) < 6:
                    continue
                alignments_seen += 1
                flag = int(fields[1])
                if flag & (0x4 | 0x100 | 0x200 | 0x400 | 0x800):
                    continue
                mapq = int(fields[4])
                if mapq < min_mapq:
                    continue
                alignments_retained += 1
                junctions = _junctions(
                    fields[2],
                    int(fields[3]),
                    fields[5],
                    min_anchor=min_anchor,
                    max_intron=max_intron,
                )
                if junctions:
                    spliced_alignments += 1
                    counts.update(junctions)
            return_code = process.wait()
    finally:
        for materialized, original in (
            (aligned_read1, read1),
            (aligned_read2, read2),
        ):
            if materialized != original:
                materialized.unlink(missing_ok=True)
        if scratch.exists() and not any(scratch.iterdir()):
            scratch.rmdir()
    if return_code != 0:
        raise RuntimeError(
            f"HISAT2 failed with return code {return_code}; inspect {log}"
        )

    partial = output.with_suffix(output.suffix + ".part")
    with gzip.open(partial, "wt", encoding="utf-8", newline="") as handle:
        handle.write("chrom\tintron_start_0based\tintron_end_0based\tread_count\n")
        for (chrom, start, end), count in sorted(counts.items()):
            handle.write(f"{chrom}\t{start}\t{end}\t{count}\n")
    partial.replace(output)
    result = {
        "aligner": "HISAT2 2.1 Windows",
        "fingerprint": source_fingerprint,
        "command": command,
        "read1": str(read1.resolve()),
        "read2": str(read2.resolve()),
        "alignments_seen": alignments_seen,
        "alignments_retained": alignments_retained,
        "spliced_alignments": spliced_alignments,
        "unique_junctions": len(counts),
        "junction_read_observations": int(sum(counts.values())),
        "elapsed_seconds": round(time.perf_counter() - started, 3),
        "filters": {
            "min_mapping_quality": min_mapq,
            "min_anchor_bp": min_anchor,
            "max_intron_length": max_intron,
            "excluded_flags": ["unmapped", "secondary", "qcfail", "duplicate", "supplementary"],
        },
    }
    summary.write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Stream HISAT2 SAM output directly into sample-level junction counts."
    )
    parser.add_argument("--hisat2", type=Path, required=True)
    parser.add_argument("--index-prefix", type=Path, required=True)
    parser.add_argument("--read1", type=Path, required=True)
    parser.add_argument("--read2", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=12)
    parser.add_argument("--min-mapq", type=int, default=1)
    parser.add_argument("--min-anchor", type=int, default=8)
    parser.add_argument("--max-intron", type=int, default=500_000)
    args = parser.parse_args(argv)
    result = align_and_count(
        executable=args.hisat2.resolve(),
        index_prefix=args.index_prefix.resolve(),
        read1=args.read1.resolve(),
        read2=args.read2.resolve(),
        output=args.output.resolve(),
        summary=args.summary.resolve(),
        log=args.log.resolve(),
        threads=args.threads,
        min_mapq=args.min_mapq,
        min_anchor=args.min_anchor,
        max_intron=args.max_intron,
    )
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
