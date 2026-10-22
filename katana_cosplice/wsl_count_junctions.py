"""Standalone pysam junction counter used by Katana's Windows WSL backend."""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
import tempfile
import time
from collections import Counter
from pathlib import Path

import pysam


def _junctions(record: object, min_anchor: int, max_intron: int) -> list[tuple[str, int, int]]:
    chrom = str(record.reference_name)
    position = int(record.reference_start)
    operations = list(record.cigartuples or [])
    output: list[tuple[str, int, int]] = []
    for index, (operation, length) in enumerate(operations):
        length = int(length)
        if operation == 3:
            left_anchor = 0
            for i in range(index - 1, -1, -1):
                if operations[i][0] in {0, 7, 8}:
                    left_anchor = int(operations[i][1])
                    break
                if operations[i][0] in {3, 4, 5}:
                    break
            right_anchor = 0
            for i in range(index + 1, len(operations)):
                if operations[i][0] in {0, 7, 8}:
                    right_anchor = int(operations[i][1])
                    break
                if operations[i][0] in {3, 4, 5}:
                    break
            if min(left_anchor, right_anchor) >= min_anchor and 0 < length <= max_intron:
                output.append((chrom, position, position + length))
            position += length
        elif operation in {0, 2, 7, 8}:
            position += length
    return output


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--min-mapq", type=int, default=1)
    parser.add_argument("--min-anchor", type=int, default=8)
    parser.add_argument("--max-intron", type=int, default=500_000)
    parser.add_argument("--include-duplicates", action="store_true")
    parser.add_argument("--scratch-copy", action="store_true")
    args = parser.parse_args()

    source = Path(args.bam)
    working_bam = source
    scratch_dir: Path | None = None
    started = time.time()
    if args.scratch_copy:
        scratch_dir = Path(tempfile.mkdtemp(prefix="katana-cosplice-"))
        working_bam = scratch_dir / source.name
        shutil.copy2(source, working_bam)

    counts: Counter[tuple[str, int, int]] = Counter()
    alignments_seen = 0
    spliced_alignments = 0
    try:
        with pysam.AlignmentFile(str(working_bam), "rb", check_sq=False) as handle:
            for record in handle.fetch(until_eof=True):
                alignments_seen += 1
                if (
                    record.is_unmapped
                    or record.is_secondary
                    or record.is_supplementary
                    or record.is_qcfail
                    or (record.is_duplicate and not args.include_duplicates)
                    or int(record.mapping_quality) < args.min_mapq
                ):
                    continue
                junctions = _junctions(record, args.min_anchor, args.max_intron)
                if junctions:
                    spliced_alignments += 1
                    counts.update(junctions)
    finally:
        if scratch_dir is not None:
            shutil.rmtree(scratch_dir, ignore_errors=True)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output, "wt", encoding="utf-8", newline="") as handle:
        handle.write("chrom\tintron_start_0based\tintron_end_0based\tread_count\n")
        for (chrom, start, end), count in sorted(counts.items()):
            handle.write(f"{chrom}\t{start}\t{end}\t{count}\n")
    summary = {
        "alignments_seen": alignments_seen,
        "spliced_alignments": spliced_alignments,
        "unique_junctions": len(counts),
        "junction_read_observations": int(sum(counts.values())),
        "elapsed_seconds": time.time() - started,
        "scratch_copy": bool(args.scratch_copy),
    }
    Path(args.summary).write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
