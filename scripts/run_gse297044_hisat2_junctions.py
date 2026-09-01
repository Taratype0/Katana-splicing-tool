from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path

import pandas as pd

try:
    from scripts.align_hisat2_stream_junctions import (
        align_and_count,
        alignment_fingerprint,
    )
except ModuleNotFoundError:  # Direct `python scripts/...py` execution.
    from align_hisat2_stream_junctions import align_and_count, alignment_fingerprint


def _valid_cache(output: Path, summary: Path, expected_fingerprint: str) -> bool:
    if not output.exists() or output.stat().st_size == 0 or not summary.exists():
        return False
    parameters = json.loads(summary.read_text(encoding="utf-8"))
    return bool(
        parameters.get("alignments_seen", 0) > 0
        and parameters.get("unique_junctions", 0) > 0
        and parameters.get("fingerprint") == expected_fingerprint
    )


def run_batch(work_root: Path, workers: int, threads_per_sample: int) -> pd.DataFrame:
    metadata = pd.read_csv(
        work_root / "inputs" / "sample_metadata_primary.tsv", sep="\t"
    )
    hisat2 = (
        work_root
        / "tools"
        / "hisat2_windows"
        / "hisat2.1"
        / "hisat2-align-s.exe"
    )
    index_prefix = work_root / "reference" / "hisat2_grch38" / "grch38" / "genome"
    output_root = work_root / "junction_counts"
    log_root = work_root / "logs" / "hisat2_windows"
    output_root.mkdir(parents=True, exist_ok=True)
    log_root.mkdir(parents=True, exist_ok=True)

    def run_sample(sample_id: str) -> dict[str, object]:
        output = output_root / f"{sample_id}.junction_counts.tsv.gz"
        summary = output_root / f"{sample_id}.alignment_summary.json"
        read1 = work_root / "fastq" / f"{sample_id}_1.fastq.gz"
        read2 = work_root / "fastq" / f"{sample_id}_2.fastq.gz"
        expected_fingerprint = alignment_fingerprint(
            executable=hisat2,
            index_prefix=index_prefix,
            read1=read1,
            read2=read2,
            threads=threads_per_sample,
            min_mapq=1,
            min_anchor=8,
            max_intron=500_000,
        )
        if _valid_cache(output, summary, expected_fingerprint):
            cached = json.loads(summary.read_text(encoding="utf-8"))
            cached["sample_id"] = sample_id
            cached["cache_status"] = "cache_hit"
            return cached
        result = align_and_count(
            executable=hisat2,
            index_prefix=index_prefix,
            read1=read1,
            read2=read2,
            output=output,
            summary=summary,
            log=log_root / f"{sample_id}.log",
            threads=threads_per_sample,
            min_mapq=1,
            min_anchor=8,
            max_intron=500_000,
        )
        result["sample_id"] = sample_id
        result["cache_status"] = "rebuilt"
        return result

    rows: list[dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=max(1, workers)) as executor:
        futures = {
            executor.submit(run_sample, str(sample_id)): str(sample_id)
            for sample_id in metadata["sample_id"]
        }
        for future in as_completed(futures):
            sample_id = futures[future]
            result = future.result()
            rows.append(result)
            print(
                f"{sample_id}: {result['cache_status']}; "
                f"junctions={result['unique_junctions']:,}; "
                f"elapsed={result['elapsed_seconds']:.1f}s",
                flush=True,
            )
    summary = pd.DataFrame(rows).sort_values("sample_id")
    summary.to_csv(output_root / "batch_alignment_summary.tsv", sep="\t", index=False)
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run resumable Windows HISAT2 junction counting for GSE297044."
    )
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--threads-per-sample", type=int, default=6)
    args = parser.parse_args(argv)
    summary = run_batch(
        args.work_root.resolve(),
        workers=args.workers,
        threads_per_sample=args.threads_per_sample,
    )
    print(f"Completed {len(summary)} samples.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
