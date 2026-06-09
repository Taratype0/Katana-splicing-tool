from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from src.adapters.native_sashimi_adapter import NativeSashimiAdapter


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run Katana-native sashimi plots from a selected-event manifest."
    )
    parser.add_argument("--manifest", required=True, help="TSV manifest with event_file, b1_txt, b2_txt, label1, label2 and outdir columns.")
    parser.add_argument("--output-root", required=True, help="Katana output root or standalone output root.")
    parser.add_argument("--gtf", default="", help="Optional GTF used to enrich exon labels in the event region.")
    parser.add_argument("--samtools", default="samtools", help="samtools executable path or command name.")
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest, sep="\t", low_memory=False)
    adapter = NativeSashimiAdapter(args.samtools)
    logs, failures = adapter.run_manifest(
        manifest,
        Path(args.output_root),
        gtf_path=Path(args.gtf) if args.gtf else None,
    )
    if logs:
        print(logs)
    if not failures.empty:
        failed_path = Path(args.output_root) / "06_sashimi" / "failed_native_sashimi.tsv"
        failed_path.parent.mkdir(parents=True, exist_ok=True)
        failures.to_csv(failed_path, sep="\t", index=False)
        print(f"Failed native sashimi events: {len(failures)}")
        print(failed_path)
        raise SystemExit(1)
    print("Katana native sashimi finished.")


if __name__ == "__main__":
    main()
