from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from katana_cosplice.benchmark import (
    render_benchmark_plots,
    render_benchmark_report,
    run_method_benchmark,
)
from katana_cosplice.config import BenchmarkConfig, NetworkConfig
from katana_cosplice.hybrid_benchmark import (
    build_hybrid_benchmark,
    write_hybrid_benchmark,
)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run same-input CoSplice and published-method baselines on a public event matrix."
    )
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--events", type=int, default=1000)
    parser.add_argument("--bootstraps", type=int, default=20)
    parser.add_argument("--min-module-size", type=int, default=20)
    args = parser.parse_args()

    matrix = pd.read_csv(args.matrix, sep="\t", index_col=0, low_memory=False)
    matrix = matrix.apply(pd.to_numeric, errors="coerce")
    matrix = matrix.loc[matrix.notna().mean(axis=1).ge(0.80)]
    matrix = matrix.apply(lambda row: row.fillna(row.median()), axis=1)
    matrix = matrix.loc[matrix.std(axis=1).gt(1e-8)]
    if len(matrix) > args.events:
        matrix = matrix.loc[matrix.var(axis=1).nlargest(args.events).index]
    samples = pd.read_csv(args.metadata, sep="\t")
    missing_samples = set(matrix.columns).difference(samples["sample_id"].astype(str))
    if missing_samples:
        raise ValueError(f"Metadata is missing {len(missing_samples)} matrix samples.")
    samples = samples.set_index(samples["sample_id"].astype(str)).loc[matrix.columns].reset_index(
        drop=True
    )
    qc = pd.DataFrame(
        {
            "event_key": matrix.index.astype(str),
            "stability_score": matrix.var(axis=1).to_numpy(),
            "selection_note": "Discovery-dataset variance ranking shared by every method.",
        }
    )
    network = NetworkConfig(
        min_module_size=args.min_module_size,
        knn_neighbors=min(50, max(10, len(matrix) // 20)),
        community_resolution=1.0,
        module_refinement=True,
        min_module_membership=0.45,
        module_merge_correlation=0.85,
        bootstrap_iterations=0,
    )
    benchmark = BenchmarkConfig(
        methods=[
            "katana_sparse",
            "pearson_hierarchical",
            "graphical_lasso",
            "wgcna_exact",
        ],
        max_events=len(matrix),
        bootstrap_iterations=args.bootstraps,
        random_seed=41,
    )
    args.output.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    result = run_method_benchmark(matrix, qc, samples, network, benchmark, progress=print)
    elapsed = time.perf_counter() - started
    outputs = {
        "common_event_universe.tsv": result.universe,
        "event_method_assignments.tsv": result.assignments,
        "method_benchmark_summary.tsv": result.method_summary,
        "bootstrap_partition_stability.tsv": result.bootstrap_stability,
        "module_matching.tsv": result.module_matching,
        "method_availability.tsv": result.method_audit,
        "method_metric_ranking.tsv": result.metric_ranking,
        "method_strengths_limitations.tsv": result.method_interpretation,
        "sample_metadata.tsv": samples,
        "external_event_qc.tsv": qc,
    }
    for name, frame in outputs.items():
        frame.to_csv(args.output / name, sep="\t", index=False)
    render_benchmark_plots(result, args.output / "plots")
    render_benchmark_report(result, args.output / "method_benchmark_report.html")
    hybrid = build_hybrid_benchmark(result.assignments, result.method_summary)
    hybrid_files = write_hybrid_benchmark(hybrid, args.output)
    (args.output / "parameters.json").write_text(
        json.dumps(
            {
                "dataset_id": args.dataset_id,
                "matrix": str(args.matrix.resolve()),
                "metadata": str(args.metadata.resolve()),
                "selection_scope": "one shared variance-ranked event universe",
                "n_events": len(matrix),
                "n_samples": matrix.shape[1],
                "network": network.__dict__,
                "benchmark": benchmark.__dict__,
                "hybrid_benchmark": {
                    "tight_branch": "katana_sparse",
                    "diffuse_branch": "wgcna_exact",
                    "fusion_policy": "parallel views; no artificial single hybrid partition",
                    "generated_files": [str(path) for path in hybrid_files],
                },
                "wall_clock_seconds": elapsed,
                "claim_limit": (
                    "Same-input computational benchmark. Single-partition metrics compare the "
                    "tight branch with published baselines. The hybrid audit reports complementary "
                    "tight/diffuse coverage and does not invent a fused partition. Biological "
                    "validation and independent cohort preservation remain separate requirements."
                ),
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
