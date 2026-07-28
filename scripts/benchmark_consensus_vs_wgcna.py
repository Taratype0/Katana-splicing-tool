from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

from katana_cosplice.benchmark import (
    _fit_method,
    _method_metrics,
    _neutral_graph,
    _rscript_with_wgcna,
    select_benchmark_universe,
)
from katana_cosplice.config import BenchmarkConfig, ConsensusConfig, load_config
from katana_cosplice.consensus import _bootstrap_indices, build_bootstrap_coassignment_consensus
from katana_cosplice.network import discover_module_labels


def _within_coassignment(labels: np.ndarray, coassignment: np.ndarray) -> float:
    weighted_sum = 0.0
    pair_count = 0
    for module in sorted(set(labels).difference({"grey"})):
        indices = np.flatnonzero(labels == module)
        if len(indices) < 2:
            continue
        values = coassignment[np.ix_(indices, indices)][np.triu_indices(len(indices), k=1)]
        weighted_sum += float(values.sum())
        pair_count += len(values)
    return weighted_sum / pair_count if pair_count else float("nan")


def _coassignment(run_labels: list[np.ndarray]) -> np.ndarray:
    n_events = len(run_labels[0])
    counts = np.zeros((n_events, n_events), dtype=np.uint16)
    for labels in run_labels:
        for module in set(labels).difference({"grey"}):
            indices = np.flatnonzero(labels == module)
            counts[np.ix_(indices, indices)] += 1
    result = counts.astype(np.float32) / len(run_labels)
    np.fill_diagonal(result, 1.0)
    return result


def _plot(summary: pd.DataFrame, output: Path) -> None:
    metrics = [
        ("mean_bootstrap_ARI", "Mean bootstrap ARI", True),
        ("weighted_mean_within_module_coassignment", "Within-module co-assignment", True),
        ("within_between_separation", "Within-between correlation separation", True),
        ("grey_fraction", "Grey fraction", False),
        ("largest_module_fraction", "Largest module fraction", False),
        ("runtime_seconds_including_bootstraps", "Runtime (seconds)", False),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.5))
    colors = ["#087E8B", "#D1495B"]
    for ax, (column, title, higher_better) in zip(axes.flat, metrics, strict=True):
        values = summary[column].to_numpy(float)
        bars = ax.bar(summary["method"], values, color=colors)
        for bar, value in zip(bars, values, strict=True):
            ax.text(bar.get_x() + bar.get_width() / 2, value, f"{value:.3g}", ha="center", va="bottom")
        ax.set_title(f"{title}\n({'higher' if higher_better else 'lower'} is better)")
        ax.tick_params(axis="x", rotation=18)
    fig.suptitle("Matched CoSplice consensus vs exact WGCNA", fontsize=16, fontweight="bold")
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Benchmark upgraded CoSplice bootstrap consensus against exact R WGCNA."
    )
    parser.add_argument("config", type=Path)
    parser.add_argument("--events", type=int, default=150)
    parser.add_argument("--bootstraps", type=int, default=20)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    config = load_config(args.config)
    config.validate()
    stable_path = config.output_dir / "03_stable_as_qc" / "stable_psi_matrix.tsv"
    metrics_path = config.output_dir / "03_stable_as_qc" / "stable_as_event_qc.tsv"
    samples_path = config.output_dir / "02_matrix" / "sample_metadata_resolved.tsv"
    if not metrics_path.exists():
        metrics_path = config.output_dir / "03_stable_as_qc" / "stable_as_qc_metrics.tsv"
    matrix = pd.read_csv(stable_path, sep="\t", index_col=0)
    qc = pd.read_csv(metrics_path, sep="\t", low_memory=False)
    samples = pd.read_csv(samples_path, sep="\t", low_memory=False)
    selected, universe = select_benchmark_universe(matrix, qc, args.events)
    observed_path = config.output_dir / "04_network" / "observed" / "event_module_assignments.tsv"
    observed = pd.read_csv(observed_path, sep="\t", low_memory=False)
    observed = observed.loc[observed["event_key"].astype(str).isin(selected.index.astype(str))]

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    benchmark = replace(config.benchmark, max_events=args.events, bootstrap_iterations=args.bootstraps)
    consensus = replace(
        config.consensus,
        bootstrap_iterations=args.bootstraps,
        dense_max_events=max(config.consensus.dense_max_events, len(selected)),
    )
    rscript, audit = _rscript_with_wgcna()
    if rscript is None:
        raise RuntimeError(audit)

    start = time.perf_counter()
    cosplice = build_bootstrap_coassignment_consensus(
        selected, observed, samples, config.network, consensus
    )
    cosplice_seconds = time.perf_counter() - start
    cosplice_labels = cosplice.assignments.set_index("event_key").reindex(selected.index)[
        "consensus_module"
    ].to_numpy(object)

    start = time.perf_counter()
    wgcna_labels, _ = _fit_method("wgcna_exact", selected, config.network, benchmark, rscript)
    wgcna_runs: list[np.ndarray] = []
    cosplice_runs: list[np.ndarray] = []
    stability_rows: list[dict[str, float | int]] = []
    rng = np.random.default_rng(consensus.random_seed)
    partition_settings = replace(
        config.network,
        consensus_bootstrap_iterations=0,
        community_consensus_runs=1,
        bootstrap_iterations=0,
    )
    for iteration in range(1, args.bootstraps + 1):
        indices = _bootstrap_indices(samples, selected.columns, rng)
        bootstrap = selected.iloc[:, indices].copy()
        bootstrap.columns = [f"matched_{iteration:03d}_{position:03d}" for position in range(len(indices))]
        c_labels = discover_module_labels(bootstrap, partition_settings)
        w_labels, _ = _fit_method("wgcna_exact", bootstrap, config.network, benchmark, rscript)
        cosplice_runs.append(c_labels)
        wgcna_runs.append(w_labels)
        stability_rows.append(
            {
                "iteration": iteration,
                "cosplice_consensus_ARI": adjusted_rand_score(cosplice_labels, c_labels),
                "wgcna_exact_ARI": adjusted_rand_score(wgcna_labels, w_labels),
            }
        )
    wgcna_seconds = time.perf_counter() - start
    stability = pd.DataFrame(stability_rows)
    wgcna_coassignment = _coassignment(wgcna_runs)
    neutral = _neutral_graph(selected, benchmark.neutral_knn_neighbors)
    rows: list[dict[str, object]] = []
    for method, labels, coassignment, runtime, ari_column in (
        (
            "Upgraded CoSplice Leiden consensus",
            cosplice_labels,
            cosplice.coassignment_matrix,
            cosplice_seconds,
            "cosplice_consensus_ARI",
        ),
        ("Exact WGCNA", wgcna_labels, wgcna_coassignment, wgcna_seconds, "wgcna_exact_ARI"),
    ):
        row = {"method": method, "n_events": len(selected), "n_samples": selected.shape[1]}
        row.update(_method_metrics(selected, labels, samples, neutral))
        row.update(
            {
                "mean_bootstrap_ARI": float(stability[ari_column].mean()),
                "median_bootstrap_ARI": float(stability[ari_column].median()),
                "min_bootstrap_ARI": float(stability[ari_column].min()),
                "weighted_mean_within_module_coassignment": _within_coassignment(labels, coassignment),
                "runtime_seconds_including_bootstraps": runtime,
            }
        )
        rows.append(row)
    summary = pd.DataFrame(rows)
    assignments = universe[["event_key", "benchmark_rank"]].copy()
    assignments["cosplice_module"] = cosplice_labels
    assignments["wgcna_module"] = wgcna_labels
    summary.to_csv(output / "cosplice_consensus_vs_exact_wgcna_summary.tsv", sep="\t", index=False)
    assignments.to_csv(output / "event_module_assignments.tsv", sep="\t", index=False)
    stability.to_csv(output / "matched_bootstrap_stability.tsv", sep="\t", index=False)
    cosplice.parameter_search.to_csv(output / "cosplice_parameter_search.tsv", sep="\t", index=False)
    cosplice.assignments.to_csv(output / "cosplice_soft_membership.tsv", sep="\t", index=False)
    np.savez_compressed(
        output / "matched_coassignment_matrices.npz",
        event_keys=selected.index.astype(str).to_numpy(),
        cosplice=cosplice.coassignment_matrix,
        wgcna=wgcna_coassignment,
    )
    _plot(summary, output / "cosplice_consensus_vs_exact_wgcna.png")
    print(summary.to_string(index=False))
    print(f"\nExact WGCNA audit: {audit}")
    print(f"Output: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
