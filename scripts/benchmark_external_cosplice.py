from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.metrics import adjusted_rand_score

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from katana_cosplice.benchmark import (
    render_benchmark_plots,
    render_benchmark_report,
    run_method_benchmark,
)
from katana_cosplice.config import BenchmarkConfig, NetworkConfig


METHODS = ["katana_sparse", "pearson_hierarchical", "graphical_lasso", "wgcna_exact"]


def _event_qc(matrix_with_missing: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    observed_fraction = matrix_with_missing.notna().mean(axis=1)
    variability = matrix_with_missing.quantile(0.95, axis=1) - matrix_with_missing.quantile(0.05, axis=1)
    variance = matrix_with_missing.var(axis=1, skipna=True)
    score = observed_fraction * variability * np.sqrt(variance.clip(lower=0.0))
    qc = pd.DataFrame(
        {
            "event_key": matrix_with_missing.index,
            "observed_fraction": observed_fraction,
            "robust_psi_range": variability,
            "psi_variance": variance,
            "stability_score": score,
            "external_selection_note": (
                "External PSI-only benchmark score; read-support Stable-AS QC was unavailable. "
                "Identical event universe used by every method."
            ),
        }
    )
    matrix = matrix_with_missing.T.fillna(matrix_with_missing.median(axis=1)).T
    matrix = matrix.loc[np.isfinite(matrix.to_numpy(dtype=float)).all(axis=1)]
    qc = qc.loc[qc["event_key"].isin(matrix.index)].copy()
    return matrix, qc


def load_gse141633(path: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, object]]:
    frame = pd.read_csv(path, sep="\t", low_memory=False)
    sample_columns = [column for column in frame.columns[8:] if str(column).startswith("TCGA-")]
    values = frame[sample_columns].apply(pd.to_numeric, errors="coerce")
    event_key = (
        frame["AC"].astype(str)
        + "|"
        + frame["chr"].astype(str)
        + ":"
        + frame["exonStart"].astype(str)
        + "-"
        + frame["exonEnd"].astype(str)
        + ":"
        + frame["strand"].astype(str)
    )
    values.index = event_key
    values = values.loc[~values.index.duplicated(keep="first")]
    values = values.loc[values.notna().mean(axis=1).ge(0.80)]
    samples = []
    for sample in sample_columns:
        fields = sample.split("-")
        sample_type = fields[3][:2] if len(fields) > 3 else ""
        condition = "normal" if sample_type.isdigit() and 10 <= int(sample_type) <= 19 else "tumor"
        samples.append({"sample_id": sample, "condition": condition, "source": "GSE141633_TCGA_PRAD"})
    metadata = pd.DataFrame(samples)
    matrix, qc = _event_qc(values)
    return matrix, metadata, qc, {
        "dataset": "GSE141633_TCGA_PRAD",
        "data_type": "bulk_RNA_seq_PSI",
        "n_events_before_qc": len(frame),
        "n_events_after_observation_filter": len(matrix),
        "n_samples": matrix.shape[1],
        "missing_fraction_before_imputation": float(values.isna().mean().mean()),
        "observation_filter": "TCGA samples; event observed in >=80% of selected samples",
    }


def load_gse85908(
    splicing_path: Path,
    metadata_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, object]]:
    sample_by_event = pd.read_csv(splicing_path, index_col="sample_id", low_memory=False)
    metadata = pd.read_csv(metadata_path, index_col=0, low_memory=False)
    common_samples = sample_by_event.index.intersection(metadata.index)
    sample_by_event = sample_by_event.loc[common_samples]
    metadata = metadata.loc[common_samples]
    matrix_with_missing = sample_by_event.T
    matrix_with_missing = matrix_with_missing.loc[matrix_with_missing.notna().mean(axis=1).ge(0.50)]
    matrix_with_missing = matrix_with_missing.loc[~matrix_with_missing.index.duplicated(keep="first")]
    matrix, qc = _event_qc(matrix_with_missing)
    samples = pd.DataFrame(
        {
            "sample_id": common_samples,
            "condition": metadata["phenotype"].astype(str).to_numpy(),
            "source": "GSE85908_single_cell_neuron_differentiation",
        }
    )
    return matrix, samples, qc, {
        "dataset": "GSE85908",
        "data_type": "single_cell_RNA_seq_PSI",
        "n_events_before_qc": sample_by_event.shape[1],
        "n_events_after_observation_filter": len(matrix),
        "n_samples": matrix.shape[1],
        "missing_fraction_before_imputation": float(matrix_with_missing.isna().mean().mean()),
        "observation_filter": "event observed in >=50% of cells; event-wise median imputation",
    }


def simulate_known_modules(seed: int = 20260716) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, object]]:
    rng = np.random.default_rng(seed)
    n_samples = 60
    modules = 3
    events_per_module = 50
    conditions = np.repeat(["control", "treatment"], n_samples // 2)
    latent = rng.normal(size=(modules, n_samples))
    latent[0, conditions == "treatment"] += 0.8
    latent[1, conditions == "treatment"] -= 0.7
    values: list[np.ndarray] = []
    truth: list[str] = []
    event_keys: list[str] = []
    for module_index in range(modules):
        for event_index in range(events_per_module):
            loading = rng.uniform(0.65, 0.95)
            signal = loading * latent[module_index] + rng.normal(scale=0.35, size=n_samples)
            psi = 1.0 / (1.0 + np.exp(-signal))
            values.append(psi)
            truth.append(f"M{module_index + 1:02d}")
            event_keys.append(f"sim_M{module_index + 1:02d}_event_{event_index + 1:03d}")
    matrix = pd.DataFrame(values, index=event_keys, columns=[f"sim_sample_{index + 1:03d}" for index in range(n_samples)])
    qc = pd.DataFrame({"event_key": event_keys, "stability_score": matrix.var(axis=1).to_numpy()})
    samples = pd.DataFrame({"sample_id": matrix.columns, "condition": conditions, "source": "known_truth_simulation"})
    truth_frame = pd.DataFrame({"event_key": event_keys, "truth_module": truth})
    return matrix, samples, qc, {
        "dataset": "known_truth_simulation",
        "data_type": "simulated_PSI_with_three_known_modules",
        "n_events_before_qc": len(matrix),
        "n_events_after_observation_filter": len(matrix),
        "n_samples": n_samples,
        "missing_fraction_before_imputation": 0.0,
        "observation_filter": "none; all simulated events retained",
        "truth": truth_frame,
    }


def write_result(
    dataset: str,
    output_root: Path,
    matrix: pd.DataFrame,
    samples: pd.DataFrame,
    qc: pd.DataFrame,
    provenance: dict[str, object],
    network: NetworkConfig,
    benchmark: BenchmarkConfig,
) -> pd.DataFrame:
    output_dir = output_root / dataset
    output_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    result = run_method_benchmark(matrix, qc, samples, network, benchmark, progress=print)
    elapsed = time.perf_counter() - started
    result.universe.to_csv(output_dir / "common_event_universe.tsv", sep="\t", index=False)
    result.assignments.to_csv(output_dir / "event_method_assignments.tsv", sep="\t", index=False)
    result.method_summary.to_csv(output_dir / "method_benchmark_summary.tsv", sep="\t", index=False)
    result.bootstrap_stability.to_csv(output_dir / "bootstrap_partition_stability.tsv", sep="\t", index=False)
    result.module_matching.to_csv(output_dir / "module_matching.tsv", sep="\t", index=False)
    result.method_audit.to_csv(output_dir / "method_availability.tsv", sep="\t", index=False)
    result.metric_ranking.to_csv(output_dir / "method_metric_ranking.tsv", sep="\t", index=False)
    result.method_interpretation.to_csv(output_dir / "method_strengths_limitations.tsv", sep="\t", index=False)
    qc.to_csv(output_dir / "external_event_qc.tsv", sep="\t", index=False)
    samples.to_csv(output_dir / "sample_metadata.tsv", sep="\t", index=False)
    render_benchmark_plots(result, output_dir / "plots")
    render_benchmark_report(result, output_dir / "method_benchmark_report.html")

    truth = provenance.pop("truth", None)
    truth_metrics: list[dict[str, object]] = []
    if isinstance(truth, pd.DataFrame):
        joined = result.assignments.merge(truth, on="event_key", how="left")
        for method in METHODS:
            column = f"module_{method}"
            if column in joined:
                truth_metrics.append(
                    {
                        "dataset": dataset,
                        "method": method,
                        "truth_adjusted_rand_index": adjusted_rand_score(joined["truth_module"], joined[column]),
                    }
                )
        pd.DataFrame(truth_metrics).to_csv(output_dir / "known_truth_recovery.tsv", sep="\t", index=False)
    (output_dir / "parameters.json").write_text(
        json.dumps(
            {
                **provenance,
                "network": network.__dict__,
                "benchmark": benchmark.__dict__,
                "wall_clock_seconds": elapsed,
                "selection_scope": "same preselected event universe and samples for every method",
                "claim_limit": "Computational generalization benchmark; not biological validation of inferred modules.",
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    summary = result.method_summary.copy()
    summary.insert(0, "dataset", dataset)
    if truth_metrics:
        summary = summary.merge(pd.DataFrame(truth_metrics), on=["dataset", "method"], how="left")
    return summary


def render_aggregate(summary: pd.DataFrame, output_root: Path) -> None:
    summary.to_csv(output_root / "aggregate_method_benchmark_summary.tsv", sep="\t", index=False)
    primary = summary.loc[summary["method"].isin(["katana_sparse", "wgcna_exact"]) & summary["status"].eq("finished")].copy()
    metrics = [
        "silhouette_cosine",
        "neutral_graph_modularity",
        "within_between_separation",
        "mean_bootstrap_ARI",
    ]
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    for ax, metric in zip(axes.flat, metrics, strict=True):
        sns.barplot(data=primary, x="dataset", y=metric, hue="method", ax=ax)
        ax.set_title(metric.replace("_", " ").title(), fontweight="bold")
        ax.tick_params(axis="x", rotation=20)
        ax.legend(fontsize=8)
    fig.suptitle("External CoSplice benchmark: identical event universe per dataset", fontweight="bold", fontsize=16)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(output_root / "external_benchmark_scorecard.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    win_rows: list[dict[str, object]] = []
    lower_is_better = {"runtime_seconds"}
    compare_metrics = [*metrics, "runtime_seconds"]
    if "truth_adjusted_rand_index" in primary.columns:
        compare_metrics.append("truth_adjusted_rand_index")
    for dataset, frame in primary.groupby("dataset", observed=True):
        indexed = frame.set_index("method")
        if not {"katana_sparse", "wgcna_exact"}.issubset(indexed.index):
            continue
        for metric in compare_metrics:
            left = indexed.at["katana_sparse", metric]
            right = indexed.at["wgcna_exact", metric]
            if pd.isna(left) or pd.isna(right):
                continue
            if np.isclose(left, right):
                winner = "tie"
            elif metric in lower_is_better:
                winner = "katana_sparse" if left < right else "wgcna_exact"
            else:
                winner = "katana_sparse" if left > right else "wgcna_exact"
            win_rows.append({"dataset": dataset, "metric": metric, "katana": left, "wgcna": right, "winner": winner})
    pd.DataFrame(win_rows).to_csv(output_root / "katana_vs_wgcna_metric_winners.tsv", sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run reproducible external Katana-CoSplice method benchmarks.")
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--extended-events", type=int, default=1000)
    args = parser.parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)
    network = NetworkConfig(
        min_module_size=20,
        knn_neighbors=50,
        community_resolution=1.0,
        module_refinement=True,
        min_module_membership=0.45,
        module_merge_correlation=0.85,
        bootstrap_iterations=0,
    )
    benchmark = BenchmarkConfig(methods=METHODS, max_events=150, bootstrap_iterations=5, random_seed=41)

    datasets: list[tuple[str, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, object]]] = []
    matrix, samples, qc, provenance = load_gse141633(args.input_root / "GSE141633_PSI.txt.gz")
    datasets.append(("GSE141633_TCGA_PRAD", matrix, samples, qc, provenance))
    matrix, samples, qc, provenance = load_gse85908(
        args.input_root / "GSE85908_splicing.csv.gz",
        args.input_root / "GSE85908_metadata.csv.gz",
    )
    datasets.append(("GSE85908", matrix, samples, qc, provenance))
    matrix, samples, qc, provenance = simulate_known_modules()
    datasets.append(("known_truth_simulation", matrix, samples, qc, provenance))

    summaries = []
    for dataset, matrix, samples, qc, provenance in datasets:
        print(f"\n=== {dataset}: {matrix.shape[0]} events x {matrix.shape[1]} samples ===")
        summaries.append(write_result(dataset, args.output_root, matrix, samples, qc, dict(provenance), network, benchmark))
        if dataset != "known_truth_simulation" and len(matrix) >= args.extended_events:
            extended = BenchmarkConfig(
                methods=["katana_sparse", "wgcna_exact"],
                max_events=args.extended_events,
                bootstrap_iterations=5,
                random_seed=41,
            )
            scaled_name = f"{dataset}_extended_{args.extended_events}"
            print(f"\n=== {scaled_name}: scale-sensitivity benchmark ===")
            summaries.append(
                write_result(
                    scaled_name,
                    args.output_root,
                    matrix,
                    samples,
                    qc,
                    {**provenance, "scale_sensitivity_test": True},
                    network,
                    extended,
                )
            )
    render_aggregate(pd.concat(summaries, ignore_index=True), args.output_root)


if __name__ == "__main__":
    main()
