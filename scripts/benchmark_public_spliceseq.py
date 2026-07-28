from __future__ import annotations

import argparse
import json
import time
from dataclasses import asdict, replace
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.metrics import adjusted_rand_score

from katana_cosplice.benchmark import (
    _fit_method,
    _method_metrics,
    _neutral_graph,
    _rscript_with_wgcna,
)
from katana_cosplice.config import (
    BenchmarkConfig,
    ConsensusConfig,
    ExternalCohortConfig,
    NetworkConfig,
    PreservationConfig,
)
from katana_cosplice.consensus import (
    _bootstrap_indices,
    build_bootstrap_coassignment_consensus,
    render_consensus_plots,
)
from katana_cosplice.network import discover_module_labels
from katana_cosplice.preservation import (
    evaluate_external_module_preservation,
    render_preservation_plots,
)


def _coassignment(labels: list[np.ndarray]) -> np.ndarray:
    count = np.zeros((len(labels[0]), len(labels[0])), dtype=np.float32)
    for run in labels:
        for module in set(run).difference({"grey"}):
            indices = np.flatnonzero(run == module)
            count[np.ix_(indices, indices)] += 1.0
    count /= len(labels)
    np.fill_diagonal(count, 1.0)
    return count


def _within_coassignment(labels: np.ndarray, matrix: np.ndarray) -> float:
    values: list[float] = []
    weights: list[int] = []
    for module in set(labels).difference({"grey"}):
        indices = np.flatnonzero(labels == module)
        if len(indices) < 2:
            continue
        upper = matrix[np.ix_(indices, indices)][np.triu_indices(len(indices), k=1)]
        values.append(float(upper.mean()))
        weights.append(len(upper))
    return float(np.average(values, weights=weights)) if values else 0.0


def _samples(columns: pd.Index, cohort: str) -> pd.DataFrame:
    return pd.DataFrame(
        {"sample_id": columns.astype(str), "condition": cohort, "cohort": cohort}
    )


def _assignments(event_keys: pd.Index, labels: np.ndarray, column: str) -> pd.DataFrame:
    return pd.DataFrame({"event_key": event_keys.astype(str), column: labels.astype(str)})


def _aggregate_preservation(method: str, result) -> dict[str, object]:
    table = result.module_preservation
    audit = result.target_audit.iloc[0].to_dict() if not result.target_audit.empty else {}
    if table.empty:
        return {
            "method": method,
            "status": audit.get("status", "no_preservable_modules"),
            "mapped_events": int(result.event_mapping["event_key"].nunique())
            if not result.event_mapping.empty
            else 0,
            "tested_modules": 0,
        }
    counts = table["preservation_class"].value_counts()
    return {
        "method": method,
        "status": audit.get("status", "finished"),
        "mapped_events": int(result.event_mapping["event_key"].nunique()),
        "tested_modules": int(len(table)),
        "preserved_modules": int(counts.get("preserved", 0)),
        "partially_preserved_modules": int(counts.get("partially_preserved", 0)),
        "not_preserved_modules": int(counts.get("not_preserved", 0)),
        "median_Zsummary": float(table["Zsummary_wgcna_like"].median()),
        "median_density_ratio": float(table["density_ratio"].median()),
        "median_connectivity_correlation": float(table["connectivity_correlation"].median()),
        "median_membership_profile_correlation": float(
            table["membership_profile_correlation"].median()
        ),
        "median_program_score_preservation": float(table["program_score_preservation"].median()),
    }


def _dashboard(
    internal: pd.DataFrame,
    external: pd.DataFrame,
    output: Path,
    *,
    discovery_name: str,
    external_name: str,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    internal_plot = internal.melt(
        id_vars="method",
        value_vars=["mean_bootstrap_ARI", "min_bootstrap_ARI"],
        var_name="metric",
        value_name="value",
    )
    sns.barplot(data=internal_plot, x="metric", y="value", hue="method", ax=axes[0])
    axes[0].set_ylim(0, 1)
    axes[0].set_title(f"{discovery_name} internal bootstrap stability")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("Adjusted Rand index")

    if not external.empty and "median_Zsummary" in external:
        sns.barplot(data=external, x="method", y="median_Zsummary", ax=axes[1])
        axes[1].axhline(2, color="#d97706", linestyle="--", linewidth=1, label="partial: 2")
        axes[1].axhline(10, color="#15803d", linestyle="--", linewidth=1, label="preserved: 10")
        axes[1].legend(frameon=False)
        axes[1].set_title(f"Independent {external_name} module preservation")
        axes[1].set_ylabel("Median WGCNA-like Zsummary")
        axes[1].tick_params(axis="x", rotation=20)
    else:
        axes[1].axis("off")
        axes[1].text(0.5, 0.5, "No preservable external modules", ha="center", va="center")
    fig.suptitle(
        f"Public co-splicing benchmark: {discovery_name} discovery to {external_name} validation",
        fontsize=14,
        fontweight="bold",
    )
    fig.tight_layout()
    fig.savefig(output, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _write_report(
    output: Path,
    internal: pd.DataFrame,
    external: pd.DataFrame,
    *,
    discovery_name: str,
    external_name: str,
    parameters: dict[str, object],
) -> Path:
    def markdown_table(frame: pd.DataFrame) -> str:
        if frame.empty:
            return "No rows."
        display = frame.copy()
        for column in display.select_dtypes(include=["number"]).columns:
            display[column] = display[column].map(
                lambda value: "" if pd.isna(value) else f"{float(value):.4f}"
            )
        headers = [str(column).replace("|", "\\|") for column in display.columns]
        lines = [
            "| " + " | ".join(headers) + " |",
            "| " + " | ".join(["---"] * len(headers)) + " |",
        ]
        for row in display.itertuples(index=False, name=None):
            lines.append(
                "| "
                + " | ".join(str(value).replace("|", "\\|") for value in row)
                + " |"
            )
        return "\n".join(lines)

    internal_rows = markdown_table(internal)
    external_rows = (
        markdown_table(external)
        if not external.empty
        else "No module passed the minimum event-overlap requirement for external preservation."
    )
    report = f"""# Public CoSplice vs exact WGCNA validation

## Design

- Discovery cohort: `{discovery_name}`
- Independent validation cohort: `{external_name}`
- Discovery events: {parameters['n_discovery_events']}
- Discovery samples: {parameters['n_discovery_samples']}
- Validation samples: {parameters['n_external_samples']}
- Exact shared event keys: {parameters['n_exact_shared_events']}
- Matched bootstrap iterations: {parameters['bootstraps']}
- Preservation permutations per module: {parameters['permutations_per_module']}

Event selection is performed from the discovery cohort only. The external cohort is used only after
module discovery, so it cannot influence the selected event universe or module labels.

## Internal matched-bootstrap benchmark

{internal_rows}

## Independent-cohort preservation

{external_rows}

## Interpretation boundary

This benchmark compares partition stability, module granularity and independent-cohort preservation
on one matched event universe. A higher value for one metric does not establish universal method
superiority. `Zsummary_wgcna_like` is Katana's permutation-calibrated density/connectivity composite,
not the exact `WGCNA::modulePreservation` statistic. Biological claims require replication across
additional independent cohorts and mechanism-specific validation.
"""
    path = output / "PUBLIC_VALIDATION_REPORT.md"
    path.write_text(report, encoding="utf-8")
    return path


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Benchmark CoSplice against exact WGCNA on public SpliceSeq PSI cohorts."
    )
    parser.add_argument("--discovery", type=Path, required=True)
    parser.add_argument("--external", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--discovery-name", default="TCGA-KIRC")
    parser.add_argument("--external-name", default="TCGA-KIRP")
    parser.add_argument("--events", type=int, default=150)
    parser.add_argument("--bootstraps", type=int, default=10)
    parser.add_argument("--permutations", type=int, default=500)
    parser.add_argument("--min-module-size", type=int, default=15)
    args = parser.parse_args()

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    discovery = pd.read_csv(args.discovery, sep="\t", index_col=0).apply(
        pd.to_numeric, errors="coerce"
    )
    external = pd.read_csv(args.external, sep="\t", index_col=0).apply(
        pd.to_numeric, errors="coerce"
    )
    discovery = discovery.loc[
        discovery.notna().mean(axis=1).ge(0.75) & discovery.std(axis=1).gt(1e-8)
    ]
    if len(discovery) > args.events:
        discovery = discovery.loc[discovery.var(axis=1).nlargest(args.events).index]
    discovery = discovery.apply(lambda row: row.fillna(row.mean()), axis=1)
    discovery.index = discovery.index.astype(str)
    external.index = external.index.astype(str)
    discovery_samples = _samples(discovery.columns, args.discovery_name)
    external_samples = _samples(external.columns, args.external_name)
    external_samples.to_csv(output / "external_sample_metadata.tsv", sep="\t", index=False)

    network = NetworkConfig(
        min_module_size=args.min_module_size,
        knn_neighbors=min(30, max(5, len(discovery) // 5)),
        consensus_knn_neighbors=[15, 25, 35],
        consensus_bootstrap_iterations=0,
        community_consensus_runs=1,
        bootstrap_iterations=0,
        edge_threshold=0.05,
        random_seed=17,
    )
    consensus_settings = ConsensusConfig(
        bootstrap_iterations=args.bootstraps,
        min_module_size=args.min_module_size,
        graph_top_k=min(35, max(5, len(discovery) - 1)),
        dense_max_events=max(200, len(discovery)),
        random_seed=37,
    )
    benchmark = BenchmarkConfig(
        methods=["wgcna_exact"],
        max_events=len(discovery),
        bootstrap_iterations=args.bootstraps,
        random_seed=41,
    )
    preservation = PreservationConfig(
        permutations=args.permutations,
        min_shared_events=max(5, min(10, args.min_module_size)),
        dense_max_events=max(200, len(discovery)),
        random_seed=43,
    )

    rscript, wgcna_audit = _rscript_with_wgcna()
    if rscript is None:
        raise RuntimeError(wgcna_audit)

    partition_network = replace(
        network,
        consensus_bootstrap_iterations=0,
        community_consensus_runs=1,
        bootstrap_iterations=0,
    )
    observed_labels = discover_module_labels(discovery, partition_network)
    observed = _assignments(discovery.index, observed_labels, "module")

    start = time.perf_counter()
    cosplice = build_bootstrap_coassignment_consensus(
        discovery,
        observed,
        discovery_samples,
        network,
        consensus_settings,
        progress=print,
    )
    cosplice_seconds = time.perf_counter() - start
    cosplice_labels = (
        cosplice.assignments.set_index("event_key")
        .reindex(discovery.index)["consensus_module"]
        .to_numpy(object)
    )

    start = time.perf_counter()
    wgcna_labels, _ = _fit_method("wgcna_exact", discovery, network, benchmark, rscript)
    wgcna_seconds = time.perf_counter() - start

    rng = np.random.default_rng(59)
    cosplice_runs: list[np.ndarray] = []
    wgcna_runs: list[np.ndarray] = []
    stability_rows: list[dict[str, object]] = []
    for iteration in range(1, args.bootstraps + 1):
        indices = _bootstrap_indices(discovery_samples, discovery.columns, rng)
        bootstrap = discovery.iloc[:, indices].copy()
        bootstrap.columns = [f"bootstrap_{iteration:03d}_{i:03d}" for i in range(len(indices))]
        c_labels = discover_module_labels(bootstrap, partition_network)
        w_labels, _ = _fit_method("wgcna_exact", bootstrap, network, benchmark, rscript)
        cosplice_runs.append(c_labels)
        wgcna_runs.append(w_labels)
        stability_rows.append(
            {
                "iteration": iteration,
                "cosplice_ARI": adjusted_rand_score(cosplice_labels, c_labels),
                "wgcna_ARI": adjusted_rand_score(wgcna_labels, w_labels),
            }
        )
        print(f"Matched bootstrap {iteration}/{args.bootstraps}")
    stability = pd.DataFrame(stability_rows)
    neutral = _neutral_graph(discovery, benchmark.neutral_knn_neighbors)
    internal_rows: list[dict[str, object]] = []
    for method, labels, coassign, seconds, column in (
        (
            "CoSplice",
            cosplice_labels,
            cosplice.coassignment_matrix,
            cosplice_seconds,
            "cosplice_ARI",
        ),
        (
            "Exact WGCNA",
            wgcna_labels,
            _coassignment(wgcna_runs),
            wgcna_seconds,
            "wgcna_ARI",
        ),
    ):
        row: dict[str, object] = {
            "method": method,
            "n_events": len(discovery),
            "n_samples": discovery.shape[1],
            "runtime_seconds_initial_fit": seconds,
            "mean_bootstrap_ARI": float(stability[column].mean()),
            "median_bootstrap_ARI": float(stability[column].median()),
            "min_bootstrap_ARI": float(stability[column].min()),
            "weighted_mean_within_module_coassignment": _within_coassignment(labels, coassign),
        }
        row.update(_method_metrics(discovery, labels, discovery_samples, neutral))
        internal_rows.append(row)
    internal = pd.DataFrame(internal_rows)

    cohort = ExternalCohortConfig(
        args.external_name,
        args.external.resolve(),
        output / "external_sample_metadata.tsv",
    )
    preservation_rows: list[pd.DataFrame] = []
    external_summary_rows: list[dict[str, object]] = []
    for method, labels in (("CoSplice", cosplice_labels), ("Exact WGCNA", wgcna_labels)):
        method_dir = output / method.lower().replace(" ", "_")
        method_dir.mkdir(parents=True, exist_ok=True)
        assignment = _assignments(discovery.index, labels, "consensus_module")
        result = evaluate_external_module_preservation(
            discovery,
            assignment,
            [cohort],
            preservation,
            progress=print,
        )
        detail = result.module_preservation.copy()
        detail.insert(0, "method", method)
        preservation_rows.append(detail)
        result.module_preservation.to_csv(
            method_dir / "external_module_preservation.tsv", sep="\t", index=False
        )
        result.event_mapping.to_csv(method_dir / "external_event_mapping.tsv", sep="\t", index=False)
        result.target_audit.to_csv(method_dir / "external_target_audit.tsv", sep="\t", index=False)
        render_preservation_plots(result, method_dir / "plots")
        external_summary_rows.append(_aggregate_preservation(method, result))

    external_detail = pd.concat(preservation_rows, ignore_index=True)
    external_summary = pd.DataFrame(external_summary_rows)
    internal.to_csv(output / "internal_matched_bootstrap_benchmark.tsv", sep="\t", index=False)
    stability.to_csv(output / "matched_bootstrap_stability.tsv", sep="\t", index=False)
    external_detail.to_csv(output / "external_module_preservation_head_to_head.tsv", sep="\t", index=False)
    external_summary.to_csv(output / "external_preservation_summary.tsv", sep="\t", index=False)
    cosplice.assignments.to_csv(output / "cosplice_assignments.tsv", sep="\t", index=False)
    _assignments(discovery.index, wgcna_labels, "wgcna_module").to_csv(
        output / "wgcna_assignments.tsv", sep="\t", index=False
    )
    cosplice.parameter_search.to_csv(output / "cosplice_parameter_search.tsv", sep="\t", index=False)
    render_consensus_plots(cosplice, output / "cosplice_plots")
    parameters = {
        "discovery": str(args.discovery.resolve()),
        "external": str(args.external.resolve()),
        "discovery_name": args.discovery_name,
        "external_name": args.external_name,
        "selection_scope": "discovery_only",
        "event_matching": "exact event key",
        "n_discovery_events": len(discovery),
        "n_discovery_samples": discovery.shape[1],
        "n_external_samples": external.shape[1],
        "n_exact_shared_events": int(discovery.index.intersection(external.index).size),
        "bootstraps": args.bootstraps,
        "permutations_per_module": args.permutations,
        "network": asdict(network),
        "consensus": asdict(consensus_settings),
        "preservation": asdict(preservation),
        "exact_wgcna": wgcna_audit,
    }
    (output / "parameters.json").write_text(json.dumps(parameters, indent=2), encoding="utf-8")
    _dashboard(
        internal,
        external_summary,
        output / "public_external_cosplice_vs_wgcna.png",
        discovery_name=args.discovery_name,
        external_name=args.external_name,
    )
    _write_report(
        output,
        internal,
        external_summary,
        discovery_name=args.discovery_name,
        external_name=args.external_name,
        parameters=parameters,
    )
    print("\nInternal benchmark")
    print(internal.to_string(index=False))
    print("\nExternal preservation")
    print(external_summary.to_string(index=False))
    print(f"\nOutput: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
