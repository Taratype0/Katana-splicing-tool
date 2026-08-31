from __future__ import annotations

import argparse
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

from katana_cosplice.benchmark import _fit_method, _rscript_with_wgcna
from katana_cosplice.config import (
    BenchmarkConfig,
    ConsensusConfig,
    ExternalCohortConfig,
    MultiObjectiveSelectionConfig,
    NetworkConfig,
    PreservationConfig,
)
from katana_cosplice.consensus import _bootstrap_indices, build_bootstrap_coassignment_consensus
from katana_cosplice.locked_ablation import (
    combine_gmt_files,
    file_record,
    functional_replication,
    load_gmt,
    local_module_enrichment,
    read_psi_matrix,
    sample_frame,
    select_discovery_events,
    sha256_file,
    utc_now,
    write_json,
)
from katana_cosplice.multiobjective_selection import build_multiobjective_candidate_scorer
from katana_cosplice.network import discover_module_labels
from katana_cosplice.preservation import evaluate_external_module_preservation


METHODS = ["stability_only", "three_objective_cosplice", "event_level_wgcna"]


def _cohort_paths(
    panel: Path, cohort: str, event_type: str = "ES"
) -> dict[str, Path]:
    prepared = panel / cohort / "prepared"
    event_type = event_type.upper()
    return {
        "psi": prepared / f"{cohort}_{event_type}_PSI_tumor_only.tsv",
        "events": prepared / f"{cohort}_{event_type}_event_metadata.tsv",
        "samples": prepared / f"{cohort}_sample_metadata.tsv",
    }


def _assignment(index: pd.Index, labels: np.ndarray, column: str = "consensus_module") -> pd.DataFrame:
    return pd.DataFrame({"event_key": index.astype(str), column: labels.astype(str)})


def _labels(frame: pd.DataFrame) -> np.ndarray:
    column = next(name for name in ("consensus_module", "wgcna_module", "module") if name in frame)
    return frame[column].astype(str).to_numpy(object)


def _selected_parameters(result) -> tuple[float, float]:
    selected = result.parameter_search.loc[result.parameter_search["selected"]].iloc[0]
    return float(selected["edge_threshold"]), float(selected["resolution"])


def _wgcna_stability(
    matrix: pd.DataFrame,
    samples: pd.DataFrame,
    labels: np.ndarray,
    network: NetworkConfig,
    benchmark: BenchmarkConfig,
    rscript: str,
    bootstraps: int,
) -> dict[str, float]:
    rng = np.random.default_rng(97)
    scores: list[float] = []
    for iteration in range(1, bootstraps + 1):
        indices = _bootstrap_indices(samples, matrix.columns, rng)
        current = matrix.iloc[:, indices].copy()
        current.columns = [f"wgcna_{iteration:03d}_{index:03d}" for index in range(len(indices))]
        current_labels, _ = _fit_method("wgcna_exact", current, network, benchmark, rscript)
        scores.append(float(adjusted_rand_score(labels, current_labels)))
        print(f"WGCNA matched bootstrap {iteration}/{bootstraps}", flush=True)
    values = np.asarray(scores, dtype=float)
    return {
        "mean_bootstrap_ARI": float(values.mean()),
        "median_bootstrap_ARI": float(np.median(values)),
        "minimum_bootstrap_ARI": float(values.min()),
    }


def _method_internal_row(method: str, assignments: pd.DataFrame, result=None, extra=None) -> dict[str, object]:
    labels = _labels(assignments)
    counts = pd.Series(labels[labels != "grey"]).value_counts()
    row: dict[str, object] = {
        "method": method,
        "n_modules": len(counts),
        "assigned_event_fraction": float(np.mean(labels != "grey")),
        "largest_module_fraction": float(counts.max() / len(labels)) if not counts.empty else 0.0,
    }
    if result is not None:
        selected = result.parameter_search.loc[result.parameter_search["selected"]].iloc[0]
        row.update(
            {
                "mean_bootstrap_ARI": float(selected["mean_bootstrap_ari"]),
                "median_bootstrap_ARI": float(selected["median_bootstrap_ari"]),
                "minimum_bootstrap_ARI": float(selected["min_bootstrap_ari"]),
                "robust_bootstrap_stability": float(selected["robust_bootstrap_stability"]),
            }
        )
    if extra:
        row.update(extra)
    return row


def _render_dashboard(
    summary: pd.DataFrame,
    output: Path,
    *,
    title: str = "Frozen three-arm CoSplice ablation",
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))
    metrics = [
        ("mean_bootstrap_ARI", "Internal bootstrap stability", "Adjusted Rand index"),
        ("median_holdout_Zsummary", "Locked LUAD preservation", "Median Zsummary-like"),
        ("median_best_term_jaccard", "Independent functional replication", "Term Jaccard"),
    ]
    palette = ["#2563EB", "#0F766E", "#D97706"]
    for axis, (column, title, ylabel), color in zip(axes, metrics, palette, strict=True):
        sns.barplot(data=summary, x="method", y=column, color=color, ax=axis)
        axis.set_title(title, fontweight="bold")
        axis.set_xlabel("")
        axis.set_ylabel(ylabel)
        axis.tick_params(axis="x", rotation=25)
        for patch, value in zip(axis.patches, summary[column], strict=False):
            axis.text(patch.get_x() + patch.get_width() / 2, patch.get_height(), f"{value:.3f}", ha="center", va="bottom", fontsize=8)
        sns.despine(ax=axis)
    fig.suptitle(title, fontsize=15, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _write_report(
    output: Path,
    summary: pd.DataFrame,
    design: dict[str, object],
) -> Path:
    display = summary.copy()
    for column in display.select_dtypes(include=["number"]).columns:
        display[column] = display[column].map(
            lambda value: "" if pd.isna(value) else f"{float(value):.4f}"
        )
    headers = [str(column) for column in display.columns]
    rows = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    rows.extend(
        "| " + " | ".join(map(str, row)) + " |"
        for row in display.itertuples(index=False, name=None)
    )
    functional_active = False
    search_path = (
        output
        / "01_prelock"
        / "three_objective_cosplice"
        / "parameter_search.tsv"
    )
    if search_path.exists():
        search = pd.read_csv(search_path, sep="\t")
        selected = search.loc[search["selected"].astype(bool)]
        if not selected.empty and "functional_objective_active" in selected:
            functional_active = bool(selected.iloc[0]["functional_objective_active"])
    report = output / "LOCKED_ABLATION_REPORT.md"
    report.write_text(
        "\n".join(
            [
                "# Frozen three-arm CoSplice ablation",
                "",
                "## Registered design",
                "",
                f"- Discovery: `{design['discovery']}`",
                f"- Tuning cohorts: `{', '.join(design['tuning_cohorts'])}`",
                f"- Holdout: `{design['holdout']}`",
                f"- Event type: `{design['event_type']}`",
                "- Direct arms: internal-quality CoSplice, adaptive multi-objective CoSplice, exact R WGCNA.",
                "- The holdout was opened only after discovery assignments and parameters were fingerprinted.",
                "",
                "## Results",
                "",
                *rows,
                "",
                "## Functional-objective gate",
                "",
                (
                    "The functional objective was active because the predeclared tuning GMT reached the formal enrichment gate."
                    if functional_active
                    else "The functional objective was automatically disabled because no candidate reached the predeclared enrichment FDR gate; its weight was redistributed over internal quality and tuning-cohort preservation."
                ),
                "",
                "## Interpretation boundary",
                "",
                "Internal ARI is interpreted only together with module count, assigned-event fraction and largest-module fraction. A coarse single-module partition is not treated as biologically superior solely because it is stable. Holdout preservation and functional replication are confirmatory endpoints; absent replication is reported as a negative result, not omitted.",
                "",
                "## Reproducibility",
                "",
                "See `analysis_plan.json`, `freeze_manifest.json`, `holdout_access_audit.json`, method-specific assignments and parameter-search tables in this directory.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Leakage-safe three-arm CoSplice ablation with a locked holdout.")
    parser.add_argument("--panel-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--discovery", default="KIRC")
    parser.add_argument("--tuning", nargs="+", default=["KIRP", "KICH"])
    parser.add_argument("--holdout", default="LUAD")
    parser.add_argument("--event-type", default="ES")
    parser.add_argument("--gmt", type=Path, action="append", required=True)
    parser.add_argument(
        "--rscript",
        type=Path,
        help="Explicit Rscript executable from an environment containing WGCNA.",
    )
    parser.add_argument("--events", type=int, default=500)
    parser.add_argument("--bootstraps", type=int, default=20)
    parser.add_argument("--permutations", type=int, default=500)
    parser.add_argument("--min-module-size", type=int, default=15)
    args = parser.parse_args(argv)

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    prelock = output / "01_prelock"
    postlock = output / "02_holdout_evaluation"
    prelock.mkdir(exist_ok=True)
    postlock.mkdir(exist_ok=True)
    event_type = args.event_type.upper()
    paths = {
        name: _cohort_paths(args.panel_root.resolve(), name, event_type)
        for name in [args.discovery, *args.tuning, args.holdout]
    }
    design = {
        "registered_at_utc": utc_now(),
        "design": "frozen three-arm ablation",
        "discovery": args.discovery,
        "tuning_cohorts": args.tuning,
        "holdout": args.holdout,
        "event_type": event_type,
        "holdout_boundary": "No holdout file is opened, fingerprinted, or summarized before freeze_manifest.json is written.",
        "historical_use_caveat": "LUAD is untouched by this ablation runner before freezing, but was used in earlier project-level validation and is not globally novel to the project.",
        "arms": METHODS,
        "primary_endpoints": ["holdout module preservation", "independent GO/KEGG term replication"],
        "holdout_planned_paths": {key: str(value.resolve()) for key, value in paths[args.holdout].items()},
    }
    write_json(output / "analysis_plan.json", design)

    merged_gmt = combine_gmt_files(args.gmt, prelock / "frozen_GO_KEGG.gmt")
    discovery_path = paths[args.discovery]
    discovery_raw = read_psi_matrix(discovery_path["psi"])
    discovery = select_discovery_events(discovery_raw, args.events)
    samples = sample_frame(discovery_path["samples"], discovery.columns, args.discovery)
    event_metadata = pd.read_csv(discovery_path["events"], sep="\t")
    event_metadata = event_metadata.loc[event_metadata["event_key"].astype(str).isin(discovery.index)]
    event_metadata.to_csv(prelock / "frozen_event_metadata.tsv", sep="\t", index=False)
    discovery.to_csv(prelock / "frozen_discovery_matrix.tsv", sep="\t")

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
    consensus = ConsensusConfig(
        bootstrap_iterations=args.bootstraps,
        min_module_size=args.min_module_size,
        graph_top_k=min(50, max(5, len(discovery) - 1)),
        dense_max_events=max(500, len(discovery)),
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
        dense_max_events=max(500, len(discovery)),
        random_seed=43,
    )
    partition_network = replace(network, consensus_bootstrap_iterations=0, community_consensus_runs=1, bootstrap_iterations=0)
    observed = _assignment(discovery.index, discover_module_labels(discovery, partition_network), "module")

    print("Fitting stability-only CoSplice arm", flush=True)
    stability_result = build_bootstrap_coassignment_consensus(discovery, observed, samples, network, consensus, progress=print)
    stability_assignment = stability_result.assignments[["event_key", "consensus_module"]].copy()

    tuning_cohorts = [ExternalCohortConfig(name, paths[name]["psi"], paths[name]["samples"]) for name in args.tuning]
    objective = MultiObjectiveSelectionConfig(
        enabled=True,
        tuning_cohorts=list(args.tuning),
        holdout_cohorts=[args.holdout],
        functional_gene_sets_gmt=merged_gmt,
        stability_weight=0.40,
        preservation_weight=0.35,
        functional_weight=0.25,
    )
    scorer = build_multiobjective_candidate_scorer(discovery, event_metadata, tuning_cohorts, objective)
    print("Fitting three-objective CoSplice arm", flush=True)
    multi_result = build_bootstrap_coassignment_consensus(discovery, observed, samples, network, consensus, candidate_scorer=scorer, progress=print)
    multi_assignment = multi_result.assignments[["event_key", "consensus_module"]].copy()

    if args.rscript is not None:
        if not args.rscript.exists():
            raise FileNotFoundError(f"Explicit Rscript does not exist: {args.rscript}")
        rscript = str(args.rscript.resolve())
        wgcna_audit = f"Explicit Rscript with WGCNA validation deferred to fit: {rscript}"
    else:
        rscript, wgcna_audit = _rscript_with_wgcna()
        if rscript is None:
            raise RuntimeError(wgcna_audit)
    print("Fitting exact event-level WGCNA arm", flush=True)
    start = time.perf_counter()
    wgcna_labels, _ = _fit_method("wgcna_exact", discovery, network, benchmark, rscript)
    wgcna_runtime = time.perf_counter() - start
    wgcna_assignment = _assignment(discovery.index, wgcna_labels, "wgcna_module")
    wgcna_stability = _wgcna_stability(discovery, samples, wgcna_labels, network, benchmark, rscript, args.bootstraps)

    assignments = {
        "stability_only": stability_assignment,
        "three_objective_cosplice": multi_assignment,
        "event_level_wgcna": wgcna_assignment,
    }
    for method, frame in assignments.items():
        method_dir = prelock / method
        method_dir.mkdir(exist_ok=True)
        frame.to_csv(method_dir / "discovery_assignments.tsv", sep="\t", index=False)
    stability_result.parameter_search.to_csv(prelock / "stability_only" / "parameter_search.tsv", sep="\t", index=False)
    multi_result.parameter_search.to_csv(prelock / "three_objective_cosplice" / "parameter_search.tsv", sep="\t", index=False)
    multi_result.selection_audit.to_csv(prelock / "three_objective_cosplice" / "selection_audit.tsv", sep="\t", index=False)

    internal = pd.DataFrame(
        [
            _method_internal_row("stability_only", stability_assignment, stability_result),
            _method_internal_row("three_objective_cosplice", multi_assignment, multi_result),
            _method_internal_row("event_level_wgcna", wgcna_assignment, extra={**wgcna_stability, "robust_bootstrap_stability": 0.5 * wgcna_stability["mean_bootstrap_ARI"] + 0.3 * wgcna_stability["median_bootstrap_ARI"] + 0.2 * wgcna_stability["minimum_bootstrap_ARI"], "runtime_seconds": wgcna_runtime}),
        ]
    )
    internal.to_csv(prelock / "internal_stability_summary.tsv", sep="\t", index=False)
    prelock_inputs = [discovery_path["psi"], discovery_path["events"], discovery_path["samples"], *args.gmt]
    for name in args.tuning:
        prelock_inputs.extend([paths[name]["psi"], paths[name]["samples"]])
    freeze_manifest = {
        "freeze_timestamp_utc": utc_now(),
        "phase": "discovery_and_tuning_frozen",
        "holdout_accessed": False,
        "discovery_event_selection": ">=75% finite PSI, non-zero variance, top variance selected in discovery only",
        "n_frozen_events": len(discovery),
        "inputs": [file_record(path, role="prelock") for path in prelock_inputs],
        "settings": {"network": asdict(network), "consensus": asdict(consensus), "benchmark": asdict(benchmark), "preservation": asdict(preservation), "multiobjective": asdict(objective)},
        "selected_parameters": {
            "stability_only": dict(zip(["edge_threshold", "resolution"], _selected_parameters(stability_result), strict=True)),
            "three_objective_cosplice": dict(zip(["edge_threshold", "resolution"], _selected_parameters(multi_result), strict=True)),
            "event_level_wgcna": {"soft_power": network.soft_power, "min_module_size": network.min_module_size, "merge_cut_height": benchmark.wgcna_merge_cut_height},
        },
        "assignment_sha256": {method: sha256_file(prelock / method / "discovery_assignments.tsv") for method in METHODS},
        "exact_wgcna_environment": wgcna_audit,
    }
    write_json(output / "freeze_manifest.json", freeze_manifest)
    print("Freeze manifest written. Opening holdout now.", flush=True)

    holdout_opened = utc_now()
    holdout_path = paths[args.holdout]
    holdout_raw = read_psi_matrix(holdout_path["psi"])
    holdout = holdout_raw.reindex(discovery.index)
    holdout = holdout.loc[holdout.notna().mean(axis=1).ge(0.75) & holdout.std(axis=1, skipna=True).gt(1e-8)]
    holdout = holdout.apply(lambda row: row.fillna(row.mean()), axis=1)
    holdout_samples = sample_frame(holdout_path["samples"], holdout.columns, args.holdout)
    holdout_samples.to_csv(postlock / "holdout_sample_metadata.tsv", sep="\t", index=False)
    holdout.to_csv(postlock / "holdout_locked_event_matrix.tsv", sep="\t")
    cohort = ExternalCohortConfig(args.holdout, holdout_path["psi"], holdout_path["samples"])

    preservation_tables: list[pd.DataFrame] = []
    preservation_summary_rows: list[dict[str, object]] = []
    for method, frame in assignments.items():
        result = evaluate_external_module_preservation(discovery, frame.rename(columns={"wgcna_module": "consensus_module"}), [cohort], preservation, progress=print)
        detail = result.module_preservation.copy()
        detail.insert(0, "method", method)
        preservation_tables.append(detail)
        summary = {
            "method": method,
            "tested_modules": len(detail),
            "median_holdout_Zsummary": float(detail["Zsummary_wgcna_like"].median()) if not detail.empty else np.nan,
            "median_holdout_program_score_preservation": float(detail["program_score_preservation"].median()) if not detail.empty else np.nan,
            "preserved_modules": int(detail["preservation_class"].eq("preserved").sum()) if not detail.empty else 0,
            "partially_preserved_modules": int(detail["preservation_class"].eq("partially_preserved").sum()) if not detail.empty else 0,
        }
        preservation_summary_rows.append(summary)
        method_dir = postlock / method
        method_dir.mkdir(exist_ok=True)
        detail.to_csv(method_dir / "holdout_module_preservation.tsv", sep="\t", index=False)

    holdout_observed = _assignment(holdout.index, discover_module_labels(holdout, partition_network), "module")
    holdout_assignments: dict[str, pd.DataFrame] = {}
    for method, discovery_result in (("stability_only", stability_result), ("three_objective_cosplice", multi_result)):
        threshold, resolution = _selected_parameters(discovery_result)
        fixed = replace(consensus, graph_edge_thresholds=[threshold], graph_resolutions=[resolution])
        print(f"Fitting frozen holdout partition for {method}", flush=True)
        result = build_bootstrap_coassignment_consensus(holdout, holdout_observed, holdout_samples, network, fixed, progress=print)
        holdout_assignments[method] = result.assignments[["event_key", "consensus_module"]].copy()
    holdout_wgcna, _ = _fit_method("wgcna_exact", holdout, network, benchmark, rscript)
    holdout_assignments["event_level_wgcna"] = _assignment(holdout.index, holdout_wgcna, "wgcna_module")
    for method, frame in holdout_assignments.items():
        frame.to_csv(postlock / method / "independent_holdout_assignments.tsv", sep="\t", index=False)

    gene_sets = load_gmt(merged_gmt)
    enrichment_tables: list[pd.DataFrame] = []
    for method in METHODS:
        enrichment_tables.append(local_module_enrichment("discovery", method, assignments[method], gene_sets))
        enrichment_tables.append(local_module_enrichment("holdout", method, holdout_assignments[method], gene_sets))
    enrichment = pd.concat(enrichment_tables, ignore_index=True)
    replication = functional_replication(enrichment)
    enrichment.to_csv(postlock / "local_GO_KEGG_enrichment.tsv", sep="\t", index=False)
    replication.to_csv(postlock / "independent_functional_replication.tsv", sep="\t", index=False)

    external_summary = pd.DataFrame(preservation_summary_rows)
    external_summary.to_csv(postlock / "holdout_preservation_summary.tsv", sep="\t", index=False)
    pd.concat(preservation_tables, ignore_index=True).to_csv(postlock / "holdout_module_preservation_all_methods.tsv", sep="\t", index=False)
    functional_summary = replication.groupby("method", observed=True).agg(
        median_best_term_jaccard=("best_term_jaccard", "median"),
        maximum_best_term_jaccard=("best_term_jaccard", "max"),
        discovery_modules_tested=("discovery_module", "nunique"),
        modules_with_replicated_terms=("shared_significant_terms", lambda values: int((values > 0).sum())),
    ).reset_index()
    summary = internal.merge(external_summary, on="method", how="left").merge(functional_summary, on="method", how="left")
    summary.to_csv(output / "locked_three_arm_summary.tsv", sep="\t", index=False)
    _render_dashboard(
        summary,
        output / "locked_three_arm_ablation.png",
        title=(
            f"Frozen three-arm {event_type} ablation: "
            f"{args.discovery} discovery to {args.holdout} holdout"
        ),
    )
    _write_report(output, summary, design)
    write_json(
        output / "holdout_access_audit.json",
        {
            "freeze_manifest_written_before_holdout_open": True,
            "holdout_opened_at_utc": holdout_opened,
            "evaluation_completed_at_utc": utc_now(),
            "holdout_inputs": [file_record(path, role="postlock_holdout") for path in holdout_path.values()],
            "n_exact_frozen_events_available_in_holdout": len(holdout),
            "generated_summary": str((output / "locked_three_arm_summary.tsv").resolve()),
        },
    )
    print(summary.to_string(index=False), flush=True)
    print(f"Output: {output}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
