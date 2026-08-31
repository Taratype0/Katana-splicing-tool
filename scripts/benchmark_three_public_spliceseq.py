from __future__ import annotations

import argparse
import json
import time
from dataclasses import asdict, replace
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.metrics import adjusted_rand_score

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from katana_cosplice.biological_interpretability import (
    evaluate_biological_interpretability,
    render_biological_interpretability,
)
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
        upper = matrix[np.ix_(indices, indices)][
            np.triu_indices(len(indices), k=1)
        ]
        values.append(float(upper.mean()))
        weights.append(len(upper))
    return float(np.average(values, weights=weights)) if values else 0.0


def _samples(columns: pd.Index, cohort: str) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_id": columns.astype(str),
            "condition": cohort,
            "cohort": cohort,
        }
    )


def _assignments(
    event_keys: pd.Index,
    labels: np.ndarray,
    column: str,
) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "event_key": event_keys.astype(str),
            column: labels.astype(str),
        }
    )


def _aggregate_preservation(method: str, result) -> dict[str, object]:
    table = result.module_preservation
    audit = (
        result.target_audit.iloc[0].to_dict()
        if not result.target_audit.empty
        else {}
    )
    if table.empty:
        return {
            "method": method,
            "status": audit.get("status", "no_preservable_modules"),
            "mapped_events": (
                int(result.event_mapping["event_key"].nunique())
                if not result.event_mapping.empty
                else 0
            ),
            "tested_modules": 0,
        }
    counts = table["preservation_class"].value_counts()
    return {
        "method": method,
        "status": audit.get("status", "finished"),
        "mapped_events": int(result.event_mapping["event_key"].nunique()),
        "tested_modules": int(len(table)),
        "preserved_modules": int(counts.get("preserved", 0)),
        "partially_preserved_modules": int(
            counts.get("partially_preserved", 0)
        ),
        "not_preserved_modules": int(counts.get("not_preserved", 0)),
        "median_Zsummary": float(table["Zsummary_wgcna_like"].median()),
        "median_density_ratio": float(table["density_ratio"].median()),
        "median_connectivity_correlation": float(
            table["connectivity_correlation"].median()
        ),
        "median_membership_profile_correlation": float(
            table["membership_profile_correlation"].median()
        ),
        "median_program_score_preservation": float(
            table["program_score_preservation"].median()
        ),
    }


def _load_matrix(path: Path, event_limit: int | None = None) -> pd.DataFrame:
    matrix = pd.read_csv(path, sep="\t", index_col=0).apply(
        pd.to_numeric,
        errors="coerce",
    )
    matrix.index = matrix.index.astype(str)
    matrix = matrix.loc[
        matrix.notna().mean(axis=1).ge(0.75)
        & matrix.std(axis=1).gt(1e-8)
    ]
    if event_limit is not None and len(matrix) > event_limit:
        matrix = matrix.loc[matrix.var(axis=1).nlargest(event_limit).index]
    return matrix.apply(lambda row: row.fillna(row.mean()), axis=1)


def _matrix_columns(path: Path) -> pd.Index:
    return pd.read_csv(path, sep="\t", nrows=0, index_col=0).columns


def _fit_discovery(
    name: str,
    matrix: pd.DataFrame,
    output: Path,
    *,
    network: NetworkConfig,
    consensus_settings: ConsensusConfig,
    benchmark: BenchmarkConfig,
    rscript: str,
    bootstraps: int,
) -> tuple[dict[str, np.ndarray], pd.DataFrame]:
    output.mkdir(parents=True, exist_ok=True)
    samples = _samples(matrix.columns, name)
    partition_network = replace(
        network,
        consensus_bootstrap_iterations=0,
        community_consensus_runs=1,
        bootstrap_iterations=0,
    )
    observed = _assignments(
        matrix.index,
        discover_module_labels(matrix, partition_network),
        "module",
    )

    start = time.perf_counter()
    cosplice = build_bootstrap_coassignment_consensus(
        matrix,
        observed,
        samples,
        network,
        consensus_settings,
        progress=print,
    )
    cosplice_seconds = time.perf_counter() - start
    labels = {
        "CoSplice": (
            cosplice.assignments.set_index("event_key")
            .reindex(matrix.index)["consensus_module"]
            .to_numpy(object)
        )
    }

    start = time.perf_counter()
    labels["Exact WGCNA"], _ = _fit_method(
        "wgcna_exact",
        matrix,
        network,
        benchmark,
        rscript,
    )
    wgcna_seconds = time.perf_counter() - start

    rng = np.random.default_rng(59)
    bootstrap_labels: dict[str, list[np.ndarray]] = {
        "CoSplice": [],
        "Exact WGCNA": [],
    }
    stability_rows: list[dict[str, object]] = []
    for iteration in range(1, bootstraps + 1):
        indices = _bootstrap_indices(samples, matrix.columns, rng)
        bootstrap = matrix.iloc[:, indices].copy()
        bootstrap.columns = [
            f"bootstrap_{iteration:03d}_{position:03d}"
            for position in range(len(indices))
        ]
        current = {
            "CoSplice": discover_module_labels(bootstrap, partition_network),
            "Exact WGCNA": _fit_method(
                "wgcna_exact",
                bootstrap,
                network,
                benchmark,
                rscript,
            )[0],
        }
        row: dict[str, object] = {"iteration": iteration}
        for method, values in current.items():
            bootstrap_labels[method].append(values)
            row[method] = adjusted_rand_score(labels[method], values)
        stability_rows.append(row)
        print(f"{name}: matched bootstrap {iteration}/{bootstraps}")

    stability = pd.DataFrame(stability_rows)
    neutral = _neutral_graph(matrix, benchmark.neutral_knn_neighbors)
    rows: list[dict[str, object]] = []
    for method, seconds in (
        ("CoSplice", cosplice_seconds),
        ("Exact WGCNA", wgcna_seconds),
    ):
        values = labels[method]
        coassignment = (
            cosplice.coassignment_matrix
            if method == "CoSplice"
            else _coassignment(bootstrap_labels[method])
        )
        row: dict[str, object] = {
            "cohort": name,
            "method": method,
            "n_events": len(matrix),
            "n_samples": matrix.shape[1],
            "runtime_seconds_initial_fit": seconds,
            "mean_bootstrap_ARI": float(stability[method].mean()),
            "median_bootstrap_ARI": float(stability[method].median()),
            "min_bootstrap_ARI": float(stability[method].min()),
            "weighted_mean_within_module_coassignment": _within_coassignment(
                values,
                coassignment,
            ),
        }
        row.update(_method_metrics(matrix, values, samples, neutral))
        row["valid_partition"] = bool(
            int(row["n_modules"]) >= 2
            and float(row["largest_module_fraction"]) <= 0.75
        )
        rows.append(row)

    internal = pd.DataFrame(rows)
    internal.to_csv(output / "internal_matched_bootstrap_benchmark.tsv", sep="\t", index=False)
    stability.to_csv(output / "matched_bootstrap_stability.tsv", sep="\t", index=False)
    cosplice.assignments.to_csv(output / "cosplice_assignments.tsv", sep="\t", index=False)
    _assignments(matrix.index, labels["Exact WGCNA"], "wgcna_module").to_csv(
        output / "wgcna_assignments.tsv",
        sep="\t",
        index=False,
    )
    cosplice.parameter_search.to_csv(
        output / "cosplice_parameter_search.tsv",
        sep="\t",
        index=False,
    )
    render_consensus_plots(cosplice, output / "cosplice_plots")
    return labels, internal


def _load_discovery_cache(
    matrix: pd.DataFrame,
    output: Path,
) -> tuple[dict[str, np.ndarray], pd.DataFrame] | None:
    cosplice_path = output / "cosplice_assignments.tsv"
    wgcna_path = output / "wgcna_assignments.tsv"
    internal_path = output / "internal_matched_bootstrap_benchmark.tsv"
    if not all(path.exists() for path in (cosplice_path, wgcna_path, internal_path)):
        return None
    cosplice = pd.read_csv(cosplice_path, sep="\t").set_index("event_key")
    wgcna = pd.read_csv(wgcna_path, sep="\t").set_index("event_key")
    if not matrix.index.equals(cosplice.index) or not matrix.index.equals(wgcna.index):
        return None
    internal = pd.read_csv(internal_path, sep="\t")
    if internal.empty or not pd.to_numeric(
        internal["n_events"],
        errors="coerce",
    ).eq(len(matrix)).all():
        return None
    return (
        {
            "CoSplice": cosplice["consensus_module"].astype(str).to_numpy(object),
            "Exact WGCNA": wgcna["wgcna_module"].astype(str).to_numpy(object),
        },
        internal,
    )


def _selection_summary(
    internal: pd.DataFrame,
    external: pd.DataFrame,
    interpretability: pd.DataFrame | None = None,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for method, group in internal.groupby("method", observed=True):
        target = external.loc[external["method"].eq(method)]
        best_z = pd.to_numeric(target["best_module_Zsummary"], errors="coerce")
        rows.append(
            {
                "method": method,
                "valid_discovery_fraction": float(group["valid_partition"].mean()),
                "median_min_bootstrap_ARI": float(group["min_bootstrap_ARI"].median()),
                "worst_min_bootstrap_ARI": float(group["min_bootstrap_ARI"].min()),
                "median_silhouette_cosine": float(group["silhouette_cosine"].median()),
                "median_best_external_Zsummary": float(best_z.median()),
                "worst_best_external_Zsummary": float(best_z.min()),
                "external_projection_count": int(len(target)),
            }
        )
    result = pd.DataFrame(rows)
    if interpretability is not None and not interpretability.empty:
        result = result.merge(
            interpretability,
            on="method",
            how="left",
            validate="one_to_one",
        )
    for column in (
        "fraction_modules_not_gene_dominated",
        "functional_annotated_module_fraction",
        "median_best_cross_cohort_term_jaccard",
        "median_best_cross_cohort_gene_jaccard",
    ):
        if column not in result:
            result[column] = np.nan
    result = result.sort_values(
        [
            "valid_discovery_fraction",
            "worst_best_external_Zsummary",
            "median_best_external_Zsummary",
            "fraction_modules_not_gene_dominated",
            "functional_annotated_module_fraction",
            "median_best_cross_cohort_gene_jaccard",
            "median_best_cross_cohort_term_jaccard",
            "worst_min_bootstrap_ARI",
        ],
        ascending=False,
        na_position="last",
    ).reset_index(drop=True)
    result["recommended_production_method"] = False
    if not result.empty:
        result.loc[0, "recommended_production_method"] = True
    result["selection_rule"] = (
        "Lexicographic: valid multi-module discovery fraction, worst external best-module "
        "Zsummary, median external best-module Zsummary, resistance to single-gene event "
        "domination, functional annotation/replication, then worst bootstrap ARI. Functional "
        "evidence validates but never defines modules."
    )
    return result


def _render_dashboard(
    internal: pd.DataFrame,
    external: pd.DataFrame,
    selection: pd.DataFrame,
    output: Path,
    *,
    cohort_count: int,
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    sns.barplot(
        data=internal,
        x="cohort",
        y="min_bootstrap_ARI",
        hue="method",
        ax=axes[0],
    )
    axes[0].set_ylim(0, 1)
    axes[0].set_title("Worst-case internal stability")
    axes[0].set_ylabel("Minimum bootstrap ARI")
    axes[0].set_xlabel("")

    heatmap = external.pivot_table(
        index="discovery_cohort",
        columns="external_cohort",
        values="best_module_Zsummary",
        aggfunc="max",
    )
    sns.heatmap(
        heatmap,
        annot=True,
        fmt=".2f",
        cmap="YlGnBu",
        ax=axes[1],
        cbar_kws={"label": "Best module Zsummary"},
    )
    axes[1].set_title("CoSplice external preservation")
    axes[1].set_xlabel("Independent validation cohort")
    axes[1].set_ylabel("Discovery cohort")

    plot = selection.melt(
        id_vars="method",
        value_vars=[
            "median_best_external_Zsummary",
            "worst_best_external_Zsummary",
        ],
        var_name="metric",
        value_name="value",
    )
    sns.barplot(data=plot, x="method", y="value", hue="metric", ax=axes[2])
    axes[2].axhline(2, color="#D97706", linestyle="--", linewidth=1)
    axes[2].set_title("Cross-cohort robustness")
    axes[2].set_ylabel("Best-module Zsummary")
    axes[2].set_xlabel("")
    axes[2].tick_params(axis="x", rotation=15)
    fig.suptitle(
        f"{cohort_count}-cohort public CoSplice benchmark",
        fontsize=16,
        fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Fit each public PSI cohort once and evaluate all directed external-cohort "
            "preservation projections."
        )
    )
    parser.add_argument(
        "--cohort",
        action="append",
        nargs=2,
        metavar=("NAME", "PSI_MATRIX"),
        required=True,
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--events", type=int, default=500)
    parser.add_argument("--bootstraps", type=int, default=20)
    parser.add_argument("--permutations", type=int, default=500)
    parser.add_argument("--min-module-size", type=int, default=20)
    parser.add_argument("--organism", default="hsapiens")
    parser.add_argument(
        "--skip-functional-enrichment",
        action="store_true",
        help="Run the local gene-diversity/RBP audit without querying g:Profiler.",
    )
    args = parser.parse_args()
    if len(args.cohort) < 3:
        raise ValueError("At least three independent cohorts are required.")

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    cohort_paths = {
        str(name): Path(path).resolve()
        for name, path in args.cohort
    }
    # Event selection is discovery-only. External projection must retain every
    # eligible event so preservation is not underestimated by intersecting two
    # independently selected top-variance lists.
    matrices = {
        name: _load_matrix(path, args.events)
        for name, path in cohort_paths.items()
    }
    network = NetworkConfig(
        min_module_size=args.min_module_size,
        knn_neighbors=min(50, max(10, args.events // 10)),
        consensus_knn_neighbors=[25, 50, 75],
        consensus_bootstrap_iterations=0,
        community_consensus_runs=1,
        bootstrap_iterations=0,
        edge_threshold=0.05,
        random_seed=17,
    )
    consensus_settings = ConsensusConfig(
        bootstrap_iterations=args.bootstraps,
        min_module_size=args.min_module_size,
        graph_top_k=min(50, max(10, args.events - 1)),
        dense_max_events=max(1_000, args.events),
        random_seed=37,
    )
    benchmark = BenchmarkConfig(
        methods=["wgcna_exact"],
        max_events=args.events,
        bootstrap_iterations=args.bootstraps,
        random_seed=41,
    )
    preservation = PreservationConfig(
        permutations=args.permutations,
        min_shared_events=max(10, args.min_module_size),
        dense_max_events=max(1_000, args.events),
        random_seed=43,
    )
    rscript, wgcna_audit = _rscript_with_wgcna()
    if rscript is None:
        raise RuntimeError(wgcna_audit)

    labels_by_cohort: dict[str, dict[str, np.ndarray]] = {}
    internal_frames: list[pd.DataFrame] = []
    for name, matrix in matrices.items():
        discovery_dir = output / "discovery" / name
        cached = _load_discovery_cache(matrix, discovery_dir)
        if cached is not None:
            print(f"{name}: discovery cache hit")
            labels, internal = cached
        else:
            labels, internal = _fit_discovery(
                name,
                matrix,
                discovery_dir,
                network=network,
                consensus_settings=consensus_settings,
                benchmark=benchmark,
                rscript=rscript,
                bootstraps=args.bootstraps,
            )
        labels_by_cohort[name] = labels
        internal_frames.append(internal)

    external_rows: list[dict[str, object]] = []
    detail_frames: list[pd.DataFrame] = []
    for external_name in matrices:
        target_columns = _matrix_columns(cohort_paths[external_name])
        for discovery_name in matrices:
            if discovery_name == external_name:
                continue
            discovery = matrices[discovery_name]
            pair_dir = output / "preservation" / f"{discovery_name}_to_{external_name}"
            pair_dir.mkdir(parents=True, exist_ok=True)
            sample_path = pair_dir / "external_sample_metadata.tsv"
            _samples(target_columns, external_name).to_csv(
                sample_path,
                sep="\t",
                index=False,
            )
            cohort = ExternalCohortConfig(
                external_name,
                cohort_paths[external_name],
                sample_path,
            )
            for method, labels in labels_by_cohort[discovery_name].items():
                method_dir = pair_dir / method.lower().replace(" ", "_")
                method_dir.mkdir(parents=True, exist_ok=True)
                result = evaluate_external_module_preservation(
                    discovery,
                    _assignments(discovery.index, labels, "consensus_module"),
                    [cohort],
                    preservation,
                    progress=print,
                )
                detail = result.module_preservation.copy()
                detail.insert(0, "method", method)
                detail.insert(0, "external_cohort", external_name)
                detail.insert(0, "discovery_cohort", discovery_name)
                detail_frames.append(detail)
                summary = _aggregate_preservation(method, result)
                tested = result.module_preservation
                summary.update(
                    {
                        "discovery_cohort": discovery_name,
                        "external_cohort": external_name,
                        "best_module_Zsummary": (
                            float(tested["Zsummary_wgcna_like"].max())
                            if not tested.empty
                            else np.nan
                        ),
                    }
                )
                external_rows.append(summary)
                result.module_preservation.to_csv(
                    method_dir / "external_module_preservation.tsv",
                    sep="\t",
                    index=False,
                )
                result.event_mapping.to_csv(
                    method_dir / "external_event_mapping.tsv",
                    sep="\t",
                    index=False,
                )
                result.target_audit.to_csv(
                    method_dir / "external_target_audit.tsv",
                    sep="\t",
                    index=False,
                )
                render_preservation_plots(result, method_dir / "plots")

    internal = pd.concat(internal_frames, ignore_index=True)
    external = pd.DataFrame(external_rows)
    detail = (
        pd.concat(detail_frames, ignore_index=True)
        if detail_frames
        else pd.DataFrame()
    )
    assignment_frames: dict[str, dict[str, pd.DataFrame]] = {}
    for cohort_name, matrix in matrices.items():
        assignment_frames[cohort_name] = {
            "CoSplice": _assignments(
                matrix.index,
                labels_by_cohort[cohort_name]["CoSplice"],
                "consensus_module",
            ),
            "Exact WGCNA": _assignments(
                matrix.index,
                labels_by_cohort[cohort_name]["Exact WGCNA"],
                "wgcna_module",
            ),
        }
    biology = evaluate_biological_interpretability(
        assignment_frames,
        organism=args.organism,
        run_functional_enrichment=not args.skip_functional_enrichment,
    )
    biology.module_gene_audit.to_csv(
        output / "module_biological_interpretability.tsv",
        sep="\t",
        index=False,
    )
    biology.functional_enrichment.to_csv(
        output / "module_go_kegg_enrichment.tsv",
        sep="\t",
        index=False,
    )
    biology.enrichment_sources.to_csv(
        output / "biological_interpretability_sources.tsv",
        sep="\t",
        index=False,
    )
    biology.functional_theme_replication.to_csv(
        output / "cross_cohort_functional_theme_replication.tsv",
        sep="\t",
        index=False,
    )
    biology.module_gene_replication.to_csv(
        output / "cross_cohort_module_gene_replication.tsv",
        sep="\t",
        index=False,
    )
    biology.method_summary.to_csv(
        output / "biological_interpretability_summary.tsv",
        sep="\t",
        index=False,
    )
    render_biological_interpretability(
        biology,
        output / "biological_interpretability_benchmark.png",
    )
    selection = _selection_summary(internal, external, biology.method_summary)
    internal.to_csv(output / "cohort_internal_benchmark.tsv", sep="\t", index=False)
    external.to_csv(
        output / "directed_projection_preservation_summary.tsv",
        sep="\t",
        index=False,
    )
    detail.to_csv(
        output / "directed_projection_module_preservation.tsv",
        sep="\t",
        index=False,
    )
    if len(matrices) == 3:
        internal.to_csv(
            output / "three_cohort_internal_benchmark.tsv",
            sep="\t",
            index=False,
        )
        external.to_csv(
            output / "six_projection_preservation_summary.tsv",
            sep="\t",
            index=False,
        )
        detail.to_csv(
            output / "six_projection_module_preservation.tsv",
            sep="\t",
            index=False,
        )
    selection.to_csv(output / "recommended_single_method.tsv", sep="\t", index=False)
    parameters = {
        "cohorts": {name: str(path) for name, path in cohort_paths.items()},
        "design": (
            f"{len(matrices)}-cohort rotating discovery with "
            f"{len(matrices) * (len(matrices) - 1)} directed external projections"
        ),
        "development_boundary": (
            "Rotating discovery cohorts are used for method-development stress testing. "
            "A separately locked cohort is still required for final confirmatory validation."
        ),
        "event_selection": (
            f"Within-discovery variance ranking after >=75% non-missing filtering; top {args.events}."
        ),
        "bootstraps": args.bootstraps,
        "permutations_per_module": args.permutations,
        "network": asdict(network),
        "consensus": asdict(consensus_settings),
        "preservation": asdict(preservation),
        "exact_wgcna": wgcna_audit,
        "biological_interpretability": {
            "organism": args.organism,
            "functional_enrichment_enabled": not args.skip_functional_enrichment,
            "principle": (
                "Gene diversity, RBP/splicing-factor support, GO/KEGG enrichment and "
                "cross-cohort functional-theme replication are validation layers. They are "
                "not used to discover or tune modules."
            ),
            "single_gene_domination_cutoff": 0.20,
        },
    }
    (output / "parameters.json").write_text(
        json.dumps(parameters, indent=2),
        encoding="utf-8",
    )
    _render_dashboard(
        internal,
        external,
        selection,
        output / "public_cohort_benchmark.png",
        cohort_count=len(matrices),
    )
    if len(matrices) == 3:
        _render_dashboard(
            internal,
            external,
            selection,
            output / "three_cohort_public_benchmark.png",
            cohort_count=3,
        )
    print("\nRecommended single production method")
    print(selection.to_string(index=False))
    print(f"\nOutput: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
