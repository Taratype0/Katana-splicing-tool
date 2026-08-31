from __future__ import annotations

import argparse
import importlib.util
import json
from dataclasses import asdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from katana_cosplice.biological_interpretability import (
    evaluate_biological_interpretability,
)
from katana_cosplice.benchmark import _rscript_with_wgcna
from katana_cosplice.config import BenchmarkConfig, ConsensusConfig, NetworkConfig


def _load_script_module(name: str, filename: str):
    path = Path(__file__).with_name(filename)
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load benchmark helper: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_external = _load_script_module(
    "katana_benchmark_external_cosplice",
    "benchmark_external_cosplice.py",
)
_three_cohort = _load_script_module(
    "katana_benchmark_three_public_spliceseq",
    "benchmark_three_public_spliceseq.py",
)
load_gse141633 = _external.load_gse141633
load_gse85908 = _external.load_gse85908
_fit_discovery = _three_cohort._fit_discovery
_load_discovery_cache = _three_cohort._load_discovery_cache


def _select_event_universe(matrix: pd.DataFrame, event_limit: int) -> pd.DataFrame:
    numeric = matrix.apply(pd.to_numeric, errors="coerce")
    numeric = numeric.loc[
        numeric.notna().mean(axis=1).ge(0.50)
        & numeric.std(axis=1, skipna=True).gt(1e-8)
    ]
    if len(numeric) > event_limit:
        numeric = numeric.loc[numeric.var(axis=1, skipna=True).nlargest(event_limit).index]
    return numeric.apply(lambda row: row.fillna(row.median()), axis=1)


def _load_prepared_matrix(path: Path, event_limit: int) -> pd.DataFrame:
    matrix = pd.read_csv(path, sep="\t", index_col=0)
    return _select_event_universe(matrix, event_limit)


def _interpretability(
    cohort: str,
    matrix: pd.DataFrame,
    labels: dict[str, np.ndarray],
    *,
    organism: str | None,
    run_functional_enrichment: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if organism is None:
        rows = []
        for method, values in labels.items():
            modules = set(values.astype(str)).difference({"grey"})
            rows.append(
                {
                    "cohort": cohort,
                    "method": method,
                    "n_modules_audited": len(modules),
                    "median_unique_genes_per_module": np.nan,
                    "median_effective_gene_fraction": np.nan,
                    "worst_top_gene_event_fraction": np.nan,
                    "fraction_modules_not_gene_dominated": np.nan,
                    "rbp_enriched_module_count_FDR_0.05": np.nan,
                    "significant_GO_KEGG_term_count": np.nan,
                    "functional_annotated_module_fraction": np.nan,
                    "median_best_cross_cohort_term_jaccard": np.nan,
                    "median_best_cross_cohort_gene_jaccard": np.nan,
                    "functional_mapping_status": "disabled:no_gene_annotation",
                }
            )
        return pd.DataFrame(rows), pd.DataFrame(), pd.DataFrame(
            [
                {
                    "cohort": cohort,
                    "method": method,
                    "status": "disabled:no_gene_annotation",
                }
                for method in labels
            ]
        )
    assignments = {
        cohort: {
            "CoSplice": pd.DataFrame(
                {
                    "event_key": matrix.index.astype(str),
                    "consensus_module": labels["CoSplice"].astype(str),
                }
            ),
            "Exact WGCNA": pd.DataFrame(
                {
                    "event_key": matrix.index.astype(str),
                    "wgcna_module": labels["Exact WGCNA"].astype(str),
                }
            ),
        }
    }
    result = evaluate_biological_interpretability(
        assignments,
        organism=organism,
        run_functional_enrichment=run_functional_enrichment,
    )
    summary = result.method_summary.copy()
    if not summary.empty:
        summary.insert(0, "cohort", cohort)
        summary["functional_mapping_status"] = (
            (
                f"enabled:{organism}"
                if run_functional_enrichment
                else f"gene_audit_only:{organism}"
            )
        )
    enrichment = result.functional_enrichment.copy()
    source = result.enrichment_sources.copy()
    return summary, enrichment, source


def _render_scorecard(internal: pd.DataFrame, biology: pd.DataFrame, output: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    sns.barplot(
        data=internal,
        x="cohort",
        y="min_bootstrap_ARI",
        hue="method",
        ax=axes[0],
    )
    axes[0].set_ylim(0, 1)
    axes[0].set_title("Worst matched-bootstrap stability")
    axes[0].set_ylabel("Minimum adjusted Rand index")
    axes[0].set_xlabel("")
    axes[0].tick_params(axis="x", rotation=18)

    sns.barplot(
        data=internal,
        x="cohort",
        y="valid_partition",
        hue="method",
        ax=axes[1],
    )
    axes[1].set_ylim(0, 1.05)
    axes[1].set_title("Scientifically usable partition")
    axes[1].set_ylabel("Valid multi-module partition")
    axes[1].set_xlabel("")
    axes[1].tick_params(axis="x", rotation=18)

    if biology.empty:
        axes[2].text(0.5, 0.5, "No biological mapping available", ha="center", va="center")
        axes[2].axis("off")
    else:
        sns.barplot(
            data=biology,
            x="cohort",
            y="fraction_modules_not_gene_dominated",
            hue="method",
            ax=axes[2],
        )
        axes[2].set_ylim(0, 1.05)
        axes[2].set_title("Resistance to single-gene pseudo-modules")
        axes[2].set_ylabel("Fraction of modules passing audit")
        axes[2].set_xlabel("")
        axes[2].tick_params(axis="x", rotation=18)
    fig.suptitle(
        "CoSplice external challenge set: matched input, no cross-catalogue claims",
        fontsize=15,
        fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Benchmark final CoSplice and exact WGCNA on three independent public "
            "splicing datasets without treating incompatible event catalogues as shared events."
        )
    )
    parser.add_argument("--gse62669-matrix", type=Path, required=True)
    parser.add_argument("--gse141633-psi", type=Path, required=True)
    parser.add_argument("--gse85908-splicing", type=Path, required=True)
    parser.add_argument("--gse85908-metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--events", type=int, default=500)
    parser.add_argument("--bootstraps", type=int, default=20)
    parser.add_argument("--min-module-size", type=int, default=20)
    parser.add_argument(
        "--skip-functional-enrichment",
        action="store_true",
        help="Keep gene-diversity audits but skip GO/KEGG web queries.",
    )
    args = parser.parse_args()

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    gse141633, _, _, _ = load_gse141633(args.gse141633_psi.resolve())
    gse85908, _, _, _ = load_gse85908(
        args.gse85908_splicing.resolve(),
        args.gse85908_metadata.resolve(),
    )
    datasets = {
        "GSE62669_mouse_brain_exon_usage": {
            "matrix": _load_prepared_matrix(args.gse62669_matrix.resolve(), args.events),
            "organism": "mmusculus",
            "biological_scope": "gene IDs available; exon-usage proxy, not rMATS PSI",
        },
        "GSE141633_TCGA_PRAD_bulk_PSI": {
            "matrix": _select_event_universe(gse141633, args.events),
            "organism": "hsapiens",
            "biological_scope": "Ensembl gene IDs available",
        },
        "GSE85908_single_cell_PSI_boundary": {
            "matrix": _select_event_universe(gse85908, args.events),
            "organism": None,
            "biological_scope": "coordinate-only event IDs; functional enrichment disabled",
        },
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
    consensus = ConsensusConfig(
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
    rscript, wgcna_audit = _rscript_with_wgcna()
    if rscript is None:
        raise RuntimeError(wgcna_audit)

    internal_frames: list[pd.DataFrame] = []
    biology_frames: list[pd.DataFrame] = []
    enrichment_frames: list[pd.DataFrame] = []
    source_frames: list[pd.DataFrame] = []
    dataset_audit: list[dict[str, object]] = []
    for name, settings in datasets.items():
        matrix = settings["matrix"]
        assert isinstance(matrix, pd.DataFrame)
        cohort_dir = output / name
        cached = _load_discovery_cache(matrix, cohort_dir)
        if cached is None:
            labels, internal = _fit_discovery(
                name,
                matrix,
                cohort_dir,
                network=network,
                consensus_settings=consensus,
                benchmark=benchmark,
                rscript=rscript,
                bootstraps=args.bootstraps,
            )
        else:
            labels, internal = cached
            print(f"{name}: discovery cache hit")
        internal_frames.append(internal)
        organism = settings["organism"]
        summary, enrichment, sources = _interpretability(
            name,
            matrix,
            labels,
            organism=organism if isinstance(organism, str) else None,
            run_functional_enrichment=not args.skip_functional_enrichment,
        )
        biology_frames.append(summary)
        if not enrichment.empty:
            enrichment_frames.append(enrichment)
        source_frames.append(sources)
        dataset_audit.append(
            {
                "dataset": name,
                "n_events": len(matrix),
                "n_samples": matrix.shape[1],
                "biological_scope": settings["biological_scope"],
                "cross_dataset_preservation_eligible": False,
                "reason": (
                    "Independent heterogeneous event catalogues; evaluate each dataset "
                    "head-to-head but do not fabricate exact cross-catalogue event matching."
                ),
            }
        )

    internal = pd.concat(internal_frames, ignore_index=True)
    biology = pd.concat(biology_frames, ignore_index=True)
    enrichment = (
        pd.concat(enrichment_frames, ignore_index=True)
        if enrichment_frames
        else pd.DataFrame()
    )
    sources = pd.concat(source_frames, ignore_index=True)
    internal.to_csv(output / "three_dataset_internal_head_to_head.tsv", sep="\t", index=False)
    biology.to_csv(output / "three_dataset_biological_audit.tsv", sep="\t", index=False)
    enrichment.to_csv(output / "three_dataset_go_kegg_enrichment.tsv", sep="\t", index=False)
    sources.to_csv(output / "three_dataset_enrichment_sources.tsv", sep="\t", index=False)
    pd.DataFrame(dataset_audit).to_csv(
        output / "three_dataset_catalogue_compatibility.tsv",
        sep="\t",
        index=False,
    )
    _render_scorecard(
        internal,
        biology,
        output / "three_dataset_external_challenge_scorecard.png",
    )
    (output / "parameters.json").write_text(
        json.dumps(
            {
                "datasets": {
                    "GSE62669": str(args.gse62669_matrix.resolve()),
                    "GSE141633": str(args.gse141633_psi.resolve()),
                    "GSE85908_splicing": str(args.gse85908_splicing.resolve()),
                    "GSE85908_metadata": str(args.gse85908_metadata.resolve()),
                },
                "events": args.events,
                "bootstraps": args.bootstraps,
                "network": asdict(network),
                "consensus": asdict(consensus),
                "exact_wgcna": wgcna_audit,
                "interpretation_boundary": (
                    "Three independent public datasets challenge general computational "
                    "behaviour. They do not provide exact cross-dataset module preservation "
                    "because their event definitions and platforms are not interchangeable."
                ),
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(internal.to_string(index=False))
    print(f"\nOutput: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
