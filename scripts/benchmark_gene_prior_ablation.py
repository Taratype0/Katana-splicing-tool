from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from katana_cosplice.config import (
    ConsensusConfig,
    GenePriorConfig,
    NetworkConfig,
)
from katana_cosplice.consensus import build_bootstrap_coassignment_consensus
from katana_cosplice.gene_prior import build_gene_mantel_prior, event_gene_id
from katana_cosplice.network import fit_network


def _event_universe(path: Path, limit: int) -> pd.DataFrame:
    matrix = pd.read_csv(path, sep="\t", index_col=0)
    numeric = matrix.apply(pd.to_numeric, errors="coerce")
    numeric = numeric.loc[
        numeric.notna().mean(axis=1).ge(0.50)
        & numeric.std(axis=1, skipna=True).gt(1e-8)
    ]
    if len(numeric) > limit:
        numeric = numeric.loc[numeric.var(axis=1, skipna=True).nlargest(limit).index]
    return numeric.apply(lambda row: row.fillna(row.median()), axis=1)


def _gene_diversity(assignments: pd.DataFrame) -> tuple[float, float]:
    current = assignments.loc[assignments["consensus_module"].ne("grey")].copy()
    if current.empty:
        return np.nan, np.nan
    current["gene_id"] = current["event_key"].map(event_gene_id)
    module = current.groupby("consensus_module", observed=True)
    unique = module["gene_id"].nunique()
    top_fraction = module["gene_id"].apply(
        lambda values: values.value_counts().iloc[0] / len(values)
    )
    return float(unique.median()), float(top_fraction.max())


def _plot(summary: pd.DataFrame, output: Path) -> None:
    metrics = [
        ("median_bootstrap_ARI", "Median bootstrap ARI", (0, 1)),
        ("grey_fraction", "Grey fraction", (0, 1)),
        ("mean_consensus_confidence", "Mean event consensus", (0, 1)),
        ("median_unique_genes_per_module", "Median genes / module", None),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.5))
    labels = [f"{value:.2f}" for value in summary["gene_prior_weight"]]
    for axis, (column, title, ylim) in zip(axes, metrics, strict=True):
        values = pd.to_numeric(summary[column], errors="coerce")
        bars = axis.bar(labels, values, color="#087E8B", edgecolor="#17324D")
        axis.set_title(title, loc="left", fontweight="bold")
        axis.set_xlabel("Gene-prior weight")
        axis.spines[["top", "right"]].set_visible(False)
        if ylim:
            axis.set_ylim(*ylim)
        for bar, value in zip(bars, values, strict=True):
            axis.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height(),
                "NA" if not np.isfinite(value) else f"{value:.3g}",
                ha="center",
                va="bottom",
                fontsize=9,
            )
    fig.suptitle(
        "GSE62669 locked event universe: AS-only versus Mantel-prior CoSplice",
        fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Ablate bounded gene-level Mantel reinforcement in CoSplice."
    )
    parser.add_argument("--event-matrix", type=Path, required=True)
    parser.add_argument("--sample-metadata", type=Path, required=True)
    parser.add_argument("--exon-counts", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--events", type=int, default=500)
    parser.add_argument("--bootstraps", type=int, default=20)
    parser.add_argument("--weights", default="0,0.1,0.2,0.3")
    args = parser.parse_args()

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    matrix = _event_universe(args.event_matrix.resolve(), args.events)
    samples = pd.read_csv(args.sample_metadata.resolve(), sep="\t")
    samples = samples.loc[samples["sample_id"].isin(matrix.columns)].copy()
    samples["sample_id"] = pd.Categorical(
        samples["sample_id"], categories=matrix.columns, ordered=True
    )
    samples = samples.sort_values("sample_id").reset_index(drop=True)
    samples["sample_id"] = samples["sample_id"].astype("object")
    if samples["sample_id"].tolist() != matrix.columns.tolist():
        raise ValueError("Sample metadata could not be aligned to the event matrix.")

    prior_settings = GenePriorConfig(
        enabled=True,
        exon_count_matrix=args.exon_counts.resolve(),
        weight=0.2,
    )
    prior = build_gene_mantel_prior(matrix.index, matrix.columns, prior_settings)
    qc = pd.DataFrame(
        {
            "stability_score": matrix.var(axis=1).rank(pct=True),
        },
        index=matrix.index,
    )
    network = NetworkConfig(
        min_module_size=20,
        knn_neighbors=50,
        consensus_knn_neighbors=[25, 50, 75],
        consensus_bootstrap_iterations=0,
        community_consensus_runs=1,
        bootstrap_iterations=0,
        edge_threshold=0.05,
        random_seed=17,
    )
    consensus_settings = ConsensusConfig(
        bootstrap_iterations=args.bootstraps,
        min_module_size=20,
        graph_top_k=50,
        dense_max_events=max(1_000, len(matrix)),
        random_seed=37,
    )
    rows: list[dict[str, object]] = []
    for weight in [float(value) for value in args.weights.split(",")]:
        active_prior = prior if weight > 0 else None
        fitted = fit_network(
            matrix,
            qc,
            network,
            samples=samples,
            gene_prior=active_prior,
            gene_prior_weight=weight,
        )
        consensus = build_bootstrap_coassignment_consensus(
            fitted.selected_matrix,
            fitted.assignments,
            samples,
            network,
            consensus_settings,
            gene_prior=active_prior,
            gene_prior_weight=weight,
        )
        label = f"weight_{weight:.2f}".replace(".", "p")
        consensus.assignments.to_csv(
            output / f"{label}_consensus_assignments.tsv",
            sep="\t",
            index=False,
        )
        fitted.edges.to_csv(
            output / f"{label}_network_edges.tsv",
            sep="\t",
            index=False,
        )
        median_genes, max_top_gene = _gene_diversity(consensus.assignments)
        non_grey = consensus.assignments["consensus_module"].ne("grey")
        rows.append(
            {
                "gene_prior_weight": weight,
                "n_events": len(consensus.assignments),
                "n_modules": int(
                    consensus.assignments.loc[
                        non_grey, "consensus_module"
                    ].nunique()
                ),
                "grey_fraction": float((~non_grey).mean()),
                "median_bootstrap_ARI": float(
                    consensus.bootstrap_runs[
                        "adjusted_rand_index_vs_observed"
                    ].median()
                ),
                "mean_consensus_confidence": float(
                    consensus.assignments["consensus_membership_confidence"].mean()
                ),
                "median_unique_genes_per_module": median_genes,
                "worst_top_gene_event_fraction": max_top_gene,
                "gene_prior_event_coverage": float(
                    fitted.assignments["gene_prior_available"].mean()
                ),
                "gene_prior_edge_coverage": float(
                    pd.to_numeric(
                        fitted.edges.get("gene_mantel_prior"), errors="coerce"
                    ).gt(0).mean()
                )
                if weight > 0
                else 0.0,
            }
        )
    summary = pd.DataFrame(rows)
    summary.to_csv(output / "gene_prior_ablation_summary.tsv", sep="\t", index=False)
    prior.gene_audit.to_csv(output / "gene_prior_gene_audit.tsv", sep="\t", index=False)
    _plot(summary, output / "gene_prior_ablation_scorecard.png")
    best = summary.sort_values(
        ["median_bootstrap_ARI", "mean_consensus_confidence", "gene_prior_weight"],
        ascending=[False, False, True],
        kind="mergesort",
    ).iloc[0]
    baseline = summary.loc[summary["gene_prior_weight"].eq(0.0)].iloc[0]
    interpretation = (
        "# Gene-level Mantel scaffold ablation\n\n"
        "The event universe and sample bootstraps are locked across weights. "
        "GO terms and phenotype labels are not used for weight selection.\n\n"
        f"- Best internal-stability weight: {best['gene_prior_weight']:.2f}\n"
        f"- Median bootstrap ARI: {baseline['median_bootstrap_ARI']:.3f} "
        f"(AS only) to {best['median_bootstrap_ARI']:.3f} (best weight)\n"
        f"- Mean consensus confidence: {baseline['mean_consensus_confidence']:.3f} "
        f"to {best['mean_consensus_confidence']:.3f}\n"
        f"- Grey fraction: {baseline['grey_fraction']:.3f} "
        f"to {best['grey_fraction']:.3f}\n\n"
        "Interpretation: retain the scaffold as optional unless improvement is "
        "reproduced in an independent cohort. A higher internal score alone does "
        "not establish biological validity or universal superiority.\n"
    )
    (output / "gene_prior_ablation_interpretation.md").write_text(
        interpretation,
        encoding="utf-8",
    )
    (output / "parameters.json").write_text(
        json.dumps(
            {
                "event_matrix": str(args.event_matrix.resolve()),
                "sample_metadata": str(args.sample_metadata.resolve()),
                "exon_counts": str(args.exon_counts.resolve()),
                "event_limit": args.events,
                "bootstraps": args.bootstraps,
                "weights": args.weights,
                "selection": "top variance after >=50% observed and non-constant filtering",
                "selection_tuning_boundary": (
                    "Weights are compared on internal stability and structural audit only; "
                    "GO/phenotype labels are not used to select a weight."
                ),
                "prior_mode": prior.prior_mode,
                "best_internal_stability_weight": float(best["gene_prior_weight"]),
                "generated_files": [
                    "gene_prior_ablation_summary.tsv",
                    "gene_prior_gene_audit.tsv",
                    "gene_prior_ablation_scorecard.png",
                    "gene_prior_ablation_interpretation.md",
                ],
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
