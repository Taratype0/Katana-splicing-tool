from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.metrics import adjusted_rand_score

from katana_cosplice.config import ConsensusConfig, NetworkConfig, UncertaintyConfig
from katana_cosplice.consensus import _confidence_tables, _select_final_modules
from katana_cosplice.network import discover_module_labels
from katana_cosplice.uncertainty import estimate_beta_binomial_uncertainty


ARM_RAW = "raw_PSI_WGCNA"
ARM_UNCERTAINTY = "WGCNA_plus_read_uncertainty"
ARM_FULL = "full_CoSplice"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_matrix(path: Path, index: pd.Index | None = None) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", index_col=0, low_memory=False)
    frame.index = frame.index.astype(str)
    if index is not None:
        frame = frame.reindex(index)
    return frame.apply(pd.to_numeric, errors="coerce")


def _bootstrap_indices(
    samples: pd.DataFrame,
    columns: pd.Index,
    rng: np.random.Generator,
) -> np.ndarray:
    ordered = samples.set_index("sample_id").reindex(columns)
    groups = ordered.groupby("condition", sort=False, observed=True).indices
    return np.concatenate(
        [rng.choice(indices, size=len(indices), replace=True) for indices in groups.values()]
    )


def _expanded_labels(
    universe: pd.Index,
    fitted_events: pd.Index,
    labels: np.ndarray,
) -> pd.Series:
    result = pd.Series("grey", index=universe, dtype=str)
    result.loc[fitted_events] = pd.Series(labels, index=fitted_events, dtype=str)
    return result


def _arm_summary(
    arm: str,
    labels: pd.Series,
    universe_size: int,
    bootstrap_rows: pd.DataFrame,
) -> dict[str, object]:
    non_grey = labels.loc[labels.ne("grey")]
    sizes = non_grey.value_counts()
    subset = bootstrap_rows.loc[bootstrap_rows["arm"].eq(arm), "adjusted_rand_index"]
    return {
        "arm": arm,
        "n_input_events": universe_size,
        "n_modules": int(len(sizes)),
        "n_assigned_events": int(len(non_grey)),
        "assigned_event_fraction": float(len(non_grey) / universe_size),
        "largest_module_fraction": float(sizes.max() / universe_size) if not sizes.empty else 0.0,
        "mean_bootstrap_ARI": float(subset.mean()),
        "median_bootstrap_ARI": float(subset.median()),
        "minimum_bootstrap_ARI": float(subset.min()),
    }


def _render_plot(
    summary: pd.DataFrame,
    bootstrap: pd.DataFrame,
    reliability: pd.Series,
    output: Path,
) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    palette = ["#6B7280", "#2D6CDF", "#0F766E"]

    sns.boxplot(
        data=bootstrap,
        x="arm",
        y="adjusted_rand_index",
        hue="arm",
        legend=False,
        palette=palette,
        ax=axes[0],
    )
    sns.stripplot(
        data=bootstrap,
        x="arm",
        y="adjusted_rand_index",
        color="#17212B",
        alpha=0.55,
        size=4,
        ax=axes[0],
    )
    axes[0].set_title("Matched-bootstrap module stability", fontweight="bold")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("ARI vs full-data assignment")
    axes[0].tick_params(axis="x", rotation=22)

    axes[1].bar(summary["arm"], summary["assigned_event_fraction"], color=palette)
    axes[1].set_ylim(0, 1.05)
    axes[1].set_title("Confident event coverage", fontweight="bold")
    axes[1].set_ylabel("Assigned event fraction")
    axes[1].tick_params(axis="x", rotation=22)
    for index, value in enumerate(summary["assigned_event_fraction"]):
        axes[1].text(index, value + 0.025, f"{value:.1%}", ha="center", fontsize=11)

    sns.histplot(reliability, bins=25, color="#2D6CDF", edgecolor="white", ax=axes[2])
    axes[2].axvline(0.15, color="#C44536", linestyle="--", linewidth=2, label="production cutoff")
    axes[2].set_title("Read/replicate reliability", fontweight="bold")
    axes[2].set_xlabel("network reliability")
    axes[2].legend(frameon=False)

    fig.suptitle(
        "Setd2 real-count three-arm ablation",
        fontsize=20,
        fontweight="bold",
        x=0.02,
        ha="left",
    )
    fig.tight_layout()
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Matched-bootstrap WGCNA/read-uncertainty/full-CoSplice ablation."
    )
    parser.add_argument("--cosplice-output", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--rscript", type=Path, required=True)
    parser.add_argument("--bootstraps", type=int, default=20)
    parser.add_argument("--seed", type=int, default=37)
    parser.add_argument("--min-module-size", type=int, default=20)
    parser.add_argument("--max-events", type=int, default=3_000)
    args = parser.parse_args(argv)

    source = args.cosplice_output.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    paths = {
        "raw": source / "03_stable_as_qc" / "stable_psi_matrix.tsv",
        "qc": source / "03_stable_as_qc" / "stable_as_event_qc.tsv",
        "ijc": source / "02_matrix" / "sample_level_ijc_matrix.tsv",
        "sjc": source / "02_matrix" / "sample_level_sjc_matrix.tsv",
        "samples": source / "02_matrix" / "sample_metadata_resolved.tsv",
    }
    missing = [str(path) for path in paths.values() if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing real-count ablation inputs: {missing}")

    raw = _load_matrix(paths["raw"])
    qc = pd.read_csv(paths["qc"], sep="\t", index_col=0, low_memory=False)
    qc.index = qc.index.astype(str)
    qc = qc.reindex(raw.index)
    if args.max_events > 0 and len(raw) > args.max_events:
        ranking = pd.to_numeric(qc.get("stability_score"), errors="coerce")
        if ranking.notna().sum() == 0:
            ranking = raw.var(axis=1, skipna=True)
        selected = ranking.sort_values(ascending=False).head(args.max_events).index
        raw = raw.reindex(selected)
        qc = qc.reindex(selected)
    ijc = _load_matrix(paths["ijc"], raw.index).reindex(columns=raw.columns)
    sjc = _load_matrix(paths["sjc"], raw.index).reindex(columns=raw.columns)
    samples = pd.read_csv(paths["samples"], sep="\t", low_memory=False)
    samples = samples.loc[samples["sample_id"].astype(str).isin(raw.columns)].copy()
    samples["sample_id"] = samples["sample_id"].astype(str)
    samples = samples.set_index("sample_id").reindex(raw.columns).reset_index()
    samples = samples.rename(columns={samples.columns[0]: "sample_id"})
    if samples["condition"].isna().any():
        raise ValueError("Setd2 sample metadata does not cover every PSI column.")

    uncertainty_cfg = UncertaintyConfig()
    uncertainty = estimate_beta_binomial_uncertainty(
        ijc,
        sjc,
        samples,
        qc,
        prior_strength=uncertainty_cfg.prior_strength,
        reliability_read_scale=uncertainty_cfg.reliability_read_scale,
    )
    reliability = uncertainty.event_metrics["network_reliability"]
    retained = reliability.loc[reliability.ge(uncertainty_cfg.min_network_reliability)].index
    if len(retained) < args.min_module_size:
        raise ValueError("Too few events passed the production uncertainty cutoff.")

    network = NetworkConfig(
        engine="tom_hierarchical",
        module_detection="wgcna_blockwise",
        min_module_size=args.min_module_size,
        max_events=None,
        tom_dense_max_events=max(3_000, len(raw)),
        wgcna_rscript=str(args.rscript.resolve()),
        consensus_bootstrap_iterations=0,
        community_consensus_runs=1,
        bootstrap_iterations=0,
        random_seed=args.seed,
    )
    consensus_cfg = ConsensusConfig(
        bootstrap_iterations=args.bootstraps,
        min_module_size=args.min_module_size,
        clustering_method="backbone_refined_consensus",
        dense_max_events=max(3_000, len(retained)),
        random_seed=args.seed,
    )

    raw_labels = discover_module_labels(raw, network)
    weighted_matrix = uncertainty.posterior_psi.reindex(retained)
    retained_reliability = reliability.reindex(retained)
    uncertainty_labels = discover_module_labels(
        weighted_matrix,
        network,
    )
    universe = raw.index
    references = {
        ARM_RAW: _expanded_labels(universe, universe, raw_labels),
        ARM_UNCERTAINTY: _expanded_labels(universe, retained, uncertainty_labels),
    }
    rng = np.random.default_rng(args.seed)
    bootstrap_rows: list[dict[str, object]] = []
    weighted_run_labels: list[np.ndarray] = []
    coassignment_counts = np.zeros((len(retained), len(retained)), dtype=np.uint16)
    for iteration in range(1, args.bootstraps + 1):
        indices = _bootstrap_indices(samples, raw.columns, rng)
        raw_boot = raw.iloc[:, indices].copy()
        raw_boot.columns = [f"boot_{iteration:03d}_{position:03d}" for position in range(len(indices))]
        weighted_boot = uncertainty.posterior_psi.reindex(retained).iloc[:, indices].copy()
        weighted_boot.columns = raw_boot.columns
        boot_raw_labels = discover_module_labels(raw_boot, network)
        boot_weighted_labels = discover_module_labels(
            weighted_boot,
            network,
        )
        weighted_run_labels.append(boot_weighted_labels.copy())
        for module in set(boot_weighted_labels).difference({"grey"}):
            module_indices = np.flatnonzero(boot_weighted_labels == module)
            coassignment_counts[np.ix_(module_indices, module_indices)] += 1
        boot_assignments = {
            ARM_RAW: _expanded_labels(universe, universe, boot_raw_labels),
            ARM_UNCERTAINTY: _expanded_labels(universe, retained, boot_weighted_labels),
        }
        for arm, labels in boot_assignments.items():
            bootstrap_rows.append(
                {
                    "bootstrap_iteration": iteration,
                    "arm": arm,
                    "adjusted_rand_index": adjusted_rand_score(references[arm], labels),
                }
            )
        print(f"Matched bootstrap {iteration}/{args.bootstraps}", flush=True)

    coassignment = coassignment_counts.astype(np.float32) / float(args.bootstraps)
    np.fill_diagonal(coassignment, 1.0)
    (
        full_label_values,
        membership_role,
        confidence,
        second_best,
        confidence_margin,
        candidate_labels,
        parameter_search,
    ) = _select_final_modules(
        coassignment,
        consensus_cfg,
        args.min_module_size,
        weighted_run_labels,
        backbone_labels=np.asarray(uncertainty_labels, dtype=object),
    )
    baseline = pd.Series(uncertainty_labels, index=retained, dtype=str)
    full_assignments, full_module_summary = _confidence_tables(
        coassignment,
        retained.astype(str).to_numpy(),
        full_label_values,
        baseline,
        membership_role,
        candidate_labels,
        confidence,
        second_best,
        confidence_margin,
    )
    references[ARM_FULL] = _expanded_labels(
        universe, retained, np.asarray(full_label_values, dtype=object)
    )
    for iteration, boot_weighted_labels in enumerate(weighted_run_labels, start=1):
        labels = _expanded_labels(universe, retained, boot_weighted_labels)
        bootstrap_rows.append(
            {
                "bootstrap_iteration": iteration,
                "arm": ARM_FULL,
                "adjusted_rand_index": adjusted_rand_score(references[ARM_FULL], labels),
            }
        )

    bootstrap = pd.DataFrame(bootstrap_rows)
    summary = pd.DataFrame(
        [
            _arm_summary(arm, references[arm], len(universe), bootstrap)
            for arm in (ARM_RAW, ARM_UNCERTAINTY, ARM_FULL)
        ]
    )
    assignments = pd.DataFrame(index=universe)
    assignments.index.name = "event_key"
    assignments["network_reliability"] = reliability
    assignments["passes_uncertainty_cutoff"] = assignments["network_reliability"].ge(
        uncertainty_cfg.min_network_reliability
    )
    for arm, labels in references.items():
        assignments[arm] = labels

    pair_rows = []
    arms = list(references)
    for left in range(len(arms) - 1):
        for right in range(left + 1, len(arms)):
            pair_rows.append(
                {
                    "arm_a": arms[left],
                    "arm_b": arms[right],
                    "ARI_all_starting_events": adjusted_rand_score(
                        references[arms[left]], references[arms[right]]
                    ),
                    "ARI_uncertainty_pass_events": adjusted_rand_score(
                        references[arms[left]].reindex(retained),
                        references[arms[right]].reindex(retained),
                    ),
                }
            )
    concordance = pd.DataFrame(pair_rows)

    summary.to_csv(output / "setd2_three_arm_summary.tsv", sep="\t", index=False)
    bootstrap.to_csv(output / "setd2_matched_bootstrap_stability.tsv", sep="\t", index=False)
    assignments.to_csv(output / "setd2_event_assignments.tsv", sep="\t")
    concordance.to_csv(output / "setd2_assignment_concordance.tsv", sep="\t", index=False)
    uncertainty.event_metrics.to_csv(output / "setd2_uncertainty_event_metrics.tsv", sep="\t")
    full_assignments.to_csv(output / "setd2_full_cosplice_confidence.tsv", sep="\t", index=False)
    full_module_summary.to_csv(output / "setd2_full_cosplice_module_summary.tsv", sep="\t", index=False)
    parameter_search.to_csv(output / "setd2_consensus_parameter_search.tsv", sep="\t", index=False)
    _render_plot(summary, bootstrap, reliability, output / "setd2_three_arm_ablation.png")

    parameters = {
        "design": "real replicate IJC/SJC three-arm ablation",
        "seed": args.seed,
        "bootstrap_iterations": args.bootstraps,
        "source": str(source),
        "inputs": {key: {"path": str(path.resolve()), "sha256": _sha256(path)} for key, path in paths.items()},
        "network": asdict(network),
        "consensus": asdict(consensus_cfg),
        "uncertainty": asdict(uncertainty_cfg),
        "n_starting_events": len(universe),
        "n_uncertainty_pass_events": len(retained),
        "max_events_predeclared": args.max_events,
    }
    (output / "parameters.json").write_text(
        json.dumps(parameters, indent=2, default=str), encoding="utf-8"
    )
    report = [
        "# Setd2 real-count three-arm ablation",
        "",
        "## Design",
        "",
        f"- Starting Stable-AS universe: {len(universe)} events across {raw.shape[1]} samples.",
        f"- Production read/replicate reliability cutoff: {uncertainty_cfg.min_network_reliability:.2f}; retained {len(retained)} events.",
        f"- All arms used seed {args.seed} and the same {args.bootstraps} condition-stratified bootstrap samples.",
        "- Raw arm: signed bicor WGCNA/TOM on raw PSI.",
        "- Uncertainty arm: the same WGCNA/TOM backbone on beta-binomial posterior PSI with low-reliability rejection; reliability is not multiplied into TOM adjacency.",
        "- Full CoSplice arm: uncertainty arm plus bootstrap co-assignment confidence and conservative boundary rejection.",
        "",
        "## Results",
        "",
        summary.to_markdown(index=False, floatfmt=".4f"),
        "",
        "## Assignment concordance",
        "",
        concordance.to_markdown(index=False, floatfmt=".4f"),
        "",
        "## Interpretation boundary",
        "",
        "This ablation tests whether read-count uncertainty and consensus confidence improve reproducibility or confidence on the same Setd2 experiment. It does not constitute external biological replication.",
    ]
    (output / "SETD2_THREE_ARM_ABLATION_REPORT.md").write_text("\n".join(report), encoding="utf-8")
    print(summary.to_string(index=False), flush=True)
    print(f"Output: {output}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
