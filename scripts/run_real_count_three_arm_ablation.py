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
from katana_cosplice.consensus import (
    _bootstrap_kme_values,
    _confidence_tables,
    _event_kme,
    _integrate_kme_membership_evidence,
    _kme_stability_table,
    _select_final_modules,
    _split_half_membership_confirmation,
)
from katana_cosplice.ensemble_network import posterior_ensemble_correlation
from katana_cosplice.network import discover_module_labels
from katana_cosplice.network import _program_scores
from katana_cosplice.membership_validation import (
    add_continuous_membership_confidence,
    build_membership_arms,
)
from katana_cosplice.r_wgcna import exact_wgcna_labels
from katana_cosplice.uncertainty import estimate_beta_binomial_uncertainty


ARM_RAW = "raw_PSI_WGCNA"
ARM_UNCERTAINTY = "posterior_mean_WGCNA"
ARM_ENSEMBLE = "posterior_ensemble_WGCNA"
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
    bootstrap_scope = {
        ARM_RAW: "condition-stratified sample bootstrap",
        ARM_UNCERTAINTY: "condition-stratified sample bootstrap",
        ARM_ENSEMBLE: "sample bootstrap plus posterior correlation ensemble",
        ARM_FULL: "sample bootstrap plus one beta-posterior PSI draw",
    }[arm]
    return {
        "arm": arm,
        "bootstrap_scope": bootstrap_scope,
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
    fig, axes = plt.subplots(1, 3, figsize=(20, 6.5))
    palette = ["#6B7280", "#4C78A8", "#E07A5F", "#0F766E"]
    display_labels = {
        ARM_RAW: "Raw PSI\nWGCNA",
        ARM_UNCERTAINTY: "Posterior\nmean WGCNA",
        ARM_ENSEMBLE: "Posterior\nensemble\n(experimental)",
        ARM_FULL: "Full CoSplice\nconfidence\nstress test",
    }

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
    axes[0].set_title("Partition recovery / uncertainty stress", fontweight="bold")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("ARI vs full-data assignment")
    axes[0].set_xticks(
        axes[0].get_xticks(),
        [display_labels.get(value, value) for value in summary["arm"]],
    )
    axes[0].tick_params(axis="x", rotation=0, labelsize=10)

    axes[1].bar(summary["arm"], summary["assigned_event_fraction"], color=palette)
    axes[1].set_ylim(0, 1.05)
    axes[1].set_title("Fixed-partition event assignment", fontweight="bold")
    axes[1].set_ylabel("Assigned event fraction")
    axes[1].set_xticks(
        axes[1].get_xticks(),
        [display_labels.get(value, value) for value in summary["arm"]],
    )
    axes[1].tick_params(axis="x", rotation=0, labelsize=10)
    for index, value in enumerate(summary["assigned_event_fraction"]):
        axes[1].text(index, value + 0.025, f"{value:.1%}", ha="center", fontsize=11)

    sns.histplot(reliability, bins=25, color="#2D6CDF", edgecolor="white", ax=axes[2])
    axes[2].axvline(0.15, color="#C44536", linestyle="--", linewidth=2, label="production cutoff")
    axes[2].set_title("Read/replicate reliability", fontweight="bold")
    axes[2].set_xlabel("network reliability")
    axes[2].legend(frameon=False)

    fig.suptitle(
        "Setd2 real-count four-arm ablation",
        fontsize=20,
        fontweight="bold",
        x=0.02,
        ha="left",
    )
    fig.text(
        0.02,
        0.01,
        "Raw/posterior mean: sample-bootstrap recovery. Ensemble/full CoSplice: "
        "additional posterior-measurement uncertainty stress. High assignment alone is not superiority.",
        fontsize=10,
        color="#4B5563",
    )
    fig.tight_layout(rect=[0, 0.06, 1, 0.94])
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _render_membership_confidence(
    assignments: pd.DataFrame,
    output: Path,
) -> None:
    visible = assignments.loc[assignments["consensus_module"].ne("grey")].copy()
    module_order = sorted(visible["consensus_module"].astype(str).unique())
    role_order = ["core", "peripheral", "uncertain"]
    role_colors = {
        "core": "#0072B2",
        "peripheral": "#56B4E9",
        "uncertain": "#B8B8B8",
    }
    counts = (
        visible.groupby(["consensus_module", "membership_role"], observed=True)
        .size()
        .unstack(fill_value=0)
        .reindex(index=module_order, columns=role_order, fill_value=0)
    )
    sns.set_theme(style="ticks", context="talk")
    fig, axes = plt.subplots(1, 3, figsize=(20, 6.2))
    bottom = np.zeros(len(counts), dtype=float)
    for role in role_order:
        values = counts[role].to_numpy(dtype=float)
        axes[0].bar(
            counts.index,
            values,
            bottom=bottom,
            color=role_colors[role],
            label=role,
        )
        bottom += values
    axes[0].set_title("Fixed WGCNA module membership", loc="left", fontweight="bold")
    axes[0].set_ylabel("AS events")
    axes[0].set_xlabel("")
    axes[0].legend(frameon=False)

    for role in role_order:
        subset = visible.loc[visible["membership_role"].eq(role)]
        axes[1].scatter(
            subset["consensus_membership_confidence"],
            subset["bootstrap_mean_abs_kME"],
            s=20,
            alpha=0.65,
            color=role_colors[role],
            edgecolors="none",
            label=role,
        )
    axes[1].axvline(0.70, color="#333333", linestyle="--", linewidth=1)
    axes[1].axhline(0.70, color="#333333", linestyle="--", linewidth=1)
    axes[1].set_xlim(-0.02, 1.02)
    axes[1].set_ylim(-0.02, 1.02)
    axes[1].set_title("Co-assignment vs bootstrap kME", loc="left", fontweight="bold")
    axes[1].set_xlabel("Bootstrap co-assignment to own module")
    axes[1].set_ylabel("Mean absolute bootstrap kME")
    axes[1].legend(frameon=False, loc="lower right")

    sns.boxplot(
        data=visible,
        x="membership_role",
        y="continuous_membership_confidence",
        order=role_order,
        palette=role_colors,
        hue="membership_role",
        legend=False,
        ax=axes[2],
    )
    sns.stripplot(
        data=visible,
        x="membership_role",
        y="continuous_membership_confidence",
        order=role_order,
        color="#111827",
        alpha=0.35,
        size=2.5,
        ax=axes[2],
    )
    axes[2].set_ylim(-0.02, 1.02)
    axes[2].set_title("Continuous four-domain confidence", loc="left", fontweight="bold")
    axes[2].set_xlabel("")
    axes[2].set_ylabel("Continuous membership confidence")
    sns.despine(fig=fig)
    fig.suptitle(
        "Setd2 CoSplice membership confidence",
        x=0.02,
        ha="left",
        fontsize=20,
        fontweight="bold",
    )
    fig.text(
        0.02,
        0.01,
        "Core membership requires confirmation in two disjoint bootstrap halves. "
        "The WGCNA partition is fixed; confidence labels do not move events between modules.",
        fontsize=10,
        color="#4B5563",
    )
    fig.tight_layout(rect=[0, 0.06, 1, 0.93])
    fig.savefig(output, dpi=240, bbox_inches="tight")
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Matched-bootstrap WGCNA/posterior-ensemble/full-CoSplice ablation."
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
        clustering_method="backbone_confidence_only",
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
    ensemble = posterior_ensemble_correlation(
        uncertainty.posterior_alpha.reindex(retained),
        uncertainty.posterior_beta.reindex(retained),
        n_draws=uncertainty_cfg.posterior_draws,
        method=network.method,
        random_seed=uncertainty_cfg.posterior_random_seed,
    )
    ensemble_labels = exact_wgcna_labels(
        weighted_matrix,
        soft_power=network.soft_power,
        min_module_size=network.min_module_size,
        merge_cut_height=network.wgcna_merge_cut_height,
        deep_split=network.dynamic_deep_split,
        pam_respects_dendro=network.dynamic_pam_respects_dendro,
        rscript=network.wgcna_rscript,
        precomputed_correlation=ensemble.mean_correlation,
    )
    universe = raw.index
    references = {
        ARM_RAW: _expanded_labels(universe, universe, raw_labels),
        ARM_UNCERTAINTY: _expanded_labels(universe, retained, uncertainty_labels),
        ARM_ENSEMBLE: _expanded_labels(universe, retained, ensemble_labels),
    }
    rng = np.random.default_rng(args.seed)
    bootstrap_rows: list[dict[str, object]] = []
    ensemble_run_labels: list[np.ndarray] = []
    full_run_labels: list[np.ndarray] = []
    full_bootstrap_kme: list[np.ndarray] = []
    full_reference_scores = _program_scores(
        weighted_matrix,
        uncertainty_labels,
    )
    full_baseline_kme = _event_kme(
        weighted_matrix,
        uncertainty_labels,
        full_reference_scores,
    )
    coassignment_counts = np.zeros((len(retained), len(retained)), dtype=np.uint16)
    for iteration in range(1, args.bootstraps + 1):
        indices = _bootstrap_indices(samples, raw.columns, rng)
        raw_boot = raw.iloc[:, indices].copy()
        raw_boot.columns = [f"boot_{iteration:03d}_{position:03d}" for position in range(len(indices))]
        weighted_boot = uncertainty.posterior_psi.reindex(retained).iloc[:, indices].copy()
        weighted_boot.columns = raw_boot.columns
        alpha_boot = uncertainty.posterior_alpha.reindex(retained).iloc[:, indices].copy()
        beta_boot = uncertainty.posterior_beta.reindex(retained).iloc[:, indices].copy()
        alpha_boot.columns = raw_boot.columns
        beta_boot.columns = raw_boot.columns
        boot_raw_labels = discover_module_labels(raw_boot, network)
        boot_weighted_labels = discover_module_labels(
            weighted_boot,
            network,
        )
        boot_ensemble = posterior_ensemble_correlation(
            alpha_boot,
            beta_boot,
            n_draws=uncertainty_cfg.posterior_draws,
            method=network.method,
            random_seed=uncertainty_cfg.posterior_random_seed + iteration,
        )
        boot_ensemble_labels = exact_wgcna_labels(
            weighted_boot,
            soft_power=network.soft_power,
            min_module_size=network.min_module_size,
            merge_cut_height=network.wgcna_merge_cut_height,
            deep_split=network.dynamic_deep_split,
            pam_respects_dendro=network.dynamic_pam_respects_dendro,
            rscript=network.wgcna_rscript,
            precomputed_correlation=boot_ensemble.mean_correlation,
        )
        ensemble_run_labels.append(boot_ensemble_labels.copy())
        posterior_draw_boot = pd.DataFrame(
            rng.beta(
                alpha_boot.to_numpy(dtype=float),
                beta_boot.to_numpy(dtype=float),
            ),
            index=retained,
            columns=raw_boot.columns,
        )
        boot_full_labels = discover_module_labels(posterior_draw_boot, network)
        full_run_labels.append(boot_full_labels.copy())
        full_bootstrap_kme.append(
            _bootstrap_kme_values(
                posterior_draw_boot,
                uncertainty_labels,
                full_reference_scores,
                indices,
            )
        )
        for module in set(boot_full_labels).difference({"grey"}):
            module_indices = np.flatnonzero(boot_full_labels == module)
            coassignment_counts[np.ix_(module_indices, module_indices)] += 1
        boot_assignments = {
            ARM_RAW: _expanded_labels(universe, universe, boot_raw_labels),
            ARM_UNCERTAINTY: _expanded_labels(universe, retained, boot_weighted_labels),
            ARM_ENSEMBLE: _expanded_labels(universe, retained, boot_ensemble_labels),
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
        full_run_labels,
        backbone_labels=np.asarray(uncertainty_labels, dtype=object),
    )
    full_kme_stability = _kme_stability_table(
        retained.astype(str).to_numpy(),
        full_baseline_kme,
        full_bootstrap_kme,
    )
    membership_role, membership_evidence = _integrate_kme_membership_evidence(
        full_label_values,
        membership_role,
        confidence_margin,
        full_kme_stability,
        consensus_cfg,
    )
    (
        membership_role,
        membership_evidence,
        split_half_role_a,
        split_half_role_b,
    ) = _split_half_membership_confirmation(
        np.asarray(uncertainty_labels, dtype=object),
        membership_role,
        membership_evidence,
        full_run_labels,
        full_baseline_kme,
        full_bootstrap_kme,
        retained.astype(str).to_numpy(),
        args.min_module_size,
        consensus_cfg,
    )
    full_kme_stability["split_half_membership_role_a"] = split_half_role_a
    full_kme_stability["split_half_membership_role_b"] = split_half_role_b
    full_kme_stability["split_half_core_confirmed"] = (
        (split_half_role_a == "core") & (split_half_role_b == "core")
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
        membership_evidence=membership_evidence,
        kme_stability=full_kme_stability,
    )
    full_assignments = add_continuous_membership_confidence(
        full_assignments,
        event_reliability=retained_reliability,
    )
    membership_arms, membership_count_audit = build_membership_arms(
        full_assignments
    )
    membership_selection_rows: list[dict[str, object]] = []
    indexed_confidence = full_assignments.set_index("event_key")
    for membership_arm, selected in membership_arms.items():
        selected_evidence = indexed_confidence.reindex(
            selected["event_key"].astype(str)
        )
        membership_selection_rows.append(
            {
                "membership_arm": membership_arm,
                "selected_events": len(selected),
                "median_continuous_confidence": pd.to_numeric(
                    selected_evidence["continuous_membership_confidence"],
                    errors="coerce",
                ).median(),
                "median_read_reliability": pd.to_numeric(
                    selected_evidence["confidence_read_component"],
                    errors="coerce",
                ).median(),
                "median_coassignment_component": pd.to_numeric(
                    selected_evidence["confidence_coassignment_component"],
                    errors="coerce",
                ).median(),
                "median_kme_component": pd.to_numeric(
                    selected_evidence["confidence_kme_component"],
                    errors="coerce",
                ).median(),
                "split_half_core_fraction": selected_evidence[
                    "membership_role"
                ].astype(str).eq("core").mean(),
            }
        )
    membership_selection = pd.DataFrame(membership_selection_rows)
    references[ARM_FULL] = _expanded_labels(
        universe, retained, np.asarray(full_label_values, dtype=object)
    )
    for iteration, boot_full_labels in enumerate(full_run_labels, start=1):
        labels = _expanded_labels(universe, retained, boot_full_labels)
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
            for arm in (ARM_RAW, ARM_UNCERTAINTY, ARM_ENSEMBLE, ARM_FULL)
        ]
    )
    role_counts = full_assignments["membership_role"].value_counts()
    full_mask = summary["arm"].eq(ARM_FULL)
    summary["n_core_events"] = np.nan
    summary["n_peripheral_events"] = np.nan
    summary["n_uncertain_assigned_events"] = np.nan
    summary["core_fraction_of_assigned"] = np.nan
    summary.loc[full_mask, "n_core_events"] = int(role_counts.get("core", 0))
    summary.loc[full_mask, "n_peripheral_events"] = int(
        role_counts.get("peripheral", 0)
    )
    summary.loc[full_mask, "n_uncertain_assigned_events"] = int(
        role_counts.get("uncertain", 0)
    )
    summary.loc[full_mask, "core_fraction_of_assigned"] = float(
        role_counts.get("core", 0) / max(np.count_nonzero(full_label_values != "grey"), 1)
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
    decision = summary.loc[
        :,
        [
            "arm",
            "n_modules",
            "assigned_event_fraction",
            "largest_module_fraction",
            "mean_bootstrap_ARI",
        ],
    ].copy()
    decision["production_role"] = decision["arm"].map(
        {
            ARM_RAW: "sensitivity comparator",
            ARM_UNCERTAINTY: "production partition",
            ARM_ENSEMBLE: "rejected experimental partition",
            ARM_FULL: "production confidence stress test",
        }
    )
    decision["decision_reason"] = decision["arm"].map(
        {
            ARM_RAW: "Established raw-PSI WGCNA reference arm.",
            ARM_UNCERTAINTY: "Best balance of module recovery, coverage and largest-module fraction.",
            ARM_ENSEMBLE: "High ARI accompanied by an over-dominant module and near-universal assignment.",
            ARM_FULL: "Does not redefine membership; quantifies sample plus read-measurement uncertainty.",
        }
    )

    summary.to_csv(output / "setd2_four_arm_summary.tsv", sep="\t", index=False)
    bootstrap.to_csv(output / "setd2_matched_bootstrap_stability.tsv", sep="\t", index=False)
    assignments.to_csv(output / "setd2_event_assignments.tsv", sep="\t")
    concordance.to_csv(output / "setd2_assignment_concordance.tsv", sep="\t", index=False)
    decision.to_csv(output / "setd2_arm_decision.tsv", sep="\t", index=False)
    uncertainty.event_metrics.to_csv(output / "setd2_uncertainty_event_metrics.tsv", sep="\t")
    full_assignments.to_csv(output / "setd2_full_cosplice_confidence.tsv", sep="\t", index=False)
    full_module_summary.to_csv(output / "setd2_full_cosplice_module_summary.tsv", sep="\t", index=False)
    membership_count_audit.to_csv(
        output / "setd2_membership_arm_count_matching.tsv",
        sep="\t",
        index=False,
    )
    membership_selection.to_csv(
        output / "setd2_membership_selection_audit.tsv",
        sep="\t",
        index=False,
    )
    for membership_arm, selected in membership_arms.items():
        selected.to_csv(
            output / f"setd2_{membership_arm}_assignments.tsv",
            sep="\t",
            index=False,
        )
    evidence_summary = (
        full_assignments.groupby(
            ["consensus_module", "membership_role", "membership_evidence"],
            observed=True,
        )
        .size()
        .rename("n_events")
        .reset_index()
    )
    evidence_summary.to_csv(
        output / "setd2_membership_confidence_summary.tsv",
        sep="\t",
        index=False,
    )
    parameter_search.to_csv(output / "setd2_consensus_parameter_search.tsv", sep="\t", index=False)
    _render_plot(summary, bootstrap, reliability, output / "setd2_four_arm_ablation.png")
    _render_membership_confidence(
        full_assignments,
        output / "setd2_membership_confidence.png",
    )

    parameters = {
        "design": "real replicate IJC/SJC four-arm ablation",
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
        "membership_confidence": {
            "partition": "fixed posterior-mean WGCNA",
            "evidence": "bootstrap co-assignment OR fixed-module bootstrap kME",
            "confirmation": "both disjoint bootstrap halves",
            "continuous_score": "equal-domain geometric mean of co-assignment, kME, split-half confirmation and read reliability",
            "thresholds": {
                key: value
                for key, value in asdict(consensus_cfg).items()
                if "membership" in key or "kme" in key
            },
        },
    }
    (output / "parameters.json").write_text(
        json.dumps(parameters, indent=2, default=str), encoding="utf-8"
    )
    report = [
        "# Setd2 real-count four-arm ablation",
        "",
        "## Design",
        "",
        f"- Starting Stable-AS universe: {len(universe)} events across {raw.shape[1]} samples.",
        f"- Production read/replicate reliability cutoff: {uncertainty_cfg.min_network_reliability:.2f}; retained {len(retained)} events.",
        f"- All arms used seed {args.seed} and the same {args.bootstraps} condition-stratified bootstrap samples.",
        "- Raw arm: signed bicor WGCNA/TOM on raw PSI.",
        "- Uncertainty arm: the same WGCNA/TOM backbone on beta-binomial posterior PSI with low-reliability rejection; reliability is not multiplied into TOM adjacency.",
        "- Posterior-ensemble arm: Fisher-z mean correlations across beta-posterior PSI draws feed the same signed WGCNA/TOM backbone.",
        "- Full CoSplice arm: posterior-mean WGCNA backbone plus condition-stratified beta-posterior-draw co-assignment confidence; confidence does not alter membership.",
        "",
        "## Results",
        "",
        summary.to_markdown(index=False, floatfmt=".4f"),
        "",
        "Raw and posterior-mean ARI compare sample-bootstrap partition recovery. "
        "The ensemble and full-CoSplice ARI columns additionally perturb PSI measurement "
        "uncertainty and are therefore stricter uncertainty stress tests, not directly "
        "interchangeable performance scores.",
        "",
        "## Assignment concordance",
        "",
        concordance.to_markdown(index=False, floatfmt=".4f"),
        "",
        "## Production decision",
        "",
        decision.to_markdown(index=False, floatfmt=".4f"),
        "",
        "## Membership confidence",
        "",
        f"- Split-half confirmed core events: {int(role_counts.get('core', 0))} "
        f"of {int(np.count_nonzero(full_label_values != 'grey'))} fixed-partition events "
        f"({float(role_counts.get('core', 0) / max(np.count_nonzero(full_label_values != 'grey'), 1)):.1%}).",
        f"- Split-half supported peripheral events: {int(role_counts.get('peripheral', 0))}.",
        f"- Assigned but unconfirmed events: {int(role_counts.get('uncertain', 0))}.",
        "- A core label requires either strong label-invariant co-assignment or strong "
        "fixed-module kME, plus confirmation in both disjoint bootstrap halves. "
        "This confidence layer never changes the WGCNA module assignment.",
        "",
        "### Count-matched member-selection audit",
        "",
        membership_selection.to_markdown(index=False, floatfmt=".4f"),
        "",
        "## Interpretation boundary",
        "",
        "This ablation tests whether read-count uncertainty and consensus confidence improve reproducibility or confidence on the same Setd2 experiment. It does not constitute external biological replication.",
    ]
    (output / "SETD2_FOUR_ARM_ABLATION_REPORT.md").write_text("\n".join(report), encoding="utf-8")
    print(summary.to_string(index=False), flush=True)
    print(f"Output: {output}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
