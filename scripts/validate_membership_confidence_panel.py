from __future__ import annotations

import argparse
import json
import math
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import hypergeom, ttest_ind, wilcoxon

from katana_cosplice.benchmark import _fit_method, _rscript_with_wgcna
from katana_cosplice.biological_interpretability import gene_from_event_key
from katana_cosplice.config import (
    BenchmarkConfig,
    ConsensusConfig,
    ExternalCohortConfig,
    NetworkConfig,
    PreservationConfig,
)
from katana_cosplice.consensus import build_bootstrap_coassignment_consensus
from katana_cosplice.formal_preservation import run_formal_preservation
from katana_cosplice.locked_ablation import (
    combine_gmt_files,
    file_record,
    load_gmt,
    read_psi_matrix,
    sample_frame,
    select_discovery_events,
    utc_now,
    write_json,
)
from katana_cosplice.membership_validation import (
    ARM_CONFIDENCE,
    ARM_CORE,
    ARM_FULL,
    ARM_MATCHED,
    ARMS,
    build_membership_arms,
)
from katana_cosplice.preservation import evaluate_external_module_preservation
from katana_cosplice.utils import bh_adjust


PALETTE = {
    ARM_FULL: "#9CA3AF",
    ARM_MATCHED: "#D97706",
    ARM_CORE: "#087E8B",
    ARM_CONFIDENCE: "#7C3AED",
}
ARM_DISPLAY = {
    ARM_FULL: "All WGCNA members",
    ARM_MATCHED: "Top-kME matched",
    ARM_CORE: "Split-half core",
    ARM_CONFIDENCE: "Continuous confidence",
}


def cohort_paths(panel: Path, cohort: str, event_type: str) -> dict[str, Path]:
    prepared = panel / cohort / "prepared"
    return {
        "psi_full": prepared / f"{cohort}_{event_type}_PSI_matrix.tsv",
        "psi_tumor": prepared / f"{cohort}_{event_type}_PSI_tumor_only.tsv",
        "events": prepared / f"{cohort}_{event_type}_event_metadata.tsv",
        "samples": prepared / f"{cohort}_sample_metadata.tsv",
    }


def assignment_frame(event_keys: pd.Index, labels: np.ndarray, column: str = "module") -> pd.DataFrame:
    return pd.DataFrame({"event_key": event_keys.astype(str), column: labels.astype(str)})


def _event_enrichment(
    cohort: str,
    arm: str,
    assignments: pd.DataFrame,
    gene_sets: dict[str, set[str]],
) -> pd.DataFrame:
    frame = assignments.copy()
    frame["gene"] = frame["event_key"].map(gene_from_event_key).str.casefold()
    frame = frame.loc[frame["consensus_module"].astype(str).ne("grey") & frame["gene"].ne("")]
    universe = set(frame["gene"])
    rows: list[dict[str, object]] = []
    for module, group in frame.groupby("consensus_module", observed=True, sort=True):
        genes = set(group["gene"])
        for term, raw_members in gene_sets.items():
            members = raw_members & universe
            overlap = genes & members
            pvalue = (
                float(hypergeom.sf(len(overlap) - 1, len(universe), len(members), len(genes)))
                if overlap and universe and members
                else 1.0
            )
            rows.append(
                {
                    "cohort": cohort,
                    "arm": arm,
                    "module": str(module),
                    "term": term,
                    "module_genes": len(genes),
                    "overlap_count": len(overlap),
                    "overlap_genes": ";".join(sorted(overlap)),
                    "pvalue": pvalue,
                }
            )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["FDR"] = result.groupby(["cohort", "arm", "module"], observed=True)["pvalue"].transform(bh_adjust)
    return result


def functional_replication_by_overlap(
    cohort: str,
    arms: dict[str, pd.DataFrame],
    holdout_assignment: pd.DataFrame,
    gene_sets: dict[str, set[str]],
    fdr: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    enrichment_parts: list[pd.DataFrame] = []
    rows: list[dict[str, object]] = []
    holdout = holdout_assignment.rename(columns={"module": "consensus_module"}).copy()
    holdout_sets = {
        str(module): set(group["event_key"].astype(str))
        for module, group in holdout.loc[holdout["consensus_module"].astype(str).ne("grey")].groupby(
            "consensus_module", observed=True
        )
    }
    for arm, discovery in arms.items():
        discovery_enrichment = _event_enrichment("discovery", arm, discovery, gene_sets)
        holdout_enrichment = _event_enrichment(cohort, arm, holdout, gene_sets)
        enrichment_parts.extend([discovery_enrichment, holdout_enrichment])
        discovery_terms = {
            str(module): set(group.loc[group["FDR"].le(fdr), "term"].astype(str))
            for module, group in discovery_enrichment.groupby("module", observed=True)
        }
        holdout_terms = {
            str(module): set(group.loc[group["FDR"].le(fdr), "term"].astype(str))
            for module, group in holdout_enrichment.groupby("module", observed=True)
        }
        for module, group in discovery.groupby("consensus_module", observed=True, sort=True):
            source_events = set(group["event_key"].astype(str))
            best_module = ""
            best_event_jaccard = 0.0
            for target_module, target_events in holdout_sets.items():
                union = source_events | target_events
                score = len(source_events & target_events) / len(union) if union else 0.0
                if score > best_event_jaccard:
                    best_module = target_module
                    best_event_jaccard = score
            source_terms = discovery_terms.get(str(module), set())
            target_terms = holdout_terms.get(best_module, set())
            shared_terms = source_terms & target_terms
            term_union = source_terms | target_terms
            rows.append(
                {
                    "cohort": cohort,
                    "arm": arm,
                    "discovery_module": str(module),
                    "best_holdout_module": best_module,
                    "event_jaccard": best_event_jaccard,
                    "discovery_significant_terms": len(source_terms),
                    "holdout_significant_terms": len(target_terms),
                    "shared_significant_terms": len(shared_terms),
                    "term_jaccard": len(shared_terms) / len(term_union) if term_union else 0.0,
                    "shared_terms": ";".join(sorted(shared_terms)),
                }
            )
    enrichment = pd.concat(enrichment_parts, ignore_index=True) if enrichment_parts else pd.DataFrame()
    return enrichment, pd.DataFrame(rows)


def _projected_program_scores(
    discovery: pd.DataFrame,
    target: pd.DataFrame,
) -> np.ndarray:
    reference = discovery.to_numpy(dtype=float)
    means = reference.mean(axis=1, keepdims=True)
    scales = reference.std(axis=1, ddof=1, keepdims=True)
    scales = np.where(scales < 1e-8, 1.0, scales)
    reference_z = (reference - means) / scales
    u, _s, _vt = np.linalg.svd(reference_z, full_matrices=False)
    weights = u[:, 0]
    target_values = target.to_numpy(dtype=float)
    target_values = np.where(np.isfinite(target_values), target_values, means)
    target_z = (target_values - means) / scales
    scores = weights @ target_z
    if np.corrcoef(scores, np.nanmean(target_values, axis=0))[0, 1] < 0:
        scores = -scores
    return scores


def phenotype_association(
    cohort: str,
    discovery: pd.DataFrame,
    target_full: pd.DataFrame,
    metadata_path: Path,
    arms: dict[str, pd.DataFrame],
    min_events: int,
) -> pd.DataFrame:
    metadata = pd.read_csv(metadata_path, sep="\t").set_index("sample_id")
    metadata = metadata.reindex(target_full.columns)
    phenotype = metadata["phenotype"].astype(str).str.casefold()
    rows: list[dict[str, object]] = []
    for arm, assignment in arms.items():
        for module, group in assignment.groupby("consensus_module", observed=True, sort=True):
            keys = [
                key for key in group["event_key"].astype(str)
                if key in discovery.index and key in target_full.index
            ]
            quality = [
                key for key in keys
                if target_full.loc[key].notna().sum() >= 6
                and float(target_full.loc[key].std(skipna=True)) > 1e-8
            ]
            if len(quality) < min_events:
                continue
            scores = _projected_program_scores(discovery.loc[quality], target_full.loc[quality])
            tumor = scores[phenotype.eq("tumor").to_numpy()]
            normal = scores[phenotype.eq("normal").to_numpy()]
            if len(tumor) < 3 or len(normal) < 3:
                continue
            statistic, pvalue = ttest_ind(tumor, normal, equal_var=False, nan_policy="omit")
            pooled = math.sqrt(
                ((len(tumor) - 1) * np.var(tumor, ddof=1) + (len(normal) - 1) * np.var(normal, ddof=1))
                / max(len(tumor) + len(normal) - 2, 1)
            )
            rows.append(
                {
                    "cohort": cohort,
                    "arm": arm,
                    "module": str(module),
                    "shared_events": len(quality),
                    "n_tumor": len(tumor),
                    "n_normal": len(normal),
                    "mean_tumor_score": float(np.mean(tumor)),
                    "mean_normal_score": float(np.mean(normal)),
                    "tumor_minus_normal": float(np.mean(tumor) - np.mean(normal)),
                    "cohens_d": float((np.mean(tumor) - np.mean(normal)) / pooled) if pooled > 0 else 0.0,
                    "welch_t": float(statistic),
                    "pvalue": float(pvalue),
                }
            )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["FDR"] = result.groupby(["cohort", "arm"], observed=True)["pvalue"].transform(bh_adjust)
    return result


def paired_confidence_test(table: pd.DataFrame, metric: str) -> dict[str, object]:
    subset = table.loc[table["arm"].isin([ARM_MATCHED, ARM_CONFIDENCE])]
    pivot = subset.pivot_table(index=["cohort", "consensus_module"], columns="arm", values=metric)
    pivot = pivot.dropna(subset=[ARM_MATCHED, ARM_CONFIDENCE])
    if pivot.empty:
        return {
            "metric": metric,
            "n_pairs": 0,
            "median_confidence_minus_matched": np.nan,
            "mean_confidence_minus_matched": np.nan,
            "confidence_wins": 0,
            "ties": 0,
            "matched_wins": 0,
            "wilcoxon_pvalue": np.nan,
        }
    delta = pivot[ARM_CONFIDENCE] - pivot[ARM_MATCHED]
    try:
        pvalue = float(wilcoxon(delta, zero_method="zsplit").pvalue)
    except ValueError:
        pvalue = 1.0
    return {
        "metric": metric,
        "n_pairs": len(delta),
        "median_confidence_minus_matched": float(delta.median()),
        "mean_confidence_minus_matched": float(delta.mean()),
        "confidence_wins": int((delta > 0).sum()),
        "ties": int((delta == 0).sum()),
        "matched_wins": int((delta < 0).sum()),
        "wilcoxon_pvalue": pvalue,
    }


def render_plots(
    output: Path,
    preservation: pd.DataFrame,
    functional: pd.DataFrame,
    phenotype: pd.DataFrame,
) -> list[Path]:
    sns.set_theme(style="ticks", context="paper", font="Arial")
    plot_dir = output / "plots"
    plot_dir.mkdir(exist_ok=True)
    generated: list[Path] = []

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2))
    sns.boxplot(data=preservation, x="arm", y="Zsummary_wgcna_like", hue="arm", palette=PALETTE, legend=False, ax=axes[0])
    sns.stripplot(data=preservation, x="arm", y="Zsummary_wgcna_like", color="#111827", alpha=0.55, size=3, ax=axes[0])
    axes[0].axhline(2.0, color="#6B7280", linestyle="--", linewidth=1)
    axes[0].set_title("External module preservation")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("WGCNA-like Zsummary")
    axes[0].tick_params(axis="x", rotation=20)
    sns.boxplot(data=preservation, x="arm", y="program_score_preservation", hue="arm", palette=PALETTE, legend=False, ax=axes[1])
    sns.stripplot(data=preservation, x="arm", y="program_score_preservation", color="#111827", alpha=0.55, size=3, ax=axes[1])
    axes[1].set_title("Program-score preservation")
    axes[1].set_xlabel("")
    axes[1].set_ylabel("Membership x variance preservation")
    axes[1].tick_params(axis="x", rotation=20)
    sns.despine(fig=fig)
    fig.tight_layout()
    path = plot_dir / "external_preservation_four_arm.png"
    fig.savefig(path, dpi=240, bbox_inches="tight")
    plt.close(fig)
    generated.append(path)

    functional_summary = functional.groupby(["cohort", "arm"], observed=True).agg(
        median_event_jaccard=("event_jaccard", "median"),
    ).reset_index()
    phenotype_summary = phenotype.assign(significant=phenotype["FDR"].le(0.05)).groupby(
        ["cohort", "arm"], observed=True
    ).agg(significant_modules=("significant", "sum"), median_abs_effect=("cohens_d", lambda x: float(np.median(np.abs(x))))).reset_index()
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.0))
    sns.barplot(data=functional_summary, x="cohort", y="median_event_jaccard", hue="arm", palette=PALETTE, ax=axes[0])
    axes[0].set_title("Independent holdout-module recovery")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("Median best event-overlap Jaccard")
    sns.barplot(data=phenotype_summary, x="cohort", y="median_abs_effect", hue="arm", palette=PALETTE, ax=axes[1])
    axes[1].set_title("Projected tumor-normal effect size")
    axes[1].set_xlabel("")
    axes[1].set_ylabel("Median absolute Cohen's d")
    for axis in axes:
        axis.tick_params(axis="x", rotation=30)
        axis.legend(title="", frameon=False, fontsize=8)
    sns.despine(fig=fig)
    fig.tight_layout()
    path = plot_dir / "functional_and_phenotype_replication.png"
    fig.savefig(path, dpi=240, bbox_inches="tight")
    plt.close(fig)
    generated.append(path)
    return generated


def render_formal_preservation_plots(
    output: Path,
    wgcna: pd.DataFrame,
    netrep: pd.DataFrame,
) -> list[Path]:
    plot_dir = output / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    wgcna = wgcna.assign(arm_display=wgcna["arm"].map(ARM_DISPLAY))
    netrep = netrep.assign(arm_display=netrep["arm"].map(ARM_DISPLAY))
    display_palette = {
        ARM_DISPLAY[arm]: color for arm, color in PALETTE.items()
    }
    sns.set_theme(style="ticks", context="talk")
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.6))
    sns.boxplot(
        data=wgcna,
        x="arm_display",
        y="Zsummary.pres",
        hue="arm_display",
        palette=display_palette,
        legend=False,
        ax=axes[0],
    )
    sns.stripplot(
        data=wgcna,
        x="arm_display",
        y="Zsummary.pres",
        color="#111827",
        alpha=0.55,
        size=3,
        ax=axes[0],
    )
    axes[0].axhline(2, color="#6B7280", linestyle="--", linewidth=1)
    axes[0].axhline(10, color="#111827", linestyle="--", linewidth=1)
    axes[0].set_title("Official WGCNA preservation", loc="left", fontweight="bold")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("WGCNA Zsummary.pres")

    sns.boxplot(
        data=netrep,
        x="arm_display",
        y="significant_statistics_FDR05",
        hue="arm_display",
        palette=display_palette,
        legend=False,
        ax=axes[1],
    )
    sns.stripplot(
        data=netrep,
        x="arm_display",
        y="significant_statistics_FDR05",
        color="#111827",
        alpha=0.55,
        size=3,
        ax=axes[1],
    )
    axes[1].set_ylim(-0.3, 7.3)
    axes[1].set_yticks(range(8))
    axes[1].set_title("NetRep seven-statistic support", loc="left", fontweight="bold")
    axes[1].set_xlabel("")
    axes[1].set_ylabel("Statistics significant at FDR <= 0.05")
    for axis in axes:
        axis.tick_params(axis="x", rotation=20)
    sns.despine(fig=fig)
    fig.suptitle(
        "Count-matched external module preservation",
        x=0.02,
        ha="left",
        fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    path = plot_dir / "official_wgcna_netrep_preservation.png"
    fig.savefig(path, dpi=240, bbox_inches="tight")
    plt.close(fig)
    return [path]


def write_report(
    output: Path,
    summary: pd.DataFrame,
    paired: pd.DataFrame,
    formal_summary: pd.DataFrame,
    formal_paired: pd.DataFrame,
    design: dict[str, object],
) -> Path:
    preservation_test = paired.loc[
        paired["metric"].eq("Zsummary_wgcna_like")
    ].iloc[0]
    score_test = paired.loc[
        paired["metric"].eq("program_score_preservation")
    ].iloc[0]
    preservation_delta = float(
        preservation_test["median_confidence_minus_matched"]
    )
    score_delta = float(score_test["median_confidence_minus_matched"])
    preservation_fdr = float(preservation_test["FDR"])
    score_fdr = float(score_test["FDR"])
    official_test = formal_paired.loc[
        formal_paired["metric"].eq("Zsummary.pres")
    ].iloc[0]
    netrep_test = formal_paired.loc[
        formal_paired["metric"].eq("significant_statistics_FDR05")
    ].iloc[0]
    official_delta = float(official_test["median_confidence_minus_matched"])
    official_fdr = float(official_test["FDR"])
    netrep_delta = float(netrep_test["median_confidence_minus_matched"])
    netrep_fdr = float(netrep_test["FDR"])
    verdict = (
        "Continuous CoSplice confidence adds externally replicated preservation value beyond count-matched kME selection."
        if official_delta > 0
        and netrep_delta > 0
        and official_fdr <= 0.05
        and netrep_fdr <= 0.05
        else "Official WGCNA and NetRep tests do not establish consistent superiority over count-matched kME selection; the confidence score remains an auditable prioritization layer."
    )
    report = output / "V122_MEMBERSHIP_EXTERNAL_VALIDATION.md"
    report.write_text(
        "\n".join(
            [
                "# CoSplice continuous membership-confidence external validation",
                "",
                "## Conclusion",
                "",
                verdict,
                "",
                f"- Discovery cohort: `{design['discovery']}`; independent evaluation cohorts: `{', '.join(design['holdouts'])}`.",
                f"- Official WGCNA confidence-minus-kME median Zsummary: `{official_delta:.3f}` (FDR `{official_fdr:.3g}`).",
                f"- NetRep confidence-minus-kME median number of supported statistics: `{netrep_delta:.3f}` (FDR `{netrep_fdr:.3g}`).",
                "- The official WGCNA/NetRep endpoints determine the verdict; the custom endpoints below are secondary diagnostics.",
                f"- Paired confidence minus matched top-kME median Zsummary: `{preservation_delta:.3f}` (FDR `{preservation_fdr:.3g}`).",
                f"- Paired confidence minus matched top-kME median program-score preservation: `{score_delta:.3f}` (FDR `{score_fdr:.3g}`).",
                "- Public SpliceSeq supplies PSI but not replicate IJC/SJC. This run validates sample-bootstrap membership confidence, not read-count uncertainty.",
                "",
                "## Four-arm control",
                "",
                "1. `wgcna_all_members`: the complete frozen WGCNA partition.",
                "2. `wgcna_top_kme_matched`: per-module top-|kME| events, matched exactly to the CoSplice core count.",
                "3. `cosplice_split_half_core`: events confirmed as core in both split halves after bootstrap kME/co-assignment assessment.",
                "4. `cosplice_confidence_matched`: the same per-module count, ranked by continuous co-assignment, kME, split-half and optional read-reliability evidence.",
                "",
                summary.to_markdown(index=False),
                "",
                "## Paired confidence-versus-kME tests",
                "",
                paired.to_markdown(index=False),
                "",
                "## Official WGCNA and NetRep preservation",
                "",
                formal_summary.to_markdown(index=False),
                "",
                formal_paired.to_markdown(index=False),
                "",
                "## Interpretation boundary",
                "",
                "The matched-count arm prevents a smaller event set from receiving automatic credit. Functional replication uses independently refitted holdout WGCNA modules and compares significant GO/KEGG themes after event-overlap matching. Phenotype association projects discovery-derived event weights into full tumor/normal matrices. Cross-cancer preservation is evidence of transportability, not disease-specific biological validation.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Validate continuous CoSplice membership confidence across independent public cohorts.")
    parser.add_argument("--panel-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--discovery", default="KIRC")
    parser.add_argument("--holdouts", nargs="+", default=["KIRP", "KICH", "LUAD", "LUSC"])
    parser.add_argument("--event-type", default="ES")
    parser.add_argument("--gmt", type=Path, action="append", required=True)
    parser.add_argument("--rscript", type=Path)
    parser.add_argument("--events", type=int, default=500)
    parser.add_argument("--bootstraps", type=int, default=20)
    parser.add_argument("--permutations", type=int, default=200)
    parser.add_argument(
        "--formal-permutations",
        type=int,
        default=1_000,
        help="Permutations for official WGCNA modulePreservation and NetRep.",
    )
    parser.add_argument(
        "--formal-workers",
        type=int,
        default=4,
        help="Parallel formal-preservation arms within one external cohort.",
    )
    parser.add_argument(
        "--formal-threads-per-arm",
        type=int,
        default=1,
        help="NetRep worker threads used by each formal-preservation arm.",
    )
    parser.add_argument(
        "--formal-cohort-workers",
        type=int,
        default=1,
        help="Independent external cohorts evaluated in parallel.",
    )
    parser.add_argument("--min-module-size", type=int, default=15)
    parser.add_argument("--functional-fdr", type=float, default=0.10)
    parser.add_argument("--random-seed", type=int, default=20260901)
    args = parser.parse_args(argv)

    output = args.output.resolve()
    prelock = output / "01_frozen_discovery"
    evaluation = output / "02_external_evaluation"
    prelock.mkdir(parents=True, exist_ok=True)
    evaluation.mkdir(parents=True, exist_ok=True)
    event_type = args.event_type.upper()
    discovery_paths = cohort_paths(args.panel_root.resolve(), args.discovery, event_type)
    design = {
        "registered_at_utc": utc_now(),
        "software_claim": "Katana-CoSplice continuous membership-confidence validation",
        "discovery": args.discovery,
        "holdouts": args.holdouts,
        "event_type": event_type,
        "arms": list(ARMS),
        "primary_comparison": f"{ARM_CONFIDENCE} versus {ARM_MATCHED}",
        "historical_use_caveat": "These public cohorts existed in prior development. This is a frozen rerun, not a globally untouched preregistered holdout experiment.",
        "read_uncertainty_boundary": "TCGA SpliceSeq PSI matrices do not contain replicate IJC/SJC; read-count uncertainty is not evaluated here.",
    }
    write_json(output / "analysis_plan.json", design)
    merged_gmt = combine_gmt_files(args.gmt, prelock / "frozen_GO_KEGG.gmt")
    discovery_raw = read_psi_matrix(discovery_paths["psi_tumor"])
    discovery = select_discovery_events(discovery_raw, args.events)
    samples = sample_frame(discovery_paths["samples"], discovery.columns, args.discovery)
    discovery.to_csv(prelock / "frozen_discovery_PSI.tsv", sep="\t")
    event_metadata = pd.read_csv(discovery_paths["events"], sep="\t")
    event_metadata.loc[event_metadata["event_key"].astype(str).isin(discovery.index)].to_csv(
        prelock / "frozen_event_metadata.tsv", sep="\t", index=False
    )

    if args.rscript:
        if not args.rscript.exists():
            raise FileNotFoundError(args.rscript)
        rscript = str(args.rscript.resolve())
        r_audit = f"Explicit Rscript: {rscript}"
    else:
        rscript, r_audit = _rscript_with_wgcna()
        if rscript is None:
            raise RuntimeError(r_audit)
    network = NetworkConfig(
        engine="tom_hierarchical",
        module_detection="wgcna_blockwise",
        method="bicor",
        signed=True,
        soft_power=6,
        min_module_size=args.min_module_size,
        max_events=None,
        tom_dense_max_events=max(args.events, 500),
        consensus_bootstrap_iterations=0,
        community_consensus_runs=1,
        bootstrap_iterations=0,
        wgcna_rscript=rscript,
        random_seed=args.random_seed,
    )
    benchmark = BenchmarkConfig(
        methods=["wgcna_exact"],
        max_events=len(discovery),
        bootstrap_iterations=args.bootstraps,
        random_seed=args.random_seed,
    )
    consensus = ConsensusConfig(
        bootstrap_iterations=args.bootstraps,
        clustering_method="backbone_confidence_only",
        min_module_size=args.min_module_size,
        dense_max_events=max(args.events, 500),
        random_seed=args.random_seed,
    )
    preservation_cfg = PreservationConfig(
        permutations=args.permutations,
        min_shared_events=max(5, min(10, args.min_module_size)),
        dense_max_events=max(args.events, 500),
        random_seed=args.random_seed,
    )
    print("Fitting frozen exact WGCNA discovery partition", flush=True)
    started = time.perf_counter()
    wgcna_labels, _ = _fit_method("wgcna_exact", discovery, network, benchmark, rscript)
    observed = assignment_frame(discovery.index, wgcna_labels)
    print("Estimating bootstrap and split-half membership confidence", flush=True)
    confidence = build_bootstrap_coassignment_consensus(
        discovery, observed, samples, network, consensus, progress=print
    )
    confidence.assignments.to_csv(prelock / "membership_confidence_assignments.tsv", sep="\t", index=False)
    arms, count_audit = build_membership_arms(confidence.assignments)
    count_audit.to_csv(prelock / "module_member_count_matching.tsv", sep="\t", index=False)
    for arm, frame in arms.items():
        frame.to_csv(prelock / f"{arm}_assignments.tsv", sep="\t", index=False)
    freeze_manifest = {
        "freeze_timestamp_utc": utc_now(),
        "holdout_files_opened_by_this_runner": False,
        "n_frozen_events": len(discovery),
        "runtime_to_freeze_seconds": time.perf_counter() - started,
        "settings": {
            "network": asdict(network),
            "benchmark": asdict(benchmark),
            "consensus": asdict(consensus),
            "preservation": asdict(preservation_cfg),
            "functional_fdr": args.functional_fdr,
        },
        "inputs": [
            file_record(discovery_paths[key], role="discovery")
            for key in ("psi_tumor", "events", "samples")
        ] + [file_record(path, role="functional_annotation") for path in args.gmt],
        "outputs": [file_record(prelock / f"{arm}_assignments.tsv", role="frozen_assignment") for arm in ARMS],
        "exact_wgcna_environment": r_audit,
    }
    write_json(output / "freeze_manifest.json", freeze_manifest)
    print("Discovery frozen; opening external cohorts", flush=True)

    gene_sets = load_gmt(merged_gmt)
    preservation_parts: list[pd.DataFrame] = []
    enrichment_parts: list[pd.DataFrame] = []
    functional_parts: list[pd.DataFrame] = []
    phenotype_parts: list[pd.DataFrame] = []
    formal_wgcna_parts: list[pd.DataFrame] = []
    formal_netrep_parts: list[pd.DataFrame] = []
    formal_audit_parts: list[pd.DataFrame] = []
    holdout_records: list[dict[str, object]] = []

    def evaluate_cohort(cohort: str) -> dict[str, object]:
        paths = cohort_paths(args.panel_root.resolve(), cohort, event_type)
        print(f"Evaluating {cohort}", flush=True)
        target_tumor_raw = read_psi_matrix(paths["psi_tumor"])
        target_tumor = target_tumor_raw.reindex(discovery.index)
        target_tumor = target_tumor.loc[
            target_tumor.notna().mean(axis=1).ge(0.75)
            & target_tumor.std(axis=1, skipna=True).gt(1e-8)
        ]
        target_tumor = target_tumor.apply(lambda row: row.fillna(row.mean()), axis=1)
        if len(target_tumor) < args.min_module_size:
            raise ValueError(f"{cohort} has only {len(target_tumor)} usable frozen events")
        holdout_labels, _ = _fit_method("wgcna_exact", target_tumor, network, benchmark, rscript)
        holdout_assignment = assignment_frame(target_tumor.index, holdout_labels)
        holdout_assignment.to_csv(evaluation / f"{cohort}_independent_wgcna_assignments.tsv", sep="\t", index=False)
        cohort_preservation: list[pd.DataFrame] = []
        for arm, frame in arms.items():
            external = ExternalCohortConfig(cohort, paths["psi_tumor"], paths["samples"])
            result = evaluate_external_module_preservation(
                discovery, frame, [external], preservation_cfg, progress=print
            )
            detail = result.module_preservation.copy()
            detail.insert(0, "arm", arm)
            cohort_preservation.append(detail)

        def formal_for_arm(item: tuple[str, pd.DataFrame]):
            arm, frame = item
            checkpoint = evaluation / "formal_checkpoints" / cohort / arm
            checkpoint.mkdir(parents=True, exist_ok=True)
            result = run_formal_preservation(
                discovery,
                target_tumor_raw,
                frame,
                cohort=cohort,
                arm=arm,
                rscript=rscript,
                permutations=args.formal_permutations,
                soft_power=network.soft_power,
                min_module_size=preservation_cfg.min_shared_events,
                random_seed=args.random_seed,
                n_threads=args.formal_threads_per_arm,
            )
            records = []
            for name in ("wgcna", "netrep", "audit"):
                path = checkpoint / f"{name}.tsv"
                temporary = path.with_suffix(".tsv.partial")
                getattr(result, name).to_csv(temporary, sep="\t", index=False)
                temporary.replace(path)
                records.append(file_record(path, role=f"formal_{name}"))
            write_json(checkpoint / "completed.json", {
                "completed_at": utc_now(), "cohort": cohort, "arm": arm,
                "permutations": args.formal_permutations, "seed": args.random_seed,
                "discovery_freeze": file_record(output / "freeze_manifest.json", role="freeze"),
                "target": file_record(paths["psi_tumor"], role="external_input"),
                "outputs": records,
            })
            print(f"Saved formal preservation: {cohort} / {arm}", flush=True)
            return result

        workers = max(1, min(args.formal_workers, len(arms)))
        print(
            f"Running formal WGCNA/NetRep preservation for {cohort} "
            f"with {workers} parallel arm(s)",
            flush=True,
        )
        with ThreadPoolExecutor(max_workers=workers) as executor:
            formal_results = executor.map(formal_for_arm, arms.items())
            formal_results = list(formal_results)
        enrichment, functional = functional_replication_by_overlap(
            cohort, arms, holdout_assignment, gene_sets, args.functional_fdr
        )
        target_full = read_psi_matrix(paths["psi_full"]).reindex(discovery.index)
        phenotype = phenotype_association(
            cohort,
            discovery,
            target_full,
            paths["samples"],
            arms,
            preservation_cfg.min_shared_events,
        )
        records = [
            file_record(paths[key], role=f"external_{cohort}")
            for key in ("psi_tumor", "psi_full", "events", "samples")
        ]
        return {
            "preservation": cohort_preservation,
            "enrichment": enrichment,
            "functional": functional,
            "phenotype": phenotype,
            "formal": formal_results,
            "records": records,
        }

    cohort_workers = max(
        1,
        min(args.formal_cohort_workers, len(args.holdouts)),
    )
    print(
        f"Evaluating {len(args.holdouts)} external cohort(s) with "
        f"{cohort_workers} cohort worker(s)",
        flush=True,
    )
    with ThreadPoolExecutor(max_workers=cohort_workers) as cohort_executor:
        cohort_results = cohort_executor.map(evaluate_cohort, args.holdouts)
        for cohort_result in cohort_results:
            preservation_parts.extend(cohort_result["preservation"])
            enrichment_parts.append(cohort_result["enrichment"])
            functional_parts.append(cohort_result["functional"])
            phenotype_parts.append(cohort_result["phenotype"])
            for formal in cohort_result["formal"]:
                formal_wgcna_parts.append(formal.wgcna)
                formal_netrep_parts.append(formal.netrep)
                formal_audit_parts.append(formal.audit)
            holdout_records.extend(cohort_result["records"])

    preservation = pd.concat(preservation_parts, ignore_index=True)
    enrichment = pd.concat(enrichment_parts, ignore_index=True)
    functional = pd.concat(functional_parts, ignore_index=True)
    phenotype = pd.concat(phenotype_parts, ignore_index=True)
    formal_wgcna = pd.concat(formal_wgcna_parts, ignore_index=True)
    formal_netrep = pd.concat(formal_netrep_parts, ignore_index=True)
    formal_audit = pd.concat(formal_audit_parts, ignore_index=True)
    preservation.to_csv(evaluation / "external_module_preservation.tsv", sep="\t", index=False)
    enrichment.to_csv(evaluation / "independent_GO_KEGG_enrichment.tsv", sep="\t", index=False)
    functional.to_csv(evaluation / "functional_replication.tsv", sep="\t", index=False)
    phenotype.to_csv(evaluation / "projected_phenotype_association.tsv", sep="\t", index=False)
    formal_wgcna.to_csv(
        evaluation / "official_wgcna_module_preservation.tsv",
        sep="\t",
        index=False,
    )
    formal_netrep.to_csv(
        evaluation / "netrep_seven_statistic_preservation.tsv",
        sep="\t",
        index=False,
    )
    formal_audit.to_csv(
        evaluation / "formal_preservation_audit.tsv",
        sep="\t",
        index=False,
    )

    summary = preservation.groupby("arm", observed=True).agg(
        tested_module_cohort_pairs=("consensus_module", "size"),
        median_Zsummary=("Zsummary_wgcna_like", "median"),
        median_program_score_preservation=("program_score_preservation", "median"),
        preserved_pairs=("preservation_class", lambda x: int((x == "preserved").sum())),
        at_least_partially_preserved_pairs=("preservation_class", lambda x: int(x.isin(["preserved", "partially_preserved"]).sum())),
    ).reset_index()
    functional_summary = functional.groupby("arm", observed=True).agg(
        modules_with_replicated_terms=("shared_significant_terms", lambda x: int((x > 0).sum())),
        median_functional_term_jaccard=("term_jaccard", "median"),
    ).reset_index()
    phenotype_summary = phenotype.assign(significant=phenotype["FDR"].le(0.05)).groupby("arm", observed=True).agg(
        significant_phenotype_associations=("significant", "sum"),
        median_abs_phenotype_effect=("cohens_d", lambda x: float(np.median(np.abs(x)))),
    ).reset_index()
    summary = summary.merge(functional_summary, on="arm", how="left").merge(phenotype_summary, on="arm", how="left")
    summary.to_csv(output / "four_arm_external_validation_summary.tsv", sep="\t", index=False)
    paired = pd.DataFrame([
        paired_confidence_test(preservation, "Zsummary_wgcna_like"),
        paired_confidence_test(preservation, "program_score_preservation"),
    ])
    paired["FDR"] = bh_adjust(paired["wilcoxon_pvalue"])
    paired.to_csv(output / "paired_confidence_vs_top_kme.tsv", sep="\t", index=False)
    formal_summary = formal_wgcna.groupby("arm", observed=True).agg(
        official_module_cohort_pairs=("module", "size"),
        median_official_Zsummary=("Zsummary.pres", "median"),
        strong_official_preservation=("Zsummary.pres", lambda x: int((x >= 10).sum())),
        moderate_or_strong_official_preservation=("Zsummary.pres", lambda x: int((x >= 2).sum())),
    ).reset_index()
    netrep_summary = formal_netrep.groupby("arm", observed=True).agg(
        median_significant_NetRep_statistics=("significant_statistics_FDR05", "median"),
        modules_all_seven_NetRep_FDR05=("all_seven_statistics_FDR05", "sum"),
    ).reset_index()
    formal_summary = formal_summary.merge(netrep_summary, on="arm", how="left")
    formal_summary.to_csv(
        output / "formal_preservation_summary.tsv", sep="\t", index=False
    )
    formal_wgcna_for_pair = formal_wgcna.rename(
        columns={"module": "consensus_module"}
    )
    formal_netrep_for_pair = formal_netrep.rename(
        columns={"module": "consensus_module"}
    )
    formal_paired = pd.DataFrame(
        [
            paired_confidence_test(
                formal_wgcna_for_pair, "Zsummary.pres"
            ),
            paired_confidence_test(
                formal_netrep_for_pair,
                "significant_statistics_FDR05",
            ),
        ]
    )
    formal_paired["FDR"] = bh_adjust(formal_paired["wilcoxon_pvalue"])
    formal_paired.to_csv(
        output / "formal_paired_confidence_vs_top_kme.tsv",
        sep="\t",
        index=False,
    )
    generated = render_plots(output, preservation, functional, phenotype)
    generated.extend(
        render_formal_preservation_plots(
            output,
            formal_wgcna,
            formal_netrep,
        )
    )
    report = write_report(
        output,
        summary,
        paired,
        formal_summary,
        formal_paired,
        design,
    )
    write_json(
        output / "parameters.json",
        {
            **design,
            "completed_at_utc": utc_now(),
            "holdout_inputs": holdout_records,
            "generated_plots": [str(path.resolve()) for path in generated],
            "formal_preservation_permutations": args.formal_permutations,
            "formal_preservation_workers": args.formal_workers,
            "formal_preservation_threads_per_arm": args.formal_threads_per_arm,
            "formal_preservation_cohort_workers": args.formal_cohort_workers,
            "report": str(report.resolve()),
        },
    )
    print(summary.to_string(index=False), flush=True)
    print(f"Output: {output}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
