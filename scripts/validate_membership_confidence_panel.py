from __future__ import annotations

import argparse
import json
import math
import time
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


def paired_core_test(table: pd.DataFrame, metric: str) -> dict[str, object]:
    subset = table.loc[table["arm"].isin([ARM_MATCHED, ARM_CORE])]
    pivot = subset.pivot_table(index=["cohort", "consensus_module"], columns="arm", values=metric)
    pivot = pivot.dropna(subset=[ARM_MATCHED, ARM_CORE])
    if pivot.empty:
        return {"metric": metric, "n_pairs": 0, "median_core_minus_matched": np.nan, "wilcoxon_pvalue": np.nan}
    delta = pivot[ARM_CORE] - pivot[ARM_MATCHED]
    try:
        pvalue = float(wilcoxon(delta, zero_method="zsplit").pvalue)
    except ValueError:
        pvalue = 1.0
    return {
        "metric": metric,
        "n_pairs": len(delta),
        "median_core_minus_matched": float(delta.median()),
        "mean_core_minus_matched": float(delta.mean()),
        "core_wins": int((delta > 0).sum()),
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
    path = plot_dir / "external_preservation_three_arm.png"
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


def write_report(
    output: Path,
    summary: pd.DataFrame,
    paired: pd.DataFrame,
    design: dict[str, object],
) -> Path:
    preservation_test = paired.loc[
        paired["metric"].eq("Zsummary_wgcna_like")
    ].iloc[0]
    score_test = paired.loc[
        paired["metric"].eq("program_score_preservation")
    ].iloc[0]
    preservation_delta = float(
        preservation_test["median_core_minus_matched"]
    )
    score_delta = float(score_test["median_core_minus_matched"])
    preservation_fdr = float(preservation_test["FDR"])
    score_fdr = float(score_test["FDR"])
    verdict = (
        "CoSplice core selection adds external value beyond count-matched kME selection."
        if preservation_delta > 0
        and score_delta > 0
        and preservation_fdr <= 0.05
        and score_fdr <= 0.05
        else "CoSplice core selection does not consistently exceed the count-matched kME baseline; publication claims must remain limited."
    )
    report = output / "V121_MEMBERSHIP_EXTERNAL_VALIDATION.md"
    report.write_text(
        "\n".join(
            [
                "# CoSplice v1.2.1 membership-confidence external validation",
                "",
                "## Conclusion",
                "",
                verdict,
                "",
                f"- Discovery cohort: `{design['discovery']}`; independent evaluation cohorts: `{', '.join(design['holdouts'])}`.",
                f"- Paired core minus matched top-kME median Zsummary: `{preservation_delta:.3f}` (FDR `{preservation_fdr:.3g}`).",
                f"- Paired core minus matched top-kME median program-score preservation: `{score_delta:.3f}` (FDR `{score_fdr:.3g}`).",
                "- Public SpliceSeq supplies PSI but not replicate IJC/SJC. This run validates sample-bootstrap membership confidence, not read-count uncertainty.",
                "",
                "## Three-arm control",
                "",
                "1. `wgcna_all_members`: the complete frozen WGCNA partition.",
                "2. `wgcna_top_kme_matched`: per-module top-|kME| events, matched exactly to the CoSplice core count.",
                "3. `cosplice_split_half_core`: events confirmed as core in both split halves after bootstrap kME/co-assignment assessment.",
                "",
                summary.to_markdown(index=False),
                "",
                "## Paired core-versus-matched tests",
                "",
                paired.to_markdown(index=False),
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
    parser = argparse.ArgumentParser(description="Validate v1.2.1 split-half membership confidence across independent public cohorts.")
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
        "software_claim": "Katana-CoSplice v1.2.1 membership-confidence validation",
        "discovery": args.discovery,
        "holdouts": args.holdouts,
        "event_type": event_type,
        "arms": list(ARMS),
        "primary_comparison": f"{ARM_CORE} versus {ARM_MATCHED}",
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
    holdout_records: list[dict[str, object]] = []
    for cohort in args.holdouts:
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
        for arm, frame in arms.items():
            external = ExternalCohortConfig(cohort, paths["psi_tumor"], paths["samples"])
            result = evaluate_external_module_preservation(
                discovery, frame, [external], preservation_cfg, progress=print
            )
            detail = result.module_preservation.copy()
            detail.insert(0, "arm", arm)
            preservation_parts.append(detail)
        enrichment, functional = functional_replication_by_overlap(
            cohort, arms, holdout_assignment, gene_sets, args.functional_fdr
        )
        enrichment_parts.append(enrichment)
        functional_parts.append(functional)
        target_full = read_psi_matrix(paths["psi_full"]).reindex(discovery.index)
        phenotype_parts.append(
            phenotype_association(
                cohort,
                discovery,
                target_full,
                paths["samples"],
                arms,
                preservation_cfg.min_shared_events,
            )
        )
        holdout_records.extend(file_record(paths[key], role=f"external_{cohort}") for key in ("psi_tumor", "psi_full", "events", "samples"))

    preservation = pd.concat(preservation_parts, ignore_index=True)
    enrichment = pd.concat(enrichment_parts, ignore_index=True)
    functional = pd.concat(functional_parts, ignore_index=True)
    phenotype = pd.concat(phenotype_parts, ignore_index=True)
    preservation.to_csv(evaluation / "external_module_preservation.tsv", sep="\t", index=False)
    enrichment.to_csv(evaluation / "independent_GO_KEGG_enrichment.tsv", sep="\t", index=False)
    functional.to_csv(evaluation / "functional_replication.tsv", sep="\t", index=False)
    phenotype.to_csv(evaluation / "projected_phenotype_association.tsv", sep="\t", index=False)

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
    summary.to_csv(output / "three_arm_external_validation_summary.tsv", sep="\t", index=False)
    paired = pd.DataFrame([
        paired_core_test(preservation, "Zsummary_wgcna_like"),
        paired_core_test(preservation, "program_score_preservation"),
    ])
    paired["FDR"] = bh_adjust(paired["wilcoxon_pvalue"])
    paired.to_csv(output / "paired_core_vs_top_kme.tsv", sep="\t", index=False)
    generated = render_plots(output, preservation, functional, phenotype)
    report = write_report(output, summary, paired, design)
    write_json(
        output / "parameters.json",
        {
            **design,
            "completed_at_utc": utc_now(),
            "holdout_inputs": holdout_records,
            "generated_plots": [str(path.resolve()) for path in generated],
            "report": str(report.resolve()),
        },
    )
    print(summary.to_string(index=False), flush=True)
    print(f"Output: {output}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
