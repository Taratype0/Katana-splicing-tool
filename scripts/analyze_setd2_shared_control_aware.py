from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

from katana_cosplice.matrix import EVENT_TYPES, _event_frame
from katana_cosplice.shared_control import (
    blocked_shared_control_permutation,
    fit_shared_control_beta_binomial,
)
from katana_cosplice.utils import bh_adjust


COLORS = {
    "observed": "#d1495b",
    "null": "#4c78a8",
    "similar": "#218380",
    "different": "#d1495b",
    "uncertain": "#9aa5b1",
    "ink": "#172033",
}


def _read_tsv(path: Path, **kwargs: object) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t", low_memory=False, **kwargs)


def _save(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, dpi=240, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _load_rmats_comparison(
    rmats_root: Path,
    comparison: str,
    mode: str,
) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for event_type in EVENT_TYPES:
        path = rmats_root / comparison / "rmats_post" / f"{event_type}.MATS.{mode}.txt"
        raw = _read_tsv(path)
        identity = _event_frame(raw, event_type)
        values = identity.merge(
            raw.assign(_row=np.arange(len(raw)))[
                [
                    "_row",
                    "FDR",
                    "IncLevelDifference",
                    "IncLevel1",
                    "IncLevel2",
                ]
            ],
            left_index=True,
            right_on="_row",
            how="left",
        ).drop(columns="_row")
        values["FDR"] = pd.to_numeric(values["FDR"], errors="coerce")
        values["dPSI"] = pd.to_numeric(
            values["IncLevelDifference"],
            errors="coerce",
        )
        values["comparison_id"] = comparison
        frames.append(values)
    result = pd.concat(frames, ignore_index=True)
    result = (
        result.sort_values(["event_key", "FDR"], na_position="last")
        .drop_duplicates("event_key")
        .reset_index(drop=True)
    )
    result["significant"] = result["FDR"].le(0.05) & result["dPSI"].abs().ge(0.10)
    return result


def _sample_ids(context: str, inhibitor: str) -> dict[str, list[str]]:
    suffix = "DP" if context == "DP" else "NO_DP"
    return {
        "WT": [f"WT_{suffix}_rep_{replicate}" for replicate in range(1, 4)],
        "MUT": [f"MUT_{suffix}_rep_{replicate}" for replicate in range(1, 4)],
        "inhibitor": [
            f"WT_{suffix}_{inhibitor}_rep_{replicate}"
            for replicate in range(1, 4)
        ],
    }


def _permutation_column_order(sample_ids: dict[str, list[str]]) -> list[str]:
    return [
        sample_ids[group][replicate]
        for replicate in range(3)
        for group in ("WT", "MUT", "inhibitor")
    ]


def _model_column_order(sample_ids: dict[str, list[str]]) -> list[str]:
    return [
        sample_id
        for group in ("WT", "MUT", "inhibitor")
        for sample_id in sample_ids[group]
    ]


def _load_sample_matrix(
    path: Path,
    columns: list[str],
) -> pd.DataFrame:
    frame = _read_tsv(path, usecols=["event_key", *columns], index_col="event_key")
    frame.index = frame.index.astype(str)
    return frame.apply(pd.to_numeric, errors="coerce")


def _build_shared_events(
    mut: pd.DataFrame,
    inhibitor: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    tested = mut[
        ["event_key", "gene_id", "gene_symbol", "event_type", "FDR", "dPSI", "significant"]
    ].merge(
        inhibitor[["event_key", "FDR", "dPSI", "significant"]],
        on="event_key",
        suffixes=("_MUT", "_inhibitor"),
        validate="one_to_one",
    )
    shared = tested.loc[
        tested["significant_MUT"] & tested["significant_inhibitor"]
    ].copy()
    shared["same_direction_original"] = (
        shared["dPSI_MUT"] * shared["dPSI_inhibitor"] > 0
    )
    shared["original_abs_delta_dPSI"] = (
        shared["dPSI_MUT"] - shared["dPSI_inhibitor"]
    ).abs()
    return tested, shared


def _permutation_analysis(
    context: str,
    inhibitor_name: str,
    tested: pd.DataFrame,
    shared: pd.DataFrame,
    psi: pd.DataFrame,
    sample_ids: dict[str, list[str]],
) -> tuple[dict[str, object], pd.DataFrame]:
    universe = pd.Index(tested["event_key"].astype(str))
    available = universe.intersection(psi.index)
    permutation_columns = _permutation_column_order(sample_ids)
    # Complete cases keep every event exchangeable under all 216 assignments.
    values = psi.loc[available, permutation_columns].dropna()
    n_mut = int(tested["significant_MUT"].sum())
    n_inhibitor = int(tested["significant_inhibitor"].sum())
    permutation = blocked_shared_control_permutation(
        values,
        n_mut=min(n_mut, len(values)),
        n_drug=min(n_inhibitor, len(values)),
    )
    null = permutation.null_distribution.copy()
    null["context"] = context
    null["inhibitor"] = inhibitor_name
    null["rank_matched_same_direction_fraction"] = np.where(
        null["rank_matched_overlap"].gt(0),
        null["rank_matched_same_direction"] / null["rank_matched_overlap"],
        np.nan,
    )
    original_overlap = len(shared)
    original_same = int(shared["same_direction_original"].sum())
    original_fraction = original_same / original_overlap if original_overlap else np.nan
    original_overlap_empirical_pvalue = float(
        (1 + null["rank_matched_overlap"].ge(original_overlap).sum())
        / (len(null) + 1)
    )
    original_direction_empirical_pvalue = float(
        (
            1
            + null["rank_matched_same_direction_fraction"]
            .ge(original_fraction)
            .sum()
        )
        / (null["rank_matched_same_direction_fraction"].notna().sum() + 1)
    )
    expected_hypergeometric = n_mut * n_inhibitor / len(tested) if len(tested) else np.nan
    row = {
        "context": context,
        "inhibitor": inhibitor_name,
        "n_exact_tested": len(tested),
        "n_events_with_sample_PSI": len(values),
        "n_MUT_significant": n_mut,
        "n_inhibitor_significant": n_inhibitor,
        "n_original_shared_significant": original_overlap,
        "n_original_shared_same_direction": original_same,
        "original_shared_same_direction_pct": 100.0 * original_fraction,
        "naive_hypergeometric_expected_overlap": expected_hypergeometric,
        "naive_hypergeometric_enrichment_ratio": (
            original_overlap / expected_hypergeometric
            if expected_hypergeometric
            else np.nan
        ),
        "naive_hypergeometric_pvalue": (
            float(
                stats.hypergeom.sf(
                    original_overlap - 1,
                    len(tested),
                    n_inhibitor,
                    n_mut,
                )
            )
            if original_overlap
            else np.nan
        ),
        "rank_matched_observed_overlap": permutation.observed_rank_overlap,
        "rank_matched_null_overlap_median": float(
            null["rank_matched_overlap"].median()
        ),
        "rank_matched_null_overlap_q025": float(
            null["rank_matched_overlap"].quantile(0.025)
        ),
        "rank_matched_null_overlap_q975": float(
            null["rank_matched_overlap"].quantile(0.975)
        ),
        "rank_matched_empirical_overlap_pvalue": permutation.empirical_overlap_pvalue,
        "original_overlap_vs_ranked_null_empirical_pvalue": (
            original_overlap_empirical_pvalue
        ),
        "original_direction_vs_ranked_null_empirical_pvalue": (
            original_direction_empirical_pvalue
        ),
        "permutation_scheme": (
            "exact enumeration of 3!^3=216 condition-label assignments within "
            "replicate/batch blocks; marginal selected-set sizes fixed to original rMATS counts"
        ),
        "interpretation_boundary": (
            "Sensitivity analysis for shared-control dependence; it does not "
            "reproduce the rMATS hierarchical likelihood under every permutation."
        ),
    }
    return row, null


def _fit_joint_events(
    context: str,
    inhibitor_name: str,
    shared: pd.DataFrame,
    psi: pd.DataFrame,
    ijc: pd.DataFrame,
    sjc: pd.DataFrame,
    sample_ids: dict[str, list[str]],
    equivalence_margin: float,
) -> pd.DataFrame:
    columns = _model_column_order(sample_ids)
    available = (
        pd.Index(shared["event_key"].astype(str))
        .intersection(psi.index)
        .intersection(ijc.index)
        .intersection(sjc.index)
    )
    metadata = shared.set_index("event_key").loc[available].reset_index()
    rows: list[dict[str, object]] = []
    for index, event_key in enumerate(available, start=1):
        if index == 1 or index % 50 == 0 or index == len(available):
            print(
                f"[{context} {inhibitor_name}] joint model "
                f"{index}/{len(available)}"
            )
        fitted = fit_shared_control_beta_binomial(
            ijc.loc[event_key, columns].to_numpy(dtype=float),
            sjc.loc[event_key, columns].to_numpy(dtype=float),
            psi.loc[event_key, columns].to_numpy(dtype=float),
            equivalence_margin=equivalence_margin,
        )
        rows.append({"event_key": event_key, **fitted})
    models = pd.DataFrame(rows)
    result = metadata.merge(models, on="event_key", how="left", validate="one_to_one")
    result.insert(0, "context", context)
    result.insert(1, "inhibitor", inhibitor_name)
    for source, target in [
        ("joint_MUT_vs_WT_pvalue", "joint_MUT_vs_WT_FDR"),
        ("joint_inhibitor_vs_WT_pvalue", "joint_inhibitor_vs_WT_FDR"),
        ("direct_MUT_vs_inhibitor_pvalue", "direct_MUT_vs_inhibitor_FDR"),
        ("equivalence_pvalue", "equivalence_FDR"),
    ]:
        result[target] = bh_adjust(pd.to_numeric(result.get(source), errors="coerce"))
    practical_similarity = (
        result["practical_similarity"]
        .astype("boolean")
        .fillna(False)
        .astype(bool)
    )
    result["formal_equivalence_FDR_0_05"] = (
        result["equivalence_FDR"].le(0.05)
        & result["equivalence_ci90_low"].gt(-equivalence_margin)
        & result["equivalence_ci90_high"].lt(equivalence_margin)
    )
    result["direct_difference_FDR_0_05"] = result[
        "direct_MUT_vs_inhibitor_FDR"
    ].le(0.05)
    result["strict_LOF_compatible_core"] = (
        result["same_direction_original"]
        & result["formal_equivalence_FDR_0_05"]
        & ~result["direct_difference_FDR_0_05"]
    )
    result["support_class"] = np.select(
        [
            result["strict_LOF_compatible_core"],
            result["same_direction_original"]
            & practical_similarity,
            result["direct_difference_FDR_0_05"],
        ],
        [
            "formal_equivalent_same_direction",
            "practically_similar_same_direction_not_formally_equivalent",
            "directly_different_MUT_vs_inhibitor",
        ],
        default="uncertain_or_direction_discordant",
    )
    return result


def _summary_joint_models(
    models: pd.DataFrame,
    equivalence_margin: float,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (context, inhibitor), subset in models.groupby(
        ["context", "inhibitor"],
        observed=True,
    ):
        rows.append(
            {
                "context": context,
                "inhibitor": inhibitor,
                "n_original_shared_significant": len(subset),
                "n_joint_model_converged": int(
                    subset["joint_model_status"].eq("converged").sum()
                ),
                "n_same_direction_original": int(
                    subset["same_direction_original"].sum()
                ),
                "n_practically_similar_same_direction": int(
                    (
                        subset["same_direction_original"]
                        & subset["practical_similarity"]
                        .astype("boolean")
                        .fillna(False)
                        .astype(bool)
                    ).sum()
                ),
                "n_formally_equivalent_FDR_0_05": int(
                    subset["formal_equivalence_FDR_0_05"].sum()
                ),
                "n_directly_different_FDR_0_05": int(
                    subset["direct_difference_FDR_0_05"].sum()
                ),
                "n_strict_LOF_compatible_core": int(
                    subset["strict_LOF_compatible_core"].sum()
                ),
                "median_abs_direct_MUT_vs_inhibitor_delta_PSI": float(
                    pd.to_numeric(
                        subset["direct_delta_psi"],
                        errors="coerce",
                    ).abs().median()
                ),
                "equivalence_margin": equivalence_margin,
            }
        )
    return pd.DataFrame(rows)


def _plot_permutation_null(
    summary: pd.DataFrame,
    null: pd.DataFrame,
    path: Path,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    for ax, ((context, inhibitor), subset) in zip(
        axes.flat,
        null.groupby(["context", "inhibitor"], sort=False),
        strict=False,
    ):
        row = summary.loc[
            summary["context"].eq(context) & summary["inhibitor"].eq(inhibitor)
        ].iloc[0]
        sns.histplot(
            subset["rank_matched_overlap"],
            bins=24,
            color=COLORS["null"],
            ax=ax,
        )
        ax.axvline(
            row["n_original_shared_significant"],
            color=COLORS["observed"],
            linewidth=2.5,
            label=f"Original rMATS overlap = {int(row['n_original_shared_significant'])}",
        )
        ax.axvline(
            row["rank_matched_null_overlap_median"],
            color=COLORS["ink"],
            linestyle="--",
            linewidth=1.5,
            label=f"Shared-control null median = {row['rank_matched_null_overlap_median']:.0f}",
        )
        ax.set_title(f"{context}: MUT vs {inhibitor}", fontweight="bold")
        ax.set_xlabel("Rank-matched exact-event overlap")
        ax.set_ylabel("Blocked permutations")
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle(
        "Shared-control-preserving overlap sensitivity analysis",
        fontsize=16,
        fontweight="bold",
    )
    _save(fig, path)


def _plot_joint_validation(models: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13, 10), constrained_layout=True)
    for ax, ((context, inhibitor), subset) in zip(
        axes.flat,
        models.groupby(["context", "inhibitor"], sort=False),
        strict=False,
    ):
        colors = subset["support_class"].map(
            {
                "formal_equivalent_same_direction": COLORS["similar"],
                "practically_similar_same_direction_not_formally_equivalent": "#6ab0a8",
                "directly_different_MUT_vs_inhibitor": COLORS["different"],
                "uncertain_or_direction_discordant": COLORS["uncertain"],
            }
        )
        ax.scatter(
            subset["dPSI_MUT"],
            subset["dPSI_inhibitor"],
            c=colors,
            s=18,
            alpha=0.75,
            edgecolors="none",
        )
        limits = np.nanmax(
            np.abs(
                subset[["dPSI_MUT", "dPSI_inhibitor"]].to_numpy(dtype=float)
            )
        )
        limits = max(float(limits), 0.1)
        ax.plot([-limits, limits], [-limits, limits], color=COLORS["ink"], linestyle="--")
        ax.axhline(0, color="#c6cbd2", linewidth=0.8)
        ax.axvline(0, color="#c6cbd2", linewidth=0.8)
        ax.set_xlim(-limits, limits)
        ax.set_ylim(-limits, limits)
        ax.set_title(f"{context}: MUT vs {inhibitor}", fontweight="bold")
        ax.set_xlabel("rMATS dPSI: MUT vs WT")
        ax.set_ylabel(f"rMATS dPSI: WT+{inhibitor} vs WT")
        strict = int(subset["strict_LOF_compatible_core"].sum())
        direct = int(subset["direct_difference_FDR_0_05"].sum())
        ax.text(
            0.03,
            0.97,
            f"strict equivalent: {strict}\ndirectly different: {direct}",
            transform=ax.transAxes,
            va="top",
            fontsize=8,
        )
    fig.suptitle(
        "Direct shared-control-aware validation of exact shared AS events",
        fontsize=16,
        fontweight="bold",
    )
    fig.legend(
        handles=[
            Patch(
                color=COLORS["similar"],
                label="Formal equivalent, same direction",
            ),
            Patch(
                color="#6ab0a8",
                label="Practical similarity only",
            ),
            Patch(
                color=COLORS["different"],
                label="Direct MUT–inhibitor difference, FDR ≤ 0.05",
            ),
            Patch(
                color=COLORS["uncertain"],
                label="Uncertain or direction-discordant",
            ),
        ],
        loc="lower center",
        bbox_to_anchor=(0.5, -0.015),
        ncol=2,
        frameon=False,
    )
    _save(fig, path)


def _write_report(
    path: Path,
    permutation_summary: pd.DataFrame,
    joint_summary: pd.DataFrame,
    equivalence_margin: float,
) -> None:
    lines = [
        "# Setd2 shared-control-aware LOF sensitivity analysis",
        "",
        "## Why this reanalysis was required",
        "",
        "The original MUT-vs-WT and inhibitor-vs-WT contrasts reuse the same WT replicates. "
        "A marginal hypergeometric overlap test treats the two significant-event sets as "
        "independent and can therefore overstate overlap evidence. The exact event overlap "
        "remains a real descriptive observation, but its naive P value is not the primary "
        "inferential result.",
        "",
        "## Corrected analysis",
        "",
        "1. A blocked condition-label permutation preserves the three replicate/batch blocks "
        "and the shared-control design. Selected-set sizes are fixed to the original rMATS "
        "margins, so the null asks how much overlap can arise from shared controls and "
        "event-level effect propensity.",
        "2. Every original exact shared significant event is re-fit in one beta-binomial "
        "WT/MUT/inhibitor model using replicate IJC/SJC counts.",
        "3. MUT and inhibitor are compared directly. A separate Welch TOST uses replicate PSI "
        f"with a practical equivalence margin of ±{equivalence_margin:.2f}. Failure to detect "
        "a difference is not treated as equivalence.",
        "",
        "## Overlap sensitivity summary",
        "",
        permutation_summary.to_markdown(index=False),
        "",
        "## Direct event-model summary",
        "",
        joint_summary.to_markdown(index=False),
        "",
        "## Interpretation boundary",
        "",
        "- `naive_hypergeometric_pvalue` is retained for traceability, not as the main evidence.",
        "- The blocked permutation is a sample-level PSI sensitivity analysis and does not "
        "re-run the full rMATS hierarchical likelihood for every label assignment.",
        "- `practically_similar` means the observed MUT–inhibitor PSI difference is within the "
        "chosen margin. `formal_equivalence` additionally requires the 90% CI and multiplicity-"
        "adjusted TOST to support equivalence.",
        "- RNA-seq convergence is compatible with a partial LOF-like response but does not by "
        "itself establish biochemical SETD2 loss of function.",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def _write_html_report(
    path: Path,
    permutation_summary: pd.DataFrame,
    joint_summary: pd.DataFrame,
    equivalence_margin: float,
) -> None:
    overlap_table = permutation_summary[
        [
            "context",
            "inhibitor",
            "n_original_shared_significant",
            "naive_hypergeometric_expected_overlap",
            "naive_hypergeometric_pvalue",
            "rank_matched_null_overlap_median",
            "original_overlap_vs_ranked_null_empirical_pvalue",
            "original_shared_same_direction_pct",
            "original_direction_vs_ranked_null_empirical_pvalue",
        ]
    ].copy()
    model_table = joint_summary.copy()
    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Setd2 shared-control-aware LOF analysis</title>
<style>
body {{ margin: 0; font: 16px/1.55 Arial, sans-serif; color: #172033; background: #f3f5f7; }}
main {{ max-width: 1200px; margin: 0 auto; padding: 42px 28px 80px; }}
h1 {{ font-size: 34px; margin: 0 0 10px; }}
h2 {{ margin-top: 38px; }}
.lede {{ font-size: 19px; color: #4f5b6b; max-width: 920px; }}
.warning {{ background: #fff2cc; border-left: 6px solid #e5a823; padding: 18px 20px; margin: 24px 0; }}
.result {{ background: #e8f4f1; border-left: 6px solid #218380; padding: 18px 20px; margin: 24px 0; }}
.card {{ background: white; padding: 24px; margin: 20px 0; border-radius: 12px; box-shadow: 0 6px 22px #17203312; }}
table {{ border-collapse: collapse; width: 100%; font-size: 13px; }}
th, td {{ border-bottom: 1px solid #dfe3e8; padding: 8px; text-align: right; }}
th:first-child, td:first-child, th:nth-child(2), td:nth-child(2) {{ text-align: left; }}
img {{ width: 100%; height: auto; display: block; }}
code {{ background: #edf0f3; padding: 2px 5px; border-radius: 4px; }}
</style>
</head>
<body><main>
<h1>Setd2 shared-control-aware LOF sensitivity analysis</h1>
<p class="lede">This report supersedes the naive claim that exact overlap is “far above chance.” MUT and inhibitor contrasts reuse the same WT replicates, so overlap, direction agreement and effect correlation must be evaluated with that dependence preserved.</p>
<div class="warning"><b>Corrected boundary.</b> The original hypergeometric P values are retained only for traceability. They assume independent event sets and are anti-conservative for this shared-control design.</div>
<h2>1. Shared-control-preserving overlap sensitivity</h2>
<div class="card">{overlap_table.to_html(index=False, float_format=lambda value: f"{value:.4g}")}</div>
<div class="card"><img src="plots/01_shared_control_permutation_null.png" alt="Shared-control permutation null"></div>
<p>The blocked null enumerates all 216 condition-label assignments within the three replicate/batch blocks and fixes the marginal selected-set sizes. It is a sensitivity analysis, not a full rMATS re-fit under every permutation.</p>
<h2>2. Direct WT/MUT/inhibitor event model</h2>
<div class="card">{model_table.to_html(index=False, float_format=lambda value: f"{value:.4g}")}</div>
<div class="card"><img src="plots/02_direct_MUT_inhibitor_event_validation.png" alt="Direct event validation"></div>
<div class="result"><b>Strict result.</b> A strict LOF-compatible event requires original same-direction significance, direct MUT–inhibitor equivalence after multiplicity correction, and no significant direct difference. The equivalence margin is ±{equivalence_margin:.2f} PSI.</div>
<h2>3. Scientific conclusion</h2>
<p>The broad “non-random overlap” claim is not supported after accounting for the shared-control structure. A limited direct-equivalence core remains and can be prioritized for validation. This supports partial molecular convergence only; it does not establish complete inhibitor phenocopy or biochemical SETD2 loss of function.</p>
</main></body></html>"""
    path.write_text(html, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Reanalyze Setd2 MUT-inhibitor AS overlap with shared-control-aware "
            "permutations and direct replicate-count models."
        )
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(r"E:\Project\RNASEQ\Result\Setd2_results"),
    )
    parser.add_argument("--rmats-mode", choices=["JC", "JCEC"], default="JC")
    parser.add_argument("--equivalence-margin", type=float, default=0.10)
    args = parser.parse_args()

    project = args.project_root
    network_root = project / "katana_output" / "09_cosplice_network"
    matrix_root = network_root / "02_matrix"
    output = network_root / "06k_shared_control_aware_lof"
    plots = output / "plots"
    output.mkdir(parents=True, exist_ok=True)
    plots.mkdir(parents=True, exist_ok=True)

    comparisons = {
        "DP": {
            "MUT": "MUT_DP-WT_DP",
            "CLK3": "WT_DP_CLK3-WT_DP",
            "Pla": "WT_DP_Pla-WT_DP",
        },
        "noDP": {
            "MUT": "MUT_NO_DP-WT_NO_DP",
            "CLK3": "WT_NO_DP_CLK3-WT_NO_DP",
            "Pla": "WT_NO_DP_Pla-WT_NO_DP",
        },
    }
    all_sample_ids = sorted(
        {
            sample_id
            for context in comparisons
            for inhibitor in ("CLK3", "Pla")
            for sample_id in sum(_sample_ids(context, inhibitor).values(), [])
        }
    )
    print("Loading sample-level PSI/IJC/SJC matrices...")
    psi = _load_sample_matrix(
        matrix_root / "sample_level_psi_matrix.tsv",
        all_sample_ids,
    )
    ijc = _load_sample_matrix(
        matrix_root / "sample_level_ijc_matrix.tsv",
        all_sample_ids,
    )
    sjc = _load_sample_matrix(
        matrix_root / "sample_level_sjc_matrix.tsv",
        all_sample_ids,
    )
    rmats_root = project / "Setd2_Splicing_genomebam" / "star_salmon" / "rmats"
    comparison_ids = sorted(
        {comparison for mapping in comparisons.values() for comparison in mapping.values()}
    )
    rmats = {
        comparison: _load_rmats_comparison(rmats_root, comparison, args.rmats_mode)
        for comparison in comparison_ids
    }

    permutation_rows: list[dict[str, object]] = []
    null_frames: list[pd.DataFrame] = []
    model_frames: list[pd.DataFrame] = []
    for context, mapping in comparisons.items():
        mut = rmats[mapping["MUT"]]
        for inhibitor_name in ("CLK3", "Pla"):
            print(f"Analyzing {context}: MUT vs {inhibitor_name}")
            tested, shared = _build_shared_events(mut, rmats[mapping[inhibitor_name]])
            sample_ids = _sample_ids(context, inhibitor_name)
            permutation_row, null = _permutation_analysis(
                context,
                inhibitor_name,
                tested,
                shared,
                psi,
                sample_ids,
            )
            permutation_rows.append(permutation_row)
            null_frames.append(null)
            model_frames.append(
                _fit_joint_events(
                    context,
                    inhibitor_name,
                    shared,
                    psi,
                    ijc,
                    sjc,
                    sample_ids,
                    args.equivalence_margin,
                )
            )

    permutation_summary = pd.DataFrame(permutation_rows)
    permutation_summary["naive_hypergeometric_FDR"] = bh_adjust(
        permutation_summary["naive_hypergeometric_pvalue"]
    )
    null_distribution = pd.concat(null_frames, ignore_index=True)
    models = pd.concat(model_frames, ignore_index=True)
    joint_summary = _summary_joint_models(models, args.equivalence_margin)

    permutation_summary.to_csv(
        output / "shared_control_overlap_sensitivity_summary.tsv",
        sep="\t",
        index=False,
    )
    null_distribution.to_csv(
        output / "shared_control_blocked_permutation_null.tsv",
        sep="\t",
        index=False,
    )
    models.to_csv(
        output / "joint_beta_binomial_shared_event_validation.tsv",
        sep="\t",
        index=False,
    )
    joint_summary.to_csv(
        output / "joint_shared_event_validation_summary.tsv",
        sep="\t",
        index=False,
    )
    models.loc[models["strict_LOF_compatible_core"]].to_csv(
        output / "strict_shared_control_aware_LOF_compatible_core.tsv",
        sep="\t",
        index=False,
    )
    _plot_permutation_null(
        permutation_summary,
        null_distribution,
        plots / "01_shared_control_permutation_null.png",
    )
    _plot_joint_validation(
        models,
        plots / "02_direct_MUT_inhibitor_event_validation.png",
    )
    _write_report(
        output / "SHARED_CONTROL_AWARE_INTERPRETATION.md",
        permutation_summary,
        joint_summary,
        args.equivalence_margin,
    )
    _write_html_report(
        output / "SHARED_CONTROL_AWARE_REPORT.html",
        permutation_summary,
        joint_summary,
        args.equivalence_margin,
    )
    parameters = {
        "analysis": "Setd2 shared-control-aware LOF sensitivity analysis",
        "rmats_mode": args.rmats_mode,
        "original_significant_event_rule": "rMATS FDR <= 0.05 and abs(dPSI) >= 0.10",
        "event_matching": "exact canonical event_key",
        "shared_control_correction": (
            "blocked exhaustive condition-label permutation within three "
            "replicate/batch blocks plus a joint WT/MUT/inhibitor beta-binomial model"
        ),
        "permutation_count_per_pair": 216,
        "equivalence_margin_PSI": args.equivalence_margin,
        "equivalence_test": "Welch two-one-sided test, 90% CI, BH within context/inhibitor",
        "comparisons": comparisons,
        "generated_tables": sorted(path.name for path in output.glob("*.tsv")),
        "generated_plots": sorted(path.name for path in plots.glob("*.png")),
        "generated_reports": [
            "SHARED_CONTROL_AWARE_INTERPRETATION.md",
            "SHARED_CONTROL_AWARE_REPORT.html",
        ],
        "limits": [
            "The blocked permutation is a shared-control sensitivity analysis and does not rerun rMATS under every permutation.",
            "The beta-binomial model uses rMATS IJC/SJC replicate counts; event-specific isoform lengths are constant across the fitted groups.",
            "RNA-seq similarity is compatible with partial LOF-like convergence but is not biochemical proof of LOF.",
        ],
    }
    (output / "parameters.json").write_text(
        json.dumps(parameters, indent=2),
        encoding="utf-8",
    )
    print(f"Finished: {output}")
    print("\nShared-control overlap sensitivity:")
    print(permutation_summary.to_string(index=False))
    print("\nDirect joint-model validation:")
    print(joint_summary.to_string(index=False))


if __name__ == "__main__":
    main()
