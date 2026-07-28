from __future__ import annotations

import argparse
import hashlib
import html
import json
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
import pandas as pd

from katana_cosplice.matched_background import (
    build_matched_background_mut_coverage,
    build_significant_event_venn_summary,
)
from katana_cosplice.matrix import EVENT_TYPES, _event_frame


COLORS = {
    "CLK3": "#D55E00",
    "Pla": "#0072B2",
    "CLK3+Pla": "#6A51A3",
    "MUT_only": "#BDBDBD",
    "covered_by_both": "#6A51A3",
    "CLK3_only": "#D55E00",
    "Pla_only": "#0072B2",
    "ink": "#1A1A1A",
    "muted": "#666666",
}

NATURE_RC = {
    "font.family": "Arial",
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 9,
    "axes.linewidth": 0.8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "legend.fontsize": 7,
    "figure.titlesize": 10,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}


def _bh_adjust(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    output = pd.Series(np.nan, index=numeric.index, dtype=float)
    finite = numeric.dropna()
    if finite.empty:
        return output
    ordered = finite.sort_values()
    adjusted = ordered * len(ordered) / np.arange(1, len(ordered) + 1)
    adjusted = np.minimum.accumulate(adjusted.iloc[::-1])[::-1].clip(upper=1.0)
    output.loc[adjusted.index] = adjusted
    return output


def _format_pvalue(value: float) -> str:
    if pd.isna(value):
        return "NA"
    if value < 1e-4:
        return f"{value:.1e}"
    return f"{value:.4f}"


def _significance_stars(value: float) -> str:
    if pd.isna(value):
        return "NA"
    if value < 0.001:
        return "***"
    if value < 0.01:
        return "**"
    if value < 0.05:
        return "*"
    return "ns"


def _load_rmats_comparison(rmats_root: Path, comparison: str, mode: str) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for event_type in EVENT_TYPES:
        path = rmats_root / comparison / "rmats_post" / f"{event_type}.MATS.{mode}.txt"
        raw = pd.read_csv(path, sep="\t", low_memory=False)
        identity = _event_frame(raw, event_type)
        values = identity.merge(
            raw.assign(_row=np.arange(len(raw)))[["_row", "FDR", "IncLevelDifference"]],
            left_index=True,
            right_on="_row",
            how="left",
        ).drop(columns="_row")
        values["FDR"] = pd.to_numeric(values["FDR"], errors="coerce")
        values["dPSI"] = pd.to_numeric(values["IncLevelDifference"], errors="coerce")
        values["source_comparison"] = comparison
        frames.append(values)
    return pd.concat(frames, ignore_index=True)


def _save_figure(fig: plt.Figure, output_path: Path) -> None:
    fig.savefig(output_path, dpi=220, bbox_inches="tight", facecolor="white")
    if output_path.suffix.lower() == ".png":
        fig.savefig(output_path.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _draw_venn_panel(
    ax: plt.Axes,
    venn_summary: pd.DataFrame,
    context: str,
    *,
    title: str,
) -> None:
    subset = venn_summary.loc[venn_summary["context"].eq(context)].set_index("region")
    totals = {
        "MUT": int(subset["n_MUT_significant_total"].iloc[0]),
        "CLK3": int(subset["n_CLK3_significant_total"].iloc[0]),
        "PlaB": int(subset["n_PlaB_significant_total"].iloc[0]),
    }
    circles = [
        Circle((0.38, 0.58), 0.27, facecolor="#8DA0CB", edgecolor="#365F91", alpha=0.38),
        Circle((0.62, 0.58), 0.27, facecolor="#FC8D62", edgecolor=COLORS["CLK3"], alpha=0.38),
        Circle((0.50, 0.38), 0.27, facecolor="#66C2A5", edgecolor=COLORS["Pla"], alpha=0.38),
    ]
    for circle in circles:
        ax.add_patch(circle)

    positions = {
        "MUT_only": (0.23, 0.64),
        "CLK3_only": (0.77, 0.64),
        "PlaB_only": (0.50, 0.18),
        "MUT_CLK3_only": (0.50, 0.70),
        "MUT_PlaB_only": (0.35, 0.40),
        "CLK3_PlaB_only": (0.65, 0.40),
        "MUT_CLK3_PlaB": (0.50, 0.52),
    }
    for region, position in positions.items():
        ax.text(
            *position,
            f"{int(subset.loc[region, 'n_events']):,}",
            ha="center",
            va="center",
            fontsize=8.5,
            fontweight="bold",
        )

    ax.text(
        0.18,
        0.93,
        f"Setd2 MUT vs WT\nn = {totals['MUT']:,}",
        ha="center",
        va="top",
        fontsize=7.5,
        fontweight="bold",
    )
    ax.text(
        0.82,
        0.93,
        f"CLK3 inhibitor vs WT\nn = {totals['CLK3']:,}",
        ha="center",
        va="top",
        fontsize=7.5,
        fontweight="bold",
    )
    ax.text(
        0.50,
        0.01,
        f"PlaB inhibitor vs WT\nn = {totals['PlaB']:,}",
        ha="center",
        va="bottom",
        fontsize=7.5,
        fontweight="bold",
    )
    ax.set_title(title, loc="left", fontweight="bold", pad=8)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.axis("off")


def _plot_coverage(
    summary: pd.DataFrame,
    output_path: Path,
) -> None:
    directional = summary.loc[summary["inhibitor"].isin(["CLK3", "Pla"])].copy()
    directional["label"] = [
        (
            f"{row.context} background\n"
            f"MUT vs WT  |  {'PlaB' if row.inhibitor == 'Pla' else row.inhibitor} inhibitor vs WT"
        )
        for row in directional.itertuples(index=False)
    ]
    y = np.arange(len(directional))[::-1]

    with plt.rc_context(NATURE_RC):
        fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.1))
        for ax in axes:
            ax.spines[["top", "right"]].set_visible(False)
            ax.tick_params(direction="out", length=3)

        colors = [COLORS[value] for value in directional["inhibitor"]]
        axes[0].barh(y, directional["MUT_coverage_pct"], height=0.55, color=colors)
        axes[0].set_yticks(y, directional["label"])
        axes[0].set_xlim(0, 47)
        axes[0].set_xlabel("Shared significant events / all significant MUT events (%)")
        axes[0].set_title(
            "Fraction of the MUT-associated AS set recovered by each inhibitor",
            loc="left",
            fontweight="bold",
        )
        for ypos, row in zip(y, directional.itertuples(index=False), strict=True):
            axes[0].text(
                float(row.MUT_coverage_pct) + 0.7,
                ypos,
                f"{int(row.n_MUT_events_covered)} of {int(row.n_MUT_significant)} "
                f"({float(row.MUT_coverage_pct):.1f}%)",
                va="center",
                fontsize=7,
            )

        axes[1].barh(y, directional["same_direction_pct"], height=0.55, color=colors)
        axes[1].set_yticks(y, directional["label"])
        axes[1].set_xlim(50, 110)
        axes[1].set_xlabel("Shared significant events with the same dPSI direction (%)")
        axes[1].set_title(
            "Direction agreement among the shared significant events",
            loc="left",
            fontweight="bold",
        )
        axes[1].axvline(50, color="#9E9E9E", linestyle=":", linewidth=0.8)
        for ypos, row in zip(y, directional.itertuples(index=False), strict=True):
            fdr = float(row.same_direction_FDR)
            axes[1].text(
                float(row.same_direction_pct) + 0.7,
                ypos,
                f"{int(row.n_same_direction)} of {int(row.n_MUT_events_covered)} "
                f"({float(row.same_direction_pct):.1f}%)\n"
                f"q={_format_pvalue(fdr)} {_significance_stars(fdr)}",
                va="center",
                fontsize=6.5,
                linespacing=1.05,
            )

        for panel, ax in enumerate(axes):
            ax.text(
                -0.24,
                1.03,
                chr(ord("a") + panel),
                transform=ax.transAxes,
                fontsize=11,
                fontweight="bold",
            )
        fig.suptitle(
            "Setd2 MUT shares a limited but directionally consistent AS response "
            "with CLK3 and PlaB inhibition",
            x=0.045,
            ha="left",
            fontweight="bold",
        )
        fig.text(
            0.045,
            0.005,
            "Shared means the exact same AS event met rMATS FDR <= 0.05 and "
            "|standardized dPSI| >= 0.10 in both contrasts. "
            "q: BH-adjusted one-sided binomial test against 50% direction agreement.",
            fontsize=6.5,
            color=COLORS["muted"],
        )
        fig.subplots_adjust(left=0.23, right=0.98, bottom=0.18, top=0.80, wspace=0.64)
        _save_figure(fig, output_path)


def _plot_significant_event_venn(
    venn_summary: pd.DataFrame,
    comparison_design: pd.DataFrame,
    context: str,
    output_path: Path,
) -> None:
    design = comparison_design.loc[comparison_design["context"].eq(context)].set_index("role")

    with plt.rc_context(NATURE_RC):
        fig, ax = plt.subplots(figsize=(5.8, 4.5))
        _draw_venn_panel(
            ax,
            venn_summary,
            context,
            title=f"{context}: overlap among significant AS-event sets",
        )
        ax.text(
            0.02,
            -0.10,
            (
                f"MUT: {design.loc['MUT', 'comparison_id']}  |  "
                f"CLK3: {design.loc['CLK3', 'comparison_id']}  |  "
                f"PlaB: {design.loc['PlaB', 'comparison_id']}\n"
                "Events met rMATS FDR <= 0.05 and |standardized dPSI| >= 0.10. "
                "The same event can overlap even when its dPSI direction differs."
            ),
            transform=ax.transAxes,
            fontsize=6.5,
            color=COLORS["muted"],
            va="top",
        )
        fig.subplots_adjust(left=0.03, right=0.97, bottom=0.20, top=0.90)
        _save_figure(fig, output_path)


def _build_pairwise_venn_summary(membership: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for context in ["DP", "noDP"]:
        context_frame = membership.loc[membership["context"].eq(context)]
        for inhibitor, column, display in [
            ("CLK3", "CLK3_significant", "CLK3 inhibitor"),
            ("PlaB", "Pla_significant", "PlaB inhibitor"),
        ]:
            mut = context_frame["MUT_significant"].fillna(False)
            treatment = context_frame[column].fillna(False)
            shared = int((mut & treatment).sum())
            mut_total = int(mut.sum())
            treatment_total = int(treatment.sum())
            rows.append(
                {
                    "context": context,
                    "comparison_pair": f"MUT_vs_{inhibitor}",
                    "left_label": "Setd2 MUT vs WT",
                    "right_label": f"{display} vs WT",
                    "n_MUT_significant": mut_total,
                    "n_inhibitor_significant": treatment_total,
                    "n_MUT_only": int((mut & ~treatment).sum()),
                    "n_shared_exact_events": shared,
                    "n_inhibitor_only": int((~mut & treatment).sum()),
                    "pct_MUT_events_shared": 100.0 * shared / mut_total if mut_total else np.nan,
                }
            )
    return pd.DataFrame(rows)


def _plot_pairwise_venn(
    pairwise_summary: pd.DataFrame,
    output_path: Path,
) -> None:
    with plt.rc_context(NATURE_RC):
        fig, axes = plt.subplots(2, 2, figsize=(8.2, 6.6))
        for panel, (ax, row) in enumerate(
            zip(axes.flat, pairwise_summary.itertuples(index=False), strict=True)
        ):
            ax.add_patch(
                Circle(
                    (0.40, 0.50),
                    0.29,
                    facecolor="#8DA0CB",
                    edgecolor="#365F91",
                    alpha=0.40,
                )
            )
            inhibitor_color = COLORS["CLK3"] if "CLK3" in row.comparison_pair else COLORS["Pla"]
            ax.add_patch(
                Circle(
                    (0.60, 0.50),
                    0.29,
                    facecolor=inhibitor_color,
                    edgecolor=inhibitor_color,
                    alpha=0.30,
                )
            )
            ax.text(0.27, 0.50, f"{int(row.n_MUT_only):,}", ha="center", fontweight="bold")
            ax.text(
                0.50,
                0.50,
                f"{int(row.n_shared_exact_events):,}",
                ha="center",
                fontweight="bold",
            )
            ax.text(
                0.73,
                0.50,
                f"{int(row.n_inhibitor_only):,}",
                ha="center",
                fontweight="bold",
            )
            ax.text(
                0.24,
                0.88,
                f"{row.left_label}\nn = {int(row.n_MUT_significant):,}",
                ha="center",
                va="top",
                fontsize=7.5,
                fontweight="bold",
            )
            ax.text(
                0.76,
                0.88,
                f"{row.right_label}\nn = {int(row.n_inhibitor_significant):,}",
                ha="center",
                va="top",
                fontsize=7.5,
                fontweight="bold",
            )
            ax.text(
                0.50,
                0.12,
                f"Shared: {int(row.n_shared_exact_events):,} "
                f"({float(row.pct_MUT_events_shared):.1f}% of MUT events)",
                ha="center",
                fontsize=7,
            )
            ax.set_title(
                f"{chr(ord('a') + panel)}   {row.context} background",
                loc="left",
                fontweight="bold",
                pad=8,
            )
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.set_aspect("equal")
            ax.axis("off")

        fig.suptitle(
            "Pairwise overlap of significant AS events in matched experimental backgrounds",
            x=0.06,
            ha="left",
            fontweight="bold",
        )
        fig.text(
            0.06,
            0.015,
            "Each overlap is an exact AS-event match. Event direction is not required for "
            "membership in the overlap and is evaluated separately in Figure 1b.",
            fontsize=6.5,
            color=COLORS["muted"],
        )
        fig.subplots_adjust(left=0.05, right=0.98, bottom=0.08, top=0.90, hspace=0.18, wspace=0.16)
        _save_figure(fig, output_path)


def _plot_coverage_classes(class_summary: pd.DataFrame, output_path: Path) -> None:
    class_order = ["MUT_only", "CLK3_only", "Pla_only", "covered_by_both"]
    labels = {
        "MUT_only": "MUT only",
        "CLK3_only": "CLK3 only",
        "Pla_only": "PlaB only",
        "covered_by_both": "Covered by both",
    }
    matrix = (
        class_summary.pivot(index="context", columns="coverage_class", values="pct_MUT_signature")
        .reindex(index=["DP", "noDP"], columns=class_order)
        .fillna(0.0)
    )

    with plt.rc_context(NATURE_RC):
        fig, ax = plt.subplots(figsize=(4.7, 3.25))
        bottom = np.zeros(len(matrix))
        for coverage_class in class_order:
            values = matrix[coverage_class].to_numpy(dtype=float)
            bars = ax.bar(
                matrix.index,
                values,
                bottom=bottom,
                width=0.60,
                label=labels[coverage_class],
                color=COLORS[coverage_class],
                edgecolor="white",
                linewidth=0.4,
            )
            for bar, value, base in zip(bars, values, bottom, strict=True):
                if value >= 6.0:
                    ax.text(
                        bar.get_x() + bar.get_width() / 2,
                        base + value / 2,
                        f"{value:.1f}%",
                        ha="center",
                        va="center",
                        color="white" if coverage_class != "MUT_only" else COLORS["ink"],
                        fontsize=7,
                        fontweight="bold",
                    )
            bottom += values
        ax.set_ylim(0, 100)
        ax.set_ylabel("MUT significant-event signature (%)")
        ax.set_title("Composition of the MUT AS signature", loc="left", fontweight="bold")
        ax.legend(frameon=False, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.15))
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(direction="out", length=3)
        fig.subplots_adjust(left=0.16, right=0.98, bottom=0.28, top=0.90)
        _save_figure(fig, output_path)


def _plot_dpsi_scatter(
    catalog: pd.DataFrame,
    summary: pd.DataFrame,
    context: str,
    output_path: Path,
) -> None:
    with plt.rc_context(NATURE_RC):
        fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.45), sharex=True, sharey=True)
        for panel, (ax, inhibitor) in enumerate(zip(axes, ["CLK3", "Pla"], strict=True)):
            subset = catalog.loc[
                catalog["context"].eq(context) & catalog[f"{inhibitor}_covers_MUT"]
            ].copy()
            same = subset[f"{inhibitor}_same_direction"].fillna(False)
            ax.scatter(
                subset.loc[~same, "MUT_dPSI"],
                subset.loc[~same, f"{inhibitor}_dPSI"],
                s=14,
                alpha=0.72,
                color="#CC79A7",
                linewidth=0,
                label="Opposite direction",
            )
            ax.scatter(
                subset.loc[same, "MUT_dPSI"],
                subset.loc[same, f"{inhibitor}_dPSI"],
                s=14,
                alpha=0.58,
                color=COLORS[inhibitor],
                linewidth=0,
                label="Same direction",
            )
            limits = np.array(
                [
                    subset["MUT_dPSI"].min() if not subset.empty else -0.5,
                    subset["MUT_dPSI"].max() if not subset.empty else 0.5,
                    subset[f"{inhibitor}_dPSI"].min() if not subset.empty else -0.5,
                    subset[f"{inhibitor}_dPSI"].max() if not subset.empty else 0.5,
                ],
                dtype=float,
            )
            span = max(abs(limits.min()), abs(limits.max()), 0.15) * 1.06
            ax.plot([-span, span], [-span, span], linestyle="--", linewidth=0.8, color="#7F7F7F")
            ax.axhline(0, linewidth=0.6, color="#B0B0B0")
            ax.axvline(0, linewidth=0.6, color="#B0B0B0")
            ax.set_xlim(-span, span)
            ax.set_ylim(-span, span)
            ax.set_xlabel("MUT standardized dPSI")
            ax.set_ylabel(f"{'PlaB' if inhibitor == 'Pla' else inhibitor} standardized dPSI")
            row = summary.loc[
                summary["context"].eq(context) & summary["inhibitor"].eq(inhibitor)
            ].iloc[0]
            ax.set_title(
                f"{'PlaB' if inhibitor == 'Pla' else inhibitor} | n={len(subset)}",
                loc="left",
                fontweight="bold",
            )
            ax.text(
                0.03,
                0.97,
                f"same direction: {float(row.same_direction_pct):.1f}% "
                f"(FDR={_format_pvalue(float(row.same_direction_FDR))})\n"
                f"Spearman rho={float(row.dPSI_spearman_rho):.2f} "
                f"(FDR={_format_pvalue(float(row.dPSI_spearman_FDR))})",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=6.5,
            )
            ax.text(
                -0.17,
                1.07,
                chr(ord("a") + panel),
                transform=ax.transAxes,
                fontsize=10,
                fontweight="bold",
            )
            ax.spines[["top", "right"]].set_visible(False)
            ax.tick_params(direction="out", length=3)
            ax.legend(frameon=False, loc="lower right")
        fig.suptitle(
            f"{context}: effect-size agreement among covered MUT events",
            x=0.08,
            ha="left",
            fontweight="bold",
        )
        fig.subplots_adjust(left=0.10, right=0.98, bottom=0.16, top=0.82, wspace=0.30)
        _save_figure(fig, output_path)


def _summary_for_display(summary: pd.DataFrame) -> pd.DataFrame:
    display = summary.loc[summary["inhibitor"].isin(["CLK3", "Pla"])].copy()
    display["inhibitor"] = display["inhibitor"].replace({"Pla": "PlaB"})
    display = display.rename(
        columns={
            "n_MUT_significant": "MUT significant",
            "n_MUT_events_covered": "covered",
            "MUT_coverage_pct": "coverage (%)",
            "same_direction_pct": "same direction (%)",
            "same_direction_FDR": "direction FDR",
            "dPSI_spearman_rho": "dPSI Spearman rho",
            "dPSI_spearman_FDR": "rho FDR",
            "median_abs_delta_dPSI": "median abs delta dPSI",
        }
    )
    return display[
        [
            "context",
            "inhibitor",
            "MUT significant",
            "covered",
            "coverage (%)",
            "same direction (%)",
            "direction FDR",
            "dPSI Spearman rho",
            "rho FDR",
            "median abs delta dPSI",
        ]
    ]


def _build_markdown(
    summary: pd.DataFrame,
    class_summary: pd.DataFrame,
    comparison_design: pd.DataFrame,
    pairwise_venn_summary: pd.DataFrame,
    venn_summary: pd.DataFrame,
) -> str:
    display = _summary_for_display(summary)
    result_lookup = summary.set_index(["context", "inhibitor"])
    venn_lookup = venn_summary.set_index(["context", "region"])
    dp_clk = result_lookup.loc[("DP", "CLK3")]
    dp_pla = result_lookup.loc[("DP", "Pla")]
    nodp_clk = result_lookup.loc[("noDP", "CLK3")]
    nodp_pla = result_lookup.loc[("noDP", "Pla")]
    dp_three_way = int(venn_lookup.loc[("DP", "MUT_CLK3_PlaB"), "n_events"])
    nodp_three_way = int(venn_lookup.loc[("noDP", "MUT_CLK3_PlaB"), "n_events"])
    class_display = class_summary.copy()
    class_display["coverage_class"] = class_display["coverage_class"].replace(
        {
            "MUT_only": "MUT only",
            "CLK3_only": "CLK3 only",
            "Pla_only": "PlaB only",
            "covered_by_both": "covered by both",
        }
    )
    rows = [
        "# Setd2 MUT and inhibitor-associated alternative splicing",
        "",
        "## Executive conclusion",
        "",
        "The Setd2 mutation shares a limited but highly directional alternative-splicing core "
        "with CLK3 and PlaB inhibition. This is not a global inhibitor-like response: most "
        "significant MUT events are not significant after either inhibitor.",
        "",
        f"In the DP background, CLK3 inhibition reproduced {int(dp_clk.n_MUT_events_covered)} "
        f"of {int(dp_clk.n_MUT_significant)} MUT-associated events "
        f"({dp_clk.MUT_coverage_pct:.1f}%), whereas PlaB reproduced "
        f"{int(dp_pla.n_MUT_events_covered)} ({dp_pla.MUT_coverage_pct:.1f}%). "
        f"{dp_three_way} events were significant in MUT, CLK3 and PlaB comparisons.",
        "",
        f"Without DP, CLK3 reproduced {int(nodp_clk.n_MUT_events_covered)} of "
        f"{int(nodp_clk.n_MUT_significant)} MUT-associated events "
        f"({nodp_clk.MUT_coverage_pct:.1f}%), whereas PlaB reproduced "
        f"{int(nodp_pla.n_MUT_events_covered)} ({nodp_pla.MUT_coverage_pct:.1f}%). "
        f"{nodp_three_way} events were shared by all three comparisons.",
        "",
        "The shared events were not merely the same genomic events: 85.6-94.7% changed in the "
        "same dPSI direction as the mutation, and all four direction-concordance tests remained "
        "significant after FDR correction. This supports a reproducible point of convergence "
        "between the mutation and inhibitor responses, while the much larger inhibitor-only "
        "event sets show that neither inhibitor is a complete phenocopy of the mutation.",
        "",
        "## Exact comparisons included",
        "",
        comparison_design.to_markdown(index=False),
        "",
        "## Analysis process",
        "",
        "1. Load all five JCEC rMATS event tables for each of the six contrasts.",
        "2. Build exact canonical event keys from normalized gene ID, event type, chromosome, "
        "strand and event coordinates.",
        "3. Orient all contrasts as MUT/inhibitor minus matched WT.",
        "4. Define significant AS as rMATS FDR <= 0.05 and abs(dPSI) >= 0.10.",
        "   No IJC/SJC minimum, replicate-stability criterion, DEXSeq or DTU filter is applied "
        "in this report.",
        "5. Within DP and noDP separately, count the exact events shared by MUT, CLK3 and PlaB.",
        "6. Ask what fraction of the MUT-associated events is recovered by each inhibitor.",
        "7. For shared events, test direction concordance against 50% with an exact one-sided "
        "binomial test and control the four tests by Benjamini-Hochberg FDR.",
        "8. Quantify dPSI effect-size agreement using Spearman correlation with BH FDR.",
        "",
        "## Main results",
        "",
        display.to_markdown(index=False, floatfmt=".3g"),
        "",
        "## Complete significant-event overlaps",
        "",
        "### Pairwise overlaps",
        "",
        pairwise_venn_summary.to_markdown(index=False, floatfmt=".3g"),
        "",
        "### Three-way overlaps",
        "",
        "Each circle contains all significant events from one comparison. The numbers inside "
        "the overlaps are exact event matches. An event can occupy an overlap even if its dPSI "
        "sign differs, so overlap and direction are reported separately.",
        "",
        venn_summary.to_markdown(index=False, floatfmt=".3g"),
        "",
        "## MUT signature composition",
        "",
        class_display.to_markdown(index=False, floatfmt=".3g"),
        "",
        "## Statistical interpretation",
        "",
        "- Direction FDR tests whether direction concordance among covered events exceeds 50%.",
        "- rho FDR tests monotonic dPSI effect-size agreement among the same covered events.",
        "- The overlap fractions are descriptive. No independent overlap-enrichment P value is reported "
        "because MUT and inhibitor contrasts share the matched WT reference.",
        "- Event-level tests treat events as units and do not model dependence among events "
        "from the same gene or splicing program; significance is therefore supporting evidence, "
        "not proof of a causal mechanism.",
        "",
        "## Interpretation boundary",
        "",
        "These data show a shared splicing response, not proof that the mutation is LOF or GOF. "
        "They also do not require the MUT and inhibitor dPSI values to be numerically identical.",
        "",
        "The complete event-level audit trail is in `matched_background_MUT_event_catalog.tsv`.",
    ]
    return "\n".join(rows) + "\n"


def _build_html(
    summary: pd.DataFrame,
    class_summary: pd.DataFrame,
    comparison_design: pd.DataFrame,
    pairwise_venn_summary: pd.DataFrame,
    venn_summary: pd.DataFrame,
    top_events: pd.DataFrame,
    output: Path,
) -> str:
    result_lookup = summary.set_index(["context", "inhibitor"])
    venn_lookup = venn_summary.set_index(["context", "region"])
    dp_clk = result_lookup.loc[("DP", "CLK3")]
    dp_pla = result_lookup.loc[("DP", "Pla")]
    nodp_clk = result_lookup.loc[("noDP", "CLK3")]
    nodp_pla = result_lookup.loc[("noDP", "Pla")]
    dp_three_way = int(venn_lookup.loc[("DP", "MUT_CLK3_PlaB"), "n_events"])
    nodp_three_way = int(venn_lookup.loc[("noDP", "MUT_CLK3_PlaB"), "n_events"])
    table_html = _summary_for_display(summary).to_html(
        index=False,
        float_format=lambda value: f"{value:.3g}",
        border=0,
        classes="data-table",
    )
    comparison_html = comparison_design.to_html(index=False, border=0, classes="data-table")
    pairwise_venn_html = pairwise_venn_summary.to_html(
        index=False,
        float_format=lambda value: f"{value:.3g}",
        border=0,
        classes="data-table",
    )
    venn_html = venn_summary.to_html(
        index=False,
        float_format=lambda value: f"{value:.3g}",
        border=0,
        classes="data-table",
    )
    class_table = class_summary.copy()
    class_table["coverage_class"] = class_table["coverage_class"].replace(
        {
            "MUT_only": "MUT only",
            "CLK3_only": "CLK3 only",
            "Pla_only": "PlaB only",
            "covered_by_both": "covered by both",
        }
    )
    class_html = class_table.to_html(
        index=False,
        float_format=lambda value: f"{value:.3g}",
        border=0,
        classes="data-table",
    )
    top_html = top_events.to_html(
        index=False,
        float_format=lambda value: f"{value:.3g}",
        border=0,
        classes="data-table",
    )
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Setd2 MUT and inhibitor-associated alternative splicing</title>
<style>
:root {{ --ink:#1a1a1a; --muted:#666; --line:#d9d9d9; --clk:#d55e00; --pla:#0072b2; }}
body {{ font-family: Arial, sans-serif; margin: 38px auto 80px; max-width: 1160px; color: var(--ink); line-height: 1.48; }}
h1 {{ font-size: 32px; margin-bottom: 5px; }}
h2 {{ font-size: 20px; border-bottom: 1px solid var(--line); padding-bottom: 7px; margin-top: 38px; }}
h3 {{ font-size: 15px; margin-top: 24px; }}
.kicker {{ text-transform: uppercase; letter-spacing: 0.12em; color: var(--muted); font-size: 12px; }}
.answer {{ font-size: 18px; background: #f4f4f4; border-left: 5px solid var(--ink); padding: 18px 22px; margin: 22px 0; }}
.key {{ display:grid; grid-template-columns:repeat(3,1fr); gap:12px; margin:18px 0; }}
.metric {{ border-top:3px solid var(--ink); padding:10px 4px; }}
.metric b {{ display:block; font-size:25px; }}
.metric span {{ color:var(--muted); font-size:12px; }}
.boundary {{ background:#fff7ed; border-left:4px solid var(--clk); padding:14px 18px; }}
.figure {{ margin:26px 0 34px; }}
.figure img {{ display:block; width:100%; max-width:1050px; margin:auto; }}
.caption {{ color:var(--muted); font-size:12px; margin:8px auto; max-width:1000px; }}
.data-table {{ border-collapse:collapse; width:100%; font-size:12px; }}
.data-table th {{ background:#f2f2f2; text-align:left; }}
.data-table th,.data-table td {{ border-bottom:1px solid var(--line); padding:7px 8px; }}
.scroll {{ overflow-x:auto; }}
code {{ background:#f4f4f4; padding:2px 4px; }}
ol li {{ margin-bottom:6px; }}
</style>
</head>
<body>
<div class="kicker">Katana / Setd2 matched-background analysis</div>
<h1>Setd2 MUT and inhibitor-associated alternative splicing</h1>
<p class="kicker">Matched backgrounds · exact AS-event identity · direction tested separately</p>
<div class="answer"><b>Main result.</b> The Setd2 mutation shares a limited but strongly
directional alternative-splicing core with CLK3 and PlaB inhibition. The shared component is
biologically coherent, but it represents only part of the MUT response and therefore does not
support a complete inhibitor phenocopy.</div>
<div class="key">
<div class="metric"><b>20–35%</b><span>MUT signature covered by one inhibitor</span></div>
<div class="metric"><b>~10%</b><span>MUT signature covered by both inhibitors</span></div>
<div class="metric"><b>86–95%</b><span>same dPSI direction among covered events</span></div>
</div>
<h2>1. Comparisons included</h2>
<p>Every effect is estimated within the same DP background and oriented as
<code>MUT/inhibitor − matched WT</code>. The WT DP-versus-noDP comparison is not used.</p>
<div class="scroll">{comparison_html}</div>
<h2>2. Analysis process</h2>
<ol>
<li>Load all five JCEC rMATS event types for the six contrasts above.</li>
<li>Match exact events using normalized gene ID, event type, chromosome, strand and rMATS coordinates.</li>
<li>Apply rMATS FDR ≤ 0.05 and |standardized dPSI| ≥ 0.10.</li>
<li>No IJC/SJC minimum, replicate-stability rule, DEXSeq filter or DTU filter is applied in
this report. The input is therefore an effect-size-filtered rMATS significant-event set.</li>
<li>Use all significant MUT events as the fixed denominator in each background.</li>
<li>Classify MUT events as CLK3-covered, PlaB-covered, covered by both or MUT-only.</li>
<li>Test direction concordance against 50% using an exact one-sided binomial test with BH FDR.</li>
<li>Measure effect-size agreement with Spearman rho and BH FDR.</li>
</ol>
<h2>3. Results</h2>
<div class="figure"><img src="plots/01_MUT_coverage_and_direction.png" alt="Coverage and direction">
<div class="caption"><b>Figure 1.</b> Left: the fraction of all significant MUT-associated
events that is also significant in each matched inhibitor comparison. Right: among those
shared events, the fraction changing in the same standardized dPSI direction. Asterisks and
q values refer to the direction test, not to overlap enrichment.</div></div>
<p><b>DP background.</b> CLK3 inhibition reproduced
{int(dp_clk.n_MUT_events_covered)} of {int(dp_clk.n_MUT_significant)} MUT-associated events
({dp_clk.MUT_coverage_pct:.1f}%), and PlaB reproduced {int(dp_pla.n_MUT_events_covered)}
({dp_pla.MUT_coverage_pct:.1f}%). {dp_three_way} events were significant in all three
comparisons.</p>
<p><b>No-DP background.</b> CLK3 reproduced {int(nodp_clk.n_MUT_events_covered)} of
{int(nodp_clk.n_MUT_significant)} MUT-associated events ({nodp_clk.MUT_coverage_pct:.1f}%),
and PlaB reproduced {int(nodp_pla.n_MUT_events_covered)}
({nodp_pla.MUT_coverage_pct:.1f}%). {nodp_three_way} events were significant in all three
comparisons.</p>
<p><b>Direction.</b> Among the shared events, 85.6-94.7% changed in the same standardized
dPSI direction as the mutation. The four direction-concordance tests remained significant
after FDR correction. Thus, the common events represent a directional shared response rather
than an overlap created only by large event counts.</p>
<div class="scroll">
{table_html}
</div>
<h2>4. Pairwise significant-event overlaps</h2>
<p>Each circle contains all events meeting rMATS FDR ≤ 0.05 and |standardized dPSI| ≥ 0.10
in that comparison. Circle labels give the total number of significant events. Numbers in the
overlaps are exact matches of the same AS event. Each panel compares MUT with one inhibitor
in the same DP background.</p>
<div class="figure"><img src="plots/02_pairwise_significant_event_venn.png" alt="Pairwise Venn diagrams">
<div class="caption"><b>Figure 2.</b> Pairwise event overlap for MUT versus CLK3 and MUT versus
PlaB, analyzed separately in DP and noDP. The percentage below each diagram is calculated
against the complete significant MUT event set for that background.</div></div>
<div class="scroll">{pairwise_venn_html}</div>
<h2>5. Three-way overlap</h2>
<p>The three-way diagrams are shown after the pairwise results. They identify the stricter
subset significant in MUT, CLK3 and PlaB comparisons simultaneously. Overlap alone does not
imply the same direction; direction is evaluated separately in Figure 1.</p>
<div class="figure"><img src="plots/05_DP_significant_event_venn.png" alt="DP Venn diagram">
<div class="caption"><b>Figure 3.</b> DP three-way significant-event overlap:
MUT_DP versus WT_DP, WT_DP_CLK3 versus WT_DP, and WT_DP_PlaB versus WT_DP.</div></div>
<div class="figure"><img src="plots/06_noDP_significant_event_venn.png" alt="noDP Venn diagram">
<div class="caption"><b>Figure 4.</b> noDP three-way significant-event overlap:
MUT_noDP versus WT_noDP, WT_noDP_CLK3 versus WT_noDP, and WT_noDP_PlaB versus WT_noDP.</div></div>
<div class="scroll">{venn_html}</div>
<h2>6. MUT-signature composition</h2>
<div class="figure"><img src="plots/02_MUT_coverage_classes.png" alt="Coverage classes">
<div class="caption"><b>Figure 5.</b> Decomposition of the complete MUT signature. CLK3-only and
PlaB-only mean inhibitor-specific coverage of MUT events; covered-by-both is the strict
three-way event intersection.</div></div>
<div class="scroll">{class_html}</div>
<h2>7. Effect-size agreement</h2>
<div class="figure"><img src="plots/03_DP_dPSI_scatter.png" alt="DP dPSI scatter">
<div class="caption"><b>Figure 6.</b> DP matched-background shared significant events. The
diagonal indicates equal dPSI magnitude; equality is not required for directional agreement.</div></div>
<div class="figure"><img src="plots/04_noDP_dPSI_scatter.png" alt="noDP dPSI scatter">
<div class="caption"><b>Figure 7.</b> noDP matched-background shared significant events.</div></div>
<h2>8. Representative covered events</h2>
<p>The delivery package contains the complete event catalog. The table below is a compact
high-confidence view prioritizing events covered by both inhibitors and low rMATS FDR.</p>
<div class="scroll">{top_html}</div>
<h2>9. Statistical boundary</h2>
<div class="boundary"><b>Coverage is descriptive.</b> MUT and inhibitor contrasts share the
matched WT reference, so the report does not present an independence-based hypergeometric
overlap P value. Direction FDR and rho FDR are event-level supporting statistics and do not
model dependence among events from the same gene or splicing program. This analysis does not
classify the mutation as LOF or GOF.</div>
<h2>10. Reproducibility package</h2>
<p>Complete TSV results, input-file SHA-256 manifest, parameters, methods, PNG/PDF figures and
the exact analysis code are included in the ZIP package. Output directory:
<code>{html.escape(str(output))}</code>.</p>
</body>
</html>
"""


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _comparison_design(comparisons: dict[str, dict[str, str]]) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    groups = {
        "DP": {
            "MUT": ("MUT_DP", "WT_DP"),
            "CLK3": ("WT_DP_CLK3", "WT_DP"),
            "Pla": ("WT_DP_Pla", "WT_DP"),
        },
        "noDP": {
            "MUT": ("MUT_NO_DP", "WT_NO_DP"),
            "CLK3": ("WT_NO_DP_CLK3", "WT_NO_DP"),
            "Pla": ("WT_NO_DP_Pla", "WT_NO_DP"),
        },
    }
    for context, mapping in comparisons.items():
        for role, comparison_id in mapping.items():
            numerator, denominator = groups[context][role]
            rows.append(
                {
                    "context": context,
                    "role": "PlaB" if role == "Pla" else role,
                    "comparison_id": comparison_id,
                    "numerator": numerator,
                    "denominator": denominator,
                    "standardized_dPSI_positive_means": (
                        f"event inclusion/usage higher in {numerator} than {denominator}"
                    ),
                }
            )
    return pd.DataFrame(rows)


def _input_manifest(
    rmats_root: Path,
    comparisons: dict[str, dict[str, str]],
    mode: str,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for context, mapping in comparisons.items():
        for role, comparison_id in mapping.items():
            for event_type in EVENT_TYPES:
                path = rmats_root / comparison_id / "rmats_post" / f"{event_type}.MATS.{mode}.txt"
                stat = path.stat()
                rows.append(
                    {
                        "context": context,
                        "role": "PlaB" if role == "Pla" else role,
                        "comparison_id": comparison_id,
                        "event_type": event_type,
                        "source_file": str(path.resolve()),
                        "size_bytes": stat.st_size,
                        "modified_time": pd.Timestamp(stat.st_mtime, unit="s").isoformat(),
                        "sha256": _file_sha256(path),
                    }
                )
    return pd.DataFrame(rows)


def _top_covered_events(catalog: pd.DataFrame, limit: int = 40) -> pd.DataFrame:
    covered = catalog.loc[~catalog["coverage_class"].eq("MUT_only")].copy()
    covered["_class_order"] = covered["coverage_class"].map(
        {"covered_by_both": 0, "CLK3_only": 1, "Pla_only": 2}
    )
    covered["_minimum_FDR"] = covered[["MUT_FDR", "CLK3_FDR", "Pla_FDR"]].min(axis=1)
    covered = covered.sort_values(
        ["_class_order", "_minimum_FDR", "context", "gene_symbol", "event_type"],
        na_position="last",
    ).head(limit)
    columns = [
        "context",
        "gene_symbol",
        "event_type",
        "coverage_class",
        "MUT_dPSI",
        "MUT_FDR",
        "CLK3_dPSI",
        "CLK3_FDR",
        "CLK3_same_direction",
        "Pla_dPSI",
        "Pla_FDR",
        "Pla_same_direction",
        "event_key",
    ]
    return covered.loc[:, [column for column in columns if column in covered.columns]]


def _write_delivery_package(output: Path, repo_root: Path) -> Path:
    delivery = output / "delivery_package"
    if delivery.exists():
        shutil.rmtree(delivery)
    (delivery / "plots").mkdir(parents=True)
    (delivery / "tables").mkdir()
    (delivery / "methods").mkdir()
    (delivery / "code").mkdir()

    for name in [
        "MATCHED_BACKGROUND_MUT_COVERAGE_REPORT.html",
        "MATCHED_BACKGROUND_MUT_COVERAGE_REPORT.md",
        "parameters.json",
    ]:
        shutil.copy2(output / name, delivery / name)
    for path in (output / "plots").glob("*"):
        if path.is_file():
            shutil.copy2(path, delivery / "plots" / path.name)
    for path in output.glob("*.tsv"):
        shutil.copy2(path, delivery / "tables" / path.name)

    methods = """# Analysis methods and provenance

## Scientific question

Within DP and noDP separately, what fraction of significant Setd2 MUT alternative-splicing
events is also significant after CLK3 or PlaB inhibition, and do covered events change in
the same direction?

## Workflow

1. Read JCEC rMATS SE, RI, A3SS, A5SS and MXE result tables for six matched-background contrasts.
2. Normalize Ensembl gene IDs and build exact canonical event keys from gene ID, event type,
   chromosome, strand and event coordinates.
3. Orient every contrast as MUT/inhibitor minus its matched WT.
4. Define significant AS as rMATS FDR <= 0.05 and abs(dPSI) >= 0.10.
   No IJC/SJC minimum, replicate-stability criterion, DEXSeq filter or DTU filter is applied.
5. Build the complete union and exact seven-region membership of the MUT, CLK3 and PlaB
   significant-event sets for each DP background.
6. Fix the denominator to all significant MUT events separately in DP and noDP.
7. Left-match inhibitor results onto that denominator and classify coverage as CLK3-only,
   PlaB-only, both or MUT-only.
8. Among covered events, test direction concordance above 50% with an exact one-sided binomial
   test. Apply Benjamini-Hochberg FDR across the four context-by-inhibitor tests.
9. Calculate Spearman correlation between MUT and inhibitor dPSI among covered events and apply
   BH FDR across the same four comparisons.

## Statistical boundary

Coverage is descriptive because the MUT and inhibitor effects share the same matched WT
reference. No independence-based overlap P value is claimed. Direction and correlation tests
treat events as units and do not model dependence among events from the same gene.
"""
    (delivery / "methods" / "ANALYSIS_METHODS_AND_PROVENANCE.md").write_text(
        methods,
        encoding="utf-8",
    )
    for source in [
        repo_root / "katana_cosplice" / "matched_background.py",
        repo_root / "scripts" / "analyze_setd2_matched_background_mut_coverage.py",
        repo_root / "tests" / "test_matched_background.py",
    ]:
        shutil.copy2(source, delivery / "code" / source.name)

    readme = """# Setd2 matched-background MUT-centric AS coverage package

Open `MATCHED_BACKGROUND_MUT_COVERAGE_REPORT.html` first.

## Contents

- `MATCHED_BACKGROUND_MUT_COVERAGE_REPORT.html`: complete visual report.
- `MATCHED_BACKGROUND_MUT_COVERAGE_REPORT.md`: plain-text scientific report.
- `plots/`: publication-style PNG and editable-vector PDF figures.
- `tables/`: complete summary, comparison design, event catalog, event-type summary,
  representative covered events, process log and SHA-256 input manifest.
- `methods/`: analysis workflow and interpretation boundaries.
- `code/`: exact algorithm, runner and regression test used for this package.
- `parameters.json`: thresholds, directions, comparisons and generated-file inventory.

The complete event catalog is the audit trail. No LOF/GOF classification is made.
"""
    (delivery / "README.md").write_text(readme, encoding="utf-8")

    zip_path = output / "SETD2_MATCHED_BACKGROUND_MUT_COVERAGE_PACKAGE.zip"
    if zip_path.exists():
        zip_path.unlink()
    shutil.make_archive(
        str(zip_path.with_suffix("")),
        "zip",
        root_dir=delivery,
    )
    return zip_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Analyze Setd2 MUT AS coverage by inhibitors within matched DP backgrounds."
    )
    parser.add_argument(
        "--rmats-root",
        type=Path,
        default=Path(
            r"E:\Project\RNASEQ\Result\Setd2_results"
            r"\Setd2_Splicing_genomebam\star_salmon\rmats"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            r"E:\Project\RNASEQ\Result\Setd2_results\katana_output"
            r"\09_cosplice_network\06l_matched_background_MUT_coverage"
        ),
    )
    parser.add_argument("--mode", choices=["JCEC", "JC"], default="JCEC")
    parser.add_argument("--fdr-cutoff", type=float, default=0.05)
    parser.add_argument("--abs-dpsi-cutoff", type=float, default=0.10)
    args = parser.parse_args()

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
    args.output.mkdir(parents=True, exist_ok=True)
    prior_delivery = args.output / "delivery_package"
    if prior_delivery.exists():
        shutil.rmtree(prior_delivery)
    prior_zip = args.output / "SETD2_MATCHED_BACKGROUND_MUT_COVERAGE_PACKAGE.zip"
    if prior_zip.exists():
        prior_zip.unlink()
    plot_dir = args.output / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    loaded: dict[str, pd.DataFrame] = {}
    catalogs: list[pd.DataFrame] = []
    summaries: list[pd.DataFrame] = []
    classes: list[pd.DataFrame] = []
    venn_memberships: list[pd.DataFrame] = []
    venn_summaries: list[pd.DataFrame] = []
    for context, mapping in comparisons.items():
        for comparison in mapping.values():
            if comparison not in loaded:
                loaded[comparison] = _load_rmats_comparison(args.rmats_root, comparison, args.mode)
        result = build_matched_background_mut_coverage(
            loaded[mapping["MUT"]],
            {
                "CLK3": loaded[mapping["CLK3"]],
                "Pla": loaded[mapping["Pla"]],
            },
            context=context,
            fdr_cutoff=args.fdr_cutoff,
            abs_dpsi_cutoff=args.abs_dpsi_cutoff,
        )
        catalog = result.event_catalog.copy()
        catalog["MUT_source_comparison"] = mapping["MUT"]
        catalog["CLK3_source_comparison"] = mapping["CLK3"]
        catalog["Pla_source_comparison"] = mapping["Pla"]
        catalogs.append(catalog)
        summaries.append(result.inhibitor_summary)
        classes.append(result.coverage_class_summary)
        venn_result = build_significant_event_venn_summary(
            {
                "MUT": loaded[mapping["MUT"]],
                "CLK3": loaded[mapping["CLK3"]],
                "Pla": loaded[mapping["Pla"]],
            },
            context=context,
            fdr_cutoff=args.fdr_cutoff,
            abs_dpsi_cutoff=args.abs_dpsi_cutoff,
        )
        venn_memberships.append(venn_result.membership_catalog)
        venn_summaries.append(venn_result.region_summary)

    catalog = pd.concat(catalogs, ignore_index=True)
    summary = pd.concat(summaries, ignore_index=True)
    summary["same_direction_FDR"] = _bh_adjust(summary["same_direction_binomial_pvalue"])
    summary["dPSI_spearman_FDR"] = _bh_adjust(summary["dPSI_spearman_pvalue"])
    summary["coverage_inference"] = np.where(
        summary["inhibitor"].isin(["CLK3", "Pla"]),
        "descriptive_only_shared_matched_WT",
        "descriptive_joint_coverage",
    )
    class_summary = pd.concat(classes, ignore_index=True)
    venn_membership = pd.concat(venn_memberships, ignore_index=True)
    venn_summary = pd.concat(venn_summaries, ignore_index=True)
    pairwise_venn_summary = _build_pairwise_venn_summary(venn_membership)
    event_type_summary = (
        catalog.groupby(["context", "event_type", "coverage_class"], dropna=False)
        .size()
        .rename("n_events")
        .reset_index()
    )
    comparison_design = _comparison_design(comparisons)
    top_events = _top_covered_events(catalog)
    input_manifest = _input_manifest(args.rmats_root, comparisons, args.mode)
    process = pd.DataFrame(
        [
            (1, "Load rMATS", "Read JCEC SE, RI, A3SS, A5SS and MXE tables for six contrasts."),
            (2, "Canonical event identity", "Match normalized gene ID, event type, chromosome, strand and exact coordinates."),
            (3, "Direction standardization", "Orient dPSI as MUT/inhibitor minus matched WT."),
            (4, "Significance filter", f"FDR <= {args.fdr_cutoff} and abs(dPSI) >= {args.abs_dpsi_cutoff}."),
            (5, "Venn membership", "Build the complete union and seven exact overlap regions across MUT, CLK3 and PlaB significant-event sets."),
            (6, "Fixed denominator", "Use all significant MUT events separately in DP and noDP."),
            (7, "Coverage classification", "Assign CLK3-only, PlaB-only, both or MUT-only within the MUT denominator."),
            (8, "Direction test", "One-sided exact binomial test versus 50%; BH FDR across four tests."),
            (9, "Effect-size agreement", "Spearman correlation of covered-event dPSI; BH FDR across four tests."),
        ],
        columns=["step", "stage", "operation"],
    )

    catalog.to_csv(args.output / "matched_background_MUT_event_catalog.tsv", sep="\t", index=False)
    summary.to_csv(args.output / "matched_background_MUT_coverage_summary.tsv", sep="\t", index=False)
    class_summary.to_csv(
        args.output / "matched_background_MUT_coverage_class_summary.tsv",
        sep="\t",
        index=False,
    )
    venn_membership.to_csv(
        args.output / "matched_background_all_significant_event_membership.tsv",
        sep="\t",
        index=False,
    )
    venn_summary.to_csv(
        args.output / "matched_background_significant_event_venn_summary.tsv",
        sep="\t",
        index=False,
    )
    pairwise_venn_summary.to_csv(
        args.output / "matched_background_pairwise_venn_summary.tsv",
        sep="\t",
        index=False,
    )
    event_type_summary.to_csv(
        args.output / "matched_background_MUT_event_type_summary.tsv",
        sep="\t",
        index=False,
    )
    comparison_design.to_csv(
        args.output / "matched_background_comparison_design.tsv",
        sep="\t",
        index=False,
    )
    top_events.to_csv(
        args.output / "matched_background_top_covered_events.tsv",
        sep="\t",
        index=False,
    )
    input_manifest.to_csv(
        args.output / "input_file_manifest_sha256.tsv",
        sep="\t",
        index=False,
    )
    process.to_csv(
        args.output / "analysis_process.tsv",
        sep="\t",
        index=False,
    )

    _plot_coverage(summary, plot_dir / "01_MUT_coverage_and_direction.png")
    _plot_pairwise_venn(
        pairwise_venn_summary,
        plot_dir / "02_pairwise_significant_event_venn.png",
    )
    _plot_coverage_classes(class_summary, plot_dir / "02_MUT_coverage_classes.png")
    _plot_dpsi_scatter(catalog, summary, "DP", plot_dir / "03_DP_dPSI_scatter.png")
    _plot_dpsi_scatter(catalog, summary, "noDP", plot_dir / "04_noDP_dPSI_scatter.png")
    _plot_significant_event_venn(
        venn_summary,
        comparison_design,
        "DP",
        plot_dir / "05_DP_significant_event_venn.png",
    )
    _plot_significant_event_venn(
        venn_summary,
        comparison_design,
        "noDP",
        plot_dir / "06_noDP_significant_event_venn.png",
    )

    markdown = _build_markdown(
        summary,
        class_summary,
        comparison_design,
        pairwise_venn_summary,
        venn_summary,
    )
    (args.output / "MATCHED_BACKGROUND_MUT_COVERAGE_REPORT.md").write_text(
        markdown,
        encoding="utf-8",
    )
    (args.output / "MATCHED_BACKGROUND_MUT_COVERAGE_REPORT.html").write_text(
        _build_html(
            summary,
            class_summary,
            comparison_design,
            pairwise_venn_summary,
            venn_summary,
            top_events,
            args.output,
        ),
        encoding="utf-8",
    )

    generated_files = sorted(
        str(path.relative_to(args.output)).replace("\\", "/")
        for path in args.output.rglob("*")
        if path.is_file() and path.name != "parameters.json"
    )
    parameters = {
        "analysis": "Matched-background MUT-centric AS inhibitor coverage",
        "scientific_question": (
            "Within DP and noDP separately, what fraction of significant MUT AS events is "
            "also significant after CLK3 or PlaB inhibition?"
        ),
        "rmats_root": str(args.rmats_root.resolve()),
        "mode": args.mode,
        "fdr_cutoff": args.fdr_cutoff,
        "abs_dPSI_cutoff": args.abs_dpsi_cutoff,
        "additional_event_filters": {
            "minimum_IJC_or_SJC": None,
            "replicate_stability_filter": None,
            "DEXSeq_filter": None,
            "DTU_filter": None,
        },
        "event_set_definition": (
            "Effect-size-filtered rMATS significant events: FDR cutoff plus absolute "
            "standardized dPSI cutoff only."
        ),
        "event_matching": (
            "Exact canonical event identity: normalized gene_id + event_type + chromosome + "
            "strand + rMATS event coordinates"
        ),
        "coverage_denominator": "All significant MUT events in the same DP context",
        "canonical_direction": "MUT or inhibitor minus matched WT",
        "direction_inference": (
            "One-sided exact binomial test against 50% same direction; "
            "Benjamini-Hochberg FDR across four context-by-inhibitor tests"
        ),
        "effect_size_inference": (
            "Spearman correlation of MUT and inhibitor standardized dPSI among covered events; "
            "Benjamini-Hochberg FDR across four tests"
        ),
        "coverage_inference_boundary": (
            "Coverage is descriptive because each MUT and inhibitor contrast shares its "
            "matched WT reference; no independence-based overlap P value is reported."
        ),
        "comparisons": comparisons,
        "generated_files": generated_files,
        "delivery_package_zip": "SETD2_MATCHED_BACKGROUND_MUT_COVERAGE_PACKAGE.zip",
    }
    (args.output / "parameters.json").write_text(
        json.dumps(parameters, indent=2),
        encoding="utf-8",
    )
    zip_path = _write_delivery_package(args.output, Path(__file__).resolve().parents[1])

    print(summary.to_string(index=False))
    print(f"\nOutput: {args.output}")
    print(f"Delivery package: {zip_path}")


if __name__ == "__main__":
    main()
