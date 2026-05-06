from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


class CardsAdapter:
    GROUP_COLORS = {
        "WT_NO_DP": "#4C78A8",
        "WT_DP": "#F58518",
        "KO_NO_DP": "#54A24B",
        "KO_DP": "#E45756",
    }

    def render_cards(self, shortlist: pd.DataFrame, expression_support: pd.DataFrame, output_root: Path) -> list[Path]:
        output_root.mkdir(parents=True, exist_ok=True)
        generated: list[Path] = []
        for _, row in shortlist.iterrows():
            display_group = str(row.get("display_group", "shortlist"))
            outdir = output_root / self._sanitize(display_group)
            outdir.mkdir(parents=True, exist_ok=True)
            fig = plt.figure(figsize=(14, 4.8))
            ax1 = plt.subplot(1, 3, 1)
            ax2 = plt.subplot(1, 3, 2)
            ax3 = plt.subplot(1, 3, 3)

            self._plot_abs_dpsi(ax1, row)
            self._plot_signed_dpsi(ax2, row)
            self._plot_expression(ax3, expression_support, str(row.get("GeneID", "")))

            title = f"{row.get('geneSymbol', '')} | {row.get('event_type', '')} | {display_group}"
            subtitle = (
                f"class={row.get('class_interpretation', '')} | "
                f"A={row.get('comparison_A', '')} | "
                f"B={row.get('comparison_B', '')} | "
                f"rank={row.get('rank_score', '')}"
            )
            fig.suptitle(title + "\n" + subtitle, fontsize=11)
            plt.tight_layout(rect=[0, 0, 1, 0.88])

            out_png = outdir / f"{self._sanitize(str(row.get('geneSymbol', 'gene')))}__{self._sanitize(str(row.get('event_type', 'event')))}__card.png"
            plt.savefig(out_png, dpi=300)
            plt.close()
            generated.append(out_png)
        return generated

    def _plot_abs_dpsi(self, ax, row: pd.Series) -> None:
        labels = [str(row.get("comparison_A", "")), str(row.get("comparison_B", ""))]
        vals = [
            float(row.get("abs_dPSI_A", 0) or 0),
            float(row.get("abs_dPSI_B", 0) or 0),
        ]
        ax.bar([0, 1], vals, color=["#4C78A8", "#E45756"])
        ax.set_xticks([0, 1])
        ax.set_xticklabels(labels, rotation=20, ha="right")
        ax.set_ylabel("abs(dPSI)")
        ax.set_title("Comparison effect size")

    def _plot_signed_dpsi(self, ax, row: pd.Series) -> None:
        labels = [str(row.get("comparison_A", "")), str(row.get("comparison_B", ""))]
        vals = [
            float(row.get("dPSI_A", 0) or 0),
            float(row.get("dPSI_B", 0) or 0),
        ]
        ax.axhline(0, linewidth=1, color="black")
        ax.bar([0, 1], vals, color=["#4C78A8", "#E45756"])
        ax.set_xticks([0, 1])
        ax.set_xticklabels(labels, rotation=20, ha="right")
        ax.set_ylabel("dPSI")
        ax.set_title(str(row.get("class_interpretation", "")))

    def _plot_expression(self, ax, expr_df: pd.DataFrame, gene_id: str) -> None:
        sub = expr_df[expr_df["gene_id"].astype(str) == gene_id].copy() if not expr_df.empty and "gene_id" in expr_df.columns else pd.DataFrame()
        if sub.empty:
            ax.text(0.5, 0.5, f"No expression support for\n{gene_id}", ha="center", va="center")
            ax.set_axis_off()
            return
        groups = sub["group"].astype(str).tolist()
        exprs = [float(value) for value in sub["expr"].tolist()]
        colors = [self.GROUP_COLORS.get(group, "#888888") for group in groups]
        ax.bar(range(len(groups)), exprs, color=colors)
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels(groups, rotation=20, ha="right")
        ax.set_ylabel("Mean normalized expression")
        ax.set_title("Expression support")

    def _sanitize(self, value: str) -> str:
        value = re.sub(r"[^A-Za-z0-9._-]+", "_", value)
        return re.sub(r"_+", "_", value).strip("_")
