"""Fast, zero-BAM-rerun ingestion and publication Matplotlib rendering of 5'–3' gene-body coverage.

Supports direct ingestion of MultiQC outputs (`multiqc_data.json`,
`mqc_rseqc_gene_body_coverage_*.txt`, `mqc_qualimap_gene_coverage_profile_*.txt`,
`multiqc_qualimap_rnaseq.txt`) as well as existing RSeQC/Qualimap coverage tables,
rendering crisp native Matplotlib PNG/SVG exhibits immune to headless Vega font drops.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

MCK_NAVY = "#051C2C"
MCK_BLUE = "#2251FF"
MCK_CERULEAN = "#0077C8"
MCK_CYAN = "#00A9F4"
MCK_CORAL = "#B83A4B"
MCK_AMBER = "#D04A02"
MCK_SLATE = "#64748B"
MCK_PALETTE = [MCK_NAVY, MCK_BLUE, MCK_CORAL, MCK_CERULEAN, MCK_AMBER, MCK_CYAN, "#475569", "#0D9488"]


def parse_multiqc_5to3_coverage(qc_root: str | Path, sample_conditions: dict[str, str] | None = None) -> dict[str, Any] | None:
    """Discover and parse MultiQC or precomputed RSeQC/Qualimap 5'–3' coverage profiles in <50ms."""
    root = Path(qc_root).resolve()
    if not root.exists():
        return None
    sample_conditions = sample_conditions or {}
    curves_rows: list[dict[str, Any]] = []
    bias_summary: dict[str, float] = {}
    provenance_source = ""

    # 1. Prefer MultiQC multiqc_data.json if present
    mqc_jsons = sorted(root.rglob("multiqc_data.json"))
    for mqc_path in mqc_jsons:
        try:
            payload = json.loads(mqc_path.read_text(encoding="utf-8-sig"))
        except Exception:
            continue
        plot_data = payload.get("report_plot_data", {})
        for plot_key in (
            "rseqc_gene_body_coverage_plot",
            "rseqc_gene_body_coverage_plot_percentages",
            "qualimap_gene_coverage_profile",
            "qualimap_gene_coverage_profile_normalised",
        ):
            block = plot_data.get(plot_key)
            if not isinstance(block, dict):
                continue
            datasets = block.get("datasets", [])
            if datasets and isinstance(datasets[0], list):
                series_list = datasets[0]
            elif datasets and isinstance(datasets[0], dict) and "lines" in datasets[0]:
                series_list = datasets[0]["lines"]
            else:
                series_list = []
            for series in series_list:
                sid = str(series.get("name", "")).strip()
                pts = series.get("data", [])
                if not sid or not pts:
                    continue
                for pt in pts:
                    if isinstance(pt, (list, tuple)) and len(pt) >= 2:
                        curves_rows.append({"sample_id": sid, "percentile": float(pt[0]), "coverage": float(pt[1])})
                    elif isinstance(pt, dict) and "x" in pt and "y" in pt:
                        curves_rows.append({"sample_id": sid, "percentile": float(pt["x"]), "coverage": float(pt["y"])})
            if curves_rows:
                provenance_source = f"MultiQC ({mqc_path.name} :: {plot_key})"
                break
        # Extract native Qualimap 5'-3' bias from general_stats if recorded
        gen_stats = payload.get("report_general_stats_data", [])
        blocks = gen_stats if isinstance(gen_stats, list) else list(gen_stats.values()) if isinstance(gen_stats, dict) else []
        for gblock in blocks:
            if not isinstance(gblock, dict):
                continue
            for sid, metrics in gblock.items():
                if isinstance(metrics, dict):
                    for k in ("5_3_bias", "qualimap_5_3_bias", "5'-3' bias"):
                        if k in metrics and metrics[k] is not None:
                            try:
                                bias_summary[str(sid)] = float(metrics[k])
                            except (TypeError, ValueError):
                                pass
        if curves_rows:
            break

    # 2. Check MultiQC flat export text tables (`mqc_rseqc_gene_body_coverage*.txt` / `mqc_qualimap_gene_coverage*.txt`)
    if not curves_rows:
        flat_candidates = sorted(
            p for p in root.rglob("*.txt")
            if any(tag in p.name.lower() for tag in ("mqc_rseqc_gene_body_coverage", "mqc_qualimap_gene_coverage", "genebodycoverage"))
        ) + sorted(
            p for p in root.rglob("*.tsv")
            if "gene_body_coverage" in p.name.lower() or "coverage_bias" in p.name.lower()
        )
        for cand in flat_candidates:
            try:
                raw = pd.read_csv(cand, sep="\t")
            except Exception:
                continue
            cols = [str(c) for c in raw.columns]
            if {"sample_id", "percentile", "coverage"} <= set(cols):
                sub = raw[["sample_id", "percentile", "coverage"]].copy()
                for _, r in sub.iterrows():
                    curves_rows.append({"sample_id": str(r["sample_id"]), "percentile": float(r["percentile"]), "coverage": float(r["coverage"])})
                provenance_source = f"Precomputed coverage table ({cand.name})"
                break
            elif len(cols) >= 11:
                # Wide matrix: row = Sample, columns = percentiles 1..100
                num_cols = []
                for c in cols[1:]:
                    try:
                        num_cols.append((c, float(c)))
                    except ValueError:
                        pass
                if len(num_cols) >= 10:
                    for _, r in raw.iterrows():
                        sid = str(r.iloc[0]).strip()
                        for col_name, pct in num_cols:
                            val = pd.to_numeric(r[col_name], errors="coerce")
                            if pd.notna(val):
                                curves_rows.append({"sample_id": sid, "percentile": pct, "coverage": float(val)})
                    if curves_rows:
                        provenance_source = f"MultiQC / RSeQC flat table ({cand.name})"
                        break

    if not curves_rows:
        return None

    df = pd.DataFrame(curves_rows)
    df = df.dropna(subset=["sample_id", "percentile", "coverage"]).copy()
    df = df[(df["percentile"] >= 0) & (df["percentile"] <= 100) & np.isfinite(df["coverage"])]
    if df.empty:
        return None
    df = df.groupby(["sample_id", "percentile"], as_index=False)["coverage"].mean()
    means = df.groupby("sample_id")["coverage"].transform("mean").replace(0, np.nan)
    df["relative_coverage"] = df["coverage"] / means
    df["condition"] = df["sample_id"].map(lambda s: sample_conditions.get(str(s), "Cohort"))

    summary_rows = []
    for sid, grp in df.groupby("sample_id", sort=True):
        grp = grp.sort_values("percentile")
        five = float(grp.loc[grp["percentile"] <= 20, "relative_coverage"].mean())
        mid = float(grp.loc[(grp["percentile"] > 20) & (grp["percentile"] < 80), "relative_coverage"].mean())
        three = float(grp.loc[grp["percentile"] >= 80, "relative_coverage"].mean())
        ratio_3_5 = three / five if five > 0 else np.nan
        summary_rows.append({
            "sample_id": str(sid),
            "condition": str(grp["condition"].iloc[0]),
            "five_prime_rel_mean": round(five, 4),
            "mid_body_rel_mean": round(mid, 4),
            "three_prime_rel_mean": round(three, 4),
            "three_to_five_window_ratio": round(float(ratio_3_5), 4) if np.isfinite(ratio_3_5) else np.nan,
            "multiqc_native_5_3_bias": bias_summary.get(str(sid), np.nan),
            "provenance_source": provenance_source,
        })
    return {
        "curves": df.sort_values(["condition", "sample_id", "percentile"]).reset_index(drop=True),
        "summary": pd.DataFrame(summary_rows),
        "provenance_source": provenance_source,
    }


def render_mckinsey_5to3_figure(
    parsed: dict[str, Any],
    png_path: str | Path,
    svg_path: str | Path | None = None,
    title: str = "Transcript 5′ → 3′ Gene-Body Coverage & Terminal Balance",
) -> None:
    """Render a crisp 2-panel executive consulting native Matplotlib figure (guaranteed zero missing glyphs)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    curves: pd.DataFrame = parsed["curves"]
    summary: pd.DataFrame = parsed["summary"]
    conditions = list(dict.fromkeys(curves["condition"].astype(str)))
    cond_colors = {c: MCK_PALETTE[i % len(MCK_PALETTE)] for i, c in enumerate(conditions)}

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    fig, (ax_curve, ax_ratio) = plt.subplots(
        1, 2, figsize=(12.8, 4.6), gridspec_kw={"width_ratios": [2.25, 1.0]}, dpi=220
    )
    fig.patch.set_facecolor("white")

    # Left panel: 5' -> 3' curves
    ax_curve.axvspan(0, 20, color="#F1F5F9", alpha=0.85, zorder=0)
    ax_curve.axvspan(80, 100, color="#F1F5F9", alpha=0.85, zorder=0)
    ax_curve.axhline(1.0, color="#CBD5E1", linestyle="--", linewidth=0.9, zorder=1)
    ax_curve.grid(axis="y", color="#E2E8F0", linewidth=0.6, zorder=0)

    for _, sgrp in curves.groupby("sample_id", sort=False):
        cond = str(sgrp["condition"].iloc[0])
        ax_curve.plot(
            sgrp["percentile"],
            sgrp["relative_coverage"],
            color=cond_colors[cond],
            alpha=0.28 if len(summary) > 4 else 0.45,
            linewidth=1.0,
            zorder=2,
        )

    for cond in conditions:
        cgrp = curves.loc[curves["condition"] == cond].groupby("percentile", as_index=False)["relative_coverage"].mean()
        n_samp = int((summary["condition"] == cond).sum())
        ax_curve.plot(
            cgrp["percentile"],
            cgrp["relative_coverage"],
            color=cond_colors[cond],
            linewidth=2.4,
            label=f"{cond} (n={n_samp})",
            zorder=4,
        )

    ax_curve.text(10, ax_curve.get_ylim()[1] * 0.96, "5′ window\n(0–20%)", ha="center", va="top", fontsize=7.8, color=MCK_SLATE)
    ax_curve.text(90, ax_curve.get_ylim()[1] * 0.96, "3′ window\n(80–100%)", ha="center", va="top", fontsize=7.8, color=MCK_SLATE)
    ax_curve.set_xlim(0, 100)
    ax_curve.set_xlabel("Transcript body percentile (5′ → 3′)", fontsize=9.5, color=MCK_NAVY, labelpad=6)
    ax_curve.set_ylabel("Normalized coverage (× sample mean)", fontsize=9.5, color=MCK_NAVY)
    ax_curve.set_title("A. Normalized 5′ → 3′ coverage profile", loc="left", fontsize=10.5, fontweight="bold", color=MCK_NAVY, pad=8)
    ax_curve.legend(frameon=False, fontsize=8.5, loc="lower center", ncol=min(3, max(1, len(conditions))))

    # Right panel: 3'/5' balance strip
    ax_ratio.axvline(1.0, color="#94A3B8", linestyle="--", linewidth=1.0, zorder=1)
    ax_ratio.grid(axis="x", color="#E2E8F0", linewidth=0.6, zorder=0)
    y_positions = np.arange(len(conditions))
    rng = np.random.default_rng(42)

    for yi, cond in enumerate(conditions):
        vals = summary.loc[summary["condition"] == cond, "three_to_five_window_ratio"].dropna().to_numpy(dtype=float)
        if not len(vals):
            continue
        jitter = rng.uniform(-0.14, 0.14, size=len(vals)) if len(vals) > 1 else np.zeros(1)
        ax_ratio.scatter(
            vals,
            np.full(len(vals), yi) + jitter,
            s=48,
            color=cond_colors[cond],
            edgecolors="white",
            linewidths=0.7,
            alpha=0.9,
            zorder=3,
        )
        med = float(np.median(vals))
        ax_ratio.plot([med, med], [yi - 0.26, yi + 0.26], color=MCK_NAVY, linewidth=2.2, zorder=4)
        ax_ratio.text(med, yi + 0.31, f"{med:.2f}×", ha="center", va="bottom", fontsize=8, color=MCK_NAVY, fontweight="bold")

    ax_ratio.set_yticks(y_positions, conditions, fontsize=8.8, color=MCK_NAVY)
    ax_ratio.set_ylim(-0.6, max(0.6, len(conditions) - 0.4))
    ax_ratio.invert_yaxis()
    ax_ratio.set_xlabel("3′ / 5′ terminal coverage ratio (80–100% ÷ 0–20%)", fontsize=9.2, color=MCK_NAVY, labelpad=6)
    ax_ratio.set_title("B. 3′ / 5′ terminal balance", loc="left", fontsize=10.5, fontweight="bold", color=MCK_NAVY, pad=8)

    fig.suptitle(title, x=0.01, ha="left", fontsize=12.5, fontweight="bold", color=MCK_NAVY, y=1.02)
    fig.tight_layout()
    png_path = Path(png_path)
    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=220, bbox_inches="tight", facecolor="white")
    if svg_path:
        fig.savefig(svg_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
