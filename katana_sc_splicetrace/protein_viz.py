"""Publication-grade visualization suite (5 distinct styles) for Splicing -> ORF / NMD / UniProt Domain Consequences."""
from __future__ import annotations

import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.path import Path as MplPath
from matplotlib.patches import FancyBboxPatch, PathPatch, Polygon, Rectangle
import numpy as np
import pandas as pd

NAVY = "#1e2630"
BLUE = "#2863a5"
RED = "#c63f4a"
TEAL = "#2d7e9b"
GREEN = "#238b45"
GOLD = "#d49a25"
PURPLE = "#6b4c9a"
GRAY = "#5b6470"
LIGHT = "#e3e7eb"

DOMAIN_PALETTE = {
    "DOMAIN": "#2863a5",
    "REPEAT": "#2d7e9b",
    "ZN_FING": "#6b4c9a",
    "DNA_BIND": "#238b45",
    "TRANSMEM": "#d49a25",
    "MOTIF": "#c63f4a",
    "COILED": "#7a8694",
    "TOPO_DOM": "#8bb9d7",
    "REGION": "#9aa5b1",
}


def _save_dual(fig: Figure, base: Path, dpi: int = 180) -> tuple[Path, Path]:
    base.parent.mkdir(parents=True, exist_ok=True)
    FigureCanvasAgg(fig)
    png_p = base.with_suffix(".png")
    svg_p = base.with_suffix(".svg")
    fig.savefig(png_p, format="png", dpi=dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(svg_p, format="svg", bbox_inches="tight", facecolor="white")
    fig.clear()
    return png_p, svg_p


def classify_proteome_fate(row: pd.Series | dict) -> str:
    nmd = str(row["nmd_switch_class"])
    fst = str(row["frame_status"])
    dom = str(row["disrupted_protein_features"])
    dpsi = float(row["dpsi"])
    if nmd == "inclusion_triggers_NMD_poison_exon":
        return "NMD Decay (Poison Exon Included)" if dpsi > 0 else "NMD Rescue (Poison Exon Skipped)"
    if nmd == "skipping_triggers_NMD_decay":
        return "NMD Decay (Essential Exon Skipped)" if dpsi < 0 else "ORF Preserved (Essential Exon Included)"
    if fst in ("utr_cis_regulatory_only", "noncoding_transcript") or nmd in ("coding_silent_utr_switch", "noncoding"):
        return "5'/3' UTR or Non-Coding Switch"
    if dom != "none":
        return "Productive UniProt/Pfam Domain Altered"
    return "Productive Linker / Tail Isoform Switch"


# ---------------------------------------------------------------------------
# STYLE 1: Global Proteome Consequence Breakdown + ΔΨ × ΔORF Quadrant Map
# ---------------------------------------------------------------------------
def plot_consequence_landscape_v24(
    contrast_dfs: dict[str, pd.DataFrame],
    focal_contrast: str,
    out_base: Path,
    dataset_title: str = "Study",
) -> tuple[Path, Path]:
    cats = [
        ("Productive UniProt/Pfam Domain Altered", BLUE),
        ("Productive Linker / Tail Isoform Switch", "#8bb9d7"),
        ("NMD Decay (Poison Exon Included)", RED),
        ("NMD Decay (Essential Exon Skipped)", "#e06d53"),
        ("NMD Rescue (Poison Exon Skipped)", GREEN),
        ("ORF Preserved (Essential Exon Included)", TEAL),
        ("5'/3' UTR or Non-Coding Switch", "#b8c0c8"),
    ]
    fig = Figure(figsize=(13.6, 4.9), layout="constrained")
    gs = fig.add_gridspec(1, 2, width_ratios=[1.18, 1.15])
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])

    cids = list(contrast_dfs.keys())
    y_pos = np.arange(len(cids))[::-1]
    lefts = np.zeros(len(cids), dtype=float)
    for cat_label, col in cats:
        vals = []
        for cid in cids:
            df_c = contrast_dfs[cid]
            if df_c.empty:
                vals.append(0.0)
            else:
                fates = df_c.apply(classify_proteome_fate, axis=1)
                vals.append(100.0 * float((fates == cat_label).mean()))
        vals_arr = np.asarray(vals, dtype=float)
        if np.any(vals_arr > 0):
            ax0.barh(y_pos, vals_arr, left=lefts, height=0.58, color=col, edgecolor="white", linewidth=0.7, label=cat_label)
            for yi, (lv, vv) in enumerate(zip(lefts, vals_arr)):
                if vv >= 9.0:
                    ax0.text(
                        lv + 0.5 * vv,
                        y_pos[yi],
                        f"{vv:.0f}%",
                        ha="center",
                        va="center",
                        fontsize=7.8,
                        color="white" if col in (BLUE, RED, GREEN, TEAL, "#e06d53") else NAVY,
                        fontweight="bold",
                    )
            lefts += vals_arr

    ax0.set_yticks(y_pos)
    ax0.set_yticklabels(cids, fontsize=8.8, fontweight="bold", color=NAVY)
    ax0.set_xlim(0, 100)
    ax0.set_xlabel("% of Regulated Alternative Splicing Events", fontsize=9, color=NAVY)
    ax0.spines[["top", "right"]].set_visible(False)
    ax0.set_title(f"{dataset_title} | Proteome & NMD Composition", fontsize=9.8, fontweight="bold", color=NAVY)
    ax0.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, frameon=False, fontsize=7.6)

    # Right panel: ΔΨ vs Signed % ORF Length Change on focal_contrast
    df_f = contrast_dfs[focal_contrast].copy()
    df_f = df_f[df_f["aa_len_inclusion"].gt(0) | df_f["aa_len_exclusion"].gt(0)].copy()
    max_aa = np.maximum(np.maximum(df_f["aa_len_inclusion"], df_f["aa_len_exclusion"]), 1)
    # Signed effective ΔORF % in numerator condition: sign(dpsi) * (len_inc - len_exc)/max_aa * 100
    eff_pct = np.sign(df_f["dpsi"]) * (df_f["aa_len_inclusion"] - df_f["aa_len_exclusion"]) / max_aa * 100.0
    df_f["eff_pct"] = eff_pct
    fates_f = df_f.apply(classify_proteome_fate, axis=1)
    col_map = dict(cats)
    pt_cols = [col_map.get(ft, GRAY) for ft in fates_f]

    ax1.axhline(0, color=LIGHT, linewidth=1.0, zorder=1)
    ax1.axvline(0, color=LIGHT, linewidth=1.0, zorder=1)
    ax1.scatter(df_f["dpsi"], df_f["eff_pct"], c=pt_cols, s=36, alpha=0.82, edgecolors=NAVY, linewidths=0.4, zorder=3)

    # Label top dramatic NMD or domain-altered targets with y-repulsion
    df_f["rank_score"] = df_f["dpsi"].abs() * (1.0 + np.abs(df_f["eff_pct"]) / 35.0) + np.where(df_f["disrupted_protein_features"].ne("none"), 0.35, 0.0)
    top_lbl = df_f.sort_values("rank_score", ascending=False).drop_duplicates("gene_symbol").head(8)
    placed_side: dict[str, list[float]] = {"L": [], "R": []}
    for r in top_lbl.itertuples(index=False):
        side = "R" if r.dpsi >= 0 else "L"
        y_cand = float(r.eff_pct)
        for prev_y in placed_side[side]:
            if abs(y_cand - prev_y) < 12.0:
                y_cand = prev_y - 12.5 if y_cand <= prev_y else prev_y + 12.5
        placed_side[side].append(y_cand)
        ax1.annotate(
            f"{r.gene_symbol} ({r.event_type})",
            (r.dpsi, r.eff_pct),
            xytext=(r.dpsi + (0.03 if side == "R" else -0.03), np.clip(y_cand, -98, 98)),
            ha="left" if side == "R" else "right",
            va="center",
            fontsize=7.6,
            fontweight="bold",
            color=NAVY,
            arrowprops={"arrowstyle": "-", "color": GRAY, "lw": 0.6},
        )
    max_x = max(0.56, float(df_f["dpsi"].abs().max()) + 0.18)
    ax1.set_xlim(-max_x, max_x)
    ax1.set_ylim(-108, 108)
    ax1.spines[["top", "right"]].set_visible(False)
    ax1.set_xlabel(f"Splicing Effect Size ΔΨ ({focal_contrast})", fontsize=9, color=NAVY)
    ax1.set_ylabel("Condition-Favored ΔORF Length (% of Canonical Protein)", fontsize=9, color=NAVY)
    ax1.set_title(f"ΔΨ × ORF Truncation / Indel Landscape | {focal_contrast}", fontsize=11, fontweight="bold", color=NAVY)
    return _save_dual(fig, out_base)


# ---------------------------------------------------------------------------
# STYLE 2: Dual-Track mRNA (Last EJC + NMD Zone) & UniProt Domain + PTC Card
# ---------------------------------------------------------------------------
def plot_dual_track_nmd_architecture_v24(
    prot_df: pd.DataFrame,
    motif_map: dict[str, str],
    out_base: Path,
    contrast_id: str,
    max_genes: int = 3,
) -> tuple[Path, Path]:
    # Pick top NMD-switching genes that also have annotated UniProt domains!
    nmd_sub = prot_df[
        prot_df["nmd_switch_class"].isin(["inclusion_triggers_NMD_poison_exon", "skipping_triggers_NMD_decay"])
        & prot_df["disrupted_protein_features"].ne("none")
    ].copy()
    if len(nmd_sub) < max_genes:
        extra = prot_df[prot_df["nmd_switch_class"].isin(["inclusion_triggers_NMD_poison_exon", "skipping_triggers_NMD_decay"])]
        nmd_sub = pd.concat([nmd_sub, extra], ignore_index=True)
    nmd_sub["_abs"] = nmd_sub["dpsi"].abs()
    picks = nmd_sub.sort_values("_abs", ascending=False).drop_duplicates("gene_symbol").head(max_genes)

    fig = Figure(figsize=(12.8, max(3.6, 2.05 * len(picks) + 0.6)), layout="constrained")
    axes = fig.subplots(max(len(picks), 1), 1, squeeze=False)[:, 0]
    if picks.empty:
        axes[0].text(0.5, 0.5, "No NMD switch events in selection", ha="center", va="center")
        axes[0].axis("off")
        return _save_dual(fig, out_base)

    for ax, row in zip(axes, picks.itertuples(index=False)):
        arch = json.loads(str(row.architecture_json)) if getattr(row, "architecture_json", "{}") else {}
        is_poison = row.nmd_switch_class == "inclusion_triggers_NMD_poison_exon"
        mrna_len = max(int(arch.get("mrna_nt_inc" if is_poison else "mrna_nt_exc", 2000) or 2000), 200)
        last_ejc = arch.get("last_ejc_inc" if is_poison else "last_ejc_exc")
        stop_nt = arch.get("stop_nt_inc" if is_poison else "stop_nt_exc")
        full_aa = max(int(row.aa_len_inclusion), int(row.aa_len_exclusion), 1)
        trunc_aa = min(int(row.aa_len_inclusion), int(row.aa_len_exclusion))

        ax.set_xlim(-35, 117)
        ax.set_ylim(-0.25, 2.38)
        ax.axis("off")

        # Header label cleanly in left margin column [-34.5 .. -10]
        m_tag = motif_map.get(str(row.event_id), "")
        if len(m_tag) > 30:
            m_tag = m_tag[:28] + ".."
        m_line = f"\nMotif: {m_tag}" if m_tag else ""
        ax.text(
            -34.5,
            1.15,
            f"{row.gene_symbol} ({row.event_type}) · {row.uniprot_acc}\nΔΨ = {row.dpsi:+.2f} ({row.frame_status}){m_line}",
            va="center",
            ha="left",
            fontsize=8.3,
            fontweight="bold",
            color=NAVY,
        )

        # --- UPPER TRACK (y = 1.68): Mature mRNA + Last EJC + 50-nt NMD Trigger Zone ---
        ax.text(-1.2, 1.68, "mRNA", ha="right", va="center", fontsize=8, color=GRAY, fontweight="bold")
        ax.plot([0, 100], [1.68, 1.68], color=NAVY, linewidth=2.2, zorder=1)
        ax.add_patch(FancyBboxPatch((0, 1.54), 100, 0.28, boxstyle="round,pad=0.01", facecolor="#edf2f7", edgecolor=NAVY, linewidth=0.8, zorder=2))

        if last_ejc and stop_nt and mrna_len > 0:
            x_ejc = 100.0 * float(last_ejc) / mrna_len
            x_nmd_lim = max(0.0, 100.0 * float(last_ejc - 50) / mrna_len)
            x_ptc_nt = 100.0 * float(stop_nt) / mrna_len
            ax.add_patch(Rectangle((0, 1.54), x_nmd_lim, 0.28, facecolor="#fdecee", edgecolor="none", zorder=2.2))
            ax.plot([x_ejc, x_ejc], [1.44, 1.92], color=BLUE, linestyle="--", linewidth=1.4, zorder=4)
            ax.text(x_ejc, 1.96, f"Last EJC ({last_ejc} nt)", ha="center", va="bottom", fontsize=7.5, color=BLUE, fontweight="bold")
            ax.scatter([x_ptc_nt], [1.68], color=RED, s=75, marker="X", zorder=5)
            ax.text(x_ptc_nt, 1.40, f"PTC ({row.ptc_to_last_ejc_nt} nt upstream ≥50 nt → NMD!)", ha="center", va="top", fontsize=7.5, color=RED, fontweight="bold")

        # --- LOWER TRACK (y = 0.55): Full-Length Protein Domains + Truncated Ghost Region ---
        ax.text(-1.2, 0.58, f"{full_aa} aa", ha="right", va="center", fontsize=8, color=GRAY, fontweight="bold")
        x_trunc = 100.0 * float(trunc_aa) / full_aa
        ax.add_patch(Rectangle((0, 0.44), x_trunc, 0.28, facecolor="#dce5ef", edgecolor=NAVY, linewidth=0.9, zorder=2))
        if x_trunc < 99.5:
            ax.add_patch(Rectangle((x_trunc, 0.44), 100.0 - x_trunc, 0.28, facecolor="#f8f9fa", edgecolor=RED, linestyle="--", linewidth=1.0, hatch="///", alpha=0.65, zorder=2))

        raw_doms = arch.get("domains", [])
        struct_doms = [d for d in raw_doms if d["kind"] != "REGION"]
        doms = struct_doms if struct_doms else [d for d in raw_doms if "Disordered" not in d["label"]]
        last_lbl_x = [-99.0, -99.0]
        for d in doms[:8]:
            x0 = 100.0 * (d["start"] - 1) / full_aa
            w_d = max(1.6, 100.0 * (d["end"] - d["start"] + 1) / full_aa)
            d_col = DOMAIN_PALETTE.get(d["kind"], BLUE)
            lost = d["state"] == "lost_by_PTC_truncation" or d["start"] > trunc_aa
            ax.add_patch(
                FancyBboxPatch(
                    (x0, 0.40),
                    w_d,
                    0.36,
                    boxstyle="round,pad=0.01",
                    facecolor="white" if lost else d_col,
                    edgecolor=RED if lost else NAVY,
                    linestyle="--" if lost else "-",
                    linewidth=1.0,
                    alpha=0.72 if lost else 0.92,
                    zorder=3,
                )
            )
            xm = x0 + 0.5 * w_d
            short_lbl = d["label"].split(":", 1)[-1][:18]
            for slot_i, y_slot in enumerate([0.12, -0.10]):
                if abs(xm - last_lbl_x[slot_i]) >= 15.0:
                    ax.text(xm, y_slot, f"{short_lbl}{'*LOST' if lost else ''}", ha="center", va="center", fontsize=7.1, color=RED if lost else NAVY, fontweight="bold")
                    last_lbl_x[slot_i] = xm
                    break

        ax.plot([x_trunc, x_trunc], [0.32, 0.86], color=RED, linewidth=2.0, zorder=5)
        ax.text(x_trunc, 0.92, f"PTC Stop ({trunc_aa}* aa)", ha="center", va="bottom", fontsize=7.8, color=RED, fontweight="bold")
        mode_txt = "Inclusion → Poison NMD" if is_poison else "Skipping → NMD Decay"
        ax.text(101.2, 0.58, f"{mode_txt}\n({trunc_aa}* vs {full_aa} aa)", va="center", ha="left", fontsize=7.6, color=RED, fontweight="bold")

    fig.suptitle(f"Dual-Track Mature mRNA (50-nt EJC NMD Rule) & UniProt Protein Truncation Architecture | {contrast_id}", fontsize=11.2, fontweight="bold", color=NAVY)
    return _save_dual(fig, out_base)


# ---------------------------------------------------------------------------
# STYLE 3: Paired Inclusion vs Exclusion Isoform In-Frame Domain Surgery Plot
# ---------------------------------------------------------------------------
def plot_side_by_side_inframe_domains_v24(
    prot_df: pd.DataFrame,
    out_base: Path,
    contrast_id: str,
    max_genes: int = 4,
) -> tuple[Path, Path]:
    sub = prot_df[
        prot_df["frame_status"].eq("in_frame_indel")
        & prot_df["nmd_switch_class"].eq("productive_protein_isoform_switch")
        & prot_df["disrupted_protein_features"].ne("none")
        & (prot_df["delta_aa"].abs() >= 4)
    ].copy()
    sub["_abs"] = sub["dpsi"].abs()
    picks = sub.sort_values("_abs", ascending=False).drop_duplicates("gene_symbol").head(max_genes)

    fig = Figure(figsize=(12.6, max(3.4, 1.65 * max(len(picks), 1) + 0.6)), layout="constrained")
    axes = fig.subplots(max(len(picks), 1), 1, squeeze=False)[:, 0]
    if picks.empty:
        axes[0].text(0.5, 0.5, "No in-frame domain-remodeling events in selection", ha="center", va="center")
        axes[0].axis("off")
        return _save_dual(fig, out_base)

    for ax, row in zip(axes, picks.itertuples(index=False)):
        arch = json.loads(str(row.architecture_json)) if getattr(row, "architecture_json", "{}") else {}
        l_inc, l_exc = int(row.aa_len_inclusion), int(row.aa_len_exclusion)
        max_l = max(l_inc, l_exc, 1)
        junc_aa = int(arch.get("splice_junction_aa", 1))
        span_aa = abs(l_inc - l_exc)

        ax.set_xlim(-20, 112)
        ax.set_ylim(-0.1, 1.85)
        ax.axis("off")

        ax.text(
            -19.5,
            0.90,
            f"{row.gene_symbol} ({row.event_type})\nUniProt: {row.uniprot_acc}\nΔΨ = {row.dpsi:+.2f} ({row.delta_aa:+d} aa)",
            va="center",
            ha="left",
            fontsize=8.6,
            fontweight="bold",
            color=NAVY,
        )

        for y_bar, iso_lbl, iso_len, is_longer in [
            (1.22, f"Inclusion ({l_inc} aa)", l_inc, l_inc >= l_exc),
            (0.42, f"Exclusion ({l_exc} aa)", l_exc, l_exc > l_inc),
        ]:
            w_bar = 100.0 * iso_len / max_l
            ax.text(-1.0, y_bar + 0.12, iso_lbl, ha="right", va="center", fontsize=7.8, color=NAVY)
            ax.add_patch(Rectangle((0, y_bar), w_bar, 0.24, facecolor="#e5ebf1", edgecolor=NAVY, linewidth=0.8, zorder=2))
            if is_longer and span_aa > 0:
                xj0 = 100.0 * (junc_aa - 1) / max_l
                wj = max(1.8, 100.0 * span_aa / max_l)
                ax.add_patch(Rectangle((xj0, y_bar - 0.04), wj, 0.32, facecolor=RED if row.dpsi > 0 else BLUE, edgecolor=NAVY, linewidth=0.9, zorder=4))

        # Overlay UniProt structural domains on reference coordinates without label overlap
        raw_d = arch.get("domains", [])
        struct_d = [d for d in raw_d if d["kind"] != "REGION"]
        doms = struct_d if struct_d else [d for d in raw_d if "Disordered" not in d["label"]]
        prev_xm = -99.0
        for d in sorted(doms[:7], key=lambda item: item["start"]):
            x0 = 100.0 * (d["start"] - 1) / max_l
            wd = max(1.8, 100.0 * (d["end"] - d["start"] + 1) / max_l)
            col = DOMAIN_PALETTE.get(d["kind"], TEAL)
            ax.add_patch(FancyBboxPatch((x0, 1.20), wd, 0.28, boxstyle="round,pad=0.01", facecolor=col, edgecolor=NAVY, linewidth=0.7, alpha=0.82, zorder=3))
            xm = x0 + 0.5 * wd
            if xm - prev_xm >= 17.5:
                lbl = d["label"].split(":", 1)[-1][:17]
                ax.text(xm, 1.58, lbl, ha="center", va="bottom", fontsize=7.2, color=RED if d["state"] != "intact" else NAVY, fontweight="bold")
                prev_xm = xm

        ax.text(101.2, 0.85, str(row.disrupted_protein_features).replace("; ", "\n")[:68], va="center", ha="left", fontsize=7.5, color=NAVY)

    fig.suptitle(f"Paired Isoform In-Frame UniProt Domain Remodeling (Highlighted Segment = Spliced Cassette Peptide) | {contrast_id}", fontsize=11.0, fontweight="bold", color=NAVY)
    return _save_dual(fig, out_base)


# ---------------------------------------------------------------------------
# STYLE 4: Single-Protein Multi-Event Lollipop Splicing Map (cBioPortal Style)
# ---------------------------------------------------------------------------
def plot_lollipop_protein_splicing_map_v24(
    prot_df: pd.DataFrame,
    out_base: Path,
    contrast_id: str,
    target_gene: str | None = None,
) -> tuple[Path, Path]:
    if target_gene is None:
        vc = prot_df[prot_df["disrupted_protein_features"].ne("none")]["gene_symbol"].value_counts()
        if vc.empty:
            vc = prot_df["gene_symbol"].value_counts()
        target_gene = str(vc.index[0])

    g_rows = prot_df[prot_df["gene_symbol"].eq(target_gene)].copy()
    # Pick row with longest architecture domain list for backbone display
    best_idx = g_rows["aa_len_inclusion"].idxmax()
    first = g_rows.loc[best_idx]
    arch = json.loads(str(first["architecture_json"]))
    prot_len = max(int(g_rows["aa_len_inclusion"].max()), int(g_rows["aa_len_exclusion"].max()), 1)

    fig = Figure(figsize=(12.4, 3.7), layout="constrained")
    ax = fig.subplots()
    ax.axhline(0.0, color=NAVY, linewidth=1.2, zorder=1)
    ax.add_patch(FancyBboxPatch((1, -0.07), prot_len, 0.14, boxstyle="round,pad=0.01", facecolor="#dde5ed", edgecolor=NAVY, linewidth=1.0, zorder=2))

    doms = [d for d in arch.get("domains", []) if d["kind"] != "REGION"]
    prev_dx = -prot_len
    for d in doms[:14]:
        col = DOMAIN_PALETTE.get(d["kind"], BLUE)
        w = max(prot_len * 0.012, d["end"] - d["start"] + 1)
        ax.add_patch(FancyBboxPatch((d["start"], -0.09), w, 0.18, boxstyle="round,pad=0.01", facecolor=col, edgecolor=NAVY, linewidth=0.7, alpha=0.9, zorder=3))
        xm = d["start"] + 0.5 * w
        if xm - prev_dx >= prot_len * 0.10:
            ax.text(xm, -0.14, d["label"].split(":", 1)[-1][:16], ha="center", va="top", fontsize=7.1, color=NAVY, fontweight="bold")
            prev_dx = xm

    for r in g_rows.itertuples(index=False):
        r_arch = json.loads(str(r.architecture_json))
        aa_pos = min(prot_len, max(1, int(r_arch.get("splice_junction_aa", 1))))
        dpsi = float(r.dpsi)
        col = RED if dpsi >= 0 else BLUE
        marker = "D" if r.frame_status == "frameshift" else "o"
        ax.plot([aa_pos, aa_pos], [0.09 if dpsi >= 0 else -0.09, dpsi], color=col, linewidth=1.8, zorder=4)
        ax.scatter([aa_pos], [dpsi], color=col, s=75, marker=marker, edgecolors=NAVY, linewidths=0.8, zorder=5)
        ha_lbl = "left" if aa_pos < prot_len * 0.08 else ("right" if aa_pos > prot_len * 0.92 else "center")
        x_lbl = aa_pos + (prot_len * 0.012 if ha_lbl == "left" else (-prot_len * 0.012 if ha_lbl == "right" else 0.0))
        ax.text(
            x_lbl,
            dpsi + (0.04 if dpsi >= 0 else -0.05),
            f"{r.event_type} @ aa {aa_pos} (ΔΨ={dpsi:+.2f}, {r.delta_aa:+d} aa)",
            ha=ha_lbl,
            va="bottom" if dpsi >= 0 else "top",
            fontsize=7.4,
            fontweight="bold",
            color=NAVY,
        )

    lim_y = max(0.58, float(g_rows["dpsi"].abs().max()) + 0.24)
    ax.set_ylim(-lim_y, lim_y)
    ax.set_xlim(-prot_len * 0.04, prot_len * 1.05)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlabel(f"{target_gene} Canonical Protein Coordinate (1 .. {prot_len:,} aa | UniProt {first['uniprot_acc']})", fontsize=9.2, color=NAVY)
    ax.set_ylabel("Splicing ΔΨ (+ Inclusion / - Skipping)", fontsize=9.2, color=NAVY)
    ax.set_title(f"Single-Protein Multi-Event Splicing Lollipop Domain Map: {target_gene} (Circle = In-frame, Diamond = Frameshift) | {contrast_id}", fontsize=10.8, fontweight="bold", color=NAVY)
    return _save_dual(fig, out_base)


# ---------------------------------------------------------------------------
# STYLE 5: 4-Stage Causal Alluvial Ribbon Flow (RBP -> AS Type -> Frame -> Proteome Fate)
# ---------------------------------------------------------------------------
def _clean_rbp_label(v: object) -> str:
    if isinstance(v, (list, tuple)):
        return str(v[0]) if v else "Other / Cis-Flank"
    s = str(v).strip()
    if s.startswith("[") and s.endswith("]"):
        parts = [p.strip(" '\"") for p in s[1:-1].split(",") if p.strip(" '\"")]
        return parts[0] if parts else "Other / Cis-Flank"
    return s if s else "Other / Cis-Flank"


def plot_causal_alluvial_flow_v24(
    prot_df: pd.DataFrame,
    conseq_df: pd.DataFrame,
    out_base: Path,
    contrast_id: str,
) -> tuple[Path, Path]:
    df = prot_df.copy()
    if not conseq_df.empty and "event_id" in conseq_df.columns and "RBP" in conseq_df.columns:
        c_tmp = conseq_df[["event_id", "RBP"]].copy()
        c_tmp["clean_rbp"] = c_tmp["RBP"].map(_clean_rbp_label)
        rbp_per_ev = c_tmp.groupby("event_id")["clean_rbp"].first().to_dict()
        df["rbp_group"] = df["event_id"].astype(str).map(rbp_per_ev).fillna("Other / Cis-Flank")
    else:
        df["rbp_group"] = "Regulated AS"

    # Keep top 5 RBPs + Other
    top_r = df.loc[df["rbp_group"].ne("Other / Cis-Flank"), "rbp_group"].value_counts().head(5).index.tolist()
    df["col0"] = np.where(df["rbp_group"].isin(top_r), df["rbp_group"], "Other Cis-Regulated")
    df["col1"] = df["event_type"].astype(str)
    df["col2"] = df["frame_status"].replace({
        "in_frame_indel": "In-Frame",
        "frameshift": "Frameshift",
        "utr_cis_regulatory_only": "UTR / Silent",
        "noncoding_transcript": "Non-coding",
        "translation_start_altered": "Start Altered",
    })
    df["col3"] = df.apply(classify_proteome_fate, axis=1).replace({
        "Productive UniProt/Pfam Domain Altered": "UniProt Domain Altered",
        "Productive Linker / Tail Isoform Switch": "Productive Linker Indel",
        "NMD Decay (Poison Exon Included)": "NMD Decay Switch",
        "NMD Decay (Essential Exon Skipped)": "NMD Decay Switch",
        "NMD Rescue (Poison Exon Skipped)": "NMD Protein Rescue",
        "ORF Preserved (Essential Exon Included)": "Productive Linker Indel",
        "5'/3' UTR or Non-Coding Switch": "UTR / Non-Coding",
    })

    fig = Figure(figsize=(13.2, 4.8), layout="constrained")
    ax = fig.subplots()
    ax.set_xlim(-0.25, 3.35)
    ax.set_ylim(-0.04, 1.08)
    ax.axis("off")

    cols = ["col0", "col1", "col2", "col3"]
    col_titles = ["1. Linked RBP Motif", "2. AS Event Type", "3. Reading Frame", "4. Proteome / NMD Fate"]
    node_colors = {
        "UniProt Domain Altered": BLUE,
        "Productive Linker Indel": "#7ea8cc",
        "NMD Decay Switch": RED,
        "NMD Protein Rescue": GREEN,
        "UTR / Non-Coding": "#9aa5b1",
        "In-Frame": TEAL,
        "Frameshift": RED,
        "UTR / Silent": "#9aa5b1",
        "Start Altered": GOLD,
        "Non-coding": "#b0b8c0",
    }
    palette_cycle = [BLUE, RED, TEAL, GREEN, PURPLE, GOLD, "#7a8694"]

    n_tot = max(len(df), 1)
    gap = 0.022
    node_boxes: list[dict[str, tuple[float, float, str]]] = []

    for ci, c_name in enumerate(cols):
        vc = df[c_name].value_counts()
        k_nodes = len(vc)
        avail = 1.0 - gap * max(k_nodes - 1, 0)
        y_cur = 1.0
        boxes_c: dict[str, tuple[float, float, str]] = {}
        for ni, (lbl, cnt) in enumerate(vc.items()):
            h = avail * (cnt / float(n_tot))
            y_bot = y_cur - h
            c_hex = node_colors.get(str(lbl), palette_cycle[ni % len(palette_cycle)])
            ax.add_patch(Rectangle((ci - 0.055, y_bot), 0.11, h, facecolor=c_hex, edgecolor=NAVY, linewidth=0.8, zorder=4))
            ha_dir = "right" if ci == 0 else ("left" if ci == 3 else "center")
            x_txt = ci - 0.075 if ci == 0 else (ci + 0.075 if ci == 3 else ci)
            if h >= 0.035 or ci in (0, 3):
                ax.text(
                    x_txt,
                    y_bot + 0.5 * h,
                    f"{lbl} ({cnt})",
                    ha=ha_dir,
                    va="center",
                    fontsize=7.8,
                    fontweight="bold",
                    color=NAVY,
                    bbox={"boxstyle": "round,pad=0.12", "facecolor": "white", "edgecolor": "none", "alpha": 0.78} if ci in (1, 2) else None,
                    zorder=6,
                )
            boxes_c[str(lbl)] = (y_bot, y_cur, c_hex)
            y_cur = y_bot - gap
        node_boxes.append(boxes_c)
        ax.text(ci, 1.04, col_titles[ci], ha="center", va="bottom", fontsize=9.5, fontweight="bold", color=NAVY)

    # Draw smooth cubic Bezier ribbons between adjacent columns
    for ci in range(3):
        c_left, c_right = cols[ci], cols[ci + 1]
        avail_l = 1.0 - gap * max(len(node_boxes[ci]) - 1, 0)
        avail_r = 1.0 - gap * max(len(node_boxes[ci + 1]) - 1, 0)
        cur_top_l = {k: v[1] for k, v in node_boxes[ci].items()}
        cur_top_r = {k: v[1] for k, v in node_boxes[ci + 1].items()}

        flows = df.groupby([c_left, c_right]).size().reset_index(name="n").sort_values("n", ascending=False)
        for r in flows.itertuples(index=False):
            kl, kr, cnt = str(r[0]), str(r[1]), int(r[2])
            hl = avail_l * (cnt / float(n_tot))
            hr = avail_r * (cnt / float(n_tot))
            yl1, yl0 = cur_top_l[kl], cur_top_l[kl] - hl
            yr1, yr0 = cur_top_r[kr], cur_top_r[kr] - hr
            cur_top_l[kl] = yl0
            cur_top_r[kr] = yr0

            x0, x1 = ci + 0.055, ci + 1 - 0.055
            xm = 0.5 * (x0 + x1)
            verts = [
                (x0, yl1), (xm, yl1), (xm, yr1), (x1, yr1),
                (x1, yr0), (xm, yr0), (xm, yl0), (x0, yl0), (x0, yl1),
            ]
            codes = [
                MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
                MplPath.LINETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4, MplPath.CLOSEPOLY,
            ]
            rib_col = node_boxes[ci + 1][kr][2] if ci == 2 else node_boxes[ci][kl][2]
            ax.add_patch(PathPatch(MplPath(verts, codes), facecolor=rib_col, edgecolor="none", alpha=0.30, zorder=2))

    fig.suptitle(f"End-to-End Regulatory-to-Proteome Causal Alluvial Flow | {contrast_id} (n = {len(df)} Regulated Events)", fontsize=11.2, fontweight="bold", color=NAVY)
    return _save_dual(fig, out_base)
