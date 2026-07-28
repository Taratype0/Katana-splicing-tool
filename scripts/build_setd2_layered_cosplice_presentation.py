from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from patsy import build_design_matrices
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests


COLORS = {
    "navy": "#183153",
    "blue": "#2F6BFF",
    "teal": "#148F77",
    "orange": "#E67E22",
    "red": "#C0392B",
    "gold": "#D4A017",
    "grey": "#667085",
    "light": "#F4F7FA",
    "line": "#C7D0DA",
}


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t", low_memory=False)


def save_figure(fig: plt.Figure, path: Path) -> Path:
    fig.savefig(path, dpi=240, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def bool_series(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.fillna(False)
    return series.astype(str).str.lower().isin({"true", "1", "yes"})


def add_panel_label(ax: plt.Axes, label: str, title: str) -> None:
    ax.set_title(f"{label}  {title}", loc="left", fontsize=13, fontweight="bold", pad=10)


def build_architecture_figure(output_dir: Path) -> Path:
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_axes([0.04, 0.48, 0.92, 0.43])
    ax.axis("off")
    stages = [
        ("1  Stable AS input", "rMATS events\nsample PSI + IJC/SJC\nreplicate/read QC"),
        ("2  Network discovery", "WGCNA soft adjacency + TOM\nplus sparse bootstrap consensus\ntight and diffuse programs"),
        ("3  Statistical association", "module eigengene/program score\nfactorial genotype x DP x treatment\nplanned contrasts + FDR"),
        ("4  Biological validation", "preservation and benchmark\nRBP/motif + DEG/DTU/isoform\nhub events + Native Sashimi"),
    ]
    xs = [0.02, 0.27, 0.52, 0.77]
    widths = [0.20] * 4
    for index in range(3):
        ax.annotate(
            "",
            xy=(xs[index + 1] - 0.015, 0.52),
            xytext=(xs[index] + widths[index] + 0.005, 0.52),
            arrowprops={"arrowstyle": "-|>", "lw": 2, "color": COLORS["grey"]},
        )
    for (title, body), xpos, width in zip(stages, xs, widths, strict=True):
        ax.text(
            xpos + width / 2,
            0.52,
            f"{title}\n\n{body}",
            ha="center",
            va="center",
            fontsize=11,
            color=COLORS["navy"],
            linespacing=1.35,
            bbox={
                "boxstyle": "round,pad=0.8",
                "facecolor": COLORS["light"],
                "edgecolor": COLORS["teal"],
                "linewidth": 1.8,
            },
        )
    ax.set_title("CoSplice architecture: WGCNA-compatible core with AS-specific evidence", fontsize=20, fontweight="bold")

    ax2 = fig.add_axes([0.05, 0.08, 0.90, 0.33])
    ax2.axis("off")
    columns = ["Event DAS", "Sample PSI QC", "Network modules", "Bootstrap stability", "Factorial model", "Functional support", "Preservation"]
    rows = ["rMATS", "Jutils", "WGCNA", "CoSplice"]
    values = [
        ["Yes", "Limited", "No", "No", "Event-level", "No", "No"],
        ["No", "Yes", "No", "No", "No", "No", "No"],
        ["No", "Input-dependent", "Yes", "Limited", "Trait association", "Generic", "Yes"],
        ["Yes", "Yes", "Yes", "Yes", "Yes", "AS-specific", "Yes"],
    ]
    table = ax2.table(cellText=values, rowLabels=rows, colLabels=columns, cellLoc="center", loc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(9.5)
    table.scale(1, 2.0)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("white")
        if row == 0:
            cell.set_facecolor(COLORS["navy"])
            cell.get_text().set_color("white")
            cell.get_text().set_fontweight("bold")
        elif col == -1:
            cell.set_facecolor("#DDE8F2" if row < 4 else "#CFEDE6")
            cell.get_text().set_fontweight("bold")
        elif row == 4:
            cell.set_facecolor("#E4F4EF")
        else:
            cell.set_facecolor("#F4F6F8")
    fig.text(
        0.5,
        0.025,
        "CoSplice does not replace rMATS significance testing; it organizes stable AS events into validated sample-level programs and integrates orthogonal support.",
        ha="center",
        fontsize=11,
        color=COLORS["grey"],
    )
    return save_figure(fig, output_dir / "01_CoSplice_architecture_and_method_boundary.png")


def build_event_overview_figure(overlap: pd.DataFrame, output_dir: Path) -> Path:
    fig, axes = plt.subplots(2, 2, figsize=(15, 9))
    conditions = overlap["condition"].astype(str).tolist()
    x = np.arange(len(conditions))

    width = 0.24
    count_specs = [
        ("MUT_n", "MUT vs WT", COLORS["blue"]),
        ("CLK3_n", "CLK3 vs baseline", COLORS["orange"]),
        ("Pla_n", "Pla vs baseline", COLORS["teal"]),
    ]
    for index, (column, label, color) in enumerate(count_specs):
        values = pd.to_numeric(overlap[column], errors="coerce")
        positions = x + (index - 1) * width
        axes[0, 0].bar(positions, values, width, label=label, color=color)
        for xpos, value in zip(positions, values, strict=True):
            axes[0, 0].text(xpos, value * 1.08, f"{int(value):,}", ha="center", fontsize=8)
    axes[0, 0].set_yscale("log")
    axes[0, 0].set_xticks(x, conditions)
    axes[0, 0].set_ylabel("Significant rMATS events (log scale)")
    axes[0, 0].legend(frameon=False, fontsize=8)
    axes[0, 0].grid(axis="y", alpha=0.2)
    add_panel_label(axes[0, 0], "A", "MUT signature is narrower than inhibitor responses")

    containment_specs = [
        ("MUT_in_CLK3_pct_of_MUT", "Contained in CLK3", COLORS["orange"]),
        ("MUT_in_Pla_pct_of_MUT", "Contained in Pla", COLORS["teal"]),
        ("MUT_in_CLK3_or_Pla_pct_of_MUT", "Contained in either", COLORS["blue"]),
    ]
    for index, (column, label, color) in enumerate(containment_specs):
        values = pd.to_numeric(overlap[column], errors="coerce")
        positions = x + (index - 1) * width
        axes[0, 1].bar(positions, values, width, label=label, color=color)
        for xpos, value in zip(positions, values, strict=True):
            axes[0, 1].text(xpos, value + 1.5, f"{value:.1f}%", ha="center", fontsize=8)
    axes[0, 1].set_ylim(0, 55)
    axes[0, 1].set_xticks(x, conditions)
    axes[0, 1].set_ylabel("Percent of MUT-significant events")
    axes[0, 1].legend(frameon=False, fontsize=8)
    axes[0, 1].grid(axis="y", alpha=0.2)
    add_panel_label(axes[0, 1], "B", "A coherent subset of MUT events is inhibitor-associated")

    direction_specs = [
        ("MUT_CLK3_same_direction_pct_of_overlap", "MUT–CLK3", COLORS["orange"]),
        ("MUT_Pla_same_direction_pct_of_overlap", "MUT–Pla", COLORS["teal"]),
    ]
    for index, (column, label, color) in enumerate(direction_specs):
        values = pd.to_numeric(overlap[column], errors="coerce")
        positions = x + (index - 0.5) * 0.34
        axes[1, 0].bar(positions, values, 0.34, label=label, color=color)
        for xpos, value in zip(positions, values, strict=True):
            axes[1, 0].text(xpos, value + 1.5, f"{value:.1f}%", ha="center", fontsize=9, fontweight="bold")
    axes[1, 0].axhline(50, color=COLORS["grey"], linestyle="--", linewidth=1, label="50% reference")
    axes[1, 0].set_ylim(0, 105)
    axes[1, 0].set_xticks(x, conditions)
    axes[1, 0].set_ylabel("Same direction among shared significant events")
    axes[1, 0].legend(frameon=False, fontsize=8)
    axes[1, 0].grid(axis="y", alpha=0.2)
    add_panel_label(axes[1, 0], "C", "Shared significant events have strong direction agreement")

    labels = []
    correlations = []
    deltas = []
    colors = []
    for _, row in overlap.iterrows():
        condition = str(row["condition"])
        for drug, color in [("CLK3", COLORS["orange"]), ("Pla", COLORS["teal"])]:
            labels.append(f"{condition}\n{drug}")
            correlations.append(float(row[f"MUT_{drug}_dPSI_pearson_overlap"]))
            deltas.append(float(row[f"MUT_{drug}_median_abs_delta_dPSI_overlap"]))
            colors.append(color)
    positions = np.arange(len(labels))
    axes[1, 1].bar(positions, correlations, color=colors)
    for xpos, corr, delta in zip(positions, correlations, deltas, strict=True):
        axes[1, 1].text(xpos, corr + 0.025, f"r={corr:.2f}\nΔ={delta:.3f}", ha="center", fontsize=8)
    axes[1, 1].set_ylim(0, 1)
    axes[1, 1].set_xticks(positions, labels)
    axes[1, 1].set_ylabel("Pearson r of dPSI in shared significant events")
    axes[1, 1].grid(axis="y", alpha=0.2)
    add_panel_label(axes[1, 1], "D", "Shared-core effect sizes are strongly correlated")

    fig.suptitle("Layer 1 — Primary event-level evidence", fontsize=20, fontweight="bold", color=COLORS["navy"])
    fig.text(
        0.5,
        0.012,
        "Interpretation: MUT produces a restricted signature embedded within broader inhibitor responses; high agreement applies to exact events significant in both contrasts.",
        ha="center",
        fontsize=10.5,
        color=COLORS["grey"],
    )
    fig.tight_layout(rect=[0.02, 0.04, 0.98, 0.95])
    return save_figure(fig, output_dir / "02_event_level_evidence_overview.png")


def load_shared_pair(path: Path) -> pd.DataFrame:
    frame = read_tsv(path)
    frame["significant_A"] = bool_series(frame["significant_A"])
    frame["significant_B"] = bool_series(frame["significant_B"])
    frame["dPSI_A"] = pd.to_numeric(frame["dPSI_A"], errors="coerce")
    frame["dPSI_B"] = pd.to_numeric(frame["dPSI_B"], errors="coerce")
    return frame.loc[frame["significant_A"] & frame["significant_B"]].dropna(subset=["dPSI_A", "dPSI_B"])


def build_shared_core_scatter_figure(workdir: Path, output_dir: Path) -> Path:
    specs = [
        ("Mutation Effect vs. CLK3 Effect (DP).tsv", "DP", "CLK3"),
        ("Mutation Effect vs. Pla Effect (DP).tsv", "DP", "Pla"),
        ("Mutation Effect vs. CLK3 Effect (noDP).tsv", "noDP", "CLK3"),
        ("Mutation Effect vs. Pla Effect (noDP).tsv", "noDP", "Pla"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 11), sharex=True, sharey=True)
    event_colors = {"SE": "#2F6BFF", "RI": "#E67E22", "A3SS": "#148F77", "A5SS": "#C0392B", "MXE": "#7D5FFF"}
    summary_rows = []
    for ax, (filename, context, drug) in zip(axes.flat, specs, strict=True):
        frame = load_shared_pair(workdir / filename)
        same = np.sign(frame["dPSI_A"]).eq(np.sign(frame["dPSI_B"]))
        corr, pvalue = pearsonr(frame["dPSI_A"], frame["dPSI_B"])
        for event_type, subset in frame.groupby("event_type"):
            ax.scatter(
                subset["dPSI_A"], subset["dPSI_B"], s=24, alpha=0.65,
                color=event_colors.get(str(event_type), COLORS["grey"]), label=str(event_type),
                edgecolor="white", linewidth=0.25,
            )
        ax.axline((0, 0), slope=1, color=COLORS["grey"], linestyle="--", linewidth=1)
        ax.axhline(0, color=COLORS["line"], linewidth=0.8)
        ax.axvline(0, color=COLORS["line"], linewidth=0.8)
        ax.set_xlim(-1.05, 1.05)
        ax.set_ylim(-1.05, 1.05)
        ax.set_aspect("equal", adjustable="box")
        ax.grid(alpha=0.12)
        ax.text(
            0.04, 0.96,
            f"n = {len(frame):,}\nsame direction = {100 * same.mean():.1f}%\nPearson r = {corr:.3f}\nmedian |ΔdPSI| = {(frame['dPSI_A'] - frame['dPSI_B']).abs().median():.3f}",
            transform=ax.transAxes, ha="left", va="top", fontsize=9.5,
            bbox={"boxstyle": "round,pad=0.45", "facecolor": "white", "edgecolor": COLORS["line"], "alpha": 0.92},
        )
        add_panel_label(ax, context, f"MUT vs {drug} shared significant events")
        summary_rows.append({
            "context": context, "inhibitor": drug, "n_shared_significant": len(frame),
            "same_direction_pct": 100 * same.mean(), "pearson_r": corr, "pearson_pvalue": pvalue,
            "median_abs_delta_dPSI": (frame["dPSI_A"] - frame["dPSI_B"]).abs().median(),
        })
    axes[1, 0].set_xlabel("MUT vs WT standardized dPSI")
    axes[1, 1].set_xlabel("MUT vs WT standardized dPSI")
    axes[0, 0].set_ylabel("Inhibitor vs baseline standardized dPSI")
    axes[1, 0].set_ylabel("Inhibitor vs baseline standardized dPSI")
    handles = [plt.Line2D([0], [0], marker="o", linestyle="", color=color, label=event_type) for event_type, color in event_colors.items()]
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, bbox_to_anchor=(0.5, 0.005))
    fig.suptitle("Layer 2 — Exact shared significant AS-event core", fontsize=20, fontweight="bold", color=COLORS["navy"])
    fig.tight_layout(rect=[0.02, 0.045, 0.98, 0.95])
    pd.DataFrame(summary_rows).to_csv(output_dir / "shared_significant_event_core_summary.tsv", sep="\t", index=False)
    return save_figure(fig, output_dir / "03_shared_core_event_effect_scatter.png")


def _design_mean(model, genotype: str, dp_status: str, batches: list[str]) -> np.ndarray:
    new = pd.DataFrame({"genotype": genotype, "DP_status": dp_status, "batch": batches})
    design = build_design_matrices([model.model.data.design_info], new)[0]
    return np.asarray(design).mean(axis=0)


def formal_baseline_module_contrasts(root: Path) -> pd.DataFrame:
    wgcna_dir = root / "exact_wgcna_full_stable_mut_effect"
    scores = read_tsv(wgcna_dir / "module_program_scores.tsv")
    metadata = read_tsv(root / "02_matrix" / "sample_metadata_resolved.tsv")
    frame = metadata.merge(scores, on="sample_id", how="inner")
    frame = frame.loc[frame["treatment"].eq("Baseline")].copy()
    frame["batch"] = frame["batch"].astype(str)
    batches = sorted(frame["batch"].unique().tolist())
    rows = []
    for module in ["M03", "M06"]:
        model = smf.ols(f"{module} ~ C(genotype) * C(DP_status) + C(batch)", data=frame).fit()
        vectors = {
            (genotype, dp): _design_mean(model, genotype, dp, batches)
            for genotype in ["WT", "MUT"] for dp in ["noDP", "DP"]
        }
        contrasts = {
            "Baseline noDP: MUT - WT": vectors[("MUT", "noDP")] - vectors[("WT", "noDP")],
            "Baseline DP: MUT - WT": vectors[("MUT", "DP")] - vectors[("WT", "DP")],
            "Baseline genotype x DP": (vectors[("MUT", "DP")] - vectors[("WT", "DP")]) - (vectors[("MUT", "noDP")] - vectors[("WT", "noDP")]),
        }
        for name, vector in contrasts.items():
            result = model.t_test(vector)
            ci = np.asarray(result.conf_int()).ravel()
            estimate = np.asarray(result.effect).item()
            standard_error = np.asarray(result.sd).item()
            pvalue = np.asarray(result.pvalue).item()
            rows.append({
                "module": module, "contrast": name, "estimate": float(estimate),
                "standard_error": float(standard_error), "ci_low": float(ci[0]), "ci_high": float(ci[1]),
                "pvalue": float(pvalue), "n_samples": int(model.nobs), "residual_df": float(model.df_resid),
            })
    table = pd.DataFrame(rows)
    table["FDR"] = multipletests(table["pvalue"], method="fdr_bh")[1]
    return table


def build_module_supplement_figure(root: Path, output_dir: Path) -> Path:
    wgcna = read_tsv(root / "exact_wgcna_full_stable_mut_effect" / "factorial_module_association.tsv")
    scores = read_tsv(root / "exact_wgcna_full_stable_mut_effect" / "module_program_scores.tsv")
    metadata = read_tsv(root / "02_matrix" / "sample_metadata_resolved.tsv")
    merged = metadata.merge(scores, on="sample_id", how="inner")
    formal = formal_baseline_module_contrasts(root)
    formal.to_csv(output_dir / "module_baseline_planned_contrasts.tsv", sep="\t", index=False)

    genotype_term = 'C(Q("genotype"), Sum)'
    interaction_term = 'C(Q("genotype"), Sum):C(Q("DP_status"), Sum)'
    result_rows = []
    for module in ["M03", "M06"]:
        for label, term in [("MUT main effect", genotype_term), ("MUT x DP", interaction_term)]:
            row = wgcna.loc[wgcna["module"].eq(module) & wgcna["model_term"].eq(term)].iloc[0]
            result_rows.append({"module": module, "test": label, "FDR": float(row["FDR"]), "partial_eta_squared": float(row["partial_eta_squared"])})
    factorial = pd.DataFrame(result_rows)
    factorial.to_csv(output_dir / "module_factorial_key_results.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(2, 2, figsize=(15, 9))
    x = np.arange(2)
    for index, (test, color) in enumerate([("MUT main effect", COLORS["grey"]), ("MUT x DP", COLORS["blue"]) ]):
        values = factorial.loc[factorial["test"].eq(test)].set_index("module").loc[["M03", "M06"], "FDR"]
        score = -np.log10(values.clip(lower=1e-300))
        positions = x + (index - 0.5) * 0.34
        axes[0, 0].bar(positions, score, 0.34, color=color, label=test)
        for xpos, fdr in zip(positions, values, strict=True):
            axes[0, 0].text(xpos, -np.log10(fdr) + 0.08, f"FDR {fdr:.3f}", ha="center", fontsize=8)
    axes[0, 0].axhline(-np.log10(0.05), color=COLORS["red"], linestyle="--", linewidth=1, label="FDR 0.05")
    axes[0, 0].set_xticks(x, ["M03", "M06"])
    axes[0, 0].set_ylabel("-log10(FDR)")
    axes[0, 0].legend(frameon=False, fontsize=8)
    axes[0, 0].grid(axis="y", alpha=0.2)
    add_panel_label(axes[0, 0], "A", "No global MUT main effect; significant MUT x DP interaction")

    offsets = {"Baseline noDP: MUT - WT": -0.12, "Baseline DP: MUT - WT": 0.12}
    colors = {"Baseline noDP: MUT - WT": COLORS["orange"], "Baseline DP: MUT - WT": COLORS["teal"]}
    for contrast, offset in offsets.items():
        subset = formal.loc[formal["contrast"].eq(contrast)].set_index("module").loc[["M03", "M06"]]
        ypos = x + offset
        axes[0, 1].errorbar(
            subset["estimate"], ypos,
            xerr=[subset["estimate"] - subset["ci_low"], subset["ci_high"] - subset["estimate"]],
            fmt="o", capsize=4, color=colors[contrast], label=contrast.replace("Baseline ", ""),
        )
        for estimate, y, fdr in zip(subset["estimate"], ypos, subset["FDR"], strict=True):
            axes[0, 1].text(estimate, y + 0.06, f"FDR {fdr:.2f}", ha="center", fontsize=8)
    axes[0, 1].axvline(0, color=COLORS["grey"], linewidth=1)
    axes[0, 1].set_yticks(x, ["M03", "M06"])
    axes[0, 1].set_xlabel("Estimated baseline MUT - WT module score (95% CI)")
    axes[0, 1].legend(frameon=False, fontsize=8)
    axes[0, 1].grid(axis="x", alpha=0.2)
    add_panel_label(axes[0, 1], "B", "Baseline simple effects are not significant")

    grouped = merged.groupby(["genotype", "DP_status", "treatment"])[["M03", "M06"]].mean()
    for ax, module in zip(axes[1], ["M03", "M06"], strict=True):
        labels = []
        wt_effects = []
        mut_effects = []
        for dp in ["noDP", "DP"]:
            for drug in ["CLK3", "Pla"]:
                labels.append(f"{dp}\n{drug}")
                wt_effects.append(grouped.loc[("WT", dp, drug), module] - grouped.loc[("WT", dp, "Baseline"), module])
                mut_effects.append(grouped.loc[("MUT", dp, drug), module] - grouped.loc[("MUT", dp, "Baseline"), module])
        pos = np.arange(len(labels))
        ax.bar(pos - 0.18, wt_effects, 0.36, color="#8EA4BA", label="Effect in WT")
        ax.bar(pos + 0.18, mut_effects, 0.36, color=COLORS["blue"], label="Effect in MUT")
        ax.axhline(0, color=COLORS["grey"], linewidth=0.8)
        ax.set_xticks(pos, labels)
        ax.set_ylabel("Inhibitor - matched baseline module score")
        ax.legend(frameon=False, fontsize=8)
        ax.grid(axis="y", alpha=0.2)
        add_panel_label(ax, "C" if module == "M03" else "D", f"{module} inhibitor response in WT and MUT")

    fig.suptitle("Layer 3 — Module-level supplementary evidence", fontsize=20, fontweight="bold", color=COLORS["navy"])
    fig.text(
        0.5, 0.012,
        "Modules summarize coordinated PSI programs. They support context dependence but do not replace event-level MUT–inhibitor overlap or prove a baseline MUT phenotype.",
        ha="center", fontsize=10.5, color=COLORS["grey"],
    )
    fig.tight_layout(rect=[0.02, 0.04, 0.98, 0.95])
    return save_figure(fig, output_dir / "04_module_level_supplementary_evidence.png")


def build_evidence_ladder_figure(output_dir: Path) -> Path:
    fig, ax = plt.subplots(figsize=(16, 9))
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_title("Layered conclusion: what is supported, and what remains unproven", fontsize=20, fontweight="bold", color=COLORS["navy"], pad=20)
    boxes = [
        (0.035, 0.53, 0.285, 0.34, "PRIMARY EVIDENCE", COLORS["teal"],
         "Exact shared significant AS events\n\n"
         "• 35.5–43.5% of MUT events occur\n  under at least one inhibitor\n"
         "• 86.7–96.4% same direction\n"
         "• shared-core dPSI r = 0.745–0.909\n\n"
         "Conclusion: coherent event-level overlap"),
        (0.3575, 0.53, 0.285, 0.34, "SUPPLEMENTARY EVIDENCE", COLORS["blue"],
         "Co-splicing programs\n\n"
         "• M03/M06: significant genotype × DP\n"
         "• no global MUT main effect\n"
         "• baseline DP/noDP effects not significant\n"
         "• inhibitor program responses remain strong\n\n"
         "Conclusion: context-dependent support"),
        (0.68, 0.53, 0.285, 0.34, "BOUNDARY", COLORS["red"],
         "Not established by RNA-seq alone\n\n"
         "• complete inhibitor phenocopy\n"
         "• same-pathway occlusion\n"
         "• SETD2 molecular loss of function\n\n"
         "Needed: H3K36me3/activity, WT rescue,\n"
         "and mutant failure-to-rescue"),
    ]
    for x, y, width, height, title, color, body in boxes:
        ax.add_patch(plt.Rectangle((x, y), width, height, facecolor="white", edgecolor=color, linewidth=2.2))
        ax.add_patch(plt.Rectangle((x, y + height - 0.055), width, 0.055, facecolor=color, edgecolor=color))
        ax.text(x + 0.015, y + height - 0.028, title, va="center", ha="left", fontsize=12, fontweight="bold", color="white")
        ax.text(x + 0.018, y + height - 0.085, body, va="top", ha="left", fontsize=10.2, color=COLORS["navy"], linespacing=1.3)
    ax.text(
        0.5, 0.315,
        "Best current statement",
        ha="center", fontsize=14, fontweight="bold", color=COLORS["navy"],
    )
    ax.text(
        0.5, 0.19,
        "Setd2 MUT produces a restricted but highly coherent AS signature embedded within broader CLK3/Pla responses.\n"
        "This supports partial molecular convergence, while the module layer indicates DP-dependent context rather than a strong global MUT effect.",
        ha="center", va="center", fontsize=14, color=COLORS["navy"], linespacing=1.5,
        bbox={"boxstyle": "round,pad=0.9", "facecolor": COLORS["light"], "edgecolor": COLORS["navy"], "linewidth": 1.5},
    )
    return save_figure(fig, output_dir / "05_evidence_ladder_and_conclusion.png")


def write_report(output_dir: Path) -> Path:
    path = output_dir / "SETD2_COSPLICE_LAYERED_RESULTS.md"
    path.write_text(
        """# Setd2 CoSplice layered result presentation

## Presentation order

1. `01_CoSplice_architecture_and_method_boundary.png`: method architecture and boundaries versus rMATS, Jutils and WGCNA.
2. `02_event_level_evidence_overview.png`: primary event-level evidence.
3. `03_shared_core_event_effect_scatter.png`: exact shared significant event effects and direction.
4. `04_module_level_supplementary_evidence.png`: module-level statistics presented explicitly as supplementary evidence.
5. `05_evidence_ladder_and_conclusion.png`: evidence hierarchy and safe conclusion.

## Recommended presentation narrative

1. **Method first:** CoSplice starts from rMATS event-level evidence and adds sample-level PSI program analysis; it does not replace the rMATS significance test.
2. **Primary result:** quantify how much of the MUT-significant AS signature is contained in CLK3/Pla responses.
3. **Direction and effect:** among exact events significant in both contrasts, show sign agreement and standardized-dPSI correlation.
4. **Supplementary result:** use M03/M06 module scores to ask whether coordinated PSI programs depend on DP context.
5. **Boundary:** report partial molecular convergence, not complete inhibitor phenocopy or proven SETD2 loss of function.

## Key Setd2 results

| Context | Inhibitor | Exact shared significant events | Same direction | Pearson r of standardized dPSI |
|---|---:|---:|---:|---:|
| DP | CLK3 | 256 | 86.7% | 0.745 |
| DP | Pla | 218 | 92.2% | 0.822 |
| noDP | CLK3 | 191 | 93.7% | 0.863 |
| noDP | Pla | 110 | 96.4% | 0.909 |

- The MUT signature is much narrower than the inhibitor signatures.
- 35.5% of DP MUT-significant events and 43.5% of noDP MUT-significant events are significant under at least one inhibitor.
- High directional agreement applies to the **exact shared significant subset**, not to all tested events.
- M03/M06 show significant genotype-by-DP terms in the full factorial model (FDR 0.018 and 0.030), while their global genotype main effects are not significant.
- Baseline DP/noDP MUT-vs-WT planned contrasts are not significant after FDR correction.

## Statistical hierarchy

### Primary result: exact event-level evidence

- Event identity is matched by exact rMATS event key, not gene alone.
- Significant events use the configured rMATS FDR and absolute dPSI criteria.
- Containment asks what fraction of MUT-significant events are also significant under CLK3 and/or Pla.
- Direction agreement and dPSI correlation are calculated only among exact events significant in both contrasts.
- These selected-overlap statistics must not be described as transcriptome-wide similarity.

### Supplementary result: module-level evidence

- Modules summarize sample-level PSI patterns using a WGCNA-compatible eigengene/program score.
- The full factorial model tests genotype, DP, treatment and their interactions while adjusting for batch.
- M03/M06 have significant genotype-by-DP terms, not significant global genotype main effects.
- Baseline MUT-vs-WT planned contrasts within DP/noDP are shown separately and are not significant.
- Therefore M03/M06 support context dependence but are not evidence of a strong baseline MUT module phenotype.

## Scientific conclusion

Setd2 MUT produces a relatively restricted AS signature. A coherent subset is shared with broad CLK3/Pla responses and has strong direction and effect-size agreement. This supports partial pathway-level convergence. It does not establish complete inhibitor phenocopy, pathway occlusion or SETD2 molecular loss of function. Direct H3K36me3/activity and rescue experiments remain necessary.

## CoSplice functional scope

- stable AS event selection using PSI plus IJC/SJC and replicate/read QC
- WGCNA soft-adjacency/TOM diffuse modules
- sparse bootstrap consensus programs
- module eigengene/program scores
- factorial association and planned contrasts
- module preservation and method benchmarking
- RBP/motif, DEG, DTU and isoform support
- hub-event and Native Sashimi follow-up

CoSplice extends rather than replaces rMATS. rMATS supplies event-level DAS evidence; CoSplice supplies sample-level program discovery, stability, association and integrated validation.

## Difference from existing methods

| Method | Primary role | What CoSplice adds or preserves |
|---|---|---|
| rMATS | Event-level differential splicing and read-count evidence | rMATS remains the primary DAS test; CoSplice organizes stable events into sample-level programs and follow-up layers. |
| Jutils | PSI-oriented visualization, heatmap and PCA utilities | CoSplice adds explicit event stability/read QC, modules, bootstrap consensus, factorial association, preservation and orthogonal support. |
| WGCNA | General soft-adjacency/TOM network and module analysis | CoSplice preserves a WGCNA-compatible diffuse-module path and adds AS-event identity, PSI/read uncertainty, sparse consensus, and splicing-specific validation. |

This is a functional boundary comparison, not a claim that one method is universally superior. CoSplice combines established components into a splicing-specific workflow and must be benchmarked and externally validated as a complete method.

## Method references

- Shen et al. rMATS, *PNAS* (2014), DOI: 10.1073/pnas.1419161111.
- Langfelder and Horvath. WGCNA, *BMC Bioinformatics* (2008), DOI: 10.1186/1471-2105-9-559.
- Love et al. DESeq2, *Genome Biology* (2014), DOI: 10.1186/s13059-014-0550-8.
- Jutils workflow article, PMCID: PMC9502236.
""",
        encoding="utf-8",
    )
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a layered Setd2 CoSplice result presentation.")
    parser.add_argument("--cosplice-root", type=Path, required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--overlap-summary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    overlap = read_tsv(args.overlap_summary)
    paths = [
        build_architecture_figure(args.output_dir),
        build_event_overview_figure(overlap, args.output_dir),
        build_shared_core_scatter_figure(args.workdir, args.output_dir),
        build_module_supplement_figure(args.cosplice_root, args.output_dir),
        build_evidence_ladder_figure(args.output_dir),
        write_report(args.output_dir),
    ]
    for path in paths:
        print(path)


if __name__ == "__main__":
    main()
