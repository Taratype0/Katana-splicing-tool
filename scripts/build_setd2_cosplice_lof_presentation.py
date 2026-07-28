from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def _read(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t", low_memory=False)


def _short_term(value: str) -> str:
    text = str(value)
    if 'C(Q("genotype"), Sum):C(Q("DP_status"), Sum)' in text:
        return "genotype x DP"
    if 'C(Q("genotype"), Sum):C(Q("treatment"), Sum)' in text:
        return "genotype x treatment"
    if 'C(Q("genotype"), Sum)' in text:
        return "genotype"
    return text


def build_method_figure(root: Path, overlap: pd.DataFrame, output_dir: Path) -> Path:
    metric = _read(root / "06c_method_benchmark" / "method_metric_ranking.tsv")
    programs = _read(root / "04d_hybrid_programs" / "hybrid_program_summary.tsv")
    association = _read(root / "05_association" / "residual" / "factorial_module_association.tsv")
    wgcna_association = _read(
        root / "exact_wgcna_full_stable_mut_effect" / "factorial_module_association.tsv"
    )

    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    methods = ["katana_sparse", "wgcna_exact"]
    metric_names = ["mean_bootstrap_ARI", "within_between_separation"]
    metric_labels = ["Bootstrap ARI", "Within/between separation"]
    x = np.arange(len(metric_names))
    width = 0.36
    colors = ["#176B87", "#C46B27"]
    for index, method in enumerate(methods):
        values = []
        for name in metric_names:
            row = metric.loc[metric["method"].eq(method) & metric["metric"].eq(name)]
            values.append(float(row["raw_value"].iloc[0]))
        axes[0, 0].bar(x + (index - 0.5) * width, values, width, label=method, color=colors[index])
        for xpos, value in zip(x + (index - 0.5) * width, values, strict=True):
            axes[0, 0].text(xpos, value + 0.012, f"{value:.3f}", ha="center", fontsize=9)
    axes[0, 0].set_xticks(x, metric_labels)
    axes[0, 0].set_ylim(0, 0.85)
    axes[0, 0].set_title("A. Same-input benchmark against exact WGCNA", loc="left", fontweight="bold")
    axes[0, 0].legend(frameon=False)
    axes[0, 0].grid(axis="y", alpha=0.2)

    tight = programs.loc[programs["program_type"].eq("tight")].copy()
    tight = tight.sort_values("overall_confidence", ascending=True)
    axes[0, 1].barh(tight["program"], tight["overall_confidence"], color="#2A9D8F")
    axes[0, 1].axvline(0.75, color="#D1495B", linestyle="--", linewidth=1, label="robust threshold")
    axes[0, 1].set_xlim(0, 1)
    axes[0, 1].set_xlabel("Internal confidence")
    axes[0, 1].set_title("B. Tight CoSplice program confidence", loc="left", fontweight="bold")
    axes[0, 1].legend(frameon=False, fontsize=8)
    axes[0, 1].grid(axis="x", alpha=0.2)

    genotype = association.loc[association["model_term"].astype(str).str.contains("genotype")].copy()
    genotype["FDR"] = pd.to_numeric(genotype["FDR"], errors="coerce")
    genotype = genotype.nsmallest(5, "FDR")
    genotype["source"] = "CoSplice tight residual"

    wgcna_genotype = wgcna_association.loc[
        wgcna_association["model_term"].astype(str).str.contains("genotype")
    ].copy()
    wgcna_genotype["FDR"] = pd.to_numeric(wgcna_genotype["FDR"], errors="coerce")
    wgcna_genotype = wgcna_genotype.loc[
        wgcna_genotype["FDR"].le(0.05)
        & wgcna_genotype["model_term"].astype(str).str.contains("DP_status")
        & ~wgcna_genotype["model_term"].astype(str).str.contains("treatment")
    ].nsmallest(2, "FDR")
    wgcna_genotype["source"] = "Exact WGCNA stable full"

    genotype = pd.concat([genotype, wgcna_genotype], ignore_index=True)
    genotype = genotype.sort_values("FDR", ascending=False)
    genotype["label"] = (
        genotype["source"].str.replace("Exact ", "", regex=False)
        + " | "
        + genotype["module"].astype(str)
        + " | "
        + genotype["model_term"].map(_short_term)
    )
    score = -np.log10(genotype["FDR"].clip(lower=1e-300))
    bar_colors = genotype["source"].map(
        {"CoSplice tight residual": "#6C5B7B", "Exact WGCNA stable full": "#C46B27"}
    )
    axes[1, 0].barh(genotype["label"], score, color=bar_colors)
    axes[1, 0].axvline(-np.log10(0.05), color="#D1495B", linestyle="--", linewidth=1, label="FDR 0.05")
    axes[1, 0].set_xlabel("-log10(FDR)")
    axes[1, 0].set_title("C. MUT terms depend on the network layer", loc="left", fontweight="bold")
    axes[1, 0].grid(axis="x", alpha=0.2)

    conditions = overlap["condition"].astype(str).tolist()
    mut_counts = pd.to_numeric(overlap["MUT_n"], errors="coerce")
    inhibitor_counts = pd.to_numeric(overlap["CLK3_n"], errors="coerce")
    pla_counts = pd.to_numeric(overlap["Pla_n"], errors="coerce")
    x = np.arange(len(conditions))
    width = 0.25
    axes[1, 1].bar(x - width, mut_counts, width, label="MUT vs WT", color="#3A6EA5")
    axes[1, 1].bar(x, inhibitor_counts, width, label="CLK3 vs baseline", color="#E38B2C")
    axes[1, 1].bar(x + width, pla_counts, width, label="Pla vs baseline", color="#2A9D8F")
    axes[1, 1].set_xticks(x, conditions)
    axes[1, 1].set_yscale("log")
    axes[1, 1].set_ylabel("Significant rMATS events (log scale)")
    axes[1, 1].set_title("D. MUT AS effect is narrower than inhibitor effects", loc="left", fontweight="bold")
    axes[1, 1].legend(frameon=False, fontsize=8)
    axes[1, 1].grid(axis="y", alpha=0.2)

    fig.suptitle("Setd2 CoSplice: method readiness and MUT effect", fontsize=17, fontweight="bold")
    fig.text(
        0.5,
        0.012,
        "Interpretation: exact WGCNA stable-full detects two genotype x DP modules; CoSplice tight-residual terms are suggestive. Both support a context-dependent, not global, MUT effect.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=[0.02, 0.04, 0.98, 0.95])
    path = output_dir / "Setd2_CoSplice_method_and_MUT_effect.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return path


def build_methods_only_figure(root: Path, output_dir: Path) -> Path:
    metric = _read(root / "06c_method_benchmark" / "method_metric_ranking.tsv")
    programs = _read(root / "04d_hybrid_programs" / "hybrid_program_summary.tsv")
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))

    axes[0, 0].axis("off")
    stages = [
        ("Stable AS input", "PSI + IJC/SJC\nreplicate/read QC"),
        ("WGCNA core", "soft adjacency + TOM\nmodule program scores"),
        ("CoSplice validation", "bootstrap co-assignment\ntight + diffuse programs"),
        ("Biological follow-up", "trait models + preservation\nRBP/DTU/isoform/sashimi"),
    ]
    stage_positions = [(0.25, 0.72), (0.75, 0.72), (0.25, 0.28), (0.75, 0.28)]
    for index, ((title, body), (xpos, ypos)) in enumerate(zip(stages, stage_positions, strict=True)):
        axes[0, 0].text(
            xpos, ypos, f"{title}\n{body}", ha="center", va="center", fontsize=9,
            bbox={"boxstyle": "round,pad=0.55", "facecolor": "#EEF5F3", "edgecolor": "#2A9D8F"},
        )
    axes[0, 0].annotate("", xy=(0.59, 0.72), xytext=(0.41, 0.72), arrowprops={"arrowstyle": "->", "color": "#52616B"})
    axes[0, 0].annotate("", xy=(0.25, 0.39), xytext=(0.75, 0.61), arrowprops={"arrowstyle": "->", "color": "#52616B"})
    axes[0, 0].annotate("", xy=(0.59, 0.28), xytext=(0.41, 0.28), arrowprops={"arrowstyle": "->", "color": "#52616B"})
    axes[0, 0].set_title("A. CoSplice analysis route", loc="left", fontweight="bold")

    methods = ["katana_sparse", "wgcna_exact"]
    names = ["mean_bootstrap_ARI", "within_between_separation"]
    labels = ["Bootstrap ARI", "Within/between separation"]
    x = np.arange(len(names))
    for index, (method, color) in enumerate(zip(methods, ["#176B87", "#C46B27"], strict=True)):
        values = [float(metric.loc[metric["method"].eq(method) & metric["metric"].eq(name), "raw_value"].iloc[0]) for name in names]
        pos = x + (index - 0.5) * 0.36
        axes[0, 1].bar(pos, values, 0.36, label=method, color=color)
        for xpos, value in zip(pos, values, strict=True):
            axes[0, 1].text(xpos, value + 0.012, f"{value:.3f}", ha="center", fontsize=9)
    axes[0, 1].set_xticks(x, labels)
    axes[0, 1].set_ylim(0, 0.85)
    axes[0, 1].set_title("B. Same-input benchmark against exact WGCNA", loc="left", fontweight="bold")
    axes[0, 1].legend(frameon=False)
    axes[0, 1].grid(axis="y", alpha=0.2)

    tight = programs.loc[programs["program_type"].eq("tight")].sort_values("overall_confidence")
    axes[1, 0].barh(tight["program"], tight["overall_confidence"], color="#2A9D8F")
    axes[1, 0].axvline(0.75, color="#D1495B", linestyle="--", linewidth=1, label="robust threshold")
    axes[1, 0].set_xlim(0, 1)
    axes[1, 0].set_xlabel("Internal confidence")
    axes[1, 0].set_title("C. Tight-program internal validation", loc="left", fontweight="bold")
    axes[1, 0].legend(frameon=False, fontsize=8)
    axes[1, 0].grid(axis="x", alpha=0.2)

    axes[1, 1].axis("off")
    axes[1, 1].set_title("D. Added value beyond standard WGCNA", loc="left", fontweight="bold")
    axes[1, 1].text(
        0.02, 0.96,
        "Splicing-specific evidence layers\n\n"
        "• exact rMATS event identity and event type\n"
        "• sample-level PSI plus IJC/SJC support\n"
        "• replicate-stability and measurement QC\n"
        "• bootstrap event co-assignment confidence\n"
        "• tight and diffuse program reporting\n"
        "• genotype/context/perturbation interaction models\n"
        "• module preservation and method benchmark\n"
        "• RBP, motif, DEG, DTU and isoform support\n"
        "• hub-event Native Sashimi validation",
        ha="left", va="top", fontsize=11, linespacing=1.35,
        bbox={"boxstyle": "round,pad=0.7", "facecolor": "#F3F5F7", "edgecolor": "#A8B3BD"},
    )
    fig.suptitle("CoSplice method: WGCNA core plus splicing-specific validation", fontsize=17, fontweight="bold")
    fig.tight_layout(rect=[0.02, 0.03, 0.98, 0.95])
    path = output_dir / "Setd2_CoSplice_methods.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return path


def build_mut_module_results_figure(root: Path, output_dir: Path) -> Path:
    wgcna_dir = root / "exact_wgcna_full_stable_mut_effect"
    assignments = _read(wgcna_dir / "event_module_assignments.tsv")
    event_meta = _read(root / "02_matrix" / "event_metadata.tsv")
    scores = _read(wgcna_dir / "module_program_scores.tsv")
    sample_meta = _read(root / "02_matrix" / "sample_metadata_resolved.tsv")
    psi = _read(root / "02_matrix" / "sample_level_psi_matrix.tsv").set_index("event_key")
    enrichment_path = wgcna_dir / "module_go_kegg_enrichment.tsv"
    enrichment = _read(enrichment_path) if enrichment_path.exists() else pd.DataFrame()
    modules = ["M03", "M06"]

    merged_scores = scores.merge(sample_meta, on="sample_id", how="left")
    condition_order = [
        f"{genotype}_{condition_token}_{treatment}" if treatment != "Baseline" else f"{genotype}_{condition_token}"
        for _, condition_token in [("noDP", "NO_DP"), ("DP", "DP")]
        for treatment in ["Baseline", "CLK3", "Pla"]
        for genotype in ["WT", "MUT"]
    ]
    grouped = merged_scores.groupby("condition", sort=False)[modules].mean()
    condition_order = [condition for condition in condition_order if condition in grouped.index]
    score_matrix = grouped.loc[condition_order, modules].T

    event_info = assignments.merge(
        event_meta[["event_key", "gene_id", "gene_symbol", "event_type"]], on="event_key", how="left"
    )
    component_rows = []
    for module in modules:
        frame = event_info.loc[event_info["wgcna_module"].eq(module)]
        counts = frame["event_type"].value_counts()
        component_rows.append({
            "module": module, "n_events": len(frame), "n_genes": frame["gene_id"].nunique(),
            **{event_type: int(counts.get(event_type, 0)) for event_type in ["SE", "RI", "A3SS", "A5SS", "MXE"]},
        })
    components = pd.DataFrame(component_rows)

    correlation_rows = []
    psi_meta = sample_meta.set_index("sample_id")
    for module in modules:
        event_ids = assignments.loc[assignments["wgcna_module"].eq(module), "event_key"]
        event_ids = psi.index.intersection(event_ids)
        for context in ["noDP", "DP"]:
            def group_mean(genotype: str, treatment: str) -> pd.Series:
                mask = psi_meta["genotype"].eq(genotype) & psi_meta["DP_status"].eq(context) & psi_meta["treatment"].eq(treatment)
                columns = psi_meta.index[mask].intersection(psi.columns)
                return psi.loc[event_ids, columns].mean(axis=1)

            wt = group_mean("WT", "Baseline")
            mut_effect = group_mean("MUT", "Baseline") - wt
            for perturbation in ["CLK3", "Pla"]:
                inhibitor_effect = group_mean("WT", perturbation) - wt
                valid = mut_effect.notna() & inhibitor_effect.notna()
                rho, pvalue = spearmanr(mut_effect[valid], inhibitor_effect[valid])
                directional = valid & mut_effect.abs().ge(0.05) & inhibitor_effect.abs().ge(0.05)
                same = np.sign(mut_effect[directional]).eq(np.sign(inhibitor_effect[directional]))
                correlation_rows.append({
                    "module": module, "context": context, "perturbation": perturbation,
                    "n_events": int(valid.sum()), "spearman_rho": float(rho), "spearman_pvalue": float(pvalue),
                    "direction_tested_n": int(directional.sum()), "same_direction_n": int(same.sum()),
                    "same_direction_pct": float(100 * same.mean()) if len(same) else np.nan,
                })
    correlations = pd.DataFrame(correlation_rows)

    score_matrix.to_csv(output_dir / "MUT_significant_module_condition_scores.tsv", sep="\t")
    components.to_csv(output_dir / "MUT_significant_module_components.tsv", sep="\t", index=False)
    correlations.to_csv(output_dir / "MUT_module_inhibitor_similarity.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(2, 2, figsize=(15, 9))
    image = axes[0, 0].imshow(score_matrix, cmap="RdBu_r", aspect="auto", vmin=-0.3, vmax=0.3)
    axes[0, 0].set_xticks(range(len(condition_order)), [x.replace("_", " ") for x in condition_order], rotation=55, ha="right", fontsize=8)
    axes[0, 0].set_yticks(range(len(modules)), modules)
    for row in range(len(modules)):
        for col in range(len(condition_order)):
            axes[0, 0].text(col, row, f"{score_matrix.iloc[row, col]:.2f}", ha="center", va="center", fontsize=7)
    fig.colorbar(image, ax=axes[0, 0], fraction=0.046, pad=0.04, label="Mean module score")
    axes[0, 0].set_title("A. M03/M06 program scores across all conditions", loc="left", fontweight="bold")

    corr_matrix = correlations.pivot(index="module", columns=["context", "perturbation"], values="spearman_rho")
    corr_image = axes[0, 1].imshow(corr_matrix, cmap="BrBG", aspect="auto", vmin=-0.5, vmax=0.5)
    corr_labels = [f"{context} {drug}" for context, drug in corr_matrix.columns]
    axes[0, 1].set_xticks(range(len(corr_labels)), corr_labels, rotation=35, ha="right")
    axes[0, 1].set_yticks(range(len(corr_matrix.index)), corr_matrix.index)
    for row in range(len(corr_matrix.index)):
        for col in range(len(corr_matrix.columns)):
            subset = correlations.loc[
                correlations["module"].eq(corr_matrix.index[row])
                & correlations["context"].eq(corr_matrix.columns[col][0])
                & correlations["perturbation"].eq(corr_matrix.columns[col][1])
            ].iloc[0]
            axes[0, 1].text(col, row, f"ρ {subset['spearman_rho']:.2f}\n{(subset['same_direction_pct']):.0f}% dir", ha="center", va="center", fontsize=8)
    fig.colorbar(corr_image, ax=axes[0, 1], fraction=0.046, pad=0.04, label="Event-level Spearman ρ")
    axes[0, 1].set_title("B. MUT module-event similarity to WT inhibitor responses", loc="left", fontweight="bold")

    bottom = np.zeros(len(components))
    colors = {"SE": "#4C78A8", "RI": "#F58518", "A3SS": "#54A24B", "A5SS": "#E45756", "MXE": "#8E6CBA"}
    for event_type in ["SE", "RI", "A3SS", "A5SS", "MXE"]:
        axes[1, 0].bar(components["module"], components[event_type], bottom=bottom, label=event_type, color=colors[event_type])
        bottom += components[event_type].to_numpy()
    for index, row in components.iterrows():
        axes[1, 0].text(index, row["n_events"] + 8, f"{row['n_events']} events\n{row['n_genes']} genes", ha="center", fontsize=9)
    axes[1, 0].set_ylim(0, float(components["n_events"].max()) * 1.18)
    axes[1, 0].set_ylabel("Stable AS events")
    axes[1, 0].set_title("C. Significant-module composition", loc="left", fontweight="bold")
    axes[1, 0].legend(frameon=False, ncol=5, fontsize=8)

    axes[1, 1].axis("off")
    m03_term = "regulation of RNA metabolic process (FDR 2.8e-12)"
    m06_term = "regulation of RNA splicing (FDR 2.6e-4)"
    if not enrichment.empty:
        m03 = enrichment.loc[enrichment["module"].eq("M03") & enrichment["term_id"].eq("GO:0051252")]
        m06 = enrichment.loc[enrichment["module"].eq("M06") & enrichment["term_id"].eq("GO:0043484")]
        if not m03.empty:
            m03_term = f"{m03.iloc[0]['term_name']} (FDR {float(m03.iloc[0]['FDR']):.1e})"
        if not m06.empty:
            m06_term = f"{m06.iloc[0]['term_name']} (FDR {float(m06.iloc[0]['FDR']):.1e})"
    axes[1, 1].set_title("D. Components, function and interpretation", loc="left", fontweight="bold")
    axes[1, 1].text(
        0.02, 0.97,
        f"M03\n• 402 events / 304 genes; RI-dominant\n• {m03_term}\n• recurrent genes include Snhg17, Mapk7 and RNA-regulatory genes\n\n"
        f"M06\n• 117 events / 104 genes; SE-dominant\n• {m06_term}\n• includes Clk1, Srsf3 and Srsf7\n\n"
        "Program-level conclusion\n"
        "• both modules have significant genotype × DP terms in exact WGCNA\n"
        "• inhibitor responses remain present in both WT and MUT; no simple occlusion\n"
        "• event-level MUT–inhibitor similarity is weak/moderate and context dependent\n"
        "• high 87–96% direction agreement applies only to shared rMATS-significant events",
        ha="left", va="top", fontsize=10.5, linespacing=1.3,
        bbox={"boxstyle": "round,pad=0.65", "facecolor": "#F3F5F7", "edgecolor": "#A8B3BD"},
    )
    fig.suptitle("Setd2 MUT program results: composition, function and inhibitor relationship", fontsize=17, fontweight="bold")
    fig.text(0.5, 0.012, "Direction percentage in panel B uses events with |mean ΔPSI| ≥ 0.05 on both sides; it is descriptive, not a causal test.", ha="center", fontsize=9)
    fig.tight_layout(rect=[0.02, 0.04, 0.98, 0.95])
    path = output_dir / "Setd2_MUT_significant_module_results.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return path


def build_lof_figure(overlap: pd.DataFrame, factor: pd.DataFrame, output_dir: Path) -> Path:
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    conditions = overlap["condition"].astype(str).tolist()
    x = np.arange(len(conditions))
    width = 0.24
    containment_columns = [
        ("MUT_in_CLK3_pct_of_MUT", "Contained in CLK3", "#E38B2C"),
        ("MUT_in_Pla_pct_of_MUT", "Contained in Pla", "#2A9D8F"),
        ("MUT_in_CLK3_or_Pla_pct_of_MUT", "Contained in either", "#3A6EA5"),
    ]
    for index, (column, label, color) in enumerate(containment_columns):
        values = pd.to_numeric(overlap[column], errors="coerce")
        axes[0, 0].bar(x + (index - 1) * width, values, width, label=label, color=color)
    axes[0, 0].set_xticks(x, conditions)
    axes[0, 0].set_ylim(0, 100)
    axes[0, 0].set_ylabel("Percent of MUT-significant events")
    axes[0, 0].set_title("A. rMATS-selected MUT signature containment", loc="left", fontweight="bold")
    axes[0, 0].legend(frameon=False, fontsize=8)
    axes[0, 0].grid(axis="y", alpha=0.2)

    for index, (column, label, color) in enumerate(
        [
            ("MUT_CLK3_same_direction_pct_of_overlap", "MUT-CLK3", "#E38B2C"),
            ("MUT_Pla_same_direction_pct_of_overlap", "MUT-Pla", "#2A9D8F"),
        ]
    ):
        values = pd.to_numeric(overlap[column], errors="coerce")
        axes[0, 1].bar(x + (index - 0.5) * 0.34, values, 0.34, label=label, color=color)
    axes[0, 1].set_xticks(x, conditions)
    axes[0, 1].set_ylim(0, 100)
    axes[0, 1].set_ylabel("Same direction within shared significant events (%)")
    axes[0, 1].set_title("B. Direction agreement is high after overlap selection", loc="left", fontweight="bold")
    axes[0, 1].legend(frameon=False)
    axes[0, 1].grid(axis="y", alpha=0.2)

    factor = factor.copy()
    factor["label"] = factor["context"].str.replace("DP_status=", "", regex=False) + " | " + factor["perturbation"]
    axes[1, 0].bar(factor["label"], factor["broad_signature_spearman"], color="#217A7A")
    axes[1, 0].axhline(0, color="#263238", linewidth=0.8)
    axes[1, 0].set_ylim(0, 0.35)
    axes[1, 0].set_ylabel("Spearman correlation across stable events")
    axes[1, 0].set_title("C. Transcriptome-wide similarity is weak", loc="left", fontweight="bold")
    axes[1, 0].tick_params(axis="x", rotation=35)
    axes[1, 0].grid(axis="y", alpha=0.2)
    for index, value in enumerate(factor["broad_signature_spearman"]):
        axes[1, 0].text(index, float(value) + 0.01, f"{float(value):.2f}", ha="center", fontsize=9)

    axes[1, 1].axis("off")
    conclusion = (
        "CURRENT EVIDENCE\n\n"
        "Supports:\n"
        "- a weak, context-dependent MUT splicing phenotype\n"
        "- partial overlap with broad CLK3/Pla responses in rMATS-selected events\n\n"
        "Does not support yet:\n"
        "- a robust transcriptome-wide inhibitor phenocopy\n"
        "- same-pathway occlusion (the one FDR-significant program interaction is opposite-response)\n"
        "- LOF classification from RNA-seq alone\n\n"
        "Decisive next evidence:\n"
        "H3K36me3 / SETD2 activity + WT rescue + mutant construct failure-to-rescue."
    )
    axes[1, 1].text(
        0.02,
        0.98,
        conclusion,
        ha="left",
        va="top",
        fontsize=11,
        linespacing=1.35,
        bbox={"boxstyle": "round,pad=0.8", "facecolor": "#F3F5F7", "edgecolor": "#A8B3BD"},
    )
    axes[1, 1].set_title("D. LOF evidence statement", loc="left", fontweight="bold")

    fig.suptitle("Does Setd2 MUT phenocopy CLK3/Pla splicing inhibition?", fontsize=17, fontweight="bold")
    fig.text(
        0.5,
        0.012,
        "Selected-overlap direction agreement is conditional on significance in both contrasts; it must be interpreted together with the weak all-event correlations.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=[0.02, 0.04, 0.98, 0.95])
    path = output_dir / "Setd2_MUT_CLK3_Pla_LOF_evidence.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description="Build slide-ready Setd2 CoSplice/LOF evidence figures.")
    parser.add_argument("--cosplice-root", type=Path, required=True)
    parser.add_argument("--overlap-summary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    overlap = _read(args.overlap_summary)
    factor = _read(args.cosplice_root / "06h_factor_contrasts" / "factor_evidence_summary.tsv")
    paths = [
        build_methods_only_figure(args.cosplice_root, args.output_dir),
        build_mut_module_results_figure(args.cosplice_root, args.output_dir),
        build_lof_figure(overlap, factor, args.output_dir),
    ]
    for path in paths:
        print(path)


if __name__ == "__main__":
    main()
