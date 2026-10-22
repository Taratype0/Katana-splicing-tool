from __future__ import annotations

import argparse
from dataclasses import dataclass
from itertools import combinations
import html
import json
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import binomtest, spearmanr
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests
import yaml


CONTRASTS = (
    ("KO_DP", "WT_DP"),
    ("KO_NODP", "WT_NODP"),
    ("KO_DP", "KO_NODP"),
    ("WT_DP", "WT_NODP"),
)


@dataclass(frozen=True)
class MatrixResult:
    label: str
    root: Path
    matrix: pd.DataFrame
    event_qc: pd.DataFrame
    metadata: pd.DataFrame


def _read_tsv(path: Path, **kwargs: object) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t", low_memory=False, **kwargs)


def _load_matrix_result(label: str, root: Path) -> MatrixResult:
    matrix = _read_tsv(root / "03_stable_as_qc" / "stable_psi_matrix.tsv")
    matrix = matrix.set_index("event_key")
    matrix = matrix.apply(pd.to_numeric, errors="coerce")
    event_qc = _read_tsv(root / "03_stable_as_qc" / "stable_as_event_qc.tsv")
    metadata = _read_tsv(root / "02_matrix" / "sample_metadata_resolved.tsv")
    metadata["sample_id"] = metadata["sample_id"].astype(str)
    metadata["condition"] = metadata["condition"].astype(str)
    missing = set(metadata["sample_id"]).difference(matrix.columns)
    if missing:
        raise ValueError(f"{label} matrix is missing metadata samples: {sorted(missing)}")
    return MatrixResult(label, root, matrix, event_qc, metadata)


def _safe_corr(values: pd.DataFrame, method: str = "spearman") -> pd.DataFrame:
    return values.corr(method=method, min_periods=max(3, len(values) // 10))


def _sample_similarity(result: MatrixResult) -> tuple[pd.DataFrame, pd.DataFrame]:
    samples = list(result.metadata["sample_id"])
    sample_by_event = result.matrix.loc[:, samples].T
    usable = sample_by_event.loc[:, sample_by_event.notna().sum(axis=0).ge(3)]
    usable = usable.fillna(usable.median(axis=0))
    variable = usable.loc[:, usable.var(axis=0).gt(0)]
    scaled = StandardScaler().fit_transform(variable)
    coords = PCA(n_components=2).fit_transform(scaled)
    pca = result.metadata[["sample_id", "condition"]].copy()
    pca["PC1"] = coords[:, 0]
    pca["PC2"] = coords[:, 1]
    pca["matrix_source"] = result.label

    corr = _safe_corr(variable.T)
    rows: list[dict[str, object]] = []
    conditions = dict(zip(result.metadata["sample_id"], result.metadata["condition"]))
    for sample_a, sample_b in combinations(samples, 2):
        rows.append(
            {
                "matrix_source": result.label,
                "sample_A": sample_a,
                "sample_B": sample_b,
                "condition_A": conditions[sample_a],
                "condition_B": conditions[sample_b],
                "same_condition": conditions[sample_a] == conditions[sample_b],
                "spearman_correlation": corr.loc[sample_a, sample_b],
            }
        )
    return pca, pd.DataFrame(rows)


def _condition_columns(result: MatrixResult, condition: str) -> list[str]:
    rows = result.metadata.loc[result.metadata["condition"].eq(condition), "sample_id"]
    return [sample for sample in rows.astype(str) if sample in result.matrix.columns]


def _effect_series(result: MatrixResult, numerator: str, denominator: str) -> pd.Series:
    num = _condition_columns(result, numerator)
    den = _condition_columns(result, denominator)
    if not num or not den:
        raise ValueError(f"Missing samples for contrast {numerator} vs {denominator}")
    return result.matrix[num].mean(axis=1) - result.matrix[den].mean(axis=1)


def _loo_effects(
    result: MatrixResult,
    numerator: str,
    denominator: str,
) -> Iterable[pd.Series]:
    num = _condition_columns(result, numerator)
    den = _condition_columns(result, denominator)
    n_runs = min(len(num), len(den))
    for index in range(n_runs):
        keep_num = [sample for position, sample in enumerate(num) if position != index]
        keep_den = [sample for position, sample in enumerate(den) if position != index]
        yield result.matrix[keep_num].mean(axis=1) - result.matrix[keep_den].mean(axis=1)


def _contrast_metrics(
    result: MatrixResult,
    effect_cutoff: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    summaries: list[dict[str, object]] = []
    effects: list[pd.DataFrame] = []
    for numerator, denominator in CONTRASTS:
        full = _effect_series(result, numerator, denominator)
        selected = full.abs().ge(effect_cutoff)
        loo_correlations: list[float] = []
        loo_sign_agreements: list[float] = []
        loo_jaccards: list[float] = []
        full_set = set(full.index[selected])
        for loo in _loo_effects(result, numerator, denominator):
            valid = full.notna() & loo.notna()
            rho = spearmanr(full[valid], loo[valid]).statistic if valid.sum() >= 2 else np.nan
            loo_correlations.append(float(rho))
            if selected.sum():
                loo_sign_agreements.append(
                    float((np.sign(full[selected]) == np.sign(loo[selected])).mean())
                )
            loo_set = set(loo.index[loo.abs().ge(effect_cutoff)])
            union = full_set | loo_set
            loo_jaccards.append(len(full_set & loo_set) / len(union) if union else 1.0)
        summaries.append(
            {
                "matrix_source": result.label,
                "contrast": f"{numerator} vs {denominator}",
                "numerator": numerator,
                "denominator": denominator,
                "n_stable_events": len(full),
                "n_abs_delta_psi_ge_cutoff": int(selected.sum()),
                "effect_cutoff": effect_cutoff,
                "median_abs_delta_psi": float(full.abs().median()),
                "median_LOO_effect_correlation": float(np.nanmedian(loo_correlations)),
                "median_LOO_sign_agreement": float(np.nanmedian(loo_sign_agreements)),
                "median_LOO_effect_set_jaccard": float(np.nanmedian(loo_jaccards)),
            }
        )
        effects.append(
            pd.DataFrame(
                {
                    "event_key": full.index,
                    "matrix_source": result.label,
                    "contrast": f"{numerator} vs {denominator}",
                    "delta_psi": full.to_numpy(),
                    "abs_delta_psi": full.abs().to_numpy(),
                    "passes_effect_cutoff": selected.to_numpy(),
                }
            )
        )
    return pd.DataFrame(summaries), pd.concat(effects, ignore_index=True)


def _matrix_qc_summary(result: MatrixResult) -> pd.DataFrame:
    stable = result.event_qc.loc[result.event_qc["stable_as"].astype(str).str.lower().eq("true")]
    all_events = result.event_qc
    row = {
        "matrix_source": result.label,
        "event_unit": (
            "compound canonical AS event"
            if result.label == "rMATS events"
            else "individual competing splice junction"
        ),
        "n_input_events": len(all_events),
        "n_stable_events": len(stable),
        "stable_fraction": len(stable) / len(all_events) if len(all_events) else np.nan,
        "median_replicate_repeatability_all": pd.to_numeric(
            all_events["replicate_repeatability"], errors="coerce"
        ).median(),
        "median_replicate_repeatability_stable": pd.to_numeric(
            stable["replicate_repeatability"], errors="coerce"
        ).median(),
        "median_condition_mean_psi_range_stable": pd.to_numeric(
            stable["condition_mean_psi_range"], errors="coerce"
        ).median(),
        "median_mean_total_junction_reads_stable": pd.to_numeric(
            stable["mean_total_junction_reads"], errors="coerce"
        ).median(),
        "psi_missing_fraction_stable": result.matrix.isna().to_numpy().mean(),
    }
    return pd.DataFrame([row])


def _bam_annotation_summary(bam: MatrixResult) -> tuple[pd.DataFrame, pd.DataFrame]:
    metadata = _read_tsv(bam.root / "02_matrix" / "event_metadata.tsv")
    stable_keys = set(bam.matrix.index.astype(str))
    stable = metadata.loc[metadata["event_key"].astype(str).isin(stable_keys)].copy()
    linked = pd.to_numeric(stable["rmats_event_count"], errors="coerce").fillna(0).gt(0)
    direct = stable["annotated_direct_junction"].astype(str).str.lower().eq("true")
    gene = stable["gene_symbol"].fillna("").astype(str).str.strip().ne("")
    summary = pd.DataFrame(
        [
            {
                "stable_bam_junctions": len(stable),
                "rmats_linked": int(linked.sum()),
                "bam_only": int((~linked).sum()),
                "gtf_direct_junction": int(direct.sum()),
                "gene_annotated": int(gene.sum()),
            }
        ]
    )
    return summary, stable


def _comparison_samples(metadata: pd.DataFrame, comparison: str) -> tuple[list[str], list[str]]:
    parts = comparison.split("-", 1)
    if len(parts) != 2:
        return [], []
    left, right = parts
    return (
        metadata.loc[metadata["condition"].eq(left), "sample_id"].astype(str).tolist(),
        metadata.loc[metadata["condition"].eq(right), "sample_id"].astype(str).tolist(),
    )


def _rmats_se_bam_concordance(
    bam: MatrixResult,
    fdr_cutoff: float,
    dpsi_cutoff: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    annotation_path = bam.root / "02_matrix" / "rmats_junction_annotation.tsv"
    usecols = [
        "event_key",
        "comparison_id",
        "rmats_event_type",
        "rmats_event_id",
        "rmats_gene_id",
        "rmats_gene_symbol",
        "junction_role",
        "rMATS_FDR",
        "raw_dPSI",
    ]
    chunks: list[pd.DataFrame] = []
    for chunk in pd.read_csv(
        annotation_path,
        sep="\t",
        usecols=usecols,
        chunksize=250_000,
        low_memory=False,
    ):
        chunk["rMATS_FDR"] = pd.to_numeric(chunk["rMATS_FDR"], errors="coerce")
        chunk["raw_dPSI"] = pd.to_numeric(chunk["raw_dPSI"], errors="coerce")
        keep = (
            chunk["rmats_event_type"].eq("SE")
            & chunk["rMATS_FDR"].le(fdr_cutoff)
            & chunk["raw_dPSI"].abs().ge(dpsi_cutoff)
        )
        if keep.any():
            chunks.append(chunk.loc[keep].copy())
    if not chunks:
        return pd.DataFrame(), pd.DataFrame()
    annotation = pd.concat(chunks, ignore_index=True)
    psi = _read_tsv(bam.root / "02_matrix" / "sample_level_psi_matrix.tsv").set_index(
        "event_key"
    )
    psi = psi.apply(pd.to_numeric, errors="coerce")
    rows: list[pd.DataFrame] = []
    for comparison, subset in annotation.groupby("comparison_id", sort=False):
        left, right = _comparison_samples(bam.metadata, str(comparison))
        left = [sample for sample in left if sample in psi.columns]
        right = [sample for sample in right if sample in psi.columns]
        if not left or not right:
            continue
        keys = subset["event_key"].astype(str)
        available = keys.isin(psi.index)
        subset = subset.loc[available].copy()
        junction_delta = (
            psi.loc[subset["event_key"], left].mean(axis=1).to_numpy()
            - psi.loc[subset["event_key"], right].mean(axis=1).to_numpy()
        )
        subset["bam_junction_delta_psi"] = junction_delta
        subset["role_sign"] = np.where(subset["junction_role"].eq("skipping"), -1.0, 1.0)
        subset["role_adjusted_bam_delta_psi"] = (
            subset["bam_junction_delta_psi"] * subset["role_sign"]
        )
        rows.append(subset)
    if not rows:
        return pd.DataFrame(), pd.DataFrame()
    junctions = pd.concat(rows, ignore_index=True)
    grouped = (
        junctions.groupby(
            [
                "comparison_id",
                "rmats_event_id",
                "rmats_gene_id",
                "rmats_gene_symbol",
            ],
            dropna=False,
            as_index=False,
        )
        .agg(
            rMATS_FDR=("rMATS_FDR", "min"),
            rMATS_dPSI=("raw_dPSI", "first"),
            mean_role_adjusted_bam_delta_psi=("role_adjusted_bam_delta_psi", "mean"),
            median_role_adjusted_bam_delta_psi=("role_adjusted_bam_delta_psi", "median"),
            n_mapped_junctions=("event_key", "nunique"),
        )
    )
    grouped["same_direction"] = np.sign(grouped["rMATS_dPSI"]) == np.sign(
        grouped["mean_role_adjusted_bam_delta_psi"]
    )
    raw_rmats_root: Path | None = None
    config_path = bam.root / "cosplice_config.yaml"
    if config_path.exists():
        config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
        raw_value = str(config.get("rmats_root", "")).strip()
        if raw_value:
            raw_rmats_root = Path(raw_value)
    total_significant: dict[str, int] = {}
    if raw_rmats_root is not None and raw_rmats_root.exists():
        for comparison_dir in raw_rmats_root.iterdir():
            if not comparison_dir.is_dir():
                continue
            candidates = list(comparison_dir.rglob("SE.MATS.JCEC.txt"))
            if not candidates:
                continue
            raw = pd.read_csv(
                candidates[0],
                sep="\t",
                usecols=lambda column: column in {"FDR", "IncLevelDifference"},
                low_memory=False,
            )
            raw_fdr = pd.to_numeric(raw.get("FDR"), errors="coerce")
            raw_dpsi = pd.to_numeric(raw.get("IncLevelDifference"), errors="coerce")
            total_significant[comparison_dir.name] = int(
                (raw_fdr.le(fdr_cutoff) & raw_dpsi.abs().ge(dpsi_cutoff)).sum()
            )

    summary_rows: list[dict[str, object]] = []
    for comparison, subset in grouped.groupby("comparison_id"):
        valid = subset[
            ["rMATS_dPSI", "mean_role_adjusted_bam_delta_psi"]
        ].dropna()
        rho = (
            spearmanr(
                valid["rMATS_dPSI"],
                valid["mean_role_adjusted_bam_delta_psi"],
            ).statistic
            if len(valid) >= 2
            else np.nan
        )
        same_count = int(subset["same_direction"].sum())
        direction_pvalue = binomtest(
            same_count,
            len(subset),
            p=0.5,
            alternative="greater",
        ).pvalue
        total = total_significant.get(str(comparison), np.nan)
        summary_rows.append(
            {
                "comparison": comparison,
                "n_significant_SE_events_total": total,
                "n_significant_SE_events_with_bam_mapping": len(subset),
                "bam_mapping_fraction": (
                    len(subset) / total if pd.notna(total) and total else np.nan
                ),
                "direction_concordance": subset["same_direction"].mean(),
                "direction_concordant_events": same_count,
                "direction_concordance_binomial_pvalue": direction_pvalue,
                "spearman_rMATS_dPSI_vs_BAM_delta": rho,
                "median_mapped_junctions_per_event": subset["n_mapped_junctions"].median(),
                "rMATS_FDR_cutoff": fdr_cutoff,
                "abs_rMATS_dPSI_cutoff": dpsi_cutoff,
            }
        )
    summary = pd.DataFrame(summary_rows)
    if not summary.empty:
        summary["direction_concordance_FDR"] = multipletests(
            summary["direction_concordance_binomial_pvalue"],
            method="fdr_bh",
        )[1]
    return grouped, summary


def _bam_only_effects(
    bam: MatrixResult,
    stable_metadata: pd.DataFrame,
    effect_cutoff: float,
) -> pd.DataFrame:
    linked = pd.to_numeric(
        stable_metadata["rmats_event_count"], errors="coerce"
    ).fillna(0).gt(0)
    bam_only = stable_metadata.loc[~linked].copy()
    qc = bam.event_qc.set_index("event_key")
    rows: list[pd.DataFrame] = []
    for numerator, denominator in CONTRASTS:
        effect = _effect_series(bam, numerator, denominator)
        subset = bam_only.loc[bam_only["event_key"].isin(effect.index)].copy()
        subset["contrast"] = f"{numerator} vs {denominator}"
        subset["delta_psi"] = subset["event_key"].map(effect)
        subset["abs_delta_psi"] = subset["delta_psi"].abs()
        subset["passes_effect_cutoff"] = subset["abs_delta_psi"].ge(effect_cutoff)
        subset["replicate_repeatability"] = subset["event_key"].map(
            pd.to_numeric(qc["replicate_repeatability"], errors="coerce")
        )
        rows.append(subset)
    result = pd.concat(rows, ignore_index=True)
    return result.sort_values(
        ["passes_effect_cutoff", "abs_delta_psi", "replicate_repeatability"],
        ascending=[False, False, False],
    )


def _plot_pca(pca: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    palette = dict(
        zip(
            sorted(pca["condition"].unique()),
            sns.color_palette("colorblind", n_colors=pca["condition"].nunique()),
        )
    )
    for axis, (source, subset) in zip(axes, pca.groupby("matrix_source", sort=False)):
        sns.scatterplot(
            data=subset,
            x="PC1",
            y="PC2",
            hue="condition",
            palette=palette,
            s=85,
            edgecolor="white",
            linewidth=0.8,
            ax=axis,
        )
        for row in subset.itertuples():
            axis.annotate(str(row.sample_id).replace("_rep_", "\nr"), (row.PC1, row.PC2), fontsize=7)
        axis.set_title(source)
        axis.legend(loc="best", fontsize=8)
    fig.suptitle("Sample structure from stable splicing features", fontweight="bold")
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _plot_similarity(similarity: pd.DataFrame, path: Path) -> None:
    fig, axis = plt.subplots(figsize=(8, 5), constrained_layout=True)
    plot = similarity.copy()
    plot["pair_type"] = np.where(plot["same_condition"], "within condition", "between conditions")
    sns.boxplot(
        data=plot,
        x="matrix_source",
        y="spearman_correlation",
        hue="pair_type",
        palette=["#177E89", "#D95F59"],
        ax=axis,
    )
    sns.stripplot(
        data=plot,
        x="matrix_source",
        y="spearman_correlation",
        hue="pair_type",
        dodge=True,
        palette=["#173B45", "#6F2B2B"],
        alpha=0.35,
        size=3,
        ax=axis,
        legend=False,
    )
    axis.set_title("Replicate similarity versus between-condition similarity", fontweight="bold")
    axis.set_xlabel("")
    axis.set_ylabel("Spearman correlation across stable features")
    axis.axhline(0, color="#999999", linewidth=0.8)
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _plot_contrast_robustness(summary: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
    metrics = (
        ("median_LOO_effect_correlation", "LOO effect correlation"),
        ("median_LOO_sign_agreement", "LOO sign agreement"),
        ("median_LOO_effect_set_jaccard", "LOO effect-set Jaccard"),
    )
    for axis, (metric, title) in zip(axes, metrics):
        sns.barplot(
            data=summary,
            x=metric,
            y="contrast",
            hue="matrix_source",
            palette="colorblind",
            ax=axis,
        )
        axis.set_xlim(0, 1.02)
        axis.set_title(title, fontweight="bold")
        axis.set_xlabel("")
        axis.set_ylabel("")
    fig.suptitle("Condition-effect robustness under leave-one-replicate-out", fontweight="bold")
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _plot_se_concordance(events: pd.DataFrame, path: Path) -> None:
    if events.empty:
        return
    comparisons = sorted(events["comparison_id"].unique())
    fig, axes = plt.subplots(
        1,
        len(comparisons),
        figsize=(5 * len(comparisons), 4.5),
        squeeze=False,
        constrained_layout=True,
    )
    for axis, comparison in zip(axes[0], comparisons):
        subset = events.loc[events["comparison_id"].eq(comparison)]
        axis.scatter(
            subset["rMATS_dPSI"],
            subset["mean_role_adjusted_bam_delta_psi"],
            s=14,
            alpha=0.55,
            color="#177E89",
            edgecolors="none",
        )
        axis.axhline(0, color="#888888", linewidth=0.8)
        axis.axvline(0, color="#888888", linewidth=0.8)
        axis.set_title(str(comparison), fontsize=10, fontweight="bold")
        axis.set_xlabel("rMATS dPSI")
        axis.set_ylabel("Role-adjusted BAM junction delta")
    fig.suptitle(
        "BAM-derived junction support for significant rMATS SE events",
        fontweight="bold",
    )
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _plot_bam_only_heatmap(
    bam: MatrixResult,
    bam_only: pd.DataFrame,
    path: Path,
    top_n: int = 30,
) -> None:
    chosen = (
        bam_only.loc[bam_only["passes_effect_cutoff"]]
        .drop_duplicates("event_key")
        .nlargest(top_n, "abs_delta_psi")
    )
    if chosen.empty:
        return
    values = bam.matrix.loc[chosen["event_key"]]
    labels = chosen.set_index("event_key")["gene_symbol"].fillna("")
    row_labels = [
        f"{labels.get(key, '') or 'unannotated'} | {key.split('|', 1)[-1]}"
        for key in values.index
    ]
    fig, axis = plt.subplots(figsize=(11, max(6, len(values) * 0.25)), constrained_layout=True)
    sns.heatmap(
        values,
        cmap="vlag",
        center=0.5,
        vmin=0,
        vmax=1,
        yticklabels=row_labels,
        cbar_kws={"label": "BAM junction usage"},
        ax=axis,
    )
    axis.set_title("Top stable BAM-only competing junctions", fontweight="bold")
    axis.set_xlabel("Sample")
    axis.set_ylabel("")
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _setd2_benchmark(setd2_summary: Path, output_dir: Path) -> pd.DataFrame:
    if not setd2_summary.exists():
        return pd.DataFrame()
    table = _read_tsv(setd2_summary)
    table.to_csv(output_dir / "setd2_cosplice_vs_wgcna.tsv", sep="\t", index=False)
    methods = table.loc[table["method"].isin(["katana_sparse", "wgcna_exact"])].copy()
    if methods.empty:
        return table
    metrics = [
        "silhouette_cosine",
        "neutral_graph_modularity",
        "within_between_separation",
        "weighted_condition_R2",
        "mean_bootstrap_ARI",
        "min_bootstrap_ARI",
    ]
    long = methods.melt(
        id_vars=["method"],
        value_vars=[metric for metric in metrics if metric in methods],
        var_name="metric",
        value_name="value",
    )
    fig, axis = plt.subplots(figsize=(11, 5.5), constrained_layout=True)
    sns.barplot(
        data=long,
        x="metric",
        y="value",
        hue="method",
        palette=["#177E89", "#D95F59"],
        ax=axis,
    )
    axis.set_ylim(0, 1.05)
    axis.tick_params(axis="x", rotation=30)
    axis.set_title("Same-input network benchmark", fontweight="bold")
    axis.set_xlabel("")
    axis.set_ylabel("Metric value (higher is better)")
    fig.savefig(output_dir / "setd2_cosplice_vs_wgcna.png", dpi=220, bbox_inches="tight")
    plt.close(fig)
    return table


def _plot_integrated_result(
    qc: pd.DataFrame,
    similarity: pd.DataFrame,
    contrasts: pd.DataFrame,
    se_summary: pd.DataFrame,
    setd2: pd.DataFrame,
    path: Path,
) -> None:
    sns.set_theme(style="ticks", font="Arial", context="paper")
    fig = plt.figure(figsize=(15, 10), constrained_layout=True)
    grid = fig.add_gridspec(2, 3)
    axes = [fig.add_subplot(grid[row, column]) for row in range(2) for column in range(3)]
    colors = {
        "rMATS events": "#4477AA",
        "BAM competing junctions": "#228833",
        "katana_sparse": "#228833",
        "wgcna_exact": "#CC6677",
    }

    repeatability = qc[
        ["matrix_source", "median_replicate_repeatability_stable"]
    ].copy()
    sns.barplot(
        data=repeatability,
        x="matrix_source",
        y="median_replicate_repeatability_stable",
        hue="matrix_source",
        palette=colors,
        legend=False,
        ax=axes[0],
    )
    axes[0].set_title("a  Stable-feature reproducibility", loc="left", fontweight="bold")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("Median replicate repeatability")
    axes[0].set_ylim(0, 0.8)
    axes[0].tick_params(axis="x", rotation=15)

    similarity_summary = (
        similarity.groupby(["matrix_source", "same_condition"], as_index=False)[
            "spearman_correlation"
        ]
        .median()
        .assign(
            pair_type=lambda frame: np.where(
                frame["same_condition"],
                "Within condition",
                "Between conditions",
            )
        )
    )
    sns.barplot(
        data=similarity_summary,
        x="matrix_source",
        y="spearman_correlation",
        hue="pair_type",
        palette=["#66C2A5", "#FC8D62"],
        ax=axes[1],
    )
    axes[1].set_title("b  Sample-level similarity", loc="left", fontweight="bold")
    axes[1].set_xlabel("")
    axes[1].set_ylabel("Median Spearman correlation")
    axes[1].set_ylim(0, 1)
    axes[1].tick_params(axis="x", rotation=15)
    axes[1].legend(frameon=False, fontsize=8)

    compact_contrasts = contrasts.copy()
    compact_contrasts["contrast"] = (
        compact_contrasts["contrast"]
        .str.replace("_NODP", " noDP", regex=False)
        .str.replace("_DP", " DP", regex=False)
    )
    sns.barplot(
        data=compact_contrasts,
        x="median_LOO_effect_set_jaccard",
        y="contrast",
        hue="matrix_source",
        palette=colors,
        ax=axes[2],
    )
    axes[2].set_title("c  Replicate-robust condition effects", loc="left", fontweight="bold")
    axes[2].set_xlabel("Median LOO effect-set Jaccard")
    axes[2].set_ylabel("")
    axes[2].set_xlim(0, 1)
    axes[2].legend(frameon=False, fontsize=8)

    se_plot = se_summary.copy()
    se_plot["comparison"] = (
        se_plot["comparison"]
        .str.replace("_NODP", " noDP", regex=False)
        .str.replace("_DP", " DP", regex=False)
    )
    se_long = se_plot.melt(
        id_vars=["comparison"],
        value_vars=["bam_mapping_fraction", "direction_concordance"],
        var_name="metric",
        value_name="fraction",
    )
    se_long["metric"] = se_long["metric"].map(
        {
            "bam_mapping_fraction": "Mapped to BAM junctions",
            "direction_concordance": "Same dPSI direction",
        }
    )
    sns.barplot(
        data=se_long,
        x="fraction",
        y="comparison",
        hue="metric",
        palette=["#88CCEE", "#117733"],
        ax=axes[3],
    )
    axes[3].set_title("d  rMATS significant SE support", loc="left", fontweight="bold")
    axes[3].set_xlabel("Fraction")
    axes[3].set_ylabel("")
    axes[3].set_xlim(0, 1)
    axes[3].legend(frameon=False, fontsize=8)

    sns.barplot(
        data=se_plot,
        x="spearman_rMATS_dPSI_vs_BAM_delta",
        y="comparison",
        color="#44AA99",
        ax=axes[4],
    )
    axes[4].set_title("e  Effect-size agreement", loc="left", fontweight="bold")
    axes[4].set_xlabel("Spearman: rMATS dPSI vs BAM delta")
    axes[4].set_ylabel("")
    axes[4].set_xlim(0, 0.8)

    methods = setd2.loc[setd2["method"].isin(["katana_sparse", "wgcna_exact"])].copy()
    if not methods.empty:
        methods["retained_event_fraction"] = 1 - pd.to_numeric(
            methods["grey_fraction"], errors="coerce"
        )
        network_metrics = [
            "silhouette_cosine",
            "neutral_graph_modularity",
            "mean_bootstrap_ARI",
            "min_bootstrap_ARI",
            "retained_event_fraction",
        ]
        network_long = methods.melt(
            id_vars=["method"],
            value_vars=network_metrics,
            var_name="metric",
            value_name="value",
        )
        network_long["metric"] = network_long["metric"].map(
            {
                "silhouette_cosine": "Silhouette",
                "neutral_graph_modularity": "Modularity",
                "mean_bootstrap_ARI": "Mean ARI",
                "min_bootstrap_ARI": "Min ARI",
                "retained_event_fraction": "Retained",
            }
        )
        sns.barplot(
            data=network_long,
            x="metric",
            y="value",
            hue="method",
            palette=colors,
            ax=axes[5],
        )
        axes[5].set_ylim(0, 1)
        axes[5].set_xlabel("")
        axes[5].set_ylabel("Metric value")
        axes[5].legend(
            labels=["CoSplice tight", "Exact WGCNA"],
            frameon=False,
            fontsize=8,
        )
    axes[5].set_title("f  Same-input network benchmark", loc="left", fontweight="bold")

    for axis in axes:
        sns.despine(ax=axis)
    fig.suptitle(
        "Functional validation of BAM-junction input and CoSplice network analysis",
        fontsize=17,
        fontweight="bold",
    )
    fig.text(
        0.5,
        -0.01,
        "Input-layer validation compares matched samples and condition effects. "
        "Network-method validation compares CoSplice and exact WGCNA on the same event universe. "
        "rMATS compound events and BAM junctions are different statistical units.",
        ha="center",
        fontsize=9,
        color="#4A5560",
    )
    fig.savefig(path, dpi=250, bbox_inches="tight")
    plt.close(fig)


def _write_result_summary(
    output_dir: Path,
    qc: pd.DataFrame,
    similarity: pd.DataFrame,
    contrasts: pd.DataFrame,
    bam_summary: pd.DataFrame,
    se_summary: pd.DataFrame,
    setd2: pd.DataFrame,
) -> Path:
    qc_index = qc.set_index("matrix_source")
    within = (
        similarity.loc[similarity["same_condition"]]
        .groupby("matrix_source")["spearman_correlation"]
        .median()
    )
    bam_loo = contrasts.loc[
        contrasts["matrix_source"].eq("BAM competing junctions"),
        "median_LOO_effect_set_jaccard",
    ]
    rmats_loo = contrasts.loc[
        contrasts["matrix_source"].eq("rMATS events"),
        "median_LOO_effect_set_jaccard",
    ]
    method_index = setd2.set_index("method") if not setd2.empty else pd.DataFrame()
    lines = [
        "# CoSplice functional-validation result",
        "",
        "## Main result",
        "",
        "BAM-derived competing-junction usage is a useful discovery layer, but it does not "
        "replace rMATS. The recommended workflow is BAM-junction discovery plus rMATS "
        "canonical event/FDR annotation, followed by parallel CoSplice-tight and "
        "WGCNA-diffuse network views when sample size is sufficient.",
        "",
        "## Input-layer validation",
        "",
        f"- Stable-feature repeatability: BAM "
        f"`{qc_index.loc['BAM competing junctions', 'median_replicate_repeatability_stable']:.3f}` "
        f"versus rMATS `{qc_index.loc['rMATS events', 'median_replicate_repeatability_stable']:.3f}`.",
        f"- Median within-condition correlation: BAM "
        f"`{within['BAM competing junctions']:.3f}` versus rMATS "
        f"`{within['rMATS events']:.3f}`.",
        f"- Median leave-one-replicate-out effect-set Jaccard across the declared contrasts: BAM "
        f"`{bam_loo.median():.3f}` versus rMATS `{rmats_loo.median():.3f}`.",
        f"- Significant SE BAM mapping range: "
        f"`{se_summary['bam_mapping_fraction'].min():.1%}`–"
        f"`{se_summary['bam_mapping_fraction'].max():.1%}`.",
        f"- Same-direction range among mapped significant SE events: "
        f"`{se_summary['direction_concordance'].min():.1%}`–"
        f"`{se_summary['direction_concordance'].max():.1%}`; BH-corrected direction-test "
        f"results are retained in the complete result table.",
        f"- Stable BAM-only junctions: `{int(bam_summary.loc[0, 'bam_only'])}`. These are "
        "follow-up candidates, not automatically significant differential AS events.",
        "",
        "## Network-method validation",
        "",
    ]
    if not method_index.empty:
        cosplice = method_index.loc["katana_sparse"]
        wgcna = method_index.loc["wgcna_exact"]
        lines.extend(
            [
                f"- Silhouette, CoSplice versus exact WGCNA: "
                f"(`{cosplice['silhouette_cosine']:.3f}` vs "
                f"`{wgcna['silhouette_cosine']:.3f}`); graph modularity: "
                f"(`{cosplice['neutral_graph_modularity']:.3f}` vs "
                f"`{wgcna['neutral_graph_modularity']:.3f}`).",
                f"- Worst-case bootstrap ARI, exact WGCNA versus CoSplice: "
                f"(`{wgcna['min_bootstrap_ARI']:.3f}` vs "
                f"`{cosplice['min_bootstrap_ARI']:.3f}`); retained event fractions: "
                f"(`{1-wgcna['grey_fraction']:.1%}` vs "
                f"`{1-cosplice['grey_fraction']:.1%}`).",
                "- Neither method is universally superior on this cohort. The production "
                "hybrid preserves both tight and diffuse programs.",
            ]
        )
    path = output_dir / "RESULT_SUMMARY.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _html_table(frame: pd.DataFrame, max_rows: int = 20) -> str:
    if frame.empty:
        return "<p>No result available.</p>"
    return frame.head(max_rows).to_html(index=False, border=0, classes="data")


def _write_report(
    output_dir: Path,
    qc: pd.DataFrame,
    similarity: pd.DataFrame,
    contrasts: pd.DataFrame,
    bam_summary: pd.DataFrame,
    se_summary: pd.DataFrame,
    bam_only: pd.DataFrame,
    setd2: pd.DataFrame,
    parameters: dict[str, object],
) -> Path:
    within = (
        similarity.loc[similarity["same_condition"]]
        .groupby("matrix_source")["spearman_correlation"]
        .median()
    )
    report = output_dir / "COSPLICE_INPUT_MODE_FUNCTIONAL_BENCHMARK.html"
    images = [
        "ythdc1_stable_feature_pca.png",
        "ythdc1_sample_similarity.png",
        "ythdc1_condition_effect_robustness.png",
        "ythdc1_significant_se_bam_concordance.png",
        "ythdc1_bam_only_stable_junctions.png",
        "setd2_cosplice_vs_wgcna.png",
    ]
    image_html = "".join(
        f'<figure><img src="{html.escape(image)}"><figcaption>{html.escape(image)}</figcaption></figure>'
        for image in images
        if (output_dir / image).exists()
    )
    conclusion = (
        "<li><b>Input-layer assessment:</b> compare observed BAM-junction catalogue coverage "
        "with canonical rMATS events in the recorded QC tables. rMATS retains its canonical "
        "event unit and formal differential-event statistics.</li>"
        "<li><b>Orthogonal computational support:</b> SE concordance tests use rMATS dPSI "
        "against role-adjusted junction usage re-extracted from the same BAM alignments. "
        "This is not an independent biological cohort.</li>"
        "<li><b>Network assessment:</b> the same-input benchmark compares CoSplice and exact "
        "WGCNA on the supplied cohort. Sample counts and cohort eligibility follow the "
        "recorded input metadata; these comparisons do not establish universal superiority.</li>"
        "<li><b>Recommended production mode:</b> hybrid evidence: BAM-derived junction "
        "matrix for broad network discovery, rMATS overlay for canonical AS interpretation "
        "and FDR-supported differential evidence.</li>"
    )
    text = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>CoSplice functional benchmark</title>
<style>
body{{font-family:Arial,sans-serif;margin:32px;color:#18212b;max-width:1500px}}
h1,h2{{color:#123c4a}} .note{{background:#edf6f5;border-left:5px solid #177e89;padding:14px}}
table{{border-collapse:collapse;font-size:12px;width:100%;margin:12px 0}}
th,td{{border:1px solid #d5dde3;padding:6px;text-align:left}} th{{background:#eef2f4}}
figure{{margin:24px 0}} img{{max-width:100%;height:auto;border:1px solid #d5dde3}}
code{{background:#f2f4f5;padding:2px 4px}}
</style></head><body>
<h1>CoSplice input-mode functional benchmark</h1>
<p class="note">This report separates two questions. Input-layer validation tests BAM-junction versus rMATS
input behavior on matched samples. Network-method validation tests CoSplice versus exact WGCNA network
behavior on the same event universe. Event counts are not treated as directly interchangeable:
rMATS uses compound AS events; BAM mode uses individual competing splice junctions.</p>
<h2>Executive conclusion</h2><ul>{conclusion}</ul>
<h2>Matrix and stable-event QC</h2>{_html_table(qc)}
<p>Median within-condition sample correlations: {html.escape(within.to_dict().__str__())}</p>
<h2>BAM annotation gain</h2>{_html_table(bam_summary)}
<h2>Leave-one-replicate-out condition effects</h2>{_html_table(contrasts, 50)}
<h2>Significant rMATS SE events supported by BAM-derived junction usage</h2>
{_html_table(se_summary, 50)}
<p>Only SE events are used for the formal sign test because inclusion and skipping junction
roles have an unambiguous sign transformation. Other AS types remain available as annotations
but are not forced into an invalid common sign rule.</p>
<h2>Stable BAM-only follow-up</h2>
<p>{int(bam_only['event_key'].nunique()) if not bam_only.empty else 0} stable BAM-only
junctions were evaluated. These are discovery candidates, not automatically significant
differential AS events.</p>
{_html_table(bam_only.loc[bam_only.get('passes_effect_cutoff', False)], 20)}
<h2>CoSplice versus exact WGCNA</h2>{_html_table(setd2, 20)}
<p>Compare observed bootstrap stability, event retention, silhouette, graph modularity and
runtime in the supplied result table. Their relative performance is dataset-specific;
one cohort does not justify a universal superiority claim.</p>
<h2>Figures</h2>{image_html}
<h2>Reproducibility</h2><pre>{html.escape(json.dumps(parameters, indent=2))}</pre>
</body></html>"""
    report.write_text(text, encoding="utf-8")
    return report


def run_benchmark(
    bam_root: Path,
    rmats_root: Path,
    output_dir: Path,
    setd2_summary: Path | None = None,
    effect_cutoff: float = 0.10,
    rmats_fdr_cutoff: float = 0.05,
    rmats_dpsi_cutoff: float = 0.10,
) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    bam = _load_matrix_result("BAM competing junctions", bam_root)
    rmats = _load_matrix_result("rMATS events", rmats_root)
    if list(bam.metadata["sample_id"]) != list(rmats.metadata["sample_id"]):
        raise ValueError("BAM and rMATS benchmark inputs do not contain the same ordered samples.")

    qc = pd.concat([_matrix_qc_summary(rmats), _matrix_qc_summary(bam)], ignore_index=True)
    pca_parts: list[pd.DataFrame] = []
    similarity_parts: list[pd.DataFrame] = []
    contrast_parts: list[pd.DataFrame] = []
    effect_parts: list[pd.DataFrame] = []
    for result in (rmats, bam):
        pca, similarity = _sample_similarity(result)
        contrast, effect = _contrast_metrics(result, effect_cutoff)
        pca_parts.append(pca)
        similarity_parts.append(similarity)
        contrast_parts.append(contrast)
        effect_parts.append(effect)
    pca = pd.concat(pca_parts, ignore_index=True)
    similarity = pd.concat(similarity_parts, ignore_index=True)
    contrasts = pd.concat(contrast_parts, ignore_index=True)
    effects = pd.concat(effect_parts, ignore_index=True)

    bam_summary, stable_metadata = _bam_annotation_summary(bam)
    se_events, se_summary = _rmats_se_bam_concordance(
        bam,
        rmats_fdr_cutoff,
        rmats_dpsi_cutoff,
    )
    bam_only = _bam_only_effects(bam, stable_metadata, effect_cutoff)
    setd2 = (
        _setd2_benchmark(setd2_summary, output_dir)
        if setd2_summary is not None
        else pd.DataFrame()
    )

    tables = {
        "matrix_qc_summary.tsv": qc,
        "sample_pca.tsv": pca,
        "sample_pair_similarity.tsv": similarity,
        "condition_effect_robustness.tsv": contrasts,
        "condition_event_effects.tsv": effects,
        "bam_stable_annotation_summary.tsv": bam_summary,
        "significant_rmats_se_bam_support_events.tsv": se_events,
        "significant_rmats_se_bam_support_summary.tsv": se_summary,
        "bam_only_stable_junction_effects.tsv": bam_only,
    }
    for name, frame in tables.items():
        frame.to_csv(output_dir / name, sep="\t", index=False)

    _plot_pca(pca, output_dir / "ythdc1_stable_feature_pca.png")
    _plot_similarity(similarity, output_dir / "ythdc1_sample_similarity.png")
    _plot_contrast_robustness(
        contrasts,
        output_dir / "ythdc1_condition_effect_robustness.png",
    )
    _plot_se_concordance(
        se_events,
        output_dir / "ythdc1_significant_se_bam_concordance.png",
    )
    _plot_bam_only_heatmap(
        bam,
        bam_only,
        output_dir / "ythdc1_bam_only_stable_junctions.png",
    )
    _plot_integrated_result(
        qc,
        similarity,
        contrasts,
        se_summary,
        setd2,
        output_dir / "cosplice_functional_validation_summary.png",
    )
    _write_result_summary(
        output_dir,
        qc,
        similarity,
        contrasts,
        bam_summary,
        se_summary,
        setd2,
    )

    parameters = {
        "bam_root": str(bam_root.resolve()),
        "rmats_root": str(rmats_root.resolve()),
        "setd2_method_summary": str(setd2_summary.resolve()) if setd2_summary else None,
        "effect_cutoff": effect_cutoff,
        "rmats_fdr_cutoff": rmats_fdr_cutoff,
        "rmats_abs_dpsi_cutoff": rmats_dpsi_cutoff,
        "contrasts": [f"{left} vs {right}" for left, right in CONTRASTS],
        "event_unit_warning": (
            "rMATS compound AS events and BAM individual competing junctions are not "
            "directly interchangeable event counts."
        ),
    }
    (output_dir / "parameters.json").write_text(
        json.dumps(parameters, indent=2),
        encoding="utf-8",
    )
    return _write_report(
        output_dir,
        qc,
        similarity,
        contrasts,
        bam_summary,
        se_summary,
        bam_only,
        setd2,
        parameters,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Benchmark CoSplice BAM-junction and rMATS event inputs."
    )
    parser.add_argument("--bam-root", type=Path, required=True)
    parser.add_argument("--rmats-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--setd2-method-summary", type=Path)
    parser.add_argument("--effect-cutoff", type=float, default=0.10)
    parser.add_argument("--rmats-fdr-cutoff", type=float, default=0.05)
    parser.add_argument("--rmats-dpsi-cutoff", type=float, default=0.10)
    args = parser.parse_args()
    report = run_benchmark(
        args.bam_root,
        args.rmats_root,
        args.output_dir,
        args.setd2_method_summary,
        args.effect_cutoff,
        args.rmats_fdr_cutoff,
        args.rmats_dpsi_cutoff,
    )
    print(report)


if __name__ == "__main__":
    main()
