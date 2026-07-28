from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

from katana_cosplice.config import FactorContrastConfig
from katana_cosplice.utils import bh_adjust


@dataclass
class FactorContrastResult:
    event_contrasts: pd.DataFrame
    module_contrasts: pd.DataFrame
    phenocopy_summary: pd.DataFrame
    evidence_summary: pd.DataFrame


def _aligned_metadata(matrix: pd.DataFrame, samples: pd.DataFrame) -> pd.DataFrame:
    metadata = samples.copy()
    if "sample_id" in metadata.columns:
        metadata = metadata.set_index("sample_id")
    metadata.index = metadata.index.astype(str)
    missing = matrix.columns.difference(metadata.index)
    if len(missing):
        raise ValueError(f"Sample metadata is missing matrix samples: {missing[:5].tolist()}")
    return metadata.loc[matrix.columns].copy()


def _context_definitions(metadata: pd.DataFrame, factors: list[str]) -> list[tuple[str, dict[str, str]]]:
    if not factors:
        return [("all", {})]
    levels: list[list[str]] = []
    for factor in factors:
        if factor not in metadata.columns:
            raise ValueError(f"Configured context factor is missing from sample metadata: {factor}")
        values = [str(value) for value in metadata[factor].dropna().unique()]
        levels.append(values)
    contexts: list[tuple[str, dict[str, str]]] = []
    for combination in product(*levels):
        filters = dict(zip(factors, combination, strict=True))
        label = " | ".join(f"{key}={value}" for key, value in filters.items())
        contexts.append((label, filters))
    return contexts


def _sample_ids(metadata: pd.DataFrame, filters: dict[str, str]) -> list[str]:
    mask = pd.Series(True, index=metadata.index)
    for column, value in filters.items():
        if column not in metadata.columns:
            raise ValueError(f"Configured factor is missing from sample metadata: {column}")
        mask &= metadata[column].astype(str).eq(str(value))
    return metadata.index[mask].tolist()


def _welch_contrast(
    matrix: pd.DataFrame,
    target_ids: list[str],
    reference_ids: list[str],
    min_samples: int,
) -> pd.DataFrame:
    target = matrix.loc[:, target_ids].apply(pd.to_numeric, errors="coerce")
    reference = matrix.loc[:, reference_ids].apply(pd.to_numeric, errors="coerce")
    n_target = target.notna().sum(axis=1).astype(float)
    n_reference = reference.notna().sum(axis=1).astype(float)
    mean_target = target.mean(axis=1)
    mean_reference = reference.mean(axis=1)
    var_target = target.var(axis=1, ddof=1)
    var_reference = reference.var(axis=1, ddof=1)
    effect = mean_target - mean_reference
    term_target = var_target / n_target
    term_reference = var_reference / n_reference
    se = np.sqrt(term_target + term_reference)
    denominator = (term_target.pow(2) / (n_target - 1.0)) + (
        term_reference.pow(2) / (n_reference - 1.0)
    )
    dof = (term_target + term_reference).pow(2) / denominator
    statistic = effect / se
    pvalue = pd.Series(2.0 * stats.t.sf(np.abs(statistic), dof), index=matrix.index)
    valid = n_target.ge(min_samples) & n_reference.ge(min_samples) & se.gt(0) & dof.gt(0)
    pvalue = pvalue.where(valid, np.nan)
    result = pd.DataFrame(
        {
            "mean_target": mean_target,
            "mean_reference": mean_reference,
            "effect": effect,
            "standard_error": se,
            "statistic": statistic.where(valid, np.nan),
            "degrees_freedom": dof.where(valid, np.nan),
            "pvalue": pvalue,
            "n_target": n_target.astype(int),
            "n_reference": n_reference.astype(int),
        }
    )
    result["FDR"] = bh_adjust(result["pvalue"])
    return result


def _difference_in_differences(
    matrix: pd.DataFrame,
    target_treated: list[str],
    target_baseline: list[str],
    reference_treated: list[str],
    reference_baseline: list[str],
    min_samples: int,
) -> pd.DataFrame:
    groups = []
    for ids in (target_treated, target_baseline, reference_treated, reference_baseline):
        group = matrix.loc[:, ids].apply(pd.to_numeric, errors="coerce")
        groups.append(
            (
                group.mean(axis=1),
                group.var(axis=1, ddof=1),
                group.notna().sum(axis=1).astype(float),
            )
        )
    (mt, vt, nt), (mb, vb, nb), (rt, vr, nr), (rb, vbr, nbr) = groups
    target_effect = mt - mb
    reference_effect = rt - rb
    interaction = target_effect - reference_effect
    se = np.sqrt(vt / nt + vb / nb + vr / nr + vbr / nbr)
    z = interaction / se
    valid = nt.ge(min_samples) & nb.ge(min_samples) & nr.ge(min_samples) & nbr.ge(min_samples) & se.gt(0)
    pvalue = pd.Series(2.0 * stats.norm.sf(np.abs(z)), index=matrix.index).where(valid, np.nan)
    result = pd.DataFrame(
        {
            "target_perturbation_effect": target_effect,
            "reference_perturbation_effect": reference_effect,
            "effect": interaction,
            "standard_error": se,
            "statistic": z.where(valid, np.nan),
            "pvalue": pvalue,
            "n_target_treated": nt.astype(int),
            "n_target_baseline": nb.astype(int),
            "n_reference_treated": nr.astype(int),
            "n_reference_baseline": nbr.astype(int),
        }
    )
    same_direction = target_effect.mul(reference_effect).gt(0)
    attenuation = 1.0 - target_effect.abs().div(reference_effect.abs().replace(0, np.nan))
    result["same_perturbation_direction"] = same_direction
    result["perturbation_attenuation_fraction"] = attenuation.where(same_direction)
    result["interaction_pattern"] = np.select(
        [
            target_effect.mul(reference_effect).lt(0),
            same_direction & attenuation.gt(0),
            same_direction & attenuation.lt(0),
        ],
        ["opposite_response", "attenuated_same_direction", "enhanced_same_direction"],
        default="undetermined",
    )
    result["FDR"] = bh_adjust(result["pvalue"])
    return result


def _add_identity(table: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    result = table.reset_index(names="event_key")
    available = [column for column in ("event_key", "gene_id", "gene_symbol", "event_type") if column in metadata]
    if len(available) > 1:
        result = result.merge(metadata.loc[:, available].drop_duplicates("event_key"), on="event_key", how="left")
    return result


def analyze_factor_contrasts(
    psi: pd.DataFrame,
    program_scores: pd.DataFrame,
    samples: pd.DataFrame,
    event_metadata: pd.DataFrame,
    settings: FactorContrastConfig,
) -> FactorContrastResult:
    if not settings.enabled:
        raise ValueError("factor_contrasts is not enabled in the CoSplice config")
    metadata = _aligned_metadata(psi, samples)
    score_matrix = program_scores.T.copy()
    score_matrix.columns = score_matrix.columns.astype(str)
    missing_scores = psi.columns.difference(score_matrix.columns)
    if len(missing_scores):
        raise ValueError(f"Program scores are missing samples: {missing_scores[:5].tolist()}")
    score_matrix = score_matrix.loc[:, psi.columns]
    contexts = _context_definitions(metadata, settings.context_factors)
    event_tables: list[pd.DataFrame] = []
    module_tables: list[pd.DataFrame] = []

    for context_label, context_filters in contexts:
        baseline = {**context_filters, settings.perturbation_factor: settings.baseline_level}
        target_baseline = _sample_ids(metadata, {**baseline, settings.primary_factor: settings.target_level})
        reference_baseline = _sample_ids(metadata, {**baseline, settings.primary_factor: settings.reference_level})
        if min(len(target_baseline), len(reference_baseline)) < settings.min_samples_per_group:
            continue
        for matrix, sink, item_name in (
            (psi, event_tables, "event"),
            (score_matrix, module_tables, "module"),
        ):
            direct = _welch_contrast(matrix, target_baseline, reference_baseline, settings.min_samples_per_group)
            direct["item_type"] = item_name
            direct["contrast_type"] = "primary_baseline"
            direct["context"] = context_label
            direct["perturbation"] = settings.baseline_level
            direct["target_group"] = settings.target_level
            direct["reference_group"] = settings.reference_level
            sink.append(direct.reset_index(names="item_id"))

        for perturbation in settings.perturbation_levels:
            treated = {**context_filters, settings.perturbation_factor: perturbation}
            reference_treated = _sample_ids(
                metadata, {**treated, settings.primary_factor: settings.reference_level}
            )
            target_treated = _sample_ids(metadata, {**treated, settings.primary_factor: settings.target_level})
            if min(len(reference_treated), len(target_treated)) < settings.min_samples_per_group:
                continue
            for matrix, sink, item_name in (
                (psi, event_tables, "event"),
                (score_matrix, module_tables, "module"),
            ):
                reference_effect = _welch_contrast(
                    matrix, reference_treated, reference_baseline, settings.min_samples_per_group
                )
                reference_effect["item_type"] = item_name
                reference_effect["contrast_type"] = "perturbation_reference"
                reference_effect["context"] = context_label
                reference_effect["perturbation"] = perturbation
                reference_effect["target_group"] = settings.reference_level
                reference_effect["reference_group"] = settings.baseline_level
                sink.append(reference_effect.reset_index(names="item_id"))

                target_effect = _welch_contrast(
                    matrix, target_treated, target_baseline, settings.min_samples_per_group
                )
                target_effect["item_type"] = item_name
                target_effect["contrast_type"] = "perturbation_target"
                target_effect["context"] = context_label
                target_effect["perturbation"] = perturbation
                target_effect["target_group"] = settings.target_level
                target_effect["reference_group"] = settings.baseline_level
                sink.append(target_effect.reset_index(names="item_id"))

                interaction = _difference_in_differences(
                    matrix,
                    target_treated,
                    target_baseline,
                    reference_treated,
                    reference_baseline,
                    settings.min_samples_per_group,
                )
                interaction["item_type"] = item_name
                interaction["contrast_type"] = "perturbation_occlusion_interaction"
                interaction["context"] = context_label
                interaction["perturbation"] = perturbation
                interaction["target_group"] = settings.target_level
                interaction["reference_group"] = settings.reference_level
                sink.append(interaction.reset_index(names="item_id"))

    events = pd.concat(event_tables, ignore_index=True) if event_tables else pd.DataFrame()
    modules = pd.concat(module_tables, ignore_index=True) if module_tables else pd.DataFrame()
    if events.empty:
        raise ValueError("No configured factor contrast had enough samples to test")
    events["significant"] = (
        pd.to_numeric(events["FDR"], errors="coerce").le(settings.fdr_cutoff)
        & pd.to_numeric(events["effect"], errors="coerce").abs().ge(settings.min_abs_delta_psi)
    )
    modules["significant"] = pd.to_numeric(modules["FDR"], errors="coerce").le(settings.fdr_cutoff)
    identity = event_metadata.copy()
    if "event_key" not in identity.columns:
        identity = identity.reset_index(names="event_key")
    identity = identity.drop_duplicates("event_key")
    events = events.rename(columns={"item_id": "event_key"})
    events = events.merge(
        identity.loc[:, [column for column in ("event_key", "gene_id", "gene_symbol", "event_type") if column in identity]],
        on="event_key",
        how="left",
    )
    modules = modules.rename(columns={"item_id": "module"})

    summaries: list[dict[str, object]] = []
    for context_label, _ in contexts:
        primary = events.loc[
            events["contrast_type"].eq("primary_baseline") & events["context"].eq(context_label),
            ["event_key", "effect", "significant"],
        ].rename(columns={"effect": "primary_effect", "significant": "primary_significant"})
        for perturbation in settings.perturbation_levels:
            perturb = events.loc[
                events["contrast_type"].eq("perturbation_reference")
                & events["context"].eq(context_label)
                & events["perturbation"].eq(perturbation),
                ["event_key", "effect", "significant"],
            ].rename(columns={"effect": "perturbation_effect", "significant": "perturbation_significant"})
            target_perturb = events.loc[
                events["contrast_type"].eq("perturbation_target")
                & events["context"].eq(context_label)
                & events["perturbation"].eq(perturbation),
                ["event_key", "effect"],
            ].rename(columns={"effect": "target_perturbation_effect"})
            merged = primary.merge(perturb, on="event_key", how="inner").merge(
                target_perturb, on="event_key", how="left"
            )
            valid = merged[["primary_effect", "perturbation_effect"]].dropna()
            shared = merged.loc[merged["primary_significant"] & merged["perturbation_significant"]].copy()
            primary_count = int(merged["primary_significant"].sum())
            perturbation_count = int(merged["perturbation_significant"].sum())
            overlap_count = len(shared)
            same_count = int((shared["primary_effect"] * shared["perturbation_effect"] > 0).sum())
            universe = len(merged)
            enrichment_p = (
                float(stats.hypergeom.sf(overlap_count - 1, universe, perturbation_count, primary_count))
                if universe and primary_count and perturbation_count
                else np.nan
            )
            reference_changed = merged.loc[merged["perturbation_significant"]].copy()
            same_perturbation_direction = (
                reference_changed["target_perturbation_effect"]
                * reference_changed["perturbation_effect"]
            ).gt(0)
            attenuation = 1.0 - (
                reference_changed["target_perturbation_effect"].abs()
                / reference_changed["perturbation_effect"].abs().replace(0, np.nan)
            )
            summaries.append(
                {
                    "context": context_label,
                    "perturbation": perturbation,
                    "n_tested_events": universe,
                    "primary_significant_events": primary_count,
                    "perturbation_significant_events": perturbation_count,
                    "shared_significant_events": overlap_count,
                    "primary_signature_contained_pct": 100.0 * overlap_count / primary_count if primary_count else 0.0,
                    "shared_same_direction_events": same_count,
                    "shared_same_direction_pct": 100.0 * same_count / overlap_count if overlap_count else np.nan,
                    "broad_signature_spearman": valid["primary_effect"].corr(valid["perturbation_effect"], method="spearman"),
                    "broad_signature_pearson": valid["primary_effect"].corr(valid["perturbation_effect"], method="pearson"),
                    "overlap_enrichment_pvalue": enrichment_p,
                    "same_direction_perturbation_pct": (
                        100.0 * float(same_perturbation_direction.mean()) if len(reference_changed) else np.nan
                    ),
                    "opposite_response_pct": (
                        100.0 * float((~same_perturbation_direction).mean()) if len(reference_changed) else np.nan
                    ),
                    "median_same_direction_attenuation_in_target": float(
                        attenuation.where(same_perturbation_direction)
                        .replace([np.inf, -np.inf], np.nan)
                        .median()
                    ),
                }
            )
    phenocopy = pd.DataFrame(summaries)
    phenocopy["overlap_enrichment_FDR"] = bh_adjust(phenocopy["overlap_enrichment_pvalue"])
    evidence = phenocopy.copy()
    evidence["evidence_scope"] = "transcriptome-wide splicing signature; prioritization, not causal proof"
    evidence["interpretation"] = np.select(
        [
            evidence["broad_signature_spearman"].ge(0.30)
            & evidence["shared_same_direction_pct"].ge(80)
            & evidence["overlap_enrichment_FDR"].le(settings.fdr_cutoff),
            evidence["broad_signature_spearman"].ge(0.30),
            evidence["broad_signature_spearman"].gt(0),
        ],
        [
            "partial concordant phenocopy support",
            "broad association only; no significant-overlap support",
            "weak positive association; insufficient phenocopy support",
        ],
        default="no coherent phenocopy support",
    )
    evidence["functional_state_conclusion"] = (
        "insufficient for LOF/GOF classification without orthogonal activity and rescue evidence"
    )
    return FactorContrastResult(events, modules, phenocopy, evidence)


def render_factor_contrast_plots(
    result: FactorContrastResult,
    output_dir: Path,
    settings: FactorContrastConfig | None = None,
) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    summary = result.phenocopy_summary.copy()
    if not summary.empty:
        labels = (summary["context"] + "\n" + summary["perturbation"]).tolist()
        fig, axes = plt.subplots(1, 3, figsize=(16, 5.4))
        axes[0].bar(labels, summary["broad_signature_spearman"], color="#217A7A")
        axes[0].axhline(0, color="#263238", linewidth=0.8)
        axes[0].set_title("Broad stable-AS signature similarity")
        axes[0].set_ylabel("Spearman correlation")
        for index, value in enumerate(summary["broad_signature_spearman"]):
            axes[0].text(index, float(value) + 0.01, f"{float(value):.2f}", ha="center", fontsize=8)

        x = np.arange(len(labels))
        width = 0.38
        primary_counts = pd.to_numeric(summary["primary_significant_events"], errors="coerce").fillna(0)
        perturbation_counts = pd.to_numeric(summary["perturbation_significant_events"], errors="coerce").fillna(0)
        axes[1].bar(x - width / 2, primary_counts, width, label="Primary target-reference", color="#3A6EA5")
        axes[1].bar(x + width / 2, perturbation_counts, width, label="Reference perturbation", color="#E38B2C")
        axes[1].set_xticks(x, labels)
        axes[1].set_yscale("symlog", linthresh=1)
        axes[1].set_title("Significant stable-AS event counts")
        axes[1].set_ylabel("Events (symlog scale)")
        axes[1].legend(frameon=False, fontsize=8)
        for positions, values in ((x - width / 2, primary_counts), (x + width / 2, perturbation_counts)):
            for position, value in zip(positions, values, strict=True):
                axes[1].text(position, max(float(value), 0.08), f"{int(value)}", ha="center", va="bottom", fontsize=8)

        contained = pd.to_numeric(summary["primary_signature_contained_pct"], errors="coerce").fillna(0)
        axes[2].bar(labels, contained, color="#D1495B")
        axes[2].set_title("Primary significant-event containment")
        axes[2].set_ylabel("Primary events also significant after perturbation (%)")
        axes[2].set_ylim(0, 100)
        for index, row in summary.reset_index(drop=True).iterrows():
            shared = int(row["shared_significant_events"])
            same_direction = row["shared_same_direction_pct"]
            annotation = f"shared n={shared}"
            if shared and pd.notna(same_direction):
                annotation += f"\n{float(same_direction):.0f}% same direction"
            axes[2].text(index, max(float(contained.iloc[index]), 2), annotation, ha="center", va="bottom", fontsize=8)
        for ax in axes:
            ax.tick_params(axis="x", rotation=35)
            ax.grid(axis="y", alpha=0.2)
        if settings is None:
            title = "Configured factor phenocopy evidence"
        else:
            title = (
                f"Configured factor evidence: {settings.target_level} vs {settings.reference_level} "
                f"relative to {', '.join(settings.perturbation_levels)}"
            )
        fig.suptitle(f"{title} (descriptive, not LOF proof)", fontweight="bold")
        fig.tight_layout()
        path = output_dir / "factor_phenocopy_dashboard.png"
        fig.savefig(path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        paths.append(path)

    events = result.event_contrasts
    pairs: list[tuple[str, str, pd.DataFrame]] = []
    for row in summary.itertuples(index=False):
        primary = events.loc[
            events["contrast_type"].eq("primary_baseline") & events["context"].eq(row.context),
            ["event_key", "effect"],
        ].rename(columns={"effect": "primary_effect"})
        perturb = events.loc[
            events["contrast_type"].eq("perturbation_reference")
            & events["context"].eq(row.context)
            & events["perturbation"].eq(row.perturbation),
            ["event_key", "effect"],
        ].rename(columns={"effect": "perturbation_effect"})
        pairs.append((row.context, row.perturbation, primary.merge(perturb, on="event_key").dropna()))
    if pairs:
        ncols = 2
        nrows = int(np.ceil(len(pairs) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(12, 5.4 * nrows), squeeze=False)
        for ax, (context, perturbation, table) in zip(axes.flat, pairs, strict=False):
            ax.hexbin(table["primary_effect"], table["perturbation_effect"], gridsize=45, mincnt=1, cmap="viridis")
            limit = max(0.1, float(np.nanmax(np.abs(table[["primary_effect", "perturbation_effect"]].to_numpy()))))
            ax.plot([-limit, limit], [-limit, limit], color="#D1495B", linewidth=1.0, linestyle="--")
            ax.axhline(0, color="#263238", linewidth=0.7)
            ax.axvline(0, color="#263238", linewidth=0.7)
            rho = table["primary_effect"].corr(table["perturbation_effect"], method="spearman")
            ax.set_title(f"{context}: {perturbation} (rho={rho:.2f})")
            ax.set_xlabel("Primary target-reference delta PSI")
            ax.set_ylabel("Reference perturbation-baseline delta PSI")
        for ax in axes.flat[len(pairs):]:
            ax.axis("off")
        comparison_label = (
            f"{settings.target_level} vs {settings.reference_level}"
            if settings is not None
            else "primary target-reference"
        )
        fig.suptitle(
            f"Stable-AS signature comparison across all matched events: {comparison_label}",
            fontweight="bold",
        )
        fig.tight_layout()
        path = output_dir / "factor_phenocopy_signature_scatter.png"
        fig.savefig(path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        paths.append(path)

    modules = result.module_contrasts
    interactions = modules.loc[modules["contrast_type"].eq("perturbation_occlusion_interaction")].copy()
    if not interactions.empty:
        interactions["comparison"] = (
            interactions["context"].str.replace(r"^[^=]+=", "", regex=True)
            + " | "
            + interactions["perturbation"]
        )
        matrix = interactions.pivot_table(index="module", columns="comparison", values="effect", aggfunc="first")
        fdr = interactions.pivot_table(index="module", columns="comparison", values="FDR", aggfunc="first")
        annotations = matrix.copy().astype(object)
        for module in matrix.index:
            for comparison in matrix.columns:
                value = matrix.loc[module, comparison]
                suffix = "*" if pd.notna(fdr.loc[module, comparison]) and fdr.loc[module, comparison] <= 0.05 else ""
                annotations.loc[module, comparison] = f"{value:.2f}{suffix}" if pd.notna(value) else ""
        fig, ax = plt.subplots(figsize=(max(8.0, 2.0 * matrix.shape[1]), max(5.0, 0.55 * matrix.shape[0] + 2)))
        sns.heatmap(matrix, cmap="vlag", center=0, annot=annotations, fmt="", ax=ax, cbar_kws={"label": "Difference-in-differences"})
        ax.set_title(
            "Perturbation x primary-factor interaction on CoSplice program scores\n"
            "* FDR <= 0.05; inspect the table for attenuation versus opposite response",
            fontweight="bold",
            fontsize=13,
            pad=12,
        )
        ax.set_xlabel("Context and perturbation")
        ax.set_ylabel("CoSplice module")
        fig.tight_layout()
        path = output_dir / "module_perturbation_occlusion_heatmap.png"
        fig.savefig(path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        paths.append(path)
    return paths
