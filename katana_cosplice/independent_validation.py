from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, ttest_ind, wilcoxon

from katana_cosplice.utils import bh_adjust


@dataclass
class ProjectionModel:
    events: list[str]
    loadings: np.ndarray


def fit_module_projection_models(
    psi: pd.DataFrame,
    assignments: pd.DataFrame,
    *,
    module_column: str = "consensus_module",
) -> dict[str, ProjectionModel]:
    """Fit locked PC1 loadings for non-grey discovery modules."""
    models: dict[str, ProjectionModel] = {}
    indexed = assignments.copy()
    indexed["event_key"] = indexed["event_key"].astype(str)
    for module, group in indexed.loc[
        indexed[module_column].astype(str).ne("grey")
    ].groupby(module_column, observed=True):
        events = psi.index.intersection(group["event_key"])
        if len(events) < 3:
            continue
        values = psi.loc[events].T.apply(pd.to_numeric, errors="coerce")
        standardized = (
            (values - values.mean(axis=0))
            / values.std(axis=0).replace(0, np.nan)
        ).fillna(0.0)
        _left, _singular, right = np.linalg.svd(
            standardized.to_numpy(dtype=float),
            full_matrices=False,
        )
        loadings = right[0]
        score = standardized.to_numpy(dtype=float) @ loadings
        mean_profile = standardized.mean(axis=1).to_numpy(dtype=float)
        if np.std(score) > 0 and np.std(mean_profile) > 0:
            if np.corrcoef(score, mean_profile)[0, 1] < 0:
                loadings = -loadings
        models[str(module)] = ProjectionModel(events.tolist(), loadings)
    return models


def project_module_scores(
    psi: pd.DataFrame,
    models: dict[str, ProjectionModel],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Project locked discovery loadings onto exact events in a target cohort."""
    scores: dict[str, pd.Series] = {}
    audit: list[dict[str, object]] = []
    for module, model in models.items():
        events = pd.Index(model.events).intersection(psi.index)
        positions = [model.events.index(event) for event in events]
        if len(events) < 3:
            audit.append(
                {
                    "module": module,
                    "reference_events": len(model.events),
                    "mapped_events": len(events),
                    "mapped_fraction": len(events) / max(len(model.events), 1),
                    "status": "insufficient_exact_events",
                }
            )
            continue
        loadings = np.asarray(model.loadings, dtype=float)[positions]
        norm = float(np.linalg.norm(loadings))
        if norm <= 0:
            continue
        loadings /= norm
        values = psi.loc[events].T.apply(pd.to_numeric, errors="coerce")
        standardized = (
            (values - values.mean(axis=0))
            / values.std(axis=0).replace(0, np.nan)
        ).fillna(0.0)
        scores[module] = pd.Series(
            standardized.to_numpy(dtype=float) @ loadings,
            index=values.index,
        )
        audit.append(
            {
                "module": module,
                "reference_events": len(model.events),
                "mapped_events": len(events),
                "mapped_fraction": len(events) / max(len(model.events), 1),
                "status": "projected",
            }
        )
    return pd.DataFrame(scores), pd.DataFrame(audit)


def _hedges_g(left: pd.Series, right: pd.Series) -> float:
    pooled_df = len(left) + len(right) - 2
    if pooled_df <= 0:
        return np.nan
    pooled = np.sqrt(
        (
            (len(left) - 1) * left.var(ddof=1)
            + (len(right) - 1) * right.var(ddof=1)
        )
        / pooled_df
    )
    if not np.isfinite(pooled) or pooled <= 0:
        return np.nan
    correction = 1.0 - 3.0 / (4.0 * pooled_df - 1.0)
    return float(correction * (left.mean() - right.mean()) / pooled)


def _paired_effect(left: pd.Series, right: pd.Series) -> float:
    difference = left.to_numpy(dtype=float) - right.to_numpy(dtype=float)
    scale = float(np.std(difference, ddof=1))
    return float(np.mean(difference) / scale) if scale > 0 else np.nan


def tumor_normal_module_association(
    scores: pd.DataFrame,
    sample_metadata: pd.DataFrame,
    *,
    tumor_label: str = "tumor",
    normal_label: str = "normal",
) -> pd.DataFrame:
    """Test projected module scores using all samples and matched patient pairs."""
    metadata = sample_metadata.copy()
    metadata["sample_id"] = metadata["sample_id"].astype(str)
    metadata = metadata.set_index("sample_id").reindex(scores.index)
    phenotype = metadata["phenotype"].astype(str)
    patient = pd.Series(
        scores.index.astype(str).str.replace(
            r"_Norm$",
            "",
            regex=True,
            case=False,
        ),
        index=scores.index,
    )
    rows: list[dict[str, object]] = []
    for module in scores:
        values = pd.to_numeric(scores[module], errors="coerce")
        tumor = values.loc[phenotype.eq(tumor_label)].dropna()
        normal = values.loc[phenotype.eq(normal_label)].dropna()
        if len(tumor) < 3 or len(normal) < 3:
            continue
        welch = ttest_ind(tumor, normal, equal_var=False, nan_policy="omit")
        mann_whitney = mannwhitneyu(
            tumor,
            normal,
            alternative="two-sided",
        )
        tumor_by_patient = pd.Series(
            tumor.to_numpy(),
            index=patient.reindex(tumor.index).to_numpy(),
        )
        normal_by_patient = pd.Series(
            normal.to_numpy(),
            index=patient.reindex(normal.index).to_numpy(),
        )
        paired_ids = tumor_by_patient.index.intersection(normal_by_patient.index)
        paired_tumor = tumor_by_patient.reindex(paired_ids)
        paired_normal = normal_by_patient.reindex(paired_ids)
        paired_test = (
            wilcoxon(paired_tumor, paired_normal, zero_method="wilcox")
            if len(paired_ids) >= 3
            and not np.allclose(paired_tumor, paired_normal)
            else None
        )
        rows.append(
            {
                "module": module,
                "contrast": "tumor_vs_adjacent_normal",
                "n_tumor": len(tumor),
                "n_normal": len(normal),
                "tumor_mean_score": float(tumor.mean()),
                "normal_mean_score": float(normal.mean()),
                "mean_score_difference": float(tumor.mean() - normal.mean()),
                "hedges_g": _hedges_g(tumor, normal),
                "welch_t": float(welch.statistic),
                "welch_pvalue": float(welch.pvalue),
                "mann_whitney_U": float(mann_whitney.statistic),
                "mann_whitney_pvalue": float(mann_whitney.pvalue),
                "n_matched_pairs": len(paired_ids),
                "paired_mean_difference": (
                    float((paired_tumor - paired_normal).mean())
                    if len(paired_ids)
                    else np.nan
                ),
                "paired_effect_dz": (
                    _paired_effect(paired_tumor, paired_normal)
                    if len(paired_ids) >= 3
                    else np.nan
                ),
                "paired_wilcoxon_pvalue": (
                    float(paired_test.pvalue) if paired_test is not None else np.nan
                ),
                "score_sign_note": (
                    "PC1 sign is oriented to the discovery mean profile but remains "
                    "mathematically arbitrary; interpret magnitude and within-module "
                    "tumor-normal direction, not cross-module score signs."
                ),
            }
        )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["welch_FDR"] = bh_adjust(result["welch_pvalue"])
        result["mann_whitney_FDR"] = bh_adjust(
            result["mann_whitney_pvalue"]
        )
        result["paired_wilcoxon_FDR"] = bh_adjust(
            result["paired_wilcoxon_pvalue"]
        )
        result["all_sample_association_FDR_0.05"] = result[
            "welch_FDR"
        ].le(0.05) & result["mann_whitney_FDR"].le(0.05)
        result["paired_association_FDR_0.05"] = result[
            "paired_wilcoxon_FDR"
        ].le(0.05)
        result["robust_all_and_paired_association"] = result[
            "all_sample_association_FDR_0.05"
        ] & result["paired_association_FDR_0.05"]
        result["association_evidence_class"] = np.select(
            [
                result["robust_all_and_paired_association"],
                result["all_sample_association_FDR_0.05"],
                result["paired_association_FDR_0.05"],
            ],
            [
                "robust_all_samples_and_matched_pairs",
                "all_samples_only_not_paired",
                "matched_pairs_only",
            ],
            default="not_significant",
        )
        result = result.sort_values(
            ["welch_FDR", "mann_whitney_FDR"],
            na_position="last",
        ).reset_index(drop=True)
    return result
