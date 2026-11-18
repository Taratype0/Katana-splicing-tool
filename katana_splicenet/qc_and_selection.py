"""Donor-equitable QC and Excess Biological Variance feature selection for CoSplice v2.

Replaces raw sample variance ranking with count-aware excess biological variance:
    sigma^2_bio,j = max(0, s^2_j - mean_w[ u_{jm} * (1 - u_{jm}) / t_{jm} ])
which subtracts exact binomial measurement noise before selecting top-K junctions,
preventing low-denominator (t ~ 10-30) coin-flip noise from entering WGCNA.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def unit_robust_audit_v2(
    measured: pd.DataFrame,
    meta: pd.DataFrame,
    read_support: pd.DataFrame | None = None,
    *,
    min_reads: float = 10.0,
    min_observed_fraction: float = 0.80,
    max_unit_variance_share: float = 0.50,
    require_positive_excess_variance: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Screen events for missingness, read depth, and single-donor dominance, and compute excess bio variance."""
    meta_idx = meta.set_index("sample_id").loc[measured.columns]
    x = measured.copy()
    t_mat: np.ndarray | None = None
    if read_support is not None:
        if not x.index.equals(read_support.index) or not x.columns.equals(read_support.columns):
            raise ValueError("Read support matrix index/columns do not match usage matrix")
        sup_vals = read_support.to_numpy(dtype=float)
        obs_sup = sup_vals[np.isfinite(sup_vals)]
        if np.isinf(sup_vals).any() or (obs_sup < 0).any():
            raise ValueError("Read support must contain finite non-negative counts or NA")
        x = x.where(read_support.ge(min_reads))
        t_mat = read_support.where(read_support.ge(min_reads)).to_numpy(dtype=float)

    # Donor-equitable observation weights: each biological unit has total weight 1
    unit_counts = meta_idx["biological_sample_id"].map(meta_idx["biological_sample_id"].value_counts()).to_numpy(dtype=float)
    weights = 1.0 / unit_counts

    a = x.to_numpy(dtype=float)
    ok = np.isfinite(a)
    den = (ok * weights).sum(axis=1)
    mean_w = np.divide(np.nansum(a * weights, axis=1), den, out=np.full(len(a), np.nan), where=den > 0)
    sq_w = np.nan_to_num((a - mean_w[:, None]) ** 2, nan=0.0) * weights

    units = meta_idx.groupby("biological_sample_id", sort=True).indices
    per_unit = np.stack([sq_w[:, idx].sum(axis=1) for idx in units.values()], axis=1)
    total_sq = per_unit.sum(axis=1)
    share = np.divide(per_unit.max(axis=1), total_sq, out=np.full(len(a), np.nan), where=total_sq > 1e-12)

    obs_var = x.var(axis=1, ddof=1).to_numpy(dtype=float)
    if t_mat is not None:
        with np.errstate(divide="ignore", invalid="ignore"):
            b_var_cell = a * (1.0 - a) / t_mat
        ok_b = ok & np.isfinite(b_var_cell) & (t_mat > 0)
        den_b = (ok_b * weights).sum(axis=1)
        binom_noise = np.divide(
            np.nansum(np.where(ok_b, b_var_cell, 0.0) * weights, axis=1),
            den_b,
            out=np.zeros(len(a), dtype=float),
            where=den_b > 0,
        )
        excess_var = np.maximum(0.0, np.nan_to_num(obs_var, nan=0.0) - binom_noise)
        median_depth = np.full(len(a), np.nan)
        any_ok = ok.any(axis=1)
        if any_ok.any():
            median_depth[any_ok] = np.nanmedian(np.where(ok[any_ok], t_mat[any_ok], np.nan), axis=1)
    else:
        binom_noise = np.full(len(a), np.nan)
        excess_var = np.nan_to_num(obs_var, nan=0.0)
        median_depth = np.full(len(a), np.nan)

    obs_frac = ok.mean(axis=1)
    eligible = (obs_frac >= min_observed_fraction) & (obs_var > 1e-12) & (share < max_unit_variance_share)
    if t_mat is not None and require_positive_excess_variance:
        eligible = eligible & (excess_var > 1e-12)

    reasons = []
    for o, v, ev, s in zip(obs_frac, obs_var, excess_var, share):
        r = []
        if o < min_observed_fraction:
            r.append("missingness_or_read_support")
        if not (v > 1e-12):
            r.append("constant")
        elif t_mat is not None and require_positive_excess_variance and not (ev > 1e-12):
            r.append("pure_binomial_sampling_noise")
        if not (s < max_unit_variance_share):
            r.append("unit_dominated_variation")
        reasons.append(";".join(r) if r else "passed")

    audit = pd.DataFrame(
        {
            "observed_fraction": obs_frac,
            "variance": obs_var,
            "binomial_noise_variance": binom_noise,
            "excess_bio_variance": excess_var,
            "median_read_depth": median_depth,
            "maximum_unit_variance_fraction": share,
            "highest_contribution_unit": np.asarray(list(units))[per_unit.argmax(axis=1)],
            "masked_low_support_cells": measured.notna().sum(axis=1).to_numpy() - ok.sum(axis=1),
            "eligible": eligible,
            "reason": reasons,
        },
        index=x.index,
    )
    return x, audit


def select_events_v2(
    measured: pd.DataFrame,
    audit: pd.DataFrame,
    feature_cap: int | None = 5000,
    *,
    ranking_metric: str = "excess_bio_variance",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Select top-K QC-eligible events by excess biological variance (or raw variance for ablation)."""
    if feature_cap is not None and (not isinstance(feature_cap, (int, np.integer)) or feature_cap < 20):
        raise ValueError("feature_cap must be an integer >= 20 or None")
    if ranking_metric not in {"excess_bio_variance", "variance"}:
        raise ValueError(f"Unsupported ranking_metric: {ranking_metric!r}")

    aud = audit.copy()
    eligible_idx = aud.index[aud["eligible"]]
    score_col = ranking_metric if aud[ranking_metric].notna().any() else "variance"
    ranks = pd.DataFrame(
        {"score": aud.loc[eligible_idx, score_col].to_numpy(dtype=float), "_tie_key": eligible_idx.astype(str)},
        index=eligible_idx,
    ).sort_values(["score", "_tie_key"], ascending=[False, True], kind="stable")

    selected = ranks.index if feature_cap is None else ranks.index[:feature_cap]
    aud["qc_eligible"] = aud["eligible"]
    aud["ranking_metric"] = score_col
    aud["selection_rank"] = pd.Series(np.arange(1, len(ranks) + 1), index=ranks.index).reindex(aud.index)
    aud["network_selected"] = aud.index.isin(selected)
    aud["selection_reason"] = np.where(
        ~aud["eligible"],
        "failed_qc",
        np.where(aud["network_selected"], "selected", "outside_user_event_limit"),
    )
    return measured.loc[measured.index.isin(selected)], aud
