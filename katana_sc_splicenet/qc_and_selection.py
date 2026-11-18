"""Single-cell Excess Biological Variance audit and diversity-constrained feature selection."""
from __future__ import annotations

import numpy as np
import pandas as pd


def estimate_sc_observation_weights_v2(
    usage_values: np.ndarray,
    support_values: np.ndarray | None,
    excess_bio_var: np.ndarray | None = None,
    *,
    prior_bio_var: float = 0.02,
) -> np.ndarray:
    """Compute Empirical-Bayes inverse-variance observation weights w_{ce} = t_{ce} / (1 + t_{ce} * sigma2_bio_e).

    When `support_values` is None, returns 1.0 on finite entries and 0.0 on missing entries.
    Crucially, `w_{ce}` depends ONLY on `support_values` and training-estimated `excess_bio_var`,
    never on `usage_values` themselves, preserving exact zero bias under varying cell library depths.
    """
    u = np.asarray(usage_values, dtype=float)
    finite = np.isfinite(u)
    if support_values is None:
        return finite.astype(float)
    t = np.where(finite, np.maximum(np.asarray(support_values, dtype=float), 0.0), 0.0)
    if excess_bio_var is None:
        sig2 = np.full(u.shape[1], max(float(prior_bio_var), 1e-4), dtype=float)
    else:
        sig2 = np.maximum(np.asarray(excess_bio_var, dtype=float), 1e-4)
    w = np.divide(t, 1.0 + t * sig2[None, :], out=np.zeros_like(t, dtype=float), where=(t > 0.0) & finite)
    return w


def audit_sc_events_v2(
    usage: pd.DataFrame,
    event_genes: pd.Series,
    *,
    support: pd.DataFrame | None = None,
    min_cells: int = 15,
    min_total_support: float = 10.0,
) -> pd.DataFrame:
    """Audit single-cell events/isoforms using count-aware Excess Biological Variance on training cells."""
    u = usage.to_numpy(dtype=float)
    finite = np.isfinite(u)
    obs_cells = finite.sum(axis=0)
    genes = event_genes.reindex(usage.columns).astype(str)
    gene_mult = genes.map(genes.value_counts()).to_numpy(dtype=int)

    if support is not None:
        t = np.where(finite, np.maximum(support.reindex(index=usage.index, columns=usage.columns).to_numpy(dtype=float), 0.0), 0.0)
        tot_sup = t.sum(axis=0)
        # First-pass inverse binomial variance weights (t_ce / (1 + 0.02 * t_ce)) independent of u_ce
        w0 = np.divide(t, 1.0 + 0.02 * t, out=np.zeros_like(t), where=(t > 0) & finite)
    else:
        t = None
        tot_sup = obs_cells.astype(float)
        w0 = finite.astype(float)

    w_sum = w0.sum(axis=0)
    w_sq_sum = (w0 * w0).sum(axis=0)
    n_eff = np.divide(w_sum * w_sum, w_sq_sum, out=np.zeros_like(w_sum), where=w_sq_sum > 0)

    u_zero = np.where(finite, u, 0.0)
    mu = np.divide((w0 * u_zero).sum(axis=0), w_sum, out=np.full(u.shape[1], np.nan), where=w_sum > 0)
    diff = np.where(finite, u - mu[None, :], 0.0)
    raw_wvar = np.divide((w0 * diff * diff).sum(axis=0), w_sum, out=np.full(u.shape[1], np.nan), where=w_sum > 0)
    bessel = np.divide(n_eff, np.maximum(n_eff - 1.0, 1.0), out=np.ones_like(n_eff), where=n_eff > 1.0)
    sample_var = np.where(obs_cells >= 2, raw_wvar * bessel, 0.0)

    unw_var = np.zeros(u.shape[1], dtype=float)
    valid_unw = obs_cells >= 2
    if valid_unw.any():
        u_sub = u[:, valid_unw]
        unw_var[valid_unw] = np.nanvar(u_sub, axis=0, ddof=0)

    if t is not None:
        # Expected binomial sampling variance E[u(1-u)/(t-1) or mu(1-mu)/t]
        # Using plug-in mu_j(1 - mu_j) * sum(w0_cj / t_cj) / sum(w0_cj) avoids downward bias at extreme 0/1 cells!
        inv_t_avg = np.divide(
            np.where((t > 0) & finite, w0 / np.maximum(t, 1.0), 0.0).sum(axis=0),
            w_sum,
            out=np.zeros_like(w_sum),
            where=w_sum > 0,
        )
        binom_floor = np.where(np.isfinite(mu), np.clip(mu * (1.0 - mu), 0.0, 0.25) * inv_t_avg, np.nan)
        excess_var = np.where(np.isfinite(sample_var) & np.isfinite(binom_floor), np.maximum(0.0, sample_var - binom_floor), 0.0)
        noise_model = "binomial_support_exact"
    else:
        binom_floor = np.zeros(u.shape[1], dtype=float)
        excess_var = np.where(np.isfinite(sample_var), np.maximum(0.0, sample_var), 0.0)
        noise_model = "continuous_no_denominator"

    snr = np.divide(excess_var, excess_var + binom_floor + 1e-8)
    qc_pass = (
        (obs_cells >= int(min_cells))
        & (tot_sup >= float(min_total_support))
        & np.isfinite(mu)
        & (unw_var > 1e-8)
    )

    return pd.DataFrame({
        "event_id": usage.columns.astype(str),
        "gene_id": genes.to_numpy(),
        "gene_isoform_count": gene_mult,
        "observed_cells": obs_cells.astype(int),
        "effective_cells": n_eff,
        "total_support": tot_sup,
        "precision_mean": mu,
        "unweighted_var": unw_var,
        "weighted_sample_var": sample_var,
        "binomial_noise_floor": binom_floor,
        "excess_bio_var": excess_var,
        "signal_to_noise": snr,
        "noise_model": noise_model,
        "qc_eligible": qc_pass,
    })


def select_sc_events_v2(
    audit_df: pd.DataFrame,
    *,
    budget: int = 200,
    max_per_gene: int = 2,
    ranking_metric: str = "excess_bio_var",
    require_multi_isoform_gene: bool = False,
    min_excess_bio_var: float = 0.0,
) -> pd.DataFrame:
    """Select diverse, biologically variable single-cell splicing targets strictly from training audit."""
    allowed_metrics = {
        "excess_bio_var",
        "reliability_excess_bio_var",
        "hybrid_manifold_excess_bio_var",
        "psix_knn_score",
        "observed_cells",
        "unweighted_var",
        "signal_to_noise",
    }
    if ranking_metric not in allowed_metrics:
        raise ValueError(f"Unsupported ranking_metric: {ranking_metric}")
    df = audit_df.loc[audit_df["qc_eligible"].eq(True)].copy()
    if require_multi_isoform_gene:
        df = df.loc[df["gene_isoform_count"] >= 2].copy()
    if min_excess_bio_var > 0.0:
        df = df.loc[df["excess_bio_var"] >= float(min_excess_bio_var)].copy()
    if df.empty:
        return df

    df["_tie_id"] = df["event_id"].astype(str)
    df = df.sort_values(
        by=[ranking_metric, "observed_cells", "_tie_id"],
        ascending=[False, False, True],
        kind="mergesort",
    )

    chosen_rows: list[int] = []
    gene_counts: dict[str, int] = {}
    for idx, row in df.iterrows():
        g = str(row["gene_id"])
        if max_per_gene > 0 and gene_counts.get(g, 0) >= max_per_gene:
            continue
        chosen_rows.append(idx)
        gene_counts[g] = gene_counts.get(g, 0) + 1
        if len(chosen_rows) >= int(budget):
            break

    out = df.loc[chosen_rows].drop(columns=["_tie_id"]).reset_index(drop=True)
    return out


def attach_split_state_reliability_and_psix_v2(
    audit_df: pd.DataFrame,
    usage: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    support: pd.DataFrame | None = None,
    host_log_expr: pd.DataFrame | None = None,
    repeats: int = 8,
    knn_k: int = 12,
    seed: int = 20261117,
) -> pd.DataFrame:
    """Compute training-only split-half state-profile repeatability (`split_state_r`),
    coverage-regularized `reliability_excess_bio_var` (suppressing kallisto low-coverage EM noise!),
    and a head-to-head `Psix`-style KNN cell-manifold autocorrelation score (`psix_knn_score`).
    """
    out = audit_df.copy()
    ev_ids = out["event_id"].astype(str).tolist()
    u = usage.reindex(columns=ev_ids).to_numpy(dtype=float)
    finite = np.isfinite(u)
    n_cells, n_ev = u.shape
    obs_frac = finite.mean(axis=0)

    if support is not None:
        t = np.where(finite, np.maximum(support.reindex(index=usage.index, columns=ev_ids).to_numpy(dtype=float), 0.0), 0.0)
        w = np.divide(t, 1.0 + 0.02 * t, out=np.zeros_like(t), where=(t > 0) & finite)
    else:
        w = finite.astype(float)

    # 1. Stratified split-half state/micro-state repeatability on training cells ONLY
    strata = (
        metadata["condition"].astype(str) + "|" + metadata["cell_type"].astype(str)
    ).to_numpy(dtype=str)
    uniq_s = sorted(np.unique(strata))
    rng = np.random.default_rng(int(seed))
    r_splits = np.zeros((repeats, n_ev), dtype=float)
    valid_splits = np.zeros((repeats, n_ev), dtype=bool)

    if len(uniq_s) >= 3:
        for rep in range(repeats):
            h = np.zeros(n_cells, dtype=bool)
            for s_lbl in uniq_s:
                ix = rng.permutation(np.flatnonzero(strata == s_lbl))
                h[ix[::2]] = True
            p1 = np.full((len(uniq_s), n_ev), np.nan)
            p2 = np.full((len(uniq_s), n_ev), np.nan)
            for s_i, s_lbl in enumerate(uniq_s):
                m1 = h & (strata == s_lbl)
                m2 = (~h) & (strata == s_lbl)
                w1 = w[m1].sum(axis=0)
                w2 = w[m2].sum(axis=0)
                ok1 = (finite[m1].sum(axis=0) >= 2) & (w1 > 0)
                ok2 = (finite[m2].sum(axis=0) >= 2) & (w2 > 0)
                p1[s_i] = np.divide((w[m1] * np.where(finite[m1], u[m1], 0.0)).sum(axis=0), w1, out=np.full(n_ev, np.nan), where=ok1)
                p2[s_i] = np.divide((w[m2] * np.where(finite[m2], u[m2], 0.0)).sum(axis=0), w2, out=np.full(n_ev, np.nan), where=ok2)
            both = np.isfinite(p1) & np.isfinite(p2)
            n_b = both.sum(axis=0)
            p1c = np.where(both, p1 - np.divide(np.where(both, p1, 0.0).sum(0), np.maximum(n_b, 1)), 0.0)
            p2c = np.where(both, p2 - np.divide(np.where(both, p2, 0.0).sum(0), np.maximum(n_b, 1)), 0.0)
            den = np.sqrt((p1c * p1c).sum(0) * (p2c * p2c).sum(0))
            ok_r = (n_b >= 3) & (den > 1e-12)
            r_splits[rep] = np.divide((p1c * p2c).sum(0), den, out=np.zeros(n_ev), where=ok_r)
            valid_splits[rep] = ok_r
        v_cnt = valid_splits.sum(axis=0)
        split_r = np.divide((r_splits * valid_splits).sum(axis=0), np.maximum(v_cnt, 1), out=np.zeros(n_ev), where=v_cnt > 0)
    else:
        split_r = np.ones(n_ev, dtype=float) * 0.5

    # Coverage^3 * positive split repeatability * ExcessBioVar selects well-measured, reproducibly modulated isoforms!
    ebv = out["excess_bio_var"].to_numpy(dtype=float)
    rel_ebv = ebv * np.maximum(split_r, 0.02) * (obs_frac ** 3)
    out["observed_fraction"] = obs_frac
    out["split_state_r"] = split_r
    out["reliability_excess_bio_var"] = rel_ebv

    # 2. Psix-style KNN cell-manifold autocorrelation comparator (Buen Abad Najar et al., Genome Biol 2022)
    # Build cell-cell similarity graph from host_log_expr (or usage PCA fallback) on training cells ONLY
    if host_log_expr is not None:
        mat_cell = np.nan_to_num(host_log_expr.reindex(index=usage.index).to_numpy(dtype=float), nan=0.0)
    else:
        col_m = np.nanmean(u, axis=0)
        mat_cell = np.where(finite, u - col_m[None, :], 0.0)
    c_std = mat_cell - mat_cell.mean(axis=1, keepdims=True)
    c_nrm = np.maximum(np.linalg.norm(c_std, axis=1, keepdims=True), 1e-12)
    sim = (c_std / c_nrm) @ (c_std / c_nrm).T
    np.fill_diagonal(sim, -np.inf)
    k_use = min(int(knn_k), max(1, n_cells - 1))
    knn_w = np.zeros((n_cells, n_cells), dtype=float)
    for ci in range(n_cells):
        nn_ix = np.argsort(sim[ci])[-k_use:]
        knn_w[ci, nn_ix] = np.exp(np.maximum(sim[ci, nn_ix], 0.0))
    knn_w = 0.5 * (knn_w + knn_w.T)

    # Standardize cell PSI residuals by binomial SD and compute Moran's I / Psix statistic across KNN graph
    mu = out["precision_mean"].to_numpy(dtype=float)
    z_u = np.where(finite, (u - mu[None, :]) * np.sqrt(w), 0.0)
    num_psix = np.sum(z_u * (knn_w @ z_u), axis=0)
    den_psix = np.maximum(np.sum(z_u * z_u, axis=0), 1e-12) * max(float(knn_w.sum()) / n_cells, 1e-12)
    psix_score = num_psix / den_psix
    out["psix_knn_score"] = psix_score
    # SC CoSplice v2.1 Hybrid Manifold x ExcessBioVar x Coverage score (combines continuous state smoothness with binomial noise subtraction)
    out["hybrid_manifold_excess_bio_var"] = (
        np.maximum(psix_score, 0.0)
        * np.power(np.maximum(ebv, 1e-6), 0.35)
        * np.power(np.maximum(obs_frac, 0.05), 1.5)
    )
    return out
