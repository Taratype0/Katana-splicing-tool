"""Publication-grade Splicing Landscape Algorithms & Executive Renderers:

1. Beyond-Jutils Precision-Shrunk PSI Z-Score Clustermap & PERMANOVA Confidence-Ellipse PCA:
   - Addresses the 4 methodological flaws of Jutils (Yang et al., Bioinformatics 2021):
     (a) Avoids complete-case `dropna()` truncation via coverage-filtered Empirical-Bayes prior shrinkage
         `PSI_tilde_{e,s} = (N_{e,s} * PSI_{e,s} + kappa * mu_e) / (N_{e,s} + kappa)`.
     (b) Selects features unsupervised via `ExcessBioVar = max(0, Var_obs(PSI) - E[PSI*(1-PSI)/(N+1)])`
         so low-depth binomial count noise is never mistaken for biological variance.
     (c) Computes precision-stabilized row Z-scores (`[-z_clip, +z_clip]`) with Ward hierarchical clustering,
         condition + coverage top strips, and event-type left strip.
     (d) Evaluates condition separation via permutation PERMANOVA pseudo-F / R^2 and 95% covariance ellipses
         without supervised pre-filtering leakage (while also offering an explicit Supervised Contrast Signature mode).

2. Precision-Kernel Event-PSI UMAP & Cell-Type Modality Divergence (Single-Cell / Large Cohort):
   - Computes pairwise cell/sample splicing distance exclusively across co-measured high-ExcessBioVar events
     weighted by joint Binomial Fisher precision `w_{ij,e} = N_{i,e} * N_{j,e} / (N_{i,e} + N_{j,e} + kappa)`
     (eliminating MARVEL / naive mean-imputation dropout blobs).
   - Aligns Expression UMAP and Precision Event-PSI UMAP via Procrustes rotation and computes per-cell
     Cross-Modality Neighborhood Divergence `Delta_kNN(i) = 1 - |N_k^expr(i) cap N_k^psi(i)| / k`,
     Cell-Type Silhouette, and Donor Inverse Simpson Index (LISI) batch-mixing score.
"""
from __future__ import annotations

import math
from pathlib import Path
import re
from typing import Any

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
from scipy.linalg import orthogonal_procrustes
from scipy.spatial.distance import pdist, squareform
from scipy.stats import chi2

MCK_NAVY = "#051C2C"
MCK_BLUE = "#2251FF"
MCK_CERULEAN = "#0077C8"
MCK_CORAL = "#B83A4B"
MCK_TEAL = "#0D9488"
MCK_AMBER = "#D97706"
MCK_VIOLET = "#4F46E5"
MCK_SLATE = "#64748B"

CONDITION_PALETTE = [
    "#051C2C", "#2251FF", "#C83E4D", "#0085A1",
    "#D97706", "#4F46E5", "#0F766E", "#475569",
]
# Deliberately distinct editorial jewel/slate palette so Event Subtypes NEVER share hex colors with Sample Conditions
EVENT_TYPE_COLORS = {
    "SE": "#64748B",
    "RI": "#8B5CF6",
    "A3SS": "#10B981",
    "A5SS": "#F59E0B",
    "MXE": "#EC4899",
    "AS": "#94A3B8",
}


def infer_event_type(eid: str) -> str:
    s = str(eid or "")
    for et in ("A3SS", "A5SS", "MXE", "RI", "SE"):
        if et in s:
            return et
    return "SE"


def permanova_r2(dist_mat: np.ndarray, groups: list[str], n_perm: int = 499, seed: int = 42) -> tuple[float, float]:
    """Non-parametric Anderson (2001) PERMANOVA R^2 and permutation p-value on distance matrix."""
    n = len(groups)
    if n < 4 or len(set(groups)) < 2:
        return 0.0, 1.0
    D2 = np.asarray(dist_mat, dtype=float) ** 2
    H = np.eye(n) - np.ones((n, n)) / n
    G = -0.5 * (H @ D2 @ H)
    ss_total = float(np.trace(G))
    if ss_total <= 1e-12:
        return 0.0, 1.0

    uniq = sorted(set(groups))
    g_arr = np.asarray(groups)

    def _ss_within(labels: np.ndarray) -> float:
        ss_w = 0.0
        for u in uniq:
            idx = np.where(labels == u)[0]
            if len(idx) > 1:
                sub_d2 = D2[np.ix_(idx, idx)]
                ss_w += float(sub_d2.sum()) / (2.0 * len(idx))
        return ss_w

    ss_t = float(D2.sum()) / (2.0 * n)
    ss_w_obs = _ss_within(g_arr)
    ss_b_obs = max(0.0, ss_t - ss_w_obs)
    r2_obs = ss_b_obs / max(1e-12, ss_t)

    rng = np.random.default_rng(seed)
    ge_count = 1
    for _ in range(n_perm):
        perm = rng.permutation(g_arr)
        ss_b_p = max(0.0, ss_t - _ss_within(perm))
        if ss_b_p >= ss_b_obs - 1e-12:
            ge_count += 1
    return round(r2_obs, 4), round(ge_count / (n_perm + 1), 4)


def compute_beyond_jutils_landscape(
    psi_matrix: pd.DataFrame,
    sample_meta: list[dict[str, Any]],
    *,
    depth_matrix: pd.DataFrame | None = None,
    event_meta: dict[str, dict[str, Any]] | None = None,
    supervised_event_ids: set[str] | None = None,
    top_n: int = 150,
    manifold_top_n: int = 2000,
    min_present_fraction: float = 0.65,
    prior_kappa: float = 10.0,
    z_clip: float = 2.5,
    metric: str = "euclidean",
    event_types: list[str] | None = None,
    run_umap: bool = True,
    umap_neighbors: int = 12,
    umap_min_dist: float = 0.18,
) -> dict[str, Any]:
    """Compute Empirical-Bayes precision-shrunk PSI Z-scores, hierarchical clustermap, PCA, and Precision-Kernel UMAP."""
    s_ids = [str(s["sample_id"]) for s in sample_meta if str(s["sample_id"]) in psi_matrix.columns]
    if len(s_ids) < 3 or psi_matrix.empty:
        return {"available": False, "reason": "At least 3 measured samples/cells are required for landscape analysis."}

    meta_by_id = {str(s["sample_id"]): s for s in sample_meta}
    P_raw = psi_matrix[s_ids].apply(pd.to_numeric, errors="coerce")
    P_raw.index = psi_matrix.index.astype(str)

    if event_types:
        allowed_et = {e.upper() for e in event_types}
        et_series = pd.Series([
            str((event_meta or {}).get(eid, {}).get("event_type") or infer_event_type(eid)).upper()
            for eid in P_raw.index
        ], index=P_raw.index)
        P_raw = P_raw.loc[et_series.isin(allowed_et)]

    present_frac = P_raw.notna().mean(axis=1)
    P_obs = P_raw.loc[present_frac >= min_present_fraction].copy()
    if len(P_obs) < 5:
        P_obs = P_raw.loc[present_frac >= 0.25].copy()
    if len(P_obs) < 3:
        return {"available": False, "reason": "Insufficient multi-sample splicing events after coverage filtering."}

    # Align junction read depth matrix N_{e,s}; if absent in compact matrix, use conservative effective depth N=25 where observed, 0 where NA
    if depth_matrix is not None and not depth_matrix.empty:
        N_mat = depth_matrix.reindex(index=P_obs.index, columns=s_ids).fillna(0.0).to_numpy(dtype=float)
    else:
        N_mat = np.where(P_obs.notna().to_numpy(), 25.0, 0.0)

    P_np = P_obs.to_numpy(dtype=float)
    obs_mask = np.isfinite(P_np) & (N_mat > 0)

    # Weighted event prior mean mu_e
    w_sum = np.where(obs_mask, N_mat, 0.0).sum(axis=1, keepdims=True)
    mu_e = np.divide(
        np.where(obs_mask, N_mat * np.nan_to_num(P_np, nan=0.0), 0.0).sum(axis=1, keepdims=True),
        np.maximum(1.0, w_sum),
    )
    mu_e = np.clip(mu_e, 0.01, 0.99)

    # Empirical-Bayes precision-shrunk PSI_tilde
    P_shrunk = np.where(
        obs_mask,
        (N_mat * np.nan_to_num(P_np, nan=0.0) + prior_kappa * mu_e) / (N_mat + prior_kappa),
        mu_e,
    )

    # ExcessBioVar = max(0, Var_obs - E[mu*(1-mu)/(N+1)])
    obs_var = np.nanvar(np.where(obs_mask, P_np, np.nan), axis=1)
    binom_noise = np.nanmean(
        np.where(obs_mask, (mu_e * (1.0 - mu_e)) / np.maximum(1.0, N_mat + 1.0), np.nan),
        axis=1,
    )
    excess_biovar = np.maximum(0.0, np.nan_to_num(obs_var - binom_noise, nan=0.0))

    eids_all = P_obs.index.to_numpy()
    man_cap = max(int(top_n), int(manifold_top_n))
    if supervised_event_ids:
        sup_mask = np.isin(eids_all, list(supervised_event_ids))
        if sup_mask.sum() >= 5:
            cand_idx = np.where(sup_mask)[0]
            man_idx = cand_idx[np.argsort(-excess_biovar[cand_idx])[: max(5, man_cap)]]
            hm_idx = man_idx[: max(5, int(top_n))]
            selection_mode = "supervised_ps_das"
        else:
            man_idx = np.argsort(-excess_biovar)[: max(5, man_cap)]
            hm_idx = man_idx[: max(5, int(top_n))]
            selection_mode = "unsupervised_excess_biovar"
    else:
        man_idx = np.argsort(-excess_biovar)[: max(5, man_cap)]
        hm_idx = man_idx[: max(5, int(top_n))]
        selection_mode = "unsupervised_excess_biovar"

    # Global Manifold (PCA, PERMANOVA & Precision-Kernel UMAP) uses `man_idx` (e.g. 2,000 HVEs or all PS-DAS events!)
    P_man = P_shrunk[man_idx, :]
    N_man = N_mat[man_idx, :]
    mask_man = obs_mask[man_idx, :]

    # Heatmap Display uses `hm_idx` (top `top_n` events for readable row resolution)
    P_sel = P_shrunk[hm_idx, :]
    N_sel = N_mat[hm_idx, :]
    mask_sel = obs_mask[hm_idx, :]
    sel_eids = eids_all[hm_idx].tolist()

    # Row-wise precision-stabilized Z-score for Heatmap & Sample Dendrogram
    man_mean = P_man.mean(axis=1, keepdims=True)
    man_sd = np.maximum(0.025, P_man.std(axis=1, keepdims=True))
    Z_man = np.clip((P_man - man_mean) / man_sd, -float(z_clip), float(z_clip))

    row_mean = P_sel.mean(axis=1, keepdims=True)
    row_sd = np.maximum(0.025, P_sel.std(axis=1, keepdims=True))
    Z_mat = np.clip((P_sel - row_mean) / row_sd, -float(z_clip), float(z_clip))

    # Hierarchical clustering of samples (columns on global manifold) and heatmap events (rows)
    dist_metric = metric if metric in ("euclidean", "cosine", "correlation") else "euclidean"
    try:
        col_dist_vec = pdist(Z_man.T, metric=dist_metric)
        col_dist_vec = np.nan_to_num(col_dist_vec, nan=0.0)
        col_link = linkage(col_dist_vec, method="ward" if dist_metric == "euclidean" else "average")
        col_order = leaves_list(col_link).tolist()
    except Exception:
        col_link = None
        col_order = list(range(len(s_ids)))

    try:
        row_dist_vec = pdist(Z_mat, metric=dist_metric)
        row_dist_vec = np.nan_to_num(row_dist_vec, nan=0.0)
        row_link = linkage(row_dist_vec, method="ward" if dist_metric == "euclidean" else "average")
        row_order = leaves_list(row_link).tolist()
    except Exception:
        row_link = None
        row_order = list(range(len(sel_eids)))

    # PCA on centered shrunk global PSI (samples x manifold_events)
    X_pca = P_man.T - P_man.T.mean(axis=0, keepdims=True)
    u_s, s_vals, _ = np.linalg.svd(X_pca, full_matrices=False)
    coords_pca = u_s * s_vals
    var_tot = float((s_vals ** 2).sum()) or 1.0
    var_exp = [round(float((sv ** 2) / var_tot * 100.0), 2) for sv in s_vals[:5]]
    while len(var_exp) < 3:
        var_exp.append(0.0)

    # Precision-Kernel Pairwise Distance Matrix D_PSI on global manifold events
    n_s = len(s_ids)
    D_prec = np.zeros((n_s, n_s), dtype=float)
    W_eff = np.where(mask_man, N_man, 0.0)  # (E_man, S)
    for i in range(n_s):
        wi = W_eff[:, i : i + 1]
        wj = W_eff[:, i:]
        joint_w = (wi * wj) / np.maximum(1.0, wi + wj + prior_kappa)
        diff2 = (P_man[:, i : i + 1] - P_man[:, i:]) ** 2
        denom = joint_w.sum(axis=0)
        d_row = np.sqrt(np.divide((joint_w * diff2).sum(axis=0), np.maximum(1e-8, denom)))
        fallback = np.sqrt(np.mean((P_man[:, i : i + 1] - P_man[:, i:]) ** 2, axis=0))
        d_final = np.where(denom > 1e-6, d_row, fallback)
        D_prec[i, i:] = d_final
        D_prec[i:, i] = d_final
    np.fill_diagonal(D_prec, 0.0)

    conditions = [str(meta_by_id[sid].get("condition", "Unknown")) for sid in s_ids]
    perm_r2, perm_p = permanova_r2(D_prec, conditions)

    # Precision-Kernel UMAP (if requested and n_samples >= 6)
    umap_coords = np.zeros((n_s, 2), dtype=float)
    umap_method = "PCA_2D_fallback"
    if run_umap and n_s >= 6:
        try:
            import umap
            nn = max(3, min(int(umap_neighbors), n_s - 1))
            reducer = umap.UMAP(
                n_neighbors=nn,
                min_dist=float(umap_min_dist),
                metric="precomputed",
                random_state=42,
            )
            umap_coords = reducer.fit_transform(D_prec)
            umap_method = f"Precision-Kernel UMAP (n_neighbors={nn}, min_dist={umap_min_dist})"
        except Exception:
            umap_coords = coords_pca[:, :2]

    # Sample records
    sample_out = []
    mean_cov_per_sample = mask_man.mean(axis=0)
    for idx, sid in enumerate(s_ids):
        m = meta_by_id[sid]
        sample_out.append({
            "sample_id": sid,
            "condition": str(m.get("condition", "")),
            "donor_id": str(m.get("donor_id", m.get("biological_sample_id", ""))),
            "cell_type": str(m.get("cell_type", "")),
            "observed_event_fraction": round(float(mean_cov_per_sample[idx]), 4),
            "PC1": round(float(coords_pca[idx, 0]), 4),
            "PC2": round(float(coords_pca[idx, 1]) if coords_pca.shape[1] > 1 else 0.0, 4),
            "PC3": round(float(coords_pca[idx, 2]) if coords_pca.shape[1] > 2 else 0.0, 4),
            "psi_UMAP1": round(float(umap_coords[idx, 0]), 4),
            "psi_UMAP2": round(float(umap_coords[idx, 1]), 4),
        })

    # Ordered heatmap payload (capped to 160 rows in JSON for instant SVG rendering; full Z_mat used for PNG)
    ordered_samples = [s_ids[i] for i in col_order]
    ordered_conditions = [conditions[i] for i in col_order]
    ordered_rows = []
    for r_i in row_order[:160]:
        eid = sel_eids[r_i]
        em = (event_meta or {}).get(eid, {})
        ordered_rows.append({
            "event_id": eid,
            "gene_symbol": str(em.get("gene_symbol") or eid.split(":")[0]),
            "event_type": str(em.get("event_type") or infer_event_type(eid)),
            "excess_biovar": round(float(excess_biovar[hm_idx[r_i]]), 5),
            "z_values": [round(float(Z_mat[r_i, c_i]), 3) for c_i in col_order],
            "psi_values": [round(float(P_sel[r_i, c_i]), 3) for c_i in col_order],
        })

    return {
        "available": True,
        "selection_mode": selection_mode,
        "metric": dist_metric,
        "prior_kappa": prior_kappa,
        "z_clip": z_clip,
        "umap_method": umap_method,
        "total_candidate_events": int(len(P_obs)),
        "selected_events": int(len(man_idx)),
        "manifold_events": int(len(man_idx)),
        "heatmap_events": int(len(sel_eids)),
        "permanova_r2": perm_r2,
        "permanova_p": perm_p,
        "pca_variance_pct": var_exp,
        "samples": sample_out,
        "ordered_samples": ordered_samples,
        "ordered_conditions": ordered_conditions,
        "heatmap_rows": ordered_rows,
        "_internal": {
            "Z_ordered": Z_mat[np.ix_(row_order, col_order)],
            "col_link": col_link,
            "col_order": col_order,
            "row_order": row_order,
            "sel_eids": [sel_eids[i] for i in row_order],
            "obs_frac_ordered": mean_cov_per_sample[col_order],
            "D_prec": D_prec,
        },
    }


def compute_dual_expression_vs_psi_umap(
    expr_df: pd.DataFrame,
    psi_landscape: dict[str, Any],
    k_neighbors: int = 12,
) -> dict[str, Any]:
    """Align Gene-Expression UMAP with Precision-Kernel Event-PSI UMAP via Procrustes & compute per-cell Divergence."""
    from sklearn.metrics import silhouette_score

    samples = psi_landscape["samples"]
    s_ids = [s["sample_id"] for s in samples]
    n_s = len(s_ids)
    expr_num = (
        expr_df.drop(columns=[c for c in ("row_id", "Unnamed: 0") if c in expr_df.columns])
        .apply(pd.to_numeric, errors="coerce")
        .dropna(axis=1, how="all")
    )
    X_e = expr_num.reindex(s_ids).fillna(0.0).to_numpy(dtype=float)
    X_e = X_e - X_e.mean(axis=0, keepdims=True)

    if n_s >= 6:
        try:
            import umap
            nn = max(4, min(int(k_neighbors), n_s - 1))
            expr_umap = umap.UMAP(n_neighbors=nn, min_dist=0.22, metric="cosine", random_state=42).fit_transform(X_e)
        except Exception:
            u, sv, _ = np.linalg.svd(X_e, full_matrices=False)
            expr_umap = (u * sv)[:, :2]
    else:
        u, sv, _ = np.linalg.svd(X_e, full_matrices=False)
        expr_umap = (u * sv)[:, :2]

    psi_umap = np.array([[s["psi_UMAP1"], s["psi_UMAP2"]] for s in samples], dtype=float)

    # Standardize & Procrustes-align PSI UMAP onto Expression UMAP orientation for side-by-side visual comparison
    e_std = (expr_umap - expr_umap.mean(axis=0)) / np.maximum(1e-6, expr_umap.std(axis=0))
    p_std = (psi_umap - psi_umap.mean(axis=0)) / np.maximum(1e-6, psi_umap.std(axis=0))
    R_proc, _ = orthogonal_procrustes(p_std, e_std)
    p_aligned = p_std @ R_proc

    # Cross-modality kNN Divergence Score Delta_kNN(i) = 1 - |N_k^expr(i) cap N_k^psi(i)| / k
    k = max(3, min(int(k_neighbors), n_s - 1))
    D_e = squareform(pdist(X_e, metric="euclidean"))
    D_p = psi_landscape["_internal"]["D_prec"]

    divergences = []
    donor_lisi = []
    donors = [s.get("donor_id", "") for s in samples]
    for i in range(n_s):
        nn_e = set(np.argsort(D_e[i])[1 : k + 1].tolist())
        nn_p = set(np.argsort(D_p[i])[1 : k + 1].tolist())
        div = 1.0 - (len(nn_e & nn_p) / float(k))
        divergences.append(round(div, 4))
        # Local Inverse Simpson Index (LISI) on PSI neighborhood for Donor batch mixing (ideal for 2 equal donors = 2.0)
        d_neigh = [donors[j] for j in nn_p if donors[j]]
        if d_neigh:
            probs = pd.Series(d_neigh).value_counts(normalize=True).to_numpy()
            donor_lisi.append(round(float(1.0 / np.sum(probs ** 2)), 3))
        else:
            donor_lisi.append(1.0)

    cell_types = [s.get("cell_type") or s.get("condition") or "Cell" for s in samples]
    sil_psi = round(float(silhouette_score(D_p, cell_types, metric="precomputed")), 4) if len(set(cell_types)) >= 2 else 0.0
    sil_expr = round(float(silhouette_score(D_e, cell_types, metric="precomputed")), 4) if len(set(cell_types)) >= 2 else 0.0

    enriched_cells = []
    for idx, s in enumerate(samples):
        enriched_cells.append({
            **s,
            "expr_UMAP1": round(float(e_std[idx, 0]), 4),
            "expr_UMAP2": round(float(e_std[idx, 1]), 4),
            "psi_UMAP1": round(float(p_aligned[idx, 0]), 4),
            "psi_UMAP2": round(float(p_aligned[idx, 1]), 4),
            "modality_divergence": divergences[idx],
            "donor_lisi": donor_lisi[idx],
        })

    return {
        "cells": enriched_cells,
        "k_neighbors": k,
        "mean_divergence": round(float(np.mean(divergences)), 4),
        "median_donor_lisi": round(float(np.median(donor_lisi)), 3),
        "cell_type_silhouette_psi": sil_psi,
        "cell_type_silhouette_expr": sil_expr,
    }


def render_mckinsey_psi_clustermap(
    landscape: dict[str, Any],
    png_path: str | Path,
    svg_path: str | Path | None = None,
    title: str = "Depth-Stabilized Splicing Z-Score Clustermap",
    event_meta: dict[str, dict[str, Any]] | None = None,
) -> None:
    """Render publication-grade PSI Z-Score Clustermap with distinct Condition vs Subtype legends & gene callouts."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.patches import Patch

    png_path = Path(png_path)
    png_path.parent.mkdir(parents=True, exist_ok=True)

    internal = landscape["_internal"]
    Z = internal["Z_ordered"]
    col_link = internal["col_link"]
    sel_eids = internal["sel_eids"]
    obs_frac = internal["obs_frac_ordered"]
    samples = landscape["ordered_samples"]
    conds = landscape["ordered_conditions"]
    z_clip = float(landscape.get("z_clip", 2.5))
    show_coverage_strip = bool(len(obs_frac) > 0 and float(np.min(obs_frac)) < 0.97)

    uniq_conds = list(dict.fromkeys(conds))
    cond_color = {c: CONDITION_PALETTE[i % len(CONDITION_PALETTE)] for i, c in enumerate(uniq_conds)}

    fig = plt.figure(figsize=(13.2, 8.4), dpi=220)
    fig.patch.set_facecolor("white")
    if show_coverage_strip:
        gs = fig.add_gridspec(
            4, 4,
            height_ratios=[0.90, 0.24, 0.18, 6.2],
            width_ratios=[0.22, 7.5, 1.55, 0.25],
            hspace=0.04, wspace=0.05,
        )
        hm_row_idx = 3
    else:
        gs = fig.add_gridspec(
            3, 4,
            height_ratios=[0.92, 0.26, 6.3],
            width_ratios=[0.22, 7.5, 1.55, 0.25],
            hspace=0.04, wspace=0.05,
        )
        hm_row_idx = 2

    # 1. Top Sample Dendrogram
    ax_dend = fig.add_subplot(gs[0, 1])
    if col_link is not None and len(samples) >= 3:
        dendrogram(col_link, ax=ax_dend, color_threshold=0, above_threshold_color=MCK_NAVY, no_labels=True)
    ax_dend.axis("off")

    # 2. Sample Group Color Bar
    ax_cond = fig.add_subplot(gs[1, 1])
    cond_rgb = np.array([[matplotlib.colors.to_rgb(cond_color[c]) for c in conds]])
    ax_cond.imshow(cond_rgb, aspect="auto", interpolation="nearest")
    ax_cond.set_xticks([])
    ax_cond.set_yticks([0], ["Sample Group"], fontsize=8.0, fontweight="bold", color=MCK_NAVY)
    for sp in ax_cond.spines.values():
        sp.set_visible(False)

    # 3. Optional Cell Completeness Bar (only shown in sparse single-cell cohorts where <97% complete)
    if show_coverage_strip:
        ax_cov = fig.add_subplot(gs[2, 1])
        cov_cmap = LinearSegmentedColormap.from_list("cov_mck", ["#E2E8F0", "#93C5FD", "#051C2C"])
        ax_cov.imshow(obs_frac[None, :], aspect="auto", cmap=cov_cmap, vmin=0.35, vmax=1.0, interpolation="nearest")
        ax_cov.set_xticks([])
        ax_cov.set_yticks([0], ["Detected %"], fontsize=7.8, color=MCK_SLATE)
        for sp in ax_cov.spines.values():
            sp.set_visible(False)

    # 4. Left Event-Subtype Strip (distinct jewel/slate palette!)
    ax_et = fig.add_subplot(gs[hm_row_idx, 0])
    etypes = [str((event_meta or {}).get(e, {}).get("event_type") or infer_event_type(e)) for e in sel_eids]
    et_rgb = np.array([[matplotlib.colors.to_rgb(EVENT_TYPE_COLORS.get(et, "#94A3B8"))] for et in etypes])
    ax_et.imshow(et_rgb, aspect="auto", interpolation="nearest")
    ax_et.set_xticks([0], ["Subtype"], fontsize=7.6, rotation=90, color=MCK_SLATE)
    ax_et.set_yticks([])
    for sp in ax_et.spines.values():
        sp.set_visible(False)

    # 5. Main Diverging PSI Z-Score Heatmap
    ax_hm = fig.add_subplot(gs[hm_row_idx, 1])
    div_cmap = LinearSegmentedColormap.from_list("mck_z", [MCK_CERULEAN, "#E0F2FE", "#FFFFFF", "#FECDD3", MCK_CORAL])
    im = ax_hm.imshow(Z, aspect="auto", cmap=div_cmap, vmin=-z_clip, vmax=z_clip, interpolation="nearest")
    if len(samples) <= 32:
        ax_hm.set_xticks(np.arange(len(samples)))
        ax_hm.set_xticklabels(samples, rotation=90, fontsize=8.2, color=MCK_NAVY)
    else:
        ax_hm.set_xticks([])
        ax_hm.set_xlabel(f"{len(samples)} cells / samples ordered by Ward linkage", fontsize=8.8, color=MCK_NAVY)

    # Always display clean Gene Symbol labels on the right axis (all rows if <=32, or ~20 evenly-spaced representative genes!)
    n_rows = len(sel_eids)
    if n_rows <= 32:
        tick_idx = list(range(n_rows))
    else:
        step = max(1, n_rows // 20)
        tick_idx = []
        seen_genes: set[str] = set()
        for r_i in range(0, n_rows, step):
            g_sym = str((event_meta or {}).get(sel_eids[r_i], {}).get("gene_symbol", "")).strip()
            if g_sym and not g_sym.startswith("ENS") and g_sym not in seen_genes:
                seen_genes.add(g_sym)
                tick_idx.append(r_i)
    row_lbls = [
        f"{(event_meta or {}).get(sel_eids[i], {}).get('gene_symbol', sel_eids[i][:12])} ({etypes[i]})"
        for i in tick_idx
    ]
    ax_hm.set_yticks(tick_idx, row_lbls, fontsize=7.4, color=MCK_NAVY)
    ax_hm.yaxis.tick_right()

    # 6. Colorbar (placed in column 3 after dedicated column-2 label gutter so right-axis gene names NEVER collide!)
    ax_cb = fig.add_subplot(gs[hm_row_idx, 3])
    cb = fig.colorbar(im, cax=ax_cb)
    cb.set_label(f"Depth-Stabilized PSI Row Z-Score (±{z_clip:.1f} SD)", fontsize=8.4, color=MCK_NAVY)
    cb.ax.tick_params(labelsize=7.6)
    cb.outline.set_edgecolor("#CBD5E1")

    # 7. TWO SEPARATE, CLEARLY TITLED LEGENDS (Sample Group vs Splicing Subtype — zero color collision!)
    cond_patches = [Patch(facecolor=cond_color[c], label=c) for c in uniq_conds[:8]]
    present_ets = [et for et in ("SE", "RI", "A3SS", "A5SS", "MXE") if et in set(etypes)]
    et_patches = [Patch(facecolor=EVENT_TYPE_COLORS[et], label=et) for et in present_ets]

    leg1 = fig.legend(
        handles=cond_patches,
        title="Sample Group (Top Bar)",
        title_fontsize=8.0,
        loc="upper left",
        bbox_to_anchor=(0.04, 0.925),
        ncol=min(6, max(1, len(cond_patches))),
        frameon=True, facecolor="#F8FAFC", edgecolor="#E2E8F0", fontsize=7.8,
    )
    leg1.get_title().set_fontweight("bold")
    fig.add_artist(leg1)

    if et_patches:
        leg2 = fig.legend(
            handles=et_patches,
            title="Splicing Subtype (Left Strip)",
            title_fontsize=8.0,
            loc="upper right",
            bbox_to_anchor=(0.94, 0.925),
            ncol=len(et_patches),
            frameon=True, facecolor="#F8FAFC", edgecolor="#E2E8F0", fontsize=7.8,
        )
        leg2.get_title().set_fontweight("bold")

    man_n = landscape.get("manifold_events", landscape["selected_events"])
    hm_n = landscape.get("heatmap_events", len(sel_eids))
    mode_lbl = "Unsupervised biological variance" if "unsupervised" in landscape["selection_mode"] else "Supervised PS-DAS signature"
    fig.text(0.04, 0.978, title, fontsize=11.8, fontweight="bold", color=MCK_NAVY, ha="left", va="top")
    fig.text(
        0.04, 0.948,
        f"Top {hm_n} events shown (sample clustering & PERMANOVA R² = {landscape['permanova_r2']*100:.1f}%, P = {landscape['permanova_p']:.3g} computed across {man_n:,} {mode_lbl} events)",
        fontsize=8.3, color=MCK_SLATE, ha="left", va="top",
    )
    fig.subplots_adjust(left=0.05, right=0.88, top=0.84, bottom=0.11)
    fig.savefig(png_path, dpi=220, bbox_inches="tight", facecolor="white")
    if svg_path:
        fig.savefig(svg_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _draw_confidence_ellipse(ax, xs: np.ndarray, ys: np.ndarray, color: str, level: float = 0.90) -> None:
    from matplotlib.patches import Ellipse
    if len(xs) < 3:
        return
    cov = np.cov(xs, ys)
    if not np.all(np.isfinite(cov)):
        return
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    if vals[0] <= 1e-10:
        return
    scale = math.sqrt(chi2.ppf(level, df=2))
    width, height = 2 * scale * np.sqrt(np.maximum(1e-10, vals))
    angle = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
    ell = Ellipse(
        xy=(float(np.mean(xs)), float(np.mean(ys))),
        width=float(width), height=float(height), angle=float(angle),
        facecolor=color, alpha=0.09, edgecolor=color, linewidth=1.1, linestyle="--", zorder=1,
    )
    ax.add_patch(ell)


def render_mckinsey_psi_pca(
    landscape: dict[str, Any],
    png_path: str | Path,
    svg_path: str | Path | None = None,
    title: str = "Precision-Shrunk Splicing Principal Components & Condition Separation",
) -> None:
    """Render publication 2-panel PSI PCA (PC1 vs PC2 with 90% confidence ellipses + PERMANOVA & Scree spectrum)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    png_path = Path(png_path)
    png_path.parent.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9.0, "axes.spines.top": False, "axes.spines.right": False})
    samples = landscape["samples"]
    conds = list(dict.fromkeys(s["condition"] for s in samples))
    cond_color = {c: CONDITION_PALETTE[i % len(CONDITION_PALETTE)] for i, c in enumerate(conds)}
    markers = ["o", "s", "^", "D", "P", "X", "v", "<"]

    fig, (ax_pca, ax_scree) = plt.subplots(
        1, 2, figsize=(12.2, 4.9),
        gridspec_kw={"width_ratios": [2.35, 1.0], "wspace": 0.24}, dpi=220,
    )
    fig.patch.set_facecolor("white")

    ax_pca.axhline(0, color="#E2E8F0", linewidth=0.8, zorder=0)
    ax_pca.axvline(0, color="#E2E8F0", linewidth=0.8, zorder=0)
    ax_pca.grid(True, color="#F1F5F9", linewidth=0.6, zorder=0)

    for idx, c in enumerate(conds):
        pts = [s for s in samples if s["condition"] == c]
        xs = np.array([p["PC1"] for p in pts], dtype=float)
        ys = np.array([p["PC2"] for p in pts], dtype=float)
        col = cond_color[c]
        mk = markers[idx % len(markers)]
        _draw_confidence_ellipse(ax_pca, xs, ys, col, level=0.90)
        ax_pca.scatter(xs, ys, s=78, color=col, marker=mk, edgecolors="white", linewidths=0.9, label=f"{c} (n={len(pts)})", zorder=3)
        if len(samples) <= 30:
            for s_i, p in enumerate(pts):
                ax_pca.annotate(
                    p["sample_id"], (p["PC1"], p["PC2"]),
                    xytext=(5, 4 if s_i % 2 == 0 else -8), textcoords="offset points",
                    fontsize=7.8, color=MCK_NAVY, fontweight="bold",
                )

    vexp = landscape["pca_variance_pct"]
    ax_pca.set_xlabel(f"PC1 ({vexp[0]:.1f}% splicing variance)", fontweight="bold", color=MCK_NAVY)
    ax_pca.set_ylabel(f"PC2 ({vexp[1]:.1f}% splicing variance)", fontweight="bold", color=MCK_NAVY)
    ax_pca.set_title(
        f"A. Sample splicing manifold (PERMANOVA R² = {landscape['permanova_r2']*100:.1f}%, P = {landscape['permanova_p']:.3g})",
        loc="left", fontsize=10.2, fontweight="bold", color=MCK_NAVY,
    )
    ax_pca.legend(frameon=False, fontsize=8.2, loc="best")

    # Panel B: Variance spectrum + Beyond-Jutils QC summary
    pcs = [f"PC{i+1}" for i in range(len(vexp))]
    cum = np.cumsum(vexp)
    ax_scree.bar(np.arange(len(vexp)), vexp, width=0.52, color=MCK_NAVY, label="Individual PC %")
    ax_scree.plot(np.arange(len(vexp)), cum, color=MCK_BLUE, marker="o", linewidth=1.8, label="Cumulative %")
    for xi, v in enumerate(vexp):
        ax_scree.text(xi, v + max(vexp) * 0.04, f"{v:.1f}%", ha="center", fontsize=7.8, color=MCK_NAVY)
    ax_scree.set_xticks(np.arange(len(vexp)), pcs)
    ax_scree.set_ylabel("Explained variance (%)")
    ax_scree.set_ylim(0, max(cum[-1] * 1.16, 25))
    ax_scree.set_title("B. Variance spectrum", loc="left", fontsize=10.2, fontweight="bold", color=MCK_NAVY)
    ax_scree.legend(frameon=False, fontsize=8.0, loc="upper left")

    fig.suptitle(title, x=0.01, ha="left", fontsize=12.0, fontweight="bold", color=MCK_NAVY, y=1.02)
    fig.subplots_adjust(wspace=0.24)
    fig.savefig(png_path, dpi=220, bbox_inches="tight", facecolor="white")
    if svg_path:
        fig.savefig(svg_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def render_mckinsey_dual_psi_umap(
    dual: dict[str, Any],
    png_path: str | Path,
    svg_path: str | Path | None = None,
    title: str = "Precision-Kernel Event-PSI UMAP & Cell-Type Modality Divergence",
) -> None:
    """Render 3-panel Single-Cell / Cohort UMAP:
    Panel A: Gene Expression UMAP (colored by Cell Type, shaped by Activation State)
    Panel B: Precision-Kernel Event-PSI UMAP (colored by Cell Type, shaped by Activation State)
    Panel C: Cross-Modality Splicing-vs-Expression Divergence Score Delta_kNN per cell
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    png_path = Path(png_path)
    png_path.parent.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.8, "axes.spines.top": False, "axes.spines.right": False})

    cells = dual["cells"]
    ctypes = sorted(set(c.get("cell_type") or c.get("condition") for c in cells))
    ct_col = {ct: CONDITION_PALETTE[i % len(CONDITION_PALETTE)] for i, ct in enumerate(ctypes)}
    conds = sorted(set(c.get("condition", "") for c in cells))
    mk_map = {cond: ["o", "^", "s", "D"][i % 4] for i, cond in enumerate(conds)}

    fig, (ax1, ax2, ax3) = plt.subplots(
        1, 3, figsize=(16.2, 4.8),
        gridspec_kw={"wspace": 0.27}, dpi=220,
    )
    fig.patch.set_facecolor("white")

    for ct in ctypes:
        for cond in conds:
            sub = [c for c in cells if (c.get("cell_type") or c.get("condition")) == ct and c.get("condition") == cond]
            if not sub:
                continue
            lbl = f"{ct} · {cond}" if ct != cond else ct
            ax1.scatter(
                [c["expr_UMAP1"] for c in sub], [c["expr_UMAP2"] for c in sub],
                s=34, color=ct_col[ct], marker=mk_map[cond], alpha=0.88, edgecolors="white", linewidths=0.4, label=lbl,
            )
            ax2.scatter(
                [c["psi_UMAP1"] for c in sub], [c["psi_UMAP2"] for c in sub],
                s=34, color=ct_col[ct], marker=mk_map[cond], alpha=0.88, edgecolors="white", linewidths=0.4, label=lbl,
            )

    ax1.set_title(f"A. Expression UMAP (Sil = {dual['cell_type_silhouette_expr']:.2f})", loc="left", fontsize=10.0, fontweight="bold", color=MCK_NAVY)
    ax1.set_xlabel("Expression UMAP 1")
    ax1.set_ylabel("Expression UMAP 2")
    ax1.legend(frameon=False, fontsize=7.5, loc="best")

    ax2.set_title(f"B. Precision Event-PSI UMAP (Sil = {dual['cell_type_silhouette_psi']:.2f})", loc="left", fontsize=10.0, fontweight="bold", color=MCK_NAVY)
    ax2.set_xlabel("Precision PSI UMAP 1 (Procrustes-aligned)")
    ax2.set_ylabel("Precision PSI UMAP 2")

    div_cmap = LinearSegmentedColormap.from_list("div_mck", ["#E2E8F0", "#60A5FA", "#2251FF", MCK_CORAL])
    divs = np.array([c["modality_divergence"] for c in cells], dtype=float)
    sc3 = ax3.scatter(
        [c["psi_UMAP1"] for c in cells], [c["psi_UMAP2"] for c in cells],
        s=38, c=divs, cmap=div_cmap, vmin=0.0, vmax=1.0, edgecolors="white", linewidths=0.4,
    )
    cb = fig.colorbar(sc3, ax=ax3, fraction=0.046, pad=0.03)
    cb.set_label("Cross-modality divergence Δ_kNN", fontsize=8.0)
    ax3.set_title(f"C. Modality Divergence (Donor LISI = {dual['median_donor_lisi']:.2f})", loc="left", fontsize=10.0, fontweight="bold", color=MCK_NAVY)
    ax3.set_xlabel("Precision PSI UMAP 1")
    ax3.set_ylabel("Precision PSI UMAP 2")

    fig.suptitle(title, x=0.01, ha="left", fontsize=12.0, fontweight="bold", color=MCK_NAVY, y=1.03)
    fig.subplots_adjust(wspace=0.24)
    fig.savefig(png_path, dpi=220, bbox_inches="tight", facecolor="white")
    if svg_path:
        fig.savefig(svg_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
