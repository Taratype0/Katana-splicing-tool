"""Task 1: Threshold-free single-cell splicing contrast & dual donor/gene-block calibrated RBP discovery."""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata


@dataclass(frozen=True)
class SCMotifContrastResultV2:
    """Output of SC MotifTrace v2 de novo RBP & mature-mRNA motif discovery."""

    rbp_ranking: pd.DataFrame
    region_tests: pd.DataFrame
    event_contrasts: pd.DataFrame
    mature_consequences: pd.DataFrame
    diagnostics: dict


def compute_sc_splicing_contrast_v2(
    usage: pd.DataFrame,
    metadata: pd.DataFrame,
    event_genes: pd.Series,
    *,
    support: pd.DataFrame | None = None,
    group_col: str = "condition",
    numerator_group: str,
    denominator_group: str,
    donor_col: str = "donor_id",
    min_cells_per_group: int = 8,
    min_support_per_group: float = 10.0,
    n_cell_blocks_if_few_donors: int = 8,
    seed: int = 20261117,
) -> tuple[pd.DataFrame, pd.DataFrame | None, dict]:
    """Compute precision-weighted continuous splicing shifts (`dpsi`, `z_shift`) without hard FDR filtering.

    - When >= 5 paired donors carry both `numerator_group` and `denominator_group`, also returns a
      paired `(donor_id x event_id)` delta table for exact population-level donor sign-symmetry inference!
    - When < 5 donors exist (e.g., 2-donor Smart-seq3 or single-library CRISPRi perturbation), estimates
      between-replicate overdispersion via stratified cell-block pseudobulks and tags `inference_unit="gene_block_contrast"`.
    """
    md = metadata.copy()
    if "qc_pass" in md.columns:
        md = md.loc[md["qc_pass"].eq(True)].copy()
    if donor_col not in md.columns and "donor" in md.columns:
        md[donor_col] = md["donor"]
    if donor_col not in md.columns:
        md[donor_col] = "single_cohort"
    md["cell_id"] = md["cell_id"].astype(str)
    md = md.loc[md[group_col].astype(str).isin([str(numerator_group), str(denominator_group)])].copy()
    if md.empty:
        raise ValueError("Zero cells match numerator_group / denominator_group")

    u_df = usage.copy()
    u_df.index = u_df.index.astype(str)
    u = u_df.loc[md["cell_id"]].to_numpy(dtype=float)
    finite = np.isfinite(u)
    if support is not None:
        s_df = support.copy()
        s_df.index = s_df.index.astype(str)
        t = np.where(finite, np.maximum(s_df.loc[md["cell_id"], u_df.columns].to_numpy(dtype=float), 0.0), 0.0)
    else:
        t = finite.astype(float) * 10.0

    g_vec = md[group_col].astype(str).to_numpy()
    d_vec = md[donor_col].astype(str).to_numpy()
    is_num = g_vec == str(numerator_group)
    is_den = g_vec == str(denominator_group)

    def _pooled_stats(mask: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        f_sub = finite[mask]
        t_sub = np.where(f_sub, t[mask], 0.0)
        u_sub = np.where(f_sub, u[mask], 0.0)
        w_sub = np.divide(t_sub, 1.0 + 0.02 * t_sub, out=np.zeros_like(t_sub), where=(t_sub > 0) & f_sub)
        w_tot = w_sub.sum(axis=0)
        t_tot = t_sub.sum(axis=0)
        n_c = f_sub.sum(axis=0)
        psi_hat = np.divide((w_sub * u_sub).sum(axis=0), w_tot, out=np.full(u.shape[1], np.nan), where=w_tot > 0)
        # Binomial + finite-cell sampling variance
        p_clip = np.clip(np.nan_to_num(psi_hat, nan=0.5), 0.02, 0.98)
        se2_bin = np.divide(p_clip * (1.0 - p_clip), np.maximum(t_tot, 1.0), out=np.full(u.shape[1], np.nan), where=t_tot > 0)
        return psi_hat, se2_bin, n_c, t_tot

    psi_num, se2_num, nc_num, sup_num = _pooled_stats(is_num)
    psi_den, se2_den, nc_den, sup_den = _pooled_stats(is_den)

    # Estimate replicate/block empirical variance across donors or stratified cell blocks
    donors_num = set(d_vec[is_num])
    donors_den = set(d_vec[is_den])
    paired_donors = sorted(donors_num & donors_den)

    donor_delta_df: pd.DataFrame | None = None
    if len(paired_donors) >= 5:
        inference_unit = "paired_donor_population"
        d_deltas = []
        d_vars = []
        for d in paired_donors:
            p1, v1, c1, s1 = _pooled_stats(is_num & (d_vec == d))
            p0, v0, c0, s0 = _pooled_stats(is_den & (d_vec == d))
            ok_d = (c1 >= 2) & (c0 >= 2) & (s1 >= 3.0) & (s0 >= 3.0) & np.isfinite(p1) & np.isfinite(p0)
            d_deltas.append(np.where(ok_d, p1 - p0, np.nan))
            d_vars.append(np.where(ok_d, v1 + v0, np.nan))
        d_mat = np.stack(d_deltas, axis=0)  # D x E
        n_d = np.isfinite(d_mat).sum(axis=0)
        mean_d = np.divide(np.nansum(d_mat, axis=0), np.maximum(n_d, 1), out=np.full(u.shape[1], np.nan), where=n_d >= 3)
        var_between = np.divide(
            np.nansum((d_mat - mean_d[None, :]) ** 2, axis=0),
            np.maximum(n_d - 1, 1) * np.maximum(n_d, 1),
            out=np.full(u.shape[1], np.nan),
            where=n_d >= 3,
        )
        donor_delta_df = pd.DataFrame(d_mat, index=paired_donors, columns=u_df.columns.astype(str))
    else:
        inference_unit = "gene_block_contrast"
        rng = np.random.default_rng(int(seed))
        K = int(max(4, n_cell_blocks_if_few_donors))
        num_idx = rng.permutation(np.flatnonzero(is_num))
        den_idx = rng.permutation(np.flatnonzero(is_den))
        num_splits = np.array_split(num_idx, K)
        den_splits = np.array_split(den_idx, K)
        b_deltas = []
        for k in range(K):
            m1 = np.zeros(len(md), dtype=bool)
            m0 = np.zeros(len(md), dtype=bool)
            m1[num_splits[k]] = True
            m0[den_splits[k]] = True
            p1, _, _, s1 = _pooled_stats(m1)
            p0, _, _, s0 = _pooled_stats(m0)
            ok_b = (s1 >= 1.0) & (s0 >= 1.0) & np.isfinite(p1) & np.isfinite(p0)
            b_deltas.append(np.where(ok_b, p1 - p0, np.nan))
        b_mat = np.stack(b_deltas, axis=0)
        n_b = np.isfinite(b_mat).sum(axis=0)
        bm = np.divide(np.nansum(b_mat, axis=0), np.maximum(n_b, 1), out=np.full(u.shape[1], np.nan), where=n_b >= 2)
        var_between = np.divide(
            np.nansum((b_mat - bm[None, :]) ** 2, axis=0),
            np.maximum(n_b - 1, 1) * np.maximum(n_b, 1),
            out=np.zeros(u.shape[1], dtype=float),
            where=n_b >= 2,
        )

    eligible = (
        (nc_num >= int(min_cells_per_group))
        & (nc_den >= int(min_cells_per_group))
        & (sup_num >= float(min_support_per_group))
        & (sup_den >= float(min_support_per_group))
        & np.isfinite(psi_num)
        & np.isfinite(psi_den)
    )
    raw_dpsi = np.where(eligible, psi_num - psi_den, np.nan)
    baseline_psi = np.where(eligible, psi_den, np.nan)
    se2_tot = np.maximum(np.nan_to_num(se2_num + se2_den, nan=0.05), np.nan_to_num(var_between, nan=0.0)) + 0.0016
    se_tot = np.sqrt(se2_tot)
    z_shift = np.where(eligible, raw_dpsi / se_tot, np.nan)
    # Empirical-Bayes shrunken dPSI for downstream locus tracing
    tau2_dpsi = max(float(np.nanmean(np.maximum(raw_dpsi[eligible] ** 2 - se2_tot[eligible], 0.0))) if np.any(eligible) else 0.005, 0.001)
    shrunk_dpsi = np.where(eligible, raw_dpsi * (tau2_dpsi / (tau2_dpsi + se2_tot)), np.nan)

    genes = event_genes.reindex(u_df.columns).astype(str).to_numpy()
    contrast_df = pd.DataFrame({
        "event_id": u_df.columns.astype(str),
        "gene_id": genes,
        "eligible": eligible,
        "psi_numerator": psi_num,
        "psi_denominator": psi_den,
        "baseline_PSI": baseline_psi,
        "raw_dpsi": raw_dpsi,
        "dpsi": shrunk_dpsi,
        "se_dpsi": np.where(eligible, se_tot, np.nan),
        "z_shift": z_shift,
        "cells_numerator": nc_num.astype(int),
        "cells_denominator": nc_den.astype(int),
        "support_numerator": sup_num,
        "support_denominator": sup_den,
    })
    info = {
        "inference_unit": inference_unit,
        "paired_donors": len(paired_donors),
        "eligible_events": int(eligible.sum()),
        "eligible_genes": int(len(set(genes[eligible]))),
    }
    return contrast_df, donor_delta_df, info


def _build_nuisance_design_v2(
    cov_df: pd.DataFrame,
    baseline_psi: np.ndarray,
    eligible: np.ndarray,
) -> np.ndarray:
    """Build 5-covariate Freedman-Lane nuisance design [1, exon_GC, flank_GC, log_length, baseline_PSI, psi*(1-psi), event_type]."""
    sub = cov_df.loc[eligible]
    b_psi = np.nan_to_num(baseline_psi[eligible], nan=0.5)
    cols = [
        np.ones(len(sub), dtype=float),
        sub["exon_GC"].to_numpy(dtype=float),
        sub["flank_GC"].to_numpy(dtype=float),
        sub["log_length"].to_numpy(dtype=float),
        b_psi,
        b_psi * (1.0 - b_psi),
    ]
    etypes = sub["event_type"].astype(str).to_numpy()
    uniq_et = sorted(set(etypes))
    for et in uniq_et[1:]:
        cols.append((etypes == et).astype(float))
    X = np.column_stack(cols)
    # Guard against NaN/Inf in external sequence covariate tables by column-median / zero imputation
    for col_i in range(1, X.shape[1]):
        col = X[:, col_i]
        bad = ~np.isfinite(col)
        if np.any(bad):
            med = float(np.median(col[~bad])) if np.any(~bad) else 0.0
            X[bad, col_i] = med
    return X


def _weighted_fwl_residuals_v2(Y: np.ndarray, X: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Weighted Frisch-Waugh-Lovell residualization for 1D or 2D Y."""
    w_root = np.sqrt(np.maximum(weights, 1e-12))
    Xw = X * w_root[:, None]
    if Y.ndim == 1:
        yw = Y * w_root
        beta, _, _, _ = np.linalg.lstsq(Xw, yw, rcond=None)
        return Y - X @ beta
    Yw = Y * w_root[:, None]
    beta, _, _, _ = np.linalg.lstsq(Xw, Yw, rcond=None)
    return Y - X @ beta


def discover_sc_rbp_programs_v2(
    contrast_df: pd.DataFrame,
    rbp_tensor: np.ndarray,
    rbps: np.ndarray,
    region_names: list[str],
    covariates_df: pd.DataFrame,
    *,
    donor_delta_df: pd.DataFrame | None = None,
    mature_consequences_df: pd.DataFrame | None = None,
    rbp_expr_z: pd.Series | None = None,
    n_permutations: int = 1999,
    seed: int = 20261117,
    deconfound_correlated_motifs: bool = True,
) -> SCMotifContrastResultV2:
    """Discover candidate RBPs and region-specific motif mechanisms from single-cell splicing contrasts.

    - Dual statistical calibration:
      1. When `donor_delta_df` has >= 5 paired donors (`paired_donor_population`), computes exact/randomization
         variance-moderated donor sign-flip max-T FWER (strictly valid under `heterogeneous_population_null`!)
         with gene-block cross-event Z tie-breaking and optional ridge de-spillover of correlated bystander motifs.
      2. When < 5 donors exist (`gene_block_contrast`), runs 5-covariate Freedman-Lane orthogonalized
         Self-Normalized Rademacher Gene-Block Westfall-Young Max-T FWER across `{0..B}` draws!
    """
    elig = contrast_df["eligible"].to_numpy(dtype=bool) & np.isfinite(contrast_df["z_shift"].to_numpy(dtype=float))
    n_elig = int(elig.sum())
    if n_elig < 25:
        empty_rank = pd.DataFrame(columns=["RBP", "best_region", "direction", "stat_Z", "rbp_P_FWER", "nominal_P"])
        return SCMotifContrastResultV2(
            rbp_ranking=empty_rank,
            region_tests=pd.DataFrame(),
            event_contrasts=contrast_df,
            mature_consequences=mature_consequences_df if mature_consequences_df is not None else pd.DataFrame(),
            diagnostics={"status": "insufficient_eligible_events", "eligible_events": n_elig},
        )

    ev_sub = contrast_df.loc[elig].reset_index(drop=True)
    genes_sub = ev_sub["gene_id"].astype(str).to_numpy()
    uniq_genes, gene_codes = np.unique(genes_sub, return_inverse=True)
    n_genes = len(uniq_genes)
    gene_mult = np.bincount(gene_codes)[gene_codes].astype(float)
    ev_weights = 1.0 / gene_mult
    ev_weights = ev_weights / ev_weights.sum()

    cov_aligned = covariates_df.set_index("event_id").reindex(contrast_df["event_id"].astype(str)).reset_index()
    X_nuis = _build_nuisance_design_v2(cov_aligned, contrast_df["baseline_PSI"].to_numpy(dtype=float), elig)

    # Winsorize extreme z_shift lightly at +-8 to prevent a single ultra-deep outlier gene from dominating
    z_raw = np.clip(ev_sub["z_shift"].to_numpy(dtype=float), -8.0, 8.0)
    z_res = _weighted_fwl_residuals_v2(z_raw, X_nuis, ev_weights)

    n_reg = len(region_names)
    n_rbp = len(rbps)
    F_sub = rbp_tensor[elig]  # n_elig x n_reg x n_rbp
    F_flat = F_sub.reshape(n_elig, n_reg * n_rbp)
    F_res = _weighted_fwl_residuals_v2(F_flat, X_nuis, ev_weights)
    # Standardize residualized motif features with gene-equal weights
    f_var = np.sum((F_res * F_res) * ev_weights[:, None], axis=0)
    f_sd = np.sqrt(np.maximum(f_var, 1e-12))
    F_std = np.where(f_var[None, :] > 1e-11, F_res / f_sd[None, :], 0.0)

    # Gene-block self-normalized scores S_{g, c} for c in (n_reg * n_rbp)
    contrib = (ev_weights * z_res)[:, None] * F_std  # n_elig x C
    S_gene = np.zeros((n_genes, n_reg * n_rbp), dtype=float)
    np.add.at(S_gene, gene_codes, contrib)
    denom_gene = np.sqrt(np.sum(S_gene * S_gene, axis=0))
    valid_col = (f_var > 1e-11) & (denom_gene > 1e-12)
    Z_gene_obs = np.divide(S_gene.sum(axis=0), denom_gene, out=np.zeros(n_reg * n_rbp), where=valid_col)

    use_donor_pop = donor_delta_df is not None and len(donor_delta_df) >= 5
    rng = np.random.default_rng(int(seed))

    if use_donor_pop:
        inference_mode = "paired_donor_moderated_maxT"
        d_sub = donor_delta_df.reindex(columns=contrast_df.loc[elig, "event_id"].astype(str)).to_numpy(dtype=float)  # D x n_elig
        n_donors = d_sub.shape[0]
        # Compute nuisance-adjusted slope per donor x (region, RBP)
        d_filled = np.where(np.isfinite(d_sub), d_sub, 0.0)
        d_res = _weighted_fwl_residuals_v2(d_filled.T, X_nuis, ev_weights).T  # D x n_elig

        if deconfound_correlated_motifs and n_reg * n_rbp > 1:
            # Ridge-decorrelated donor slopes suppress collinear bystander spillover while preserving sign symmetry!
            G_cov = (F_std.T * ev_weights) @ F_std
            ridge_lam = 0.35 * float(np.median(np.diag(G_cov)[valid_col])) if np.any(valid_col) else 0.35
            rhs_d = (d_res * ev_weights[None, :]) @ F_std  # D x C
            B_donor = np.linalg.solve(G_cov + max(ridge_lam, 0.05) * np.eye(G_cov.shape[0]), rhs_d.T).T
            B_donor[:, ~valid_col] = 0.0
        else:
            B_donor = (d_res * ev_weights[None, :]) @ F_std

        if n_donors <= 14:
            signs = np.array(list(product([-1.0, 1.0], repeat=n_donors)), dtype=float)
            exact_signs = True
        else:
            signs = rng.choice([-1.0, 1.0], size=(int(n_permutations), n_donors))
            exact_signs = False

        # Variance-moderated donor t-denominator (invariant to donor sign flips b_d -> +-b_d!)
        # Note: sum_d b_d^2 IS strictly invariant under arbitrary sign flips s_d in {-1, +1}!
        ss_b = np.sum(B_donor * B_donor, axis=0)
        pos_ss = ss_b[valid_col & (ss_b > 0)]
        s0_sq = float(np.quantile(pos_ss, 0.25)) * 0.15 if len(pos_ss) else 1e-8
        mod_denom = np.sqrt(ss_b + s0_sq)
        Z_obs = np.divide(B_donor.sum(axis=0), mod_denom, out=np.zeros(n_reg * n_rbp), where=valid_col & (mod_denom > 0))
        null_Z = np.divide(signs @ B_donor, mod_denom[None, :], out=np.zeros((len(signs), n_reg * n_rbp)), where=valid_col[None, :] & (mod_denom[None, :] > 0))

        # Per-RBP max(|Z|) across its regions for exact subset-compatible familywise calibration
        Z_obs_2d = Z_obs.reshape(n_reg, n_rbp)
        null_3d = null_Z.reshape(len(signs), n_reg, n_rbp)
        rbp_obs_max = np.max(np.abs(Z_obs_2d), axis=0)
        rbp_null_max = np.max(np.abs(null_3d), axis=1)  # n_signs x n_rbp
        global_null_max = np.max(rbp_null_max, axis=1) if n_rbp > 0 else np.zeros(len(signs))

        corr_add = 0.0 if exact_signs else 1.0
        denom_p = float(len(signs)) + corr_add
        nom_p = np.ones(n_reg * n_rbp, dtype=float)
        reg_fwer = np.ones(n_reg * n_rbp, dtype=float)
        for c in np.flatnonzero(valid_col):
            thr = abs(Z_obs[c]) - 1e-12
            nom_p[c] = (corr_add + float(np.sum(np.abs(null_Z[:, c]) >= thr))) / denom_p
            reg_fwer[c] = (corr_add + float(np.sum(global_null_max >= thr))) / denom_p

        rbp_nom_p = np.ones(n_rbp, dtype=float)
        rbp_fwer_p = np.ones(n_rbp, dtype=float)
        for r_i in range(n_rbp):
            thr_r = rbp_obs_max[r_i] - 1e-12
            if thr_r > 0:
                rbp_nom_p[r_i] = (corr_add + float(np.sum(rbp_null_max[:, r_i] >= thr_r))) / denom_p
                rbp_fwer_p[r_i] = (corr_add + float(np.sum(global_null_max >= thr_r))) / denom_p

    else:
        inference_mode = "gene_block_rademacher_maxT"
        B = int(n_permutations)
        signs = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, n_genes))
        S_norm = np.divide(S_gene, denom_gene[None, :], out=np.zeros_like(S_gene), where=valid_col[None, :])
        Z_obs = Z_gene_obs.copy()
        null_Z = (signs @ np.ascontiguousarray(S_norm, dtype=np.float32)).astype(np.float64)

        Z_obs_2d = Z_obs.reshape(n_reg, n_rbp)
        null_3d = null_Z.reshape(B, n_reg, n_rbp)
        rbp_obs_max = np.max(np.abs(Z_obs_2d), axis=0)
        rbp_null_max = np.max(np.abs(null_3d), axis=1)
        global_null_max = np.max(rbp_null_max, axis=1) if n_rbp > 0 else np.zeros(B)

        nom_p = (1.0 + np.sum(np.abs(null_Z) >= np.abs(Z_obs)[None, :] - 1e-12, axis=0)) / (B + 1.0)
        reg_fwer = (1.0 + np.sum(global_null_max[:, None] >= np.abs(Z_obs)[None, :] - 1e-12, axis=0)) / (B + 1.0)
        nom_p[~valid_col] = 1.0
        reg_fwer[~valid_col] = 1.0

        rbp_nom_p = (1.0 + np.sum(rbp_null_max >= rbp_obs_max[None, :] - 1e-12, axis=0)) / (B + 1.0)
        rbp_fwer_p = (1.0 + np.sum(global_null_max[:, None] >= rbp_obs_max[None, :] - 1e-12, axis=0)) / (B + 1.0)

    # Assemble region-level and RBP-level tables
    Z_gene_2d = Z_gene_obs.reshape(n_reg, n_rbp)
    Z_obs_2d = Z_obs.reshape(n_reg, n_rbp)
    nom_2d = nom_p.reshape(n_reg, n_rbp)
    fwer_2d = reg_fwer.reshape(n_reg, n_rbp)

    reg_rows: list[dict] = []
    for r_idx, reg_name in enumerate(region_names):
        for m_idx, rbp_name in enumerate(rbps):
            zv = float(Z_obs_2d[r_idx, m_idx])
            reg_rows.append({
                "RBP": str(rbp_name),
                "region": str(reg_name),
                "layer": "affected_sequence" if reg_name == "mature_consequence_delta" else "regulatory_pre_mRNA",
                "direction": "positive" if zv >= 0 else "negative",
                "stat_Z": zv,
                "gene_block_Z": float(Z_gene_2d[r_idx, m_idx]),
                "nominal_P": float(nom_2d[r_idx, m_idx]),
                "region_P_FWER": float(fwer_2d[r_idx, m_idx]),
            })
    reg_df = pd.DataFrame(reg_rows)

    rbp_rows: list[dict] = []
    for m_idx, rbp_name in enumerate(rbps):
        best_r = int(np.argmax(np.abs(Z_obs_2d[:, m_idx])))
        zv = float(Z_obs_2d[best_r, m_idx])
        gzv = float(Z_gene_2d[best_r, m_idx])
        mature_z = float(Z_obs_2d[0, m_idx]) if region_names[0] == "mature_consequence_delta" else np.nan
        mature_p = float(fwer_2d[0, m_idx]) if region_names[0] == "mature_consequence_delta" else np.nan
        # Continuous tie-free evidence score: analytical tail survival + gene-block Z tie-breaker
        seq_score = float(-norm.logsf(abs(zv)) / np.log(10.0) + 0.02 * abs(gzv))
        expr_z = float(rbp_expr_z.get(str(rbp_name), np.nan)) if rbp_expr_z is not None else np.nan
        concordant_score = seq_score + (0.35 * min(abs(expr_z), 8.0) if np.isfinite(expr_z) else 0.0)
        rbp_rows.append({
            "RBP": str(rbp_name),
            "best_region": str(region_names[best_r]),
            "best_layer": "affected_sequence" if region_names[best_r] == "mature_consequence_delta" else "regulatory_pre_mRNA",
            "direction": "inclusion_promoted" if zv >= 0 else "skipping_promoted",
            "stat_Z": zv,
            "abs_Z": abs(zv),
            "gene_block_Z": gzv,
            "mature_consequence_Z": mature_z,
            "mature_consequence_P_FWER": mature_p,
            "nominal_P": float(rbp_nom_p[m_idx]),
            "rbp_P_FWER": float(rbp_fwer_p[m_idx]),
            "sequence_evidence_score": seq_score,
            "rbp_expr_z": expr_z,
            "concordant_activity_score": concordant_score,
        })

    rbp_df = pd.DataFrame(rbp_rows).sort_values(
        by=["rbp_P_FWER", "sequence_evidence_score", "RBP"],
        ascending=[True, False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    rbp_df["sequence_rank"] = np.arange(1, len(rbp_df) + 1, dtype=int)

    return SCMotifContrastResultV2(
        rbp_ranking=rbp_df,
        region_tests=reg_df,
        event_contrasts=contrast_df,
        mature_consequences=mature_consequences_df if mature_consequences_df is not None else pd.DataFrame(),
        diagnostics={
            "status": "complete",
            "inference_mode": inference_mode,
            "eligible_events": n_elig,
            "eligible_genes": n_genes,
            "evaluated_rbps": int(n_rbp),
            "significant_rbps_fwer05": int((rbp_df["rbp_P_FWER"] <= 0.05).sum()),
        },
    )


def compute_sc_stratified_meta_contrast_v2(
    usage: pd.DataFrame,
    metadata: pd.DataFrame,
    event_genes: pd.Series,
    *,
    support: pd.DataFrame | None = None,
    group_col: str = "condition",
    numerator_group: str,
    denominator_group: str,
    donor_col: str = "donor_id",
    min_cells_per_donor_group: int = 4,
    min_support_per_donor_group: float = 6.0,
    seed: int = 20261117,
) -> tuple[pd.DataFrame, dict[str, pd.DataFrame], dict]:
    """2-4 Donor Stratified Fixed-Effect Meta-Analysis Contrast Engine.

    Eliminates between-donor baseline batch offsets by estimating precision-weighted `(dpsi_d, se_d)`
    STRICTLY WITHIN each donor first, then combining across donors via inverse-variance meta-analysis:
        w_{de} = 1 / se_{de}^2,   dpsi_meta_e = sum_d(w_{de} dpsi_{de}) / sum_d(w_{de}),
        z_meta_e = sum_d(w_{de} dpsi_{de}) / sqrt(sum_d(w_{de})).
    Also returns the per-donor contrast dictionary so downstream discovery can enforce directional sign concordance!
    """
    md = metadata.copy()
    if "qc_pass" in md.columns:
        md = md.loc[md["qc_pass"].eq(True)].copy()
    if donor_col not in md.columns and "donor" in md.columns:
        md[donor_col] = md["donor"]
    donors = sorted(md[donor_col].astype(str).unique())

    per_donor: dict[str, pd.DataFrame] = {}
    for d_idx, d_id in enumerate(donors):
        md_d = md.loc[md[donor_col].astype(str).eq(d_id)].copy()
        if md_d[group_col].astype(str).nunique() < 2:
            continue
        c_d, _, _ = compute_sc_splicing_contrast_v2(
            usage,
            md_d,
            event_genes,
            support=support,
            group_col=group_col,
            numerator_group=numerator_group,
            denominator_group=denominator_group,
            donor_col=donor_col,
            min_cells_per_group=min_cells_per_donor_group,
            min_support_per_group=min_support_per_donor_group,
            seed=int(seed) + 31 * d_idx,
        )
        per_donor[d_id] = c_d

    if not per_donor:
        raise ValueError("Zero donors have both numerator and denominator cells")

    ev_ids = usage.columns.astype(str).to_numpy()
    genes = event_genes.reindex(usage.columns).astype(str).to_numpy()
    n_ev = len(ev_ids)
    w_sum = np.zeros(n_ev, dtype=float)
    wd_sum = np.zeros(n_ev, dtype=float)
    wb_sum = np.zeros(n_ev, dtype=float)
    d_cnt = np.zeros(n_ev, dtype=int)
    c_num_tot = np.zeros(n_ev, dtype=int)
    c_den_tot = np.zeros(n_ev, dtype=int)
    s_num_tot = np.zeros(n_ev, dtype=float)
    s_den_tot = np.zeros(n_ev, dtype=float)

    pos_d = np.zeros(n_ev, dtype=int)
    neg_d = np.zeros(n_ev, dtype=int)
    for c_d in per_donor.values():
        el = c_d["eligible"].to_numpy(dtype=bool)
        se = c_d["se_dpsi"].to_numpy(dtype=float)
        dp = c_d["raw_dpsi"].to_numpy(dtype=float)
        bp = c_d["baseline_PSI"].to_numpy(dtype=float)
        ok = el & np.isfinite(se) & (se > 0) & np.isfinite(dp)
        w = np.where(ok, 1.0 / (se * se), 0.0)
        w_sum += w
        wd_sum += np.where(ok, w * dp, 0.0)
        wb_sum += np.where(ok & np.isfinite(bp), w * bp, 0.0)
        d_cnt += ok.astype(int)
        pos_d += (ok & (dp > 0)).astype(int)
        neg_d += (ok & (dp < 0)).astype(int)
        c_num_tot += c_d["cells_numerator"].to_numpy(dtype=int)
        c_den_tot += c_d["cells_denominator"].to_numpy(dtype=int)
        s_num_tot += c_d["support_numerator"].to_numpy(dtype=float)
        s_den_tot += c_d["support_denominator"].to_numpy(dtype=float)

    elig_meta = w_sum > 0.0
    sign_agree = np.divide(np.maximum(pos_d, neg_d), np.maximum(d_cnt, 1), out=np.zeros(n_ev, dtype=float), where=d_cnt > 0)
    raw_meta = np.divide(wd_sum, w_sum, out=np.full(n_ev, np.nan), where=elig_meta)
    base_meta = np.divide(wb_sum, w_sum, out=np.full(n_ev, 0.5), where=elig_meta)
    se_meta = np.divide(1.0, np.sqrt(np.maximum(w_sum, 1e-12)), out=np.full(n_ev, np.nan), where=elig_meta)
    z_meta = np.divide(wd_sum, np.sqrt(np.maximum(w_sum, 1e-12)), out=np.full(n_ev, np.nan), where=elig_meta)
    se2_m = np.where(elig_meta, se_meta * se_meta, 0.05)
    tau2_m = max(float(np.nanmean(np.maximum(raw_meta[elig_meta] ** 2 - se2_m[elig_meta], 0.0))) if np.any(elig_meta) else 0.005, 0.001)
    shrunk_meta = np.where(elig_meta, raw_meta * (tau2_m / (tau2_m + se2_m)), np.nan)

    meta_df = pd.DataFrame({
        "event_id": ev_ids,
        "gene_id": genes,
        "eligible": elig_meta,
        "contributing_donors": d_cnt,
        "donor_sign_agreement": sign_agree,
        "baseline_PSI": base_meta,
        "raw_dpsi": raw_meta,
        "dpsi": shrunk_meta,
        "se_dpsi": se_meta,
        "z_shift": z_meta,
        "cells_numerator": c_num_tot,
        "cells_denominator": c_den_tot,
        "support_numerator": s_num_tot,
        "support_denominator": s_den_tot,
    })
    return meta_df, per_donor, {
        "inference_unit": "stratified_multi_donor_meta",
        "n_donors_meta": len(per_donor),
        "donors_evaluated": list(per_donor.keys()),
        "eligible_events": int(elig_meta.sum()),
        "multi_donor_events": int((d_cnt >= 2).sum()),
    }
