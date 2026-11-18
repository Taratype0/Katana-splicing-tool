"""Task 2 (Conditional Target Event Ranking) & Task 3 (Precision-Weighted Single-Cell Continuous PSI Prediction)."""
from __future__ import annotations

import inspect
import numpy as np
import pandas as pd
from scipy.sparse.linalg import LinearOperator, cg


def estimate_sc_psi_eb_baseline_v2(
    psi_train: np.ndarray,
    support_train: np.ndarray | None,
    *,
    prior_kappa: float = 3.0,
    prior_bio_var: float = 0.025,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute Empirical-Bayes precision-weighted baseline PSI mu_e and observation precision weights w_{ce}
    strictly on training cells (`psi_train`, `support_train`).

    Zero query-cell support or query-cell PSI is ever accessed!
    """
    y = np.asarray(psi_train, dtype=float)
    finite = np.isfinite(y)
    if support_train is None:
        w0 = finite.astype(float)
    else:
        t = np.where(finite, np.maximum(np.asarray(support_train, dtype=float), 0.0), 0.0)
        w0 = np.divide(t, 1.0 + float(prior_bio_var) * t, out=np.zeros_like(t), where=(t > 0.0) & finite)

    w_sum = w0.sum(axis=0)
    y_safe = np.where(finite, y, 0.0)
    global_mu = float((w0 * y_safe).sum() / max(w0.sum(), 1e-12)) if w0.sum() > 0 else 0.5

    # Shrink low-support training events smoothly toward global_mu to reduce small-sample binomial variance
    kappa = max(float(prior_kappa), 0.0)
    mu_eb = np.divide((w0 * y_safe).sum(axis=0) + kappa * global_mu, w_sum + kappa, out=np.full(y.shape[1], global_mu), where=(w_sum + kappa) > 0)

    # Second-pass event-specific biological dispersion estimate for observation weights
    if support_train is not None:
        diff = np.where(finite, y - mu_eb[None, :], 0.0)
        w_var = np.divide((w0 * diff * diff).sum(axis=0), np.maximum(w_sum, 1e-12))
        inv_t = np.divide(
            np.where((t > 0) & finite, w0 / np.maximum(t, 1.0), 0.0).sum(axis=0),
            np.maximum(w_sum, 1e-12),
        )
        binom_v = np.clip(mu_eb * (1.0 - mu_eb), 0.01, 0.25) * inv_t
        sig2_bio = np.clip(w_var - binom_v, 0.002, 0.15)
        w_final = np.divide(t, 1.0 + t * sig2_bio[None, :], out=np.zeros_like(t), where=(t > 0.0) & finite)
        w_sum2 = w_final.sum(axis=0)
        mu_eb = np.divide(
            (w_final * y_safe).sum(axis=0) + kappa * global_mu,
            w_sum2 + kappa,
            out=np.full(y.shape[1], global_mu),
            where=(w_sum2 + kappa) > 0,
        )
    else:
        sig2_bio = np.full(y.shape[1], float(prior_bio_var))
        w_final = w0

    return mu_eb, w_final, sig2_bio


def predict_sc_psi_precision_v2(
    expression: np.ndarray,
    psi: np.ndarray,
    training_mask: np.ndarray,
    features: np.ndarray,
    genes: np.ndarray,
    alpha: float,
    *,
    support_train_full: np.ndarray | None = None,
    prior_kappa: float = 0.0,
    unbiased_cell_intercept: bool = False,
) -> tuple[np.ndarray, dict]:
    """Inductive precision-weighted factored Gaussian MAP / Preconditioned-CG positional RBP x motif predictor.

    - When `unbiased_cell_intercept=False` (default when `Cov(support, PSI) == 0`), uses precision-weighted EB
      intercept AND precision-weighted GLS slope.
    - When `unbiased_cell_intercept=True` (for cell-equal evaluation under transcription-coupled splicing where
      high-depth cells differ slightly in mean PSI from low-depth cells), keeps the unbiased cell-equal training
      mean `E[y_ce]` as the intercept while using `support` precision weights strictly in the GLS RBP x motif slope!
    """
    x = np.asarray(expression, dtype=float)
    y = np.asarray(psi, dtype=float)
    fit = np.asarray(training_mask, dtype=bool)
    f = np.asarray(features, dtype=float)
    g = np.asarray(genes, dtype=str)

    if x.ndim != 2 or y.ndim != 2 or len(x) != len(y) or fit.shape != (len(x),):
        raise ValueError("Aligned expression, psi, and training_mask required")
    if f.ndim != 3 or f.shape[0] != y.shape[1] or f.shape[2] != x.shape[1] or g.shape != (y.shape[1],):
        raise ValueError("Aligned (E x P x R) motif feature tensor and gene array required")
    if not fit.any() or alpha < 0.0:
        raise ValueError("Non-empty training partition and alpha >= 0 required")

    sup_fit = None
    if support_train_full is not None:
        # Strictly slice training rows only so query support is NEVER touched!
        sup_fit = np.asarray(support_train_full, dtype=float)[fit]

    means_eb, w_fit_mat, sig2_bio = estimate_sc_psi_eb_baseline_v2(
        y[fit], sup_fit, prior_kappa=prior_kappa
    )
    if unbiased_cell_intercept:
        fin_fit = np.isfinite(y[fit])
        cnt_fit = fin_fit.sum(axis=0)
        unw_means = np.divide(np.where(fin_fit, y[fit], 0.0).sum(axis=0), cnt_fit, out=means_eb.copy(), where=cnt_fit > 0)
        if prior_kappa > 0.0:
            g_mu = float(np.mean(unw_means[cnt_fit > 0])) if np.any(cnt_fit > 0) else 0.5
            means = (cnt_fit * unw_means + prior_kappa * g_mu) / np.maximum(cnt_fit + prior_kappa, 1e-12)
        else:
            means = unw_means
    else:
        means = means_eb

    x_fit = x[fit]
    sd = x_fit.std(axis=0)
    sd[sd < 1e-8] = 1.0

    # Center expression on the matching intercept weighting so E[pred_query] == means when E[x_query] == E[x_fit]
    if unbiased_cell_intercept:
        fin_fit = np.isfinite(y[fit]).astype(float)
        f_sum = fin_fit.sum(axis=0)
        obs_x = np.where(
            (f_sum > 0.0)[:, None],
            (fin_fit.T @ x_fit) / np.maximum(f_sum[:, None], 1.0),
            x_fit.mean(axis=0)[None, :],
        )
    else:
        w_fit_sum = w_fit_mat.sum(axis=0)
        obs_x = np.where(
            (w_fit_sum > 0.0)[:, None],
            (w_fit_mat.T @ x_fit) / np.maximum(w_fit_sum[:, None], 1e-12),
            x_fit.mean(axis=0)[None, :],
        )

    n_ev, n_pos, n_rbp = f.shape
    dimension = n_pos * n_rbp
    beta = np.zeros(dimension, dtype=float)
    iterations = 0
    rel_res = 0.0

    fit_indices = np.flatnonzero(fit)
    ci_local, ei = np.where(np.isfinite(y[fit]) & (w_fit_mat > 0.0))
    if alpha > 0.0 and len(ci_local) > 0 and np.any(np.abs(f) > 1e-12):
        ci_global = fit_indices[ci_local]
        z = (x[ci_global] - obs_x[ei]) / sd[None, :]

        # Normalize observation precision weights so each gene contributes equal total weight 1 / N_active_genes
        obs_prec = w_fit_mat[ci_local, ei]
        ev_prec_tot = np.bincount(ei, weights=obs_prec, minlength=n_ev)
        uniq_g, g_counts = np.unique(g[ev_prec_tot > 0], return_counts=True)
        g_mult = dict(zip(uniq_g, g_counts))
        w = np.array(
            [obs_prec[k] / (max(ev_prec_tot[ei[k]], 1e-12) * g_mult.get(g[ei[k]], 1)) for k in range(len(ei))],
            dtype=float,
        )
        w /= max(w.sum(), 1e-12)

        def transpose(val: np.ndarray) -> np.ndarray:
            grouped = np.zeros((n_ev, n_rbp), dtype=float)
            np.add.at(grouped, ei, val[:, None] * z)
            return np.einsum("epr,er->pr", f, grouped).ravel()

        def matvec(v: np.ndarray) -> np.ndarray:
            wt = np.einsum("epr,pr->er", f, v.reshape(n_pos, n_rbp))
            fwd = np.einsum("nr,nr->n", z, wt[ei])
            return transpose(w * fwd) + float(alpha) * v

        rhs = transpose(w * (y[ci_global, ei] - means[ei]))
        if np.linalg.norm(rhs) > 1e-14:
            op = LinearOperator((dimension, dimension), matvec=matvec, dtype=float)
            moments = np.zeros((n_ev, n_rbp), dtype=float)
            np.add.at(moments, ei, w[:, None] * z * z)
            diag = np.einsum("epr,er->pr", f * f, moments).ravel() + float(alpha)
            diag = np.maximum(diag, 1e-12)
            precond = LinearOperator((dimension, dimension), matvec=lambda v: v / diag, dtype=float)

            def _cb(_: np.ndarray) -> None:
                nonlocal iterations
                iterations += 1

            tol_kw = {"rtol" if "rtol" in inspect.signature(cg).parameters else "tol": 1e-12}
            beta, info = cg(op, rhs, M=precond, atol=1e-14, maxiter=3000, callback=_cb, **tol_kw)
            res_norm = float(np.linalg.norm(op @ beta - rhs))
            rel_res = res_norm / max(float(np.linalg.norm(rhs)), 1e-12)
            if info != 0:
                raise RuntimeError(f"Preconditioned CG did not converge: info={info}, rel_res={rel_res}")

    weights_er = np.einsum("epr,pr->er", f, beta.reshape(n_pos, n_rbp))
    pred = means[None, :] + np.einsum("cer,er->ce", (x[:, None, :] - obs_x[None, :, :]) / sd[None, None, :], weights_er)
    return np.clip(pred, 0.0, 1.0), {
        "coefficient": beta,
        "event_mean": means,
        "sig2_bio": sig2_bio,
        "expression_scale": sd,
        "observed_expression_mean": obs_x,
        "iterations": iterations,
        "relative_residual": rel_res,
    }


def rank_rbp_target_events_v2(
    query_rbp: str,
    event_ids: np.ndarray,
    genes: np.ndarray,
    rbp_tensor: np.ndarray,
    rbps: np.ndarray,
    region_names: list[str],
    *,
    control_excess_bio_var: np.ndarray | None = None,
    region_weights: np.ndarray | None = None,
    spliceability_blend: float = 0.5,
) -> pd.DataFrame:
    """Task 2: Rank candidate target events and predicted consequence directions for a queried RBP.

    Explicitly separates:
      - `motif_susceptibility`: sequence-only RBP motif strength across mature (`affected_sequence`) & cis-intronic regions
      - `control_spliceability`: optional control-population Excess Biological Variance prior
      - `combined_score`: transparent convex rank/z blend (never claiming motif gain without reporting sequence-only score!)
    """
    r_list = [str(r) for r in rbps]
    if str(query_rbp) not in r_list:
        return pd.DataFrame({
            "event_id": event_ids.astype(str),
            "gene_id": genes.astype(str),
            "query_RBP": str(query_rbp),
            "rbp_in_library": False,
            "motif_susceptibility": 0.0,
            "mature_motif_delta": 0.0,
            "control_spliceability": np.zeros(len(event_ids)) if control_excess_bio_var is None else np.asarray(control_excess_bio_var, dtype=float),
            "combined_score": 0.0,
        })

    r_idx = r_list.index(str(query_rbp))
    mat_er = rbp_tensor[:, :, r_idx]  # E x n_reg
    if region_weights is not None:
        rw = np.asarray(region_weights, dtype=float)
    else:
        # Emphasize mature consequence + proximal flanks slightly over distal flanks
        rw = np.array([1.0 if "distal" not in rn else 0.6 for rn in region_names], dtype=float)
    rw = rw / max(np.linalg.norm(rw), 1e-12)

    motif_susc = np.sqrt(np.sum((mat_er * rw[None, :]) ** 2, axis=1))
    mature_delta = mat_er[:, 0] if (len(region_names) > 0 and region_names[0] == "mature_consequence_delta") else np.zeros(len(event_ids))

    def _zscore(v: np.ndarray) -> np.ndarray:
        s = float(np.std(v))
        return (v - float(np.mean(v))) / s if s > 1e-12 else np.zeros_like(v)

    z_m = _zscore(motif_susc)
    if control_excess_bio_var is not None:
        sp = np.maximum(np.asarray(control_excess_bio_var, dtype=float), 0.0)
        z_s = _zscore(np.sqrt(sp))
        lam = float(np.clip(spliceability_blend, 0.0, 1.0))
        comb = (1.0 - lam) * z_m + lam * z_s
    else:
        sp = np.zeros(len(event_ids), dtype=float)
        comb = z_m

    df = pd.DataFrame({
        "event_id": event_ids.astype(str),
        "gene_id": genes.astype(str),
        "query_RBP": str(query_rbp),
        "rbp_in_library": True,
        "motif_susceptibility": motif_susc,
        "mature_motif_delta": mature_delta,
        "predicted_mature_consequence": np.where(
            mature_delta > 0.1,
            "inclusion_carries_motif",
            np.where(mature_delta < -0.1, "exclusion_junction_carries_motif", "intronic_or_neutral"),
        ),
        "control_spliceability": sp,
        "combined_score": comb,
    })
    return df.sort_values(by=["combined_score", "motif_susceptibility", "event_id"], ascending=[False, False, True], kind="mergesort").reset_index(drop=True)


def gene_equal_psi_mse_v2(pred: np.ndarray, truth: np.ndarray, genes: np.ndarray) -> tuple[float, pd.Series]:
    """Compute event-mean then gene-equal MSE on observed query labels."""
    n = np.isfinite(truth).sum(axis=0)
    ss = np.nansum((pred - truth) ** 2, axis=0)
    ev_mse = np.divide(ss, n, out=np.full(len(n), np.nan), where=n > 0)
    gm = pd.DataFrame({"gene": np.asarray(genes, dtype=str), "mse": ev_mse}).groupby("gene")["mse"].mean()
    return float(gm.mean()), gm
