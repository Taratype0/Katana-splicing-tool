"""Evaluation endpoints, host-expression orthogonalization, pathway enrichment, and inductive usage prediction."""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score, average_precision_score, normalized_mutual_info_score

from katana_sc_splicenet.network import SCCoSpliceNetworkV2, precision_group_profiles_v2, precision_pairwise_correlation_v2
from katana_sc_splicenet.qc_and_selection import estimate_sc_observation_weights_v2


def orthogonalize_sc_host_expression_v2(
    train_usage: np.ndarray,
    train_host_log_expr: np.ndarray,
    test_usage: np.ndarray,
    test_host_log_expr: np.ndarray,
    *,
    min_train_cells: int = 10,
) -> tuple[np.ndarray, np.ndarray]:
    """Fit per-event usage ~ host_log_expr linear projection on training cells ONLY and residualize held-out cells."""
    tu = np.asarray(train_usage, dtype=float)
    th = np.asarray(train_host_log_expr, dtype=float)
    qu = np.asarray(test_usage, dtype=float)
    qh = np.asarray(test_host_log_expr, dtype=float)

    tr_res = tu.copy()
    te_res = qu.copy()
    for j in range(tu.shape[1]):
        ok = np.isfinite(tu[:, j]) & np.isfinite(th[:, j])
        if int(ok.sum()) < min_train_cells:
            continue
        x = th[ok, j]
        y = tu[ok, j]
        vx = float(np.var(x))
        if vx <= 1e-10:
            continue
        mx = float(np.mean(x))
        my = float(np.mean(y))
        slope = float(np.mean((x - mx) * (y - my)) / vx)
        ok_tr = np.isfinite(tu[:, j]) & np.isfinite(th[:, j])
        tr_res[ok_tr, j] = tu[ok_tr, j] - slope * (th[ok_tr, j] - mx)
        ok_te = np.isfinite(qu[:, j]) & np.isfinite(qh[:, j])
        te_res[ok_te, j] = qu[ok_te, j] - slope * (qh[ok_te, j] - mx)
    return tr_res, te_res


def evaluate_heldout_edges_and_modules_v2(
    net: SCCoSpliceNetworkV2,
    heldout_usage: pd.DataFrame,
    heldout_metadata: pd.DataFrame,
    *,
    heldout_support: pd.DataFrame | None = None,
    train_usage: pd.DataFrame | None = None,
    train_host_log_expr: pd.DataFrame | None = None,
    heldout_host_log_expr: pd.DataFrame | None = None,
    eval_level: str = "cell",
    min_pairs: int = 5,
) -> dict:
    """Evaluate primary co-splicing network and module endpoints on unshrunk held-out observed_usage."""
    event_ids = net.modules["event_id"].astype(str).tolist()
    genes = net.modules["gene_id"].astype(str).to_numpy()
    u_te = heldout_usage.reindex(columns=event_ids).to_numpy(dtype=float)

    if eval_level == "donor":
        donors = sorted(heldout_metadata["donor_id"].astype(str).unique())
        w_te = np.isfinite(u_te).astype(float)
        if heldout_support is not None:
            s_te = heldout_support.reindex(index=heldout_usage.index, columns=event_ids).to_numpy(dtype=float)
            w_te = estimate_sc_observation_weights_v2(u_te, s_te, net.modules["excess_bio_var"].to_numpy(dtype=float))
        prof_te, wt_te = precision_group_profiles_v2(
            u_te, w_te, heldout_metadata["donor_id"].astype(str).to_numpy(), donors, min_cells=3
        )
        # Evaluate held-out correlation on unshrunk observed profiles (unweighted Pearson across donors)
        r_te, _, _ = precision_pairwise_correlation_v2(prof_te, None, min_pairs=min_pairs, min_eff_n=float(min_pairs))
    else:
        # Unweighted observed_usage correlation across held-out cells so all methods share the exact same target metric
        r_te, _, _ = precision_pairwise_correlation_v2(u_te, None, min_pairs=min_pairs, min_eff_n=float(min_pairs))

    idx_map = {eid: i for i, eid in enumerate(event_ids)}
    heldout_edge_r = np.nan
    evaluable_edges = 0
    if not net.edges.empty:
        si = np.array([idx_map[s] for s in net.edges["source_event"].astype(str)], dtype=int)
        ti = np.array([idx_map[t] for t in net.edges["target_event"].astype(str)], dtype=int)
        vals = r_te[si, ti]
        fin = np.isfinite(vals)
        evaluable_edges = int(fin.sum())
        if evaluable_edges > 0:
            heldout_edge_r = float(np.mean(vals[fin]))

    def _module_coherence_from_corr(corr_mat: np.ndarray) -> tuple[float, float, float, int]:
        mods = net.modules["module"].astype(str).to_numpy()
        active_mods = sorted(set(mods) - {"M0"})
        if not active_mods:
            return np.nan, np.nan, np.nan, 0
        ia, ib = np.triu_indices(len(event_ids), k=1)
        cross_g = (genes[ia] != genes[ib]) & np.isfinite(corr_mat[ia, ib])
        same_m = (mods[ia] == mods[ib]) & (mods[ia] != "M0") & cross_g
        diff_m = (mods[ia] != mods[ib]) & ((mods[ia] != "M0") | (mods[ib] != "M0")) & cross_g
        w_in = float(np.mean(corr_mat[ia[same_m], ib[same_m]])) if np.any(same_m) else np.nan
        w_out = float(np.mean(corr_mat[ia[diff_m], ib[diff_m]])) if np.any(diff_m) else np.nan
        gap = (w_in - w_out) if (np.isfinite(w_in) and np.isfinite(w_out)) else np.nan
        return w_in, w_out, gap, len(active_mods)

    within_r, outside_r, coh_gap, n_mods = _module_coherence_from_corr(r_te)

    orth_within_r = np.nan
    orth_coh_gap = np.nan
    if (
        train_usage is not None
        and train_host_log_expr is not None
        and heldout_host_log_expr is not None
        and eval_level == "cell"
    ):
        _, u_te_orth = orthogonalize_sc_host_expression_v2(
            train_usage.reindex(columns=event_ids).to_numpy(dtype=float),
            train_host_log_expr.reindex(columns=event_ids).to_numpy(dtype=float),
            u_te,
            heldout_host_log_expr.reindex(columns=event_ids).to_numpy(dtype=float),
        )
        r_te_orth, _, _ = precision_pairwise_correlation_v2(u_te_orth, None, min_pairs=min_pairs, min_eff_n=float(min_pairs))
        orth_within_r, _, orth_coh_gap, _ = _module_coherence_from_corr(r_te_orth)

    return {
        "selected_edges": int(len(net.edges)),
        "evaluable_edges": evaluable_edges,
        "heldout_edge_observed_r": heldout_edge_r,
        "discovered_modules": n_mods,
        "heldout_module_within_r": within_r,
        "heldout_module_outside_r": outside_r,
        "heldout_module_coherence_gap": coh_gap,
        "host_orth_within_r": orth_within_r,
        "host_orth_coherence_gap": orth_coh_gap,
    }


def evaluate_module_reproducibility_v2(mod_a: pd.DataFrame, mod_b: pd.DataFrame) -> dict:
    """Measure cross-donor / split module reproducibility on shared assigned events."""
    a = mod_a.set_index("event_id")["module"].astype(str)
    b = mod_b.set_index("event_id")["module"].astype(str)
    shared = a.index.intersection(b.index)
    la = a.loc[shared].to_numpy()
    lb = b.loc[shared].to_numpy()
    both_assigned = (la != "M0") & (lb != "M0")
    if int(both_assigned.sum()) < 5:
        return {"shared_events": int(len(shared)), "both_assigned": int(both_assigned.sum()), "ari": np.nan, "nmi": np.nan}
    return {
        "shared_events": int(len(shared)),
        "both_assigned": int(both_assigned.sum()),
        "ari": float(adjusted_rand_score(la[both_assigned], lb[both_assigned])),
        "nmi": float(normalized_mutual_info_score(la[both_assigned], lb[both_assigned])),
    }


def evaluate_planted_truth_edges_v2(
    shrunk_corr: np.ndarray,
    truth_adj: np.ndarray,
    genes: np.ndarray,
    *,
    budget: int | None = None,
) -> dict:
    """Compute AUPRC and top-budget precision against planted cross-gene co-splicing truth."""
    ia, ib = np.triu_indices(len(genes), k=1)
    elig = (genes[ia] != genes[ib]) & np.isfinite(shrunk_corr[ia, ib])
    if not np.any(elig):
        return {"auprc": np.nan, "precision_at_budget": np.nan}
    scores = np.nan_to_num(shrunk_corr[ia[elig], ib[elig]], nan=-1.0)
    labels = (truth_adj[ia[elig], ib[elig]] > 0).astype(int)
    if labels.sum() == 0 or labels.sum() == len(labels):
        return {"auprc": np.nan, "precision_at_budget": np.nan}
    auprc = float(average_precision_score(labels, scores))
    k = int(labels.sum()) if budget is None else min(int(budget), len(scores))
    top_ix = np.argsort(scores, kind="mergesort")[-k:]
    prec = float(np.mean(labels[top_ix])) if k > 0 else np.nan
    return {"auprc": auprc, "precision_at_budget": prec}


def enrich_sc_module_genes_v2(
    modules_df: pd.DataFrame,
    gmt_path: Path | str,
    *,
    gene_symbol_map: pd.Series | None = None,
    min_set_size: int = 10,
    max_set_size: int = 500,
    min_overlap: int = 3,
) -> pd.DataFrame:
    """Self-contained hypergeometric pathway enrichment per single-cell co-splicing module (unique network genes background)."""
    from scipy.stats import hypergeom

    df = modules_df.copy()
    if gene_symbol_map is not None:
        g_ser = df["event_id"].astype(str).map(gene_symbol_map)
    else:
        g_ser = df["gene_id"].astype(str)
    df["_g"] = g_ser.fillna("").astype(str)
    df = df[df["_g"].map(lambda s: s not in {"", "nan", "NA", "UNKNOWN"} and not any(c in s for c in ";,|"))].copy()

    background = set(df["_g"].unique())
    if len(background) < 15:
        return pd.DataFrame()

    gene_sets: dict[str, set[str]] = {}
    for line in Path(gmt_path).read_text(encoding="utf-8").splitlines():
        parts = line.strip().split("\t")
        if len(parts) >= 3:
            tg = set(parts[2:]) & background
            if min_set_size <= len(tg) <= max_set_size:
                gene_sets[parts[0]] = tg

    rows: list[dict] = []
    active_mods = sorted(m for m in df["module"].astype(str).unique() if m not in {"0", "M0", ""})
    for mod in active_mods:
        mod_genes = set(df.loc[df["module"].astype(str).eq(mod), "_g"].unique()) & background
        if len(mod_genes) < 3:
            continue
        p_vals, mod_rows = [], []
        for term, t_genes in gene_sets.items():
            hits = sorted(mod_genes & t_genes)
            pv = float(hypergeom.sf(len(hits) - 1, len(background), len(t_genes), len(mod_genes)))
            p_vals.append(pv)
            mod_rows.append({
                "module": mod if mod.startswith("M") else f"M{mod}",
                "pathway": term,
                "module_unique_genes": len(mod_genes),
                "network_background_genes": len(background),
                "pathway_genes_in_bg": len(t_genes),
                "overlap": len(hits),
                "p_value": pv,
                "hit_genes": ";".join(hits),
            })
        if mod_rows:
            p_arr = np.asarray(p_vals, dtype=float)
            ord_p = np.argsort(p_arr, kind="mergesort")
            rk = p_arr[ord_p] * len(p_arr) / np.arange(1, len(p_arr) + 1, dtype=float)
            q_s = np.minimum.accumulate(rk[::-1])[::-1]
            q_arr = np.empty_like(p_arr)
            q_arr[ord_p] = np.clip(q_s, 0.0, 1.0)
            for rd, qv in zip(mod_rows, q_arr):
                rd["q_BH"] = float(qv)
            rows.extend(rd for rd in mod_rows if rd["overlap"] >= min_overlap)

    out = pd.DataFrame(rows)
    return out.sort_values(["module", "q_BH", "p_value"]).reset_index(drop=True) if not out.empty else out


def gene_equal_usage_mse_v2(pred: np.ndarray, truth: np.ndarray, genes: np.ndarray) -> tuple[float, pd.Series]:
    """Compute event-mean then gene-equal MSE on observed query labels (matching N117/N118 exact loss)."""
    n = np.isfinite(truth).sum(axis=0)
    ss = np.nansum((pred - truth) ** 2, axis=0)
    ev_mse = np.divide(ss, n, out=np.full(len(n), np.nan), where=n > 0)
    gm = pd.DataFrame({"gene": np.asarray(genes, dtype=str), "mse": ev_mse}).groupby("gene")["mse"].mean()
    return float(gm.mean()), gm


def predict_sc_usage_inductive_v2(
    train_usage: np.ndarray,
    query_usage: np.ndarray,
    genes: np.ndarray,
    *,
    train_support: np.ndarray | None = None,
    train_states: np.ndarray | None = None,
    query_states: np.ndarray | None = None,
    shrunk_corr: np.ndarray | None = None,
    modules: np.ndarray | None = None,
    train_aux_repr: np.ndarray | None = None,
    query_aux_repr: np.ndarray | None = None,
    state_shrink_kappa: float = 0.0,
    ridge: float = 1.0,
    use_module_agg: bool = True,
    top_k_neighbors: int = 25,
) -> np.ndarray:
    """Strictly inductive, target-gene-masked single-cell usage predictor.

    - Baseline anchor uses training-only precision-weighted global/state usage mean (zero query support needed).
    - Network residual predictor uses leave-target-gene-out co-splicing module profiles or top-k cross-gene
      soft-weighted neighbors on `usage` (or explicitly matched auxiliary representation `train_aux_repr` /
      `query_aux_repr` when testing matched-representation adapters), never reading ANY column of `target_gene`
      in `query_usage` or `query_aux_repr`!
    """
    tu = np.asarray(train_usage, dtype=float)
    qu = np.asarray(query_usage, dtype=float)
    g = np.asarray(genes, dtype=str)
    n_tr, n_ev = tu.shape
    n_qu = qu.shape[0]

    tr_w = estimate_sc_observation_weights_v2(tu, train_support, None, prior_bio_var=0.02)
    w_sum = tr_w.sum(axis=0)
    global_mean = np.divide(
        (tr_w * np.where(np.isfinite(tu), tu, 0.0)).sum(axis=0),
        w_sum,
        out=np.full(n_ev, 0.5, dtype=float),
        where=w_sum > 0.0,
    )

    if train_states is not None and query_states is not None and state_shrink_kappa > 0.0:
        ts = np.asarray(train_states, dtype=str)
        qs = np.asarray(query_states, dtype=str)
        tr_anchor = np.broadcast_to(global_mean, (n_tr, n_ev)).copy()
        qu_anchor = np.broadcast_to(global_mean, (n_qu, n_ev)).copy()
        for st in np.unique(ts):
            sel_tr = ts == st
            sw = tr_w[sel_tr].sum(axis=0)
            s_mu = np.divide(
                (tr_w[sel_tr] * np.where(np.isfinite(tu[sel_tr]), tu[sel_tr], 0.0)).sum(axis=0),
                sw,
                out=global_mean,
                where=sw > 0.0,
            )
            lam = sw / (sw + float(state_shrink_kappa))
            shrunk_st_mu = (1.0 - lam) * global_mean + lam * s_mu
            tr_anchor[sel_tr] = shrunk_st_mu
            sel_qu = qs == st
            if np.any(sel_qu):
                qu_anchor[sel_qu] = shrunk_st_mu
    else:
        tr_anchor = np.broadcast_to(global_mean, (n_tr, n_ev)).copy()
        qu_anchor = np.broadcast_to(global_mean, (n_qu, n_ev)).copy()

    if shrunk_corr is None or ridge <= 0.0:
        return np.clip(qu_anchor, 0.0, 1.0)

    # Feature matrix for neighbors (strictly centered & scaled on training cells ONLY)
    x_tr_raw = tu if train_aux_repr is None else np.asarray(train_aux_repr, dtype=float)
    x_qu_raw = qu if query_aux_repr is None else np.asarray(query_aux_repr, dtype=float)
    f_tr = np.isfinite(x_tr_raw)
    xm = np.divide(np.nansum(x_tr_raw, axis=0), np.maximum(f_tr.sum(axis=0), 1), out=np.zeros(n_ev), where=f_tr.sum(axis=0) > 0)
    xz_tr = np.where(f_tr, x_tr_raw - xm[None, :], 0.0)
    xsd = np.sqrt(np.divide((xz_tr * xz_tr).sum(axis=0), np.maximum(f_tr.sum(axis=0), 1)))
    xsd = np.where(xsd > 1e-8, xsd, 1.0)
    xz_tr = xz_tr / xsd[None, :]
    xz_qu = np.where(np.isfinite(x_qu_raw), (x_qu_raw - xm[None, :]) / xsd[None, :], 0.0)

    y_res = tu - tr_anchor
    pred_qu = qu_anchor.copy()
    same_gene = g[:, None] == g[None, :]
    w_mat = np.where(np.isfinite(shrunk_corr) & (~same_gene), shrunk_corr, 0.0)
    mod_arr = np.asarray(modules, dtype=str) if modules is not None else np.full(n_ev, "M0", dtype=str)
    active_mods = [m for m in sorted(set(mod_arr)) if m != "M0"]

    for j in range(n_ev):
        valid_y = np.isfinite(y_res[:, j]) & (tr_w[:, j] > 0.0)
        if int(valid_y.sum()) < 6:
            continue
        allowed = ~same_gene[j]
        # Construct compact low-noise feature vector Z_j (module centroids + top cross-gene neighbors, all excluding gene[j]!)
        feat_tr_cols: list[np.ndarray] = []
        feat_qu_cols: list[np.ndarray] = []

        if use_module_agg and active_mods:
            for m_name in active_mods:
                members = np.flatnonzero((mod_arr == m_name) & allowed)
                if len(members) < 2:
                    continue
                # Confidence-preserving module summary
                mw = np.maximum(np.abs(w_mat[j, members]), 0.05)
                conf = float(np.sqrt(np.mean(np.abs(w_mat[j, members])))) if np.any(w_mat[j, members] != 0) else 0.25
                mw = (mw / mw.sum()) * max(conf, 0.15)
                feat_tr_cols.append(xz_tr[:, members] @ mw)
                feat_qu_cols.append(xz_qu[:, members] @ mw)

        # Direct soft-weighted top cross-gene neighbor summary
        nb_scores = np.abs(w_mat[j]) * allowed
        pos_nb = np.flatnonzero(nb_scores > 1e-6)
        if len(pos_nb) > 0:
            if len(pos_nb) > top_k_neighbors:
                pos_nb = pos_nb[np.argsort(nb_scores[pos_nb], kind="mergesort")[-top_k_neighbors:]]
            nw = w_mat[j, pos_nb]
            scale_c = float(np.sqrt(np.mean(np.abs(nw))))
            nw_norm = (nw / np.sum(np.abs(nw))) * scale_c
            feat_tr_cols.append(xz_tr[:, pos_nb] @ nw_norm)
            feat_qu_cols.append(xz_qu[:, pos_nb] @ nw_norm)

        if not feat_tr_cols:
            continue

        Z_tr = np.column_stack(feat_tr_cols)
        Z_qu = np.column_stack(feat_qu_cols)
        # Scale predictor columns on training cells so ridge is invariant to global correlation shrinkage factor
        z_col_sd = np.maximum(np.std(Z_tr[valid_y], axis=0), 1e-4)
        Z_tr = Z_tr / z_col_sd[None, :]
        Z_qu = Z_qu / z_col_sd[None, :]
        wv = tr_w[valid_y, j]
        wv = wv / np.mean(wv)
        Zv = Z_tr[valid_y]
        yv = y_res[valid_y, j]
        A = (Zv.T * wv) @ Zv + float(ridge) * Z_tr.shape[1] * np.eye(Z_tr.shape[1])
        b = (Zv.T * wv) @ yv
        beta = np.linalg.solve(A, b)
        pred_qu[:, j] = qu_anchor[:, j] + Z_qu @ beta

    return np.clip(pred_qu, 0.0, 1.0)
