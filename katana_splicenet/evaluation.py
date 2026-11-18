"""Held-out donor PC1 projection, cross-gene coherence, and host-depth orthogonalization."""
from __future__ import annotations

import math
import numpy as np
import pandas as pd


def score_modules_v2(
    measured: pd.DataFrame, assignments: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Compute training-standardized PC1 module scores and frozen basis loadings."""
    scores = pd.DataFrame(index=measured.columns)
    coverage = pd.DataFrame(index=measured.columns)
    rows: list[dict] = []
    basis: list[dict] = []

    labels = assignments.set_index("event_key").reindex(measured.index)["module"]
    for mod in sorted(labels.dropna().unique()):
        if int(mod) == 0:
            continue
        block = measured.loc[labels.eq(mod)]
        name = f"M{int(mod)}"
        means = block.mean(axis=1)
        scales = block.std(axis=1, ddof=1).clip(lower=1e-8)
        z = ((block.T.fillna(means) - means) / scales).to_numpy(dtype=float)
        _, s_vals, vt = np.linalg.svd(z, full_matrices=False)
        load = vt[0]
        raw = z @ load
        sign = 1.0 if float(load.sum()) >= 0.0 else -1.0
        scale = max(float(np.std(raw, ddof=1)), 1e-12)
        load = load * sign / scale
        coverage[name] = block.notna().mean(axis=0)
        scores[name] = pd.Series(z @ load, index=measured.columns).where(coverage[name].ge(0.80))
        rows.append({
            "module": name,
            "events": len(block),
            "pc1_variance_fraction": float(s_vals[0] ** 2 / max(float(s_vals @ s_vals), 1e-12)),
            "evaluable_samples": int(scores[name].notna().sum()),
        })
        for j, key in enumerate(block.index):
            basis.append({
                "event_key": str(key),
                "module": name,
                "center": float(means.iloc[j]),
                "scale": float(scales.iloc[j]),
                "loading": float(load[j]),
            })
    return scores, coverage, pd.DataFrame(rows), pd.DataFrame(basis)


def visit_correlation(observed: np.ndarray, expected: np.ndarray, min_visits: int = 6) -> tuple[float, int]:
    x = np.asarray(observed, dtype=float)
    y = np.asarray(expected, dtype=float)
    usable = np.isfinite(x) & np.isfinite(y)
    count = int(usable.sum())
    if count < min_visits:
        return np.nan, count
    x = x[usable] - x[usable].mean()
    y = y[usable] - y[usable].mean()
    denom = float(np.linalg.norm(x) * np.linalg.norm(y))
    return (float(x @ y / denom) if denom > 1e-12 else np.nan), count


def pair_correlations(values: np.ndarray, idx_a: np.ndarray, idx_b: np.ndarray, min_obs: int = 6) -> np.ndarray:
    """Pairwise Pearson r across columns for specified row pairs, requiring >= min_obs shared finite samples."""
    va = values[idx_a]
    vb = values[idx_b]
    mask = np.isfinite(va) & np.isfinite(vb)
    counts = mask.sum(axis=1)
    out = np.full(len(idx_a), np.nan, dtype=float)
    ok = counts >= min_obs
    if not ok.any():
        return out
    va_ok = np.where(mask[ok], va[ok], 0.0)
    vb_ok = np.where(mask[ok], vb[ok], 0.0)
    m_ok = mask[ok]
    n_ok = counts[ok].astype(float)
    ma = va_ok.sum(axis=1, keepdims=True) / n_ok[:, None]
    mb = vb_ok.sum(axis=1, keepdims=True) / n_ok[:, None]
    da = np.where(m_ok, va_ok - ma, 0.0)
    db = np.where(m_ok, vb_ok - mb, 0.0)
    denom = np.sqrt((da * da).sum(axis=1) * (db * db).sum(axis=1))
    valid_den = denom > 1e-12
    sub_out = np.full(len(n_ok), np.nan, dtype=float)
    sub_out[valid_den] = (da[valid_den] * db[valid_den]).sum(axis=1) / denom[valid_den]
    out[ok] = np.clip(sub_out, -1.0, 1.0)
    return out


def orthogonalize_host_depth(
    train_usage: pd.DataFrame,
    test_usage: pd.DataFrame,
    train_support: pd.DataFrame,
    test_support: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Regress out log2 cluster read depth on training samples ONLY, applying frozen slopes to test samples."""
    tr_u = train_usage.to_numpy(dtype=float)
    te_u = test_usage.to_numpy(dtype=float)
    tr_x = np.log2(np.maximum(train_support.to_numpy(dtype=float), 1.0))
    te_x = np.log2(np.maximum(test_support.to_numpy(dtype=float), 1.0))

    tr_res = np.full_like(tr_u, np.nan)
    te_res = np.full_like(te_u, np.nan)
    for i in range(len(tr_u)):
        ok = np.isfinite(tr_u[i]) & np.isfinite(tr_x[i])
        if ok.sum() < 10:
            continue
        xm = float(tr_x[i, ok].mean())
        ym = float(tr_u[i, ok].mean())
        dx = tr_x[i, ok] - xm
        dy = tr_u[i, ok] - ym
        v_x = float(dx @ dx)
        beta = float(dx @ dy / v_x) if v_x > 1e-12 else 0.0
        # Residual + training mean preserves baseline scale
        tr_ok = np.isfinite(tr_u[i]) & np.isfinite(tr_x[i])
        te_ok = np.isfinite(te_u[i]) & np.isfinite(te_x[i])
        tr_res[i, tr_ok] = tr_u[i, tr_ok] - beta * (tr_x[i, tr_ok] - xm)
        te_res[i, te_ok] = te_u[i, te_ok] - beta * (te_x[i, te_ok] - xm)

    return (
        pd.DataFrame(tr_res, index=train_usage.index, columns=train_usage.columns),
        pd.DataFrame(te_res, index=test_usage.index, columns=test_usage.columns),
    )


def select_best_training_module(
    train_usage: pd.DataFrame,
    train_meta: pd.DataFrame,
    assignments: pd.DataFrame,
    mapping: pd.DataFrame,
) -> tuple[str, pd.DataFrame, pd.Series] | None:
    scores, coverage, modules, basis = score_modules_v2(train_usage, assignments)
    if basis.empty:
        return None
    observed = train_meta.set_index("sample_id").loc[train_usage.columns, ["biological_sample_id", "time_day"]]
    donors = observed["biological_sample_id"].unique()
    candidates: list[dict] = []
    for mod, loadings in basis.groupby("module", sort=False):
        vals = observed.assign(score=scores[mod])
        tr_r = []
        for donor in donors:
            held = vals[vals["biological_sample_id"].eq(donor)]
            prof = vals[vals["biological_sample_id"].ne(donor)].groupby("time_day")["score"].mean()
            r, _ = visit_correlation(held["score"].to_numpy(), held["time_day"].map(prof).to_numpy())
            if np.isfinite(r):
                tr_r.append(r)
        n_genes = int(mapping.loc[loadings["event_key"], "gene_id"].nunique())
        if len(loadings) >= 20 and n_genes >= 2 and len(tr_r) >= math.ceil(0.8 * len(donors)):
            candidates.append({"module": mod, "median_r": float(np.median(tr_r))})
    if not candidates:
        return None
    best = min(candidates, key=lambda r: (-r["median_r"], int(r["module"][1:])))
    label = best["module"]
    load_df = basis[basis["module"].eq(label)].set_index("event_key")
    profile = scores[label].groupby(observed["time_day"]).mean()
    return label, load_df, profile


def heldout_trajectory_and_coherence(
    test_usage: pd.DataFrame,
    test_meta: pd.DataFrame,
    mapping: pd.DataFrame,
    selection: tuple[str, pd.DataFrame, pd.Series] | None,
    controls: list[str],
) -> list[dict]:
    if selection is None:
        return [
            {
                "donor": str(donor),
                "module": "",
                "trajectory_r": np.nan,
                "coherence_gap": np.nan,
            }
            for donor, _ in test_meta.groupby("biological_sample_id", sort=False)
        ]

    label, basis, profile = selection
    keys = basis.index.tolist()
    block = test_usage.loc[keys]
    z = (block.T.fillna(basis["center"]) - basis["center"]) / basis["scale"]
    cov = block.notna().mean(axis=0)
    projected = pd.Series(z.to_numpy(dtype=float) @ basis["loading"].to_numpy(dtype=float), index=test_usage.columns).where(
        cov.ge(0.80)
    )

    genes = mapping.loc[keys, "gene_id"]
    mult = genes.map(genes.value_counts()).to_numpy(dtype=float)
    ia, ib = np.triu_indices(len(keys), k=1)
    cross = genes.to_numpy()[ia] != genes.to_numpy()[ib]
    ia, ib = ia[cross], ib[cross]
    w_in = 1.0 / (mult[ia] * mult[ib])

    ca = np.repeat(np.arange(len(keys)), len(controls))
    cb = np.tile(np.arange(len(controls)) + len(keys), len(keys))
    w_out = 1.0 / mult[ca]

    joined = test_usage.loc[keys + list(controls)]
    rows: list[dict] = []
    for donor, obs in test_meta.groupby("biological_sample_id", sort=False):
        s_names = obs["sample_id"].tolist()
        pred = projected.loc[s_names].to_numpy(dtype=float)
        exp = obs["time_day"].map(profile).to_numpy(dtype=float)
        r_val, n_vis = visit_correlation(pred, exp)
        mat = joined[s_names].to_numpy(dtype=float)

        def _w_coh(idx1: np.ndarray, idx2: np.ndarray, w: np.ndarray) -> float:
            if not len(idx1):
                return np.nan
            ws = 0.0
            wt = 0.0
            for st in range(0, len(idx1), 40000):
                sl = slice(st, st + 40000)
                rr = pair_correlations(mat, idx1[sl], idx2[sl])
                kp = np.isfinite(rr)
                ww = w[sl][kp]
                ws += float(rr[kp] @ ww)
                wt += float(ww.sum())
            return ws / wt if wt > 0 else np.nan

        c_in = _w_coh(ia, ib, w_in)
        c_out = _w_coh(ca, cb, w_out)
        rows.append({
            "donor": str(donor),
            "module": label,
            "paired_visits": n_vis,
            "trajectory_r": r_val,
            "within_signed_r": c_in,
            "outside_signed_r": c_out,
            "coherence_gap": c_in - c_out if (np.isfinite(c_in) and np.isfinite(c_out)) else np.nan,
        })
    return rows


def enrich_module_genes_v2(
    assignments: pd.DataFrame,
    event_to_gene: pd.Series,
    gmt_path: Path,
    *,
    min_set_size: int = 10,
    max_set_size: int = 500,
    min_overlap: int = 3,
) -> pd.DataFrame:
    """Hypergeometric pathway enrichment per module using unique deduplicated network genes as background."""
    from pathlib import Path as _P
    from scipy.stats import hypergeom
    def _bh_fdr(p_values: np.ndarray) -> np.ndarray:
        p = np.asarray(p_values, dtype=float)
        n = len(p)
        if n == 0:
            return p
        order = np.argsort(p, kind="mergesort")
        ranked = p[order] * n / np.arange(1, n + 1, dtype=float)
        q_sorted = np.minimum.accumulate(ranked[::-1])[::-1]
        out = np.empty(n, dtype=float)
        out[order] = np.clip(q_sorted, 0.0, 1.0)
        return out

    asg = assignments.set_index("event_key")["module"]
    genes = event_to_gene.reindex(asg.index).dropna().astype(str)
    valid_gene = genes.map(lambda s: s not in {"", "nan", "NA", "UNKNOWN"} and not any(c in s for c in ";,|"))
    genes = genes[valid_gene]
    asg = asg.reindex(genes.index)

    # Background is strictly unique genes across ALL network-entered events (including grey module 0), counted ONCE per gene
    background = set(genes.unique())
    if len(background) < 20:
        return pd.DataFrame()

    gene_sets: dict[str, set[str]] = {}
    for line in _P(gmt_path).read_text(encoding="utf-8").splitlines():
        parts = line.strip().split("\t")
        if len(parts) < 3:
            continue
        term_genes = set(parts[2:]) & background
        if min_set_size <= len(term_genes) <= max_set_size:
            gene_sets[parts[0]] = term_genes

    rows: list[dict] = []
    for mod in sorted(set(asg.dropna().astype(int)) - {0}):
        mod_genes = set(genes[asg.eq(mod)].unique()) & background
        if len(mod_genes) < 5:
            continue
        p_vals = []
        mod_rows = []
        for term, t_genes in gene_sets.items():
            hits = sorted(mod_genes & t_genes)
            p_val = float(hypergeom.sf(len(hits) - 1, len(background), len(t_genes), len(mod_genes)))
            p_vals.append(p_val)
            mod_rows.append({
                "module": f"M{mod}",
                "pathway": term,
                "module_unique_genes": len(mod_genes),
                "network_background_genes": len(background),
                "pathway_genes_in_bg": len(t_genes),
                "overlap": len(hits),
                "p_value": p_val,
                "hit_genes": ";".join(hits),
            })
        if mod_rows:
            q_vals = _bh_fdr(np.asarray(p_vals, dtype=float))
            for r_dict, qv in zip(mod_rows, q_vals):
                r_dict["q_BH"] = float(qv)
            rows.extend(r_dict for r_dict in mod_rows if r_dict["overlap"] >= min_overlap)

    out = pd.DataFrame(rows)
    return out.sort_values(["module", "q_BH", "p_value"]).reset_index(drop=True) if not out.empty else out
