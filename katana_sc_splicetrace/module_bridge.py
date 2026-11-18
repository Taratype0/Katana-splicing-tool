"""Cross-Package Regulatory Bridge: Connect SpliceNet Co-Splicing Modules -> SpliceTrace RBP Motifs.

Lives inside `katana-splicetrace` (Bulk) and `katana-sc-splicetrace` (Single-Cell) because PWM sequence scanning
requires `--fasta` and the bundled v3.0 RBP library (`621` PWMs), keeping `splicenet` pure matrix-only.
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd


def bridge_splicenet_modules_to_rbps(
    modules_df: pd.DataFrame,
    events_df: pd.DataFrame,
    fasta_path: Path | str,
    library: object,
    *,
    extract_rna_fn,
    score_odds_fn,
    min_module_events: int = 8,
    n_permutations: int = 999,
    seed: int = 20261118,
) -> pd.DataFrame:
    """Test every non-grey SpliceNet co-splicing module against all curated RBP motifs with GC/length orthogonalization
    and gene-block Rademacher max-T FWER calibration.
    """
    id_col = "event_id" if "event_id" in modules_df.columns else "event_key"
    if id_col not in modules_df.columns or "module" not in modules_df.columns:
        raise ValueError("modules_df must contain ('event_id' or 'event_key') and 'module'")

    mod = modules_df[[c for c in (id_col, "module", "kME") if c in modules_df.columns]].copy()
    mod = mod.rename(columns={id_col: "event_id"})
    mod["event_id"] = mod["event_id"].astype(str)

    def _norm_mod_label(v: object) -> str:
        s = str(v).strip()
        if s.lower() in ("0", "m0", "unassigned", "grey", "gray", "none", "nan", ""):
            return "M0"
        if s.isdigit():
            return f"M{int(s)}"
        return s

    mod["mod_label"] = mod["module"].map(_norm_mod_label)
    if "kME" not in mod.columns:
        mod["kME"] = np.where(mod["mod_label"].ne("M0"), 1.0, 0.0)

    ev = events_df.copy()
    ev["event_id"] = ev["event_id"].astype(str)
    merged = ev.merge(mod, on="event_id", how="inner").drop_duplicates("event_id").reset_index(drop=True)
    non_grey = sorted(
        [m for m in merged["mod_label"].unique() if m != "M0"],
        key=lambda x: (int(x[1:]) if x.startswith("M") and x[1:].isdigit() else 999, x),
    )
    if not non_grey or len(merged) < 20:
        return pd.DataFrame(columns=["module", "module_events", "module_genes", "RBP", "best_region", "stat_Z", "nominal_P", "module_P_FWER"])

    if "comparison_id" not in merged.columns:
        merged["comparison_id"] = "splicenet_bridge"
    rna_df = extract_rna_fn(merged, Path(fasta_path), acceptor_core_nt=3, feature_only=False, compute_sha256=False)

    models = library.models
    m_ids = list(models.keys())
    mid_idx = {mid: i for i, mid in enumerate(m_ids)}
    rbps = sorted(library.rbp_to_models.keys())
    rbp_cols = {r: [mid_idx[m] for m in library.rbp_to_models[r] if m in mid_idx] for r in rbps}

    ev_index = {eid: i for i, eid in enumerate(merged["event_id"])}
    n_ev = len(merged)
    reg_slots = ["mature_delta", "upstream_intron", "downstream_intron"]
    F_mod = np.zeros((n_ev, len(reg_slots), len(m_ids)), dtype=np.float64)

    mat_sub = rna_df[rna_df["layer"].eq("mature_local_form") & rna_df["region"].isin(["inclusion", "exclusion"])]
    if not mat_sub.empty:
        ms = score_odds_fn(mat_sub["sequence"].astype(str).tolist(), models)
        inc_m = np.zeros((n_ev, len(m_ids)), dtype=np.float64)
        exc_m = np.zeros((n_ev, len(m_ids)), dtype=np.float64)
        for ri, r_tup in enumerate(mat_sub[["event_id", "region"]].itertuples(index=False)):
            ei = ev_index.get(str(r_tup.event_id))
            if ei is not None:
                if r_tup.region == "inclusion":
                    inc_m[ei] = ms[ri]
                else:
                    exc_m[ei] = ms[ri]
        F_mod[:, 0, :] = inc_m - exc_m

    reg_sub = rna_df[rna_df["layer"].eq("regulatory_pre_mRNA")]
    flank_gc = np.full(n_ev, 0.5, dtype=np.float64)
    if not reg_sub.empty:
        rs = score_odds_fn(reg_sub["sequence"].astype(str).tolist(), models)
        up_cnt = np.zeros(n_ev, dtype=np.float64)
        dn_cnt = np.zeros(n_ev, dtype=np.float64)
        gc_num = np.zeros(n_ev, dtype=np.float64)
        gc_den = np.zeros(n_ev, dtype=np.float64)
        for ri, r_tup in enumerate(reg_sub[["event_id", "region", "sequence"]].itertuples(index=False)):
            ei = ev_index.get(str(r_tup.event_id))
            if ei is None:
                continue
            seq_s = str(r_tup.sequence).upper()
            if seq_s:
                gc_num[ei] += seq_s.count("G") + seq_s.count("C")
                gc_den[ei] += len(seq_s)
            if "upstream" in str(r_tup.region) or "acceptor" in str(r_tup.region) or "3ss" in str(r_tup.region):
                F_mod[ei, 1, :] += rs[ri]
                up_cnt[ei] += 1.0
            else:
                F_mod[ei, 2, :] += rs[ri]
                dn_cnt[ei] += 1.0
        F_mod[:, 1, :] /= np.maximum(up_cnt[:, None], 1.0)
        F_mod[:, 2, :] /= np.maximum(dn_cnt[:, None], 1.0)
        flank_gc = np.divide(gc_num, gc_den, out=np.full(n_ev, 0.5), where=gc_den > 0)

    F_rbp = np.zeros((n_ev, len(reg_slots), len(rbps)), dtype=np.float64)
    for r_i, r_name in enumerate(rbps):
        cols_r = rbp_cols[r_name]
        if cols_r:
            F_rbp[:, :, r_i] = np.max(F_mod[:, :, cols_r], axis=2)

    F_flat = F_rbp.reshape(n_ev, len(reg_slots) * len(rbps))
    X_nuis = np.column_stack([np.ones(n_ev), flank_gc])
    beta_f, _, _, _ = np.linalg.lstsq(X_nuis, F_flat, rcond=None)
    F_res = F_flat - X_nuis @ beta_f
    f_sd = np.maximum(F_res.std(axis=0, ddof=1), 1e-12)
    F_std = np.where(F_res.var(axis=0) > 1e-11, F_res / f_sd[None, :], 0.0)

    genes = merged["gene_id"].astype(str).to_numpy()
    uniq_g, g_codes = np.unique(genes, return_inverse=True)
    n_genes = len(uniq_g)
    rng = np.random.default_rng(int(seed))
    signs = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(int(n_permutations), n_genes))

    out_rows: list[dict] = []
    for m_lbl in non_grey:
        in_m = merged["mod_label"].to_numpy() == m_lbl
        n_m_ev = int(in_m.sum())
        if n_m_ev < int(min_module_events):
            continue
        n_m_genes = int(len(set(genes[in_m])))
        kme_vals = np.nan_to_num(merged["kME"].to_numpy(dtype=float), nan=0.0)
        if float(np.nanstd(kme_vals[in_m])) > 0.05:
            y_raw = np.where(in_m, np.maximum(np.abs(kme_vals), 0.25), 0.0)
        else:
            y_raw = in_m.astype(float)
        beta_y, _, _, _ = np.linalg.lstsq(X_nuis, y_raw, rcond=None)
        y_res = y_raw - X_nuis @ beta_y

        contrib = y_res[:, None] * F_std
        S_g = np.zeros((n_genes, F_std.shape[1]), dtype=np.float64)
        np.add.at(S_g, g_codes, contrib)
        den_g = np.sqrt(np.sum(S_g * S_g, axis=0))
        valid_c = den_g > 1e-12
        S_norm = np.divide(S_g, den_g[None, :], out=np.zeros_like(S_g), where=valid_c[None, :])
        Z_obs = S_norm.sum(axis=0)
        null_Z = (signs @ np.ascontiguousarray(S_norm, dtype=np.float32)).astype(np.float64)

        Z_2d = Z_obs.reshape(len(reg_slots), len(rbps))
        null_3d = null_Z.reshape(int(n_permutations), len(reg_slots), len(rbps))
        rbp_max = np.max(np.abs(Z_2d), axis=0)
        null_rbp_max = np.max(np.abs(null_3d), axis=1)
        glob_max = np.max(null_rbp_max, axis=1)

        nom_p = (1.0 + np.sum(null_rbp_max >= rbp_max[None, :] - 1e-6, axis=0)) / (n_permutations + 1.0)
        fwer_p = (1.0 + np.sum(glob_max[:, None] >= rbp_max[None, :] - 1e-6, axis=0)) / (n_permutations + 1.0)

        order = np.argsort(-rbp_max)[:15]
        for r_i in order:
            best_reg = int(np.argmax(np.abs(Z_2d[:, r_i])))
            out_rows.append({
                "module": m_lbl,
                "module_events": n_m_ev,
                "module_genes": n_m_genes,
                "RBP": rbps[r_i],
                "best_region": reg_slots[best_reg],
                "stat_Z": round(float(Z_2d[best_reg, r_i]), 4),
                "nominal_P": round(float(nom_p[r_i]), 5),
                "module_P_FWER": round(float(fwer_p[r_i]), 5),
            })

    return pd.DataFrame(out_rows)
