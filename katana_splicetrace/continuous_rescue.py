"""Self-contained MotifTrace v2.2 small-sample 2v2 rescue engine (11-channel PPT + ESE/ESS + CpG FWL)."""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm

from katana_splicetrace.library import CuratedMotifLibrary, load_bundled_human_rbp_library
from katana_splicetrace.rna_forms import (
    compute_flank_gc_from_rna,
    compute_sequence_gc,
    event_forms,
    extract_event_rna_v2,
    trace_motif_instances_v2,
)
from katana_splicetrace.fast_accel import configure_server_threads, fast_gene_block_maxt_permutations, fast_score_mean_log_odds


@dataclass(frozen=True)
class MotifRescueResultV22:
    rbp_ranking: pd.DataFrame
    region_tests: pd.DataFrame
    event_contrasts: pd.DataFrame
    mature_consequences: pd.DataFrame
    diagnostics: dict


def _score_mean_log_odds(sequences: list[str], models: dict[str, dict]) -> np.ndarray:
    return fast_score_mean_log_odds(sequences, models)


def build_directional_11channel_tensor_v22(
    events_df: pd.DataFrame,
    rna_df: pd.DataFrame,
    library: CuratedMotifLibrary,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    models = library.models
    m_ids = list(models.keys())
    mid_to_idx = {mid: idx for idx, mid in enumerate(m_ids)}
    rbps = sorted(library.rbp_to_models.keys())
    rbp_to_mcols = {r: [mid_to_idx[m] for m in library.rbp_to_models[r] if m in mid_to_idx] for r in rbps}

    event_ids = events_df["event_id"].astype(str).to_numpy()
    ev_index = {eid: i for i, eid in enumerate(event_ids)}
    n_ev = len(event_ids)
    region_names = [
        "mature_consequence_delta",
        "upstream_intron_0_35",
        "upstream_intron_35_85",
        "upstream_intron_85_200",
        "downstream_intron_0_35",
        "downstream_intron_35_85",
        "downstream_intron_85_200",
        "upstream_exon_3p_50",
        "alt_exon_5p_50",
        "alt_exon_3p_50",
        "downstream_exon_5p_50",
    ]
    model_tensor = np.zeros((n_ev, len(region_names), len(m_ids)), dtype=float)

    mature = rna_df[rna_df["layer"].eq("mature_local_form") & rna_df["region"].isin(["inclusion", "exclusion"])].copy()
    if not mature.empty:
        m_scores = _score_mean_log_odds(mature["sequence"].astype(str).tolist(), models)
        inc_mat = np.zeros((n_ev, len(m_ids)), dtype=float)
        exc_mat = np.zeros((n_ev, len(m_ids)), dtype=float)
        ex_seqs, ex_targets = [], []
        for row_i, rec in enumerate(mature[["event_id", "region", "sequence", "junctions_json"]].itertuples(index=False)):
            ei = ev_index.get(str(rec.event_id))
            if ei is None:
                continue
            if rec.region == "inclusion":
                inc_mat[ei] = m_scores[row_i]
                seq_inc = str(rec.sequence)
                joins = json.loads(str(rec.junctions_json)) if rec.junctions_json else []
                if len(joins) >= 2:
                    j1, j2 = int(joins[0]), int(joins[-1])
                    for piece, sl_i in (
                        (seq_inc[max(0, j1 - 50) : j1], 7),
                        (seq_inc[j1 : min(j2, j1 + 50)], 8),
                        (seq_inc[max(j1, j2 - 50) : j2], 9),
                        (seq_inc[j2 : min(len(seq_inc), j2 + 50)], 10),
                    ):
                        if len(piece) >= 6:
                            ex_seqs.append(piece)
                            ex_targets.append((ei, sl_i))
                elif len(joins) == 1:
                    j1 = int(joins[0])
                    for piece, sl_i in ((seq_inc[max(0, j1 - 50) : j1], 7), (seq_inc[j1 : min(len(seq_inc), j1 + 50)], 10)):
                        if len(piece) >= 6:
                            ex_seqs.append(piece)
                            ex_targets.append((ei, sl_i))
            elif rec.region == "exclusion":
                exc_mat[ei] = m_scores[row_i]
        model_tensor[:, 0, :] = inc_mat - exc_mat
        if ex_seqs:
            ex_scores = _score_mean_log_odds(ex_seqs, models)
            for idx_e, (ei, sl_i) in enumerate(ex_targets):
                model_tensor[ei, sl_i, :] = ex_scores[idx_e]

    reg_df = rna_df[rna_df["layer"].eq("regulatory_pre_mRNA")].copy()
    if not reg_df.empty:
        etype_map = {str(eid): str(et) for eid, et in zip(events_df["event_id"], events_df["event_type"])}
        sub_seqs, sub_targets = [], []
        for rec in reg_df[["event_id", "region", "sequence"]].itertuples(index=False):
            ei = ev_index.get(str(rec.event_id))
            if ei is None:
                continue
            base_r = str(rec.region)
            seq = str(rec.sequence)
            if not seq or base_r.endswith("_distal"):
                continue
            kind = etype_map.get(str(rec.event_id), "SE")
            route = None
            if kind == "SE":
                if base_r == "inclusion_junction_1_acceptor_intron":
                    route = (True, 1.0)
                elif base_r == "inclusion_junction_2_donor_intron":
                    route = (False, 1.0)
            elif kind == "RI":
                if base_r == "exclusion_junction_1_acceptor_intron":
                    route = (True, 1.0)
                elif base_r == "exclusion_junction_1_donor_intron":
                    route = (False, 1.0)
            elif kind == "A3SS":
                if base_r == "inclusion_junction_1_acceptor_intron":
                    route = (True, 1.0)
                elif base_r == "exclusion_junction_1_acceptor_intron":
                    route = (True, -1.0)
                elif base_r == "inclusion_junction_1_donor_intron":
                    route = (False, 1.0)
            elif kind == "A5SS":
                if base_r == "inclusion_junction_1_donor_intron":
                    route = (False, 1.0)
                elif base_r == "exclusion_junction_1_donor_intron":
                    route = (False, -1.0)
                elif base_r == "inclusion_junction_1_acceptor_intron":
                    route = (True, 1.0)
            elif kind == "MXE":
                if base_r == "inclusion_junction_1_acceptor_intron":
                    route = (True, 1.0)
                elif base_r == "exclusion_junction_1_acceptor_intron":
                    route = (True, -1.0)
                elif base_r == "inclusion_junction_2_donor_intron":
                    route = (False, 1.0)
                elif base_r == "exclusion_junction_2_donor_intron":
                    route = (False, -1.0)
            if route is None:
                continue
            is_up, sgn = route
            L = len(seq)
            for b_i, (d0, d1) in enumerate(((0, 35), (35, 85), (85, 200))):
                piece = seq[max(0, L - d1) : max(0, L - d0)] if is_up else seq[min(L, d0) : min(L, d1)]
                if len(piece) >= 6:
                    sub_seqs.append(piece)
                    sub_targets.append((ei, (1 if is_up else 4) + b_i, sgn))
        if sub_seqs:
            b_scores = _score_mean_log_odds(sub_seqs, models)
            for idx_s, (ei, slot, sgn) in enumerate(sub_targets):
                model_tensor[ei, slot, :] += sgn * b_scores[idx_s]

    rms_m = np.sqrt(np.mean(model_tensor * model_tensor, axis=0, keepdims=True))
    norm_m = np.divide(model_tensor, rms_m, out=np.zeros_like(model_tensor), where=rms_m > 1e-12)
    rbp_tensor = np.zeros((n_ev, len(region_names), len(rbps)), dtype=float)
    for r_i, r_name in enumerate(rbps):
        cols = rbp_to_mcols[r_name]
        if cols:
            rbp_tensor[:, :, r_i] = np.mean(norm_m[:, :, cols], axis=2)
    return rbp_tensor, np.asarray(rbps, dtype=str), region_names


def _fwl_residuals(Y: np.ndarray, X: np.ndarray, weights: np.ndarray) -> np.ndarray:
    w_root = np.sqrt(np.maximum(weights, 1e-12))
    Xw = X * w_root[:, None]
    if Y.ndim == 1:
        beta, _, _, _ = np.linalg.lstsq(Xw, Y * w_root, rcond=None)
        return Y - X @ beta
    beta, _, _, _ = np.linalg.lstsq(Xw, Y * w_root[:, None], rcond=None)
    return Y - X @ beta


def run_bulk_continuous_rescue_v22(
    all_measurable_events: pd.DataFrame,
    fasta_path: Path | str,
    library: CuratedMotifLibrary | None = None,
    *,
    max_events_blind_depth: int = 3500,
    n_permutations: int = 1999,
    seed: int = 20261117,
    threads: int | None = None,
) -> MotifRescueResultV22:
    """Run Katana-SpliceTrace v2.4 multi-core calibrated rescue engine."""
    configure_server_threads(threads)
    lib = load_bundled_human_rbp_library() if library is None else library
    df = all_measurable_events.copy()
    if "gene_symbol" not in df.columns:
        df["gene_symbol"] = df["gene_id"].astype(str)
    if "comparison_id" not in df.columns:
        df["comparison_id"] = "rescue_contrast"
    if "depth" not in df.columns:
        df["depth"] = 50.0

    p_ctrl = np.clip(df["baseline_PSI"].to_numpy(dtype=float), 0.0, 1.0)
    p_kd = np.clip(p_ctrl + df["dpsi"].to_numpy(dtype=float), 0.0, 1.0)
    non_const = (0.5 * (p_ctrl + p_kd) >= 0.03) & (0.5 * (p_ctrl + p_kd) <= 0.97)
    if int(non_const.sum()) >= 120:
        df = df.loc[non_const].reset_index(drop=True)

    fdr_all = pd.to_numeric(df.get("fdr", pd.Series(1.0, index=df.index)), errors="coerce").fillna(1.0).to_numpy(dtype=float)
    pval_all = np.clip(pd.to_numeric(df.get("PValue", df.get("fdr", 1.0)), errors="coerce").fillna(1.0).to_numpy(dtype=float), 1e-12, 1.0)
    dpsi_v = df["dpsi"].to_numpy(dtype=float)
    df["log_depth"] = np.log1p(np.maximum(df["depth"].to_numpy(dtype=float), 1.0))

    is_fg = (fdr_all <= 0.10) & (np.abs(dpsi_v) >= 0.05)
    if int(is_fg.sum()) < 35:
        is_fg = (pval_all <= 0.015) & (np.abs(dpsi_v) >= 0.05)
    is_bg = (fdr_all > 0.25) & (np.abs(dpsi_v) < 0.03)

    fg_df = df.loc[is_fg].copy()
    bg_pool = df.loc[is_bg & (~df["gene_id"].isin(set(fg_df["gene_id"])))].copy()
    if len(fg_df) >= 10 and len(bg_pool) >= 40:
        bg_parts = []
        for et, fg_et in fg_df.groupby("event_type"):
            bp_et = bg_pool.loc[bg_pool["event_type"].eq(et)].copy()
            if bp_et.empty:
                continue
            med_d, med_p = float(fg_et["log_depth"].median()), float(fg_et["baseline_PSI"].median())
            bp_et["_dist"] = ((bp_et["log_depth"] - med_d) / 0.85) ** 2 + ((bp_et["baseline_PSI"] - med_p) / 0.22) ** 2
            n_take = min(len(bp_et), max(250, 5 * len(fg_et)))
            bg_parts.append(bp_et.sort_values(["_dist", "event_id"], kind="mergesort").head(n_take).drop(columns=["_dist"]))
        bg_df = pd.concat(bg_parts, ignore_index=True) if bg_parts else bg_pool.head(800)
        panel_raw = pd.concat([fg_df, bg_df], ignore_index=True)
        panel_raw["_is_fg"] = [True] * len(fg_df) + [False] * len(bg_df)
    else:
        df["_h"] = df["event_id"].astype(str).map(lambda s: hashlib.md5(s.encode()).hexdigest())
        panel_raw = df.sort_values(["depth", "_h"], ascending=[False, True]).head(int(max_events_blind_depth)).drop(columns=["_h"]).copy()
        panel_raw["_is_fg"] = True

    valid_mask = []
    for rec in panel_raw.to_dict("records"):
        try:
            event_forms(rec)
            valid_mask.append(True)
        except Exception:
            valid_mask.append(False)
    cont_df = panel_raw.loc[valid_mask].reset_index(drop=True)

    rna_df = extract_event_rna_v2(cont_df, Path(fasta_path), intron_nt=200, acceptor_core_nt=3, feature_only=False, compute_sha256=False)
    rbp_tensor, rbps, region_names = build_directional_11channel_tensor_v22(cont_df, rna_df, lib)

    aff = rna_df[rna_df["layer"].eq("affected_sequence")]
    exon_gc_map, aff_len_map = {}, {}
    for eid, grp in aff.groupby("event_id", sort=False):
        joined = "".join(grp["sequence"].astype(str))
        exon_gc_map[str(eid)] = compute_sequence_gc(joined)[0]
        aff_len_map[str(eid)] = float(len(joined))
    flank_gc_ser = compute_flank_gc_from_rna(rna_df)

    ev_i = {str(e): i for i, e in enumerate(cont_df["event_id"])}
    cpg_num = np.zeros(len(cont_df), dtype=float)
    cpg_den = np.zeros(len(cont_df), dtype=float)
    for rec in rna_df[["event_id", "sequence"]].itertuples(index=False):
        idx = ev_i.get(str(rec.event_id))
        if idx is not None:
            s = str(rec.sequence)
            cpg_num[idx] += s.count("CG")
            cpg_den[idx] += max(len(s) - 1, 1)
    cpg_rate = cpg_num / np.maximum(cpg_den, 1.0)

    exon_gc = np.array([exon_gc_map.get(str(x), 0.5) for x in cont_df["event_id"]], dtype=float)
    flank_gc = np.array([float(flank_gc_ser.get(str(x), 0.5)) for x in cont_df["event_id"]], dtype=float)
    exon_gc = np.where(np.isfinite(exon_gc), exon_gc, 0.5)
    flank_gc = np.where(np.isfinite(flank_gc), flank_gc, 0.5)
    log_len = np.array([float(np.log1p(aff_len_map.get(str(x), 100.0))) for x in cont_df["event_id"]], dtype=float)
    b_psi = np.clip(cont_df["baseline_PSI"].to_numpy(dtype=float), 0.02, 0.98)
    etypes = cont_df["event_type"].astype(str).to_numpy()
    uniq_et = sorted(set(etypes))

    X_nuis = np.column_stack([
        np.ones(len(cont_df), dtype=float),
        exon_gc,
        flank_gc,
        cpg_rate,
        log_len,
        b_psi,
        b_psi * (1.0 - b_psi),
        cont_df["log_depth"].to_numpy(dtype=float),
        *[(etypes == et).astype(float) for et in uniq_et[1:]],
    ])

    dpsi_c = cont_df["dpsi"].to_numpy(dtype=float)
    pv_c = np.clip(pd.to_numeric(cont_df.get("PValue", cont_df.get("fdr", 1.0)), errors="coerce").fillna(1.0).to_numpy(dtype=float), 1e-12, 1.0)
    z_raw = np.sign(dpsi_c) * np.clip(norm.isf(0.5 * pv_c), 0.0, 8.0)
    y_sparse = np.where(cont_df["_is_fg"].to_numpy(dtype=bool), z_raw, 0.0)
    cont_df["z_shift"] = y_sparse

    genes = cont_df["gene_id"].astype(str).to_numpy()
    ug, g_codes = np.unique(genes, return_inverse=True)
    G = len(ug)
    w = 1.0 / np.bincount(g_codes)[g_codes].astype(float)
    w /= w.sum()

    n_reg, n_rbp = len(region_names), len(rbps)
    F_flat = rbp_tensor.reshape(len(cont_df), n_reg * n_rbp)
    F_res = _fwl_residuals(F_flat, X_nuis, w)
    f_sd0 = np.sqrt(np.maximum(np.sum((F_res * F_res) * w[:, None], axis=0), 1e-12))
    F_win = _fwl_residuals(np.clip(F_res / f_sd0[None, :], -3.5, 3.5), X_nuis, w)
    y_res = _fwl_residuals(y_sparse, X_nuis, w)

    F_g = np.zeros((G, n_reg * n_rbp), dtype=float)
    np.add.at(F_g, g_codes, w[:, None] * F_win)
    y_g = np.zeros(G, dtype=float)
    np.add.at(y_g, g_codes, w * y_res)

    F_gc = F_g - F_g.mean(axis=0, keepdims=True)
    yc = y_g - y_g.mean()
    sd_perm = np.sqrt(np.sum(yc * yc) * np.sum(F_gc * F_gc, axis=0) / max(G - 1, 1))
    valid_c = sd_perm > 1e-12
    Z_obs_2d, rbp_nom_p, rbp_fwer_p, reg_nom_p, reg_fwer_p = fast_gene_block_maxt_permutations(
        yc, F_gc, sd_perm, valid_c, n_reg, n_rbp, n_permutations=int(n_permutations), seed=int(seed)
    )
    reg_rows = [
        {
            "RBP": str(rbps[m_i]),
            "region": str(region_names[r_i]),
            "stat_Z": float(Z_obs_2d[r_i, m_i]),
            "nominal_P": float(reg_nom_p[r_i, m_i]),
            "region_P_FWER": float(reg_fwer_p[r_i, m_i]),
        }
        for r_i in range(n_reg)
        for m_i in range(n_rbp)
    ]
    rbp_rows = []
    for m_i, r_name in enumerate(rbps):
        best_r = int(np.argmax(np.abs(Z_obs_2d[:, m_i])))
        zv = float(Z_obs_2d[best_r, m_i])
        rbp_rows.append({
            "RBP": str(r_name),
            "best_region": str(region_names[best_r]),
            "direction": "inclusion_promoted" if zv >= 0 else "skipping_promoted",
            "stat_Z": zv,
            "abs_Z": abs(zv),
            "nominal_P": float(rbp_nom_p[m_i]),
            "rbp_P_FWER": float(rbp_fwer_p[m_i]),
            "sequence_evidence_score": float(-norm.logsf(abs(zv)) / np.log(10.0)),
        })
    rbp_df = pd.DataFrame(rbp_rows).sort_values(
        by=["rbp_P_FWER", "sequence_evidence_score", "RBP"], ascending=[True, False, True], kind="mergesort"
    ).reset_index(drop=True)
    rbp_df["sequence_rank"] = np.arange(1, len(rbp_df) + 1, dtype=int)

    sig_rbps = rbp_df.loc[rbp_df["rbp_P_FWER"] <= 0.05, "RBP"].tolist()[:10]
    if not sig_rbps:
        sig_rbps = rbp_df.loc[rbp_df["nominal_P"] <= 0.01, "RBP"].tolist()[:5]
    conseq_df = pd.DataFrame()
    if sig_rbps:
        sig_m = {m: lib.models[m] for r in sig_rbps for m in lib.rbp_to_models.get(r, ()) if m in lib.models}
        fg_only = cont_df.loc[cont_df["_is_fg"]].head(200)
        _, conseq_df = trace_motif_instances_v2(rna_df[rna_df["event_id"].isin(set(fg_only["event_id"]))], fg_only, sig_m)

    return MotifRescueResultV22(
        rbp_ranking=rbp_df,
        region_tests=pd.DataFrame(reg_rows),
        event_contrasts=cont_df,
        mature_consequences=conseq_df,
        diagnostics={
            "status": "complete",
            "inference_mode": "case_control_cpg_gene_block_maxT",
            "evaluated_events": len(cont_df),
            "foreground_events": int(cont_df["_is_fg"].sum()),
            "eligible_genes": int(G),
            "significant_rbps_fwer05": int((rbp_df["rbp_P_FWER"] <= 0.05).sum()),
        },
    )
