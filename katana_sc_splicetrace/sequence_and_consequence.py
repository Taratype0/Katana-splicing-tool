"""Single-cell RNA form reconstruction, mature-mRNA motif consequence tracing, and 11-channel positional tensors."""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import pandas as pd

from katana_sc_splicetrace.library import (
    CuratedMotifLibrary,
    load_bundled_human_rbp_library,
    load_curated_human_rbp_library,
)
from katana_sc_splicetrace.fast_accel import fast_score_mean_log_odds
from katana_sc_splicetrace.rna_forms import (
    compute_flank_gc_from_rna,
    compute_sequence_gc,
    extract_event_rna_v2,
    trace_motif_instances_v2,
)


def load_default_curated_library_v2(
    models_json: Path | None = None,
    assignments_tsv: Path | None = None,
    *,
    species: str = "human",
    extra_motif_files: list[Path | str] | None = None,
    include_v3_expansion: bool = True,
) -> CuratedMotifLibrary:
    if models_json is not None and assignments_tsv is not None:
        return load_curated_human_rbp_library(Path(models_json), Path(assignments_tsv), species=species, extra_motif_files=extra_motif_files, include_v3_expansion=include_v3_expansion)
    return load_bundled_human_rbp_library(species=species, extra_motif_files=extra_motif_files, include_v3_expansion=include_v3_expansion)


def extract_sc_event_sequences_v2(
    events_df: pd.DataFrame,
    fasta_path: Path,
    *,
    anchor_nt: int = 100,
    intron_nt: int = 200,
    distal_intron_nt: int = 500,
    acceptor_core_nt: int = 3,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Extract mature (`affected_sequence`, `mature_local_form`) and cis-intronic (`regulatory_pre_mRNA`) sequences
    plus 5 sequence composition covariates (`exon_GC`, `flank_GC`, `log_length`, `event_type`) for single-cell events.
    Defaults to `acceptor_core_nt=3` so the [-20..-3] nt 3'ss polypyrimidine tract is retained!
    """
    ev = events_df.copy()
    if "gene_symbol" not in ev.columns and "gene_name" in ev.columns:
        ev["gene_symbol"] = ev["gene_name"]
    if "comparison_id" not in ev.columns:
        ev["comparison_id"] = "sc_panel"

    rna_df = extract_event_rna_v2(
        ev,
        Path(fasta_path),
        anchor_nt=anchor_nt,
        intron_nt=intron_nt,
        distal_intron_nt=distal_intron_nt,
        acceptor_core_nt=acceptor_core_nt,
        feature_only=False, compute_sha256=False,
    )

    # Compute per-event affected_sequence length & GC + deduplicated regulatory flank GC
    aff = rna_df[rna_df["layer"].eq("affected_sequence")]
    exon_gc_map: dict[str, float] = {}
    aff_len_map: dict[str, float] = {}
    for eid, grp in aff.groupby("event_id", sort=False):
        joined = "".join(grp["sequence"].astype(str))
        gc_val, _ = compute_sequence_gc(joined)
        exon_gc_map[str(eid)] = gc_val
        aff_len_map[str(eid)] = float(len(joined))

    flank_gc_ser = compute_flank_gc_from_rna(rna_df)

    cov_df = pd.DataFrame({
        "event_id": ev["event_id"].astype(str).to_numpy(),
        "gene_id": ev["gene_id"].astype(str).to_numpy(),
        "event_type": ev["event_type"].astype(str).to_numpy(),
        "exon_GC": [exon_gc_map.get(str(x), 0.5) for x in ev["event_id"]],
        "flank_GC": [float(flank_gc_ser.get(str(x), 0.5)) for x in ev["event_id"]],
        "log_length": [float(np.log1p(aff_len_map.get(str(x), 100.0))) for x in ev["event_id"]],
    })
    cov_df["exon_GC"] = cov_df["exon_GC"].fillna(cov_df["exon_GC"].median() if cov_df["exon_GC"].notna().any() else 0.5)
    cov_df["flank_GC"] = cov_df["flank_GC"].fillna(cov_df["flank_GC"].median() if cov_df["flank_GC"].notna().any() else 0.5)
    return rna_df, cov_df


def score_mean_log_odds_per_form_v2(sequences: list[str], models: dict[str, dict]) -> np.ndarray:
    return fast_score_mean_log_odds(sequences, models)


def build_sc_directional_motif_tensor_v2(
    events_df: pd.DataFrame,
    rna_df: pd.DataFrame,
    library: CuratedMotifLibrary | None = None,
    *,
    include_exonic_boundaries: bool = True,
    include_distal_flanks: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    """Build signed directional x positional RBP motif tensor (E x R_regions x N_RBPs) combining:
      0. `mature_consequence_delta`: log2(1+odds_inclusion) - log2(1+odds_exclusion) on mature mRNA (`mature_local_form`)
      1..3. `upstream_intron_[0_35, 35_85, 85_200]`: 3'ss-adjacent intronic distance bins (including [-20..-3] nt PPT!)
      4..6. `downstream_intron_[0_35, 35_85, 85_200]`: 5'ss-adjacent intronic distance bins
      (+ 4 exonic ESE/ESS 50-nt splice-junction boundary bins when `include_exonic_boundaries=True`)
      (+ 2 distal [200..500] nt flanks when `include_distal_flanks=True`).
    """
    lib = load_default_curated_library_v2() if library is None else library
    models = lib.models
    m_ids = list(models.keys())
    mid_to_idx = {mid: idx for idx, mid in enumerate(m_ids)}
    rbps = sorted(lib.rbp_to_models.keys())
    rbp_to_mcols: dict[str, list[int]] = {
        r: [mid_to_idx[mid] for mid in lib.rbp_to_models[r] if mid in mid_to_idx]
        for r in rbps
    }

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
    ]
    exonic_offset = len(region_names)
    if include_exonic_boundaries:
        region_names.extend([
            "upstream_exon_3p_50",
            "alt_exon_5p_50",
            "alt_exon_3p_50",
            "downstream_exon_5p_50",
        ])
    distal_offset = len(region_names)
    if include_distal_flanks:
        region_names.extend(["upstream_intron_distal", "downstream_intron_distal"])

    model_tensor = np.zeros((n_ev, len(region_names), len(m_ids)), dtype=float)

    # 1. Mature local forms: inclusion minus exclusion + 50-nt ESE/ESS exonic boundaries around splice joins
    mature = rna_df[rna_df["layer"].eq("mature_local_form") & rna_df["region"].isin(["inclusion", "exclusion"])].copy()
    if not mature.empty:
        m_scores = score_mean_log_odds_per_form_v2(mature["sequence"].astype(str).tolist(), models)
        inc_mat = np.zeros((n_ev, len(m_ids)), dtype=float)
        exc_mat = np.zeros((n_ev, len(m_ids)), dtype=float)
        ex_seqs: list[str] = []
        ex_targets: list[tuple[int, int]] = []
        for row_i, rec in enumerate(mature[["event_id", "region", "sequence", "junctions_json"]].itertuples(index=False)):
            ei = ev_index.get(str(rec.event_id))
            if ei is None:
                continue
            if rec.region == "inclusion":
                inc_mat[ei] = m_scores[row_i]
                if include_exonic_boundaries:
                    seq_inc = str(rec.sequence)
                    joins = json.loads(str(rec.junctions_json)) if rec.junctions_json else []
                    if len(joins) >= 2:
                        j1, j2 = int(joins[0]), int(joins[-1])
                        segs = (
                            (seq_inc[max(0, j1 - 50) : j1], exonic_offset + 0),
                            (seq_inc[j1 : min(j2, j1 + 50)], exonic_offset + 1),
                            (seq_inc[max(j1, j2 - 50) : j2], exonic_offset + 2),
                            (seq_inc[j2 : min(len(seq_inc), j2 + 50)], exonic_offset + 3),
                        )
                        for piece, sl_i in segs:
                            if len(piece) >= 6:
                                ex_seqs.append(piece)
                                ex_targets.append((ei, sl_i))
                    elif len(joins) == 1:
                        j1 = int(joins[0])
                        for piece, sl_i in (
                            (seq_inc[max(0, j1 - 50) : j1], exonic_offset + 0),
                            (seq_inc[j1 : min(len(seq_inc), j1 + 50)], exonic_offset + 3),
                        ):
                            if len(piece) >= 6:
                                ex_seqs.append(piece)
                                ex_targets.append((ei, sl_i))
            elif rec.region == "exclusion":
                exc_mat[ei] = m_scores[row_i]
        model_tensor[:, 0, :] = inc_mat - exc_mat
        if ex_seqs:
            ex_scores = score_mean_log_odds_per_form_v2(ex_seqs, models)
            for idx_e, (ei, sl_i) in enumerate(ex_targets):
                model_tensor[ei, sl_i, :] = ex_scores[idx_e]

    # 2. Regulatory intronic flanks resolved by exact alternative splice-site geometry (no 2:1 constitutive dilution!)
    reg_df = rna_df[rna_df["layer"].eq("regulatory_pre_mRNA")].copy()
    if not reg_df.empty:
        etype_map = {str(eid): str(et) for eid, et in zip(events_df["event_id"], events_df["event_type"])}
        sub_seqs: list[str] = []
        # Each target tuple: (event_idx, slot_idx, sign_weight)
        sub_targets: list[tuple[int, int, float]] = []
        for rec in reg_df[["event_id", "region", "sequence"]].itertuples(index=False):
            ei = ev_index.get(str(rec.event_id))
            if ei is None:
                continue
            rname = str(rec.region)
            seq = str(rec.sequence)
            if not seq:
                continue
            kind = etype_map.get(str(rec.event_id), "SE")
            is_distal = rname.endswith("_distal")
            base_r = rname[: -len("_distal")] if is_distal else rname

            # Route each splice-junction flank to (is_3ss_acceptor_slot, sign_weight) without constitutive 2:1 dilution
            route: tuple[bool, float] | None = None
            if kind == "SE":
                if base_r == "inclusion_junction_1_acceptor_intron":
                    route = (True, 1.0)   # Pure cassette exon 3'ss (incl. [-20..-3] PPT!)
                elif base_r == "inclusion_junction_2_donor_intron":
                    route = (False, 1.0)  # Pure cassette exon 5'ss (RBFOX2 / TIA1 / HNRNPL / QKI!)
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
            else:
                is_up_fb = any(tok in base_r for tok in ("upstream", "acceptor_intron", "intron1"))
                route = (is_up_fb, 1.0)

            if route is None:
                continue
            is_up, sgn = route

            if is_distal:
                if include_distal_flanks:
                    sub_seqs.append(seq)
                    sub_targets.append((ei, distal_offset + (0 if is_up else 1), sgn))
                continue

            L = len(seq)
            bins_nt = ((0, 35), (35, 85), (85, 200))
            for b_i, (d0, d1) in enumerate(bins_nt):
                if is_up:
                    s_end = max(0, L - d0)
                    s_beg = max(0, L - d1)
                    piece = seq[s_beg:s_end]
                else:
                    piece = seq[min(L, d0) : min(L, d1)]
                if len(piece) >= 6:
                    sub_seqs.append(piece)
                    sub_targets.append((ei, (1 if is_up else 4) + b_i, sgn))

        if sub_seqs:
            b_scores = score_mean_log_odds_per_form_v2(sub_seqs, models)
            for idx_s, (ei, slot, sgn) in enumerate(sub_targets):
                model_tensor[ei, slot, :] += sgn * b_scores[idx_s]

    # Standardize each model column across events before averaging into archetype RBP features
    rms_m = np.sqrt(np.mean(model_tensor * model_tensor, axis=0, keepdims=True))
    norm_model_tensor = np.divide(model_tensor, rms_m, out=np.zeros_like(model_tensor), where=rms_m > 1e-12)

    rbp_tensor = np.zeros((n_ev, len(region_names), len(rbps)), dtype=float)
    for r_i, r_name in enumerate(rbps):
        cols = rbp_to_mcols[r_name]
        if cols:
            rbp_tensor[:, :, r_i] = np.mean(norm_model_tensor[:, :, cols], axis=2)

    return rbp_tensor, event_ids, np.asarray(rbps, dtype=str), region_names


def trace_sc_mature_motif_consequences_v2(
    rna_df: pd.DataFrame,
    contrast_events_df: pd.DataFrame,
    library: CuratedMotifLibrary | None = None,
    *,
    top_n_events: int = 150,
    min_abs_z: float = 1.5,
    relative_likelihood_cutoff: float = 0.80,
) -> pd.DataFrame:
    """Trace exact skipped (`inclusion_only`) and de novo splice-junction (`exclusion_only`, `crosses_splice_join=True`)
    mature-mRNA RBP motif instances on the top shifted single-cell splicing events.
    """
    lib = load_default_curated_library_v2() if library is None else library
    ev = contrast_events_df.copy()
    if "z_shift" in ev.columns:
        ev = ev.loc[ev["z_shift"].abs() >= float(min_abs_z)].copy()
        ev = ev.reindex(ev["z_shift"].abs().sort_values(ascending=False).index).head(int(top_n_events))
    elif "dpsi" in ev.columns:
        ev = ev.reindex(ev["dpsi"].abs().sort_values(ascending=False).index).head(int(top_n_events))
    if ev.empty:
        return pd.DataFrame()

    sub_rna = rna_df[rna_df["event_id"].isin(set(ev["event_id"].astype(str)))].copy()
    _, changes_df = trace_motif_instances_v2(
        sub_rna,
        ev,
        lib.models,
        relative_likelihood_cutoff=relative_likelihood_cutoff,
    )
    if changes_df.empty:
        return changes_df

    z_map = ev.set_index("event_id")["z_shift"].to_dict() if "z_shift" in ev.columns else {}
    dpsi_map = ev.set_index("event_id")["dpsi"].to_dict() if "dpsi" in ev.columns else {}
    changes_df["sc_contrast_z"] = [float(z_map.get(str(x), np.nan)) for x in changes_df["event_id"]]
    changes_df["sc_shrunk_dpsi"] = [float(dpsi_map.get(str(x), np.nan)) for x in changes_df["event_id"]]
    return changes_df
