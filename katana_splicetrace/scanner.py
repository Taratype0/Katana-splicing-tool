"""Vectorized width-batched MotifTrace v2 scan engine."""
from __future__ import annotations

import numpy as np
import pandas as pd
from katana_splicetrace.inference import BlockFamilyV2
from katana_splicetrace.library import CuratedMotifLibrary
import katana_splicetrace.fast_accel as _fa

_ENCODE_TABLE = np.full(256, 4, dtype=np.int16)
_ENCODE_TABLE[np.frombuffer(b"ACGT", dtype=np.uint8)] = np.arange(4)


def batch_pwm_profiles_for_width(
    sequences: dict[str, str],
    odds_batch: np.ndarray,
    w: int,
    widths: tuple[int, ...] = (0, 10, 20, 50),
    *,
    reverse_positional: bool = False,
) -> dict[str, dict[int, np.ndarray]]:
    """Compute exact PWM profiles simultaneously across all M_w models of width w."""
    m_count = len(odds_batch)
    result: dict[str, dict[int, np.ndarray]] = {}
    for ev_id, raw_seq in sequences.items():
        seq = raw_seq.upper().replace("U", "T")
        seq_len = len(seq)
        n = max(0, seq_len - w + 1)
        ev_out: dict[int, np.ndarray] = {}
        if n == 0:
            for width in widths:
                span = 1 if width == 0 else max(0, seq_len - width + 1)
                ev_out[width] = np.full((m_count, span), np.nan, dtype=float)
            result[ev_id] = ev_out
            continue

        encoded = _ENCODE_TABLE[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]
        invalid = np.r_[0, np.cumsum(encoded == 4)]
        valid = (invalid[w:] - invalid[:-w]) == 0
        starts = np.flatnonzero(valid)

        if _fa.HAS_NUMBA and n >= 8:
            pos_batch, valid_nb = _fa._numba_positional_odds_matrix(encoded.astype(np.uint8), 0, seq_len, np.ascontiguousarray(odds_batch, dtype=np.float64), int(w))
            valid = valid_nb.astype(bool)
            starts = np.flatnonzero(valid)
            vals = pos_batch[:, starts] if len(starts) else np.empty((m_count, 0), dtype=float)
        else:
            pos_batch = np.zeros((m_count, n), dtype=float)
            if len(starts):
                vals = np.ones((m_count, len(starts)), dtype=float)
                for k in range(w):
                    vals *= odds_batch[:, k, encoded[starts + k]]
                pos_batch[:, starts] = vals
            else:
                vals = np.empty((m_count, 0), dtype=float)

        csum = np.concatenate([np.zeros((m_count, 1), dtype=float), np.cumsum(pos_batch, axis=1)], axis=1)
        ccnt = np.concatenate([np.zeros(1, dtype=int), np.cumsum(valid.astype(int))])

        for width in widths:
            if width == 0:
                if len(starts):
                    ev_out[0] = np.log2(1.0 + np.mean(vals, axis=1))[:, None]
                else:
                    ev_out[0] = np.full((m_count, 1), np.nan, dtype=float)
                continue
            span = max(0, seq_len - width + 1)
            k = width - w + 1
            prof = np.full((m_count, span), np.nan, dtype=float)
            if k > 0 and span > 0:
                w_sum = csum[:, k : k + span] - csum[:, :span]
                w_cnt = ccnt[k : k + span] - ccnt[:span]
                ok_w = w_cnt > 0
                if ok_w.any():
                    prof[:, ok_w] = np.log2(1.0 + w_sum[:, ok_w] / w_cnt[None, ok_w])
                if reverse_positional:
                    prof = prof[:, ::-1]
            ev_out[width] = prof
        result[ev_id] = ev_out
    return result


def run_motiftrace_v2_scan(
    matched_pairs: pd.DataFrame,
    rna_df: pd.DataFrame,
    library: CuratedMotifLibrary,
    *,
    permutations: int = 999,
    seed: int = 20261015,
    minimum_blocks: int = 20,
    use_dpsi_weights: bool = True,
    scale_normalize_clusters: bool = True,
    widths: tuple[int, ...] = (0, 10, 20, 50),
) -> dict[str, pd.DataFrame]:
    """Execute full MotifTrace v2 scan returning whole_regions, clusters, and rbp_ranking tables."""
    pairs = matched_pairs[matched_pairs["status"].eq("matched")].copy()
    if pairs.empty:
        return {
            "whole_regions": pd.DataFrame(),
            "clusters": pd.DataFrame(),
            "rbp_ranking_all": pd.DataFrame(),
            "rbp_ranking_affected": pd.DataFrame(),
            "rbp_ranking_regulatory": pd.DataFrame(),
        }

    pairs["block"] = pairs["foreground_gene"].astype(str)
    family = BlockFamilyV2(
        pairs["block"].tolist(),
        permutations=permutations,
        seed=seed,
        minimum_blocks=minimum_blocks,
        scale_normalize_clusters=scale_normalize_clusters,
    )

    regions_df = rna_df[rna_df["layer"].isin(["affected_sequence", "regulatory_pre_mRNA"])]
    seq_lookup = {(str(r.event_id), str(r.region)): str(r.sequence) for r in regions_df.itertuples(index=False)}

    # Group library PWMs by motif width
    by_width: dict[int, tuple[list[str], np.ndarray]] = {}
    for w_val in sorted({len(m["matrix"]) for m in library.models.values()}):
        mids = sorted([mid for mid, m in library.models.items() if len(m["matrix"]) == w_val])
        odds = np.stack(
            [(20.0 * np.asarray(library.models[mid]["matrix"], dtype=float) + 0.025) / 20.1 / 0.25 for mid in mids],
            axis=0,
        )
        by_width[w_val] = (mids, odds)

    whole_frames: list[pd.DataFrame] = []
    all_clusters: list[dict] = []

    for kind in sorted(regions_df["event_type"].unique()):
        kind_regions = sorted(regions_df.loc[regions_df["event_type"].eq(kind), "region"].unique())
        for direction in ("inclusion_down", "inclusion_up"):
            grp = pairs[pairs["event_type"].eq(kind) & pairs["direction"].eq(direction)].copy()
            if grp["block"].nunique() < minimum_blocks:
                continue
            fg_evs = grp["foreground_event"].astype(str).tolist()
            bg_evs = grp["background_event"].astype(str).tolist()
            blocks_list = grp["block"].tolist()
            if use_dpsi_weights and "dpsi" in grp.columns and grp["dpsi"].notna().all():
                # Bounded |dPSI| weighting (clipped to [0.10, 0.80]) invariant under sign flips
                p_weights = np.clip(grp["dpsi"].abs().to_numpy(dtype=float), 0.10, 0.80)
            else:
                p_weights = None

            needed_evs = sorted(set(fg_evs) | set(bg_evs))
            for region in kind_regions:
                seqs = {ev: seq_lookup.get((ev, region), "") for ev in needed_evs}
                layer = "affected_sequence" if region.startswith("affected_") else "regulatory_pre_mRNA"
                rev_pos = ("acceptor_intron" in region)

                # Precompute estimable spans per width (depends only on sequence lengths and width, not PWM!)
                spans: dict[int, int] = {}
                for width in widths:
                    if width == 0:
                        spans[0] = 1
                    else:
                        lens = np.array([
                            min(max(0, len(seqs[fe]) - width + 1), max(0, len(seqs[be]) - width + 1))
                            for fe, be in zip(fg_evs, bg_evs)
                        ])
                        gl = pd.DataFrame({"gene": blocks_list, "length": lens}).groupby("gene")["length"].max().sort_values(ascending=False)
                        spans[width] = int(gl.iloc[minimum_blocks - 1]) if len(gl) >= minimum_blocks else 0

                for w_val, (mids, odds_batch) in by_width.items():
                    prof_map = batch_pwm_profiles_for_width(seqs, odds_batch, w_val, widths, reverse_positional=rev_pos)
                    for width in widths:
                        span = spans[width]
                        if span <= 0:
                            continue
                        diff_3d = np.full((len(mids), len(grp), span), np.nan, dtype=float)
                        for pair_i, (fe, be) in enumerate(zip(fg_evs, bg_evs)):
                            a_mat = prof_map[fe][width]
                            b_mat = prof_map[be][width]
                            n_common = min(a_mat.shape[1], b_mat.shape[1], span)
                            if n_common > 0:
                                diff_3d[:, pair_i, :n_common] = a_mat[:, :n_common] - b_mat[:, :n_common]
                        rbp_labels = [library.models[mid].get("RBPs", "") for mid in mids]
                        base_meta = {
                            "event_type": kind,
                            "direction": direction,
                            "region": region,
                            "window_nt": width,
                            "layer": layer,
                        }
                        rows, clust = family.evaluate_model_batch(
                            diff_3d, blocks_list, mids, rbp_labels, base_meta, pair_weights=p_weights
                        )
                        if width == 0 and not rows.empty:
                            valid_rows = rows[rows["raw_P"].notna()]
                            if not valid_rows.empty:
                                whole_frames.append(valid_rows)
                        elif width > 0:
                            all_clusters.extend(clust)

    whole_raw = pd.concat(whole_frames, ignore_index=True) if whole_frames else pd.DataFrame()
    whole_adj = family.adjust_whole_regions(whole_raw)
    clusters_adj = family.adjust_clusters(all_clusters)
    return {
        "whole_regions": whole_adj,
        "clusters": clusters_adj,
        "rbp_ranking_all": family.summarize_rbp_pivotal(library.rbp_to_models, layer="all"),
        "rbp_ranking_affected": family.summarize_rbp_pivotal(library.rbp_to_models, layer="affected_sequence"),
        "rbp_ranking_regulatory": family.summarize_rbp_pivotal(library.rbp_to_models, layer="regulatory_pre_mRNA"),
    }
