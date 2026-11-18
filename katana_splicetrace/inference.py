"""MotifTrace v2 gene-block randomization with dual exact FWER & symmetric multi-PWM pivotal calibration.

Key mathematical upgrades over v1.1.0:
1. Dedicated Whole-Region Exact Westfall-Young Max-Z FWER (`window_nt == 0`):
   eliminates the 1.23M sliding-window BH dilution bug and rescues diffuse
   exonic/intronic regulators (e.g. SRSF1) with exact finite-sample FWER control.
2. Scale-Normalized Positional Cluster Mass (`mass / sqrt(window_nt)`):
   puts 10-nt, 20-nt, and 50-nt sliding-window clusters on equal variance footing
   under Westfall-Young max-T randomization.
3. Optional |dPSI| effect-magnitude pair weighting (`pair_weights`):
   upweights strong biological switches while remaining strictly invariant under
   gene-block sign flips in {-1, +1}, preserving exact conditional validity.
4. Exchangeable Pooled-Survival Pivotal RBP Family Summary (`summarize_rbp_pivotal`):
   maps per-model global whole-Z and cluster-mass maxima symmetrically across all
   (B + 1) draws onto empirical log-survival scores, standardizes each RBP's
   multi-PWM maximum symmetrically across {0..B}, and computes exact finite-sample
   Westfall-Young RBP-level FWER P-values immune to sparse-cluster zero inflation.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.ndimage import label
from scipy.stats import norm, rankdata

THRESHOLD = float(norm.isf(0.01))


def bh(p: np.ndarray) -> np.ndarray:
    p_arr = np.asarray(p, dtype=float)
    out = np.full(len(p_arr), np.nan)
    valid = np.flatnonzero(np.isfinite(p_arr))
    if len(valid):
        order = valid[np.argsort(p_arr[valid], kind="stable")]
        q = p_arr[order] * len(order) / np.arange(1, len(order) + 1)
        out[order] = np.minimum(1.0, np.minimum.accumulate(q[::-1])[::-1])
    return out


def _fast_average_ranks_2d(values: np.ndarray) -> np.ndarray:
    """Exact C-vectorized columnwise average ranks omitting NaNs (matches scipy.stats.rankdata)."""
    n_rows, n_cols = values.shape
    if n_rows == 0 or n_cols == 0:
        return np.zeros_like(values, dtype=float)
    fin = np.isfinite(values)
    v_inf = np.where(fin, values, np.inf)
    order = np.argsort(v_inf, axis=0, kind="mergesort")
    sorted_v = np.take_along_axis(v_inf, order, axis=0)
    sorted_fin = np.isfinite(sorted_v)

    base_ranks = np.where(sorted_fin, (np.arange(1, n_rows + 1, dtype=float))[:, None], 0.0)
    out = np.zeros((n_rows, n_cols), dtype=float)
    np.put_along_axis(out, order, base_ranks, axis=0)

    if n_rows > 1:
        tie_cols = np.flatnonzero(((sorted_v[1:] == sorted_v[:-1]) & sorted_fin[1:]).any(axis=0))
        if len(tie_cols):
            out[:, tie_cols] = np.nan_to_num(
                rankdata(values[:, tie_cols], axis=0, method="average", nan_policy="omit")
            )
    return out


def signed_block_ranks_v2(
    differences: np.ndarray,
    blocks: list[str] | np.ndarray,
    *,
    pair_weights: np.ndarray | None = None,
    minimum_blocks: int = 20,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute standardized gene-block signed rank vectors and observed Z statistics."""
    d = np.asarray(differences, dtype=float)
    if d.ndim != 2 or len(blocks) != len(d):
        raise ValueError("Event difference matrix and block identities disagree")
    if pair_weights is not None:
        w = np.asarray(pair_weights, dtype=float)
        if w.shape != (len(d),) or not np.isfinite(w).all() or (w < 0).any():
            raise ValueError("pair_weights must be finite non-negative 1D array matching events")
    else:
        w = None

    keys, codes = np.unique(np.asarray(blocks, dtype=str), return_inverse=True)
    valid = np.isfinite(d)
    n_events = valid.sum(axis=0)
    values = np.where(valid & (d != 0.0), np.abs(d), np.nan)
    ranks = _fast_average_ranks_2d(values) * np.sign(np.nan_to_num(d))
    if w is not None:
        ranks = ranks * w[:, None]

    sums = np.zeros((len(keys), d.shape[1]), dtype=float)
    presence = np.zeros((len(keys), d.shape[1]), dtype=np.int32)
    np.add.at(sums, codes, ranks)
    np.add.at(presence, codes, valid.astype(np.int32))
    n_blocks = (presence > 0).sum(axis=0)

    denom = np.sqrt(np.square(sums).sum(axis=0))
    denom[denom == 0.0] = 1.0
    normalized = sums / denom
    observed = normalized.sum(axis=0)
    observed[n_blocks < minimum_blocks] = np.nan
    return keys, normalized, observed, n_events, n_blocks


def contiguous_runs(z: np.ndarray, threshold: float = THRESHOLD) -> list[tuple[int, int, float]]:
    idx = np.flatnonzero(np.isfinite(z) & (z > threshold))
    if not len(idx):
        return []
    splits = np.split(idx, np.flatnonzero(np.diff(idx) != 1) + 1)
    return [(int(p[0]), int(p[-1] + 1), float((z[p] - threshold).sum())) for p in splits]


def max_cluster_rows(z: np.ndarray, threshold: float = THRESHOLD) -> np.ndarray:
    mask = np.isfinite(z) & (z > threshold)
    result = np.zeros(len(z), dtype=float)
    active = np.flatnonzero(mask.any(axis=1))
    if not len(active):
        return result
    sub_mask = mask[active]
    sub_z = z[active]
    labs, n_comp = label(sub_mask, structure=np.array([[0, 0, 0], [1, 1, 1], [0, 0, 0]]))
    if n_comp:
        rr, cc = np.nonzero(sub_mask)
        if len(rr):
            masses = np.bincount(labs[rr, cc], weights=sub_z[rr, cc] - threshold, minlength=n_comp + 1)
            comp_vals = masses[labs[rr, cc]]
            row_starts = np.r_[0, np.flatnonzero(rr[1:] != rr[:-1]) + 1]
            sub_res = np.zeros(len(active), dtype=float)
            sub_res[rr[row_starts]] = np.maximum.reduceat(comp_vals, row_starts)
            result[active] = sub_res
    return result


def max_cluster_rows_3d(z3: np.ndarray, threshold: float = THRESHOLD) -> np.ndarray:
    b_len, m_len, span = z3.shape
    flat = z3.reshape(b_len * m_len, span)
    return max_cluster_rows(flat, threshold=threshold).reshape(b_len, m_len)


def _whole_z_log_survival(mat: np.ndarray) -> np.ndarray:
    """Exact tail-preserving -log10(Phi(-Z)) score without finite-sample rank ceiling compression."""
    out = np.zeros_like(mat, dtype=float)
    fin = np.isfinite(mat)
    if fin.any():
        out[fin] = -norm.logsf(mat[fin]) / np.log(10.0)
    return out


def _cluster_pareto_log_survival(mat: np.ndarray) -> np.ndarray:
    """Symmetric {0..B}-exchangeable empirical + exponential-Pareto upper-tail log-survival for cluster mass."""
    flat = mat.ravel()
    out = np.zeros_like(flat, dtype=float)
    pos = np.isfinite(flat) & (flat > 0.0)
    if not pos.any():
        return out.reshape(mat.shape)
    n_tot = int(np.isfinite(flat).sum())
    vals = flat[pos]
    ranks_desc = rankdata(-vals, method="max")
    surv = ranks_desc / float(max(1, n_tot))
    scores = -np.log10(np.clip(surv, 1e-15, 1.0))

    # Smooth symmetric exponential-Pareto continuation above 99th percentile of positive pooled entries
    if len(vals) >= 200:
        q99 = float(np.quantile(vals, 0.99))
        tail_mask = vals > q99
        if int(tail_mask.sum()) >= 10:
            beta = float(np.mean(vals[tail_mask] - q99))
            if beta > 1e-6:
                p_q99 = float(np.sum(vals >= q99)) / float(max(1, n_tot))
                base_s = -np.log10(max(p_q99, 1e-15))
                scores[tail_mask] = base_s + (vals[tail_mask] - q99) / (beta * np.log(10.0))

    out[pos] = scores
    return out.reshape(mat.shape)


class BlockFamilyV2:
    """Shared gene-block sign-flip inference tracking Whole-Region, Cluster, and Exchangeable RBP Pivotal nulls."""

    def __init__(
        self,
        blocks: list[str] | np.ndarray,
        *,
        permutations: int = 999,
        seed: int = 20261015,
        minimum_blocks: int = 20,
        scale_normalize_clusters: bool = True,
    ) -> None:
        if permutations < 99:
            raise ValueError("At least 99 randomizations required")
        self.keys = np.array(sorted(set(map(str, blocks))))
        self.lookup = {x: i for i, x in enumerate(self.keys)}
        self.signs = np.random.default_rng(seed).choice(np.array([-1.0, 1.0]), size=(permutations, len(self.keys)))
        self.B = int(permutations)
        self.minimum = int(minimum_blocks)
        self.scale_normalize = bool(scale_normalize_clusters)

        self.whole_null_max: dict[str, np.ndarray] = {
            k: np.full(self.B, -np.inf) for k in ("all", "affected_sequence", "regulatory_pre_mRNA")
        }
        self.cluster_null_max: dict[str, np.ndarray] = {
            k: np.zeros(self.B, dtype=float) for k in ("all", "affected_sequence", "regulatory_pre_mRNA")
        }
        # Global per-model maxima across strata over all (B + 1) draws (index 0 = observed, 1..B = null)
        self._model_whole_max: dict[tuple[str, str], np.ndarray] = {}
        self._model_cluster_max: dict[tuple[str, str], np.ndarray] = {}

    def _record_model_whole(self, layer: str, model_id: str, obs_val: float, null_vals: np.ndarray) -> None:
        if not np.isfinite(obs_val):
            return
        vec = np.empty(self.B + 1, dtype=float)
        vec[0] = obs_val
        vec[1:] = np.where(np.isfinite(null_vals), null_vals, -np.inf)
        for lyr in ("all", layer):
            key = (lyr, model_id)
            cur = self._model_whole_max.get(key)
            if cur is None:
                self._model_whole_max[key] = vec.copy()
            else:
                np.maximum(cur, vec, out=cur)

    def _record_model_cluster(self, layer: str, model_id: str, obs_val: float, null_vals: np.ndarray) -> None:
        if not np.isfinite(obs_val):
            return
        vec = np.empty(self.B + 1, dtype=float)
        vec[0] = obs_val
        vec[1:] = np.nan_to_num(null_vals, nan=0.0)
        for lyr in ("all", layer):
            key = (lyr, model_id)
            cur = self._model_cluster_max.get(key)
            if cur is None:
                self._model_cluster_max[key] = vec.copy()
            else:
                np.maximum(cur, vec, out=cur)

    def evaluate(
        self,
        differences: np.ndarray,
        blocks: list[str] | np.ndarray,
        meta: dict,
        *,
        pair_weights: np.ndarray | None = None,
    ) -> tuple[pd.DataFrame, list[dict]]:
        diff_3d = np.asarray(differences, dtype=float)[None, :, :]
        mid = str(meta.get("model_id", ""))
        rbp = str(meta.get("RBP", ""))
        base_meta = {k: v for k, v in meta.items() if k not in {"model_id", "RBP"}}
        return self.evaluate_model_batch(
            diff_3d, blocks, [mid], [rbp], base_meta, pair_weights=pair_weights
        )

    def evaluate_model_batch(
        self,
        diff_3d: np.ndarray,
        blocks: list[str] | np.ndarray,
        model_ids: list[str],
        rbp_labels: list[str],
        base_meta: dict,
        *,
        pair_weights: np.ndarray | None = None,
        max_models_per_subbatch: int = 36,
    ) -> tuple[pd.DataFrame, list[dict]]:
        """Evaluate M models simultaneously with L3-cache-friendly sub-batching."""
        m_total = len(model_ids)
        if m_total > max_models_per_subbatch and int(base_meta.get("window_nt", 0)) > 0:
            r_parts: list[pd.DataFrame] = []
            c_parts: list[dict] = []
            for st in range(0, m_total, max_models_per_subbatch):
                ed = min(m_total, st + max_models_per_subbatch)
                rp, cp = self.evaluate_model_batch(
                    diff_3d[st:ed],
                    blocks,
                    model_ids[st:ed],
                    rbp_labels[st:ed],
                    base_meta,
                    pair_weights=pair_weights,
                    max_models_per_subbatch=max_models_per_subbatch,
                )
                if not rp.empty:
                    r_parts.append(rp)
                c_parts.extend(cp)
            return (pd.concat(r_parts, ignore_index=True) if r_parts else pd.DataFrame()), c_parts

        m_count, n_pairs, span = diff_3d.shape
        flat_diff = diff_3d.transpose(1, 0, 2).reshape(n_pairs, m_count * span)
        keys, ranks, z_flat, n_ev_flat, n_bg_flat = signed_block_ranks_v2(
            flat_diff, blocks, pair_weights=pair_weights, minimum_blocks=self.minimum
        )
        indices = np.array([self.lookup[x] for x in keys], dtype=int)
        signs = self.signs[:, indices]
        good_flat = np.isfinite(z_flat)
        exceed_flat = np.zeros(m_count * span, dtype=np.int32)

        width = int(base_meta.get("window_nt", 0))
        layer = str(base_meta.get("layer", "affected_sequence"))
        scale_div = np.sqrt(float(width)) if (width > 0 and self.scale_normalize) else 1.0

        track_null = np.zeros((self.B, m_count), dtype=float) if width > 0 else np.full((self.B, m_count), -np.inf)
        chunk_sz = 512
        for start in range(0, self.B, chunk_sz):
            end = min(self.B, start + chunk_sz)
            null_flat = signs[start:end] @ ranks
            null_flat[:, ~good_flat] = np.nan
            exceed_flat += np.sum(null_flat >= z_flat - 1e-12, axis=0)
            null_3d = null_flat.reshape(end - start, m_count, span)

            if width == 0:
                vals_0 = null_3d[:, :, 0]
                track_null[start:end] = vals_0
                valid_m = np.isfinite(vals_0).any(axis=0)
                if valid_m.any():
                    row_max = np.nanmax(vals_0[:, valid_m], axis=1)
                    self.whole_null_max["all"][start:end] = np.maximum(self.whole_null_max["all"][start:end], row_max)
                    if layer in self.whole_null_max:
                        self.whole_null_max[layer][start:end] = np.maximum(self.whole_null_max[layer][start:end], row_max)
            else:
                cm_bm = max_cluster_rows_3d(null_3d) / scale_div
                track_null[start:end] = cm_bm
                row_max = np.max(cm_bm, axis=1)
                self.cluster_null_max["all"][start:end] = np.maximum(self.cluster_null_max["all"][start:end], row_max)
                if layer in self.cluster_null_max:
                    self.cluster_null_max[layer][start:end] = np.maximum(self.cluster_null_max[layer][start:end], row_max)

        p_flat = np.full(m_count * span, np.nan)
        p_flat[good_flat] = (1.0 + exceed_flat[good_flat]) / (self.B + 1.0)
        sums_flat = np.nansum(flat_diff, axis=0)
        means_flat = np.divide(sums_flat, n_ev_flat, out=np.full(m_count * span, np.nan), where=n_ev_flat > 0)

        z_2d = z_flat.reshape(m_count, span)
        p_2d = p_flat.reshape(m_count, span)
        nev_2d = n_ev_flat.reshape(m_count, span)
        nbg_2d = n_bg_flat.reshape(m_count, span)
        means_2d = means_flat.reshape(m_count, span)

        clusters: list[dict] = []
        row_frames: list[pd.DataFrame] = []
        pos_vec = np.arange(span) if width > 0 else -np.ones(span, dtype=int)

        for m_idx, (mid, rbp) in enumerate(zip(model_ids, rbp_labels)):
            z_m = z_2d[m_idx]
            good_m = np.isfinite(z_m)
            meta_m = {**base_meta, "model_id": mid}
            if rbp != "" or "RBP" in base_meta:
                meta_m["RBP"] = rbp

            if width == 0 or m_total == 1:
                df_m = pd.DataFrame({
                    "window_start0": pos_vec,
                    "n_matched_events": nev_2d[m_idx],
                    "n_independent_gene_blocks": nbg_2d[m_idx],
                    "paired_rank_Z": z_m,
                    "mean_paired_score_difference": means_2d[m_idx],
                    "raw_P": p_2d[m_idx],
                })
                for k, v in meta_m.items():
                    df_m[k] = v
                row_frames.append(df_m)

            if width == 0 and good_m.any():
                self._record_model_whole(layer, mid, float(np.nanmax(z_m)), track_null[:, m_idx])
            elif width > 0 and good_m.any():
                raw_runs = contiguous_runs(z_m)
                max_scaled = 0.0
                for a, b_end, raw_mass in raw_runs:
                    scaled_mass = raw_mass / scale_div
                    max_scaled = max(max_scaled, scaled_mass)
                    clusters.append({
                        **meta_m,
                        "first_window_start0": a,
                        "last_window_start0": b_end - 1,
                        "cluster_mass_raw": raw_mass,
                        "cluster_mass": scaled_mass,
                        "minimum_events": int(nev_2d[m_idx, a:b_end].min()),
                        "minimum_gene_blocks": int(nbg_2d[m_idx, a:b_end].min()),
                    })
                self._record_model_cluster(layer, mid, max_scaled, track_null[:, m_idx])

        out_rows = pd.concat(row_frames, ignore_index=True) if row_frames else pd.DataFrame()
        return out_rows, clusters

    def adjust_whole_regions(self, whole_rows: pd.DataFrame) -> pd.DataFrame:
        df = whole_rows.copy()
        if df.empty:
            df["whole_P_FWER_all"] = pd.Series(dtype=float)
            df["whole_P_FWER_layer"] = pd.Series(dtype=float)
            df["whole_q_BH"] = pd.Series(dtype=float)
            return df
        null_all = self.whole_null_max["all"]
        df["whole_P_FWER_all"] = [
            (1.0 + float(np.sum(null_all >= val - 1e-12))) / (self.B + 1.0) if np.isfinite(val) else np.nan
            for val in df["paired_rank_Z"]
        ]
        layer_p = []
        for val, lyr in zip(df["paired_rank_Z"], df["layer"]):
            if not np.isfinite(val) or lyr not in self.whole_null_max:
                layer_p.append(np.nan)
            else:
                nl = self.whole_null_max[lyr]
                layer_p.append((1.0 + float(np.sum(nl >= val - 1e-12))) / (self.B + 1.0))
        df["whole_P_FWER_layer"] = layer_p
        df["normal_upper_P"] = norm.sf(df["paired_rank_Z"].to_numpy(dtype=float))
        df["whole_q_BH"] = bh(df["normal_upper_P"].to_numpy(dtype=float))
        return df

    def adjust_clusters(self, clusters: list[dict]) -> pd.DataFrame:
        c = pd.DataFrame(clusters)
        if c.empty:
            return c
        null_all = self.cluster_null_max["all"]
        c["cluster_P_FWER"] = [
            (1.0 + float(np.sum(null_all >= val - 1e-12))) / (self.B + 1.0) for val in c["cluster_mass"]
        ]
        c["cluster_P_FWER_layer"] = [
            (1.0 + float(np.sum(self.cluster_null_max[lyr] >= val - 1e-12))) / (self.B + 1.0)
            for val, lyr in zip(c["cluster_mass"], c["layer"])
        ]
        return c

    def summarize_rbp_pivotal(self, rbp_to_models: dict[str, tuple[str, ...]], layer: str = "all") -> pd.DataFrame:
        """Exchangeable pooled-survival multi-PWM pivotal RBP ranking and exact Westfall-Young FWER."""
        all_mids = sorted(
            {
                m
                for mids in rbp_to_models.values()
                for m in mids
                if (layer, m) in self._model_whole_max or (layer, m) in self._model_cluster_max
            }
        )
        if not all_mids:
            return pd.DataFrame(columns=["RBP", "layer", "evaluable_models", "rbp_pivotal_Z", "rbp_P_FWER", "rank"])

        w_mat = np.vstack([
            self._model_whole_max.get((layer, m), np.full(self.B + 1, -np.inf)) for m in all_mids
        ])
        c_mat = np.vstack([
            self._model_cluster_max.get((layer, m), np.zeros(self.B + 1, dtype=float)) for m in all_mids
        ])
        s_w = _whole_z_log_survival(w_mat)
        s_c = _cluster_pareto_log_survival(c_mat)
        s_comb = np.maximum(s_w, s_c)
        mid_idx = {m: i for i, m in enumerate(all_mids)}

        rbp_names: list[str] = []
        rbp_piv_vecs: list[np.ndarray] = []
        n_models_used: list[int] = []

        for rbp, mids in sorted(rbp_to_models.items()):
            idxs = [mid_idx[m] for m in mids if m in mid_idx]
            if not idxs:
                continue
            u_vec = np.max(s_comb[idxs], axis=0)  # shape (B + 1,)
            # Symmetric standardization across ALL (B + 1) draws preserves exact {0..B} exchangeability!
            mu = float(np.mean(u_vec))
            sd = float(np.std(u_vec, ddof=1))
            if not np.isfinite(sd) or sd < 1e-9:
                piv = np.zeros(self.B + 1, dtype=float)
            else:
                piv = (u_vec - mu) / sd
            rbp_names.append(rbp)
            rbp_piv_vecs.append(piv)
            n_models_used.append(len(idxs))

        if not rbp_names:
            return pd.DataFrame(columns=["RBP", "layer", "evaluable_models", "rbp_pivotal_Z", "rbp_P_FWER", "rank"])

        mat = np.vstack(rbp_piv_vecs)  # shape (n_rbps, B + 1)
        obs = mat[:, 0]
        # Exact (B + 1)-exchangeable global max across RBPs (index 0 included via +1 in numerator!)
        global_null_max = np.max(mat[:, 1:], axis=0)
        p_fwer = np.array([(1.0 + float(np.sum(global_null_max >= v - 1e-12))) / (self.B + 1.0) for v in obs])

        out = pd.DataFrame({
            "RBP": rbp_names,
            "layer": layer,
            "evaluable_models": n_models_used,
            "rbp_pivotal_Z": obs,
            "rbp_P_FWER": p_fwer,
        })
        out["rank"] = out["rbp_pivotal_Z"].rank(ascending=False, method="average")
        out["rank_percentile"] = 100.0 * (len(out) - out["rank"]) / max(1, len(out) - 1)
        return out.sort_values("rank", kind="stable").reset_index(drop=True)


def fast_pwm_positions(sequence: str, model: dict) -> np.ndarray:
    if model["kind"] != "PWM":
        raise ValueError("Explicit PWM required")
    seq = sequence.upper().replace("U", "T")
    w = len(model["matrix"])
    n = max(0, len(seq) - w + 1)
    if not n:
        return np.empty(0, dtype=float)
    table = np.full(256, 4, dtype=np.int16)
    table[np.frombuffer(b"ACGT", dtype=np.uint8)] = np.arange(4)
    encoded = table[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]
    invalid = np.r_[0, np.cumsum(encoded == 4)]
    valid = (invalid[w:] - invalid[:-w]) == 0
    starts = np.flatnonzero(valid)
    odds = (20.0 * np.asarray(model["matrix"], dtype=float) + 0.025) / 20.1 / 0.25
    values = np.ones(len(starts), dtype=float)
    for k in range(w):
        values *= odds[k, encoded[starts + k]]
    out = np.full(n, np.nan, dtype=float)
    out[starts] = values
    return out


def pwm_profiles(sequence: str, model: dict, widths: tuple[int, ...] = (0, 10, 20, 50)) -> dict[int, np.ndarray]:
    pos = fast_pwm_positions(sequence, model)
    w = len(model["matrix"])
    valid = np.isfinite(pos)
    out: dict[int, np.ndarray] = {}
    for width in widths:
        if width == 0:
            out[0] = np.array([np.log2(1.0 + float(np.nanmean(pos))) if valid.any() else np.nan])
            continue
        span = max(0, len(sequence) - width + 1)
        values = np.full(span, np.nan, dtype=float)
        k = width - w + 1
        if k > 0 and span:
            view = np.lib.stride_tricks.sliding_window_view(pos, k)
            mask = np.isfinite(view).any(axis=1)
            values[mask] = np.log2(1.0 + np.nanmean(view[mask], axis=1))
        out[width] = values
    return out
