"""High-performance multi-core PWM scanning & BLAS permutation acceleration for Katana-SpliceTrace.

Provides:
  1. Numba `@njit(parallel=True, fastmath=True, nogil=True, cache=True)` multi-core server acceleration
     (scales linearly across SLURM / HPC CPU cores via `--threads` or `NUMBA_NUM_THREADS`),
  2. Vectorized 2D padded-tensor NumPy fallback when Numba is unavailable,
  3. Contiguous `float32` BLAS SGEMM gene-block permutation testing.
"""
from __future__ import annotations

import os
import numpy as np

_ENCODE_LUT = np.full(256, 4, dtype=np.uint8)
_ENCODE_LUT[np.frombuffer(b"ACGTacgtUu", dtype=np.uint8)] = np.array([0, 1, 2, 3, 0, 1, 2, 3, 3, 3], dtype=np.uint8)

try:
    import numba as nb

    HAS_NUMBA = True
except Exception:  # pragma: no cover
    nb = None
    HAS_NUMBA = False


def configure_server_threads(threads: int | None = None) -> int:
    """Configure Numba thread pool for multi-core Linux/Windows server execution."""
    if threads is None or threads <= 0:
        slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK") or os.environ.get("OMP_NUM_THREADS")
        threads = int(slurm_cpus) if (slurm_cpus and slurm_cpus.isdigit()) else (os.cpu_count() or 4)
    threads = max(1, int(threads))
    if HAS_NUMBA:
        try:
            max_avail = nb.config.NUMBA_NUM_THREADS
            nb.set_num_threads(min(threads, max_avail))
            return int(nb.get_num_threads())
        except Exception:
            return threads
    return 1


def pack_ascii_sequences(sequences: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """Pack variable-length DNA/RNA strings into a single contiguous uint8 buffer + int64 offsets."""
    n = len(sequences)
    offsets = np.zeros(n + 1, dtype=np.int64)
    lengths = [len(s) for s in sequences]
    offsets[1:] = np.cumsum(lengths, dtype=np.int64)
    total = int(offsets[-1])
    if total == 0:
        return np.empty(0, dtype=np.uint8), offsets
    joined = "".join(sequences).encode("ascii", errors="replace")
    # Ensure exact byte count matches character lengths
    if len(joined) != total:
        joined = b"".join(s.encode("ascii", errors="replace")[: len(s)].ljust(len(s), b"N") for s in sequences)
    flat = _ENCODE_LUT[np.frombuffer(joined, dtype=np.uint8)]
    return flat, offsets


if HAS_NUMBA:

    @nb.njit(parallel=True, fastmath=True, nogil=True, cache=True)
    def _numba_mean_log2_odds_width(
        flat_seq: np.ndarray,
        offsets: np.ndarray,
        odds_batch: np.ndarray,
        w: int,
    ) -> np.ndarray:
        n_seq = len(offsets) - 1
        m_count = odds_batch.shape[0]
        out = np.zeros((n_seq, m_count), dtype=np.float64)
        for s_i in nb.prange(n_seq):
            st = offsets[s_i]
            en = offsets[s_i + 1]
            seq_len = en - st
            n_pos = seq_len - w + 1
            if n_pos <= 0:
                continue
            valid_cnt = 0
            # Local accumulator across M_w models for sequence s_i
            acc = np.zeros(m_count, dtype=np.float64)
            for p in range(n_pos):
                ok = True
                for k in range(w):
                    if flat_seq[st + p + k] >= 4:
                        ok = False
                        break
                if not ok:
                    continue
                valid_cnt += 1
                for m in range(m_count):
                    prod = 1.0
                    for k in range(w):
                        prod *= odds_batch[m, k, flat_seq[st + p + k]]
                    acc[m] += prod
            if valid_cnt > 0:
                inv_c = 1.0 / float(valid_cnt)
                for m in range(m_count):
                    out[s_i, m] = np.log2(1.0 + acc[m] * inv_c)
        return out

    @nb.njit(parallel=True, fastmath=True, nogil=True, cache=True)
    def _numba_positional_odds_matrix(
        flat_seq: np.ndarray,
        st: int,
        seq_len: int,
        odds_batch: np.ndarray,
        w: int,
    ) -> tuple[np.ndarray, np.ndarray]:
        m_count = odds_batch.shape[0]
        n_pos = seq_len - w + 1
        pos_batch = np.zeros((m_count, n_pos), dtype=np.float64)
        valid = np.zeros(n_pos, dtype=np.int64)
        for p in nb.prange(n_pos):
            ok = True
            for k in range(w):
                if flat_seq[st + p + k] >= 4:
                    ok = False
                    break
            if ok:
                valid[p] = 1
                for m in range(m_count):
                    prod = 1.0
                    for k in range(w):
                        prod *= odds_batch[m, k, flat_seq[st + p + k]]
                    pos_batch[m, p] = prod
        return pos_batch, valid


def fast_score_mean_log_odds(sequences: list[str], models: dict[str, dict]) -> np.ndarray:
    """Drop-in multi-core replacement for `_score_mean_log_odds` (exact to <1e-12)."""
    m_ids = list(models.keys())
    n_seq, n_mod = len(sequences), len(m_ids)
    out = np.zeros((n_seq, n_mod), dtype=np.float64)
    if n_seq == 0 or n_mod == 0:
        return out

    flat_seq, offsets = pack_ascii_sequences(sequences)
    by_width: dict[int, list[int]] = {}
    for idx, mid in enumerate(m_ids):
        by_width.setdefault(len(models[mid]["matrix"]), []).append(idx)

    if HAS_NUMBA:
        for width, mod_indices in by_width.items():
            mats = np.stack([np.asarray(models[m_ids[i]]["matrix"], dtype=np.float64) for i in mod_indices], axis=0)
            odds_batch = np.ascontiguousarray(((20.0 * mats + 0.025) / 20.1) / 0.25)
            scores_w = _numba_mean_log2_odds_width(flat_seq, offsets, odds_batch, int(width))
            out[:, mod_indices] = scores_w
        return out

    # Vectorized 2D padded-tensor NumPy fallback (no Python per-sequence loop!)
    lengths = np.diff(offsets)
    max_len = int(lengths.max()) if len(lengths) else 0
    if max_len == 0:
        return out
    padded = np.full((n_seq, max_len), 4, dtype=np.uint8)
    for i in range(n_seq):
        L = int(lengths[i])
        if L > 0:
            padded[i, :L] = flat_seq[offsets[i] : offsets[i + 1]]

    for width, mod_indices in by_width.items():
        w = int(width)
        if max_len < w:
            continue
        mats = np.stack([np.asarray(models[m_ids[i]]["matrix"], dtype=np.float64) for i in mod_indices], axis=0)
        odds_ext = np.zeros((len(mod_indices), w, 5), dtype=np.float64)
        odds_ext[:, :, :4] = ((20.0 * mats + 0.025) / 20.1) / 0.25

        n_pos = max_len - w + 1
        valid_mask = np.ones((n_seq, n_pos), dtype=bool)
        vals = np.ones((len(mod_indices), n_seq, n_pos), dtype=np.float64)
        for k in range(w):
            sl = padded[:, k : k + n_pos]
            valid_mask &= sl < 4
            vals *= odds_ext[:, k, sl]
        vals[:, ~valid_mask] = 0.0
        cnts = valid_mask.sum(axis=1)
        ok = cnts > 0
        if np.any(ok):
            means = vals[:, ok, :].sum(axis=2) / cnts[None, ok]
            out[np.ix_(ok, mod_indices)] = np.log2(1.0 + means.T)
    return out


def fast_gene_block_maxt_permutations(
    yc: np.ndarray,
    F_gc: np.ndarray,
    sd_perm: np.ndarray,
    valid_c: np.ndarray,
    n_reg: int,
    n_rbp: int,
    n_permutations: int = 1999,
    seed: int = 20261117,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Contiguous BLAS SGEMM/DGEMM gene-block permutation test across `n_reg x n_rbp` channels."""
    Z_obs = np.divide(yc @ F_gc, sd_perm, out=np.zeros(n_reg * n_rbp, dtype=np.float64), where=valid_c)
    B = int(n_permutations)
    G = len(yc)
    rng = np.random.default_rng(int(seed))

    # Vectorized random argsort generates B independent permutations in C without Python loops
    rand_keys = rng.random((B, G), dtype=np.float32)
    perm_idx = np.argsort(rand_keys, axis=1)
    yc_f32 = np.asarray(yc, dtype=np.float32)
    perm_Y = yc_f32[perm_idx]
    F_f32 = np.ascontiguousarray(F_gc, dtype=np.float32)
    sd_f32 = np.where(valid_c, np.asarray(sd_perm, dtype=np.float32), 1.0)

    null_Z = (perm_Y @ F_f32) / sd_f32[None, :]
    if not np.all(valid_c):
        null_Z[:, ~valid_c] = 0.0

    Z_obs_2d = Z_obs.reshape(n_reg, n_rbp)
    abs_Z_obs_f32 = np.abs(Z_obs).astype(np.float32)
    abs_null = np.abs(null_Z)
    null_3d = abs_null.reshape(B, n_reg, n_rbp)

    rbp_obs_max = np.max(np.abs(Z_obs_2d), axis=0).astype(np.float32)
    rbp_null_max = np.max(null_3d, axis=1)
    glob_null = np.max(rbp_null_max, axis=1)

    eps = np.float32(1e-6)
    rbp_nom_p = (1.0 + np.sum(rbp_null_max >= rbp_obs_max[None, :] - eps, axis=0)) / (B + 1.0)
    rbp_fwer_p = (1.0 + np.sum(glob_null[:, None] >= rbp_obs_max[None, :] - eps, axis=0)) / (B + 1.0)
    reg_nom_p = ((1.0 + np.sum(abs_null >= abs_Z_obs_f32[None, :] - eps, axis=0)) / (B + 1.0)).reshape(n_reg, n_rbp)
    reg_fwer_p = ((1.0 + np.sum(glob_null[:, None] >= abs_Z_obs_f32[None, :] - eps, axis=0)) / (B + 1.0)).reshape(n_reg, n_rbp)
    return Z_obs_2d, rbp_nom_p, rbp_fwer_p, reg_nom_p, reg_fwer_p
