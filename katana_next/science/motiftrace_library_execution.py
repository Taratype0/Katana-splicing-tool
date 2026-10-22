"""Bounded all-library execution of the existing paired MotifTrace statistic.

Only arithmetic/storage are optimized. A finite motif retains count/coverage;
a native numeric PWM retains the declared continuous log2(1+mean odds) score.
"""
from __future__ import annotations

import warnings
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

from katana_next.science.motiftrace_clusters import CLUSTER_Z, cluster_runs
from katana_next.science.motiftrace_full_execution import profiles, signed_ranks


class PackedPWM:
    def __init__(self, words):
        self.offsets, pieces, offset = {}, [], 0
        for word in sorted(set(words)):
            self.offsets[word] = (offset, len(word))
            pieces.append(word+"N"*25)
            offset += len(word)+25
        text = "".join(pieces).encode("ascii")
        code = np.full(256, 4, np.uint8)
        code[np.frombuffer(b"ACGT", np.uint8)] = np.arange(4)
        self.encoded = code[np.frombuffer(text, np.uint8)]
        self.bad_prefix = np.r_[np.int32(0), np.cumsum(self.encoded == 4, dtype=np.int32)]

    def activate(self, matrix):
        matrix = np.asarray(matrix, dtype=float)
        if matrix.ndim != 2 or matrix.shape[1] != 4 or not np.isfinite(matrix).all() or (matrix < 0).any() or not np.allclose(matrix.sum(axis=1), 1):
            raise ValueError("Original normalized finite PWM required")
        self.width = len(matrix)
        smoothed = (20*matrix+.1*np.full(4, .25))/20.1
        odds = np.c_[smoothed/.25, np.ones(self.width)]
        n = max(0, len(self.encoded)-self.width+1)
        values = np.ones(n)
        for j in range(self.width):
            values *= odds[j, self.encoded[j:j+n]]
        valid = self.bad_prefix[self.width:self.width+n] == self.bad_prefix[:n]
        values[~valid] = np.nan
        self.values = values

    def positions(self, word):
        offset, length = self.offsets[word]
        return self.values[offset:offset+max(0, length-self.width+1)]

    def profile(self, word):
        values, result = self.positions(word), {}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            result[0] = (np.array([np.log2(1+np.nanmean(values))]) if len(values) else np.array([np.nan]),)
            for window in (10, 20, 50):
                span = max(0, len(word)-window+1)
                if not span or window < self.width or not len(values):
                    result[window] = (np.full(span, np.nan),)
                else:
                    # Direct row reductions match the original window kernel;
                    # prefix subtraction would lose small odds after large peaks.
                    view = sliding_window_view(values, window-self.width+1)
                    result[window] = (np.log2(1+np.nanmean(view, axis=1)),)
        return result


class DatasetPlan:
    def __init__(self, rna, pairs):
        regions = rna[rna.layer.isin(["affected_sequence", "regulatory_pre_mRNA"])]
        lookup = {(r.event_id, r.region): r.sequence for r in regions.itertuples()}
        self.specs = []
        for (kind, direction), group in pairs.groupby(["event_type", "direction"], sort=True):
            for region in sorted(set(regions.loc[regions.event_type.eq(kind), "region"])):
                words = [(lookup.get((r.foreground_event, region), ""), lookup.get((r.background_event, region), "")) for r in group.itertuples()]
                indices = group.global_pair_index.to_numpy()
                for window in (0, 10, 20, 50):
                    span = max(1, max([max(0, min(len(a), len(b))-window+1) for a, b in words], default=0)) if window else 1
                    coord = "masked_acceptor_3prime_minus20" if region.endswith("acceptor_intron") else "masked_donor_5prime_plus6" if region.endswith("donor_intron") else "affected_RNA_5prime"
                    self.specs.append({"meta": {"event_type": kind, "direction": direction, "region": region,
                        "window_nt": window, "coordinate_system": coord}, "words": words, "span": span, "pair_indices": indices})

    def datasets(self, model, packed):
        methods = ("finite_nonoverlap_counts", "overlapping_union_coverage") if model["kind"] == "consensus" else ("PWM_mean_odds",)
        cache = {}
        for spec in self.specs:
            for word in {w for pair in spec["words"] for w in pair if w}:
                if word not in cache:
                    cache[word] = profiles(word, model["pattern"]) if model["kind"] == "consensus" else packed.profile(word)
            window = spec["meta"]["window_nt"]
            for m, method in enumerate(methods):
                diff = np.full((len(spec["words"]), spec["span"]), np.nan)
                for i, (a, b) in enumerate(spec["words"]):
                    if not a or not b:
                        continue
                    x, y = cache[a][window][m], cache[b][window][m]
                    if window and spec["meta"]["region"].endswith("acceptor_intron"):
                        x, y = x[::-1], y[::-1]
                    n = min(len(x), len(y), spec["span"])
                    diff[i, :n] = x[:n]-y[:n]
                yield {"meta": {**spec["meta"], "method": method}, "differences": diff, "pair_indices": spec["pair_indices"]}


def compressed_null_max(values, lengths, threshold=CLUSTER_Z):
    positive = np.isfinite(values) & (values > threshold)
    sums = np.cumsum(np.where(positive, (values-threshold)*lengths[None, :], 0), axis=1)
    reset = np.maximum.accumulate(np.where(~positive, sums, 0), axis=1)
    return np.max(sums-reset, axis=1, initial=0)


def evaluate(diff, signs, is_scalar=False):
    """Keep original ranks and swaps; do not spend BLAS work on zero ranks."""
    n = np.isfinite(diff).sum(axis=0)
    observed = np.full(len(n), np.nan)
    p_index = np.zeros(len(n), np.int16)
    null_max = np.zeros(len(signs))
    if len(diff) < 20:
        return observed, n, p_index, null_max
    ranks, observed, n = signed_ranks(diff)
    eligible = np.isfinite(observed)
    p_index[eligible] = 1000
    if not eligible.any() or not np.any(ranks[:, eligible]):
        return observed, n, p_index, null_max
    cuts = np.r_[0, np.flatnonzero(np.any(ranks[:, 1:] != ranks[:, :-1], axis=0) | (eligible[1:] != eligible[:-1]))+1]
    lengths = np.diff(np.r_[cuts, ranks.shape[1]])
    selected = eligible[cuts] & np.any(ranks[:, cuts] != 0, axis=0)
    exceed = np.zeros(len(cuts), np.int32)
    for a in range(0, len(signs), 64):
        values = np.zeros((min(64, len(signs)-a), len(cuts)))
        values[:, ~eligible[cuts]] = np.nan
        values[:, selected] = signs[a:a+64] @ ranks[:, cuts[selected]]
        exceed += np.sum(values >= observed[cuts]-1e-12, axis=0)
        if not is_scalar:
            null_max[a:a+len(values)] = compressed_null_max(values, lengths)
    p_index[eligible] = np.repeat(exceed+1, lengths)[eligible]
    return observed, n, p_index, null_max


def histogram_bh(counts):
    """Exact BH for the frozen 1/1000 Monte Carlo grid, without a huge sort."""
    counts = np.asarray(counts, np.int64)
    total, cumulative = counts[1:].sum(), np.cumsum(counts[1:])
    values = np.ones(1000)
    occupied = counts[1:] > 0
    values[occupied] = np.arange(1, 1001)[occupied]/1000*total/cumulative[occupied]
    values = np.minimum.accumulate(values[::-1])[::-1].clip(0, 1)
    return np.r_[np.nan, values]
