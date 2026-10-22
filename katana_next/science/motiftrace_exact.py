"""Conditional signed-rank tails with ties, without a Monte Carlo P floor.

Each nonzero paired difference contributes its integer doubled midrank with
probability 1/2 under the same conditional pair-swap null used previously.
Polynomial coefficients enumerate that distribution, including tied ranks.
No asymptotic-normal or saddlepoint fallback is used.
"""
from __future__ import annotations

from collections import OrderedDict
from math import gcd
from functools import reduce

import numpy as np
from scipy.stats import rankdata


class ExactSignedRank:
    def __init__(self, maximum_cache_bytes=128 * 1024**2):
        self.cache = OrderedDict()
        self.cache_bytes = 0
        self.maximum_cache_bytes = maximum_cache_bytes
        self.hits = self.misses = 0

    def tail(self, weights):
        weights = np.sort(np.asarray(weights, dtype=np.int64))
        divisor = reduce(gcd, map(int, weights), 0) or 1
        weights = weights // divisor
        key = weights.astype('<i2').tobytes()
        if key in self.cache:
            self.hits += 1
            self.cache.move_to_end(key)
            return self.cache[key], divisor
        self.misses += 1
        # All coefficients are nonnegative: adding tiny extreme-tail masses
        # never involves subtracting two nearly equal CDF values.
        pmf = np.ones(1, dtype=np.float64)
        for weight in weights:
            weight = int(weight)
            updated = np.zeros(len(pmf) + weight, dtype=np.float64)
            updated[:len(pmf)] += .5 * pmf
            updated[weight:] += .5 * pmf
            pmf = updated
        tail = np.cumsum(pmf[::-1])[::-1]
        tail = np.minimum(tail, 1.)
        while self.cache and self.cache_bytes + tail.nbytes > self.maximum_cache_bytes:
            _, previous = self.cache.popitem(last=False)
            self.cache_bytes -= previous.nbytes
        if tail.nbytes <= self.maximum_cache_bytes:
            self.cache[key] = tail
            self.cache_bytes += tail.nbytes
        return tail, divisor

    def one(self, differences):
        d = np.asarray(differences, float)
        d = d[np.isfinite(d) & (d != 0)]
        if not len(d):
            return 1., 0, 1.
        weights = np.rint(2 * rankdata(np.abs(d), method='average')).astype(np.int64)
        tail, divisor = self.tail(weights)
        observed = int(weights[d > 0].sum()) // divisor
        return float(tail[observed]), len(d), float(np.ldexp(1., -len(d)))

    def evaluate(self, differences, minimum=20):
        d = np.asarray(differences, float)
        valid = np.isfinite(d)
        n = valid.sum(axis=0)
        absolute = np.abs(d)
        absolute[(absolute == 0) | ~valid] = np.nan
        ranks = rankdata(absolute, axis=0, method='average', nan_policy='omit')
        weights = np.rint(2 * np.nan_to_num(ranks)).astype(np.int16)
        signed = weights.astype(float) * np.sign(np.nan_to_num(d))
        denominator = np.sqrt(np.sum(weights.astype(float)**2, axis=0))
        z = np.divide(signed.sum(axis=0), denominator, out=np.zeros(len(n)), where=denominator > 0)
        z[n < minimum] = np.nan
        p = np.full(len(n), np.nan)
        nonzero = (weights > 0).sum(axis=0)
        minimum_p = np.ldexp(np.ones(len(n)), -nonzero)
        # Adjacent identical rank/sign columns share the exact same statistic.
        cuts = np.r_[0, np.flatnonzero(np.any(signed[:, 1:] != signed[:, :-1], axis=0)
                                      | (n[1:] != n[:-1])) + 1]
        lengths = np.diff(np.r_[cuts, len(n)])
        for start, length in zip(cuts, lengths):
            if n[start] < minimum:
                continue
            w = weights[:, start]
            w = w[w > 0]
            if len(w):
                tail, divisor = self.tail(w)
                observed = int(weights[d[:, start] > 0, start].sum()) // divisor
                value = float(tail[observed])
            else:
                value = 1.
            p[start:start + length] = value
        return z, n, p, nonzero, minimum_p


def full_family_adjust(p_values):
    """BH plus arbitrary-dependence BY and Holm; no family subsetting."""
    p = np.asarray(p_values, float)
    if not np.isfinite(p).all() or ((p < 0) | (p > 1)).any():
        raise ValueError('Finite full-family P values required')
    m = len(p)
    order = np.argsort(p, kind='stable')
    ranks = np.arange(1, m + 1)
    ordered = p[order]
    bh = np.minimum.accumulate((ordered * (m / ranks))[::-1])[::-1].clip(0, 1)
    holm = np.maximum.accumulate(ordered * (m - ranks + 1)).clip(0, 1)
    harmonic = float(np.sum(1. / ranks))
    outputs = []
    for values in (bh, np.minimum(1., bh * harmonic), holm):
        out = np.empty(m)
        out[order] = values
        outputs.append(out)
    return outputs
