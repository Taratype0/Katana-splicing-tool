from __future__ import annotations

import math
from collections.abc import Iterable

import numpy as np
from scipy.optimize import minimize
from scipy.special import betaln, expit, gammaln
from scipy.stats import chi2, hypergeom, spearmanr


def bh_adjust(values: Iterable[float]) -> np.ndarray:
    array = np.asarray(list(values), dtype=float)
    adjusted = np.full(array.shape, np.nan, dtype=float)
    valid = np.flatnonzero(np.isfinite(array))
    if not len(valid):
        return adjusted
    order = valid[np.argsort(array[valid])]
    ranked = array[order] * len(order) / np.arange(1, len(order) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted[order] = np.clip(ranked, 0.0, 1.0)
    return adjusted


def hypergeometric_overlap(
    universe_size: int, set_a_size: int, set_b_size: int, overlap_size: int
) -> float:
    if universe_size <= 0 or min(set_a_size, set_b_size, overlap_size) < 0:
        return math.nan
    return float(
        hypergeom.sf(overlap_size - 1, universe_size, set_a_size, set_b_size)
    )


def safe_spearman(left: Iterable[float], right: Iterable[float]) -> tuple[float, float]:
    a = np.asarray(list(left), dtype=float)
    b = np.asarray(list(right), dtype=float)
    valid = np.isfinite(a) & np.isfinite(b)
    if valid.sum() < 3 or np.nanstd(a[valid]) == 0 or np.nanstd(b[valid]) == 0:
        return math.nan, math.nan
    result = spearmanr(a[valid], b[valid])
    return float(result.statistic), float(result.pvalue)


def beta_binomial_log_likelihood(
    coefficients: np.ndarray,
    log_concentration: float,
    design: np.ndarray,
    successes: np.ndarray,
    totals: np.ndarray,
) -> float:
    probability = np.clip(expit(design @ coefficients), 1e-7, 1 - 1e-7)
    concentration = float(np.exp(np.clip(log_concentration, -5.0, 14.0)))
    alpha = probability * concentration
    beta = (1.0 - probability) * concentration
    combinations = (
        gammaln(totals + 1)
        - gammaln(successes + 1)
        - gammaln(totals - successes + 1)
    )
    values = combinations + betaln(successes + alpha, totals - successes + beta) - betaln(alpha, beta)
    return float(np.sum(values))


def beta_binomial_two_group(
    successes: Iterable[float], totals: Iterable[float], groups: Iterable[int]
) -> dict[str, float | str]:
    y = np.asarray(list(successes), dtype=float)
    n = np.asarray(list(totals), dtype=float)
    g = np.asarray(list(groups), dtype=float)
    valid = np.isfinite(y) & np.isfinite(n) & np.isfinite(g) & (n > 0) & (y >= 0) & (y <= n)
    y, n, g = y[valid], n[valid], g[valid]
    if len(y) < 4 or len(np.unique(g)) != 2 or min(np.sum(g == 0), np.sum(g == 1)) < 2:
        return {
            "status": "insufficient_replicates",
            "dpsi": math.nan,
            "log_odds": math.nan,
            "rho": math.nan,
            "pvalue": math.nan,
        }

    observed = (y + 0.5) / (n + 1.0)
    p0 = float(np.clip(np.average(observed, weights=n), 1e-4, 1 - 1e-4))
    base_logit = math.log(p0 / (1 - p0))
    null_design = np.ones((len(y), 1), dtype=float)
    full_design = np.column_stack([np.ones(len(y), dtype=float), g])

    def fit(design: np.ndarray, start: np.ndarray):
        return minimize(
            lambda parameters: -beta_binomial_log_likelihood(
                parameters[:-1], parameters[-1], design, y, n
            ),
            start,
            method="L-BFGS-B",
            bounds=[(-20.0, 20.0)] * (design.shape[1]) + [(-5.0, 14.0)],
            options={"maxiter": 500, "ftol": 1e-10},
        )

    null = fit(null_design, np.array([base_logit, math.log(20.0)]))
    group0 = float(np.clip(np.average(observed[g == 0], weights=n[g == 0]), 1e-4, 1 - 1e-4))
    group1 = float(np.clip(np.average(observed[g == 1], weights=n[g == 1]), 1e-4, 1 - 1e-4))
    start_effect = math.log(group1 / (1 - group1)) - math.log(group0 / (1 - group0))
    full = fit(full_design, np.array([math.log(group0 / (1 - group0)), start_effect, math.log(20.0)]))
    if not null.success or not full.success:
        return {
            "status": "optimizer_failed",
            "dpsi": group1 - group0,
            "log_odds": float(full.x[1]) if len(full.x) > 1 else math.nan,
            "rho": math.nan,
            "pvalue": math.nan,
        }
    likelihood_ratio = max(0.0, 2.0 * (-full.fun + null.fun))
    concentration = float(np.exp(full.x[-1]))
    return {
        "status": "fit",
        "dpsi": float(expit(full.x[0] + full.x[1]) - expit(full.x[0])),
        "log_odds": float(full.x[1]),
        "rho": float(1.0 / (concentration + 1.0)),
        "pvalue": float(chi2.sf(likelihood_ratio, 1)),
    }
