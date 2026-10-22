from __future__ import annotations

import math
from collections.abc import Iterable

import numpy as np
from scipy.optimize import minimize
from scipy.special import betaln, expit, gammaln, digamma, softmax
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
    return float(hypergeom.sf(overlap_size - 1, universe_size, set_a_size, set_b_size))


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
        gammaln(totals + 1) - gammaln(successes + 1) - gammaln(totals - successes + 1)
    )
    values = (
        combinations
        + betaln(successes + alpha, totals - successes + beta)
        - betaln(alpha, beta)
    )
    return float(np.sum(values))


def beta_binomial_two_group(
    successes: Iterable[float], totals: Iterable[float], groups: Iterable[int]
) -> dict[str, float | str]:
    return beta_binomial_group_test(successes, totals, groups)


def beta_binomial_group_test(
    successes: Iterable[float],
    totals: Iterable[float],
    groups: Iterable[int],
    covariates: np.ndarray | None = None,
    *, analytic_gradient: bool = False,
) -> dict[str, float | str]:
    """Test one binary condition effect with optional numeric nuisance terms."""

    y = np.asarray(list(successes), dtype=float)
    n = np.asarray(list(totals), dtype=float)
    g = np.asarray(list(groups), dtype=float)
    nuisance = _numeric_covariates(covariates, len(y))
    valid = (
        np.isfinite(y)
        & np.isfinite(n)
        & np.isfinite(g)
        & np.isfinite(nuisance).all(axis=1)
        & (n > 0)
        & (y >= 0)
        & (y <= n)
    )
    y, n, g, nuisance = y[valid], n[valid], g[valid], nuisance[valid]
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
    null_design = np.column_stack([np.ones(len(y), dtype=float), nuisance])
    full_design = np.column_stack([null_design, g])
    if (
        np.linalg.matrix_rank(null_design) < null_design.shape[1]
        or np.linalg.matrix_rank(full_design) < full_design.shape[1]
    ):
        return {
            "status": "confounded_covariates",
            "dpsi": math.nan,
            "log_odds": math.nan,
            "rho": math.nan,
            "pvalue": math.nan,
        }

    def fit(design: np.ndarray, start: np.ndarray):
        return minimize(
            (lambda parameters: _bb_value_gradient(parameters, design, y, n)) if analytic_gradient
            else (lambda parameters: -beta_binomial_log_likelihood(parameters[:-1], parameters[-1], design, y, n)),
            start,
            method="L-BFGS-B",
            jac=analytic_gradient,
            bounds=[(-20.0, 20.0)] * (design.shape[1]) + [(-5.0, 14.0)],
            options={"maxiter": 500, "ftol": 1e-10},
        )

    null = fit(
        null_design,
        np.concatenate([[base_logit], np.zeros(nuisance.shape[1]), [math.log(20.0)]]),
    )
    group0 = float(
        np.clip(np.average(observed[g == 0], weights=n[g == 0]), 1e-4, 1 - 1e-4)
    )
    group1 = float(
        np.clip(np.average(observed[g == 1], weights=n[g == 1]), 1e-4, 1 - 1e-4)
    )
    start_effect = math.log(group1 / (1 - group1)) - math.log(group0 / (1 - group0))
    full = fit(
        full_design,
        np.concatenate(
            [
                [math.log(group0 / (1 - group0))],
                np.zeros(nuisance.shape[1]),
                [start_effect, math.log(20.0)],
            ]
        ),
    )
    if analytic_gradient and null.success:
        nested = fit(full_design, np.r_[null.x[:-1], 0., null.x[-1]])
        successful = [candidate for candidate in (full, nested) if candidate.success and np.isfinite(candidate.fun)]
        if successful:
            full = min(successful, key=lambda candidate: candidate.fun)
    if not null.success or not full.success:
        return {
            "status": "optimizer_failed",
            "dpsi": group1 - group0,
            "log_odds": float(full.x[-2]) if len(full.x) > 1 else math.nan,
            "rho": math.nan,
            "pvalue": math.nan,
        }
    likelihood_ratio = max(0.0, 2.0 * (-full.fun + null.fun) * (len(y) if analytic_gradient else 1))
    concentration = float(np.exp(full.x[-1]))
    mean_nuisance = nuisance.mean(axis=0) if nuisance.shape[1] else np.empty(0)
    group0_design = np.concatenate([[1.0], mean_nuisance, [0.0]])
    group1_design = group0_design.copy()
    group1_design[-1] = 1.0
    coefficients = full.x[:-1]
    return {
        "status": "fit",
        "dpsi": float(
            expit(group1_design @ coefficients) - expit(group0_design @ coefficients)
        ),
        "log_odds": float(full.x[-2]),
        "rho": float(1.0 / (concentration + 1.0)),
        "pvalue": float(chi2.sf(likelihood_ratio, 1)),
        "null_log_likelihood": float(-null.fun * (len(y) if analytic_gradient else 1)),
        "full_log_likelihood": float(-full.fun * (len(y) if analytic_gradient else 1)),
        "optimizer": "analytic_gradient_nested_start" if analytic_gradient else "finite_difference",
    }


def _bb_value_gradient(parameters, design, successes, totals):
    """Mean negative BB likelihood and derivative of the existing clipped model."""
    raw = expit(design @ parameters[:-1])
    probability = np.clip(raw, 1e-7, 1-1e-7)
    concentration = np.exp(np.clip(parameters[-1], -5., 14.))
    a, b = probability*concentration, (1-probability)*concentration
    da = digamma(successes+a)-digamma(a)
    db = digamma(totals-successes+b)-digamma(b)
    derivative = raw*(1-raw)*((raw > 1e-7) & (raw < 1-1e-7))
    beta_score = design.T @ (concentration*(da-db)*derivative)
    concentration_score = concentration*(probability*da+(1-probability)*db+digamma(concentration)-digamma(totals+concentration))
    k_score = concentration_score.sum() if -5. <= parameters[-1] <= 14. else 0.
    value = -beta_binomial_log_likelihood(parameters[:-1], parameters[-1], design, successes, totals)/len(totals)
    return value, -np.r_[beta_score, k_score]/len(totals)


def _dm_value_gradient(parameters: np.ndarray, design: np.ndarray, matrix: np.ndarray):
    """Mean negative DM log likelihood and analytic gradient (last choice baseline)."""
    choices = matrix.shape[1]
    coefficients = parameters[:-1].reshape(design.shape[1], choices - 1)
    logits = np.column_stack([design @ coefficients, np.zeros(len(matrix))])
    probability = softmax(logits, axis=1)
    concentration = np.exp(parameters[-1])
    alpha = probability * concentration
    totals = matrix.sum(axis=1)
    loglik = (gammaln(totals + 1) - gammaln(matrix + 1).sum(axis=1)
              + gammaln(concentration) - gammaln(totals + concentration)
              + (gammaln(matrix + alpha) - gammaln(alpha)).sum(axis=1))
    alpha_score = digamma(matrix + alpha) - digamma(alpha)
    centered_score = alpha_score - (alpha_score * probability).sum(axis=1, keepdims=True)
    beta_score = design.T @ (alpha * centered_score)[:, :-1]
    concentration_score = concentration * (digamma(concentration) - digamma(totals + concentration)
                                           + (alpha_score * probability).sum(axis=1))
    gradient = -np.concatenate([beta_score.ravel(), [concentration_score.sum()]]) / len(matrix)
    return -float(loglik.mean()), gradient


def dirichlet_multinomial_two_group(
    counts: np.ndarray,
    groups: Iterable[int],
    covariates: np.ndarray | None = None,
) -> dict[str, object]:
    """Likelihood-ratio test for a multi-junction choice across two groups.

    Rows are independent biological units and columns are mutually exclusive
    junction choices. A single concentration parameter captures extra-
    multinomial variation, matching the parsimonious scQuint model family.
    """

    matrix = np.asarray(counts, dtype=float)
    group = np.asarray(list(groups), dtype=float)
    if matrix.ndim != 2 or matrix.shape[1] < 2 or len(group) != len(matrix):
        return _failed_dirichlet_multinomial(
            "invalid_shape", matrix.shape[1] if matrix.ndim == 2 else 0
        )
    if not np.isfinite(matrix).all() or (matrix < 0).any() or (matrix % 1).any():
        return _failed_dirichlet_multinomial("invalid_counts", matrix.shape[1])
    nuisance = _numeric_covariates(covariates, len(matrix))
    valid = (
        np.isfinite(matrix).all(axis=1)
        & (matrix >= 0).all(axis=1)
        & (matrix.sum(axis=1) > 0)
        & np.isfinite(group)
        & np.isfinite(nuisance).all(axis=1)
    )
    matrix, group, nuisance = matrix[valid], group[valid], nuisance[valid]
    choices = matrix.shape[1]
    if (
        len(matrix) < 4
        or len(np.unique(group)) != 2
        or min(np.sum(group == 0), np.sum(group == 1)) < 2
    ):
        return _failed_dirichlet_multinomial("insufficient_replicates", choices)

    totals = matrix.sum(axis=1)
    pooled = (matrix.sum(axis=0) + 0.5) / (totals.sum() + 0.5 * choices)
    reference = max(float(pooled[-1]), 1e-7)
    base_logits = np.log(np.clip(pooled[:-1], 1e-7, 1.0) / reference)

    def probabilities(parameters: np.ndarray, design: np.ndarray) -> np.ndarray:
        coefficients = parameters[:-1].reshape(design.shape[1], choices - 1)
        logits = design @ coefficients
        logits = np.column_stack([logits, np.zeros(len(logits))])
        logits -= logits.max(axis=1, keepdims=True)
        exponential = np.exp(logits)
        return exponential / exponential.sum(axis=1, keepdims=True)

    def fit(design: np.ndarray, start: np.ndarray):
        return minimize(
            _dm_value_gradient,
            start,
            args=(design, matrix),
            jac=True,
            method="L-BFGS-B",
            bounds=[(-20.0, 20.0)] * (design.shape[1] * (choices - 1)) + [(-5.0, 14.0)],
            options={"maxiter": 2000, "ftol": 1e-13, "gtol": 1e-7, "maxls": 50},
        )

    null_design = np.column_stack([np.ones(len(matrix), dtype=float), nuisance])
    full_design = np.column_stack([null_design, group])
    if (
        np.linalg.matrix_rank(null_design) < null_design.shape[1]
        or np.linalg.matrix_rank(full_design) < full_design.shape[1]
    ):
        return _failed_dirichlet_multinomial("confounded_covariates", choices)
    null_coefficients = np.zeros((null_design.shape[1], choices - 1), dtype=float)
    null_coefficients[0] = base_logits
    full_coefficients = np.zeros((full_design.shape[1], choices - 1), dtype=float)
    full_coefficients[0] = base_logits
    null_start = np.concatenate([null_coefficients.ravel(), [math.log(20.0)]])
    full_start = np.concatenate([full_coefficients.ravel(), [math.log(20.0)]])
    null_fits = [fit(null_design, np.concatenate([null_start[:-1], [math.log(c)]])) for c in (2., 20., 200.)]
    converged_null = [r for r in null_fits if r.success and np.isfinite(r.fun)]
    if not converged_null:
        return _failed_dirichlet_multinomial("optimizer_failed", choices)
    null = min(converged_null, key=lambda r: r.fun)
    nested = np.vstack([null.x[:-1].reshape(null_design.shape[1], choices - 1), np.zeros(choices - 1)])
    warm_start = np.concatenate([nested.ravel(), [null.x[-1]]])
    grouped_coefficients = full_coefficients.copy()
    for g in (0, 1):
        proportions = matrix[group == g].sum(axis=0) + .5
        logits = np.log(proportions[:-1] / proportions[-1])
        if g == 0:
            grouped_coefficients[0] = logits
        else:
            grouped_coefficients[-1] = logits - grouped_coefficients[0]
    group_start = np.concatenate([grouped_coefficients.ravel(), [null.x[-1]]])
    full_fits = [fit(full_design, start) for start in (warm_start, group_start, full_start)]
    converged_full = [r for r in full_fits if r.success and np.isfinite(r.fun)]
    if not converged_full:
        return _failed_dirichlet_multinomial("optimizer_failed", choices)
    full = min(converged_full, key=lambda r: r.fun)
    if full.fun > null.fun + 1e-7:
        return _failed_dirichlet_multinomial("non_nested_optimum", choices)

    mean_nuisance = nuisance.mean(axis=0) if nuisance.shape[1] else np.empty(0)
    group0_design = np.concatenate([[1.0], mean_nuisance, [0.0]])[None, :]
    group1_design = group0_design.copy()
    group1_design[0, -1] = 1.0
    group0 = probabilities(full.x, group0_design)[0]
    group1 = probabilities(full.x, group1_design)[0]
    likelihood_ratio = max(0.0, 2.0 * len(matrix) * (-full.fun + null.fun))
    degrees = choices - 1
    concentration = float(np.exp(full.x[-1]))
    return {
        "status": "fit",
        "delta_proportions": group1 - group0,
        "rho": float(1.0 / (concentration + 1.0)),
        "likelihood_ratio": likelihood_ratio,
        "degrees_of_freedom": degrees,
        "pvalue": float(chi2.sf(likelihood_ratio, degrees)),
        "null_log_likelihood": float(-null.fun * len(matrix)),
        "full_log_likelihood": float(-full.fun * len(matrix)),
        "optimizer": "analytic_gradient_multistart_LBFGSB",
    }


def _failed_dirichlet_multinomial(status: str, choices: int) -> dict[str, object]:
    return {
        "status": status,
        "delta_proportions": np.repeat(np.nan, max(0, choices)),
        "rho": math.nan,
        "likelihood_ratio": math.nan,
        "degrees_of_freedom": max(0, choices - 1),
        "pvalue": math.nan,
    }


def _numeric_covariates(covariates: np.ndarray | None, observations: int) -> np.ndarray:
    if covariates is None:
        return np.empty((observations, 0), dtype=float)
    matrix = np.asarray(covariates, dtype=float)
    if matrix.ndim == 1:
        matrix = matrix[:, None]
    if matrix.ndim != 2 or len(matrix) != observations:
        raise ValueError("Covariates must have one row per observation.")
    return matrix
