from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import norm


@dataclass
class BayesianPreservationResult:
    module_posteriors: pd.DataFrame
    metric_posteriors: pd.DataFrame


def _normal_hierarchical_grid(
    values: np.ndarray,
    standard_errors: np.ndarray,
    *,
    prior_sd: float,
    tau_scale: float,
    grid_size: int,
    partial_threshold: float,
    strict_threshold: float,
) -> dict[str, float]:
    tau = np.linspace(0.0, tau_scale * 4.0, grid_size)
    log_tau_prior = np.log(np.sqrt(2.0 / np.pi) / tau_scale) - (tau**2) / (
        2.0 * tau_scale**2
    )
    post_mean = np.zeros(grid_size)
    post_variance = np.zeros(grid_size)
    log_weight = np.zeros(grid_size)
    for index, tau_value in enumerate(tau):
        variance = standard_errors**2 + tau_value**2
        precision = 1.0 / prior_sd**2 + np.sum(1.0 / variance)
        post_variance[index] = 1.0 / precision
        post_mean[index] = post_variance[index] * np.sum(values / variance)
        marginal_variance = variance + prior_sd**2
        log_weight[index] = (
            log_tau_prior[index]
            + np.sum(norm.logpdf(values, loc=0.0, scale=np.sqrt(marginal_variance)))
        )
    weights = np.exp(log_weight - logsumexp(log_weight))
    mean = float(np.sum(weights * post_mean))
    second = float(np.sum(weights * (post_variance + post_mean**2)))
    variance = max(second - mean**2, 1e-12)
    draws_z = (post_mean[:, None] + np.sqrt(post_variance)[:, None] * np.linspace(-4, 4, 801))
    draw_weights = weights[:, None] * norm.pdf(np.linspace(-4, 4, 801))[None, :]
    draw_weights /= draw_weights.sum()
    flat = draws_z.ravel()
    flat_weight = draw_weights.ravel()
    order = np.argsort(flat)
    cumulative = np.cumsum(flat_weight[order])
    lower = float(flat[order][np.searchsorted(cumulative, 0.025)])
    upper = float(flat[order][np.searchsorted(cumulative, 0.975)])
    return {
        "posterior_mean": mean,
        "posterior_sd": float(np.sqrt(variance)),
        "credible_interval_low": lower,
        "credible_interval_high": upper,
        "posterior_probability_positive": float(
            np.sum(draw_weights[draws_z > 0.0])
        ),
        "posterior_probability_partial": float(
            np.sum(draw_weights[draws_z >= partial_threshold])
        ),
        "posterior_probability_strict": float(
            np.sum(draw_weights[draws_z >= strict_threshold])
        ),
        "posterior_tau_mean": float(np.sum(weights * tau)),
    }


def bayesian_multi_cohort_preservation(
    preservation: pd.DataFrame,
    *,
    module_columns: tuple[str, ...] = ("consensus_module",),
    prior_sd: float = 5.0,
    tau_scale: float = 2.0,
    grid_size: int = 301,
) -> BayesianPreservationResult:
    """Combine external preservation evidence with a hierarchical normal model.

    Density uses the permutation-calibrated WGCNA-like Z statistic. Structural
    similarity metrics use their direct correlations after Fisher
    transformation. A negative permutation Zconnectivity is not treated as a
    negative observed connectivity correlation.
    """
    connectivity_column = (
        "connectivity_correlation"
        if "connectivity_correlation" in preservation
        else "Zconnectivity"
    )
    metrics = {
        "density": ("Zdensity", "density_empirical_pvalue", "z_statistic"),
        "connectivity": (
            connectivity_column,
            (
                None
                if connectivity_column == "connectivity_correlation"
                else "connectivity_empirical_pvalue"
            ),
            (
                "correlation"
                if connectivity_column == "connectivity_correlation"
                else "z_statistic"
            ),
        ),
        "program": ("program_score_preservation", None, "correlation"),
        "membership": ("membership_profile_correlation", None, "correlation"),
    }
    missing_id = set(module_columns).difference(preservation.columns)
    if missing_id:
        raise ValueError(f"Preservation table is missing module identifiers: {sorted(missing_id)}")
    rows: list[dict[str, object]] = []
    for key, part in preservation.groupby(list(module_columns), observed=True, dropna=False):
        key_values = key if isinstance(key, tuple) else (key,)
        identity = dict(zip(module_columns, key_values, strict=True))
        for metric, (value_column, pvalue_column, value_scale) in metrics.items():
            if value_column not in part:
                continue
            values = pd.to_numeric(part[value_column], errors="coerce").to_numpy()
            valid = np.isfinite(values)
            if valid.sum() < 1:
                continue
            values = values[valid]
            if value_scale == "z_statistic" and pvalue_column and pvalue_column in part:
                pvalues = pd.to_numeric(part.loc[valid, pvalue_column], errors="coerce").clip(
                    lower=1e-12, upper=1.0
                )
                z_equivalent = np.maximum(norm.isf(pvalues.to_numpy() / 2.0), 0.25)
                standard_errors = np.maximum(np.abs(values) / z_equivalent, 0.25)
                partial_threshold = 2.0
                strict_threshold = 10.0
            else:
                standard_errors = np.full(len(values), 0.5, dtype=float)
                if value_scale == "correlation":
                    values = np.arctanh(np.clip(values, -0.999, 0.999))
                partial_threshold = 0.0
                strict_threshold = float(np.arctanh(0.50))
            posterior = _normal_hierarchical_grid(
                values,
                standard_errors,
                prior_sd=prior_sd,
                tau_scale=tau_scale,
                grid_size=grid_size,
                partial_threshold=partial_threshold,
                strict_threshold=strict_threshold,
            )
            rows.append(
                {
                    **identity,
                    "metric": metric,
                    "n_external_cohorts": len(values),
                    "observed_mean": float(np.mean(values)),
                    "input_column": value_column,
                    "modeled_scale": value_scale,
                    "partial_threshold_used": partial_threshold,
                    "strict_threshold_used": strict_threshold,
                    **posterior,
                }
            )
    metric_table = pd.DataFrame(rows)
    summary_rows: list[dict[str, object]] = []
    if not metric_table.empty:
        for key, part in metric_table.groupby(list(module_columns), observed=True):
            key_values = key if isinstance(key, tuple) else (key,)
            structural = part.loc[part["metric"].isin(["density", "connectivity"])]
            supporting = part.loc[part["metric"].isin(["program", "membership"])]
            partial_probability = float(
                structural["posterior_probability_partial"].min()
            ) if not structural.empty else np.nan
            positive_probability = float(part["posterior_probability_positive"].min())
            strict_probability = float(
                structural["posterior_probability_strict"].min()
            ) if not structural.empty else np.nan
            supporting_probability = float(
                supporting["posterior_probability_partial"].min()
            ) if not supporting.empty else 1.0
            summary_rows.append(
                {
                    **dict(zip(module_columns, key_values, strict=True)),
                    "n_metrics": len(part),
                    "conservative_probability_positive_all_metrics": positive_probability,
                    "conservative_probability_partial_preservation": partial_probability,
                    "conservative_probability_strict_preservation": strict_probability,
                    "conservative_probability_supporting_metrics_positive": supporting_probability,
                    "bayesian_preservation_class": (
                        "preserved"
                        if strict_probability >= 0.95 and supporting_probability >= 0.80
                        else "partially_preserved"
                        if partial_probability >= 0.80 and supporting_probability >= 0.80
                        else "uncertain_or_not_preserved"
                    ),
                    "decision_rule": (
                        "density/connectivity use WGCNA-like Z thresholds; program/"
                        "membership correlations are required to be positive support"
                    ),
                }
            )
    return BayesianPreservationResult(pd.DataFrame(summary_rows), metric_table)
