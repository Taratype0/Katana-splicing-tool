from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass
class UncertaintyResult:
    posterior_psi: pd.DataFrame
    posterior_alpha: pd.DataFrame
    posterior_beta: pd.DataFrame
    event_metrics: pd.DataFrame
    sample_metrics: pd.DataFrame


def _condition_rho(
    successes: np.ndarray,
    totals: np.ndarray,
    condition_values: np.ndarray,
) -> float:
    """Estimate beta-binomial overdispersion after removing condition means."""
    numerators: list[float] = []
    denominators: list[float] = []
    for condition in pd.unique(condition_values):
        mask = (condition_values == condition) & np.isfinite(successes) & (totals > 0)
        if mask.sum() < 2:
            continue
        n = totals[mask].astype(float)
        p = successes[mask] / n
        mean = float(np.average(p, weights=n))
        if mean <= 1e-6 or mean >= 1.0 - 1e-6:
            continue
        observed = float(np.average((p - mean) ** 2, weights=n))
        binomial = float(np.average(mean * (1.0 - mean) / n, weights=n))
        beta_component = float(
            np.average(mean * (1.0 - mean) * (1.0 - 1.0 / n), weights=n)
        )
        if beta_component > 1e-12:
            numerators.append(max(0.0, observed - binomial))
            denominators.append(beta_component)
    if not denominators:
        return 0.0
    return float(np.clip(np.sum(numerators) / np.sum(denominators), 0.0, 0.95))


def estimate_beta_binomial_uncertainty(
    ijc: pd.DataFrame,
    sjc: pd.DataFrame,
    samples: pd.DataFrame,
    qc_metrics: pd.DataFrame | None = None,
    *,
    prior_strength: float = 1.0,
    reliability_read_scale: float = 20.0,
) -> UncertaintyResult:
    """Estimate read-depth and replicate-aware PSI uncertainty.

    A condition-adjusted method-of-moments beta-binomial rho estimate inflates
    binomial sampling variance. The resulting effective read count is used in a
    Jeffreys-like beta posterior. This is an uncertainty model for network
    weighting, not a differential-splicing significance test.
    """
    if not ijc.index.equals(sjc.index) or not ijc.columns.equals(sjc.columns):
        raise ValueError("IJC and SJC matrices must have identical event and sample axes.")
    if not {"sample_id", "condition"}.issubset(samples.columns):
        raise ValueError("Sample metadata must contain sample_id and condition.")
    ordered = samples.set_index("sample_id").reindex(ijc.columns)
    if ordered["condition"].isna().any():
        missing = ordered.index[ordered["condition"].isna()].tolist()
        raise ValueError(f"Missing sample condition metadata: {missing[:10]}")
    if prior_strength <= 0 or reliability_read_scale <= 0:
        raise ValueError("prior_strength and reliability_read_scale must be positive.")

    inclusion = ijc.to_numpy(dtype=float)
    skipping = sjc.to_numpy(dtype=float)
    totals = inclusion + skipping
    if np.nanmin(inclusion) < 0 or np.nanmin(skipping) < 0:
        raise ValueError("IJC and SJC counts cannot be negative.")
    conditions = ordered["condition"].astype(str).to_numpy()
    rho = np.asarray(
        [_condition_rho(row, total, conditions) for row, total in zip(inclusion, totals)],
        dtype=float,
    )
    inflation = 1.0 + np.maximum(totals - 1.0, 0.0) * rho[:, None]
    effective_reads = np.divide(
        totals,
        inflation,
        out=np.zeros_like(totals, dtype=float),
        where=np.isfinite(totals) & (inflation > 0),
    )
    raw_psi = np.divide(
        inclusion,
        totals,
        out=np.full_like(inclusion, 0.5, dtype=float),
        where=totals > 0,
    )
    alpha = prior_strength / 2.0 + raw_psi * effective_reads
    beta = prior_strength / 2.0 + (1.0 - raw_psi) * effective_reads
    posterior_psi = alpha / (alpha + beta)
    posterior_variance = alpha * beta / (
        (alpha + beta) ** 2 * (alpha + beta + 1.0)
    )
    read_precision = effective_reads / (effective_reads + reliability_read_scale)

    repeatability = np.ones(len(ijc), dtype=float)
    if qc_metrics is not None:
        aligned = qc_metrics.reindex(ijc.index)
        for column in (
            "replicate_repeatability",
            "repeatability",
            "condition_repeatability",
            "stability_score",
        ):
            if column in aligned:
                candidate = pd.to_numeric(aligned[column], errors="coerce").fillna(0.0)
                repeatability = candidate.clip(0.0, 1.0).to_numpy(dtype=float)
                break
    event_reliability = (
        np.nanmedian(read_precision, axis=1)
        * (1.0 - rho)
        * np.sqrt(np.clip(repeatability, 0.0, 1.0))
    )
    event_reliability = np.clip(event_reliability, 0.0, 1.0)

    posterior_frame = pd.DataFrame(
        posterior_psi,
        index=ijc.index,
        columns=ijc.columns,
    )
    event_metrics = pd.DataFrame(
        {
            "event_key": ijc.index,
            "beta_binomial_rho": rho,
            "median_total_reads": np.nanmedian(totals, axis=1),
            "median_effective_reads": np.nanmedian(effective_reads, axis=1),
            "median_posterior_sd": np.sqrt(np.nanmedian(posterior_variance, axis=1)),
            "median_read_precision": np.nanmedian(read_precision, axis=1),
            "replicate_consistency_component": repeatability,
            "network_reliability": event_reliability,
            "uncertainty_model": "condition_adjusted_beta_binomial_moments",
        }
    ).set_index("event_key")
    sample_metrics = pd.DataFrame(
        {
            "event_key": np.repeat(ijc.index.to_numpy(), len(ijc.columns)),
            "sample_id": np.tile(ijc.columns.to_numpy(), len(ijc.index)),
            "condition": np.tile(conditions, len(ijc.index)),
            "total_reads": totals.ravel(),
            "effective_reads": effective_reads.ravel(),
            "posterior_psi": posterior_psi.ravel(),
            "posterior_sd": np.sqrt(posterior_variance).ravel(),
            "read_precision": read_precision.ravel(),
        }
    )
    alpha_frame = pd.DataFrame(alpha, index=ijc.index, columns=ijc.columns)
    beta_frame = pd.DataFrame(beta, index=ijc.index, columns=ijc.columns)
    return UncertaintyResult(
        posterior_frame,
        alpha_frame,
        beta_frame,
        event_metrics,
        sample_metrics,
    )
