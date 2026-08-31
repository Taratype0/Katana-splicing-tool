from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

from katana_cosplice.network import correlation_matrix


@dataclass
class PosteriorCorrelationResult:
    mean_correlation: np.ndarray
    correlation_sd: np.ndarray
    sign_probability: np.ndarray
    event_keys: np.ndarray
    n_draws: int
    method: str
    random_seed: int


def posterior_ensemble_correlation(
    posterior_alpha: pd.DataFrame,
    posterior_beta: pd.DataFrame,
    *,
    n_draws: int = 25,
    method: str = "bicor",
    random_seed: int = 17,
    draw_transform: Callable[[pd.DataFrame], pd.DataFrame] | None = None,
) -> PosteriorCorrelationResult:
    """Propagate beta-posterior PSI uncertainty into event correlations.

    Correlations are averaged on the Fisher-z scale. The returned standard
    deviation and sign probability quantify edge uncertainty across posterior
    PSI draws; they are not differential-splicing p-values.
    """

    if n_draws < 2:
        raise ValueError("posterior ensemble requires at least two draws")
    if not posterior_alpha.index.equals(posterior_beta.index) or not posterior_alpha.columns.equals(
        posterior_beta.columns
    ):
        raise ValueError("posterior alpha and beta matrices must have identical axes")
    alpha = posterior_alpha.to_numpy(dtype=float)
    beta = posterior_beta.to_numpy(dtype=float)
    if not np.isfinite(alpha).all() or not np.isfinite(beta).all():
        raise ValueError("posterior alpha and beta matrices must be finite")
    if np.any(alpha <= 0) or np.any(beta <= 0):
        raise ValueError("posterior alpha and beta parameters must be positive")

    rng = np.random.default_rng(random_seed)
    shape = (len(posterior_alpha), len(posterior_alpha))
    mean_z = np.zeros(shape, dtype=np.float64)
    m2_z = np.zeros(shape, dtype=np.float64)
    positive = np.zeros(shape, dtype=np.uint32)
    template = posterior_alpha.copy()
    for draw_index in range(1, n_draws + 1):
        template.iloc[:, :] = rng.beta(alpha, beta)
        draw = draw_transform(template.copy()) if draw_transform is not None else template
        correlation = correlation_matrix(draw, method)
        positive += correlation >= 0
        z = np.arctanh(np.clip(correlation, -0.999999, 0.999999))
        delta = z - mean_z
        mean_z += delta / draw_index
        m2_z += delta * (z - mean_z)

    mean_correlation = np.tanh(mean_z)
    z_sd = np.sqrt(m2_z / (n_draws - 1))
    # Delta-method conversion keeps the uncertainty scale interpretable as r.
    correlation_sd = z_sd * (1.0 - mean_correlation**2)
    positive_probability = positive.astype(np.float64) / n_draws
    sign_probability = np.maximum(positive_probability, 1.0 - positive_probability)
    for matrix, diagonal in (
        (mean_correlation, 1.0),
        (correlation_sd, 0.0),
        (sign_probability, 1.0),
    ):
        matrix[:] = (matrix + matrix.T) / 2.0
        np.fill_diagonal(matrix, diagonal)
    return PosteriorCorrelationResult(
        mean_correlation=mean_correlation.astype(np.float32),
        correlation_sd=correlation_sd.astype(np.float32),
        sign_probability=sign_probability.astype(np.float32),
        event_keys=posterior_alpha.index.to_numpy(dtype=str),
        n_draws=n_draws,
        method=method,
        random_seed=random_seed,
    )


def save_posterior_correlation(result: PosteriorCorrelationResult, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        mean_correlation=result.mean_correlation,
        correlation_sd=result.correlation_sd,
        sign_probability=result.sign_probability,
        event_keys=result.event_keys,
        n_draws=np.asarray([result.n_draws], dtype=np.int32),
        method=np.asarray([result.method]),
        random_seed=np.asarray([result.random_seed], dtype=np.int64),
    )


def load_posterior_correlation(path: Path, event_keys: pd.Index) -> PosteriorCorrelationResult:
    payload = np.load(path, allow_pickle=False)
    cached_keys = payload["event_keys"].astype(str)
    requested = event_keys.astype(str).to_numpy()
    if not np.array_equal(cached_keys, requested):
        raise ValueError("posterior correlation cache event order does not match the network matrix")
    return PosteriorCorrelationResult(
        mean_correlation=payload["mean_correlation"],
        correlation_sd=payload["correlation_sd"],
        sign_probability=payload["sign_probability"],
        event_keys=cached_keys,
        n_draws=int(payload["n_draws"][0]),
        method=str(payload["method"][0]),
        random_seed=int(payload["random_seed"][0]),
    )
