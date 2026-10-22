from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from katana_cosplice.config import StableASConfig


@dataclass
class StableASResult:
    metrics: pd.DataFrame
    stable_psi: pd.DataFrame
    condition_metrics: pd.DataFrame
    summary: pd.DataFrame


def _validate_inputs(
    psi: pd.DataFrame,
    ijc: pd.DataFrame,
    sjc: pd.DataFrame,
    samples: pd.DataFrame,
) -> None:
    if not psi.index.is_unique or not psi.columns.is_unique:
        raise ValueError("PSI matrix event and sample labels must be unique.")
    if not psi.index.equals(ijc.index) or not psi.index.equals(sjc.index):
        raise ValueError("PSI, IJC, and SJC matrices must contain the same ordered event_key rows.")
    if not psi.columns.equals(ijc.columns) or not psi.columns.equals(sjc.columns):
        raise ValueError("PSI, IJC, and SJC matrices must contain the same ordered sample columns.")
    if not {"sample_id", "condition"}.issubset(samples.columns):
        raise ValueError("Sample metadata must contain sample_id and condition columns.")
    if not samples["sample_id"].is_unique:
        raise ValueError("Sample metadata sample_id values must be unique.")
    metadata_samples = samples["sample_id"].tolist()
    if set(metadata_samples) != set(psi.columns):
        raise ValueError(
            "Sample metadata and matrix columns differ: "
            f"metadata_only={sorted(set(metadata_samples)-set(psi.columns))}, "
            f"matrix_only={sorted(set(psi.columns)-set(metadata_samples))}"
        )
    # Range checks need numeric values, not a large temporary stacked index.
    psi_values = psi.to_numpy(copy=False)
    measured_psi = psi_values[~pd.isna(psi_values)]
    if np.any(measured_psi < 0.0) or np.any(measured_psi > 1.0):
        raise ValueError("PSI values must be within [0, 1].")
    for counts in (ijc, sjc):
        values = counts.to_numpy(copy=False)
        if np.any(values[~pd.isna(values)] < 0):
            raise ValueError("IJC and SJC values cannot be negative.")


def _condition_metrics(
    psi: pd.DataFrame,
    total_reads: pd.DataFrame,
    samples: pd.DataFrame,
    min_reads: float,
    min_replicates: int,
    event_keys: pd.Index,
) -> pd.DataFrame:
    psi = psi.loc[event_keys]
    total_reads = total_reads.loc[event_keys]
    frames: list[pd.DataFrame] = []
    for condition, subset in samples.groupby("condition", sort=False, observed=True):
        columns = [sample for sample in subset["sample_id"] if sample in psi.columns]
        values = psi.loc[:, columns]
        reads = total_reads.loc[:, columns]
        measured = values.notna()
        supported = measured & reads.ge(min_reads)
        frame = pd.DataFrame(
            {
                "event_key": psi.index,
                "condition": str(condition),
                "n_replicates": len(columns),
                "n_measured": measured.sum(axis=1).to_numpy(),
                "n_read_supported": supported.sum(axis=1).to_numpy(),
                "mean_psi": values.mean(axis=1, skipna=True).to_numpy(),
                "sd_psi": values.std(axis=1, ddof=1, skipna=True).to_numpy(),
                "median_total_junction_reads": reads.where(measured).median(axis=1, skipna=True).to_numpy(),
            }
        )
        frame["passes_min_supported_replicates"] = frame["n_read_supported"].ge(min_replicates)
        frames.append(frame)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _condition_moments(matrix: pd.DataFrame, samples: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series]:
    condition_means: dict[str, pd.Series] = {}
    within_variances: list[pd.Series] = []
    for condition, subset in samples.groupby("condition", sort=False, observed=True):
        columns = [sample for sample in subset["sample_id"] if sample in matrix.columns]
        if not columns:
            continue
        values = matrix.loc[:, columns]
        condition_means[str(condition)] = values.mean(axis=1, skipna=True)
        within_variances.append(values.var(axis=1, ddof=1, skipna=True))
    means = pd.DataFrame(condition_means, index=matrix.index)
    within = pd.concat(within_variances, axis=1).mean(axis=1, skipna=True) if within_variances else pd.Series(np.nan, index=matrix.index)
    return means, within


def _condition_support_fraction(
    psi: pd.DataFrame,
    total_reads: pd.DataFrame,
    samples: pd.DataFrame,
    min_reads: float,
    min_replicates: int,
) -> pd.Series:
    supported: list[pd.Series] = []
    for _, subset in samples.groupby("condition", sort=False, observed=True):
        columns = [sample for sample in subset["sample_id"] if sample in psi.columns]
        if not columns:
            continue
        measured = psi.loc[:, columns].notna() & total_reads.loc[:, columns].ge(min_reads)
        supported.append(measured.sum(axis=1).ge(min_replicates))
    if not supported:
        return pd.Series(0.0, index=psi.index)
    return pd.concat(supported, axis=1).mean(axis=1)


def stable_as_qc(
    psi: pd.DataFrame,
    ijc: pd.DataFrame,
    sjc: pd.DataFrame,
    samples: pd.DataFrame,
    settings: StableASConfig,
) -> StableASResult:
    _validate_inputs(psi, ijc, sjc, samples)
    thresholds = settings.resolved()
    total_reads = ijc.fillna(0.0) + sjc.fillna(0.0)
    valid_fraction = psi.notna().mean(axis=1)
    psi_variance = psi.var(axis=1, ddof=1, skipna=True)
    mean_total_reads = total_reads.where(psi.notna()).mean(axis=1, skipna=True)
    read_supported_fraction = total_reads.ge(thresholds["min_total_junction_reads"]).mean(axis=1)

    condition_means, within_variance = _condition_moments(psi, samples)
    condition_supported_fraction = _condition_support_fraction(
        psi,
        total_reads,
        samples,
        thresholds["min_total_junction_reads"],
        int(thresholds["min_supported_replicates_per_condition"]),
    )
    condition_mean_psi_range = condition_means.max(axis=1, skipna=True) - condition_means.min(axis=1, skipna=True)
    alternative_psi_fraction = (
        psi.ge(thresholds["alternative_psi_low"])
        & psi.le(thresholds["alternative_psi_high"])
    ).mean(axis=1)
    sample_psi_range = psi.max(axis=1, skipna=True) - psi.min(axis=1, skipna=True)
    between_variance = condition_means.var(axis=1, ddof=1, skipna=True)
    repeatability = between_variance / (between_variance + within_variance + 1e-12)
    repeatability = repeatability.clip(0.0, 1.0).fillna(0.0)

    expected_binomial = (psi * (1.0 - psi) / total_reads.where(total_reads.gt(0))).replace([np.inf, -np.inf], np.nan)
    expected_within_variance = expected_binomial.mean(axis=1, skipna=True)
    overdispersion_ratio = (within_variance / (expected_within_variance + 1e-12)).replace([np.inf, -np.inf], np.nan)

    metrics = pd.DataFrame(
        {
            "valid_fraction": valid_fraction,
            "mean_total_junction_reads": mean_total_reads,
            "read_supported_fraction": read_supported_fraction,
            "condition_supported_fraction": condition_supported_fraction,
            "psi_variance": psi_variance,
            "condition_mean_psi_range": condition_mean_psi_range,
            "alternative_psi_fraction": alternative_psi_fraction,
            "sample_psi_range": sample_psi_range,
            "between_condition_variance": between_variance,
            "within_condition_variance": within_variance,
            "replicate_repeatability": repeatability,
            "expected_binomial_variance": expected_within_variance,
            "overdispersion_ratio": overdispersion_ratio,
        },
        index=psi.index,
    )
    metrics.index.name = "event_key"
    metrics["passes_valid_fraction"] = metrics["valid_fraction"].ge(thresholds["min_valid_fraction"])
    metrics["passes_read_support"] = metrics["read_supported_fraction"].ge(thresholds["min_read_supported_fraction"])
    metrics["passes_condition_support"] = metrics["condition_supported_fraction"].ge(thresholds["min_condition_supported_fraction"])
    metrics["passes_variance"] = metrics["psi_variance"].ge(thresholds["min_psi_variance"])
    metrics["passes_repeatability"] = metrics["replicate_repeatability"].ge(thresholds["min_repeatability"])
    metrics["passes_dynamic_range"] = (
        metrics["condition_mean_psi_range"].ge(thresholds["min_condition_mean_psi_range"])
        & (
            metrics["alternative_psi_fraction"].ge(thresholds["min_alternative_fraction"])
            | metrics["sample_psi_range"].ge(thresholds["min_sample_psi_range"])
        )
    )
    pass_columns = [
        "passes_valid_fraction",
        "passes_read_support",
        "passes_condition_support",
        "passes_variance",
        "passes_repeatability",
        "passes_dynamic_range",
    ]
    metrics["stable_as"] = metrics[pass_columns].all(axis=1)
    metrics["stability_score"] = (
        metrics["valid_fraction"].clip(0, 1)
        * metrics["read_supported_fraction"].clip(0, 1)
        * np.sqrt(metrics["replicate_repeatability"].clip(0, 1))
        * np.log1p(metrics["mean_total_junction_reads"].clip(lower=0))
        * np.sqrt(metrics["psi_variance"].clip(lower=0))
    )
    metrics["failed_rules"] = metrics.apply(
        lambda row: ",".join(
            label
            for label, column in [
                ("missingness", "passes_valid_fraction"),
                ("read_support", "passes_read_support"),
                ("condition_coverage", "passes_condition_support"),
                ("low_variance", "passes_variance"),
                ("replicate_instability", "passes_repeatability"),
                ("insufficient_dynamic_range", "passes_dynamic_range"),
            ]
            if not bool(row[column])
        ),
        axis=1,
    )

    stable_keys = metrics.index[metrics["stable_as"]]
    stable_psi = psi.loc[stable_keys].copy()
    stable_psi = stable_psi.T.fillna(stable_psi.median(axis=1)).T
    condition_metrics = _condition_metrics(
        psi,
        total_reads,
        samples,
        thresholds["min_total_junction_reads"],
        int(thresholds["min_supported_replicates_per_condition"]),
        stable_keys,
    )
    summary_rows = [
        {
            "rule": column,
            "n_passed": int(metrics[column].sum()),
            "n_failed": int((~metrics[column]).sum()),
            "pass_fraction": float(metrics[column].mean()),
        }
        for column in pass_columns
    ]
    summary_rows.append(
        {
            "rule": "stable_as",
            "n_passed": int(metrics["stable_as"].sum()),
            "n_failed": int((~metrics["stable_as"]).sum()),
            "pass_fraction": float(metrics["stable_as"].mean()),
        }
    )
    return StableASResult(
        metrics=metrics.sort_values("stability_score", ascending=False),
        stable_psi=stable_psi,
        condition_metrics=condition_metrics,
        summary=pd.DataFrame(summary_rows),
    )
