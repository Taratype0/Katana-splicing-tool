from __future__ import annotations

from dataclasses import asdict, dataclass
import json

import numpy as np
import pandas as pd

from katana_cosplice.config import StableASConfig


NETWORK_INPUT_POLICY_VERSION = "network-all-samples-counts-v1"
NETWORK_MIN_READS = 10.0
NETWORK_NONCONSTANT_RULE = "sample_max_gt_min"


def network_input_policy() -> dict[str, object]:
    """The active input policy, separate from historical diagnostic settings."""
    return {
        "policy_id": NETWORK_INPUT_POLICY_VERSION,
        "min_count": NETWORK_MIN_READS,
        "scope": "every declared sample in the complete network cohort",
        "count_definition": "IJC + SJC, both finite and nonnegative",
        "psi_definition": "finite PSI within [0, 1] for every declared sample",
        "missing_or_invalid": "exclude the event; do not impute or shrink the cohort",
        "constant_criterion": NETWORK_NONCONSTANT_RULE,
        "differential_preselection": False,
        "legacy_quality_thresholds": "diagnostic only; never determine eligibility",
        "stable_as_alias": "network_input_eligible; not replicate stability",
        "optional_cap_ranking": "cohort_PSI_variance; condition independent",
    }


@dataclass
class StableASResult:
    metrics: pd.DataFrame
    stable_psi: pd.DataFrame
    condition_metrics: pd.DataFrame
    summary: pd.DataFrame


def _validate_labels(labels: pd.Index | pd.Series, label: str) -> None:
    if pd.isna(labels).any():
        raise ValueError(f"{label} labels cannot be missing.")
    normalized = pd.Index([str(value).strip() for value in labels])
    if (normalized == "").any():
        raise ValueError(f"{label} labels cannot be empty.")
    if not normalized.is_unique:
        raise ValueError(f"{label} labels must be unique, including after whitespace normalization.")


def _validate_inputs(
    psi: pd.DataFrame,
    ijc: pd.DataFrame,
    sjc: pd.DataFrame,
    samples: pd.DataFrame,
) -> None:
    for name, frame in (("PSI", psi), ("IJC", ijc), ("SJC", sjc)):
        _validate_labels(frame.index, f"{name} event")
        _validate_labels(frame.columns, f"{name} sample")
    if not psi.index.equals(ijc.index) or not psi.index.equals(sjc.index):
        raise ValueError("PSI, IJC, and SJC matrices must contain the same ordered event_key rows.")
    if not psi.columns.equals(ijc.columns) or not psi.columns.equals(sjc.columns):
        raise ValueError("PSI, IJC, and SJC matrices must contain the same ordered sample columns.")
    if not {"sample_id", "condition"}.issubset(samples.columns):
        raise ValueError("Sample metadata must contain sample_id and condition columns.")
    _validate_labels(samples["sample_id"], "Sample metadata sample_id")
    if len(psi.columns) == 0:
        raise ValueError("The network cohort must contain declared samples.")
    metadata_samples = samples["sample_id"].tolist()
    if set(metadata_samples) != set(psi.columns):
        raise ValueError(
            "Sample metadata and matrix columns differ: "
            f"metadata_only={sorted(set(metadata_samples)-set(psi.columns), key=str)}, "
            f"matrix_only={sorted(set(psi.columns)-set(metadata_samples), key=str)}"
        )
    if samples["condition"].isna().any() or samples["condition"].astype(str).str.strip().eq("").any():
        raise ValueError("Sample metadata condition values cannot be missing or empty.")


def _numeric_frame(frame: pd.DataFrame) -> pd.DataFrame:
    # Mark malformed observations unusable without modifying the native matrix.
    numeric = frame.apply(pd.to_numeric, errors="coerce")
    if np.iscomplexobj(numeric.to_numpy()):
        complex_cells = frame.map(lambda value: isinstance(value, (complex, np.complexfloating)))
        numeric = numeric.where(~complex_cells, np.nan)
        numeric = pd.DataFrame(np.real(numeric.to_numpy()), index=frame.index, columns=frame.columns)
    return numeric.astype(float)


def _condition_metrics(
    psi: pd.DataFrame,
    total_reads: pd.DataFrame,
    usable: pd.DataFrame,
    supported: pd.DataFrame,
    samples: pd.DataFrame,
    event_keys: pd.Index,
) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for condition, subset in samples.groupby("condition", sort=False, observed=True):
        columns = subset["sample_id"].tolist()
        values = psi.loc[event_keys, columns]
        reads = total_reads.loc[event_keys, columns]
        frame = pd.DataFrame(
            {
                "event_key": event_keys,
                "condition": str(condition),
                "n_replicates": len(columns),
                "n_measured": usable.loc[event_keys, columns].sum(axis=1).to_numpy(),
                "n_read_supported": supported.loc[event_keys, columns].sum(axis=1).to_numpy(),
                "mean_psi": values.mean(axis=1, skipna=True).to_numpy(),
                "sd_psi": values.std(axis=1, ddof=1, skipna=True).to_numpy(),
                "median_total_junction_reads": reads.median(axis=1, skipna=True).to_numpy(),
            }
        )
        frame["passes_all_samples"] = frame["n_read_supported"].eq(len(columns))
        # Compatibility field: all declared samples, never a minimum of two.
        frame["passes_min_supported_replicates"] = frame["passes_all_samples"]
        frame["input_policy"] = NETWORK_INPUT_POLICY_VERSION
        frame["usage"] = "diagnostic; eligibility is evaluated over the complete cohort"
        frames.append(frame)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _condition_moments(matrix: pd.DataFrame, samples: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series]:
    condition_means: dict[str, pd.Series] = {}
    within_variances: list[pd.Series] = []
    for condition, subset in samples.groupby("condition", sort=False, observed=True):
        values = matrix.loc[:, subset["sample_id"].tolist()]
        condition_means[str(condition)] = values.mean(axis=1, skipna=True)
        within_variances.append(values.var(axis=1, ddof=1, skipna=True))
    means = pd.DataFrame(condition_means, index=matrix.index)
    within = pd.concat(within_variances, axis=1).mean(axis=1, skipna=True) if within_variances else pd.Series(np.nan, index=matrix.index)
    return means, within


def stable_as_qc(
    psi: pd.DataFrame,
    ijc: pd.DataFrame,
    sjc: pd.DataFrame,
    samples: pd.DataFrame,
    settings: StableASConfig,
) -> StableASResult:
    """Select complete, read-supported, nonconstant events without group preselection.

    The function/result names remain compatible with older projects. Historical
    preset and override values are recorded but do not control the active policy.
    """
    _validate_inputs(psi, ijc, sjc, samples)
    numeric_psi = _numeric_frame(psi)
    numeric_ijc = _numeric_frame(ijc)
    numeric_sjc = _numeric_frame(sjc)
    valid_psi = pd.DataFrame(
        np.isfinite(numeric_psi.to_numpy()) & numeric_psi.ge(0.0).to_numpy() & numeric_psi.le(1.0).to_numpy(),
        index=psi.index,
        columns=psi.columns,
    )
    valid_counts = pd.DataFrame(
        np.isfinite(numeric_ijc.to_numpy()) & np.isfinite(numeric_sjc.to_numpy())
        & numeric_ijc.ge(0.0).to_numpy() & numeric_sjc.ge(0.0).to_numpy(),
        index=psi.index,
        columns=psi.columns,
    )
    with np.errstate(over="ignore", invalid="ignore"):
        raw_totals = numeric_ijc + numeric_sjc
    valid_counts &= np.isfinite(raw_totals)
    total_reads = raw_totals.where(valid_counts)
    count_supported = valid_counts & total_reads.ge(NETWORK_MIN_READS)
    supported = valid_psi & count_supported
    usable_psi = numeric_psi.where(valid_psi)
    usable = valid_psi & valid_counts

    condition_means, within_variance = _condition_moments(usable_psi, samples)
    condition_mean_psi_range = condition_means.max(axis=1, skipna=True) - condition_means.min(axis=1, skipna=True)
    between_variance = condition_means.var(axis=1, ddof=1, skipna=True)
    repeatability = (between_variance / (between_variance + within_variance + 1e-12)).clip(0.0, 1.0).fillna(0.0)
    psi_variance = usable_psi.var(axis=1, ddof=1, skipna=True)
    sample_psi_range = usable_psi.max(axis=1, skipna=True) - usable_psi.min(axis=1, skipna=True)
    mean_total_reads = total_reads.where(valid_psi).mean(axis=1, skipna=True)
    expected_binomial = (usable_psi * (1.0 - usable_psi) / total_reads.where(total_reads.gt(0))).replace([np.inf, -np.inf], np.nan)
    expected_within_variance = expected_binomial.mean(axis=1, skipna=True)
    overdispersion_ratio = (within_variance / (expected_within_variance + 1e-12)).replace([np.inf, -np.inf], np.nan)
    condition_support = [
        supported.loc[:, subset["sample_id"].tolist()].all(axis=1)
        for _, subset in samples.groupby("condition", sort=False, observed=True)
    ]
    condition_supported_fraction = pd.concat(condition_support, axis=1).mean(axis=1)

    metrics = pd.DataFrame(
        {
            "input_policy": NETWORK_INPUT_POLICY_VERSION,
            "min_total_junction_reads": NETWORK_MIN_READS,
            "expected_samples": len(psi.columns),
            "n_valid_psi_samples": valid_psi.sum(axis=1),
            "n_valid_read_count_samples": valid_counts.sum(axis=1),
            "n_read_supported_samples": count_supported.sum(axis=1),
            "n_usable_read_supported_samples": supported.sum(axis=1),
            "valid_fraction": valid_psi.mean(axis=1),
            "mean_total_junction_reads": mean_total_reads,
            "read_supported_fraction": count_supported.mean(axis=1),
            "condition_supported_fraction": condition_supported_fraction,
            "psi_variance": psi_variance,
            "condition_mean_psi_range": condition_mean_psi_range,
            "alternative_psi_fraction": (usable_psi.ge(0.20) & usable_psi.le(0.80)).mean(axis=1),
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
    metrics["passes_complete_psi"] = valid_psi.all(axis=1)
    metrics["passes_complete_read_counts"] = valid_counts.all(axis=1)
    metrics["passes_read_support"] = count_supported.all(axis=1)
    # No variance floor or numerical epsilon: strictly unequal valid values suffice.
    metrics["passes_nonconstant"] = sample_psi_range.gt(0.0) & metrics["passes_complete_psi"]
    pass_columns = ["passes_complete_psi", "passes_complete_read_counts", "passes_read_support", "passes_nonconstant"]
    metrics["network_input_eligible"] = metrics[pass_columns].all(axis=1)
    metrics["stable_as"] = metrics["network_input_eligible"]
    metrics["legacy_stability_score_diagnostic"] = (
        metrics["valid_fraction"].clip(0, 1)
        * metrics["read_supported_fraction"].clip(0, 1)
        * np.sqrt(metrics["replicate_repeatability"].clip(0, 1))
        * np.log1p(metrics["mean_total_junction_reads"].clip(lower=0))
        * np.sqrt(metrics["psi_variance"].clip(lower=0))
    ).fillna(0.0)
    # An explicit optional event cap must not select events by known conditions.
    metrics["stability_score"] = metrics["psi_variance"].fillna(0.0)
    metrics["stable_as_score"] = metrics["stability_score"]
    metrics["diagnostic_only_metrics"] = (
        "valid_fraction,read_supported_fraction,condition_supported_fraction,psi_variance,"
        "condition_mean_psi_range,alternative_psi_fraction,sample_psi_range,"
        "replicate_repeatability,expected_binomial_variance,overdispersion_ratio,"
        "legacy_stability_score_diagnostic"
    )
    metrics["ignored_legacy_qc_settings"] = json.dumps(asdict(settings), sort_keys=True, allow_nan=False)
    metrics["failed_rules"] = [
        ",".join(
            label for label, passed in (
                ("invalid_or_missing_psi", bool(row.passes_complete_psi)),
                ("invalid_or_missing_counts", bool(row.passes_complete_read_counts)),
                ("sample_below_10_reads", bool(row.passes_read_support) or not bool(row.passes_complete_read_counts)),
                ("constant_psi", bool(row.passes_nonconstant) or not bool(row.passes_complete_psi)),
            ) if not passed
        )
        for row in metrics.itertuples()
    ]

    eligible_keys = metrics.index[metrics["network_input_eligible"]]
    # Eligible data are complete. Never impute observations or reorder samples.
    stable_psi = numeric_psi.loc[eligible_keys].copy()
    condition_metrics = _condition_metrics(
        usable_psi, total_reads, usable, supported, samples, eligible_keys,
    )
    summary_rows = [
        {
            "rule": column,
            "n_passed": int(metrics[column].sum()),
            "n_failed": int((~metrics[column]).sum()),
            "pass_fraction": float(metrics[column].mean()),
            "input_policy": NETWORK_INPUT_POLICY_VERSION,
        }
        for column in [*pass_columns, "network_input_eligible", "stable_as"]
    ]
    return StableASResult(
        metrics=metrics.sort_values("stability_score", ascending=False, kind="stable"),
        stable_psi=stable_psi,
        condition_metrics=condition_metrics,
        summary=pd.DataFrame(summary_rows),
    )
