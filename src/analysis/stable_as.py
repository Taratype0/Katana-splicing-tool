from __future__ import annotations

import math
from dataclasses import dataclass

import pandas as pd

from src.domain.models.analysis import ThresholdConfig


STABILITY_METRIC_COLUMNS = [
    "mean_junction_reads_group1",
    "mean_junction_reads_group2",
    "replicate_count_group1",
    "replicate_count_group2",
    "supported_replicate_count_group1",
    "supported_replicate_count_group2",
    "replicate_direction_consistency",
    "stable_support_status",
]


@dataclass(frozen=True)
class StableAsMetrics:
    mean_junction_reads_group1: float
    mean_junction_reads_group2: float
    replicate_count_group1: int
    replicate_count_group2: int
    supported_replicate_count_group1: int
    supported_replicate_count_group2: int
    replicate_direction_consistency: float
    stable_support_status: str

    def as_dict(self) -> dict[str, object]:
        return {
            "mean_junction_reads_group1": self.mean_junction_reads_group1,
            "mean_junction_reads_group2": self.mean_junction_reads_group2,
            "replicate_count_group1": self.replicate_count_group1,
            "replicate_count_group2": self.replicate_count_group2,
            "supported_replicate_count_group1": self.supported_replicate_count_group1,
            "supported_replicate_count_group2": self.supported_replicate_count_group2,
            "replicate_direction_consistency": self.replicate_direction_consistency,
            "stable_support_status": self.stable_support_status,
        }


def parse_float_list(value: object) -> list[float]:
    if value is None:
        return []
    try:
        if pd.isna(value):
            return []
    except (TypeError, ValueError):
        pass
    values: list[float] = []
    for token in str(value).replace(";", ",").split(","):
        text = token.strip()
        if not text or text.upper() in {"NA", "NAN", "NONE"}:
            continue
        try:
            number = float(text)
        except ValueError:
            continue
        if math.isfinite(number):
            values.append(number)
    return values


def paired_read_totals(ijc_values: list[float], sjc_values: list[float]) -> list[float]:
    n_items = max(len(ijc_values), len(sjc_values))
    totals: list[float] = []
    for index in range(n_items):
        ijc = ijc_values[index] if index < len(ijc_values) else 0.0
        sjc = sjc_values[index] if index < len(sjc_values) else 0.0
        totals.append(float(ijc) + float(sjc))
    return totals


def pairwise_direction_consistency(
    psi_group1: list[float],
    psi_group2: list[float],
    standardized_dpsi: float,
) -> float:
    if not math.isfinite(standardized_dpsi) or standardized_dpsi == 0:
        return math.nan
    expected_positive = standardized_dpsi > 0
    differences = [
        left - right
        for left in psi_group1
        for right in psi_group2
        if math.isfinite(left) and math.isfinite(right) and left != right
    ]
    if not differences:
        return math.nan
    return sum((difference > 0) == expected_positive for difference in differences) / len(differences)


def stable_event_metrics(
    row: pd.Series,
    standardized_dpsi: float,
    *,
    reverse_splicing: bool,
    min_reads_per_replicate: float,
) -> dict[str, object]:
    psi_1 = parse_float_list(row.get("IncLevel1"))
    psi_2 = parse_float_list(row.get("IncLevel2"))
    ijc_1 = parse_float_list(row.get("IJC_SAMPLE_1"))
    ijc_2 = parse_float_list(row.get("IJC_SAMPLE_2"))
    sjc_1 = parse_float_list(row.get("SJC_SAMPLE_1"))
    sjc_2 = parse_float_list(row.get("SJC_SAMPLE_2"))
    if reverse_splicing:
        psi_1, psi_2 = psi_2, psi_1
        ijc_1, ijc_2 = ijc_2, ijc_1
        sjc_1, sjc_2 = sjc_2, sjc_1

    reads_1 = paired_read_totals(ijc_1, sjc_1)
    reads_2 = paired_read_totals(ijc_2, sjc_2)
    if not psi_1 or not psi_2:
        status = "psi_unavailable"
    elif not reads_1 or not reads_2:
        status = "read_support_unavailable"
    else:
        status = "available"
    metrics = StableAsMetrics(
        mean_junction_reads_group1=float(sum(reads_1) / len(reads_1)) if reads_1 else math.nan,
        mean_junction_reads_group2=float(sum(reads_2) / len(reads_2)) if reads_2 else math.nan,
        replicate_count_group1=max(len(psi_1), len(reads_1)),
        replicate_count_group2=max(len(psi_2), len(reads_2)),
        supported_replicate_count_group1=sum(value >= min_reads_per_replicate for value in reads_1),
        supported_replicate_count_group2=sum(value >= min_reads_per_replicate for value in reads_2),
        replicate_direction_consistency=pairwise_direction_consistency(psi_1, psi_2, standardized_dpsi),
        stable_support_status=status,
    )
    return metrics.as_dict()


def apply_stable_as_filter(
    frame: pd.DataFrame,
    *,
    fdr: pd.Series,
    standardized_dpsi: pd.Series,
    reverse_splicing: bool,
    thresholds: ThresholdConfig,
) -> pd.DataFrame:
    working = frame.copy()
    for column in STABILITY_METRIC_COLUMNS:
        if column not in working.columns:
            working[column] = "not_evaluated" if column == "stable_support_status" else math.nan
    base_sig = fdr.le(thresholds.splicing_fdr) & standardized_dpsi.abs().ge(thresholds.splicing_dpsi)
    working["rmats_sig_before_stability"] = base_sig.fillna(False)
    significant_index = working.index[working["rmats_sig_before_stability"]]
    if len(significant_index):
        metric_rows = working.loc[significant_index].apply(
            lambda row: stable_event_metrics(
                row,
                float(standardized_dpsi.loc[row.name]) if pd.notna(standardized_dpsi.loc[row.name]) else math.nan,
                reverse_splicing=reverse_splicing,
                min_reads_per_replicate=float(thresholds.stable_min_junction_reads),
            ),
            axis=1,
        )
        metric_frame = pd.DataFrame(metric_rows.tolist(), index=metric_rows.index)
        for column in STABILITY_METRIC_COLUMNS:
            if column in metric_frame.columns:
                working.loc[metric_frame.index, column] = metric_frame[column]
    working["stable_as_event"] = (
        working["rmats_sig_before_stability"]
        & working["stable_support_status"].eq("available")
        & pd.to_numeric(working["supported_replicate_count_group1"], errors="coerce").ge(
            thresholds.stable_min_replicates_per_group
        )
        & pd.to_numeric(working["supported_replicate_count_group2"], errors="coerce").ge(
            thresholds.stable_min_replicates_per_group
        )
        & pd.to_numeric(working["replicate_direction_consistency"], errors="coerce").ge(
            thresholds.stable_min_pairwise_consistency
        )
    )
    return working
