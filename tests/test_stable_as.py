from __future__ import annotations

import math

import pandas as pd

from src.analysis.stable_as import apply_stable_as_filter, pairwise_direction_consistency
from src.domain.models.analysis import ThresholdConfig


def _thresholds() -> ThresholdConfig:
    return ThresholdConfig(
        splicing_fdr=0.05,
        splicing_dpsi=0.10,
        stable_min_junction_reads=20,
        stable_min_replicates_per_group=2,
        stable_min_pairwise_consistency=0.67,
    )


def test_direction_consistency_is_independent_of_replicate_order() -> None:
    group_1 = [0.8, 0.7, 0.6]
    group_2 = [0.2, 0.3, 0.4]
    expected = pairwise_direction_consistency(group_1, group_2, 0.4)
    reordered = pairwise_direction_consistency(group_1[::-1], [0.3, 0.4, 0.2], 0.4)
    assert expected == 1.0
    assert reordered == expected


def test_stable_filter_requires_read_support_in_individual_replicates() -> None:
    frame = pd.DataFrame(
        {
            "FDR": [0.001],
            "IncLevelDifference": [0.4],
            "IncLevel1": ["0.8,0.7,0.6"],
            "IncLevel2": ["0.2,0.3,0.4"],
            "IJC_SAMPLE_1": ["60,0,0"],
            "SJC_SAMPLE_1": ["0,0,0"],
            "IJC_SAMPLE_2": ["20,20,20"],
            "SJC_SAMPLE_2": ["0,0,0"],
        }
    )
    result = apply_stable_as_filter(
        frame,
        fdr=frame["FDR"],
        standardized_dpsi=frame["IncLevelDifference"],
        reverse_splicing=False,
        thresholds=_thresholds(),
    )
    assert bool(result.loc[0, "rmats_sig_before_stability"])
    assert result.loc[0, "supported_replicate_count_group1"] == 1
    assert not bool(result.loc[0, "stable_as_event"])


def test_stable_filter_marks_missing_support_instead_of_treating_it_as_zero() -> None:
    frame = pd.DataFrame(
        {
            "FDR": [0.001],
            "IncLevelDifference": [0.4],
            "IncLevel1": ["0.8,0.7,0.6"],
            "IncLevel2": ["0.2,0.3,0.4"],
        }
    )
    result = apply_stable_as_filter(
        frame,
        fdr=frame["FDR"],
        standardized_dpsi=frame["IncLevelDifference"],
        reverse_splicing=False,
        thresholds=_thresholds(),
    )
    assert result.loc[0, "stable_support_status"] == "read_support_unavailable"
    assert math.isnan(float(result.loc[0, "mean_junction_reads_group1"]))
    assert not bool(result.loc[0, "stable_as_event"])


def test_reverse_splicing_swaps_groups_but_preserves_standardized_direction() -> None:
    frame = pd.DataFrame(
        {
            "FDR": [0.001],
            "IncLevel1": ["0.2,0.3,0.4"],
            "IncLevel2": ["0.8,0.7,0.6"],
            "IJC_SAMPLE_1": ["20,20,20"],
            "SJC_SAMPLE_1": ["10,10,10"],
            "IJC_SAMPLE_2": ["20,20,20"],
            "SJC_SAMPLE_2": ["10,10,10"],
        }
    )
    standardized = pd.Series([0.4])
    result = apply_stable_as_filter(
        frame,
        fdr=frame["FDR"],
        standardized_dpsi=standardized,
        reverse_splicing=True,
        thresholds=_thresholds(),
    )
    assert result.loc[0, "replicate_direction_consistency"] == 1.0
    assert bool(result.loc[0, "stable_as_event"])
