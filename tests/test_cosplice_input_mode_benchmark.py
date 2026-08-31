from __future__ import annotations

from pathlib import Path

import pandas as pd

from katana_cosplice.input_mode_benchmark import (
    MatrixResult,
    _contrast_metrics,
    _matrix_qc_summary,
)


def _result(tmp_path: Path, label: str) -> MatrixResult:
    samples = [
        "A_rep_1",
        "A_rep_2",
        "A_rep_3",
        "B_rep_1",
        "B_rep_2",
        "B_rep_3",
    ]
    matrix = pd.DataFrame(
        {
            samples[0]: [0.8, 0.2],
            samples[1]: [0.82, 0.22],
            samples[2]: [0.78, 0.18],
            samples[3]: [0.2, 0.3],
            samples[4]: [0.22, 0.32],
            samples[5]: [0.18, 0.28],
        },
        index=["E1", "E2"],
    )
    qc = pd.DataFrame(
        {
            "event_key": ["E1", "E2"],
            "stable_as": [True, True],
            "replicate_repeatability": [0.95, 0.9],
            "condition_mean_psi_range": [0.6, 0.1],
            "mean_total_junction_reads": [100, 80],
        }
    )
    metadata = pd.DataFrame(
        {
            "sample_id": samples,
            "condition": ["A", "A", "A", "B", "B", "B"],
        }
    )
    return MatrixResult(label, tmp_path, matrix, qc, metadata)


def test_matrix_qc_keeps_event_unit_explicit(tmp_path: Path) -> None:
    table = _matrix_qc_summary(_result(tmp_path, "rMATS events"))

    assert table.loc[0, "event_unit"] == "compound canonical AS event"
    assert table.loc[0, "n_stable_events"] == 2


def test_condition_effect_loo_is_stable_for_consistent_replicates(tmp_path: Path) -> None:
    result = _result(tmp_path, "test")
    from katana_cosplice import input_mode_benchmark as benchmark

    original = benchmark.CONTRASTS
    benchmark.CONTRASTS = (("A", "B"),)
    try:
        summary, effects = _contrast_metrics(result, 0.1)
    finally:
        benchmark.CONTRASTS = original

    assert summary.loc[0, "median_LOO_effect_correlation"] > 0.99
    assert summary.loc[0, "median_LOO_sign_agreement"] == 1.0
    assert effects["passes_effect_cutoff"].any()
