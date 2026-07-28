from __future__ import annotations

import json

import pandas as pd

from src.analysis.native_splicing_patterns import (
    psi_pattern_clusters,
    replicate_group_separation,
    row_zscore_matrix,
    sample_similarity_tables,
)
from src.services.project_service import ProjectService


def _psi_matrix() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "event_matching_key": [f"event_{index}" for index in range(8)],
            "cmp:A:rep1": [0.1, 0.2, 0.8, 0.7, 0.1, 0.2, 0.8, 0.7],
            "cmp:A:rep2": [0.11, 0.21, 0.79, 0.69, 0.12, 0.22, 0.78, 0.68],
            "cmp:B:rep1": [0.8, 0.7, 0.2, 0.1, 0.8, 0.7, 0.2, 0.1],
            "cmp:B:rep2": [0.79, 0.69, 0.21, 0.11, 0.78, 0.68, 0.22, 0.12],
        }
    )


def test_sample_similarity_keeps_overlap_denominator() -> None:
    long, correlation, overlap = sample_similarity_tables(_psi_matrix())
    assert not long.empty
    row = long.loc[
        (long["sample_A"] == "cmp:A:rep1") & (long["sample_B"] == "cmp:A:rep2")
    ].iloc[0]
    assert row["shared_measured_events"] == 8
    assert row["spearman_PSI_similarity"] > 0.95
    assert correlation.shape == overlap.shape == (4, 5)


def test_replicate_separation_distinguishes_within_and_between_groups() -> None:
    matrix = _psi_matrix()
    similarity, _, _ = sample_similarity_tables(matrix)
    rows = []
    for sample in matrix.columns[1:]:
        group = sample.split(":")[1]
        for event in matrix["event_matching_key"]:
            rows.append(
                {
                    "sample_id": sample,
                    "comparison_id": "cmp",
                    "sample_group": group,
                    "event_matching_key": event,
                }
            )
    summary = replicate_group_separation(pd.DataFrame(rows), similarity)
    assert len(summary) == 1
    assert summary.iloc[0]["replicate_separation_delta"] > 1.5


def test_pattern_clusters_use_sample_psi_and_cap_cluster_count() -> None:
    matrix = _psi_matrix()
    metadata = pd.DataFrame(
        {
            "event_matching_key": matrix["event_matching_key"],
            "gene_symbol": [f"Gene{index}" for index in range(len(matrix))],
            "event_type": ["SE"] * len(matrix),
            "event_id": matrix["event_matching_key"],
        }
    )
    clusters, summary = psi_pattern_clusters(matrix, metadata)
    assert len(clusters) == len(matrix)
    assert 2 <= clusters["psi_pattern_cluster"].nunique() <= 8
    assert "cluster_method" in summary.columns
    assert "sample/replicate PSI trajectories" in summary.iloc[0]["cluster_method"]


def test_row_zscore_matrix_preserves_event_key() -> None:
    result = row_zscore_matrix(_psi_matrix())
    assert result.columns[0] == "event_matching_key"
    assert result["event_matching_key"].tolist() == _psi_matrix()["event_matching_key"].tolist()
    numeric = result.drop(columns=["event_matching_key"])
    assert (numeric.mean(axis=1).abs() < 1e-10).all()


def test_sample_pca_keeps_biological_sample_identity() -> None:
    matrix = _psi_matrix().rename(
        columns={
            "cmp:A:rep1": "A_rep1",
            "cmp:A:rep2": "A_rep2",
            "cmp:B:rep1": "B_rep1",
            "cmp:B:rep2": "B_rep2",
        }
    )
    psi_long = pd.DataFrame(
        [
            {
                "biological_sample_id": sample_id,
                "sample_group": sample_id.split("_")[0],
                "comparison_id": "cmp",
            }
            for sample_id in matrix.columns[1:]
        ]
    )
    result = ProjectService()._native_jutils_sample_pca(matrix, psi_long)
    assert set(result["sample_id"]) == set(matrix.columns[1:])
    assert set(result["sample_group"]) == {"A", "B"}


def test_native_pattern_cache_requires_current_schema(tmp_path) -> None:
    service = ProjectService()
    assert not service._native_jutils_cache_schema_matches(tmp_path)
    (tmp_path / "parameters.json").write_text(
        json.dumps({"analysis_schema_version": "old"}),
        encoding="utf-8",
    )
    assert not service._native_jutils_cache_schema_matches(tmp_path)
    (tmp_path / "parameters.json").write_text(
        json.dumps({"analysis_schema_version": service.NATIVE_PATTERN_SCHEMA_VERSION}),
        encoding="utf-8",
    )
    assert service._native_jutils_cache_schema_matches(tmp_path)


def test_native_pattern_summary_records_denominators_and_low_confidence() -> None:
    tables = {
        "unified_as_event_table": pd.DataFrame({"event": ["a", "b", "c"]}),
        "sample_level_psi_matrix": pd.DataFrame(
            {
                "event_matching_key": ["a", "b"],
                "sample_1": [0.1, 0.2],
                "sample_2": [0.3, 0.4],
            }
        ),
        "psi_pattern_cluster_summary": pd.DataFrame({"cluster": [1, 2]}),
        "pairwise_as_similarity": pd.DataFrame(
            {
                "direction_similarity": [0.8, 0.6],
                "weighted_cosine_similarity": [0.7, 0.5],
                "similarity_confidence": ["high", "low_shared_support"],
            }
        ),
    }
    summary = ProjectService._native_pattern_summary_table(
        tables,
        comparison_ids=("cmp_a", "cmp_b"),
        fdr_cutoff=0.05,
        abs_dpsi_cutoff=0.1,
        top_events=100,
        selection_mode="supervised",
    ).iloc[0]
    assert summary["comparison_count"] == 2
    assert summary["event_count"] == 3
    assert summary["sample_count"] == 2
    assert summary["pair_count"] == 2
    assert summary["low_confidence_pair_count"] == 1
    assert summary["median_direction_similarity"] == 0.7
    assert "denominator" in summary["interpretation"]
