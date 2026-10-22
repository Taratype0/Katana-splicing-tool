from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from katana_cosplice.publication_benchmark import build_publication_benchmark


def _write_method_result(path: Path, *, cosplice_ari: float, wgcna_ari: float) -> None:
    path.mkdir(parents=True)
    pd.DataFrame(
        [
            {
                "method": "katana_sparse",
                "status": "finished",
                "valid_partition": True,
                "n_events": 100,
                "n_samples": 30,
                "n_modules": 3,
                "grey_fraction": 0.1,
                "mean_bootstrap_ARI": cosplice_ari,
                "runtime_seconds": 2.0,
            },
            {
                "method": "wgcna_exact",
                "status": "finished",
                "valid_partition": True,
                "n_events": 100,
                "n_samples": 30,
                "n_modules": 2,
                "grey_fraction": 0.2,
                "mean_bootstrap_ARI": wgcna_ari,
                "runtime_seconds": 5.0,
            },
        ]
    ).to_csv(path / "method_benchmark_summary.tsv", sep="\t", index=False)


def test_publication_benchmark_builds_traceable_report(tmp_path: Path) -> None:
    result_dir = tmp_path / "dataset"
    _write_method_result(result_dir, cosplice_ari=0.8, wgcna_ari=0.7)
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "benchmark_title": "Test publication benchmark",
                "datasets": [
                    {
                        "dataset_id": "PUBLIC_A",
                        "accession": "GSE00000",
                        "result_dir": str(result_dir),
                        "data_type": "event-by-sample PSI",
                        "design": "test",
                        "validation_role": "same-input comparison",
                        "citation_url": "https://example.org",
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    output = tmp_path / "output"
    result = build_publication_benchmark(manifest, output)

    assert (output / "dataset_method_metrics.tsv").exists()
    assert (output / "paired_cosplice_tight_vs_wgcna.tsv").exists()
    assert (output / "publication_head_to_head_scorecard.png").stat().st_size > 0
    assert (output / "PUBLICATION_HEAD_TO_HEAD_REPORT.md").exists()
    assert (output / "parameters.json").exists()
    contrast = result.paired_contrasts
    row = contrast.loc[contrast["metric"].eq("mean_bootstrap_ARI")].iloc[0]
    assert row["winner"] == "cosplice_tight"
    assert row["interpretable"]


def test_degenerate_partitions_are_not_scored_as_wins(tmp_path: Path) -> None:
    result_dir = tmp_path / "dataset"
    result_dir.mkdir()
    pd.DataFrame(
        [
            {
                "method": "katana_sparse",
                "status": "degenerate",
                "valid_partition": False,
                "n_modules": 0,
                "mean_bootstrap_ARI": 1.0,
            },
            {
                "method": "wgcna_exact",
                "status": "finished",
                "valid_partition": True,
                "n_modules": 2,
                "mean_bootstrap_ARI": 0.8,
            },
        ]
    ).to_csv(result_dir / "method_benchmark_summary.tsv", sep="\t", index=False)
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "datasets": [
                    {
                        "dataset_id": "PUBLIC_B",
                        "result_dir": str(result_dir),
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    result = build_publication_benchmark(manifest, tmp_path / "output")
    row = result.paired_contrasts.loc[
        result.paired_contrasts["metric"].eq("mean_bootstrap_ARI")
    ].iloc[0]
    assert not row["interpretable"]
    assert row["winner"] == "not_interpretable"


def test_string_false_partition_flag_is_not_treated_as_truthy(tmp_path: Path) -> None:
    result_dir = tmp_path / "dataset"
    result_dir.mkdir()
    pd.DataFrame(
        [
            {
                "method": "katana_sparse",
                "status": "degenerate",
                "valid_partition": "False",
                "n_modules": 1,
                "mean_bootstrap_ARI": 1.0,
            },
            {
                "method": "wgcna_exact",
                "status": "finished",
                "valid_partition": "True",
                "n_modules": 2,
                "mean_bootstrap_ARI": 0.5,
            },
        ]
    ).to_csv(result_dir / "method_benchmark_summary.tsv", sep="\t", index=False)
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "datasets": [{"dataset_id": "PUBLIC_C", "result_dir": str(result_dir)}],
            }
        ),
        encoding="utf-8",
    )

    result = build_publication_benchmark(manifest, tmp_path / "output")

    row = result.paired_contrasts.loc[
        result.paired_contrasts["metric"].eq("mean_bootstrap_ARI")
    ].iloc[0]
    assert not row["interpretable"]
    assert row["winner"] == "not_interpretable"
