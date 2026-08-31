from __future__ import annotations

from pathlib import Path

import pandas as pd

from katana_cosplice.hybrid_benchmark import (
    build_hybrid_benchmark,
    write_hybrid_benchmark,
)


def test_hybrid_benchmark_reports_parallel_branches_and_increment() -> None:
    assignments = pd.DataFrame(
        {
            "event_key": [f"event_{index}" for index in range(8)],
            "module_katana_sparse": [
                "M01",
                "M01",
                "M01",
                "M02",
                "M02",
                "M02",
                "grey",
                "grey",
            ],
            "module_wgcna_exact": [
                "M01",
                "M01",
                "grey",
                "M02",
                "M02",
                "grey",
                "M03",
                "grey",
            ],
        }
    )
    summary = pd.DataFrame(
        [
            {
                "method": "katana_sparse",
                "mean_bootstrap_ARI": 0.7,
                "min_bootstrap_ARI": 0.5,
                "silhouette_cosine": 0.4,
                "within_between_separation": 0.6,
                "runtime_seconds": 1.0,
            },
            {
                "method": "wgcna_exact",
                "mean_bootstrap_ARI": 0.8,
                "min_bootstrap_ARI": 0.6,
                "silhouette_cosine": 0.5,
                "within_between_separation": 0.7,
                "runtime_seconds": 2.0,
            },
        ]
    )

    result = build_hybrid_benchmark(assignments, summary)

    scopes = result.event_assignments["program_scope"].value_counts().to_dict()
    assert scopes == {
        "both": 4,
        "tight_only": 2,
        "diffuse_only": 1,
        "unassigned": 1,
    }
    assert set(result.branch_summary["branch"]) == {"tight", "diffuse"}
    assert result.tight_incremental_programs["n_tight_only_events"].sum() == 2
    equivalence = result.architecture_summary.loc[
        result.architecture_summary["metric"].eq("diffuse_exact_wgcna_equivalence"),
        "value",
    ].iloc[0]
    assert equivalence == 1.0


def test_hybrid_benchmark_writes_tables_and_plot(tmp_path: Path) -> None:
    assignments = pd.DataFrame(
        {
            "event_key": [f"event_{index}" for index in range(6)],
            "module_katana_sparse": ["M01", "M01", "M01", "M02", "M02", "M02"],
            "module_wgcna_exact": ["M01", "M01", "M01", "M02", "M02", "M02"],
        }
    )
    summary = pd.DataFrame(
        [
            {"method": "katana_sparse", "mean_bootstrap_ARI": 0.8},
            {"method": "wgcna_exact", "mean_bootstrap_ARI": 0.8},
        ]
    )
    result = build_hybrid_benchmark(assignments, summary)

    files = write_hybrid_benchmark(result, tmp_path)

    assert (tmp_path / "hybrid_architecture_summary.tsv").exists()
    assert (tmp_path / "hybrid_tight_incremental_programs.tsv").exists()
    assert (tmp_path / "plots" / "hybrid_architecture_audit.png").stat().st_size > 0
    assert files
