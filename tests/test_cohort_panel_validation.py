from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from katana_cosplice.cohort_panel_validation import (
    evaluate_cohort_panel,
    write_cohort_panel_validation,
)


def test_cohort_panel_separates_primary_and_stress_projections(
    tmp_path: Path,
) -> None:
    benchmark = tmp_path / "benchmark"
    benchmark.mkdir()
    pd.DataFrame(
        [
            {
                "cohort": cohort,
                "method": method,
                "min_bootstrap_ARI": ari,
                "median_bootstrap_ARI": ari,
                "weighted_mean_within_module_coassignment": ari,
                "n_modules": 2,
                "grey_fraction": 0.1,
                "valid_partition": True,
            }
            for cohort in ("A", "B", "C")
            for method, ari in (("CoSplice", 0.8), ("Exact WGCNA", 0.7))
        ]
    ).to_csv(benchmark / "cohort_internal_benchmark.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "discovery_cohort": left,
                "external_cohort": right,
                "method": method,
                "best_module_Zsummary": z,
                "median_Zsummary": z,
                "median_density_ratio": 1.0,
                "median_connectivity_correlation": 0.8,
                "median_membership_profile_correlation": 0.8,
                "median_program_score_preservation": score,
                "mapped_events": 50,
            }
            for left, right in (("A", "B"), ("A", "C"))
            for method, z, score in (
                ("CoSplice", 3.0, 0.8),
                ("Exact WGCNA", 1.0, 0.6),
            )
        ]
    ).to_csv(
        benchmark / "directed_projection_preservation_summary.tsv",
        sep="\t",
        index=False,
    )
    config = tmp_path / "config.json"
    config.write_text(
        json.dumps(
            {
                "cohorts": [
                    {"cohort": "A", "family": "organ_1"},
                    {"cohort": "B", "family": "organ_1"},
                    {"cohort": "C", "family": "other"},
                ]
            }
        ),
        encoding="utf-8",
    )
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "cohorts": [
                    {"cohort": cohort, "n_events": 100, "n_samples": 30}
                    for cohort in ("A", "B", "C")
                ]
            }
        ),
        encoding="utf-8",
    )

    result = evaluate_cohort_panel(benchmark, config, manifest)

    scopes = result.preservation_comparison.set_index(
        ["discovery_cohort", "external_cohort"]
    )["validation_scope"]
    assert scopes.loc[("A", "B")] == "organ_matched_primary"
    assert scopes.loc[("A", "C")] == "cross_organ_stress"
    assert result.claim_audit.loc[
        result.claim_audit["claim"].str.contains("preservation"),
        "status",
    ].iloc[0] == "not_supported"
    primary = result.scope_summary.loc[
        result.scope_summary["validation_scope"].eq("organ_matched_primary")
    ].iloc[0]
    assert primary["paired_projection_count"] == 1

    generated = write_cohort_panel_validation(result, tmp_path / "summary")
    generated_names = {path.name for path in generated}
    assert "TEN_COHORT_VALIDATION_REPORT.md" in generated_names
    assert "ten_cohort_validation_summary.png" in generated_names
