from pathlib import Path

import json
import pandas as pd

from katana_cosplice.locked_panel_report import build_locked_panel_report


def test_locked_panel_report_collects_completed_runs(tmp_path: Path) -> None:
    run = tmp_path / "primary"
    search_dir = run / "01_prelock" / "three_objective_cosplice"
    search_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "method": ["stability_only", "event_level_wgcna"],
            "n_modules": [3, 1],
            "assigned_event_fraction": [0.8, 0.5],
            "largest_module_fraction": [0.4, 0.5],
            "mean_bootstrap_ARI": [0.7, 0.9],
            "median_holdout_Zsummary": [3.0, 1.0],
            "median_holdout_program_score_preservation": [0.7, 0.5],
        }
    ).to_csv(run / "locked_three_arm_summary.tsv", sep="\t", index=False)
    (run / "analysis_plan.json").write_text(
        json.dumps(
            {
                "event_type": "ES",
                "discovery": "A",
                "tuning_cohorts": ["B"],
                "holdout": "C",
            }
        ),
        encoding="utf-8",
    )
    (run / "freeze_manifest.json").write_text(
        json.dumps({"n_frozen_events": 100}), encoding="utf-8"
    )
    pd.DataFrame(
        {"selected": [True], "functional_objective_active": [False]}
    ).to_csv(search_dir / "parameter_search.tsv", sep="\t", index=False)
    paths = build_locked_panel_report(tmp_path)
    assert len(paths) == 4
    summary = pd.read_csv(tmp_path / "LOCKED_ABLATION_PANEL_SUMMARY.tsv", sep="\t")
    assert summary["valid_multi_module_partition"].tolist() == [True, False]
    assert (tmp_path / "LOCKED_ABLATION_PANEL_REPORT.md").exists()


def test_locked_panel_report_excludes_all_grey_ari(tmp_path: Path) -> None:
    run = tmp_path / "all_grey_stress_test"
    run.mkdir(parents=True)
    pd.DataFrame(
        {
            "method": ["stability_only", "event_level_wgcna"],
            "n_modules": [0, 1],
            "assigned_event_fraction": [0.0, 1.0],
            "largest_module_fraction": [0.0, 1.0],
            "mean_bootstrap_ARI": [1.0, 0.8],
            "median_holdout_Zsummary": [float("nan"), 0.0],
            "median_holdout_program_score_preservation": [float("nan"), 0.4],
        }
    ).to_csv(run / "locked_three_arm_summary.tsv", sep="\t", index=False)
    (run / "analysis_plan.json").write_text(
        json.dumps(
            {
                "event_type": "ME",
                "discovery": "A",
                "tuning_cohorts": ["B"],
                "holdout": "C",
            }
        ),
        encoding="utf-8",
    )
    (run / "freeze_manifest.json").write_text(
        json.dumps({"n_frozen_events": 68}), encoding="utf-8"
    )

    build_locked_panel_report(tmp_path)

    summary = pd.read_csv(tmp_path / "LOCKED_ABLATION_PANEL_SUMMARY.tsv", sep="\t")
    assert summary["stability_evaluable"].tolist() == [False, True]
    report = (tmp_path / "LOCKED_ABLATION_PANEL_REPORT.md").read_text(encoding="utf-8")
    assert "bootstrap ARI 0/0" in report
    assert "Zero-module/all-grey partitions are excluded" in report
