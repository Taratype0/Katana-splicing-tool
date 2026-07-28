from __future__ import annotations

from pathlib import Path

import yaml

from src.domain.models.project import ComparisonDefinition, ProjectConfig
from src.io.scanner import ProjectScanner
from src.services.project_service import ProjectService


def test_saved_comparison_id_replaces_scanner_alias(tmp_path: Path) -> None:
    config_path = tmp_path / ".katana_project.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "comparisons": [
                    {
                        "comparison_id": "stable-comparison-id",
                        "display_name": "Treatment effect",
                        "rmats_name": "Treatment_vs_Control",
                        "experiment_group": "Treatment",
                        "control_group": "Control",
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    discovered = ComparisonDefinition(
        comparison_id="Treatment_vs_Control",
        rmats_name="Treatment_vs_Control",
        rmats_path=tmp_path / "rmats",
        experiment_group="Treatment",
        control_group="Control",
    )
    project = ProjectConfig(
        project_name="test",
        project_root=tmp_path,
        config_path=config_path,
        available_comparisons=[discovered],
    )

    scanner = ProjectScanner()
    scanner._apply_saved_config(project)
    scanner._merge_alias_comparisons(project)

    assert [item.comparison_id for item in project.available_comparisons] == ["stable-comparison-id"]
    assert project.available_comparisons[0].rmats_path == tmp_path / "rmats"


def test_display_name_is_not_reversed_by_analysis_direction() -> None:
    comparison = ComparisonDefinition(
        comparison_id="stable-id",
        rmats_name="Treatment_vs_Control",
        deg_name="Treatment_vs_Control",
        experiment_group="Treatment",
        control_group="Control",
        reverse_splicing=True,
        reverse_deg=False,
    )
    assert comparison.display_resolved_name == "Treatment vs Control"
    assert comparison.final_as_groups == ("Control", "Treatment")
    assert comparison.final_deg_groups == ("Treatment", "Control")


def test_direction_payload_keeps_deg_and_splicing_directions_independent(tmp_path: Path) -> None:
    comparison = ComparisonDefinition(
        comparison_id="stable-id",
        display_name="Treatment effect",
        experiment_group="Treatment",
        control_group="Control",
        reverse_splicing=True,
        reverse_deg=False,
    )
    service = ProjectService()
    service.current_project = ProjectConfig(
        project_name="test",
        project_root=tmp_path,
        available_comparisons=[comparison],
    )

    payload = service._comparison_parameters_payload(comparison)

    assert payload["log2fc_direction_label"] == (
        "log2FC > 0 means transcript expression is higher in Treatment than Control"
    )
    assert payload["dpsi_direction_label"] == "dPSI > 0 means event usage is higher in Control than Treatment"
    assert payload["display_name"] == "Treatment effect"
