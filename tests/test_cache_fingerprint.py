from __future__ import annotations

from pathlib import Path

from src.domain.models.project import ComparisonDefinition, ProjectConfig
from src.services.project_service import ProjectService


def _service_with_project(tmp_path: Path) -> tuple[ProjectService, ComparisonDefinition]:
    comparison = ComparisonDefinition(
        comparison_id="comparison-1",
        rmats_name="Treatment_vs_Control",
        deg_name="Treatment_vs_Control",
        experiment_group="Treatment",
        control_group="Control",
    )
    service = ProjectService()
    service.current_project = ProjectConfig(
        project_name="test",
        project_root=tmp_path,
        output_root=tmp_path / "katana_output",
        available_comparisons=[comparison],
        selected_comparison_ids=[comparison.comparison_id],
        comparison_order_ids=[comparison.comparison_id],
    )
    return service, comparison


def test_candidate_signature_changes_with_thresholds_and_direction(tmp_path: Path) -> None:
    service, comparison = _service_with_project(tmp_path)
    initial = service._candidate_cache_signature()

    service.current_thresholds.splicing_dpsi = 0.2
    threshold_changed = service._candidate_cache_signature()
    assert threshold_changed != initial

    service.current_thresholds.splicing_dpsi = 0.1
    comparison.reverse_splicing = True
    direction_changed = service._candidate_cache_signature()
    assert direction_changed != initial


def test_candidate_cache_without_fingerprint_is_not_loaded(tmp_path: Path) -> None:
    service, _ = _service_with_project(tmp_path)
    legacy = tmp_path / "katana_output" / "04a_candidate_gene_screening" / "data"
    legacy.mkdir(parents=True)
    (legacy / "gene_level_integrated_candidates.tsv").write_text(
        "comparison_id\tgene_symbol\ncomparison-1\tGeneA\n",
        encoding="utf-8",
    )
    assert service._load_cached_candidate_gene_table().empty
