from __future__ import annotations

from pathlib import Path

import pandas as pd

from src.domain.models.project import ComparisonDefinition, ProjectConfig
from src.services.project_service import ProjectService


def _service(tmp_path: Path) -> ProjectService:
    deg_path = tmp_path / "Treatment_vs_Control.deseq2.results.tsv"
    if not deg_path.exists():
        pd.DataFrame(
            {
                "gene_id": ["ENSG1.1", "ENSG2.1", "ENSG3.1"],
                "baseMean": [100.0, 200.0, 50.0],
                "log2FoldChange": [2.0, -1.5, 0.1],
                "padj": [0.001, 0.01, 0.9],
            }
        ).to_csv(deg_path, sep="\t", index=False)
    comparison = ComparisonDefinition(
        comparison_id="comparison-1",
        display_name="Treatment effect",
        deg_name="Treatment_vs_Control",
        deg_path=deg_path,
        has_deg=True,
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
    return service


def test_deg_explorer_runs_independently_and_reloads_cache(tmp_path: Path) -> None:
    service = _service(tmp_path)
    result = service.preview_deg_explorer("comparison-1", allow_generate=True, force_rebuild=True)
    assert len(result) == 3
    assert int(result["DEG_significant"].sum()) == 2
    assert service.module_state("deg_explorer")["status"] == "finished"
    assert all(
        state["status"] == "not_run"
        for name, state in service._module_states.items()
        if name != "deg_explorer"
    )
    files = service.deg_explorer_file_map("comparison-1")
    assert all(files[key].exists() for key in ("all", "significant", "up", "down", "summary", "parameters", "volcano", "ma"))

    second = _service(tmp_path)
    cached = second.preview_deg_explorer("comparison-1", allow_generate=False)
    assert len(cached) == 3
    assert second.module_state("deg_explorer")["status"] == "cache_hit"


def test_functional_enrichment_cache_load_updates_status_without_running(tmp_path: Path) -> None:
    service = _service(tmp_path)
    files = service.enrichment_file_map("comparison-1")
    files["results"].parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {
            "comparison_id": ["comparison-1"],
            "source": ["GO:BP"],
            "term_id": ["GO:0000001"],
            "term_name": ["example process"],
            "FDR": [0.01],
        }
    ).to_csv(files["results"], sep="\t", index=False)

    result = service.load_cached_functional_enrichment("comparison-1")

    assert len(result) == 1
    assert service.module_state("functional_enrichment")["status"] == "cache_hit"
    assert all(
        state["status"] == "not_run"
        for name, state in service._module_states.items()
        if name != "functional_enrichment"
    )
