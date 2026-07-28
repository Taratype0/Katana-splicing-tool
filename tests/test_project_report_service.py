from __future__ import annotations

import json

import pandas as pd

from src.services.artifact_registry import AnalysisArtifactRegistry
from src.services.project_report_service import ProjectReportService


def test_project_report_uses_cached_registry_without_running_analysis(tmp_path) -> None:
    registry = AnalysisArtifactRegistry(tmp_path, "Report project")
    registry.publish_selection(
        module_id="5.3.2",
        artifact_key="candidate_gene_selection",
        entity_type="gene",
        frame=pd.DataFrame(
            {
                "comparison_id": ["cmp-a"],
                "gene_symbol": ["GeneA"],
            }
        ),
        scope="comparison",
        scope_id="cmp-a",
        comparison_ids=["cmp-a"],
        id_columns=["comparison_id", "gene_symbol"],
    )

    outputs = ProjectReportService(tmp_path, "Report project").build(
        registry=registry,
        project_summary={"project_name": "Report project", "selected_comparisons": 1},
        module_states={"candidate_selection": {"status": "cache_hit", "message": "Loaded cache."}},
    )

    assert all(path.exists() for path in outputs.values())
    inventory = pd.read_csv(outputs["inventory"], sep="\t")
    assert inventory["entity_type"].tolist() == ["gene"]
    provenance = json.loads(outputs["provenance"].read_text(encoding="utf-8"))
    assert provenance["artifact_count"] == 1
    assert provenance["module_states"]["candidate_selection"]["status"] == "cache_hit"
    assert "Katana reproducible analysis report" in outputs["report"].read_text(encoding="utf-8")
