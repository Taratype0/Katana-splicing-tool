from __future__ import annotations

import json

import pandas as pd

from src.services.artifact_registry import AnalysisArtifactRegistry


def test_registry_publishes_deduplicated_comparison_selection(tmp_path) -> None:
    registry = AnalysisArtifactRegistry(tmp_path / "katana_output", "demo")
    frame = pd.DataFrame(
        [
            {"comparison_id": "A_vs_B", "gene_symbol": "Gene1", "rank": 1},
            {"comparison_id": "A_vs_B", "gene_symbol": "Gene1", "rank": 2},
            {"comparison_id": "A_vs_B", "gene_symbol": "Gene2", "rank": 3},
        ]
    )

    artifact = registry.publish_selection(
        module_id="5.3.2",
        artifact_key="candidate_gene_selection",
        entity_type="gene",
        frame=frame,
        scope="comparison",
        scope_id="A_vs_B",
        comparison_ids=["A_vs_B"],
        id_columns=["comparison_id", "gene_symbol"],
        consumer_modules=["5.5", "7.1", "8.1"],
        selection_rule="test selection",
    )

    loaded, loaded_artifact = registry.latest_selection(
        entity_type="gene",
        scope_id="A_vs_B",
        module_id="5.3.2",
        artifact_key="candidate_gene_selection",
    )
    assert loaded_artifact == artifact
    assert loaded["gene_symbol"].tolist() == ["Gene1", "Gene2"]
    assert artifact.row_count == 2
    assert artifact.comparison_ids == ("A_vs_B",)
    assert registry.resolve_path(artifact).exists()

    payload = json.loads(registry.path.read_text(encoding="utf-8"))
    assert payload["schema_version"] == 1
    assert artifact.artifact_id in payload["artifacts"]


def test_registry_keeps_selection_types_separate(tmp_path) -> None:
    registry = AnalysisArtifactRegistry(tmp_path, "demo")
    for key, gene in (
        ("candidate_gene_selection", "Selected"),
        ("candidate_gene_blacklist", "Blocked"),
    ):
        registry.publish_selection(
            module_id="5.3.2",
            artifact_key=key,
            entity_type="gene",
            frame=pd.DataFrame({"gene_symbol": [gene]}),
            scope="comparison",
            scope_id="A_vs_B",
            id_columns=["gene_symbol"],
        )

    selected, artifact = registry.latest_selection(
        entity_type="gene",
        scope_id="A_vs_B",
        module_id="5.3.2",
        artifact_key="candidate_gene_selection",
    )
    assert artifact is not None
    assert artifact.artifact_key == "candidate_gene_selection"
    assert selected["gene_symbol"].tolist() == ["Selected"]
