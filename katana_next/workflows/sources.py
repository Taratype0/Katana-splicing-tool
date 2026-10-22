from __future__ import annotations

import json
from functools import lru_cache
from typing import Any

import pandas as pd

from katana_next.domain.models import AnalysisArtifact
from katana_next.services.project_service import ProjectContext


def source_key(analysis_id: str, name: str | None) -> str:
    return f"{analysis_id}::{name or '*'}"


def selected_source(context: ProjectContext, parameters: dict[str, Any]) -> AnalysisArtifact | None:
    selection_id = str(parameters.get("selection_id", "")).strip()
    if not selection_id:
        return None
    selection = next((s for s in context.catalog.list_selections(context.manifest.project_id)
                      if s.selection_id == selection_id), None)
    if selection is None:
        raise ValueError("The saved selection no longer exists.")
    source = context.catalog.get_artifact(selection.source_artifact_id)
    validate_source(context, source)
    return source


def validate_source(context: ProjectContext, artifact: AnalysisArtifact | None) -> None:
    if artifact is None or artifact.project_id != context.manifest.project_id or not artifact.path.is_file():
        raise ValueError("The selected source artifact is unavailable in this project.")
    if artifact.metadata.get("project_snapshot_fingerprint") != context.snapshot.manifest_fingerprint:
        raise ValueError("The selected source belongs to an older project snapshot. Rebuild explicitly before using it.")
    if artifact.analysis_id in {"quality.stable_as", "differential.as", "landscape.bulk"}:
        from katana_next.science.bulk import READ_SUPPORT_POLICY_VERSION

        if (artifact.metadata.get("read_support_policy") != READ_SUPPORT_POLICY_VERSION
                or artifact.metadata.get("workflow_version") != _workflow_versions().get(artifact.analysis_id)):
            raise ValueError("This saved result uses a historical read-support policy or method version. It is preserved, but cannot stand in for a new every-sample result; run the required module explicitly.")
    if not source_complete(context, artifact):
        raise ValueError("The selected source run has not completed. Finish or rerun that module explicitly first.")


def source_complete(context: ProjectContext, artifact: AnalysisArtifact) -> bool:
    if artifact.analysis_id in {"quality.stable_as", "differential.as", "landscape.bulk"}:
        from katana_next.science.bulk import READ_SUPPORT_POLICY_VERSION

        current = _workflow_versions().get(artifact.analysis_id)
        if (artifact.metadata.get("read_support_policy") != READ_SUPPORT_POLICY_VERSION
                or current is None or artifact.metadata.get("workflow_version") != current):
            return False
    if artifact.metadata.get("artifact_contract_version") != "3.0-correctness":
        return True
    current = _workflow_versions().get(artifact.analysis_id)
    if current is not None and artifact.metadata.get("workflow_version") != current:
        return False
    path = context.results.artifact_directory(artifact.analysis_id, artifact.run_fingerprint) / "parameters.json"
    try:
        marker = json.loads(path.read_text(encoding="utf-8"))
        return marker.get("completed") is True and artifact.artifact_id in {a["artifact_id"] for a in marker["generated_artifacts"]}
    except (OSError, ValueError, KeyError, TypeError):
        return False


@lru_cache(maxsize=1)
def _workflow_versions() -> dict[str, str]:
    from katana_next.workflows.registry import WorkflowRegistry
    registry = WorkflowRegistry()
    return {item["workflow_id"]: registry.get(item["workflow_id"]).version for item in registry.definitions()}


def resolve_table(context: ProjectContext, analysis_id: str, name: str | None = None,
                  parameters: dict[str, Any] | None = None) -> AnalysisArtifact | None:
    parameters = parameters or {}
    bound = parameters.get("source_artifact_ids", {}).get(source_key(analysis_id, name))
    artifacts = ([context.catalog.get_artifact(bound)] if bound else
                 context.catalog.list_artifacts(context.manifest.project_id, analysis_id=analysis_id, kind="table"))
    requested = set(parameters.get("contrast_ids") or [c.contrast_id for c in context.manifest.contrasts])
    contrast_scoped = analysis_id.startswith("differential.") or analysis_id in {"prioritize.evidence", "quality.stable_as"}
    for artifact in artifacts:
        if artifact is None or artifact.analysis_id != analysis_id or (name and artifact.name != name):
            continue
        if artifact.metadata.get("project_snapshot_fingerprint") != context.snapshot.manifest_fingerprint or not artifact.path.is_file():
            continue
        if not source_complete(context, artifact):
            continue
        available = set(artifact.contrast_ids or artifact.parameters.get("contrast_ids") or
                        [c.contrast_id for c in context.manifest.contrasts])
        if contrast_scoped and not requested.issubset(available):
            continue
        return artifact
    if bound:
        raise ValueError(f"Bound source for {name or analysis_id} is stale or does not cover the requested comparisons.")
    return None


def selection_rows(context: ProjectContext, parameters: dict[str, Any], source: AnalysisArtifact,
                   frame: pd.DataFrame) -> pd.DataFrame:
    selection_id = str(parameters.get("selection_id", "")).strip()
    if not selection_id:
        return frame
    selection = next(s for s in context.catalog.list_selections(context.manifest.project_id)
                     if s.selection_id == selection_id)
    if selection.source_artifact_id != source.artifact_id:
        raise ValueError("A saved selection must be resolved against its original source, not the latest table.")
    ids = set(selection.entity_ids)
    selected = frame[frame["row_id"].astype(str).isin(ids)]
    missing = ids - set(selected["row_id"].astype(str))
    if missing:
        raise ValueError(f"{len(missing)} selected rows are absent from their source artifact.")
    return selected.copy()
