from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Protocol

import pandas as pd

from katana_next.domain.models import AnalysisArtifact, ArtifactKind
from katana_next.infrastructure.fingerprints import fingerprint_file, fingerprint_value
from katana_next.infrastructure.jobs import JobContext
from katana_next.services.project_service import ProjectContext


@dataclass(frozen=True)
class ParameterDefinition:
    key: str
    label: str
    kind: str
    default: Any
    help: str
    options: tuple[str, ...] = ()
    minimum: float | None = None
    maximum: float | None = None
    step: float | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class WorkflowDefinition:
    workflow_id: str
    workspace_id: str
    label: str
    question: str
    method: str
    output_summary: str
    analysis_unit: str
    expected_artifacts: int = 1
    parameters: tuple[ParameterDefinition, ...] = ()
    required_asset_kinds: tuple[str, ...] = ()
    optional_asset_kinds: tuple[str, ...] = ()
    enabled: bool = True
    limitations: str = ""

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["parameters"] = [parameter.to_dict() for parameter in self.parameters]
        return payload


class Workflow(Protocol):
    definition: WorkflowDefinition
    version: str

    def input_paths(
        self, context: ProjectContext, parameters: dict[str, Any]
    ) -> list[Path]: ...

    def parent_artifacts(
        self, context: ProjectContext, parameters: dict[str, Any]
    ) -> list[AnalysisArtifact]: ...

    def run(
        self,
        context: ProjectContext,
        job: JobContext,
        run_fingerprint: str,
        parameters: dict[str, Any],
    ) -> list[str]: ...


class WorkflowBase:
    definition: WorkflowDefinition
    version = "1.0.0"

    def input_paths(
        self, context: ProjectContext, parameters: dict[str, Any]
    ) -> list[Path]:
        del parameters
        kinds = set(self.definition.required_asset_kinds) | set(
            self.definition.optional_asset_kinds
        )
        paths = [asset.path for asset in context.manifest.assets if asset.kind in kinds]
        return list(dict.fromkeys(path.resolve() for path in paths if path.exists()))

    def parent_artifacts(
        self, context: ProjectContext, parameters: dict[str, Any]
    ) -> list[AnalysisArtifact]:
        del context, parameters
        return []

    def normalize_parameters(self, supplied: dict[str, Any] | None) -> dict[str, Any]:
        values = {
            definition.key: definition.default
            for definition in self.definition.parameters
        }
        values.update(supplied or {})
        return values

    def validate_assets(self, context: ProjectContext) -> None:
        available = {asset.kind for asset in context.manifest.assets if asset.path.exists()}
        missing = [
            kind for kind in self.definition.required_asset_kinds if kind not in available
        ]
        if missing:
            raise ValueError(
                f"{self.definition.label} requires input asset kinds: {', '.join(missing)}"
            )


@dataclass
class ArtifactWriter:
    context: ProjectContext
    analysis_id: str
    run_fingerprint: str
    parameters: dict[str, Any]
    workflow_version: str
    common_metadata: dict[str, Any] = field(default_factory=dict)
    artifacts: list[AnalysisArtifact] = field(default_factory=list)

    @property
    def directory(self) -> Path:
        return self.context.results.artifact_directory(
            self.analysis_id, self.run_fingerprint
        )

    def table(
        self,
        name: str,
        filename: str,
        frame: pd.DataFrame,
        *,
        entity_type: str,
        metadata: dict[str, Any] | None = None,
        parents: tuple[str, ...] = (),
        warnings: tuple[str, ...] = (),
    ) -> AnalysisArtifact:
        path = self.directory / filename
        row_count = self.context.results.write_dataframe(frame, path)
        artifact = self._artifact(
            name,
            ArtifactKind.TABLE,
            entity_type,
            path,
            metadata=metadata,
            parents=parents,
            warnings=warnings,
            row_count=row_count,
        )
        return self._register(artifact)

    def plot(
        self,
        name: str,
        filename: str,
        specification: dict[str, Any],
        *,
        entity_type: str,
        metadata: dict[str, Any] | None = None,
        parents: tuple[str, ...] = (),
        warnings: tuple[str, ...] = (),
    ) -> AnalysisArtifact:
        path = self.context.results.write_plot_spec(specification, self.directory / filename)
        artifact = self._artifact(
            name,
            ArtifactKind.PLOT_SPEC,
            entity_type,
            path,
            metadata=metadata,
            parents=parents,
            warnings=warnings,
        )
        return self._register(artifact)

    def existing_file(
        self,
        name: str,
        path: Path,
        *,
        kind: ArtifactKind = ArtifactKind.FILE,
        entity_type: str,
        metadata: dict[str, Any] | None = None,
        parents: tuple[str, ...] = (),
        warnings: tuple[str, ...] = (),
    ) -> AnalysisArtifact:
        if not path.is_file():
            raise FileNotFoundError(path)
        artifact = self._artifact(
            name,
            kind,
            entity_type,
            path.resolve(),
            metadata=metadata,
            parents=parents,
            warnings=warnings,
        )
        return self._register(artifact)

    def finish(self) -> list[str]:
        self.context.results.write_plot_spec(
            {
                "schema_version": "2.0",
                "workflow_id": self.analysis_id,
                "workflow_version": self.workflow_version,
                "run_fingerprint": self.run_fingerprint,
                "parameters": self.parameters,
                "generated_artifacts": [
                    {
                        "artifact_id": artifact.artifact_id,
                        "name": artifact.name,
                        "kind": artifact.kind.value,
                        "path": str(artifact.path),
                        "fingerprint": artifact.fingerprint,
                    }
                    for artifact in self.artifacts
                ],
            },
            self.directory / "parameters.json",
        )
        self.context.catalog.export_manifest(self.context.manifest.project_id)
        return [artifact.artifact_id for artifact in self.artifacts]

    def _artifact(
        self,
        name: str,
        kind: ArtifactKind,
        entity_type: str,
        path: Path,
        *,
        metadata: dict[str, Any] | None,
        parents: tuple[str, ...],
        warnings: tuple[str, ...],
        row_count: int | None = None,
    ) -> AnalysisArtifact:
        merged_metadata = {
            "workflow_version": self.workflow_version,
            "project_snapshot_fingerprint": self.context.snapshot.manifest_fingerprint,
            **self.common_metadata,
            **(metadata or {}),
        }
        content_fingerprint = fingerprint_file(path, force_full=True)
        return AnalysisArtifact(
            project_id=self.context.manifest.project_id,
            analysis_id=self.analysis_id,
            name=name,
            kind=kind,
            entity_type=entity_type,
            fingerprint=fingerprint_value(
                {
                    "run": self.run_fingerprint,
                    "name": name,
                    "path": path.name,
                    "content": content_fingerprint.model_dump(mode="json"),
                }
            ),
            run_fingerprint=self.run_fingerprint,
            path=path,
            row_count=row_count,
            parent_artifact_ids=parents,
            contrast_ids=tuple(self.parameters.get("contrast_ids", ())),
            parameters=self.parameters,
            metadata=merged_metadata,
            warnings=warnings,
        )

    def _register(self, artifact: AnalysisArtifact) -> AnalysisArtifact:
        self.context.catalog.register_artifact(artifact)
        self.artifacts.append(artifact)
        return artifact
