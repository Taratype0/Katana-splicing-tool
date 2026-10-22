from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Protocol

import pandas as pd
from katana_cosplice.branding import display_module_name

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
        payload = asdict(self)
        for key in ("label", "help"):
            payload[key] = display_module_name(payload[key])
        return payload


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
    supports_contrast_selection: bool = False
    enabled: bool = True
    limitations: str = ""
    presentation_role: str = "primary"

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        for key in ("label", "question", "method", "output_summary", "limitations"):
            payload[key] = display_module_name(payload[key])
        payload["parameters"] = [parameter.to_dict() for parameter in self.parameters]
        from katana_next.workflows.presentation import apply_workflow_copy
        payload = apply_workflow_copy(payload)
        from katana_next.science.das_selection import workflow_das_capability
        capability = workflow_das_capability(self.workflow_id)
        payload["das_selection"] = capability
        if capability["mode"] == "switchable" and not any(p["key"] == "event_set" for p in payload["parameters"]):
            payload["parameters"].append(ParameterDefinition(
                "event_set", "DAS event set", "select", "stable_significant",
                "RS-DAS is the default: significant plus IJC + SJC >= 10 in every comparison sample. DAS uses original significance without the read-support foreground gate. Both original flags and caller statistics are preserved.",
                ("stable_significant", "significant"),
            ).to_dict())
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
        from katana_next.science.das_selection import normalize_workflow_selection
        values = normalize_workflow_selection(self.definition.workflow_id, values)
        for definition in self.definition.parameters:
            value = values.get(definition.key)
            if value is None:
                continue
            try:
                if definition.kind == "number":
                    if isinstance(definition.default, int) and not isinstance(
                        definition.default, bool
                    ):
                        numeric = float(value)
                        if not numeric.is_integer():
                            raise ValueError("must be an integer")
                        values[definition.key] = int(numeric)
                    else:
                        values[definition.key] = float(value)
                elif definition.kind == "boolean":
                    if isinstance(value, str):
                        normalized = value.strip().casefold()
                        if normalized not in {"true", "false"}:
                            raise ValueError("must be true or false")
                        values[definition.key] = normalized == "true"
                    else:
                        values[definition.key] = bool(value)
                if definition.kind == "number":
                    import math
                    number = values[definition.key]
                    if not math.isfinite(number):
                        raise ValueError("must be finite")
                    if definition.minimum is not None and number < definition.minimum:
                        raise ValueError(f"must be at least {definition.minimum}")
                    if definition.maximum is not None and number > definition.maximum:
                        raise ValueError(f"must be at most {definition.maximum}")
                if definition.kind == "select" and definition.options and values[definition.key] not in definition.options:
                    raise ValueError("must be one of the declared choices")
            except (TypeError, ValueError) as error:
                raise ValueError(
                    f"Invalid value for {definition.label}: {value!r} ({error})."
                ) from error
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
        if not len(frame.columns):
            if len(frame):
                raise ValueError("A nonempty artifact table requires named columns")
            frame = pd.DataFrame({"row_id": pd.Series(dtype="string")})
        path = self.directory / filename
        row_count = self.context.results.write_dataframe(frame, path)
        tsv = path.with_suffix(".tsv")
        frame.to_csv(tsv, sep="\t", index=False)
        metadata = {**(metadata or {}), "exports": {"tsv": str(tsv)}}
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
        from katana_next.workflows.column_details import COLUMN_DETAILS
        artifact = artifact.model_copy(update={"metadata": {**artifact.metadata, "column_definitions": {
            column: COLUMN_DETAILS[column] for column in frame if column in COLUMN_DETAILS
        }}})
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
        from katana_next.workflows.figure_exports import export_vega
        exports = export_vega(self.context.results.read_plot_spec(path), self.directory / "plots", Path(filename).stem)
        metadata = {**(metadata or {}), "exports": exports, "static_export_font": "Arial"}
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
                "completed": True,
                "workflow_id": self.analysis_id,
                "workflow_version": self.workflow_version,
                "run_fingerprint": self.run_fingerprint,
                "parameters": self.parameters,
                "generated_plot_files": [path for artifact in self.artifacts for kind, path in artifact.metadata.get("exports", {}).items() if kind in {"png", "svg"}],
                "generated_table_files": [path for artifact in self.artifacts for kind, path in artifact.metadata.get("exports", {}).items() if kind == "tsv"],
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
            "artifact_contract_version": "3.0-correctness",
            "project_snapshot_fingerprint": self.context.snapshot.manifest_fingerprint,
            **self.common_metadata,
            **(metadata or {}),
        }
        from katana_next.science.das_selection import workflow_das_capability, selection_metadata
        capability = workflow_das_capability(self.analysis_id)
        merged_metadata["das_selection_capability"] = capability
        if capability["mode"] in {"switchable", "rs_required"}:
            merged_metadata["das_selection"] = selection_metadata(self.parameters.get("event_set"))
            if self.analysis_id == "differential.as":
                merged_metadata["das_selection"]["source_table_scope"] = "all events with both original significant and stable_significant flags"
            elif capability["mode"] == "rs_required":
                merged_metadata["das_selection"]["requirement"] = "RS-required; ordinary-DAS switching is unsupported"
        if self.context.manifest.metadata.get("gene_filter") and not self.analysis_id.startswith("single_cell."):
            from katana_next.science.gene_filter import gene_filter_details, gene_filter_footnote
            merged_metadata["gene_filter"] = gene_filter_details(self.context.manifest)
            merged_metadata["gene_filter_footnote"] = gene_filter_footnote(self.context.manifest)
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
