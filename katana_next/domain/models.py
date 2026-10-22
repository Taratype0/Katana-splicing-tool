from __future__ import annotations

from datetime import datetime, timezone
from enum import StrEnum
from pathlib import Path
from typing import Any
from uuid import uuid4

from pydantic import BaseModel, ConfigDict, Field, model_validator


PLATFORM_SCHEMA_VERSION = "2.0"


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


class FrozenModel(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")


class InputAsset(FrozenModel):
    asset_id: str
    kind: str
    path: Path
    required: bool = True
    metadata: dict[str, Any] = Field(default_factory=dict)


class Sample(FrozenModel):
    sample_id: str
    biological_sample_id: str
    condition: str
    batch: str | None = None
    attributes: dict[str, str] = Field(default_factory=dict)
    assets: dict[str, Path] = Field(default_factory=dict)


class ContrastArm(FrozenModel):
    label: str
    sample_ids: tuple[str, ...]

    @model_validator(mode="after")
    def require_samples(self) -> "ContrastArm":
        if not self.sample_ids:
            raise ValueError("A contrast arm must contain at least one sample.")
        if len(self.sample_ids) != len(set(self.sample_ids)):
            raise ValueError("A contrast arm cannot contain duplicate sample IDs.")
        return self


class Contrast(FrozenModel):
    contrast_id: str
    label: str
    numerator: ContrastArm
    denominator: ContrastArm

    @model_validator(mode="after")
    def require_disjoint_arms(self) -> "Contrast":
        overlap = set(self.numerator.sample_ids) & set(self.denominator.sample_ids)
        if overlap:
            raise ValueError(f"Contrast arms overlap: {sorted(overlap)}")
        return self


class AnnotationReference(FrozenModel):
    genome_build: str
    gtf: Path | None = None
    genome_fasta: Path | None = None
    transcript_fasta: Path | None = None


class ProjectManifest(FrozenModel):
    schema_version: str = PLATFORM_SCHEMA_VERSION
    project_id: str = Field(default_factory=lambda: str(uuid4()))
    name: str
    input_root: Path
    output_root: Path
    annotation: AnnotationReference
    assets: tuple[InputAsset, ...] = ()
    samples: tuple[Sample, ...] = ()
    contrasts: tuple[Contrast, ...] = ()
    metadata: dict[str, Any] = Field(default_factory=dict)

    @model_validator(mode="after")
    def validate_references(self) -> "ProjectManifest":
        if self.schema_version != PLATFORM_SCHEMA_VERSION:
            raise ValueError(
                f"Unsupported manifest schema {self.schema_version!r}; "
                f"expected {PLATFORM_SCHEMA_VERSION!r}."
            )
        sample_ids = [sample.sample_id for sample in self.samples]
        if len(sample_ids) != len(set(sample_ids)):
            raise ValueError("Sample IDs must be unique.")
        known = set(sample_ids)
        contrast_ids = [contrast.contrast_id for contrast in self.contrasts]
        if len(contrast_ids) != len(set(contrast_ids)):
            raise ValueError("Contrast IDs must be unique; display labels do not define identity.")
        for contrast in self.contrasts:
            referenced = set(contrast.numerator.sample_ids) | set(
                contrast.denominator.sample_ids
            )
            missing = referenced - known
            if missing:
                raise ValueError(
                    f"Contrast {contrast.contrast_id!r} references unknown samples: "
                    f"{sorted(missing)}"
                )
        asset_ids = [asset.asset_id for asset in self.assets]
        if len(asset_ids) != len(set(asset_ids)):
            raise ValueError("Input asset IDs must be unique.")
        return self


class FileFingerprint(FrozenModel):
    path: Path
    exists: bool
    size: int | None = None
    modified_ns: int | None = None
    digest: str | None = None
    hash_method: str = "missing"


class ProjectSnapshot(FrozenModel):
    project_id: str
    manifest_path: Path
    name: str
    input_root: Path
    output_root: Path
    genome_build: str
    manifest_fingerprint: str
    input_fingerprints: tuple[FileFingerprint, ...]
    sample_count: int
    contrast_count: int
    opened_at: datetime = Field(default_factory=utc_now)
    warnings: tuple[str, ...] = ()


class ArtifactStatus(StrEnum):
    PENDING = "pending"
    RUNNING = "running"
    COMPLETE = "complete"
    FAILED = "failed"
    CANCELLED = "cancelled"


class ArtifactKind(StrEnum):
    TABLE = "table"
    PLOT_SPEC = "plot_spec"
    SELECTION = "selection"
    REPORT = "report"
    FILE = "file"


class AnalysisArtifact(FrozenModel):
    artifact_id: str = Field(default_factory=lambda: str(uuid4()))
    project_id: str
    analysis_id: str
    name: str
    kind: ArtifactKind
    entity_type: str
    status: ArtifactStatus = ArtifactStatus.COMPLETE
    fingerprint: str
    run_fingerprint: str
    path: Path
    schema_version: str = PLATFORM_SCHEMA_VERSION
    row_count: int | None = None
    parent_artifact_ids: tuple[str, ...] = ()
    contrast_ids: tuple[str, ...] = ()
    parameters: dict[str, Any] = Field(default_factory=dict)
    metadata: dict[str, Any] = Field(default_factory=dict)
    warnings: tuple[str, ...] = ()
    created_at: datetime = Field(default_factory=utc_now)


class SelectionSet(FrozenModel):
    selection_id: str = Field(default_factory=lambda: str(uuid4()))
    project_id: str
    name: str
    entity_type: str
    source_artifact_id: str
    entity_ids: tuple[str, ...]
    filters: dict[str, Any] = Field(default_factory=dict)
    created_at: datetime = Field(default_factory=utc_now)


class JobStatus(StrEnum):
    QUEUED = "queued"
    RUNNING = "running"
    COMPLETE = "complete"
    FAILED = "failed"
    CANCELLED = "cancelled"
    CACHE_HIT = "cache_hit"


class JobEvent(FrozenModel):
    sequence: int
    timestamp: datetime = Field(default_factory=utc_now)
    progress: float = Field(ge=0.0, le=1.0)
    message: str


class JobRecord(FrozenModel):
    job_id: str = Field(default_factory=lambda: str(uuid4()))
    project_id: str
    analysis_id: str
    status: JobStatus
    run_fingerprint: str
    parameters: dict[str, Any] = Field(default_factory=dict)
    progress: float = Field(default=0.0, ge=0.0, le=1.0)
    message: str = ""
    cache_status: str = "rebuild"
    artifact_ids: tuple[str, ...] = ()
    events: tuple[JobEvent, ...] = ()
    created_at: datetime = Field(default_factory=utc_now)
    started_at: datetime | None = None
    completed_at: datetime | None = None
    error: str | None = None
