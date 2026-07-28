from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path
from typing import Any, Protocol, runtime_checkable


class StageStatus(StrEnum):
    NOT_RUN = "not_run"
    RUNNING = "running"
    CACHE_HIT = "cache_hit"
    FINISHED = "finished"
    FAILED = "failed"


class ArtifactKind(StrEnum):
    TABLE = "table"
    MATRIX = "matrix"
    PLOT = "plot"
    PARAMETERS = "parameters"
    REPORT = "report"
    LOG = "log"
    OTHER = "other"


@dataclass(frozen=True)
class ArtifactSpec:
    key: str
    relative_path: str
    kind: ArtifactKind
    required: bool = True
    description: str = ""

    def resolve(self, output_root: Path) -> Path:
        return output_root / self.relative_path


@dataclass(frozen=True)
class StageSpec:
    key: str
    title: str
    order: int
    description: str
    dependencies: tuple[str, ...] = ()
    optional_dependencies: tuple[str, ...] = ()
    artifacts: tuple[ArtifactSpec, ...] = ()
    cli_command: str | None = None
    accepts_force: bool = True
    visible: bool = True
    implemented: bool = True
    explicit_only: bool = True


@dataclass(frozen=True)
class AnnotationProviderSpec:
    key: str
    title: str
    scope: str
    required_inputs: tuple[str, ...]
    description: str


@dataclass
class AnnotationResult:
    provider: str
    available: bool
    message: str
    tables: dict[str, Any] = field(default_factory=dict)
    artifacts: list[Path] = field(default_factory=list)


@runtime_checkable
class AnnotationProvider(Protocol):
    """Extension contract for mechanism annotation without coupling it to networking."""

    spec: AnnotationProviderSpec

    def availability(self, context: Any) -> tuple[bool, str]: ...

    def annotate(self, context: Any) -> AnnotationResult: ...
