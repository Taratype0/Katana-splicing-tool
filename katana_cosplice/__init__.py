"""Sample-level co-splicing network analysis for Katana projects."""

__version__ = "1.2.2"

from katana_cosplice.annotations import AnnotationRegistry
from katana_cosplice.architecture import (
    AnnotationProvider,
    AnnotationProviderSpec,
    AnnotationResult,
    ArtifactKind,
    ArtifactSpec,
    StageSpec,
    StageStatus,
)
from katana_cosplice.config import CoSpliceConfig, load_config
from katana_cosplice.pipeline import CoSplicePipeline
from katana_cosplice.project import CoSpliceProject
from katana_cosplice.stage_registry import stage_spec, stage_specs

__all__ = [
    "AnnotationProvider",
    "AnnotationProviderSpec",
    "AnnotationRegistry",
    "AnnotationResult",
    "ArtifactKind",
    "ArtifactSpec",
    "CoSpliceConfig",
    "CoSplicePipeline",
    "CoSpliceProject",
    "StageSpec",
    "StageStatus",
    "load_config",
    "stage_spec",
    "stage_specs",
]
