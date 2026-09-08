from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class WorkspaceDefinition:
    workspace_id: str
    label: str
    purpose: str
    icon: str
    enabled: bool
    phase: int

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


WORKSPACES = (
    WorkspaceDefinition(
        "project", "Project", "Audit inputs, samples and contrasts.", "folder", True, 1
    ),
    WorkspaceDefinition(
        "quality",
        "Quality",
        "Assess coverage, replicates and exclusions.",
        "shield",
        False,
        2,
    ),
    WorkspaceDefinition(
        "landscape",
        "Landscape",
        "Explore expression and splicing structure.",
        "activity",
        False,
        2,
    ),
    WorkspaceDefinition(
        "differential",
        "Differential",
        "Test DEG, DTU and alternative splicing.",
        "split",
        False,
        2,
    ),
    WorkspaceDefinition(
        "prioritize",
        "Prioritize",
        "Combine transparent evidence dimensions.",
        "target",
        False,
        2,
    ),
    WorkspaceDefinition(
        "validate",
        "Validate",
        "Inspect junction support with Native Sashimi.",
        "scan",
        False,
        2,
    ),
    WorkspaceDefinition(
        "programs",
        "Programs",
        "Discover and validate co-splicing programs.",
        "network",
        False,
        3,
    ),
    WorkspaceDefinition(
        "single-cell",
        "Single Cell",
        "Analyze assay-aware single-cell splicing.",
        "cells",
        False,
        4,
    ),
    WorkspaceDefinition(
        "reports",
        "Reports",
        "Build reproducible evidence packages.",
        "report",
        False,
        5,
    ),
)
