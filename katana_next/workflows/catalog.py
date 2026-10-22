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
        "project", "Samples & Design", "Check input files, sample groups and comparisons.", "folder", True, 1
    ),
    WorkspaceDefinition(
        "quality",
        "Quality Control",
        "Check read coverage and replicate consistency.",
        "shield",
        True,
        2,
    ),
    WorkspaceDefinition(
        "landscape",
        "Splicing Overview",
        "View splice-event distributions, PCA and sample clustering.",
        "activity",
        True,
        2,
    ),
    WorkspaceDefinition(
        "differential",
        "Differential Analysis",
        "Analyze DAS, DEG, transcript usage and factor interactions.",
        "split",
        True,
        2,
    ),
    WorkspaceDefinition(
        "prioritize",
        "Functional Analysis",
        "Analyze GO/KEGG, RBP motifs, evidence overlap and isoforms.",
        "target",
        True,
        2,
    ),
    WorkspaceDefinition(
        "validate",
        "Read Validation",
        "View BAM coverage, Sashimi plots and target-gene RNA evidence.",
        "scan",
        True,
        2,
    ),
    WorkspaceDefinition(
        "programs",
        "Splicing Networks",
        "Build cosplice modules and examine their functions and stability.",
        "network",
        True,
        3,
    ),
    WorkspaceDefinition(
        "single-cell",
        "Single-cell Splicing",
        "Analyze imported single-cell splicing data by cell type and donor.",
        "cells",
        True,
        4,
    ),
    WorkspaceDefinition(
        "reports",
        "Reports & Export",
        "Generate illustrated study reports and export full results.",
        "report",
        True,
        5,
    ),
)
