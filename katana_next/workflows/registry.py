from __future__ import annotations

from collections.abc import Iterable

from katana_next.workflows.base import Workflow
from katana_next.workflows.bulk import (
    CandidateEvidenceWorkflow,
    CrossComparisonWorkflow,
    DegWorkflow,
    DifferentialAsWorkflow,
    DtuWorkflow,
    LandscapeWorkflow,
    ProjectAuditWorkflow,
    StableAsWorkflow,
)
from katana_next.workflows.interpretation import (
    GeneSetEnrichmentWorkflow,
    IsoformSupportWorkflow,
    RbpSupportWorkflow,
)
from katana_next.workflows.programs import CoSpliceWorkflow
from katana_next.workflows.single_cell import (
    SingleCellAuditWorkflow,
    SingleCellDifferentialWorkflow,
    SingleCellQuantifyWorkflow,
)
from katana_next.workflows.validation import NativeSashimiWorkflow


class WorkflowRegistry:
    """Single production inventory for discoverable and runnable workflows."""

    def __init__(self, workflows: Iterable[Workflow] | None = None) -> None:
        registered = tuple(workflows or _production_workflows())
        self._workflows = {
            workflow.definition.workflow_id: workflow for workflow in registered
        }
        if len(self._workflows) != len(registered):
            raise ValueError("Workflow IDs must be unique.")

    def get(self, workflow_id: str) -> Workflow:
        try:
            return self._workflows[workflow_id]
        except KeyError as exc:
            raise KeyError(f"Unknown workflow: {workflow_id}") from exc

    def all(self) -> tuple[Workflow, ...]:
        return tuple(self._workflows.values())

    def definitions(self) -> list[dict[str, object]]:
        return [workflow.definition.to_dict() for workflow in self.all()]


def _production_workflows() -> tuple[Workflow, ...]:
    from katana_next.workflows.reporting import ReproducibilityReportWorkflow

    return (
        ProjectAuditWorkflow(),
        StableAsWorkflow(),
        LandscapeWorkflow(),
        DifferentialAsWorkflow(),
        DegWorkflow(),
        DtuWorkflow(),
        CrossComparisonWorkflow(),
        CandidateEvidenceWorkflow(),
        GeneSetEnrichmentWorkflow(),
        RbpSupportWorkflow(),
        IsoformSupportWorkflow(),
        NativeSashimiWorkflow(),
        CoSpliceWorkflow(),
        SingleCellAuditWorkflow(),
        SingleCellQuantifyWorkflow(),
        SingleCellDifferentialWorkflow(),
        ReproducibilityReportWorkflow(),
    )
