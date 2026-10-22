from __future__ import annotations

from collections.abc import Iterable

from katana_next.workflows.base import Workflow
from katana_next.workflows.joint_response import JointResponseWorkflow
from katana_next.workflows.expression_response import ExpressionResponseWorkflow
from katana_next.workflows.module_explanation import ModuleExplanationWorkflow
from katana_next.workflows.narrative_report import NarrativeReportWorkflow
from katana_next.workflows.perturbation import PerturbationWorkflow
from katana_next.workflows.joint_pathways import JointPathwayWorkflow
from katana_next.workflows.robustness_summary import RobustnessSummaryWorkflow
from katana_next.workflows.bulk_interaction import ExploratoryInteractionWorkflow
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
from katana_next.workflows.psi_network import PsiNetworkWorkflow
from katana_next.workflows.djexpress import DJExpressWorkflow
from katana_next.workflows.psi_projection import PsiProjectionWorkflow
from katana_next.workflows.module_functions import ModuleFunctionsWorkflow
from katana_next.workflows.hub_events import HubEventsWorkflow
from katana_next.workflows.single_cell import (
    SingleCellAuditWorkflow,
    SingleCellDifferentialWorkflow,
    SingleCellQuantifyWorkflow,
)
from katana_next.workflows.single_cell_extended import (
    SingleCellBenchmarkWorkflow,
    SingleCellLandscapeWorkflow,
    SingleCellProgramsWorkflow,
)
from katana_next.workflows.validation import NativeSashimiWorkflow
from katana_next.workflows.selected_gene_views import SelectedGeneViewsWorkflow
from katana_next.workflows.motif import EventMotifWorkflow
from katana_next.workflows.dexseq_interaction import DexseqInteractionWorkflow
from katana_next.workflows.isoform_consequences import IsoformConsequencesWorkflow
from katana_next.workflows.single_cell_sashimi import SingleCellSashimiWorkflow
from katana_next.workflows.reliability import ModuleDiagnosticsWorkflow, ReadCalibrationWorkflow
from katana_next.workflows.external_validation import ExternalPreservationWorkflow
from katana_next.workflows.binding import BindingEvidenceWorkflow
from katana_next.workflows.single_cell_isoforms import SingleCellIsoformWorkflow
from katana_next.workflows.single_cell_imports import SingleCellResultImportWorkflow
from katana_next.workflows.coverage_bias import CoverageBiasWorkflow
from katana_next.workflows.cell_embedding import CellEmbeddingWorkflow
from katana_next.workflows.binding_enrichment import BindingEnrichmentWorkflow
from katana_next.workflows.rbp_mechanism import RbpMechanismWorkflow
from katana_next.workflows.se_rbp_impact import SeRbpImpactWorkflow
from katana_next.workflows.rbp_sequence_results import CompletedRbpSequenceWorkflow
from katana_next.workflows.motiftrace_consequences import MotifTraceConsequencesWorkflow
from katana_next.workflows.motiftrace_usage import MotifTraceUsageWorkflow
from katana_next.workflows.motiftrace_parallel import MotifTraceParallelWorkflow


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
        CoverageBiasWorkflow(),
        ReadCalibrationWorkflow(),
        LandscapeWorkflow(),
        DifferentialAsWorkflow(),
        DegWorkflow(),
        DtuWorkflow(),
        CrossComparisonWorkflow(),
        DexseqInteractionWorkflow(),
        ExpressionResponseWorkflow(),
        ExploratoryInteractionWorkflow(),
        CandidateEvidenceWorkflow(),
        JointResponseWorkflow(),
        JointPathwayWorkflow(),
        GeneSetEnrichmentWorkflow(),
        RbpMechanismWorkflow(),
        SeRbpImpactWorkflow(),
        CompletedRbpSequenceWorkflow(),
        MotifTraceConsequencesWorkflow(),
        MotifTraceUsageWorkflow(),
        MotifTraceParallelWorkflow(),
        RbpSupportWorkflow(),
        BindingEvidenceWorkflow(),
        BindingEnrichmentWorkflow(),
        IsoformSupportWorkflow(),
        IsoformConsequencesWorkflow(),
        SelectedGeneViewsWorkflow(),
        EventMotifWorkflow(),
        NativeSashimiWorkflow(),
        PerturbationWorkflow(),
        RobustnessSummaryWorkflow(),
        DJExpressWorkflow(),
        CoSpliceWorkflow(),
        PsiNetworkWorkflow(),
        PsiProjectionWorkflow(),
        ModuleFunctionsWorkflow(),
        HubEventsWorkflow(),
        ModuleExplanationWorkflow(),
        ModuleDiagnosticsWorkflow(),
        ExternalPreservationWorkflow(),
        SingleCellAuditWorkflow(),
        SingleCellQuantifyWorkflow(),
        SingleCellLandscapeWorkflow(),
        CellEmbeddingWorkflow(),
        SingleCellDifferentialWorkflow(),
        SingleCellIsoformWorkflow(),
        SingleCellResultImportWorkflow('marvel'),
        SingleCellResultImportWorkflow('longread_isoforms'),
        SingleCellProgramsWorkflow(),
        SingleCellSashimiWorkflow(),
        SingleCellBenchmarkWorkflow(),
        ReproducibilityReportWorkflow(),
        NarrativeReportWorkflow(),
    )
