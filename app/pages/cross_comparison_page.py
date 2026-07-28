from __future__ import annotations

from app.pages.cross_comparison_controller import (
    CrossComparisonController,
    CrossComparisonPageContract,
)
from app.pages.results_page import ResultsPage
from src.services.project_service import ProjectService


class CrossComparisonResultsPage(ResultsPage):
    """Dedicated 5.4 page boundary backed by the shared results renderer.

    The class intentionally keeps the existing renderers while making the 5.4
    route contract explicit. Individual renderers can now move out of
    ``ResultsPage`` incrementally without changing main-window navigation.
    """

    def __init__(
        self,
        project_service: ProjectService,
        on_run_modules=None,
        *,
        subbranch: str,
    ) -> None:
        self.page_contract: CrossComparisonPageContract = CrossComparisonController.contract(subbranch)
        super().__init__(
            project_service,
            on_run_modules,
            fixed_branch="group",
            fixed_subbranch=self.page_contract.subbranch,
        )
