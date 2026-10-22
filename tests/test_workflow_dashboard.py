from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication, QPushButton

from app.pages.workflow_dashboard_page import WorkflowDashboardPage
from src.services.project_service import ProjectService


def test_workflow_dashboard_status_aggregation() -> None:
    assert WorkflowDashboardPage._aggregate_status([]) == ("on_demand", "On demand")
    assert WorkflowDashboardPage._aggregate_status(["not_run"]) == ("not_run", "Not run")
    assert WorkflowDashboardPage._aggregate_status(["finished", "cache_hit"]) == ("finished", "Complete")
    assert WorkflowDashboardPage._aggregate_status(["finished", "not_run"]) == ("partial", "1/2 ready")
    assert WorkflowDashboardPage._aggregate_status(["running", "not_run"]) == ("running", "Running")
    assert WorkflowDashboardPage._aggregate_status(["failed", "finished"]) == ("failed", "Needs attention")


def test_workflow_dashboard_routes_without_running_analysis() -> None:
    _app = QApplication.instance() or QApplication([])
    opened: list[str] = []
    service = ProjectService()
    page = WorkflowDashboardPage(service, opened.append)

    primary_buttons = [
        button for button in page.findChildren(QPushButton)
        if button.text() == "Overview"
    ]
    assert len(primary_buttons) == 1
    primary_buttons[0].click()
    assert opened == ["5.1 Overview"]
    assert service.module_state("splicing_landscape")["status"] == "not_run"
