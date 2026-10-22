from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication, QComboBox

from app.pages.comparison_page import ComparisonPage
from src.domain.models.project import ComparisonDefinition, ProjectConfig
from src.services.project_service import ProjectService


def _application() -> QApplication:
    return QApplication.instance() or QApplication([])


def test_pairing_shows_read_only_canonical_directions_without_direction_selectors(tmp_path: Path) -> None:
    _application()
    comparison = ComparisonDefinition(
        comparison_id="Treatment_vs_Control",
        display_name="Treatment effect",
        rmats_name="Treatment_vs_Control",
        deg_name="Treatment_vs_Control",
        experiment_group="Treatment",
        control_group="Control",
        reverse_splicing=True,
        reverse_deg=False,
        detected_rmats_modes=["JC"],
    )
    service = ProjectService()
    service.current_project = ProjectConfig(
        project_name="test",
        project_root=tmp_path,
        available_comparisons=[comparison],
        selected_comparison_ids=[comparison.comparison_id],
        confirmed=True,
    )

    page = ComparisonPage(service)
    page.refresh()

    headers = [page.table.horizontalHeaderItem(index).text() for index in range(page.table.columnCount())]
    assert headers == [
        "Use",
        "Display Name",
        "rMATS Source",
        "DEG Source",
        "Canonical AS Direction",
        "Canonical DEG Direction",
        "rMATS File",
        "DEG File",
        "Open rMATS",
        "Open DEG",
        "rMATS",
        "Remove",
    ]
    assert page.table.item(0, 4).text() == "Control vs Treatment"
    assert page.table.item(0, 5).text() == "Treatment vs Control"
    assert not page.findChildren(QComboBox)
    assert "read-only" in page.label.text()

    page.deleteLater()
