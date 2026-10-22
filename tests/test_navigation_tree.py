from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication

from app.widgets.navigation_tree import WorkflowNavigation


def _app() -> QApplication:
    return QApplication.instance() or QApplication([])


def test_workflow_navigation_builds_numbered_hierarchy() -> None:
    _app()
    navigation = WorkflowNavigation()
    for index, name in enumerate(
        [
            "1. Project",
            "5. Per-comparison Analysis",
            "5.1 Overview",
            "5.1.1 Landscape",
            "5.2 Transcript + Splicing Integration",
            "7. RBP",
            "7.1 Per-comparison RBP support",
        ]
    ):
        navigation.add_route(index, name)

    analysis = navigation.topLevelItem(1)
    assert analysis.text(0) == "5  Per-comparison Analysis"
    assert analysis.child(0).text(0) == "5.1  Overview"
    assert analysis.child(0).child(0).text(0) == "5.1.1  Landscape"

    rbp = navigation.topLevelItem(2)
    assert rbp.text(0) == "7  RBP"
    assert rbp.child(0).text(0) == "7.1  Per-comparison RBP support"


def test_workflow_navigation_preserves_legacy_route_index_api() -> None:
    app = _app()
    navigation = WorkflowNavigation()
    names = [
        "1. Project",
        "5. Per-comparison Analysis",
        "5.1 Overview",
        "5.1.1 Landscape",
    ]
    for index, name in enumerate(names):
        navigation.add_route(index, name)

    emitted: list[int] = []
    navigation.currentRowChanged.connect(emitted.append)
    navigation.setCurrentRow(3)
    app.processEvents()

    assert navigation.currentRow() == 3
    assert navigation.currentItem().text(0) == "5.1.1  Landscape"
    assert navigation.currentItem().parent().isExpanded()
    assert emitted[-1] == 3
