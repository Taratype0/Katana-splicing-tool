from __future__ import annotations

from PySide6.QtGui import QFontDatabase
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QComboBox,
    QHBoxLayout,
    QLabel,
    QPlainTextEdit,
    QPushButton,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from src.services.project_service import ProjectService


class InputCheckPage(QWidget):
    def __init__(self, project_service: ProjectService) -> None:
        super().__init__()
        self.project_service = project_service

        self.label = QLabel(
            "Review the loaded files for the currently selected comparisons. "
            "Use the selector to inspect one comparison at a time or all selected comparisons."
        )
        self.comparison_tabs = QTabWidget()
        self.comparison_tabs.currentChanged.connect(self._render_report)
        self.comparison_selector = QComboBox()
        self.comparison_selector.currentIndexChanged.connect(self._render_report)
        self.check_button = QPushButton("Check")
        self.check_button.setMaximumWidth(120)
        self.check_button.setStyleSheet(
            "QPushButton { background-color: #2563EB; color: white; border-radius: 6px; padding: 6px 12px; font-weight: 600; }"
            "QPushButton:disabled { background-color: #64748B; color: #E2E8F0; }"
        )
        self.check_button.clicked.connect(self._render_report)

        self.report = QPlainTextEdit()
        self.report.setReadOnly(True)
        self.report.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        self.report.setFont(QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont))
        self.report.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.report.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)

        layout = QVBoxLayout(self)
        layout.addWidget(self.label)
        layout.addWidget(self.comparison_tabs)
        toolbar = QHBoxLayout()
        toolbar.addWidget(self.comparison_selector, 1)
        toolbar.addWidget(self.check_button)
        layout.addLayout(toolbar)
        layout.addWidget(self.report, 1)

    def refresh(self) -> None:
        project = self.project_service.current_project
        self.comparison_selector.blockSignals(True)
        self.comparison_selector.clear()
        if project is None:
            self.check_button.setEnabled(False)
            self.report.setPlainText("No project loaded.")
            self.comparison_tabs.clear()
            self.comparison_selector.blockSignals(False)
            return

        selected = self.project_service.selected_comparisons_for_display()
        current_tab = self._selected_comparison_id()
        self.comparison_tabs.blockSignals(True)
        self.comparison_tabs.clear()
        self.comparison_tabs.addTab(QWidget(), "All selected comparisons")
        self.comparison_tabs.setTabToolTip(0, "__all__")
        for index, comparison in enumerate(selected, start=1):
            self.comparison_tabs.addTab(QWidget(), comparison.resolved_name)
            self.comparison_tabs.setTabToolTip(index, comparison.comparison_id)
        if current_tab:
            for index in range(self.comparison_tabs.count()):
                if self.comparison_tabs.tabToolTip(index) == current_tab:
                    self.comparison_tabs.setCurrentIndex(index)
                    break
        self.comparison_tabs.blockSignals(False)
        self.comparison_selector.addItem("All selected comparisons", "__all__")
        for comparison in selected:
            self.comparison_selector.addItem(comparison.resolved_name, comparison.comparison_id)
        self.check_button.setEnabled(True)
        self.comparison_selector.blockSignals(False)
        self._render_report()

    def _render_report(self) -> None:
        project = self.project_service.current_project
        if project is None:
            self.report.setPlainText("No project loaded.")
            return
        selected_id = self._selected_comparison_id()
        comparison_id = None if selected_id in (None, "__all__") else str(selected_id)
        report = self.project_service.build_selected_input_confirmation_report(comparison_id)
        self.report.setPlainText(report)
        self.report.verticalScrollBar().setValue(0)

    def _selected_comparison_id(self):
        if self.comparison_tabs.count() > 0:
            tooltip = self.comparison_tabs.tabToolTip(self.comparison_tabs.currentIndex())
            if tooltip:
                return tooltip
        return self.comparison_selector.currentData()
