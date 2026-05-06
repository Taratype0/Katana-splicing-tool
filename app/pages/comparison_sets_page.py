from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QPushButton,
    QAbstractScrollArea,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from src.services.project_service import ProjectService


class ComparisonSetsPage(QWidget):
    def __init__(self, project_service: ProjectService, on_confirm=None) -> None:
        super().__init__()
        self.project_service = project_service
        self.on_confirm = on_confirm
        self.label = QLabel("Define pairwise comparisons for program comparison analysis.")
        self.label.setWordWrap(True)
        self.confirm_button = QPushButton("Confirm")
        self.confirm_button.clicked.connect(self._confirm_comparison_sets)
        self.confirm_button.setMaximumWidth(120)
        self.confirm_button.setStyleSheet(
            "QPushButton { background-color: #2563EB; color: white; border-radius: 6px; padding: 6px 12px; font-weight: 600; }"
            "QPushButton:disabled { background-color: #64748B; color: #E2E8F0; }"
        )
        self.table = QTableWidget(0, 10)
        self.table.setHorizontalHeaderLabels(
            [
                "Use",
                "Display Name",
                "Comparison A",
                "A Meaning",
                "Comparison B",
                "B Meaning",
                "Set Experiment",
                "Set Control",
                "Status",
                "Remove",
            ]
        )
        self.table.itemChanged.connect(self._handle_item_changed)
        self.table.setHorizontalScrollMode(QAbstractItemView.ScrollMode.ScrollPerPixel)
        self.table.setVerticalScrollMode(QAbstractItemView.ScrollMode.ScrollPerPixel)
        self.table.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.table.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.table.setSizeAdjustPolicy(QAbstractScrollArea.SizeAdjustPolicy.AdjustIgnored)
        self.table.setMinimumWidth(0)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Fixed)
        self._is_refreshing = False
        self._pending_hscroll = 0
        self._pending_vscroll = 0

        button_row = QHBoxLayout()
        button_row.addWidget(self.confirm_button)
        self.add_button = QPushButton("Add Pair")
        self.add_button.clicked.connect(self._add_pair)
        button_row.addWidget(self.add_button)
        button_row.addStretch(1)

        layout = QVBoxLayout(self)
        layout.addWidget(self.label)
        layout.addLayout(button_row)
        layout.addWidget(self.table, 1)

    def refresh(self) -> None:
        project = self.project_service.current_project
        self._is_refreshing = True
        self._pending_hscroll = self.table.horizontalScrollBar().value()
        self._pending_vscroll = self.table.verticalScrollBar().value()
        if project is None:
            self.label.setText("Load a project first.")
            self.confirm_button.setEnabled(False)
            self.table.setRowCount(0)
            self._is_refreshing = False
            return

        self.label.setText(
            "Create the pairwise program-comparison definitions you want to run. "
            "Each row compares Comparison A against Comparison B. "
            "A Meaning and B Meaning show what each selected comparison represents after Pairing. "
            "Set Experiment and Set Control are the labels/semantics for this pairwise comparison result, not the raw source file orientation. "
            "Display Name, Set Experiment and Set Control auto-follow your A/B choices unless you overwrite them. "
            "The Pairing names and directions are inherited here automatically. "
            "Status shows whether both sides are ready for DEG, rMATS and isoform-aware analysis."
        )
        self.confirm_button.setText("Confirmed" if project.comparison_sets_confirmed else "Confirm")
        self.confirm_button.setEnabled(project.pairing_confirmed and not project.comparison_sets_confirmed)
        if not project.comparison_pairs:
            self.project_service.add_comparison_pair()
            project = self.project_service.current_project

        comparisons = self.project_service.selected_comparisons_for_display()
        comparison_options = [("", "None")] + [
            (comparison.comparison_id, comparison.display_resolved_name) for comparison in comparisons
        ]

        self.table.setRowCount(len(project.comparison_pairs))
        for row, pair in enumerate(project.comparison_pairs):
            use_box = QCheckBox()
            use_box.setChecked(pair.enabled)
            use_box.stateChanged.connect(
                lambda state, pair_id=pair.pair_id: self._toggle_pair(pair_id, state == Qt.CheckState.Checked.value)
            )
            self.table.setCellWidget(row, 0, use_box)

            name_item = QTableWidgetItem(pair.resolved_name)
            name_item.setData(Qt.ItemDataRole.UserRole, (pair.pair_id, "display_name"))
            self.table.setItem(row, 1, name_item)

            combo_a = QComboBox()
            combo_b = QComboBox()
            for value, label in comparison_options:
                combo_a.addItem(label, value)
                combo_b.addItem(label, value)
            resolved_a = self._resolved_comparison_id(pair.comparison_a)
            resolved_b = self._resolved_comparison_id(pair.comparison_b)
            self._set_combo_value(combo_a, resolved_a or "")
            self._set_combo_value(combo_b, resolved_b or "")
            combo_a.currentIndexChanged.connect(
                lambda _index, pair_id=pair.pair_id, combo=combo_a: self._update_pair_source(
                    pair_id, "comparison_a", combo.currentData()
                )
            )
            combo_b.currentIndexChanged.connect(
                lambda _index, pair_id=pair.pair_id, combo=combo_b: self._update_pair_source(
                    pair_id, "comparison_b", combo.currentData()
                )
            )
            self.table.setCellWidget(row, 2, combo_a)

            left_summary = QTableWidgetItem(self._comparison_meaning(resolved_a or pair.comparison_a))
            left_summary.setFlags(left_summary.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.table.setItem(row, 3, left_summary)

            self.table.setCellWidget(row, 4, combo_b)

            right_summary = QTableWidgetItem(self._comparison_meaning(resolved_b or pair.comparison_b))
            right_summary.setFlags(right_summary.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.table.setItem(row, 5, right_summary)

            experiment_item = QTableWidgetItem(pair.experiment_group or "")
            experiment_item.setData(Qt.ItemDataRole.UserRole, (pair.pair_id, "experiment_group"))
            self.table.setItem(row, 6, experiment_item)

            control_item = QTableWidgetItem(pair.control_group or "")
            control_item.setData(Qt.ItemDataRole.UserRole, (pair.pair_id, "control_group"))
            self.table.setItem(row, 7, control_item)

            status_text, status_tooltip = self._pair_status(pair)
            status_item = QTableWidgetItem(status_text)
            status_item.setFlags(status_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            status_item.setToolTip(status_tooltip)
            self.table.setItem(row, 8, status_item)

            remove_button = QPushButton("Remove")
            remove_button.clicked.connect(lambda _checked=False, pair_id=pair.pair_id: self._remove_pair(pair_id))
            self.table.setCellWidget(row, 9, remove_button)

        self.table.setColumnWidth(0, 44)
        self.table.setColumnWidth(1, 240)
        self.table.setColumnWidth(2, 220)
        self.table.setColumnWidth(3, 220)
        self.table.setColumnWidth(4, 220)
        self.table.setColumnWidth(5, 220)
        self.table.setColumnWidth(6, 180)
        self.table.setColumnWidth(7, 180)
        self.table.setColumnWidth(8, 180)
        self.table.setColumnWidth(9, 90)
        self._is_refreshing = False
        self.table.horizontalScrollBar().setValue(self._pending_hscroll)
        self.table.verticalScrollBar().setValue(self._pending_vscroll)

    def _handle_item_changed(self, item: QTableWidgetItem) -> None:
        if self._is_refreshing:
            return
        payload = item.data(Qt.ItemDataRole.UserRole)
        if not payload:
            return
        pair_id, field = payload
        self.project_service.update_comparison_pair(str(pair_id), **{field: item.text().strip() or None})
        self.refresh()

    def _toggle_pair(self, pair_id: str, enabled: bool) -> None:
        if self._is_refreshing:
            return
        self.project_service.update_comparison_pair(pair_id, enabled=enabled)
        self.refresh()

    def _update_pair_source(self, pair_id: str, field: str, comparison_id: str) -> None:
        if self._is_refreshing:
            return
        self.project_service.update_comparison_pair(pair_id, **{field: comparison_id or None})
        self.refresh()

    def _add_pair(self) -> None:
        self.project_service.add_comparison_pair()
        self.refresh()

    def _remove_pair(self, pair_id: str) -> None:
        self.project_service.remove_comparison_pair(pair_id)
        self.refresh()

    def _confirm_comparison_sets(self) -> None:
        self.project_service.confirm_comparison_sets()
        self.refresh()
        if callable(self.on_confirm):
            self.on_confirm()

    def _pair_status(self, pair) -> tuple[str, str]:
        project = self.project_service.current_project
        if project is None:
            return "No project", "Load a project first."
        left = self._resolve_comparison(pair.comparison_a)
        right = self._resolve_comparison(pair.comparison_b)
        if left is None or right is None:
            missing = []
            if left is None:
                missing.append("Comparison A")
            if right is None:
                missing.append("Comparison B")
            message = f"Missing selection: {', '.join(missing)}"
            return message, message

        issues = []
        if left.rmats_path is None:
            issues.append(f"{left.display_resolved_name}: missing rMATS")
        if right.rmats_path is None:
            issues.append(f"{right.display_resolved_name}: missing rMATS")
        if left.deg_path is None:
            issues.append(f"{left.display_resolved_name}: missing DEG")
        if right.deg_path is None:
            issues.append(f"{right.display_resolved_name}: missing DEG")
        if not left.has_quant:
            issues.append(f"{left.display_resolved_name}: no quant.sf")
        if not right.has_quant:
            issues.append(f"{right.display_resolved_name}: no quant.sf")

        if issues:
            summary = "Warning"
            if all("missing DEG" in issue for issue in issues):
                summary = "Missing DEG"
            elif all("missing rMATS" in issue for issue in issues):
                summary = "Missing rMATS"
            elif all("no quant.sf" in issue for issue in issues):
                summary = "No quant"
            return summary, "\n".join(issues)

        return "Ready", (
            f"{left.display_resolved_name} and {right.display_resolved_name} both have rMATS, DEG and quant.sf support."
        )

    def _comparison_meaning(self, comparison_id: str | None) -> str:
        comparison = self._resolve_comparison(comparison_id)
        if comparison is None:
            return ""
        experiment = comparison.experiment_group or comparison.source_experiment_group or "?"
        control = comparison.control_group or comparison.source_control_group or "?"
        return f"{experiment} vs {control}"

    def _resolved_comparison_id(self, reference: str | None) -> str | None:
        comparison = self._resolve_comparison(reference)
        return comparison.comparison_id if comparison is not None else None

    def _resolve_comparison(self, reference: str | None):
        project = self.project_service.current_project
        if project is None or not reference:
            return None
        reference = str(reference)
        for comparison in project.available_comparisons:
            if reference in {
                comparison.comparison_id,
                comparison.display_resolved_name,
                comparison.resolved_name,
                comparison.display_name or "",
            }:
                return comparison
        return None

    @staticmethod
    def _set_combo_value(combo: QComboBox, value: str) -> None:
        index = combo.findData(value)
        if index >= 0:
            combo.setCurrentIndex(index)
