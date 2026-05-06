from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QHeaderView,
    QCheckBox,
    QComboBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QAbstractScrollArea,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from src.services.project_service import ProjectService


class VisualizationGroupsPage(QWidget):
    def __init__(self, project_service: ProjectService, on_confirm=None) -> None:
        super().__init__()
        self.project_service = project_service
        self.on_confirm = on_confirm
        self._is_refreshing = False
        self._pending_hscroll = 0
        self._pending_vscroll = 0

        self.label = QLabel(
            "Create custom comparison groups for Cards and Sashimi. "
            "Each row can contain multiple comparisons, for example 1/2/3/4 as one group and 2/3 as another."
        )
        self.label.setWordWrap(True)
        self.add_button = QPushButton("Add Group")
        self.add_button.clicked.connect(self._add_group)
        self.add_button.setMaximumWidth(140)
        self.add_button.setStyleSheet(
            "QPushButton { background-color: #2563EB; color: white; border-radius: 6px; padding: 6px 12px; font-weight: 600; }"
        )
        self.confirm_button = QPushButton("Confirm")
        self.confirm_button.clicked.connect(self._confirm_groups)
        self.confirm_button.setMaximumWidth(120)
        self.confirm_button.setStyleSheet(
            "QPushButton { background-color: #2563EB; color: white; border-radius: 6px; padding: 6px 12px; font-weight: 600; }"
        )

        self.table = QTableWidget(0, 8)
        self.table.setHorizontalHeaderLabels(
            [
                "Use",
                "Display Name",
                "Comparison 1",
                "Comparison 2",
                "Comparison 3",
                "Comparison 4",
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
        self.table.verticalHeader().setVisible(False)
        header = self.table.horizontalHeader()
        header.setStretchLastSection(False)
        header.setSectionResizeMode(QHeaderView.ResizeMode.Fixed)

        header = QHBoxLayout()
        header.addWidget(self.confirm_button)
        header.addWidget(self.add_button)
        header.addStretch(1)

        layout = QVBoxLayout(self)
        layout.addWidget(self.label)
        layout.addLayout(header)
        layout.addWidget(self.table, 1)

    def refresh(self) -> None:
        project = self.project_service.current_project
        self._is_refreshing = True
        self._pending_hscroll = self.table.horizontalScrollBar().value()
        self._pending_vscroll = self.table.verticalScrollBar().value()
        if project is None:
            self.table.setRowCount(0)
            self.add_button.setEnabled(False)
            self.confirm_button.setEnabled(False)
            self._is_refreshing = False
            return

        if not project.visualization_groups:
            self.project_service.add_visualization_group()
            project = self.project_service.current_project

        selected = {
            comparison.comparison_id: comparison.display_resolved_name
            for comparison in self.project_service.selected_comparisons_for_display()
        }
        by_label = {label: comparison_id for comparison_id, label in selected.items()}
        self.add_button.setEnabled(True)
        self.confirm_button.setText("Confirmed" if project.visualization_groups_confirmed else "Confirm")
        self.confirm_button.setEnabled(project.comparison_sets_confirmed and not project.visualization_groups_confirmed)
        self.table.clearContents()
        self.table.setRowCount(len(project.visualization_groups))

        for row, group in enumerate(project.visualization_groups):
            enabled_box = QCheckBox()
            enabled_box.setChecked(group.enabled)
            enabled_box.stateChanged.connect(
                lambda state, group_id=group.group_id: self._toggle_group(
                    group_id, state == Qt.CheckState.Checked.value
                )
            )
            self.table.setCellWidget(row, 0, enabled_box)

            name_item = QTableWidgetItem(group.resolved_name)
            name_item.setData(Qt.ItemDataRole.UserRole, (group.group_id, "display_name"))
            self.table.setItem(row, 1, name_item)

            resolved_ids = [
                self._resolve_comparison_choice(comparison_id, selected, by_label)
                for comparison_id in group.comparison_ids[:4]
            ]
            missing = [
                comparison_id
                for comparison_id in group.comparison_ids[:4]
                if self._resolve_comparison_choice(comparison_id, selected, by_label) is None
            ]
            ordered_ids = [comparison_id for comparison_id in resolved_ids if comparison_id]

            for slot in range(4):
                combo = QComboBox()
                combo.addItem("", None)
                for comparison_id, label in selected.items():
                    combo.addItem(label, comparison_id)
                current_id = resolved_ids[slot] if slot < len(resolved_ids) else None
                index = combo.findData(current_id)
                combo.setCurrentIndex(index if index >= 0 else 0)
                combo.currentIndexChanged.connect(
                    lambda _index, group_id=group.group_id, row_index=row: self._update_group_slots(
                        group_id, row_index
                    )
                )
                self.table.setCellWidget(row, 2 + slot, combo)

            if not ordered_ids:
                status_text = "Empty"
                tooltip = "Select up to four comparisons for this visual group."
            elif missing:
                status_text = "Warning"
                tooltip = "Could not map these stored comparisons:\n" + "\n".join(missing)
            else:
                status_text = f"Ready ({len(ordered_ids)})"
                tooltip = "Included comparisons in order:\n" + "\n".join(
                    selected.get(comparison_id, comparison_id) for comparison_id in ordered_ids
                )
            status_item = QTableWidgetItem(status_text)
            status_item.setToolTip(tooltip)
            status_item.setFlags(status_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.table.setItem(row, 6, status_item)

            remove_button = QPushButton("Remove")
            remove_button.clicked.connect(lambda _checked=False, group_id=group.group_id: self._remove_group(group_id))
            self.table.setCellWidget(row, 7, remove_button)

        self.table.setColumnWidth(0, 44)
        self.table.setColumnWidth(1, 240)
        self.table.setColumnWidth(2, 220)
        self.table.setColumnWidth(3, 220)
        self.table.setColumnWidth(4, 220)
        self.table.setColumnWidth(5, 220)
        self.table.setColumnWidth(6, 140)
        self.table.setColumnWidth(7, 90)
        self._is_refreshing = False
        self.table.horizontalScrollBar().setValue(self._pending_hscroll)
        self.table.verticalScrollBar().setValue(self._pending_vscroll)

    def _handle_item_changed(self, item: QTableWidgetItem) -> None:
        if self._is_refreshing:
            return
        payload = item.data(Qt.ItemDataRole.UserRole)
        if not payload:
            return
        group_id, field = payload
        if field == "display_name":
            self.project_service.update_visualization_group(
                group_id, display_name=item.text().strip() or None
            )
        self.refresh()

    def _toggle_group(self, group_id: str, enabled: bool) -> None:
        if self._is_refreshing:
            return
        self.project_service.update_visualization_group(group_id, enabled=enabled)
        self.refresh()

    def _add_group(self) -> None:
        self.project_service.add_visualization_group()
        self.refresh()

    def _remove_group(self, group_id: str) -> None:
        self.project_service.remove_visualization_group(group_id)
        self.refresh()

    def _update_group_slots(self, group_id: str, row_index: int) -> None:
        if self._is_refreshing:
            return
        ordered_ids: list[str] = []
        for slot in range(4):
            combo = self.table.cellWidget(row_index, 2 + slot)
            if isinstance(combo, QComboBox):
                comparison_id = combo.currentData()
                if comparison_id:
                    ordered_ids.append(str(comparison_id))
        self.project_service.update_visualization_group(group_id, comparison_ids=ordered_ids)
        self.refresh()

    def _confirm_groups(self) -> None:
        self.project_service.confirm_visualization_groups()
        self.refresh()
        if callable(self.on_confirm):
            self.on_confirm()

    @staticmethod
    def _resolve_comparison_choice(
        value: str,
        by_id: dict[str, str],
        by_label: dict[str, str],
    ) -> str | None:
        if value in by_id:
            return value
        return by_label.get(value)
