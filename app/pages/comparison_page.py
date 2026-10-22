from __future__ import annotations

import os
import subprocess
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QHBoxLayout,
    QCheckBox,
    QLineEdit,
    QLabel,
    QPushButton,
    QAbstractScrollArea,
    QHeaderView,
    QSizePolicy,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from app.widgets.page_header import PageHeader
from src.services.project_service import ProjectService


class ComparisonPage(QWidget):
    def __init__(self, project_service: ProjectService, on_confirm=None) -> None:
        super().__init__()
        self.project_service = project_service
        self.on_confirm = on_confirm
        self.header = PageHeader(
            "Pairing and canonical direction",
            "Choose which detected comparisons continue downstream, verify source files, and preserve their analysis order.",
            "02 · Experimental design",
        )
        self.label = QLabel("Comparison pairing will appear here.")
        self.label.setWordWrap(True)
        self.label.setObjectName("contextPanel")
        self.confirm_button = QPushButton("Confirm")
        self.confirm_button.setProperty("class", "primary")
        self.confirm_button.clicked.connect(self._confirm_pairing)
        self.confirm_button.setMaximumWidth(120)
        self.move_up_button = QPushButton("Move up")
        self.move_down_button = QPushButton("Move down")
        self.move_up_button.clicked.connect(self._move_selected_up)
        self.move_down_button.clicked.connect(self._move_selected_down)
        self.move_up_button.setEnabled(False)
        self.move_down_button.setEnabled(False)
        self.table = QTableWidget(0, 12)
        self.table.setHorizontalHeaderLabels(
            [
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
        )
        self.table.itemChanged.connect(self._handle_item_changed)
        self.table.itemSelectionChanged.connect(self._update_move_buttons)
        self._is_refreshing = False
        self._pending_hscroll = 0
        self._pending_vscroll = 0
        self._pending_current_id: str | None = None

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(12)
        layout.addWidget(self.header)
        layout.addWidget(self.label)
        button_row = QHBoxLayout()
        button_row.addWidget(self.confirm_button)
        button_row.addWidget(self.move_up_button)
        button_row.addWidget(self.move_down_button)
        button_row.addStretch(1)
        layout.addLayout(button_row)
        layout.addWidget(self.table, 1)
        header = self.table.horizontalHeader()
        header.setStretchLastSection(False)
        header.setSectionResizeMode(QHeaderView.ResizeMode.Fixed)
        header.setMinimumSectionSize(42)
        self.table.verticalHeader().setVisible(False)
        self.table.setWordWrap(False)
        self.table.setHorizontalScrollMode(QAbstractItemView.ScrollMode.ScrollPerPixel)
        self.table.setVerticalScrollMode(QAbstractItemView.ScrollMode.ScrollPerPixel)
        self.table.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.table.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.table.setSizeAdjustPolicy(QAbstractScrollArea.SizeAdjustPolicy.AdjustIgnored)
        self.table.setMinimumWidth(0)
        self.table.viewport().setMinimumWidth(0)
        self.table.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Expanding)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)

    def refresh(self) -> None:
        project = self.project_service.current_project
        self._is_refreshing = True
        self._pending_hscroll = self.table.horizontalScrollBar().value()
        self._pending_vscroll = self.table.verticalScrollBar().value()
        self._pending_current_id = self._current_comparison_id()
        if project is None:
            self.label.setText("Load a project first.")
            self.confirm_button.setEnabled(False)
            self.table.setRowCount(0)
            self._is_refreshing = False
            return

        comparisons = list(self.project_service.available_comparisons_for_display())
        selected_ids = set(project.selected_comparison_ids)
        self.label.setText(
            "Pair rMATS and DEG comparisons here. Canonical AS and DEG directions are read-only and inherited from the "
            "detected comparison/configuration; downstream modules use standardized dPSI and standardized log2FC. "
            "DEXSeq and DTU inherit the canonical AS/rMATS direction. "
            "Confirm saves the selected comparisons and display order. Display names are inherited by Comparison Sets. "
            "Use controls whether a detected comparison flows into Comparison Sets and Analysis, without changing the raw Project scan. "
            "Legacy reverse settings remain readable for old projects but are no longer editable in this page.\n"
            f"Detected comparisons: {len(comparisons)} | Enabled for downstream steps: {len(selected_ids)}"
        )
        self.confirm_button.setText("Confirmed" if project.pairing_confirmed else "Confirm")
        self.confirm_button.setEnabled(project.confirmed and not project.pairing_confirmed)
        self.table.blockSignals(True)
        self.table.clearContents()
        self.table.setRowCount(len(comparisons))

        for row, comparison in enumerate(comparisons):
            use_box = QCheckBox()
            use_box.setChecked(comparison.comparison_id in selected_ids)
            use_box.stateChanged.connect(
                lambda state, comparison_id=comparison.comparison_id: self._toggle_comparison_use(
                    comparison_id, state == Qt.CheckState.Checked.value
                )
            )
            self.table.setCellWidget(row, 0, use_box)

            display_edit = QLineEdit(comparison.display_resolved_name)
            display_edit.setToolTip(
                "Display name is saved by stable comparison_id.\n"
                "Editing this row only updates this comparison, even after refresh or reordering."
            )
            display_edit.editingFinished.connect(
                lambda cid=comparison.comparison_id, editor=display_edit: self._commit_display_name_edit(cid, editor)
            )
            self.table.setCellWidget(row, 1, display_edit)

            rmats_name_item = QTableWidgetItem(comparison.display_source_name("rmats"))
            rmats_name_item.setData(Qt.ItemDataRole.UserRole, (comparison.comparison_id, "rmats_name"))
            rmats_name_item.setToolTip(
                f"Original source name: {comparison.source_name('rmats')}\n"
                f"Resolved display name: {comparison.display_source_name('rmats')}\n"
                f"Matched file: {comparison.rmats_file_label}"
            )
            self.table.setItem(row, 2, rmats_name_item)

            deg_name_item = QTableWidgetItem(comparison.display_source_name("deg"))
            deg_name_item.setData(Qt.ItemDataRole.UserRole, (comparison.comparison_id, "deg_name"))
            deg_name_item.setToolTip(
                f"Original source name: {comparison.source_name('deg')}\n"
                f"Resolved display name: {comparison.display_source_name('deg')}\n"
                f"Matched file: {comparison.deg_file_label}"
            )
            self.table.setItem(row, 3, deg_name_item)

            as_direction_item = QTableWidgetItem(self._format_groups(comparison.final_as_groups))
            as_direction_item.setFlags(as_direction_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            as_direction_item.setToolTip(
                "Read-only canonical splicing direction. standardized_dPSI > 0 means greater "
                "event usage/inclusion on the left side. DEXSeq and DTU inherit this direction."
            )
            self.table.setItem(row, 4, as_direction_item)

            deg_direction_item = QTableWidgetItem(self._format_groups(comparison.final_deg_groups))
            deg_direction_item.setFlags(deg_direction_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            deg_direction_item.setToolTip(
                "Read-only canonical transcript direction. standardized_log2FC > 0 means greater "
                "expression on the left side."
            )
            self.table.setItem(row, 5, deg_direction_item)

            rmats_file_item = QTableWidgetItem(comparison.rmats_file_label)
            rmats_file_item.setToolTip(comparison.rmats_file_label)
            rmats_file_item.setFlags(rmats_file_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.table.setItem(row, 6, rmats_file_item)

            deg_file_item = QTableWidgetItem(comparison.deg_file_label)
            deg_file_item.setToolTip(comparison.deg_file_label)
            deg_file_item.setFlags(deg_file_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.table.setItem(row, 7, deg_file_item)

            open_rmats_button = self._build_open_button(
                comparison.rmats_path,
                label="Open",
                tooltip=comparison.rmats_file_label,
            )
            self.table.setCellWidget(row, 8, open_rmats_button)

            open_deg_button = self._build_open_button(
                comparison.deg_path,
                label="Open",
                tooltip=comparison.deg_file_label,
            )
            self.table.setCellWidget(row, 9, open_deg_button)

            rmats_mode_item = QTableWidgetItem(comparison.detected_rmats_modes_label)
            rmats_mode_item.setFlags(rmats_mode_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.table.setItem(row, 10, rmats_mode_item)

            remove_button = QPushButton("Remove")
            remove_button.clicked.connect(
                lambda _checked=False, comparison_id=comparison.comparison_id: self._remove_comparison(comparison_id)
            )
            self.table.setCellWidget(row, 11, remove_button)

        self.table.setColumnWidth(0, 44)
        self.table.setColumnWidth(1, 200)
        self.table.setColumnWidth(2, 200)
        self.table.setColumnWidth(3, 200)
        self.table.setColumnWidth(4, 220)
        self.table.setColumnWidth(5, 220)
        self.table.setColumnWidth(6, 240)
        self.table.setColumnWidth(7, 240)
        self.table.setColumnWidth(8, 96)
        self.table.setColumnWidth(9, 96)
        self.table.setColumnWidth(10, 90)
        self.table.setColumnWidth(11, 90)
        self.table.blockSignals(False)
        self._is_refreshing = False
        self.table.horizontalScrollBar().setValue(self._pending_hscroll)
        self.table.verticalScrollBar().setValue(self._pending_vscroll)
        self._restore_current_selection()
        self._update_move_buttons()

    def _handle_item_changed(self, item: QTableWidgetItem) -> None:
        if self._is_refreshing:
            return

        payload = item.data(Qt.ItemDataRole.UserRole)
        if not payload:
            return
        comparison_id, field_name = payload
        comparison = self._find_comparison(str(comparison_id))
        if comparison is None:
            return

        new_text = item.text().strip() or None
        if field_name == "display_name":
            self.project_service.invalidate_pairing_confirmation()
            self.project_service.rename_comparison(str(comparison_id), item.text())
            self.refresh()
            return
        if field_name == "rmats_name":
            self.project_service.invalidate_pairing_confirmation()
            stored_text = comparison._swap_comparison_name(new_text or "", comparison.reverse_splicing) if new_text else None
            self.project_service.update_comparison_mapping(str(comparison_id), rmats_name=stored_text)
            self.refresh()
            return
        if field_name == "deg_name":
            self.project_service.invalidate_pairing_confirmation()
            stored_text = comparison._swap_comparison_name(new_text or "", comparison.reverse_deg) if new_text else None
            self.project_service.update_comparison_mapping(str(comparison_id), deg_name=stored_text)
            self.refresh()
            return

        return

    def _commit_display_name_edit(self, comparison_id: str, editor: QLineEdit) -> None:
        if self._is_refreshing:
            return
        comparison = self._find_comparison(comparison_id)
        if comparison is None:
            return
        new_text = editor.text()
        current_text = comparison.display_name or comparison.display_resolved_name
        if new_text == current_text:
            return
        self.project_service.invalidate_pairing_confirmation()
        self.project_service.rename_comparison(comparison_id, new_text)
        self.refresh()

    def _confirm_pairing(self) -> None:
        try:
            self.project_service.confirm_pairing()
        except RuntimeError as exc:
            self.label.setText(f"{self.label.text()}\n{exc}")
            self.refresh()
            return
        self.refresh()
        if callable(self.on_confirm):
            self.on_confirm()

    def _toggle_comparison_use(self, comparison_id: str, enabled: bool) -> None:
        if self._is_refreshing or self.project_service.current_project is None:
            return
        selected = list(self.project_service.current_project.selected_comparison_ids)
        if enabled:
            if comparison_id not in selected:
                selected.append(comparison_id)
        else:
            selected = [item for item in selected if item != comparison_id]
        self.project_service.invalidate_pairing_confirmation()
        self.project_service.set_selected_comparisons(selected)
        self.refresh()

    def _remove_comparison(self, comparison_id: str) -> None:
        self.project_service.remove_comparison_from_selection(comparison_id)
        self.project_service.invalidate_pairing_confirmation()
        self.refresh()

    def _move_selected_up(self) -> None:
        comparison_id = self._current_comparison_id()
        if not comparison_id:
            return
        self.project_service.move_comparison_up(comparison_id)
        self.refresh()

    def _move_selected_down(self) -> None:
        comparison_id = self._current_comparison_id()
        if not comparison_id:
            return
        self.project_service.move_comparison_down(comparison_id)
        self.refresh()

    def _build_open_button(self, path: Path | None, *, label: str, tooltip: str) -> QPushButton:
        button = QPushButton(label)
        button.setMaximumWidth(84)
        button.setToolTip(tooltip)
        if path is None or not path.exists():
            button.setEnabled(False)
            return button
        button.clicked.connect(lambda _checked=False, target=Path(path): self._open_original_path(target))
        return button

    def _open_original_path(self, path: Path) -> None:
        try:
            target = path.expanduser().resolve()
            if not target.exists():
                return
            if target.is_file():
                subprocess.Popen(["explorer.exe", f"/select,{target}"])
                return
            os.startfile(str(target))
        except Exception:
            # Keep the UI responsive even if Explorer fails to open.
            return

    def _find_comparison(self, comparison_id: str):
        project = self.project_service.current_project
        if project is None:
            return None
        return next((item for item in project.available_comparisons if item.comparison_id == comparison_id), None)

    def _current_comparison_id(self) -> str | None:
        row = self.table.currentRow()
        if row < 0:
            return None
        item = self.table.item(row, 1) or self.table.item(row, 2) or self.table.item(row, 3)
        if item is None:
            return None
        payload = item.data(Qt.ItemDataRole.UserRole)
        if isinstance(payload, tuple):
            return str(payload[0])
        return None

    def _restore_current_selection(self) -> None:
        target_id = self._pending_current_id
        if not target_id:
            if self.table.rowCount() > 0:
                self.table.selectRow(0)
            return
        for row in range(self.table.rowCount()):
            item = self.table.item(row, 1) or self.table.item(row, 2) or self.table.item(row, 3)
            if item is None:
                continue
            payload = item.data(Qt.ItemDataRole.UserRole)
            if isinstance(payload, tuple) and str(payload[0]) == target_id:
                self.table.selectRow(row)
                return
        if self.table.rowCount() > 0:
            self.table.selectRow(0)

    def _update_move_buttons(self) -> None:
        row = self.table.currentRow()
        row_count = self.table.rowCount()
        self.move_up_button.setEnabled(row > 0)
        self.move_down_button.setEnabled(0 <= row < row_count - 1)

    @staticmethod
    def _format_groups(groups: tuple[str | None, str | None]) -> str:
        left, right = groups
        if left and right:
            return f"{left} vs {right}"
        return "missing"
