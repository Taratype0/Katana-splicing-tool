from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QAbstractItemView, QAbstractScrollArea, QLabel, QLineEdit, QPushButton, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget

from src.services.project_service import ProjectService


class ProjectPage(QWidget):
    def __init__(self, project_service: ProjectService, on_confirm=None) -> None:
        super().__init__()
        self.project_service = project_service
        self.on_confirm = on_confirm
        self.summary = QLabel("Select a project directory to begin.")
        self.confirm_button = QPushButton("Confirm")
        self.confirm_button.clicked.connect(self._confirm_project)
        self.confirm_button.setMaximumWidth(120)
        self.confirm_button.setStyleSheet(
            "QPushButton { background-color: #2563EB; color: white; border-radius: 6px; padding: 6px 12px; font-weight: 600; }"
            "QPushButton:disabled { background-color: #64748B; color: #E2E8F0; }"
        )
        self.table = QTableWidget(0, 12)
        self.table.setHorizontalHeaderLabels(
            [
                "Comparison",
                "Display Name",
                "rMATS Name",
                "DEG Name",
                "SUPPA Name",
                "DEXSeq Name",
                "DTU Name",
                "Quant Name",
                "Sashimi Name",
                "BAM Group",
                "rMATS",
                "Isoform",
            ]
        )
        self.table.itemChanged.connect(self._handle_item_changed)
        self.table.setHorizontalScrollMode(QAbstractItemView.ScrollMode.ScrollPerPixel)
        self.table.setVerticalScrollMode(QAbstractItemView.ScrollMode.ScrollPerPixel)
        self.table.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.table.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.table.setSizeAdjustPolicy(QAbstractScrollArea.SizeAdjustPolicy.AdjustIgnored)
        self.table.setMinimumWidth(0)
        self._is_refreshing = False

        layout = QVBoxLayout(self)
        layout.addWidget(self.summary)
        layout.addWidget(self.confirm_button)
        layout.addWidget(self.table, 1)

    def refresh(self) -> None:
        project = self.project_service.current_project
        if project is None:
            self.summary.setText("Select a project directory to begin.")
            self.confirm_button.setEnabled(False)
            self._is_refreshing = True
            self.table.setRowCount(0)
            self._is_refreshing = False
            return

        mode = "manual input mode" if project.input_paths else "auto scan mode"
        inputs = [
            f"rMATS={project.input_paths.get('rmats_root', project.rmats_root or 'auto')}",
            f"DEG={project.input_paths.get('deg_root', project.deg_root or 'auto')}",
            f"quant={project.input_paths.get('quant_root', project.quant_root or 'auto')}",
            f"counts={project.input_paths.get('counts_path', project.counts_path or 'auto')}",
            f"contrastsheet={project.input_paths.get('contrastsheet_path', project.contrastsheet_path or 'auto')}",
        ]
        self.summary.setText(
            f"Project: {project.project_name} | Root: {project.project_root} | "
            f"Comparisons: {len(project.available_comparisons)} | Mode: {mode}\n"
            f"Using: {' | '.join(str(value) for value in inputs)}\n"
            "Step 1: confirm the detected inputs and source names here. "
            "Step 2: go to Pairing to confirm experiment/control direction. "
            "Those names and directions will carry forward into Comparison Sets, Analysis, Shortlist and Cards."
        )
        self.confirm_button.setText("Confirmed" if project.confirmed else "Confirm")
        self.confirm_button.setEnabled(not project.confirmed)
        comparisons = self.project_service.available_comparisons_for_display()
        self._is_refreshing = True
        self.table.blockSignals(True)
        self.table.setRowCount(len(comparisons))
        for row, comparison in enumerate(comparisons):
            self.table.setItem(row, 0, QTableWidgetItem(comparison.comparison_id))
            name_edit = QLineEdit(comparison.resolved_name)
            name_edit.setToolTip("Display name is saved by stable comparison_id, not by temporary row index.")
            name_edit.editingFinished.connect(
                lambda cid=comparison.comparison_id, editor=name_edit: self._commit_display_name_edit(cid, editor)
            )
            self.table.setCellWidget(row, 1, name_edit)
            self._set_mapping_item(row, 2, comparison.comparison_id, comparison.rmats_name or comparison.comparison_id)
            self._set_mapping_item(row, 3, comparison.comparison_id, comparison.deg_name or comparison.comparison_id)
            self._set_mapping_item(row, 4, comparison.comparison_id, comparison.suppa_name or comparison.comparison_id)
            self._set_mapping_item(row, 5, comparison.comparison_id, comparison.dexseq_name or comparison.comparison_id)
            self._set_mapping_item(row, 6, comparison.comparison_id, comparison.dtu_name or comparison.comparison_id)
            self._set_mapping_item(row, 7, comparison.comparison_id, comparison.quant_name or comparison.comparison_id)
            self._set_mapping_item(row, 8, comparison.comparison_id, comparison.sashimi_name or comparison.comparison_id)
            self._set_mapping_item(row, 9, comparison.comparison_id, comparison.bam_group_name or comparison.comparison_id)
            rmats_status_item = QTableWidgetItem(comparison.detected_rmats_modes_label)
            rmats_status_item.setToolTip(comparison.rmats_file_label)
            self.table.setItem(row, 10, rmats_status_item)
            quant_item = QTableWidgetItem("yes" if comparison.has_quant else "no")
            quant_item.setToolTip(str(comparison.quant_dir) if comparison.quant_dir else "missing")
            self.table.setItem(row, 11, quant_item)
        self.table.blockSignals(False)
        self._is_refreshing = False

    def _handle_item_changed(self, item: QTableWidgetItem) -> None:
        if self._is_refreshing:
            return
        comparison_id = item.data(Qt.ItemDataRole.UserRole)
        if comparison_id is None:
            return
        comparison_id = str(comparison_id)
        if item.column() == 1:
            self.project_service.invalidate_project_confirmation()
            self.project_service.rename_comparison(comparison_id, item.text())
            self.refresh()
            return
        mapping_columns = {
            2: "rmats_name",
            3: "deg_name",
            4: "suppa_name",
            5: "dexseq_name",
            6: "dtu_name",
            7: "quant_name",
            8: "sashimi_name",
            9: "bam_group_name",
        }
        field = mapping_columns.get(item.column())
        if field is None:
            return
        self.project_service.invalidate_project_confirmation()
        self.project_service.update_comparison_mapping(comparison_id, **{field: item.text()})
        self.refresh()

    def _commit_display_name_edit(self, comparison_id: str, editor: QLineEdit) -> None:
        if self._is_refreshing:
            return
        comparison = next(
            (item for item in self.project_service.current_project.available_comparisons if item.comparison_id == comparison_id),
            None,
        ) if self.project_service.current_project is not None else None
        if comparison is None:
            return
        new_text = editor.text()
        current_text = comparison.display_name or comparison.resolved_name
        if new_text == current_text:
            return
        self.project_service.invalidate_project_confirmation()
        self.project_service.rename_comparison(comparison_id, new_text)
        self.refresh()

    def _confirm_project(self) -> None:
        self.project_service.confirm_project()
        self.refresh()
        if callable(self.on_confirm):
            self.on_confirm()

    def _set_mapping_item(self, row: int, column: int, comparison_id: str, value: str) -> None:
        item = QTableWidgetItem(value)
        item.setData(Qt.ItemDataRole.UserRole, comparison_id)
        comparison = next(
            (entry for entry in self.project_service.current_project.available_comparisons if entry.comparison_id == comparison_id),
            None,
        ) if self.project_service.current_project is not None else None
        if comparison is not None:
            if column == 2:
                item.setToolTip(comparison.rmats_file_label)
            elif column == 3:
                item.setToolTip(comparison.deg_file_label)
            elif column == 7:
                item.setToolTip(str(comparison.quant_dir) if comparison.quant_dir else "missing")
        self.table.setItem(row, column, item)
