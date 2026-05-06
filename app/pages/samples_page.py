from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QFileDialog,
    QAbstractScrollArea,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from src.services.project_service import ProjectService


class SamplesPage(QWidget):
    def __init__(self, project_service: ProjectService) -> None:
        super().__init__()
        self.project_service = project_service
        self.label = QLabel("Sample mapping")
        self.import_button = QPushButton("Import Sample Metadata")
        self.import_button.clicked.connect(self._import_metadata)
        self.export_button = QPushButton("Export Sample Template")
        self.export_button.clicked.connect(self._export_template)
        self.table = QTableWidget(0, 6)
        self.table.setHorizontalHeaderLabels(
            ["Include", "Sample", "Sample Group", "Condition", "Batch", "quant.sf"]
        )
        self.table.itemChanged.connect(self._handle_item_changed)
        self.table.setHorizontalScrollMode(QAbstractItemView.ScrollMode.ScrollPerPixel)
        self.table.setVerticalScrollMode(QAbstractItemView.ScrollMode.ScrollPerPixel)
        self.table.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.table.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.table.setSizeAdjustPolicy(QAbstractScrollArea.SizeAdjustPolicy.AdjustIgnored)
        self.table.setMinimumWidth(0)
        self.details = QTextEdit()
        self.details.setReadOnly(True)
        self._is_refreshing = False

        layout = QVBoxLayout(self)
        layout.addWidget(self.label)
        button_row = QHBoxLayout()
        button_row.addWidget(self.import_button)
        button_row.addWidget(self.export_button)
        layout.addLayout(button_row)
        layout.addWidget(self.table)
        layout.addWidget(self.details)

    def refresh(self) -> None:
        project = self.project_service.current_project
        if project is None:
            self.table.setRowCount(0)
            self.details.setPlainText("")
            return

        samples = list(project.isoform_samples)
        self._is_refreshing = True
        self.table.setRowCount(len(samples))
        for row, sample in enumerate(samples):
            include_item = QTableWidgetItem()
            include_item.setFlags(
                Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsUserCheckable | Qt.ItemFlag.ItemIsSelectable
            )
            include_item.setCheckState(Qt.CheckState.Checked if sample.include else Qt.CheckState.Unchecked)
            include_item.setData(Qt.ItemDataRole.UserRole, (sample.sample_id, str(sample.quant_sf), "include"))
            self.table.setItem(row, 0, include_item)

            sample_item = QTableWidgetItem(sample.sample_id)
            sample_item.setFlags(sample_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.table.setItem(row, 1, sample_item)

            group_item = QTableWidgetItem(sample.sample_group or "")
            group_item.setData(Qt.ItemDataRole.UserRole, (sample.sample_id, str(sample.quant_sf), "sample_group"))
            self.table.setItem(row, 2, group_item)

            condition_item = QTableWidgetItem(sample.condition or "")
            condition_item.setData(Qt.ItemDataRole.UserRole, (sample.sample_id, str(sample.quant_sf), "condition"))
            self.table.setItem(row, 3, condition_item)

            batch_item = QTableWidgetItem(sample.batch or "")
            batch_item.setData(Qt.ItemDataRole.UserRole, (sample.sample_id, str(sample.quant_sf), "batch"))
            self.table.setItem(row, 4, batch_item)

            quant_item = QTableWidgetItem(str(sample.quant_sf))
            quant_item.setFlags(quant_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.table.setItem(row, 5, quant_item)

        self._is_refreshing = False

        groups = sorted({sample.sample_group for sample in samples if sample.sample_group})
        conditions = sorted({sample.condition for sample in samples if sample.condition})
        self.label.setText(f"Samples ({len(samples)} detected)")
        self.details.setPlainText(
            "\n".join(
                [
                    f"Detected quant.sf files: {len(samples)}",
                    f"Sample groups: {', '.join(groups) if groups else 'not set'}",
                    f"Conditions: {', '.join(conditions) if conditions else 'not set'}",
                    "",
                    "Use this page to normalize sample naming across projects.",
                    "Sample Group should be the biological group name without replicate suffixes.",
                ]
            )
        )

    def _handle_item_changed(self, item: QTableWidgetItem) -> None:
        if self._is_refreshing:
            return
        metadata = item.data(Qt.ItemDataRole.UserRole)
        if not metadata:
            return
        sample_id, quant_sf, field = metadata
        kwargs = {}
        if field == "include":
            kwargs["include"] = item.checkState() == Qt.CheckState.Checked
        else:
            kwargs[field] = item.text()
        self.project_service.update_isoform_sample(sample_id, quant_sf, **kwargs)
        self.refresh()

    def _import_metadata(self) -> None:
        project = self.project_service.current_project
        if project is None:
            return
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Import Sample Metadata",
            str(project.project_root),
            "Tables (*.tsv *.csv *.txt);;All Files (*)",
        )
        if not file_path:
            return
        try:
            updated, created = self.project_service.import_isoform_metadata(file_path)
        except Exception as exc:
            self.details.setPlainText(f"Sample metadata import failed:\n{exc}")
            return
        self.refresh()
        self.details.setPlainText(
            f"Sample metadata imported.\n\nUpdated rows: {updated}\nCreated rows: {created}\nSource: {file_path}"
        )

    def _export_template(self) -> None:
        project = self.project_service.current_project
        if project is None:
            return
        default_path = project.project_root / "sample_metadata_template.tsv"
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export Sample Metadata Template",
            str(default_path),
            "TSV Files (*.tsv);;All Files (*)",
        )
        if not file_path:
            return
        try:
            out = self.project_service.export_isoform_design_template(file_path)
        except Exception as exc:
            self.details.setPlainText(f"Sample metadata export failed:\n{exc}")
            return
        self.details.setPlainText(f"Sample metadata template exported:\n{out}")
