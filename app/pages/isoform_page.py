from __future__ import annotations

import pandas as pd
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
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure

from src.services.project_service import ProjectService


class IsoformPage(QWidget):
    def __init__(self, project_service: ProjectService) -> None:
        super().__init__()
        self.project_service = project_service
        self.label = QLabel("IsoformSwitchAnalyzeR sample design")
        self.explanation = QLabel(
            "Isoform is the follow-up validation layer after DEG/splicing candidate selection. "
            "Use this page to confirm sample conditions and to inspect whether the selected samples support isoform-level follow-up analysis."
        )
        self.explanation.setWordWrap(True)
        self.import_button = QPushButton("Import Metadata")
        self.import_button.clicked.connect(self._import_metadata)
        self.export_button = QPushButton("Export Design Template")
        self.export_button.clicked.connect(self._export_template)
        self.run_button = QPushButton("Run IsoformSwitchAnalyzeR")
        self.run_button.clicked.connect(self._run_isoform_switch)
        self.table = QTableWidget(0, 6)
        self.table.setHorizontalHeaderLabels(
            ["Include", "Sample", "Comparison", "Condition", "Batch", "quant.sf"]
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
        self.figure = Figure(figsize=(8, 3), tight_layout=True)
        self.canvas = FigureCanvasQTAgg(self.figure)
        self._is_refreshing = False
        self._frame = pd.DataFrame()

        layout = QVBoxLayout(self)
        layout.addWidget(self.label)
        layout.addWidget(self.explanation)
        button_row = QHBoxLayout()
        button_row.addWidget(self.import_button)
        button_row.addWidget(self.export_button)
        button_row.addWidget(self.run_button)
        layout.addLayout(button_row)
        layout.addWidget(self.canvas)
        layout.addWidget(self.table)
        layout.addWidget(self.details)

    def refresh(self) -> None:
        project = self.project_service.current_project
        if project is None:
            self.table.setRowCount(0)
            self.details.setPlainText("")
            return

        samples = self.project_service.selected_isoform_samples()
        if not samples:
            samples = list(project.isoform_samples)
        frame = self.project_service.isoform_adapter.build_import_manifest_frame(samples) if samples else pd.DataFrame()
        self._frame = frame

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

            comparison_item = QTableWidgetItem(sample.comparison_id or "")
            comparison_item.setData(Qt.ItemDataRole.UserRole, (sample.sample_id, str(sample.quant_sf), "comparison_id"))
            self.table.setItem(row, 2, comparison_item)

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

        included = sum(1 for sample in samples if sample.include)
        conditions = sorted({sample.condition for sample in samples if sample.include and sample.condition})
        self._draw_summary(samples)
        self.label.setText(f"IsoformSwitchAnalyzeR sample design ({included} included samples)")
        self.details.setPlainText(
            "\n".join(
                [
                    f"Detected transcript-quant samples: {len(project.isoform_samples)}",
                    f"Included samples: {included}",
                    f"Unique conditions: {', '.join(conditions) if conditions else 'not set'}",
                    "",
                    "Required before running:",
                    "- use the Samples page to confirm sample groups and conditions",
                    "- include the samples you want to analyze",
                    "- assign a condition for each included sample",
                    "- provide at least 2 conditions across included samples",
                    "- configure Rscript and GTF in Settings",
                    "- use isoform as follow-up validation after candidate-gene classification",
                    "",
                    "Metadata import columns:",
                    "- sample_id (required)",
                    "- quant_sf (recommended when sample names are reused)",
                    "- comparison_id, condition, batch, include (optional)",
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

    def _run_isoform_switch(self) -> None:
        try:
            output = self.project_service.run_isoform_switch_pipeline()
        except Exception as exc:
            self.details.setPlainText(f"IsoformSwitchAnalyzeR failed:\n{exc}")
            return
        self.refresh()
        self.details.setPlainText("IsoformSwitchAnalyzeR completed.\n\n" + output)

    def _draw_summary(self, samples) -> None:
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        included = [sample for sample in samples if sample.include]
        if not included:
            axis.text(0.5, 0.5, "No included samples for isoform analysis.", ha="center", va="center")
            axis.axis("off")
            self.canvas.draw_idle()
            return
        counts: dict[str, int] = {}
        for sample in included:
            label = sample.condition or sample.sample_group or sample.sample_id
            counts[label] = counts.get(label, 0) + 1
        axis.bar(list(counts.keys()), list(counts.values()), color="#5E60CE")
        axis.set_ylabel("Included samples")
        axis.set_title("Isoform follow-up overview")
        axis.tick_params(axis="x", rotation=20)
        self.canvas.draw_idle()

    def _import_metadata(self) -> None:
        project = self.project_service.current_project
        if project is None:
            return
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Import Isoform Metadata",
            str(project.project_root),
            "Tables (*.tsv *.csv *.txt);;All Files (*)",
        )
        if not file_path:
            return
        try:
            updated, created = self.project_service.import_isoform_metadata(file_path)
        except Exception as exc:
            self.details.setPlainText(f"Isoform metadata import failed:\n{exc}")
            return
        self.refresh()
        self.details.setPlainText(
            f"Isoform metadata imported.\n\nUpdated rows: {updated}\nCreated rows: {created}\nSource: {file_path}"
        )

    def _export_template(self) -> None:
        project = self.project_service.current_project
        if project is None:
            return
        default_path = project.project_root / "isoform_design_template.tsv"
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export Isoform Design Template",
            str(default_path),
            "TSV Files (*.tsv);;All Files (*)",
        )
        if not file_path:
            return
        try:
            out = self.project_service.export_isoform_design_template(file_path)
        except Exception as exc:
            self.details.setPlainText(f"Isoform design export failed:\n{exc}")
            return
        self.details.setPlainText(f"Isoform design template exported:\n{out}")
