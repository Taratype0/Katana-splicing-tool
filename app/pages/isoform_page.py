from __future__ import annotations

import os
import shutil
import pandas as pd
from pandas.api.types import is_bool_dtype, is_numeric_dtype
from PySide6.QtCore import QObject, QThread, QTimer, Signal
from PySide6.QtCore import Qt
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import (
    QAbstractItemView,
    QDialog,
    QFileDialog,
    QAbstractScrollArea,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QProgressBar,
    QScrollArea,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure

from app.pages.module_details import module_detail
from src.services.project_service import ProjectService


class IsoformSupportWorker(QObject):
    finished = Signal(object)
    failed = Signal(str)

    def __init__(self, project_service: ProjectService, comparison_id: str, top_n: int) -> None:
        super().__init__()
        self.project_service = project_service
        self.comparison_id = comparison_id
        self.top_n = top_n

    def run(self) -> None:
        try:
            tables = self.project_service.preview_isoform_candidate_support(
                self.comparison_id,
                source="selection",
                top_n=self.top_n,
                allow_generate=True,
                force_rebuild=True,
            )
        except Exception as exc:
            self.failed.emit(str(exc))
            return
        self.finished.emit(tables)


class IsoformCrossSupportWorker(QObject):
    finished = Signal(object)
    failed = Signal(str)

    def __init__(self, project_service: ProjectService, top_n: int) -> None:
        super().__init__()
        self.project_service = project_service
        self.top_n = top_n

    def run(self) -> None:
        try:
            tables = self.project_service.preview_isoform_cross_comparison_support(
                source="selection",
                top_n=self.top_n,
                allow_generate=True,
                force_rebuild=True,
            )
        except Exception as exc:
            self.failed.emit(str(exc))
            return
        self.finished.emit(tables)


class IsoformPage(QWidget):
    def __init__(self, project_service: ProjectService, initial_tab: int = 0) -> None:
        super().__init__()
        self.project_service = project_service
        self.label = QLabel("8. Isoform / transcript usage support")
        self.explanation = QLabel(
            "8 focuses on transcript/isoform evidence from quant.sf, GTF and optional DTU results. "
            "This is separate from 7 RBP support: it asks whether candidate genes have isoform usage or transcript-abundance support."
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
        self.support_title = QLabel("Candidate isoform support")
        self.support_title.setStyleSheet("font-size: 16px; font-weight: 700;")
        self.support_context = QLabel(
            "8.1 uses selected candidate genes plus quant.sf, GTF transcript annotation and optional DTU q values. "
            "It is transcript/isoform support; DEG+rMATS remain the core candidate definition."
        )
        self.support_context.setWordWrap(True)
        self.support_comparison_tabs = QTabWidget()
        self.support_comparison_tabs.setDocumentMode(True)
        self.support_comparison_tabs.setUsesScrollButtons(True)
        self.support_comparison_tabs.setMaximumHeight(36)
        self.support_comparison_tabs.currentChanged.connect(self._refresh_support_cached)
        self.support_search = QLineEdit()
        self.support_search.setPlaceholderText("Search gene / transcript / status...")
        self.support_search.textChanged.connect(self._apply_support_filter)
        self.support_top_n = QSpinBox()
        self.support_top_n.setRange(1, 5000)
        self.support_top_n.setValue(50)
        self.run_support_button = QPushButton("Run Isoform support")
        self.run_support_button.clicked.connect(self._run_isoform_support)
        self.support_detail_button = QPushButton("Detail")
        self.support_detail_button.clicked.connect(lambda: self._show_detail_dialog("8.isoform"))
        self.open_support_button = QPushButton("Open output folder")
        self.open_support_button.clicked.connect(self._open_isoform_support_folder)
        self.export_support_button = QPushButton("Export current table")
        self.export_support_button.clicked.connect(self._export_isoform_support_table)
        self.support_progress = QProgressBar()
        self.support_progress.setRange(0, 0)
        self.support_progress.setMaximumWidth(170)
        self.support_progress.setVisible(False)
        self.support_table_tabs = QTabWidget()
        self.gene_summary_table = QTableWidget(0, 0)
        self.gene_summary_table.setSortingEnabled(True)
        self.transcript_table = QTableWidget(0, 0)
        self.transcript_table.setSortingEnabled(True)
        self.support_table_tabs.addTab(self.gene_summary_table, "Gene summary")
        self.support_table_tabs.addTab(self.transcript_table, "Transcript table")
        self.support_status = QLabel("")
        self.support_status.setWordWrap(True)
        self._support_tables: dict[str, pd.DataFrame] = {"gene_summary": pd.DataFrame(), "transcript_table": pd.DataFrame()}
        self._support_thread: QThread | None = None
        self._support_worker: IsoformSupportWorker | None = None
        self.cross_title = QLabel("Cross-comparison isoform support")
        self.cross_title.setStyleSheet("font-size: 16px; font-weight: 700;")
        self.cross_context = QLabel(
            "8.2 compares transcript/isoform support patterns across selected comparisons using per-comparison quant.sf + GTF + optional DTU support. "
            "It is a support-prioritization comparison, not an IsoformSwitchAnalyzeR causal model."
        )
        self.cross_context.setWordWrap(True)
        self.cross_search = QLineEdit()
        self.cross_search.setPlaceholderText("Search gene / transcript / support pattern...")
        self.cross_search.textChanged.connect(self._apply_cross_filter)
        self.cross_top_n = QSpinBox()
        self.cross_top_n.setRange(1, 5000)
        self.cross_top_n.setValue(50)
        self.run_cross_button = QPushButton("Run Cross-comparison Isoform support")
        self.run_cross_button.clicked.connect(self._run_isoform_cross_support)
        self.cross_detail_button = QPushButton("Detail")
        self.cross_detail_button.clicked.connect(lambda: self._show_detail_dialog("8.isoform"))
        self.open_cross_button = QPushButton("Open output folder")
        self.open_cross_button.clicked.connect(self._open_isoform_cross_folder)
        self.export_cross_button = QPushButton("Export current table")
        self.export_cross_button.clicked.connect(self._export_isoform_cross_table)
        self.cross_progress = QProgressBar()
        self.cross_progress.setRange(0, 0)
        self.cross_progress.setMaximumWidth(170)
        self.cross_progress.setVisible(False)
        self.cross_table_tabs = QTabWidget()
        self.cross_gene_table = QTableWidget(0, 0)
        self.cross_gene_table.setSortingEnabled(True)
        self.cross_transcript_table = QTableWidget(0, 0)
        self.cross_transcript_table.setSortingEnabled(True)
        self.cross_switch_table = QTableWidget(0, 0)
        self.cross_switch_table.setSortingEnabled(True)
        self.cross_plot_table = QTableWidget(0, 2)
        self.cross_plot_table.setHorizontalHeaderLabels(["Plot", "Path"])
        self.cross_plot_table.setSortingEnabled(True)
        self.cross_plot_table.currentCellChanged.connect(self._preview_isoform_cross_plot)
        self.cross_preview_image = QLabel("No cached cross-comparison plot result.")
        self.cross_preview_image.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.cross_preview_image.setMinimumHeight(320)
        self.cross_preview_scroll = QScrollArea()
        self.cross_preview_scroll.setWidgetResizable(True)
        self.cross_preview_scroll.setWidget(self.cross_preview_image)
        self.cross_export_plot_button = QPushButton("Export selected plot")
        self.cross_export_plot_button.clicked.connect(self._export_isoform_cross_plot)
        self.cross_table_tabs.addTab(self.cross_gene_table, "Gene support matrix")
        self.cross_table_tabs.addTab(self.cross_transcript_table, "Transcript support matrix")
        self.cross_table_tabs.addTab(self.cross_switch_table, "Dominant isoform switches")
        cross_plot_tab = QWidget()
        cross_plot_layout = QVBoxLayout(cross_plot_tab)
        cross_plot_controls = QHBoxLayout()
        cross_plot_controls.addWidget(self.cross_export_plot_button)
        cross_plot_controls.addStretch(1)
        cross_plot_layout.addLayout(cross_plot_controls)
        cross_plot_layout.addWidget(self.cross_plot_table)
        cross_plot_layout.addWidget(self.cross_preview_scroll, 1)
        self.cross_table_tabs.addTab(cross_plot_tab, "Plots / Preview")
        self.cross_status = QLabel("")
        self.cross_status.setWordWrap(True)
        self._cross_tables: dict[str, pd.DataFrame] = {
            "gene_support_matrix": pd.DataFrame(),
            "transcript_support_matrix": pd.DataFrame(),
            "dominant_isoform_switches": pd.DataFrame(),
        }
        self._cross_plot_frame = pd.DataFrame()
        self._cross_preview_pixmap = QPixmap()
        self._cross_thread: QThread | None = None
        self._cross_worker: IsoformCrossSupportWorker | None = None

        sample_tab = QWidget()
        sample_layout = QVBoxLayout(sample_tab)
        sample_layout.addWidget(self.label)
        sample_layout.addWidget(self.explanation)
        button_row = QHBoxLayout()
        button_row.addWidget(self.import_button)
        button_row.addWidget(self.export_button)
        button_row.addWidget(self.run_button)
        sample_layout.addLayout(button_row)
        sample_layout.addWidget(self.canvas)
        sample_layout.addWidget(self.table)
        sample_layout.addWidget(self.details)

        support_tab = QWidget()
        support_layout = QVBoxLayout(support_tab)
        support_controls = QHBoxLayout()
        support_controls.addWidget(QLabel("Top N"))
        support_controls.addWidget(self.support_top_n)
        support_controls.addWidget(self.support_search, 1)
        support_controls.addWidget(self.run_support_button)
        support_controls.addWidget(self.support_detail_button)
        support_controls.addWidget(self.support_progress)
        support_controls.addWidget(self.open_support_button)
        support_controls.addWidget(self.export_support_button)
        support_layout.addWidget(self.support_title)
        support_layout.addWidget(self.support_context)
        support_layout.addWidget(self.support_comparison_tabs)
        support_layout.addLayout(support_controls)
        support_layout.addWidget(self.support_table_tabs, 1)
        support_layout.addWidget(self.support_status)

        cross_tab = QWidget()
        cross_layout = QVBoxLayout(cross_tab)
        cross_controls = QHBoxLayout()
        cross_controls.addWidget(QLabel("Top N"))
        cross_controls.addWidget(self.cross_top_n)
        cross_controls.addWidget(self.cross_search, 1)
        cross_controls.addWidget(self.run_cross_button)
        cross_controls.addWidget(self.cross_detail_button)
        cross_controls.addWidget(self.cross_progress)
        cross_controls.addWidget(self.open_cross_button)
        cross_controls.addWidget(self.export_cross_button)
        cross_layout.addWidget(self.cross_title)
        cross_layout.addWidget(self.cross_context)
        cross_layout.addLayout(cross_controls)
        cross_layout.addWidget(self.cross_table_tabs, 1)
        cross_layout.addWidget(self.cross_status)

        self.main_tabs = QTabWidget()
        self.main_tabs.addTab(sample_tab, "Sample design")
        self.main_tabs.addTab(support_tab, "8.1 Per-comparison isoform support")
        self.main_tabs.addTab(cross_tab, "8.2 Cross-comparison isoform support")
        self.main_tabs.setCurrentIndex(max(0, min(initial_tab, self.main_tabs.count() - 1)))

        layout = QVBoxLayout(self)
        layout.addWidget(self.main_tabs)

    def navigate_to(self, _branch: str | None = None, subbranch: str | None = None) -> None:
        if subbranch == "per":
            self.main_tabs.setCurrentIndex(1)
        elif subbranch == "cross":
            self.main_tabs.setCurrentIndex(2)
        elif subbranch == "design":
            self.main_tabs.setCurrentIndex(0)

    def _show_detail_dialog(self, key: str) -> None:
        text = module_detail(key)
        if not text:
            return
        dialog = QDialog(self)
        dialog.setWindowTitle("Module details")
        dialog.resize(760, 420)
        layout = QVBoxLayout(dialog)
        viewer = QTextEdit(dialog)
        viewer.setReadOnly(True)
        viewer.setPlainText(text)
        close_button = QPushButton("Close", dialog)
        close_button.clicked.connect(dialog.accept)
        layout.addWidget(viewer, 1)
        layout.addWidget(close_button)
        QTimer.singleShot(0, dialog.showMaximized)
        dialog.exec()

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
        self._refresh_support_cached()
        self._refresh_cross_cached()

    def _populate_support_comparison_tabs(self) -> None:
        comparisons = self.project_service.selected_comparisons_for_display()
        labels = [self.support_comparison_tabs.tabText(i) for i in range(self.support_comparison_tabs.count())]
        wanted = [comparison.display_resolved_name for comparison in comparisons]
        if labels == wanted:
            return
        current_id = self._current_support_comparison_id()
        self.support_comparison_tabs.blockSignals(True)
        self.support_comparison_tabs.clear()
        for comparison in comparisons:
            tab = QWidget()
            tab.setProperty("comparison_id", comparison.comparison_id)
            self.support_comparison_tabs.addTab(tab, comparison.display_resolved_name)
        if current_id:
            for index, comparison in enumerate(comparisons):
                if comparison.comparison_id == current_id:
                    self.support_comparison_tabs.setCurrentIndex(index)
                    break
        self.support_comparison_tabs.blockSignals(False)

    def _current_support_comparison_id(self) -> str:
        widget = self.support_comparison_tabs.currentWidget()
        return str(widget.property("comparison_id") or "") if widget is not None else ""

    def _refresh_support_cached(self) -> None:
        project = self.project_service.current_project
        if project is None:
            self._set_support_tables({"gene_summary": pd.DataFrame(), "transcript_table": pd.DataFrame()})
            self.support_status.setText("No project loaded.")
            return
        self._populate_support_comparison_tabs()
        comparison_id = self._current_support_comparison_id()
        if not comparison_id:
            self.support_status.setText("No comparison selected.")
            return
        tables = self.project_service.preview_isoform_candidate_support(
            comparison_id,
            source="selection",
            top_n=self.support_top_n.value(),
            allow_generate=False,
            force_rebuild=False,
        )
        self._set_support_tables(tables)
        state = self.project_service.module_state("isoform_support")
        self.support_status.setText(f"{state.get('status', 'not_run')}: {state.get('message', '')}")

    def _run_isoform_support(self) -> None:
        comparison_id = self._current_support_comparison_id()
        if not comparison_id:
            self.support_status.setText("No comparison selected.")
            return
        self.support_progress.setVisible(True)
        self.run_support_button.setEnabled(False)
        self._support_thread = QThread(self)
        self._support_worker = IsoformSupportWorker(self.project_service, comparison_id, self.support_top_n.value())
        self._support_worker.moveToThread(self._support_thread)
        self._support_thread.started.connect(self._support_worker.run)
        self._support_worker.finished.connect(self._on_isoform_support_finished)
        self._support_worker.failed.connect(self._on_isoform_support_failed)
        self._support_worker.finished.connect(self._support_thread.quit)
        self._support_worker.failed.connect(self._support_thread.quit)
        self._support_thread.finished.connect(self._support_thread.deleteLater)
        self._support_thread.start()

    def _on_isoform_support_finished(self, tables: object) -> None:
        self.support_progress.setVisible(False)
        self.run_support_button.setEnabled(True)
        self._set_support_tables(tables if isinstance(tables, dict) else {})
        state = self.project_service.module_state("isoform_support")
        self.support_status.setText(f"{state.get('status', 'finished')}: {state.get('message', '')}")

    def _on_isoform_support_failed(self, message: str) -> None:
        self.support_progress.setVisible(False)
        self.run_support_button.setEnabled(True)
        self.support_status.setText(f"Isoform support failed: {message}")

    def _refresh_cross_cached(self) -> None:
        if self.project_service.current_project is None:
            self._set_cross_tables({})
            self._set_isoform_cross_plot_files(pd.DataFrame())
            self.cross_status.setText("No project loaded.")
            return
        tables = self.project_service.preview_isoform_cross_comparison_support(
            source="selection",
            top_n=self.cross_top_n.value(),
            allow_generate=False,
            force_rebuild=False,
        )
        self._set_cross_tables(tables)
        self._refresh_cross_plots()
        state = self.project_service.module_state("isoform_cross_support")
        self.cross_status.setText(f"{state.get('status', 'not_run')}: {state.get('message', '')}")

    def _run_isoform_cross_support(self) -> None:
        if self.project_service.current_project is None:
            self.cross_status.setText("No project loaded.")
            return
        self.cross_progress.setVisible(True)
        self.run_cross_button.setEnabled(False)
        self._cross_thread = QThread(self)
        self._cross_worker = IsoformCrossSupportWorker(self.project_service, self.cross_top_n.value())
        self._cross_worker.moveToThread(self._cross_thread)
        self._cross_thread.started.connect(self._cross_worker.run)
        self._cross_worker.finished.connect(self._on_isoform_cross_finished)
        self._cross_worker.failed.connect(self._on_isoform_cross_failed)
        self._cross_worker.finished.connect(self._cross_thread.quit)
        self._cross_worker.failed.connect(self._cross_thread.quit)
        self._cross_thread.finished.connect(self._cross_thread.deleteLater)
        self._cross_thread.start()

    def _on_isoform_cross_finished(self, tables: object) -> None:
        self.cross_progress.setVisible(False)
        self.run_cross_button.setEnabled(True)
        self._set_cross_tables(tables if isinstance(tables, dict) else {})
        self._refresh_cross_plots()
        state = self.project_service.module_state("isoform_cross_support")
        self.cross_status.setText(f"{state.get('status', 'finished')}: {state.get('message', '')}")

    def _on_isoform_cross_failed(self, message: str) -> None:
        self.cross_progress.setVisible(False)
        self.run_cross_button.setEnabled(True)
        self.cross_status.setText(f"Cross-comparison isoform support failed: {message}")

    def _set_support_tables(self, tables: dict[str, pd.DataFrame]) -> None:
        self._support_tables = {
            "gene_summary": tables.get("gene_summary", pd.DataFrame()).copy() if isinstance(tables, dict) else pd.DataFrame(),
            "transcript_table": tables.get("transcript_table", pd.DataFrame()).copy() if isinstance(tables, dict) else pd.DataFrame(),
        }
        self._apply_support_filter()

    def _apply_support_filter(self) -> None:
        query = self.support_search.text().strip().casefold()
        for key, table in (("gene_summary", self.gene_summary_table), ("transcript_table", self.transcript_table)):
            frame = self._support_tables.get(key, pd.DataFrame()).copy()
            if query and not frame.empty:
                mask = pd.Series(False, index=frame.index)
                for column in frame.columns:
                    mask |= frame[column].fillna("").astype(str).str.casefold().str.contains(query, regex=False)
                frame = frame.loc[mask].copy()
            self._populate_dataframe_table(table, frame)

    def _set_cross_tables(self, tables: dict[str, pd.DataFrame]) -> None:
        self._cross_tables = {
            "gene_support_matrix": tables.get("gene_support_matrix", pd.DataFrame()).copy() if isinstance(tables, dict) else pd.DataFrame(),
            "transcript_support_matrix": tables.get("transcript_support_matrix", pd.DataFrame()).copy() if isinstance(tables, dict) else pd.DataFrame(),
            "dominant_isoform_switches": tables.get("dominant_isoform_switches", pd.DataFrame()).copy() if isinstance(tables, dict) else pd.DataFrame(),
        }
        self._apply_cross_filter()

    def _refresh_cross_plots(self) -> None:
        if self.project_service.current_project is None:
            self._set_isoform_cross_plot_files(pd.DataFrame())
            return
        self._set_isoform_cross_plot_files(self.project_service.isoform_cross_support_plot_files())

    def _set_isoform_cross_plot_files(self, frame: pd.DataFrame) -> None:
        self._cross_plot_frame = frame.copy() if isinstance(frame, pd.DataFrame) else pd.DataFrame()
        self.cross_plot_table.setSortingEnabled(False)
        self.cross_plot_table.setRowCount(0)
        if self._cross_plot_frame.empty:
            self._cross_preview_pixmap = QPixmap()
            self.cross_preview_image.setPixmap(QPixmap())
            self.cross_preview_image.setText("No cached cross-comparison plot result.")
            self.cross_plot_table.setSortingEnabled(True)
            return
        self.cross_plot_table.setRowCount(len(self._cross_plot_frame))
        for row_index, (_, row) in enumerate(self._cross_plot_frame.iterrows()):
            self.cross_plot_table.setItem(row_index, 0, QTableWidgetItem(str(row.get("plot_name", ""))))
            self.cross_plot_table.setItem(row_index, 1, QTableWidgetItem(str(row.get("absolute_path", ""))))
        self.cross_plot_table.resizeColumnsToContents()
        self.cross_plot_table.setSortingEnabled(True)
        self.cross_plot_table.selectRow(0)
        self._preview_isoform_cross_plot(0, 0, -1, -1)

    def _selected_isoform_cross_plot_path(self) -> str:
        row = self.cross_plot_table.currentRow()
        if row < 0:
            return ""
        item = self.cross_plot_table.item(row, 1)
        return item.text() if item is not None else ""

    def _preview_isoform_cross_plot(self, _current_row: int, _current_column: int, _previous_row: int, _previous_column: int) -> None:
        path = self._selected_isoform_cross_plot_path()
        if not path:
            return
        pixmap = QPixmap(path)
        if pixmap.isNull():
            self._cross_preview_pixmap = QPixmap()
            self.cross_preview_image.setPixmap(QPixmap())
            self.cross_preview_image.setText(f"Plot preview unavailable:\n{path}")
            return
        self._cross_preview_pixmap = pixmap
        self.cross_preview_image.setText("")
        self._fit_isoform_cross_preview_pixmap()

    def _fit_isoform_cross_preview_pixmap(self) -> None:
        if self._cross_preview_pixmap.isNull():
            return
        viewport = self.cross_preview_scroll.viewport().size()
        if viewport.width() <= 10 or viewport.height() <= 10:
            return
        scaled = self._cross_preview_pixmap.scaled(
            max(200, viewport.width() - 20),
            max(200, viewport.height() - 20),
            Qt.AspectRatioMode.KeepAspectRatio,
            Qt.TransformationMode.SmoothTransformation,
        )
        self.cross_preview_image.setPixmap(scaled)

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        QTimer.singleShot(0, self._fit_isoform_cross_preview_pixmap)

    def _apply_cross_filter(self) -> None:
        query = self.cross_search.text().strip().casefold()
        mapping = (
            ("gene_support_matrix", self.cross_gene_table),
            ("transcript_support_matrix", self.cross_transcript_table),
            ("dominant_isoform_switches", self.cross_switch_table),
        )
        for key, table in mapping:
            frame = self._cross_tables.get(key, pd.DataFrame()).copy()
            if query and not frame.empty:
                mask = pd.Series(False, index=frame.index)
                for column in frame.columns:
                    mask |= frame[column].fillna("").astype(str).str.casefold().str.contains(query, regex=False)
                frame = frame.loc[mask].copy()
            self._populate_dataframe_table(table, frame)

    def _open_isoform_support_folder(self) -> None:
        folder = self.project_service.isoform_support_output_dir(self._current_support_comparison_id())
        if folder is None:
            return
        folder.mkdir(parents=True, exist_ok=True)
        os.startfile(str(folder))

    def _open_isoform_cross_folder(self) -> None:
        folder = self.project_service.isoform_cross_support_output_dir()
        if folder is None:
            return
        folder.mkdir(parents=True, exist_ok=True)
        os.startfile(str(folder))

    def _export_isoform_support_table(self) -> None:
        key = "gene_summary" if self.support_table_tabs.currentIndex() == 0 else "transcript_table"
        frame = self._support_tables.get(key, pd.DataFrame())
        if frame.empty:
            return
        project = self.project_service.current_project
        start = str(project.project_root / f"{key}.tsv") if project is not None else f"{key}.tsv"
        path, _ = QFileDialog.getSaveFileName(self, "Export isoform support table", start, "TSV Files (*.tsv);;All Files (*)")
        if path:
            frame.to_csv(path, sep="\t", index=False)

    def _export_isoform_cross_table(self) -> None:
        keys = ["gene_support_matrix", "transcript_support_matrix", "dominant_isoform_switches"]
        if self.cross_table_tabs.currentIndex() >= len(keys):
            self._export_isoform_cross_plot()
            return
        key = keys[self.cross_table_tabs.currentIndex()]
        frame = self._cross_tables.get(key, pd.DataFrame())
        if frame.empty:
            return
        project = self.project_service.current_project
        start = str(project.project_root / f"{key}.tsv") if project is not None else f"{key}.tsv"
        path, _ = QFileDialog.getSaveFileName(self, "Export cross-comparison isoform table", start, "TSV Files (*.tsv);;All Files (*)")
        if path:
            frame.to_csv(path, sep="\t", index=False)

    def _export_isoform_cross_plot(self) -> None:
        source = self._selected_isoform_cross_plot_path()
        if not source or not os.path.exists(source):
            return
        project = self.project_service.current_project
        start_dir = project.project_root if project is not None else None
        default_name = os.path.basename(source)
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Export selected isoform plot",
            str(start_dir / default_name) if start_dir is not None else default_name,
            "PNG Files (*.png);;All Files (*)",
        )
        if path:
            shutil.copy2(source, path)

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

    @staticmethod
    def _populate_dataframe_table(table: QTableWidget, frame: pd.DataFrame) -> None:
        table.setSortingEnabled(False)
        table.setRowCount(0)
        table.setColumnCount(0)
        if frame is None or frame.empty:
            table.setSortingEnabled(True)
            return
        table.setColumnCount(len(frame.columns))
        table.setHorizontalHeaderLabels([str(column) for column in frame.columns])
        table.setRowCount(len(frame))
        for row_index, (_, row) in enumerate(frame.iterrows()):
            for col_index, column in enumerate(frame.columns):
                value = row[column]
                item = QTableWidgetItem("" if pd.isna(value) else str(value))
                if is_numeric_dtype(frame[column]):
                    item.setData(Qt.ItemDataRole.EditRole, float(value) if pd.notna(value) else None)
                elif is_bool_dtype(frame[column]):
                    item.setData(Qt.ItemDataRole.EditRole, bool(value))
                table.setItem(row_index, col_index, item)
        table.resizeColumnsToContents()
        table.setSortingEnabled(True)
