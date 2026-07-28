from __future__ import annotations

import os
import re

import pandas as pd
from PySide6.QtCore import QObject, QThread, Signal, Qt
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import (
    QComboBox,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QProgressBar,
    QPushButton,
    QScrollArea,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from src.services.project_service import ProjectService


class EnrichmentWorker(QObject):
    finished = Signal(object)
    failed = Signal(str)

    def __init__(self, service: ProjectService, comparison_id: str, source: str, organism: str, genes: list[str]) -> None:
        super().__init__()
        self.service = service
        self.comparison_id = comparison_id
        self.source = source
        self.organism = organism
        self.genes = genes

    def run(self) -> None:
        try:
            result = self.service.run_functional_enrichment(
                self.comparison_id,
                source=self.source,
                organism=self.organism,
                custom_genes=self.genes,
                force_rebuild=True,
            )
        except Exception as exc:
            self.failed.emit(str(exc))
            return
        self.finished.emit(result)


class EnrichmentPage(QWidget):
    def __init__(self, project_service: ProjectService) -> None:
        super().__init__()
        self.project_service = project_service
        self.comparison_tabs = QTabWidget()
        self.comparison_tabs.setDocumentMode(True)
        self.comparison_tabs.setUsesScrollButtons(True)
        self.comparison_tabs.currentChanged.connect(self.refresh)
        self.source = QComboBox()
        self.source.addItem("Significant DEG", "deg_significant")
        self.source.addItem("Upregulated DEG", "deg_up")
        self.source.addItem("Downregulated DEG", "deg_down")
        self.source.addItem("Candidate gene selection", "candidate_selection")
        self.source.addItem("Custom gene list", "custom")
        self.source.currentIndexChanged.connect(self._update_custom_visibility)
        self.organism = QComboBox()
        self.organism.addItem("Mouse", "mmusculus")
        self.organism.addItem("Human", "hsapiens")
        self.organism.addItem("Rat", "rnorvegicus")
        self.custom_genes = QLineEdit()
        self.custom_genes.setPlaceholderText("Custom gene symbols separated by comma, space, or tab")
        self.custom_genes.setVisible(False)
        self.search = QLineEdit()
        self.search.setPlaceholderText("Search term / source / gene")
        self.search.textChanged.connect(self.refresh)
        self.run_button = QPushButton("Run GO / KEGG Enrichment")
        self.run_button.clicked.connect(self._run)
        self.open_button = QPushButton("Open output folder")
        self.open_button.clicked.connect(self._open_folder)
        self.export_button = QPushButton("Export filtered table")
        self.export_button.clicked.connect(self._export)
        self.progress = QProgressBar()
        self.progress.setRange(0, 0)
        self.progress.setVisible(False)
        self.status = QLabel(
            "Not run. This module performs g:Profiler over-representation analysis on an existing DEG/candidate/custom gene list."
        )
        self.status.setWordWrap(True)
        self.table = QTableWidget(0, 0)
        self.table.setSortingEnabled(True)
        self.plot = QLabel("No cached GO/KEGG plot result")
        self.plot.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.plot.setMinimumSize(640, 420)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(self.plot)
        self.views = QTabWidget()
        self.views.addTab(self.table, "Results")
        self.views.addTab(scroll, "Plots / Preview")
        self._frame = pd.DataFrame()
        self._thread: QThread | None = None
        self._worker: EnrichmentWorker | None = None

        controls = QHBoxLayout()
        controls.addWidget(QLabel("Gene source"))
        controls.addWidget(self.source)
        controls.addWidget(QLabel("Organism"))
        controls.addWidget(self.organism)
        controls.addWidget(self.custom_genes, 1)
        controls.addWidget(self.search, 1)
        controls.addWidget(self.run_button)
        controls.addWidget(self.progress)
        controls.addWidget(self.open_button)
        controls.addWidget(self.export_button)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("GO / KEGG Enrichment", objectName="pageTitle"))
        layout.addWidget(self.comparison_tabs)
        layout.addLayout(controls)
        layout.addWidget(self.status)
        layout.addWidget(self.views, 1)

    def refresh(self) -> None:
        comparisons = self.project_service.selected_comparisons_for_display()
        current_id = self._comparison_id()
        if self.comparison_tabs.count() != len(comparisons) or any(
            self.comparison_tabs.tabToolTip(index) != item.comparison_id for index, item in enumerate(comparisons)
        ):
            self.comparison_tabs.blockSignals(True)
            self.comparison_tabs.clear()
            for comparison in comparisons:
                self.comparison_tabs.addTab(QWidget(), comparison.display_resolved_name)
                self.comparison_tabs.setTabToolTip(self.comparison_tabs.count() - 1, comparison.comparison_id)
            target = next((index for index, item in enumerate(comparisons) if item.comparison_id == current_id), 0)
            if self.comparison_tabs.count():
                self.comparison_tabs.setCurrentIndex(target)
            self.comparison_tabs.blockSignals(False)
        frame = self.project_service.load_cached_functional_enrichment(self._comparison_id())
        query = self.search.text().strip().casefold()
        if query and not frame.empty:
            mask = pd.Series(False, index=frame.index)
            for column in ("term_name", "term_id", "source", "intersection_genes"):
                if column in frame.columns:
                    mask |= frame[column].astype(str).str.casefold().str.contains(query, regex=False, na=False)
            frame = frame.loc[mask]
        self._frame = frame.copy()
        self._show_frame(frame)
        state = self.project_service.module_state("functional_enrichment")
        if not frame.empty or state.get("status") != "not_run":
            self.status.setText(str(state.get("message", "Cached enrichment result.")))
        self._refresh_plot()

    def _comparison_id(self) -> str:
        index = self.comparison_tabs.currentIndex()
        return self.comparison_tabs.tabToolTip(index) if index >= 0 else ""

    def _show_frame(self, frame: pd.DataFrame) -> None:
        self.table.setSortingEnabled(False)
        self.table.clear()
        if frame is None or frame.empty:
            self.table.setRowCount(0)
            self.table.setColumnCount(0)
            self.table.setSortingEnabled(True)
            return
        show = frame.head(2000)
        self.table.setColumnCount(len(show.columns))
        self.table.setHorizontalHeaderLabels(list(show.columns))
        self.table.setRowCount(len(show))
        for row_index, (_, row) in enumerate(show.iterrows()):
            for column_index, column in enumerate(show.columns):
                value = row[column]
                item = QTableWidgetItem("" if pd.isna(value) else str(value))
                if isinstance(value, (int, float)) and not isinstance(value, bool) and pd.notna(value):
                    item.setData(Qt.ItemDataRole.DisplayRole, float(value))
                self.table.setItem(row_index, column_index, item)
        self.table.setSortingEnabled(True)
        self.table.resizeColumnsToContents()

    def _run(self) -> None:
        comparison_id = self._comparison_id()
        if not comparison_id or self._thread is not None:
            return
        genes = [item for item in re.split(r"[\s,;|]+", self.custom_genes.text().strip()) if item]
        thread = QThread(self)
        worker = EnrichmentWorker(
            self.project_service,
            comparison_id,
            str(self.source.currentData()),
            str(self.organism.currentData()),
            genes,
        )
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.finished.connect(self._finished)
        worker.failed.connect(self._failed)
        worker.finished.connect(thread.quit)
        worker.failed.connect(thread.quit)
        thread.finished.connect(self._cleanup)
        self._thread, self._worker = thread, worker
        self.progress.setVisible(True)
        self.run_button.setEnabled(False)
        self.status.setText(f"Running GO/KEGG enrichment for {comparison_id}...")
        thread.start()

    def _finished(self, _result: object) -> None:
        self.refresh()

    def _failed(self, message: str) -> None:
        self.status.setText(f"GO/KEGG enrichment failed: {message}")

    def _cleanup(self) -> None:
        self.progress.setVisible(False)
        self.run_button.setEnabled(True)
        if self._worker is not None:
            self._worker.deleteLater()
        if self._thread is not None:
            self._thread.deleteLater()
        self._worker = None
        self._thread = None

    def _update_custom_visibility(self) -> None:
        self.custom_genes.setVisible(self.source.currentData() == "custom")

    def _refresh_plot(self) -> None:
        path = self.project_service.enrichment_file_map(self._comparison_id()).get("plot")
        if path is None or not path.exists():
            self.plot.setPixmap(QPixmap())
            self.plot.setText("No cached GO/KEGG plot result")
            return
        pixmap = QPixmap(str(path))
        self.plot.setText("")
        self.plot.setPixmap(pixmap.scaled(self.plot.size(), Qt.AspectRatioMode.KeepAspectRatio, Qt.TransformationMode.SmoothTransformation))

    def resizeEvent(self, event) -> None:
        super().resizeEvent(event)
        self._refresh_plot()

    def _open_folder(self) -> None:
        path = self.project_service.enrichment_file_map(self._comparison_id()).get("results")
        if path is not None:
            os.startfile(str(path.parent))

    def _export(self) -> None:
        if self._frame.empty:
            return
        target, _ = QFileDialog.getSaveFileName(self, "Export enrichment table", "go_kegg_filtered.tsv", "TSV (*.tsv);;CSV (*.csv)")
        if not target:
            return
        self._frame.to_csv(target, sep="," if target.lower().endswith(".csv") else "\t", index=False)
