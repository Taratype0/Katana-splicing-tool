from __future__ import annotations

import os

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


class DegWorker(QObject):
    finished = Signal(object)
    failed = Signal(str)

    def __init__(self, service: ProjectService, comparison_id: str) -> None:
        super().__init__()
        self.service = service
        self.comparison_id = comparison_id

    def run(self) -> None:
        try:
            result = self.service.preview_deg_explorer(
                self.comparison_id,
                allow_generate=True,
                force_rebuild=True,
            )
        except Exception as exc:
            self.failed.emit(str(exc))
            return
        self.finished.emit(result)


class DegPage(QWidget):
    def __init__(self, project_service: ProjectService) -> None:
        super().__init__()
        self.project_service = project_service
        self.comparison_tabs = QTabWidget()
        self.comparison_tabs.setDocumentMode(True)
        self.comparison_tabs.setUsesScrollButtons(True)
        self.comparison_tabs.currentChanged.connect(self.refresh)
        self.search = QLineEdit()
        self.search.setPlaceholderText("Search gene symbol / gene ID")
        self.search.textChanged.connect(self.refresh)
        self.class_filter = QComboBox()
        self.class_filter.addItems(["All genes", "Significant only", "Up only", "Down only"])
        self.class_filter.currentTextChanged.connect(self.refresh)
        self.run_button = QPushButton("Run DEG Explorer")
        self.run_button.clicked.connect(self._run)
        self.plot_selector = QComboBox()
        self.plot_selector.addItem("Volcano", "volcano")
        self.plot_selector.addItem("MA plot", "ma")
        self.plot_selector.currentIndexChanged.connect(self._refresh_plot)
        self.open_button = QPushButton("Open output folder")
        self.open_button.clicked.connect(self._open_folder)
        self.export_button = QPushButton("Export filtered table")
        self.export_button.clicked.connect(self._export)
        self.progress = QProgressBar()
        self.progress.setRange(0, 0)
        self.progress.setVisible(False)
        self.status = QLabel("Not run. Select a comparison and click Run DEG Explorer.")
        self.status.setWordWrap(True)
        self.table = QTableWidget(0, 0)
        self.table.setSortingEnabled(True)
        self.plot = QLabel("No cached DEG plot result")
        self.plot.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.plot.setMinimumSize(640, 420)
        plot_scroll = QScrollArea()
        plot_scroll.setWidgetResizable(True)
        plot_scroll.setWidget(self.plot)
        self.views = QTabWidget()
        self.views.addTab(self.table, "Data")
        self.views.addTab(plot_scroll, "Plots / Preview")
        self._frame = pd.DataFrame()
        self._thread: QThread | None = None
        self._worker: DegWorker | None = None

        controls = QHBoxLayout()
        controls.addWidget(self.search, 1)
        controls.addWidget(self.class_filter)
        controls.addWidget(self.plot_selector)
        controls.addWidget(self.run_button)
        controls.addWidget(self.progress)
        controls.addWidget(self.open_button)
        controls.addWidget(self.export_button)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("DEG Explorer", objectName="pageTitle"))
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
                page = QWidget()
                self.comparison_tabs.addTab(page, comparison.display_resolved_name)
                self.comparison_tabs.setTabToolTip(self.comparison_tabs.count() - 1, comparison.comparison_id)
            target = next((index for index, item in enumerate(comparisons) if item.comparison_id == current_id), 0)
            if self.comparison_tabs.count():
                self.comparison_tabs.setCurrentIndex(target)
            self.comparison_tabs.blockSignals(False)
        comparison_id = self._comparison_id()
        if not comparison_id:
            self._show_frame(pd.DataFrame())
            return
        frame = self.project_service.preview_deg_explorer(comparison_id, allow_generate=False)
        state = self.project_service.module_state("deg_explorer")
        self.status.setText(str(state.get("message", "Not run.")))
        self._frame = self._filter(frame)
        self._show_frame(self._frame)
        self._refresh_plot()

    def _comparison_id(self) -> str:
        index = self.comparison_tabs.currentIndex()
        return self.comparison_tabs.tabToolTip(index) if index >= 0 else ""

    def _filter(self, frame: pd.DataFrame) -> pd.DataFrame:
        if frame is None or frame.empty:
            return pd.DataFrame()
        result = frame.copy()
        mode = self.class_filter.currentText()
        if mode == "Significant only":
            result = result.loc[result["DEG_significant"].fillna(False)]
        elif mode == "Up only":
            result = result.loc[result["direction_class"].eq("up")]
        elif mode == "Down only":
            result = result.loc[result["direction_class"].eq("down")]
        query = self.search.text().strip().casefold()
        if query:
            mask = pd.Series(False, index=result.index)
            for column in ("gene_symbol", "gene_id"):
                if column in result.columns:
                    mask |= result[column].astype(str).str.casefold().str.contains(query, regex=False, na=False)
            result = result.loc[mask]
        return result

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
        self.table.setHorizontalHeaderLabels([str(column) for column in show.columns])
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
        self.progress.setVisible(True)
        self.run_button.setEnabled(False)
        self.status.setText(f"Running DEG Explorer for {comparison_id}...")
        thread = QThread(self)
        worker = DegWorker(self.project_service, comparison_id)
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.finished.connect(self._finished)
        worker.failed.connect(self._failed)
        worker.finished.connect(thread.quit)
        worker.failed.connect(thread.quit)
        thread.finished.connect(self._cleanup)
        self._thread, self._worker = thread, worker
        thread.start()

    def _finished(self, _result: object) -> None:
        self.status.setText("DEG Explorer finished.")
        self.refresh()

    def _failed(self, message: str) -> None:
        self.status.setText(f"DEG Explorer failed: {message}")

    def _cleanup(self) -> None:
        self.progress.setVisible(False)
        self.run_button.setEnabled(True)
        if self._worker is not None:
            self._worker.deleteLater()
        if self._thread is not None:
            self._thread.deleteLater()
        self._worker = None
        self._thread = None

    def _refresh_plot(self) -> None:
        files = self.project_service.deg_explorer_file_map(self._comparison_id())
        path = files.get(str(self.plot_selector.currentData() or "volcano"))
        if path is None or not path.exists():
            self.plot.setPixmap(QPixmap())
            self.plot.setText("No cached DEG plot result")
            return
        pixmap = QPixmap(str(path))
        target = self.plot.size()
        self.plot.setText("")
        self.plot.setPixmap(pixmap.scaled(target, Qt.AspectRatioMode.KeepAspectRatio, Qt.TransformationMode.SmoothTransformation))

    def resizeEvent(self, event) -> None:
        super().resizeEvent(event)
        self._refresh_plot()

    def _open_folder(self) -> None:
        files = self.project_service.deg_explorer_file_map(self._comparison_id())
        path = files.get("all")
        if path is not None:
            os.startfile(str(path.parent))

    def _export(self) -> None:
        if self._frame.empty:
            return
        target, _ = QFileDialog.getSaveFileName(self, "Export filtered DEG table", "deg_filtered.tsv", "TSV (*.tsv);;CSV (*.csv)")
        if not target:
            return
        sep = "," if target.lower().endswith(".csv") else "\t"
        self._frame.to_csv(target, sep=sep, index=False)
