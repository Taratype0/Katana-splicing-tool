from __future__ import annotations

import os
import shutil

import pandas as pd
from pandas.api.types import is_bool_dtype, is_numeric_dtype
from PySide6.QtCore import QObject, QThread, QTimer, Signal, Qt
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import (
    QFileDialog,
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
    QVBoxLayout,
    QWidget,
)

from src.services.project_service import ProjectService


class RbpSupportWorker(QObject):
    finished = Signal(object)
    failed = Signal(str)

    def __init__(
        self,
        project_service: ProjectService,
        comparison_id: str,
        source: str,
        top_n: int,
        custom_gene_set_path: str,
    ) -> None:
        super().__init__()
        self.project_service = project_service
        self.comparison_id = comparison_id
        self.source = source
        self.top_n = top_n
        self.custom_gene_set_path = custom_gene_set_path

    def run(self) -> None:
        try:
            tables = self.project_service.preview_rbp_support(
                self.comparison_id,
                source=self.source,
                top_n=self.top_n,
                custom_gene_set_path=self.custom_gene_set_path or None,
                allow_generate=True,
                force_rebuild=True,
            )
        except Exception as exc:
            self.failed.emit(str(exc))
            return
        self.finished.emit(tables)


class RbpCrossSupportWorker(QObject):
    finished = Signal(object)
    failed = Signal(str)

    def __init__(
        self,
        project_service: ProjectService,
        top_n: int,
        custom_gene_set_path: str,
    ) -> None:
        super().__init__()
        self.project_service = project_service
        self.top_n = top_n
        self.custom_gene_set_path = custom_gene_set_path

    def run(self) -> None:
        try:
            tables = self.project_service.preview_rbp_cross_comparison_support(
                source="selection",
                top_n=self.top_n,
                custom_gene_set_path=self.custom_gene_set_path or None,
                allow_generate=True,
                force_rebuild=True,
            )
        except Exception as exc:
            self.failed.emit(str(exc))
            return
        self.finished.emit(tables)


class RbpPage(QWidget):
    def __init__(self, project_service: ProjectService) -> None:
        super().__init__()
        self.project_service = project_service
        self.title = QLabel("RBP support / enrichment")
        self.title.setStyleSheet("font-size: 16px; font-weight: 700;")
        self.context = QLabel(
            "RBP support is a prioritization layer: it checks whether current candidate genes overlap known or custom RBP/splicing-factor gene sets. "
            "It does not claim causal RBP regulation."
        )
        self.context.setWordWrap(True)
        self.comparison_tabs = QTabWidget()
        self.comparison_tabs.setDocumentMode(True)
        self.comparison_tabs.setUsesScrollButtons(True)
        self.comparison_tabs.setMaximumHeight(36)
        self.comparison_tabs.currentChanged.connect(self.refresh)

        self.search = QLineEdit()
        self.search.setPlaceholderText("Search RBP set / gene / evidence...")
        self.search.textChanged.connect(self._apply_filters)
        self.top_n = QSpinBox()
        self.top_n.setRange(1, 5000)
        self.top_n.setValue(50)
        self.source_label = QLabel("Source: candidate selection, fallback to cached top N if no selection")
        self.custom_path = QLineEdit()
        self.custom_path.setPlaceholderText("Optional custom RBP gene set TSV/CSV")
        self.choose_button = QPushButton("Choose gene set")
        self.choose_button.clicked.connect(self._choose_gene_set)
        self.run_button = QPushButton("Run RBP support")
        self.run_button.clicked.connect(self._run)
        self.open_button = QPushButton("Open output folder")
        self.open_button.clicked.connect(self._open_output_folder)
        self.export_button = QPushButton("Export current table")
        self.export_button.clicked.connect(self._export_current_table)
        self.progress = QProgressBar()
        self.progress.setRange(0, 0)
        self.progress.setVisible(False)
        self.progress.setMaximumWidth(170)

        self.table_tabs = QTabWidget()
        self.summary_table = QTableWidget(0, 0)
        self.summary_table.setSortingEnabled(True)
        self.matches_table = QTableWidget(0, 0)
        self.matches_table.setSortingEnabled(True)
        self.table_tabs.addTab(self.summary_table, "Summary")
        self.table_tabs.addTab(self.matches_table, "Matched genes")

        self.status = QLabel("")
        self.status.setWordWrap(True)
        self._tables: dict[str, pd.DataFrame] = {"summary": pd.DataFrame(), "matched_genes": pd.DataFrame()}
        self._thread: QThread | None = None
        self._worker: RbpSupportWorker | None = None
        self.cross_search = QLineEdit()
        self.cross_search.setPlaceholderText("Search RBP set / gene / support pattern...")
        self.cross_search.textChanged.connect(self._apply_cross_filters)
        self.cross_top_n = QSpinBox()
        self.cross_top_n.setRange(1, 5000)
        self.cross_top_n.setValue(50)
        self.cross_custom_path = QLineEdit()
        self.cross_custom_path.setPlaceholderText("Optional custom RBP gene set TSV/CSV")
        self.cross_choose_button = QPushButton("Choose gene set")
        self.cross_choose_button.clicked.connect(self._choose_cross_gene_set)
        self.cross_run_button = QPushButton("Run Cross-comparison RBP support")
        self.cross_run_button.clicked.connect(self._run_cross)
        self.cross_open_button = QPushButton("Open output folder")
        self.cross_open_button.clicked.connect(self._open_cross_output_folder)
        self.cross_export_button = QPushButton("Export current table")
        self.cross_export_button.clicked.connect(self._export_cross_table)
        self.cross_progress = QProgressBar()
        self.cross_progress.setRange(0, 0)
        self.cross_progress.setVisible(False)
        self.cross_progress.setMaximumWidth(170)
        self.cross_table_tabs = QTabWidget()
        self.cross_matrix_table = QTableWidget(0, 0)
        self.cross_matrix_table.setSortingEnabled(True)
        self.cross_matches_table = QTableWidget(0, 0)
        self.cross_matches_table.setSortingEnabled(True)
        self.cross_comparison_summary_table = QTableWidget(0, 0)
        self.cross_comparison_summary_table.setSortingEnabled(True)
        self.cross_plot_table = QTableWidget(0, 2)
        self.cross_plot_table.setHorizontalHeaderLabels(["Plot", "Path"])
        self.cross_plot_table.setSortingEnabled(True)
        self.cross_plot_table.currentCellChanged.connect(self._preview_cross_plot)
        self.cross_preview_image = QLabel("No cached cross-comparison plot result.")
        self.cross_preview_image.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.cross_preview_image.setMinimumHeight(320)
        self.cross_preview_scroll = QScrollArea()
        self.cross_preview_scroll.setWidgetResizable(True)
        self.cross_preview_scroll.setWidget(self.cross_preview_image)
        self.cross_export_plot_button = QPushButton("Export selected plot")
        self.cross_export_plot_button.clicked.connect(self._export_cross_plot)
        self.cross_table_tabs.addTab(self.cross_matrix_table, "RBP enrichment matrix")
        self.cross_table_tabs.addTab(self.cross_matches_table, "Matched gene matrix")
        self.cross_table_tabs.addTab(self.cross_comparison_summary_table, "Comparison summary")
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
            "enrichment_matrix": pd.DataFrame(),
            "matched_gene_matrix": pd.DataFrame(),
            "comparison_summary": pd.DataFrame(),
        }
        self._cross_plot_frame = pd.DataFrame()
        self._cross_preview_pixmap = QPixmap()
        self._cross_thread: QThread | None = None
        self._cross_worker: RbpCrossSupportWorker | None = None

        controls = QHBoxLayout()
        controls.addWidget(QLabel("Top N"))
        controls.addWidget(self.top_n)
        controls.addWidget(self.search, 1)
        controls.addWidget(self.choose_button)
        controls.addWidget(self.run_button)
        controls.addWidget(self.progress)
        controls.addWidget(self.open_button)
        controls.addWidget(self.export_button)

        path_row = QHBoxLayout()
        path_row.addWidget(self.source_label)
        path_row.addWidget(self.custom_path, 1)

        per_tab = QWidget()
        per_layout = QVBoxLayout(per_tab)
        per_layout.addWidget(self.context)
        per_layout.addWidget(self.comparison_tabs)
        per_layout.addLayout(path_row)
        per_layout.addLayout(controls)
        per_layout.addWidget(self.table_tabs, 1)
        per_layout.addWidget(self.status)

        cross_context = QLabel(
            "Cross-comparison RBP support compares per-comparison RBP/splicing-factor overlap and enrichment patterns. "
            "It is a support-prioritization view, not causal RBP regulation."
        )
        cross_context.setWordWrap(True)
        cross_path_row = QHBoxLayout()
        cross_path_row.addWidget(QLabel("Source: candidate selection, fallback to cached top N if no selection"))
        cross_path_row.addWidget(self.cross_custom_path, 1)
        cross_controls = QHBoxLayout()
        cross_controls.addWidget(QLabel("Top N"))
        cross_controls.addWidget(self.cross_top_n)
        cross_controls.addWidget(self.cross_search, 1)
        cross_controls.addWidget(self.cross_choose_button)
        cross_controls.addWidget(self.cross_run_button)
        cross_controls.addWidget(self.cross_progress)
        cross_controls.addWidget(self.cross_open_button)
        cross_controls.addWidget(self.cross_export_button)
        cross_tab = QWidget()
        cross_layout = QVBoxLayout(cross_tab)
        cross_layout.addWidget(cross_context)
        cross_layout.addLayout(cross_path_row)
        cross_layout.addLayout(cross_controls)
        cross_layout.addWidget(self.cross_table_tabs, 1)
        cross_layout.addWidget(self.cross_status)

        self.main_tabs = QTabWidget()
        self.main_tabs.addTab(per_tab, "7.1 Per-comparison RBP support")
        self.main_tabs.addTab(cross_tab, "7.2 Cross-comparison RBP support")

        layout = QVBoxLayout(self)
        layout.addWidget(self.title)
        layout.addWidget(self.main_tabs)

    def refresh(self) -> None:
        project = self.project_service.current_project
        if project is None:
            self._set_tables({"summary": pd.DataFrame(), "matched_genes": pd.DataFrame()})
            self.status.setText("No project loaded.")
            return
        self._populate_comparison_tabs()
        comparison_id = self._current_comparison_id()
        if not comparison_id:
            self.status.setText("No comparison selected.")
            return
        tables = self.project_service.preview_rbp_support(
            comparison_id,
            top_n=self.top_n.value(),
            custom_gene_set_path=self.custom_path.text().strip() or None,
            allow_generate=False,
            force_rebuild=False,
        )
        self._set_tables(tables)
        state = self.project_service.module_state("rbp_support")
        self.status.setText(f"{state.get('status', 'not_run')}: {state.get('message', '')}")
        self._refresh_cross_cached()

    def _populate_comparison_tabs(self) -> None:
        comparisons = self.project_service.selected_comparisons_for_display()
        labels = [self.comparison_tabs.tabText(i) for i in range(self.comparison_tabs.count())]
        wanted = [comparison.display_resolved_name for comparison in comparisons]
        if labels == wanted:
            return
        current_id = self._current_comparison_id()
        self.comparison_tabs.blockSignals(True)
        self.comparison_tabs.clear()
        for comparison in comparisons:
            tab = QWidget()
            tab.setProperty("comparison_id", comparison.comparison_id)
            self.comparison_tabs.addTab(tab, comparison.display_resolved_name)
        if current_id:
            for index, comparison in enumerate(comparisons):
                if comparison.comparison_id == current_id:
                    self.comparison_tabs.setCurrentIndex(index)
                    break
        self.comparison_tabs.blockSignals(False)

    def _current_comparison_id(self) -> str:
        widget = self.comparison_tabs.currentWidget()
        return str(widget.property("comparison_id") or "") if widget is not None else ""

    def _choose_gene_set(self) -> None:
        project = self.project_service.current_project
        start = str(project.project_root) if project is not None else ""
        path, _ = QFileDialog.getOpenFileName(self, "Choose RBP gene set", start, "Tables (*.tsv *.csv *.txt);;All Files (*)")
        if path:
            self.custom_path.setText(path)
            self.refresh()

    def _choose_cross_gene_set(self) -> None:
        project = self.project_service.current_project
        start = str(project.project_root) if project is not None else ""
        path, _ = QFileDialog.getOpenFileName(self, "Choose RBP gene set", start, "Tables (*.tsv *.csv *.txt);;All Files (*)")
        if path:
            self.cross_custom_path.setText(path)
            self._refresh_cross_cached()

    def _run(self) -> None:
        comparison_id = self._current_comparison_id()
        if not comparison_id:
            self.status.setText("No comparison selected.")
            return
        self.progress.setVisible(True)
        self.run_button.setEnabled(False)
        self._thread = QThread(self)
        self._worker = RbpSupportWorker(
            self.project_service,
            comparison_id,
            "selection",
            self.top_n.value(),
            self.custom_path.text().strip(),
        )
        self._worker.moveToThread(self._thread)
        self._thread.started.connect(self._worker.run)
        self._worker.finished.connect(self._on_run_finished)
        self._worker.failed.connect(self._on_run_failed)
        self._worker.finished.connect(self._thread.quit)
        self._worker.failed.connect(self._thread.quit)
        self._thread.finished.connect(self._thread.deleteLater)
        self._thread.start()

    def _on_run_finished(self, tables: object) -> None:
        self.progress.setVisible(False)
        self.run_button.setEnabled(True)
        self._set_tables(tables if isinstance(tables, dict) else {})
        state = self.project_service.module_state("rbp_support")
        self.status.setText(f"{state.get('status', 'finished')}: {state.get('message', '')}")

    def _on_run_failed(self, message: str) -> None:
        self.progress.setVisible(False)
        self.run_button.setEnabled(True)
        self.status.setText(f"RBP support failed: {message}")

    def _refresh_cross_cached(self) -> None:
        if self.project_service.current_project is None:
            self._set_cross_tables({})
            self._set_cross_plot_files(pd.DataFrame())
            self.cross_status.setText("No project loaded.")
            return
        tables = self.project_service.preview_rbp_cross_comparison_support(
            source="selection",
            top_n=self.cross_top_n.value(),
            custom_gene_set_path=self.cross_custom_path.text().strip() or None,
            allow_generate=False,
            force_rebuild=False,
        )
        self._set_cross_tables(tables)
        self._refresh_cross_plots()
        state = self.project_service.module_state("rbp_cross_support")
        self.cross_status.setText(f"{state.get('status', 'not_run')}: {state.get('message', '')}")

    def _run_cross(self) -> None:
        if self.project_service.current_project is None:
            self.cross_status.setText("No project loaded.")
            return
        self.cross_progress.setVisible(True)
        self.cross_run_button.setEnabled(False)
        self._cross_thread = QThread(self)
        self._cross_worker = RbpCrossSupportWorker(
            self.project_service,
            self.cross_top_n.value(),
            self.cross_custom_path.text().strip(),
        )
        self._cross_worker.moveToThread(self._cross_thread)
        self._cross_thread.started.connect(self._cross_worker.run)
        self._cross_worker.finished.connect(self._on_cross_finished)
        self._cross_worker.failed.connect(self._on_cross_failed)
        self._cross_worker.finished.connect(self._cross_thread.quit)
        self._cross_worker.failed.connect(self._cross_thread.quit)
        self._cross_thread.finished.connect(self._cross_thread.deleteLater)
        self._cross_thread.start()

    def _on_cross_finished(self, tables: object) -> None:
        self.cross_progress.setVisible(False)
        self.cross_run_button.setEnabled(True)
        self._set_cross_tables(tables if isinstance(tables, dict) else {})
        self._refresh_cross_plots()
        state = self.project_service.module_state("rbp_cross_support")
        self.cross_status.setText(f"{state.get('status', 'finished')}: {state.get('message', '')}")

    def _on_cross_failed(self, message: str) -> None:
        self.cross_progress.setVisible(False)
        self.cross_run_button.setEnabled(True)
        self.cross_status.setText(f"Cross-comparison RBP support failed: {message}")

    def _set_tables(self, tables: dict[str, pd.DataFrame]) -> None:
        self._tables = {
            "summary": tables.get("summary", pd.DataFrame()).copy() if isinstance(tables, dict) else pd.DataFrame(),
            "matched_genes": tables.get("matched_genes", pd.DataFrame()).copy() if isinstance(tables, dict) else pd.DataFrame(),
        }
        self._apply_filters()

    def _apply_filters(self) -> None:
        query = self.search.text().strip().casefold()
        for key, table in (("summary", self.summary_table), ("matched_genes", self.matches_table)):
            frame = self._tables.get(key, pd.DataFrame()).copy()
            if query and not frame.empty:
                mask = pd.Series(False, index=frame.index)
                for column in frame.columns:
                    mask |= frame[column].fillna("").astype(str).str.casefold().str.contains(query, regex=False)
                frame = frame.loc[mask].copy()
            self._populate_table(table, frame)

    def _set_cross_tables(self, tables: dict[str, pd.DataFrame]) -> None:
        self._cross_tables = {
            "enrichment_matrix": tables.get("enrichment_matrix", pd.DataFrame()).copy() if isinstance(tables, dict) else pd.DataFrame(),
            "matched_gene_matrix": tables.get("matched_gene_matrix", pd.DataFrame()).copy() if isinstance(tables, dict) else pd.DataFrame(),
            "comparison_summary": tables.get("comparison_summary", pd.DataFrame()).copy() if isinstance(tables, dict) else pd.DataFrame(),
        }
        self._apply_cross_filters()

    def _refresh_cross_plots(self) -> None:
        if self.project_service.current_project is None:
            self._set_cross_plot_files(pd.DataFrame())
            return
        self._set_cross_plot_files(self.project_service.rbp_cross_support_plot_files())

    def _set_cross_plot_files(self, frame: pd.DataFrame) -> None:
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
        self._preview_cross_plot(0, 0, -1, -1)

    def _selected_cross_plot_path(self) -> str:
        row = self.cross_plot_table.currentRow()
        if row < 0:
            return ""
        item = self.cross_plot_table.item(row, 1)
        return item.text() if item is not None else ""

    def _preview_cross_plot(self, _current_row: int, _current_column: int, _previous_row: int, _previous_column: int) -> None:
        path = self._selected_cross_plot_path()
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
        self._fit_cross_preview_pixmap()

    def _fit_cross_preview_pixmap(self) -> None:
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
        QTimer.singleShot(0, self._fit_cross_preview_pixmap)

    def _apply_cross_filters(self) -> None:
        query = self.cross_search.text().strip().casefold()
        mapping = (
            ("enrichment_matrix", self.cross_matrix_table),
            ("matched_gene_matrix", self.cross_matches_table),
            ("comparison_summary", self.cross_comparison_summary_table),
        )
        for key, table in mapping:
            frame = self._cross_tables.get(key, pd.DataFrame()).copy()
            if query and not frame.empty:
                mask = pd.Series(False, index=frame.index)
                for column in frame.columns:
                    mask |= frame[column].fillna("").astype(str).str.casefold().str.contains(query, regex=False)
                frame = frame.loc[mask].copy()
            self._populate_table(table, frame)

    def _open_output_folder(self) -> None:
        comparison_id = self._current_comparison_id()
        folder = self.project_service.rbp_support_output_dir(comparison_id)
        if folder is None:
            return
        folder.mkdir(parents=True, exist_ok=True)
        os.startfile(str(folder))

    def _open_cross_output_folder(self) -> None:
        folder = self.project_service.rbp_cross_support_output_dir()
        if folder is None:
            return
        folder.mkdir(parents=True, exist_ok=True)
        os.startfile(str(folder))

    def _export_current_table(self) -> None:
        key = "summary" if self.table_tabs.currentIndex() == 0 else "matched_genes"
        frame = self._tables.get(key, pd.DataFrame())
        if frame.empty:
            return
        project = self.project_service.current_project
        default = (project.project_root if project else None) or None
        path, _ = QFileDialog.getSaveFileName(self, "Export RBP table", str(default / f"{key}.tsv") if default else f"{key}.tsv", "TSV Files (*.tsv);;All Files (*)")
        if path:
            frame.to_csv(path, sep="\t", index=False)

    def _export_cross_table(self) -> None:
        keys = ["enrichment_matrix", "matched_gene_matrix", "comparison_summary"]
        if self.cross_table_tabs.currentIndex() >= len(keys):
            self._export_cross_plot()
            return
        key = keys[self.cross_table_tabs.currentIndex()]
        frame = self._cross_tables.get(key, pd.DataFrame())
        if frame.empty:
            return
        project = self.project_service.current_project
        default = (project.project_root if project else None) or None
        path, _ = QFileDialog.getSaveFileName(self, "Export cross-comparison RBP table", str(default / f"{key}.tsv") if default else f"{key}.tsv", "TSV Files (*.tsv);;All Files (*)")
        if path:
            frame.to_csv(path, sep="\t", index=False)

    def _export_cross_plot(self) -> None:
        source = self._selected_cross_plot_path()
        if not source or not os.path.exists(source):
            return
        project = self.project_service.current_project
        start_dir = project.project_root if project is not None else None
        default_name = os.path.basename(source)
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Export selected RBP plot",
            str(start_dir / default_name) if start_dir is not None else default_name,
            "PNG Files (*.png);;All Files (*)",
        )
        if path:
            shutil.copy2(source, path)

    @staticmethod
    def _populate_table(table: QTableWidget, frame: pd.DataFrame) -> None:
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
