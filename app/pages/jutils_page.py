from __future__ import annotations

import os
import shutil
from pathlib import Path

import pandas as pd
from PySide6.QtCore import Qt, QUrl
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import (
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QScrollArea,
    QSplitter,
    QStackedWidget,
    QTableWidget,
    QTableWidgetItem,
    QTextBrowser,
    QTextEdit,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)
from PySide6.QtSvgWidgets import QSvgWidget

from src.services.project_service import ProjectService


class JutilsPage(QWidget):
    def __init__(self, project_service: ProjectService) -> None:
        super().__init__()
        self.project_service = project_service
        self.label = QLabel("Jutils output browser")
        self.status_label = QLabel("")
        self.status_label.setWordWrap(True)
        self.comparison_tabs = QTabWidget()
        self.comparison_tabs.currentChanged.connect(self.refresh)
        self.command_label = QTextEdit()
        self.command_label.setReadOnly(True)
        self.command_label.setMaximumHeight(90)
        self.output_dir_label = QLabel("")
        self.output_dir_label.setWordWrap(True)

        self.run_button = QPushButton("Run Jutils")
        self.run_button.clicked.connect(self._run_jutils)
        self.preview_button = QPushButton("Preview Selected")
        self.preview_button.clicked.connect(self._preview_selected)
        self.open_output_button = QPushButton("Open Output Folder")
        self.open_output_button.clicked.connect(self._open_output_folder)
        self.open_selected_button = QPushButton("Open Selected File")
        self.open_selected_button.clicked.connect(self._open_selected_file)
        self.download_button = QPushButton("Download Selected")
        self.download_button.clicked.connect(self._download_selected_file)

        self.table = QTableWidget(0, 0)
        self.table.setSortingEnabled(True)
        self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.table.itemSelectionChanged.connect(self._preview_selected)
        self._frame = pd.DataFrame()

        self.preview_stack = QStackedWidget()
        self.image_label = QLabel("No preview selected.")
        self.image_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.image_label.setScaledContents(False)
        self.image_scroll = QScrollArea()
        self.image_scroll.setWidgetResizable(True)
        self.image_scroll.setWidget(self.image_label)
        self.svg_widget = QSvgWidget()
        self.svg_widget.setMinimumHeight(420)
        self.html_browser = QTextBrowser()
        self.table_preview = QTableWidget(0, 0)
        self.table_preview.setSortingEnabled(True)
        self.text_preview = QTextEdit()
        self.text_preview.setReadOnly(True)
        self.preview_stack.addWidget(self.image_scroll)
        self.preview_stack.addWidget(self.svg_widget)
        self.preview_stack.addWidget(self.html_browser)
        self.preview_stack.addWidget(self.table_preview)
        self.preview_stack.addWidget(self.text_preview)

        self.details = QTextEdit()
        self.details.setReadOnly(True)

        top_buttons = QHBoxLayout()
        top_buttons.addWidget(self.run_button)
        top_buttons.addWidget(self.preview_button)
        top_buttons.addWidget(self.open_output_button)
        top_buttons.addWidget(self.open_selected_button)
        top_buttons.addWidget(self.download_button)
        top_buttons.addStretch(1)

        meta_layout = QVBoxLayout()
        meta_layout.addWidget(self.label)
        meta_layout.addLayout(top_buttons)
        meta_layout.addWidget(self.comparison_tabs)
        meta_layout.addWidget(self.status_label)
        meta_layout.addWidget(self.output_dir_label)
        meta_layout.addWidget(QLabel("Command log"))
        meta_layout.addWidget(self.command_label)

        splitter = QSplitter(Qt.Orientation.Vertical)
        splitter.addWidget(self.table)
        splitter.addWidget(self.preview_stack)
        splitter.addWidget(self.details)
        splitter.setStretchFactor(0, 2)
        splitter.setStretchFactor(1, 3)
        splitter.setStretchFactor(2, 1)
        splitter.setSizes([220, 520, 180])

        layout = QVBoxLayout(self)
        layout.addLayout(meta_layout)
        layout.addWidget(splitter, 1)

    def refresh(self) -> None:
        project = self.project_service.current_project
        if project is None:
            self.table.setRowCount(0)
            self.table.setColumnCount(0)
            self.status_label.setText("No project loaded.")
            self.comparison_tabs.clear()
            self.output_dir_label.clear()
            self.command_label.clear()
            self.details.clear()
            self._frame = pd.DataFrame()
            return

        report = self.project_service.jutils_report()
        browser = self.project_service.build_jutils_output_browser()
        self._populate_comparison_tabs(browser)
        selected_comparison = self._selected_comparison_id()
        if selected_comparison not in (None, "__all__") and not browser.empty and "comparison_id" in browser.columns:
            browser = browser.loc[browser["comparison_id"].astype(str) == str(selected_comparison)].copy()
        self._frame = browser
        self._populate_file_table(browser)
        status = str(report.get("status", "idle"))
        output_dir = str(report.get("output_directory", ""))
        self.status_label.setText(f"Jutils status: {status}")
        self.output_dir_label.setText(f"Output directory: {output_dir or 'not available'}")
        commands = report.get("command_lines", []) or []
        self.command_label.setPlainText("\n".join(str(line) for line in commands) or "No command has been run yet.")

        if status == "finished":
            generated_files = report.get("generated_files", []) or []
            lines = [
                "Jutils finished.",
                f"Output saved to: {output_dir or 'not available'}",
                "Generated files:",
            ]
            lines.extend(f"- {row.get('relative_path', '')}" for row in generated_files[:200])
            self.details.setPlainText("\n".join(lines))
        elif status == "failed":
            self.details.setPlainText(
                "\n".join(
                    [
                        "Jutils failed.",
                        "Command:",
                        "\n".join(str(line) for line in commands) or "not available",
                        "",
                        "stdout:",
                        str(report.get("stdout", "") or ""),
                        "",
                        "stderr:",
                        str(report.get("stderr", "") or ""),
                        "",
                        "error log:",
                        str(report.get("error_log", "") or ""),
                        "",
                        f"output directory: {output_dir or 'not available'}",
                    ]
                )
            )
        else:
            self.details.setPlainText("Run Jutils to populate output files and previews.")
        self._preview_selected()

    def _populate_comparison_tabs(self, frame: pd.DataFrame) -> None:
        current = self._selected_comparison_id()
        comparisons = []
        if "comparison_id" in frame.columns:
            comparisons = sorted(frame["comparison_id"].dropna().astype(str).unique().tolist())
        elif self.project_service.current_project is not None:
            comparisons = [item.comparison_id for item in self.project_service.selected_comparisons_for_display()]
        self.comparison_tabs.blockSignals(True)
        self.comparison_tabs.clear()
        self.comparison_tabs.addTab(QWidget(), "All comparisons")
        self.comparison_tabs.setTabToolTip(0, "__all__")
        for index, comparison_id in enumerate(comparisons, start=1):
            self.comparison_tabs.addTab(QWidget(), comparison_id)
            self.comparison_tabs.setTabToolTip(index, comparison_id)
        if current:
            for index in range(self.comparison_tabs.count()):
                if self.comparison_tabs.tabToolTip(index) == current:
                    self.comparison_tabs.setCurrentIndex(index)
                    break
        self.comparison_tabs.blockSignals(False)

    def _selected_comparison_id(self):
        if self.comparison_tabs.count() > 0:
            tooltip = self.comparison_tabs.tabToolTip(self.comparison_tabs.currentIndex())
            if tooltip:
                return tooltip
        return "__all__"

    def _populate_file_table(self, frame: pd.DataFrame) -> None:
        if frame.empty:
            self.table.setRowCount(0)
            self.table.setColumnCount(0)
            return
        columns = ["comparison_id", "relative_path", "file_type", "file_size", "last_modified", "preview_kind"]
        self.table.setColumnCount(len(columns))
        self.table.setHorizontalHeaderLabels(columns)
        self.table.setRowCount(len(frame))
        for row_idx, (_, row) in enumerate(frame[columns].iterrows()):
            for col_idx, column in enumerate(columns):
                item = QTableWidgetItem("" if pd.isna(row[column]) else str(row[column]))
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(row_idx, col_idx, item)

    def _run_jutils(self) -> None:
        self.status_label.setText("Jutils status: running")
        self.details.setPlainText("Jutils is running...")
        try:
            self.project_service.run_jutils_pipeline()
        except Exception:
            pass
        self.refresh()

    def _selected_row(self) -> pd.Series | None:
        if self._frame.empty:
            return None
        row_idx = self.table.currentRow()
        if row_idx < 0 or row_idx >= len(self._frame):
            return None
        return self._frame.iloc[row_idx]

    def _preview_selected(self) -> None:
        row = self._selected_row()
        if row is None:
            self.preview_stack.setCurrentWidget(self.text_preview)
            self.text_preview.setPlainText("Select a generated Jutils file to preview it.")
            return
        path = Path(str(row.get("absolute_path", "")))
        preview_kind = str(row.get("preview_kind", "binary"))
        if not path.exists():
            self.preview_stack.setCurrentWidget(self.text_preview)
            self.text_preview.setPlainText(f"Missing file:\n{path}")
            return
        if preview_kind == "image":
            if path.suffix.lower() == ".svg":
                self.svg_widget.load(str(path))
                self.preview_stack.setCurrentWidget(self.svg_widget)
            else:
                pixmap = QPixmap(str(path))
                self.image_label.setPixmap(pixmap)
                self.preview_stack.setCurrentWidget(self.image_scroll)
        elif preview_kind == "html":
            self.html_browser.setSource(QUrl.fromLocalFile(str(path)))
            self.preview_stack.setCurrentWidget(self.html_browser)
        elif preview_kind in {"table", "excel"}:
            self._load_table_preview(path)
            self.preview_stack.setCurrentWidget(self.table_preview)
        elif preview_kind == "pdf":
            self.preview_stack.setCurrentWidget(self.text_preview)
            self.text_preview.setPlainText(f"PDF preview is not embedded.\nOpen or download:\n{path}")
        else:
            self.preview_stack.setCurrentWidget(self.text_preview)
            self.text_preview.setPlainText(
                "\n".join(
                    [
                        f"File: {path}",
                        f"Type: {row.get('file_type', '')}",
                        f"Size: {row.get('file_size', '')}",
                        f"Modified: {row.get('last_modified', '')}",
                    ]
                )
            )

    def _load_table_preview(self, path: Path) -> None:
        try:
            if path.suffix.lower() in {".xlsx", ".xls"}:
                frame = pd.read_excel(path).head(500)
            elif path.suffix.lower() == ".csv":
                frame = pd.read_csv(path, low_memory=False).head(500)
            else:
                frame = pd.read_csv(path, sep="\t", low_memory=False).head(500)
        except Exception as exc:
            self.preview_stack.setCurrentWidget(self.text_preview)
            self.text_preview.setPlainText(f"Failed to preview table:\n{path}\n\n{exc}")
            return
        self.table_preview.setColumnCount(len(frame.columns))
        self.table_preview.setHorizontalHeaderLabels(frame.columns.astype(str).tolist())
        self.table_preview.setRowCount(len(frame))
        for row_idx, (_, row) in enumerate(frame.iterrows()):
            for col_idx, column in enumerate(frame.columns):
                item = QTableWidgetItem("" if pd.isna(row[column]) else str(row[column]))
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self.table_preview.setItem(row_idx, col_idx, item)

    def _open_output_folder(self) -> None:
        report = self.project_service.jutils_report()
        output_dir = Path(str(report.get("output_directory", ""))) if report.get("output_directory") else None
        if output_dir and output_dir.exists():
            os.startfile(str(output_dir))

    def _open_selected_file(self) -> None:
        row = self._selected_row()
        if row is None:
            return
        path = Path(str(row.get("absolute_path", "")))
        if path.exists():
            os.startfile(str(path))

    def _download_selected_file(self) -> None:
        row = self._selected_row()
        if row is None:
            return
        path = Path(str(row.get("absolute_path", "")))
        if not path.exists():
            return
        output_path, _ = QFileDialog.getSaveFileName(
            self,
            "Download Selected Jutils File",
            path.name,
            "All Files (*.*)",
        )
        if not output_path:
            return
        shutil.copy2(path, output_path)
