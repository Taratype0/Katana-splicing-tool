from __future__ import annotations

import os

import pandas as pd
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QDoubleSpinBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from src.services.project_service import ProjectService


class SashimiPage(QWidget):
    def __init__(self, project_service: ProjectService, view_mode: str = "run") -> None:
        super().__init__()
        self.project_service = project_service
        self.view_mode = view_mode

        self.label = QLabel("Run Sashimi")
        self.explanation = QLabel("")
        self.explanation.setWordWrap(True)

        self.comparison_tabs = QTabWidget()
        self.comparison_tabs.setDocumentMode(True)
        self.comparison_tabs.setUsesScrollButtons(True)
        self.comparison_tabs.setMaximumHeight(36)
        self.comparison_tabs.currentChanged.connect(self.refresh)

        self.search = QLineEdit()
        self.search.setPlaceholderText("Search gene / event_id / event_type / coordinates")
        self.search.textChanged.connect(self.refresh)

        self.event_type_filter = QComboBox()
        self.event_type_filter.currentTextChanged.connect(self.refresh)

        self.fdr_filter = QDoubleSpinBox()
        self.fdr_filter.setRange(0.0, 1.0)
        self.fdr_filter.setSingleStep(0.01)
        self.fdr_filter.setValue(0.05)
        self.fdr_filter.setPrefix("FDR <= ")
        self.fdr_filter.valueChanged.connect(self.refresh)

        self.dpsi_filter = QDoubleSpinBox()
        self.dpsi_filter.setRange(0.0, 5.0)
        self.dpsi_filter.setSingleStep(0.05)
        self.dpsi_filter.setValue(0.10)
        self.dpsi_filter.setPrefix("|dPSI| >= ")
        self.dpsi_filter.valueChanged.connect(self.refresh)

        self.generate_manifest_button = QPushButton("Generate selected\nmanifest")
        self.generate_manifest_button.clicked.connect(self._generate_manifest)
        self.run_button = QPushButton("Run selected\nsashimi")
        self.run_button.clicked.connect(self._run_sashimi)
        self.open_output_button = QPushButton("Open output\nfolder")
        self.open_output_button.clicked.connect(self._open_output_folder)
        for button in (self.generate_manifest_button, self.run_button, self.open_output_button):
            button.setMinimumHeight(40)
            button.setMinimumWidth(130)

        self.event_table = QTableWidget(0, 0)
        self.event_table.setSortingEnabled(True)
        self.event_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.event_table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.event_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.event_table.itemSelectionChanged.connect(self._show_details)

        self.manifest_table = QTableWidget(0, 0)
        self.manifest_table.setSortingEnabled(True)
        self.manifest_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.manifest_table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.manifest_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.manifest_table.itemSelectionChanged.connect(self._show_details)

        self.output_table = QTableWidget(0, 0)
        self.output_table.setSortingEnabled(True)
        self.output_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.output_table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.output_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.output_table.itemSelectionChanged.connect(self._show_details)

        self.failed_table = QTableWidget(0, 0)
        self.failed_table.setSortingEnabled(True)
        self.failed_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.failed_table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.failed_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.failed_table.itemSelectionChanged.connect(self._show_details)

        self.events_label = QLabel("Available splicing events")
        self.manifest_label = QLabel("Generated sashimi manifest rows")
        self.output_label = QLabel("Generated sashimi output files")
        self.failed_label = QLabel("Failed sashimi jobs")

        self.status_label = QLabel("")
        self.status_label.setWordWrap(True)
        self.status_label.setStyleSheet("font-size: 11px; color: #C9CDD3;")

        self.details = QTextEdit()
        self.details.setReadOnly(True)

        filters = QHBoxLayout()
        filters.addWidget(self.search, 2)
        filters.addWidget(self.event_type_filter)
        filters.addWidget(self.fdr_filter)
        filters.addWidget(self.dpsi_filter)
        filters.addWidget(self.generate_manifest_button)
        filters.addWidget(self.run_button)
        filters.addWidget(self.open_output_button)

        layout = QVBoxLayout(self)
        layout.addWidget(self.label)
        layout.addWidget(self.explanation)
        layout.addWidget(self.comparison_tabs)
        layout.addLayout(filters)
        layout.addWidget(self.events_label)
        layout.addWidget(self.event_table, 3)
        layout.addWidget(self.manifest_label)
        layout.addWidget(self.manifest_table, 2)
        layout.addWidget(self.output_label)
        layout.addWidget(self.output_table, 2)
        layout.addWidget(self.failed_label)
        layout.addWidget(self.failed_table, 1)
        layout.addWidget(self.status_label)
        layout.addWidget(self.details, 1)

        self._comparison_ids: list[str] = []
        self._events = pd.DataFrame()
        self._manifest = pd.DataFrame()
        self._outputs = pd.DataFrame()
        self._failures = pd.DataFrame()
        self._apply_mode_visibility()

    def refresh(self) -> None:
        project = self.project_service.current_project
        if project is None:
            self._comparison_ids = []
            self.comparison_tabs.clear()
            self._events = pd.DataFrame()
            self._manifest = pd.DataFrame()
            self._outputs = pd.DataFrame()
            self._failures = pd.DataFrame()
            self._populate_table(self.event_table, pd.DataFrame())
            self._populate_table(self.manifest_table, pd.DataFrame())
            self._populate_table(self.output_table, pd.DataFrame())
            self._populate_table(self.failed_table, pd.DataFrame())
            self.status_label.setText("Load a project first.")
            self.details.setPlainText("")
            return

        self._refresh_comparison_tabs()
        comparison_id = self._current_comparison_id()

        events = self.project_service.available_sashimi_events(comparison_id, allow_generate=False)
        events = self._apply_filters(events)
        self._events = events.copy()
        self._populate_event_type_filter(events)
        self._populate_table(
            self.event_table,
            events[
                [
                    column
                    for column in [
                        "gene_symbol",
                        "event_type",
                        "event_id",
                        "comparison_display_name",
                        "dPSI",
                        "FDR",
                        "direction",
                        "coordinates",
                    ]
                    if column in events.columns
                ]
            ].copy()
            if not events.empty
            else pd.DataFrame(),
        )

        manifest = self.project_service.run_state.results.sashimi_manifest.copy()
        if not manifest.empty and comparison_id and "comparison_id" in manifest.columns:
            manifest = manifest.loc[manifest["comparison_id"].astype(str) == comparison_id].copy()
        self._manifest = manifest.copy()
        self._populate_table(
            self.manifest_table,
            manifest[
                [
                    column
                    for column in [
                        "geneSymbol",
                        "event_type",
                        "event_uid",
                        "comparison_name",
                        "label1",
                        "label2",
                        "event_file",
                        "outdir",
                    ]
                    if column in manifest.columns
                ]
            ].copy()
            if not manifest.empty
            else pd.DataFrame(),
        )

        outputs = self.project_service.build_sashimi_output_browser(comparison_id)
        self._outputs = outputs.copy()
        self._populate_table(
            self.output_table,
            outputs[
                [
                    column
                    for column in [
                        "comparison_id",
                        "relative_path",
                        "file_type",
                        "file_size",
                        "last_modified",
                        "preview_kind",
                    ]
                    if column in outputs.columns
                ]
            ].copy()
            if not outputs.empty
            else pd.DataFrame(),
        )

        failures = self.project_service.last_sashimi_failures()
        if not failures.empty and comparison_id and "comparison_id" in failures.columns:
            failures = failures.loc[failures["comparison_id"].astype(str) == comparison_id].copy()
        self._failures = failures.copy()
        self._populate_table(
            self.failed_table,
            failures[
                [
                    column
                    for column in [
                        "gene",
                        "event_type",
                        "event_id",
                        "comparison_id",
                        "error_message",
                        "output_folder",
                    ]
                    if column in failures.columns
                ]
            ].copy()
            if not failures.empty
            else pd.DataFrame(),
        )

        manifest_state = self.project_service.module_state("sashimi_manifest")
        plot_state = self.project_service.module_state("sashimi_plot")
        if self.view_mode == "run":
            if events.empty:
                self.status_label.setText(
                    "No available splicing events found for this comparison. "
                    "Missing input: candidate event follow-up source not available. Please run 5.3.2 / 5.3.4 first."
                )
            elif plot_state.get("status") in {"finished", "failed"}:
                self.status_label.setText(str(plot_state.get("message", "")))
            elif manifest_state.get("status") in {"finished", "failed"}:
                self.status_label.setText(str(manifest_state.get("message", "")))
            else:
                self.status_label.setText(
                    f"Available events: {len(events)} | Generated manifest rows: {len(manifest)}. "
                    "Select one or more events, then click 'Generate selected manifest' or 'Run selected sashimi'."
                )
        elif self.view_mode == "preview":
            if manifest_state.get("status") == "failed":
                self.status_label.setText(str(manifest_state.get("message", "")))
            elif plot_state.get("status") == "failed":
                self.status_label.setText(str(plot_state.get("message", "")))
            elif outputs.empty:
                self.status_label.setText(
                    f"Generated manifest rows: {len(manifest)} | Generated sashimi output files: 0. "
                    "No generated sashimi plot files were found yet. Run selected events from 5.5.2 first."
                )
            else:
                self.status_label.setText(
                    f"Generated manifest rows: {len(manifest)} | Generated sashimi output files: {len(outputs)}."
                )
        else:
            self.status_label.setText(f"Failed sashimi jobs: {len(failures)}.")
        self._show_details()

    def _apply_mode_visibility(self) -> None:
        run_mode = self.view_mode == "run"
        preview_mode = self.view_mode == "preview"
        failed_mode = self.view_mode == "failed"

        self.search.setVisible(run_mode)
        self.event_type_filter.setVisible(run_mode)
        self.fdr_filter.setVisible(run_mode)
        self.dpsi_filter.setVisible(run_mode)
        self.generate_manifest_button.setVisible(run_mode)
        self.run_button.setVisible(run_mode)

        self.events_label.setVisible(run_mode)
        self.event_table.setVisible(run_mode)
        self.manifest_label.setVisible(run_mode or preview_mode)
        self.manifest_table.setVisible(run_mode or preview_mode)
        self.output_label.setVisible(preview_mode)
        self.output_table.setVisible(preview_mode)
        self.failed_label.setVisible(failed_mode)
        self.failed_table.setVisible(failed_mode)

        if run_mode:
            self.label.setText("5.5.2 Run Sashimi")
            self.explanation.setText(
                "This page only lists available splicing events for the current comparison. Select events here, generate a manifest, then run rmats2sashimi only for the selected events."
            )
        elif preview_mode:
            self.label.setText("5.5.3 Sashimi Preview")
            self.explanation.setText(
                "This page only shows generated sashimi manifest rows and any output files already written for the current comparison. It does not auto-run sashimi."
            )
        else:
            self.label.setText("5.5.4 Failed Sashimi Jobs")
            self.explanation.setText(
                "This page only shows event-level sashimi failures for the current comparison, including the error message and output folder."
            )

    def _refresh_comparison_tabs(self) -> None:
        comparisons = self.project_service.selected_comparisons_for_display()
        current = self._current_comparison_id()
        self._comparison_ids = [item.comparison_id for item in comparisons]
        self.comparison_tabs.blockSignals(True)
        self.comparison_tabs.clear()
        for comparison in comparisons:
            self.comparison_tabs.addTab(QWidget(), comparison.display_resolved_name)
        if self._comparison_ids:
            index = self._comparison_ids.index(current) if current in self._comparison_ids else 0
            self.comparison_tabs.setCurrentIndex(index)
        self.comparison_tabs.blockSignals(False)

    def _current_comparison_id(self) -> str | None:
        index = self.comparison_tabs.currentIndex()
        if index < 0 or index >= len(self._comparison_ids):
            return None
        return self._comparison_ids[index]

    def _populate_event_type_filter(self, frame: pd.DataFrame) -> None:
        current = self.event_type_filter.currentText()
        values = ["All"]
        if "event_type" in frame.columns and not frame.empty:
            values.extend(sorted(frame["event_type"].dropna().astype(str).unique().tolist()))
        self.event_type_filter.blockSignals(True)
        self.event_type_filter.clear()
        self.event_type_filter.addItems(values)
        if current in values:
            self.event_type_filter.setCurrentText(current)
        self.event_type_filter.blockSignals(False)

    def _apply_filters(self, frame: pd.DataFrame) -> pd.DataFrame:
        if frame.empty:
            return frame
        filtered = frame.copy()
        event_type = self.event_type_filter.currentText().strip()
        if event_type and event_type != "All" and "event_type" in filtered.columns:
            filtered = filtered.loc[filtered["event_type"].astype(str) == event_type].copy()
        if "FDR" in filtered.columns:
            filtered = filtered.loc[pd.to_numeric(filtered["FDR"], errors="coerce") <= self.fdr_filter.value()].copy()
        if "dPSI" in filtered.columns:
            filtered = filtered.loc[pd.to_numeric(filtered["dPSI"], errors="coerce").abs() >= self.dpsi_filter.value()].copy()
        query = self.search.text().strip().lower()
        if query:
            mask = pd.Series(False, index=filtered.index)
            for column in ("gene_symbol", "event_type", "event_id", "comparison_display_name", "coordinates"):
                if column in filtered.columns:
                    mask |= filtered[column].fillna("").astype(str).str.lower().str.contains(query, na=False)
            filtered = filtered.loc[mask].copy()
        return filtered.reset_index(drop=True)

    def _selected_event_ids(self) -> list[str]:
        rows = sorted({index.row() for index in self.event_table.selectionModel().selectedRows()})
        if not rows:
            row = self.event_table.currentRow()
            if row >= 0:
                rows = [row]
        if not rows:
            return []
        event_col = -1
        for idx in range(self.event_table.columnCount()):
            header = self.event_table.horizontalHeaderItem(idx)
            if header is not None and header.text() == "event_id":
                event_col = idx
                break
        if event_col < 0:
            return []
        ids: list[str] = []
        for row in rows:
            item = self.event_table.item(row, event_col)
            if item is None:
                continue
            event_id = item.text().strip()
            if event_id and event_id not in ids:
                ids.append(event_id)
        return ids

    def _generate_manifest(self) -> None:
        comparison_id = self._current_comparison_id()
        event_ids = self._selected_event_ids()
        if not event_ids:
            self.status_label.setText("No events selected. Please select one or more splicing events first.")
            return
        try:
            manifest = self.project_service.build_sashimi_manifest_for_events(
                comparison_id,
                event_ids,
                allow_generate=False,
            )
        except Exception as exc:
            self.status_label.setText(f"Manifest generation failed: {exc}")
            self.details.setPlainText(f"Manifest generation failed:\n{exc}")
            return
        self.status_label.setText(
            f"Generated sashimi manifest rows: {len(manifest)}. Open 5.5.3 to inspect the manifest and output files."
        )
        self.refresh()

    def _run_sashimi(self) -> None:
        comparison_id = self._current_comparison_id()
        event_ids = self._selected_event_ids()
        if not event_ids:
            self.status_label.setText("No events selected. Please select one or more splicing events first.")
            return
        try:
            manifest = self.project_service.build_sashimi_manifest_for_events(
                comparison_id,
                event_ids,
                allow_generate=False,
            )
            output = self.project_service.run_sashimi_pipeline(manifest)
        except Exception as exc:
            self.status_label.setText(f"rmats2sashimi failed: {exc}")
            self.details.setPlainText(f"rmats2sashimi failed:\n{exc}")
            self.refresh()
            return
        failures = self.project_service.last_sashimi_failures()
        if failures.empty:
            self.status_label.setText("rmats2sashimi finished for the selected events.")
            self.details.setPlainText(output or "rmats2sashimi finished with no stdout output.")
        else:
            first = failures.iloc[0]
            self.status_label.setText(
                f"rmats2sashimi finished with {len(failures)} failed event(s). Open 5.5.4 for the failed-job table."
            )
            self.details.setPlainText(
                "\n".join(
                    [
                        output or "rmats2sashimi produced no successful stdout output.",
                        "",
                        "First failed event:",
                        f"Gene: {first.get('gene', '')}",
                        f"Event type: {first.get('event_type', '')}",
                        f"Event ID: {first.get('event_id', '')}",
                        f"Comparison: {first.get('comparison_id', '')}",
                        f"Error: {first.get('error_message', '')}",
                        f"Output folder: {first.get('output_folder', '')}",
                    ]
                )
            )
        self.refresh()

    def _show_details(self) -> None:
        if self.view_mode == "failed":
            if self._failures.empty or self.failed_table.currentRow() < 0:
                self.details.setPlainText("Select a failed sashimi row to inspect its event-level error details.")
                return
            row = min(self.failed_table.currentRow(), len(self._failures) - 1)
            record = self._failures.iloc[row]
            lines = [
                f"Gene: {record.get('gene', '')}",
                f"Event type: {record.get('event_type', '')}",
                f"Event ID: {record.get('event_id', '')}",
                f"Comparison: {record.get('comparison_id', '')}",
                f"Error: {record.get('error_message', '')}",
                f"Output folder: {record.get('output_folder', '')}",
                f"Command: {record.get('command', '')}",
                f"Script path: {record.get('script_path', '')}",
                f"stderr: {record.get('stderr', '')}",
            ]
            self.details.setPlainText("\n".join(lines))
            return

        if self.view_mode == "preview":
            if not self._outputs.empty and 0 <= self.output_table.currentRow() < len(self._outputs):
                record = self._outputs.iloc[self.output_table.currentRow()]
                self.details.setPlainText(
                    "\n".join(
                        [
                            f"Comparison: {record.get('comparison_id', '')}",
                            f"Relative path: {record.get('relative_path', '')}",
                            f"Absolute path: {record.get('absolute_path', '')}",
                            f"File type: {record.get('file_type', '')}",
                            f"Preview kind: {record.get('preview_kind', '')}",
                            f"File size: {record.get('file_size', '')}",
                            f"Modified: {record.get('last_modified', '')}",
                            "",
                            "Use Open output folder to inspect or open/download the generated sashimi file.",
                        ]
                    )
                )
                return
            if not self._manifest.empty and 0 <= self.manifest_table.currentRow() < len(self._manifest):
                record = self._manifest.iloc[self.manifest_table.currentRow()]
                self.details.setPlainText(
                    "\n".join(
                        [
                            f"Gene: {record.get('geneSymbol', '')}",
                            f"Event type: {record.get('event_type', '')}",
                            f"Event ID: {record.get('event_uid', '')}",
                            f"Comparison: {record.get('comparison_name', '')}",
                            f"Label 1: {record.get('label1', '')}",
                            f"Label 2: {record.get('label2', '')}",
                            f"Event file: {record.get('event_file', '')}",
                            f"Output dir: {record.get('outdir', '')}",
                        ]
                    )
                )
                return
            self.details.setPlainText(
                "No generated sashimi plot files were found for this comparison yet. Generate and run selected events from 5.5.2 first, then return here."
            )
            return

        if self._events.empty:
            self.details.setPlainText(
                "Open this page to inspect available candidate-event follow-up rows. "
                "This page does not auto-generate sashimi plots; select events first, then click Generate or Run."
            )
            return
        row = self.event_table.currentRow()
        if row < 0 or row >= len(self._events):
            self.details.setPlainText("Select one or more events to inspect their details before generating sashimi.")
            return
        record = self._events.iloc[row]
        lines = [
            f"Comparison: {record.get('comparison_display_name', '')}",
            f"Gene: {record.get('gene_symbol', '')}",
            f"Gene ID: {record.get('gene_id', '')}",
            f"Event type: {record.get('event_type', '')}",
            f"Event ID: {record.get('event_id', '')}",
            f"dPSI: {record.get('dPSI', '')}",
            f"FDR: {record.get('FDR', '')}",
            f"Direction: {record.get('direction', '')}",
            f"Coordinates / region: {record.get('coordinates', '')}",
        ]
        self.details.setPlainText("\n".join(lines))

    def _open_output_folder(self) -> None:
        project = self.project_service.current_project
        if project is None or project.output_root is None:
            self.status_label.setText("Output directory is not configured.")
            return
        output_dir = project.output_root / "06_sashimi"
        output_dir.mkdir(parents=True, exist_ok=True)
        os.startfile(output_dir)

    @staticmethod
    def _populate_table(table: QTableWidget, frame: pd.DataFrame) -> None:
        table.setSortingEnabled(False)
        table.clear()
        if frame is None or frame.empty:
            table.setRowCount(0)
            table.setColumnCount(0)
            table.setSortingEnabled(True)
            return
        display = frame.copy()
        table.setColumnCount(len(display.columns))
        table.setHorizontalHeaderLabels(display.columns.tolist())
        table.setRowCount(len(display))
        for row_idx, (_, row) in enumerate(display.iterrows()):
            for col_idx, column in enumerate(display.columns):
                value = row[column]
                item = QTableWidgetItem("" if pd.isna(value) else str(value))
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                table.setItem(row_idx, col_idx, item)
        table.setSortingEnabled(True)
        if table.rowCount() > 0:
            table.selectRow(0)
