from __future__ import annotations

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
    def __init__(self, project_service: ProjectService) -> None:
        super().__init__()
        self.project_service = project_service
        self.label = QLabel("Sashimi event follow-up")
        self.explanation = QLabel(
            "Sashimi is event-driven. Open this page to inspect available splicing events from the current candidate/event follow-up source, "
            "filter them, select one or more events, then generate/run rmats2sashimi only for those selected events."
        )
        self.explanation.setWordWrap(True)

        self.comparison_tabs = QTabWidget()
        self.comparison_tabs.setDocumentMode(True)
        self.comparison_tabs.setUsesScrollButtons(True)
        self.comparison_tabs.setMaximumHeight(36)
        self.comparison_tabs.currentChanged.connect(self.refresh)

        self.search = QLineEdit()
        self.search.setPlaceholderText("Search gene / event_id / event_type")
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

        self.generate_manifest_button = QPushButton("Generate selected manifest")
        self.generate_manifest_button.clicked.connect(self._generate_manifest)
        self.run_button = QPushButton("Run selected sashimi")
        self.run_button.clicked.connect(self._run_sashimi)
        self.open_output_button = QPushButton("Open output folder")
        self.open_output_button.clicked.connect(self._open_output_folder)

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

        self.failed_table = QTableWidget(0, 0)
        self.failed_table.setSortingEnabled(True)
        self.failed_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.failed_table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.failed_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)

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
        layout.addWidget(QLabel("Available splicing events"))
        layout.addWidget(self.event_table, 3)
        layout.addWidget(QLabel("Generated sashimi manifest preview"))
        layout.addWidget(self.manifest_table, 2)
        layout.addWidget(QLabel("Failed sashimi jobs"))
        layout.addWidget(self.failed_table, 1)
        layout.addWidget(self.status_label)
        layout.addWidget(self.details, 1)

        self._events = pd.DataFrame()
        self._manifest = pd.DataFrame()
        self._comparison_ids: list[str] = []

    def refresh(self) -> None:
        project = self.project_service.current_project
        if project is None:
            self._comparison_ids = []
            self.comparison_tabs.clear()
            self._populate_table(self.event_table, pd.DataFrame())
            self._populate_table(self.manifest_table, pd.DataFrame())
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
                        "comparison_id",
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

        failures = self.project_service.last_sashimi_failures()
        if not failures.empty and comparison_id and "comparison_id" in failures.columns:
            failures = failures.loc[failures["comparison_id"].astype(str) == comparison_id].copy()
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

        if events.empty:
            self.status_label.setText(
                "No available splicing events found for this comparison. "
                "Missing input: candidate event follow-up source not available. Please run 5.3.2 / 5.3.4 first."
            )
        else:
            self.status_label.setText(
                f"Available events: {len(events)} | Generated manifest rows: {len(manifest)} | Failed jobs: {len(failures)}. "
                "Select one or more events, then click 'Generate selected manifest' or 'Run selected sashimi'."
            )
        self._show_details()

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
        if not rows or self._events.empty or "event_id" not in self._events.columns:
            return []
        ids: list[str] = []
        for row in rows:
            if row < len(self._events):
                event_id = str(self._events.iloc[row]["event_id"]).strip()
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
            self.project_service.build_sashimi_manifest_for_events(
                comparison_id,
                event_ids,
                allow_generate=False,
            )
        except Exception as exc:
            self.status_label.setText(f"Manifest generation failed: {exc}")
            self.details.setPlainText(f"Manifest generation failed:\n{exc}")
            return
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
        self.status_label.setText("rmats2sashimi finished for the selected events.")
        self.details.setPlainText(output or "rmats2sashimi finished with no stdout output.")
        self.refresh()

    def _show_details(self) -> None:
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
        import os

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
