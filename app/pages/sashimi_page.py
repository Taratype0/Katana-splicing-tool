from __future__ import annotations

import os
import re
from pathlib import Path

import pandas as pd
from PySide6.QtCore import QTimer, Qt
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QDoubleSpinBox,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QScrollArea,
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

        self.runner_mode = QComboBox()
        self.runner_mode.addItems(["Katana native", "rmats2sashimiplot legacy"])
        self.runner_mode.setMinimumWidth(170)

        self.generate_manifest_button = QPushButton("Generate selected\nmanifest")
        self.generate_manifest_button.clicked.connect(self._generate_manifest)
        self.run_button = QPushButton("Run selected\nsashimi")
        self.run_button.clicked.connect(self._run_sashimi)
        self.open_output_button = QPushButton("Open output\nfolder")
        self.open_output_button.clicked.connect(self._open_output_folder)
        for button in (self.generate_manifest_button, self.run_button, self.open_output_button):
            button.setMinimumHeight(40)
            button.setMinimumWidth(130)

        self.custom_label = QLabel("Custom native event (optional, generic SE-like skipped block)")
        self.custom_label.setStyleSheet("font-weight: 600;")
        self.custom_gene = QLineEdit()
        self.custom_gene.setPlaceholderText("gene symbol")
        self.custom_gene_id = QLineEdit()
        self.custom_gene_id.setPlaceholderText("gene id optional")
        self.custom_event_id = QLineEdit()
        self.custom_event_id.setPlaceholderText("custom event id optional")
        self.custom_event_type = QComboBox()
        self.custom_event_type.addItems(["SE"])
        self.custom_chrom = QLineEdit()
        self.custom_chrom.setPlaceholderText("chr")
        self.custom_strand = QComboBox()
        self.custom_strand.addItems(["+", "-"])
        self.custom_upstream = QLineEdit()
        self.custom_upstream.setPlaceholderText("upstream start-end, 1-based")
        self.custom_event = QLineEdit()
        self.custom_event.setPlaceholderText("event/skipped block start-end, 1-based")
        self.custom_downstream = QLineEdit()
        self.custom_downstream.setPlaceholderText("downstream start-end, 1-based")
        self.custom_dpsi = QLineEdit()
        self.custom_dpsi.setPlaceholderText("dPSI optional")
        self.custom_fdr = QLineEdit()
        self.custom_fdr.setPlaceholderText("FDR optional")
        self.generate_custom_manifest_button = QPushButton("Generate custom\nmanifest")
        self.generate_custom_manifest_button.clicked.connect(self._generate_custom_manifest)
        self.run_custom_button = QPushButton("Run custom\nnative sashimi")
        self.run_custom_button.clicked.connect(self._run_custom_native_sashimi)
        for button in (self.generate_custom_manifest_button, self.run_custom_button):
            button.setMinimumHeight(40)
            button.setMinimumWidth(150)

        self.detect_label = QLabel("Detect significant continuous skipped-exon junctions")
        self.detect_label.setStyleSheet("font-weight: 600;")
        self.detect_gene = QLineEdit()
        self.detect_gene.setPlaceholderText("gene symbol or gene id required")
        self.detect_fdr = QDoubleSpinBox()
        self.detect_fdr.setRange(0.0, 1.0)
        self.detect_fdr.setSingleStep(0.01)
        self.detect_fdr.setValue(0.05)
        self.detect_fdr.setPrefix("q <= ")
        self.detect_dpsi = QDoubleSpinBox()
        self.detect_dpsi.setRange(0.0, 5.0)
        self.detect_dpsi.setSingleStep(0.05)
        self.detect_dpsi.setValue(0.10)
        self.detect_dpsi.setPrefix("|dPSI| >= ")
        self.detect_button = QPushButton("Detect significant rMATS\ncontinuous skips")
        self.detect_button.clicked.connect(self._detect_continuous_skips)
        self.detect_button.setMinimumHeight(40)
        self.detect_button.setMinimumWidth(180)

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
        self.output_label = QLabel("Generated sashimi PNG previews")
        self.failed_label = QLabel("Failed sashimi jobs")
        self.preview_label = QLabel("Sashimi image preview")
        self.preview_image = QLabel("Select a generated rMATS event PNG row to preview it here.")
        self.preview_image.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.preview_image.setMinimumHeight(320)
        self.preview_image.setWordWrap(True)
        self.preview_scroll = QScrollArea()
        self.preview_scroll.setWidgetResizable(True)
        self.preview_scroll.setWidget(self.preview_image)
        self.preview_scroll.setMinimumHeight(420)

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
        filters.addWidget(self.runner_mode)
        filters.addWidget(self.generate_manifest_button)
        filters.addWidget(self.run_button)
        filters.addWidget(self.open_output_button)

        custom_grid = QGridLayout()
        custom_grid.addWidget(self.custom_label, 0, 0, 1, 4)
        custom_grid.addWidget(self.custom_gene, 1, 0)
        custom_grid.addWidget(self.custom_gene_id, 1, 1)
        custom_grid.addWidget(self.custom_event_id, 1, 2)
        custom_grid.addWidget(self.custom_event_type, 1, 3)
        custom_grid.addWidget(self.custom_chrom, 2, 0)
        custom_grid.addWidget(self.custom_strand, 2, 1)
        custom_grid.addWidget(self.custom_dpsi, 2, 2)
        custom_grid.addWidget(self.custom_fdr, 2, 3)
        custom_grid.addWidget(self.custom_upstream, 3, 0)
        custom_grid.addWidget(self.custom_event, 3, 1)
        custom_grid.addWidget(self.custom_downstream, 3, 2)
        custom_grid.addWidget(self.generate_custom_manifest_button, 3, 3)
        custom_grid.addWidget(self.run_custom_button, 3, 4)
        custom_grid.addWidget(self.detect_label, 4, 0, 1, 4)
        custom_grid.addWidget(self.detect_gene, 5, 0)
        custom_grid.addWidget(self.detect_fdr, 5, 1)
        custom_grid.addWidget(self.detect_dpsi, 5, 2)
        custom_grid.addWidget(self.detect_button, 5, 3)

        layout = QVBoxLayout(self)
        layout.addWidget(self.label)
        layout.addWidget(self.explanation)
        layout.addWidget(self.comparison_tabs)
        layout.addLayout(filters)
        layout.addLayout(custom_grid)
        layout.addWidget(self.events_label)
        layout.addWidget(self.event_table, 3)
        layout.addWidget(self.manifest_label)
        layout.addWidget(self.manifest_table, 2)
        layout.addWidget(self.output_label)
        layout.addWidget(self.output_table, 2)
        layout.addWidget(self.preview_label)
        layout.addWidget(self.preview_scroll, 3)
        layout.addWidget(self.failed_label)
        layout.addWidget(self.failed_table, 1)
        layout.addWidget(self.status_label)
        layout.addWidget(self.details, 1)

        self._comparison_ids: list[str] = []
        self._events = pd.DataFrame()
        self._manifest = pd.DataFrame()
        self._outputs = pd.DataFrame()
        self._failures = pd.DataFrame()
        self._preview_pixmap = QPixmap()
        self._detected_events_by_comparison: dict[str, pd.DataFrame] = {}
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
            self._clear_preview("Load a project and run selected sashimi events to preview outputs.")
            self.status_label.setText("Load a project first.")
            self.details.setPlainText("")
            return

        self._refresh_comparison_tabs()
        comparison_id = self._current_comparison_id()

        events = self.project_service.available_sashimi_events(comparison_id, allow_generate=False)
        detected = self._detected_events_by_comparison.get(str(comparison_id or ""), pd.DataFrame())
        if not detected.empty:
            events = pd.concat([events, detected], ignore_index=True, sort=False) if not events.empty else detected.copy()
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
                        "source",
                        "event_type",
                        "event_id",
                        "comparison_display_name",
                        "dPSI",
                        "FDR",
                        "skipped_exons",
                        "transcript_name",
                        "direction",
                        "coordinates",
                        "splicing_direction_flipped",
                    ]
                    if column in events.columns
                ]
            ].copy()
            if not events.empty
            else pd.DataFrame(),
            select_first=False,
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
                        "dPSI",
                        "FDR",
                    ]
                    if column in manifest.columns
                ]
            ].copy()
            if not manifest.empty
            else pd.DataFrame(),
        )

        outputs = self.project_service.build_sashimi_output_browser(comparison_id)
        if self.view_mode == "preview":
            outputs = self._preview_png_output_rows(outputs)
        self._outputs = outputs.copy()
        self._populate_table(
            self.output_table,
            outputs[
                [
                    column
                    for column in [
                        "gene",
                        "event_type",
                        "event_id",
                        "comparison_display_name",
                        "dPSI",
                        "FDR",
                        "comparison_id",
                        "relative_path",
                        "last_modified",
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
                    "Generated sashimi PNG previews: 0. "
                    "No generated PNG preview was found yet. Run selected rMATS events from 5.5.2 first; "
                    "use Open output folder for PDFs/SVGs/tables."
                )
            else:
                self.status_label.setText(
                    f"Generated sashimi PNG preview events: {len(outputs)}. Select a row to preview the PNG."
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
        self.runner_mode.setVisible(run_mode)
        self.generate_manifest_button.setVisible(run_mode)
        self.run_button.setVisible(run_mode)
        for widget in (
            self.custom_label,
            self.custom_gene,
            self.custom_gene_id,
            self.custom_event_id,
            self.custom_event_type,
            self.custom_chrom,
            self.custom_strand,
            self.custom_upstream,
            self.custom_event,
            self.custom_downstream,
            self.custom_dpsi,
            self.custom_fdr,
            self.generate_custom_manifest_button,
            self.run_custom_button,
            self.detect_label,
            self.detect_gene,
            self.detect_fdr,
            self.detect_dpsi,
            self.detect_button,
        ):
            widget.setVisible(run_mode)

        self.events_label.setVisible(run_mode)
        self.event_table.setVisible(run_mode)
        self.manifest_label.setVisible(run_mode)
        self.manifest_table.setVisible(run_mode)
        self.output_label.setVisible(preview_mode)
        self.output_table.setVisible(preview_mode)
        self.preview_label.setVisible(preview_mode)
        self.preview_scroll.setVisible(preview_mode)
        self.failed_label.setVisible(failed_mode)
        self.failed_table.setVisible(failed_mode)

        if run_mode:
            self.label.setText("5.5.2 Run Sashimi")
            self.explanation.setText(
                "This page lists available rMATS splicing events for the current comparison. Select events here, or define a custom SE-like skipped block by generic coordinates, then run Katana native sashimi only for the selected/custom events."
            )
        elif preview_mode:
            self.label.setText("5.5.3 Sashimi Preview")
            self.explanation.setText(
                "This page only shows generated rMATS event PNG previews for the current comparison. It does not auto-run sashimi; use Open output folder for non-preview files."
            )
        else:
            self.label.setText("Sashimi Diagnostics")
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

    def _selected_event_rows(self) -> pd.DataFrame:
        rows = sorted({index.row() for index in self.event_table.selectionModel().selectedRows()})
        source_rows = [self._source_row_for_table(self.event_table, row) for row in rows]
        valid_rows = [row for row in source_rows if row is not None and 0 <= row < len(self._events)]
        return self._events.iloc[valid_rows].copy() if valid_rows else pd.DataFrame()

    def _build_manifest_for_selected_rows(self, comparison_id: str | None, selected_rows: pd.DataFrame) -> pd.DataFrame:
        if selected_rows.empty:
            raise RuntimeError("No events selected. Please select one or more splicing events first.")
        detected_mask = selected_rows.get("source", pd.Series("", index=selected_rows.index)).fillna("").astype(str) == "continuous_skip_detector"
        detected_rows = selected_rows.loc[detected_mask].copy()
        rmats_rows = selected_rows.loc[~detected_mask].copy()
        frames: list[pd.DataFrame] = []
        if not rmats_rows.empty:
            event_ids = rmats_rows["event_id"].dropna().astype(str).drop_duplicates().tolist()
            if event_ids:
                frames.append(
                    self.project_service.build_sashimi_manifest_for_events(
                        comparison_id,
                        event_ids,
                        allow_generate=False,
                    )
                )
        if not detected_rows.empty:
            frames.append(
                self.project_service.build_custom_native_sashimi_manifest_for_rows(
                    comparison_id,
                    detected_rows,
                )
            )
        manifest = pd.concat([frame for frame in frames if frame is not None and not frame.empty], ignore_index=True) if frames else pd.DataFrame()
        if manifest.empty:
            raise RuntimeError("No sashimi manifest rows could be generated for the selected events.")
        self.project_service.run_state.results.sashimi_manifest = manifest.copy()
        return manifest

    def _generate_manifest(self) -> None:
        comparison_id = self._current_comparison_id()
        selected_rows = self._selected_event_rows()
        if selected_rows.empty:
            self.status_label.setText("No events selected. Please select one or more splicing events first.")
            return
        try:
            manifest = self._build_manifest_for_selected_rows(comparison_id, selected_rows)
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
        selected_rows = self._selected_event_rows()
        if selected_rows.empty:
            self.status_label.setText("No events selected. Please select one or more splicing events first.")
            return
        try:
            manifest = self._build_manifest_for_selected_rows(comparison_id, selected_rows)
            has_detected = (
                selected_rows.get("source", pd.Series("", index=selected_rows.index)).fillna("").astype(str)
                == "continuous_skip_detector"
            ).any()
            use_native = has_detected or self.runner_mode.currentText().startswith("Katana native")
            output = (
                self.project_service.run_native_sashimi_pipeline(manifest)
                if use_native
                else self.project_service.run_sashimi_pipeline(manifest)
            )
        except Exception as exc:
            self.status_label.setText(f"Sashimi failed: {exc}")
            self.details.setPlainText(f"Sashimi failed:\n{exc}")
            self.refresh()
            return
        failures = self.project_service.last_sashimi_failures()
        runner_name = "Katana native" if use_native else self.runner_mode.currentText()
        if failures.empty:
            self.status_label.setText(f"{runner_name} finished for the selected events.")
            self.details.setPlainText(output or f"{runner_name} finished with no stdout output.")
        else:
            first = failures.iloc[0]
            self.status_label.setText(
                f"{runner_name} finished with {len(failures)} failed event(s). Failed rows were written to the sashimi output folder."
            )
            self.details.setPlainText(
                "\n".join(
                    [
                        output or f"{runner_name} produced no successful stdout output.",
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

    def _detect_continuous_skips(self) -> None:
        comparison_id = self._current_comparison_id()
        try:
            detected = self.project_service.detect_significant_continuous_skipped_exon_junctions(
                comparison_id,
                gene_query=self.detect_gene.text().strip(),
                fdr_cutoff=self.detect_fdr.value(),
                abs_dpsi_cutoff=self.detect_dpsi.value(),
            )
        except Exception as exc:
            self.status_label.setText(f"Continuous skip detection failed: {exc}")
            self.details.setPlainText(f"Continuous skip detection failed:\n{exc}")
            return
        self._detected_events_by_comparison[str(comparison_id or "")] = detected.copy()
        if detected.empty:
            self.status_label.setText(
                "No significant continuous skipped-exon junctions found for this gene under the active rMATS FDR/dPSI thresholds."
            )
            self.details.setPlainText(
                "No events found. Detection only reports rMATS SE rows passing FDR/q cutoff and abs(dPSI) cutoff, then labels continuous skipped exons from the GTF."
            )
        else:
            self.status_label.setText(
                f"Detected {len(detected)} significant continuous skipped-exon junction(s). Select rows and run Katana native sashimi."
            )
            self.details.setPlainText(
                "Detected continuous skipped-exon junctions are rMATS-significant SE events filtered by FDR/q and abs(dPSI), then annotated with GTF exon labels. They are not auto-run."
            )
        self.refresh()

    @staticmethod
    def _parse_range(text: str, label: str) -> tuple[int, int]:
        values = re.findall(r"\d[\d,]*", text or "")
        if len(values) < 2:
            raise RuntimeError(f"{label} must contain start and end coordinates, e.g. 12345-12500.")
        start = int(values[0].replace(",", ""))
        end = int(values[1].replace(",", ""))
        if start <= 0 or end <= 0:
            raise RuntimeError(f"{label} coordinates must be positive 1-based positions.")
        return min(start, end), max(start, end)

    def _custom_event_payload(self) -> dict[str, object]:
        upstream_start, upstream_end = self._parse_range(self.custom_upstream.text(), "upstream")
        event_start, event_end = self._parse_range(self.custom_event.text(), "event/skipped block")
        downstream_start, downstream_end = self._parse_range(self.custom_downstream.text(), "downstream")
        return {
            "gene_symbol": self.custom_gene.text().strip(),
            "gene_id": self.custom_gene_id.text().strip(),
            "event_id": self.custom_event_id.text().strip(),
            "event_type": self.custom_event_type.currentText().strip(),
            "chrom": self.custom_chrom.text().strip(),
            "strand": self.custom_strand.currentText().strip(),
            "upstream_start": upstream_start,
            "upstream_end": upstream_end,
            "event_start": event_start,
            "event_end": event_end,
            "downstream_start": downstream_start,
            "downstream_end": downstream_end,
            "dpsi": self.custom_dpsi.text().strip(),
            "fdr": self.custom_fdr.text().strip(),
        }

    def _generate_custom_manifest(self) -> None:
        comparison_id = self._current_comparison_id()
        try:
            manifest = self.project_service.build_custom_native_sashimi_manifest(
                comparison_id,
                self._custom_event_payload(),
            )
        except Exception as exc:
            self.status_label.setText(f"Custom manifest generation failed: {exc}")
            self.details.setPlainText(f"Custom manifest generation failed:\n{exc}")
            return
        self.status_label.setText(
            f"Generated custom native sashimi manifest rows: {len(manifest)}. Run custom native sashimi when ready."
        )
        self.refresh()

    def _run_custom_native_sashimi(self) -> None:
        comparison_id = self._current_comparison_id()
        try:
            manifest = self.project_service.build_custom_native_sashimi_manifest(
                comparison_id,
                self._custom_event_payload(),
            )
            output = self.project_service.run_native_sashimi_pipeline(manifest)
        except Exception as exc:
            self.status_label.setText(f"Custom native sashimi failed: {exc}")
            self.details.setPlainText(f"Custom native sashimi failed:\n{exc}")
            self.refresh()
            return
        failures = self.project_service.last_sashimi_failures()
        if failures.empty:
            self.status_label.setText("Katana native sashimi finished for the custom event.")
            self.details.setPlainText(output or "Katana native sashimi finished with no stdout output.")
        else:
            first = failures.iloc[0]
            self.status_label.setText(
                f"Katana native sashimi finished with {len(failures)} failed custom event(s). Failed rows were written to the sashimi output folder."
            )
            self.details.setPlainText(
                "\n".join(
                    [
                        output or "Katana native sashimi produced no successful stdout output.",
                        "",
                        "First failed custom event:",
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
            output_source_row = self._source_row_for_table(self.output_table, self.output_table.currentRow())
            if output_source_row is not None and not self._outputs.empty and 0 <= output_source_row < len(self._outputs):
                record = self._outputs.iloc[output_source_row]
                self._preview_output_record(record)
                self.details.setPlainText(
                    "\n".join(
                        [
                            f"Comparison: {record.get('comparison_id', '')}",
                            f"Relative path: {record.get('relative_path', '')}",
                            f"Absolute path: {record.get('absolute_path', '')}",
                            f"dPSI: {record.get('dPSI', '')}",
                            f"FDR: {record.get('FDR', '')}",
                            f"Modified: {record.get('last_modified', '')}",
                            "",
                            "Only PNG previews are listed here. Use Open output folder to inspect PDF/SVG/table outputs.",
                        ]
                    )
                )
                return
            self._clear_preview("No generated PNG preview selected.")
            self.details.setPlainText(
                "No generated PNG previews were found for this comparison yet. Generate and run selected events from 5.5.2 first, then return here."
            )
            return

        if self._events.empty:
            self.details.setPlainText(
                "Open this page to inspect available candidate-event follow-up rows. "
                "This page does not auto-generate sashimi plots; select events first, then click Generate or Run."
            )
            return
        row = self.event_table.currentRow()
        source_row = self._source_row_for_table(self.event_table, row)
        if source_row is None or source_row < 0 or source_row >= len(self._events):
            self.details.setPlainText("Select one or more events to inspect their details before generating sashimi.")
            return
        record = self._events.iloc[source_row]
        lines = [
            f"Comparison: {record.get('comparison_display_name', '')}",
            f"Source: {record.get('source', 'rMATS')}",
            f"Gene: {record.get('gene_symbol', '')}",
            f"Gene ID: {record.get('gene_id', '')}",
            f"Transcript: {record.get('transcript_name', '')}",
            f"Event type: {record.get('event_type', '')}",
            f"Event ID: {record.get('event_id', '')}",
            f"Skipped exons: {record.get('skipped_exons', '')}",
            f"PSI group1 / group2: {record.get('PSI_group1', '')} / {record.get('PSI_group2', '')}",
            f"Skip reads group1 / group2: {record.get('skip_reads_group1', '')} / {record.get('skip_reads_group2', '')}",
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
    def _populate_table(table: QTableWidget, frame: pd.DataFrame, *, select_first: bool = True) -> None:
        table.blockSignals(True)
        try:
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
                    item.setData(Qt.ItemDataRole.UserRole, row_idx)
                    table.setItem(row_idx, col_idx, item)
            table.setSortingEnabled(True)
            if select_first and table.rowCount() > 0:
                table.selectRow(0)
            elif not select_first:
                table.clearSelection()
        finally:
            table.blockSignals(False)

    @staticmethod
    def _source_row_for_table(table: QTableWidget, visual_row: int) -> int | None:
        if visual_row < 0 or visual_row >= table.rowCount():
            return None
        item = table.item(visual_row, 0)
        if item is None:
            return visual_row
        value = item.data(Qt.ItemDataRole.UserRole)
        try:
            return int(value)
        except (TypeError, ValueError):
            return visual_row

    def _clear_preview(self, message: str) -> None:
        self._preview_pixmap = QPixmap()
        self.preview_image.setPixmap(QPixmap())
        self.preview_image.setText(message)

    def _preview_output_record(self, record: pd.Series) -> None:
        path = Path(str(record.get("absolute_path", "")))
        preview_kind = str(record.get("preview_kind", "")).lower()
        suffix = path.suffix.lower()
        if not path.exists():
            self._clear_preview(f"Missing output file:\n{path}")
            return
        if preview_kind != "image" and suffix not in {".png", ".jpg", ".jpeg", ".svg"}:
            self._clear_preview(
                "\n".join(
                    [
                        "No embedded image preview for this file.",
                        f"File: {path}",
                        "Use Open output folder to inspect it directly.",
                    ]
                )
            )
            return
        pixmap = QPixmap(str(path))
        if pixmap.isNull():
            self._clear_preview(f"Could not load image preview:\n{path}")
            return
        self._preview_pixmap = pixmap
        self.preview_image.setText("")
        self._fit_preview_pixmap()
        QTimer.singleShot(0, self._fit_preview_pixmap)
        QTimer.singleShot(120, self._fit_preview_pixmap)

    def _fit_preview_pixmap(self) -> None:
        if self._preview_pixmap.isNull() or not self.preview_scroll.isVisible():
            return
        viewport = self.preview_scroll.viewport()
        max_width = max(viewport.width() - 28, 360)
        max_height = max(viewport.height() - 28, 260)
        pixmap = self._preview_pixmap.scaled(
            max_width,
            max_height,
            Qt.AspectRatioMode.KeepAspectRatio,
            Qt.TransformationMode.SmoothTransformation,
        )
        self.preview_image.setPixmap(pixmap)
        self.preview_image.setMinimumSize(0, 0)
        self.preview_image.adjustSize()

    def resizeEvent(self, event) -> None:
        super().resizeEvent(event)
        self._fit_preview_pixmap()

    @staticmethod
    def _preview_png_output_rows(frame: pd.DataFrame) -> pd.DataFrame:
        if frame is None or frame.empty:
            return pd.DataFrame(columns=frame.columns if frame is not None else [])
        working = frame.copy()
        absolute = working.get("absolute_path", pd.Series("", index=working.index)).fillna("").astype(str)
        suffix = absolute.map(lambda value: Path(value).suffix.lower())
        file_type = working.get("file_type", pd.Series("", index=working.index)).fillna("").astype(str).str.lower()
        png = working.loc[(file_type == "png") | (suffix == ".png")].copy()
        if png.empty:
            return png
        relative = png.get("relative_path", pd.Series("", index=png.index)).fillna("").astype(str).str.replace("\\", "/", regex=False).str.lower()
        summary_mask = relative.str.endswith("native_summary.png") | relative.str.endswith("summary.png")
        if summary_mask.any():
            png = png.loc[summary_mask].copy()
        key_columns = [column for column in ["comparison_id", "gene", "event_type", "event_id"] if column in png.columns]
        if key_columns:
            png = png.drop_duplicates(subset=key_columns, keep="first")
        return png.reset_index(drop=True)
