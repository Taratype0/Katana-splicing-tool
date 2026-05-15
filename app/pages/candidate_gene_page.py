from __future__ import annotations

import os
import shutil

import pandas as pd
from pandas.api.types import is_bool_dtype, is_numeric_dtype
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QComboBox,
    QDoubleSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from src.services.project_service import ProjectService


class CandidateGenePage(QWidget):
    def __init__(self, project_service: ProjectService, mode: str = "ranking") -> None:
        super().__init__()
        self.project_service = project_service
        self.mode = mode
        self.label = QLabel("Candidate gene table")
        self.context_label = QLabel("")
        self.context_label.setWordWrap(True)
        self.search = QLineEdit()
        self.search.setPlaceholderText("Search gene, comparison, class...")
        self.search.textChanged.connect(self.refresh)
        self.comparison_tabs = QTabWidget()
        self.comparison_tabs.currentChanged.connect(self.refresh)
        self.comparison_tabs.setDocumentMode(True)
        self.comparison_tabs.setStyleSheet("QTabWidget::pane { border: 0; margin: 0; padding: 0; }")
        self.comparison_tabs.setUsesScrollButtons(True)
        self.comparison_tabs.setMaximumHeight(36)
        self.comparison_filter = QComboBox()
        self.comparison_filter.currentTextChanged.connect(self.refresh)
        self.comparison_filter.setVisible(False)
        self.class_filter = QComboBox()
        self.class_filter.addItems(["All", "Tier 1", "Tier 2", "Tier 3", "Tier 4", "Tier 5", "Unclassified"])
        self.class_filter.currentTextChanged.connect(self.refresh)
        self.evidence_filter = QComboBox()
        self.evidence_filter.addItems(["All", "DEG_rMATS_DEXSeq_DTU", "DEG_rMATS_DEXSeq", "DEG_rMATS_DTU", "DEG_rMATS", "rMATS_DEXSeq", "rMATS_DTU", "rMATS", "DEG", "unclassified"])
        self.evidence_filter.currentTextChanged.connect(self.refresh)
        self.direction_filter = QComboBox()
        self.direction_filter.addItems(
            [
                "All",
                "transcript_up_splicing_up",
                "transcript_up_splicing_down",
                "transcript_down_splicing_up",
                "transcript_down_splicing_down",
                "DEG_only",
                "rMATS_only",
                "rMATS_DEXSeq_or_DTU_only",
                "unclassified",
            ]
        )
        self.direction_filter.currentTextChanged.connect(self.refresh)
        self.min_abs_log2fc = QDoubleSpinBox()
        self.min_abs_log2fc.setRange(0.0, 100.0)
        self.min_abs_log2fc.setSingleStep(0.1)
        self.min_abs_log2fc.setPrefix("|log2FC| >= ")
        self.min_abs_log2fc.valueChanged.connect(self.refresh)
        self.min_abs_dpsi = QDoubleSpinBox()
        self.min_abs_dpsi.setRange(0.0, 100.0)
        self.min_abs_dpsi.setSingleStep(0.05)
        self.min_abs_dpsi.setPrefix("|dPSI| >= ")
        self.min_abs_dpsi.valueChanged.connect(self.refresh)
        self.path_label = QLineEdit()
        self.path_label.setReadOnly(True)
        self.path_label.setPlaceholderText("Current source file path")
        self.generate_button = QPushButton("Generate / Load Candidate Results")
        self.generate_button.clicked.connect(self._generate_candidate_results)
        self.open_button = QPushButton("Open output folder")
        self.open_button.clicked.connect(self.open_output_folder)
        self.download_button = QPushButton("Download file")
        self.download_button.clicked.connect(self.download_current_file)
        self.export_button = QPushButton("Export Current View")
        self.export_button.clicked.connect(self.export_current_view)
        self.table = QTableWidget(0, 0)
        self.table.setSortingEnabled(True)
        self._current_frame = pd.DataFrame()
        self._current_source_path = None

        top = QHBoxLayout()
        top.addWidget(self.label)
        top.addWidget(self.search)
        top.addWidget(self.comparison_filter)
        top.addWidget(self.class_filter)
        top.addWidget(self.evidence_filter)
        top.addWidget(self.direction_filter)
        top.addWidget(self.min_abs_log2fc)
        top.addWidget(self.min_abs_dpsi)
        top.addWidget(self.generate_button)
        top.addWidget(self.open_button)
        top.addWidget(self.download_button)
        top.addWidget(self.export_button)

        layout = QVBoxLayout(self)
        layout.addLayout(top)
        layout.addWidget(self.context_label)
        layout.addWidget(self.comparison_tabs)
        layout.addWidget(self.path_label)
        layout.addWidget(self.table)

    def refresh(self) -> None:
        project = self.project_service.current_project
        if project is None:
            self.table.setRowCount(0)
            self.table.setColumnCount(0)
            self.context_label.clear()
            return

        cached = self.project_service.run_state.results.candidate_gene_table
        if cached is None or cached.empty:
            self.table.setRowCount(0)
            self.table.setColumnCount(0)
            self.label.setText("Candidate gene table (not run)")
            self.context_label.setText("No cached candidate ranking is loaded. Click 'Generate / Load Candidate Results' to build it for the current project.")
            return
        frame = self.project_service.preview_candidate_gene_screening(allow_generate=False)
        if frame is None or frame.empty:
            self.table.setRowCount(0)
            self.table.setColumnCount(0)
            self.label.setText("Candidate gene table (no cached result)")
            self.context_label.setText("No cached candidate ranking is loaded in memory. Click 'Generate / Load Candidate Results' to compute/load it.")
            return
        self._refresh_comparison_filter(frame)
        self._refresh_comparison_tabs(frame)

        filtered = self._filter_frame(frame, self.search.text())
        self._current_frame = filtered.copy()
        self._update_context_label()
        self._current_source_path = self._resolve_current_source_path()
        self.path_label.setText(str(self._current_source_path) if self._current_source_path else "Source file: in-memory filtered view")
        if self.mode == "tiers":
            columns = [
                "comparison_id",
                "rank",
                "gene_symbol",
                "gene_id",
                "n_rMATS_significant_events",
                "significant_as_event_types",
                "significant_as_event_ids",
                "significant_as_event_list",
                "as_event_summary",
                "candidate_tier",
                "evidence_class",
                "evidence_count",
                "direction_class",
                "DE_FDR",
                "standardized_log2FC",
                "best_rMATS_FDR",
                "dominant_rMATS_standardized_dPSI",
                "dominant_rMATS_event_type",
                "candidate_reason",
            ]
        else:
            columns = [
                "comparison_id",
                "group1",
                "group2",
                "rank",
                "gene_symbol",
                "gene_id",
                "n_rMATS_significant_events",
                "significant_as_event_types",
                "significant_as_event_ids",
                "significant_as_event_list",
                "as_event_summary",
                "candidate_tier",
                "evidence_class",
                "evidence_count",
                "DEG_significant",
                "rMATS_significant",
                "dominant_rMATS_event_type",
                "best_rMATS_event_id",
                "DEXSeq_significant",
                "DTU_significant",
                "DE_FDR",
                "standardized_log2FC",
                "best_rMATS_FDR",
                "dominant_rMATS_standardized_dPSI",
                "direction_class",
                "shortlist_gene",
                "blacklist_gene",
            ]
        available_columns = [column for column in columns if column in filtered.columns]
        self.table.setSortingEnabled(False)
        self.table.setColumnCount(len(available_columns))
        self.table.setHorizontalHeaderLabels(available_columns)
        self.table.setRowCount(len(filtered))
        for row_idx, (_, row) in enumerate(filtered[available_columns].iterrows()):
            for col_idx, column in enumerate(available_columns):
                value = row[column]
                item = QTableWidgetItem()
                if pd.isna(value):
                    item.setText("")
                elif self._is_numeric_column(filtered, column):
                    numeric_value = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
                    if pd.notna(numeric_value):
                        if float(numeric_value).is_integer():
                            item.setData(Qt.ItemDataRole.EditRole, int(numeric_value))
                        else:
                            item.setData(Qt.ItemDataRole.EditRole, float(numeric_value))
                    else:
                        item.setText(str(value))
                elif self._is_boolean_column(filtered, column):
                    item.setData(Qt.ItemDataRole.EditRole, bool(value))
                else:
                    item.setText(str(value))
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(row_idx, col_idx, item)
        self.table.setSortingEnabled(True)
        label_prefix = "Candidate auto ranking" if self.mode == "ranking" else "Tier 1-5 candidate table"
        self.label.setText(f"{label_prefix} ({len(filtered)} rows)")

    def _generate_candidate_results(self) -> None:
        frame = self.project_service.preview_candidate_gene_screening(allow_generate=True, force_rebuild=True)
        if frame is None or frame.empty:
            if self.project_service.run_state.results.tx_splicing_gene_table.empty:
                self.context_label.setText(
                    "Missing input: transcript + splicing integration not run for this project. Please run 5.2.1 first."
                )
            else:
                self.context_label.setText(
                    "Candidate ranking generation finished, but no candidate rows were produced for the current project."
                )
        self.refresh()

    def _resolve_current_source_path(self):
        comparison_id = self._selected_comparison()
        if not comparison_id:
            return None
        file_map = self.project_service.candidate_screening_file_map(comparison_id)
        if self.mode == "tiers":
            selected_tier = self.class_filter.currentText()
            tier_key = {
                "Tier 1": "tier1",
                "Tier 2": "tier2",
                "Tier 3": "tier3",
                "Tier 4": "tier4",
                "Tier 5": "tier5",
            }.get(selected_tier)
            if tier_key and tier_key in file_map:
                return file_map[tier_key]
        return file_map.get("gene_level_integrated_candidates")

    def _filter_frame(self, frame: pd.DataFrame, query: str) -> pd.DataFrame:
        filtered = frame.copy()
        if not filtered.empty:
            filtered = self.project_service._dedupe_candidate_frame(filtered)
        selected_comparison = self._selected_comparison()
        if selected_comparison and "comparison_id" in filtered.columns:
            filtered = filtered.loc[filtered["comparison_id"].astype(str) == selected_comparison].copy()

        selected_class = self.class_filter.currentText()
        if selected_class != "All" and "candidate_tier" in filtered.columns:
            filtered = filtered.loc[filtered["candidate_tier"].astype(str) == selected_class].copy()

        selected_evidence = self.evidence_filter.currentText()
        if selected_evidence != "All" and "evidence_class" in filtered.columns:
            filtered = filtered.loc[filtered["evidence_class"].astype(str) == selected_evidence].copy()

        selected_direction = self.direction_filter.currentText()
        if selected_direction != "All" and "direction_class" in filtered.columns:
            filtered = filtered.loc[filtered["direction_class"].astype(str) == selected_direction].copy()

        if "standardized_log2FC" in filtered.columns:
            log2fc = pd.to_numeric(filtered["standardized_log2FC"], errors="coerce").abs()
            filtered = filtered.loc[log2fc >= self.min_abs_log2fc.value()].copy()

        if "dominant_rMATS_standardized_dPSI" in filtered.columns:
            dpsi = pd.to_numeric(filtered["dominant_rMATS_standardized_dPSI"], errors="coerce").abs()
            filtered = filtered.loc[dpsi >= self.min_abs_dpsi.value()].copy()

        query = query.strip().lower()
        if not query:
            return filtered
        mask = pd.Series(False, index=filtered.index)
        for column in ["comparison_id", "gene_symbol", "gene_id", "candidate_tier", "evidence_class", "direction_class", "candidate_reason"]:
            if column in filtered.columns:
                mask = mask | filtered[column].astype(str).str.lower().str.contains(query, na=False)
        return filtered.loc[mask].copy()

    def _refresh_comparison_filter(self, frame: pd.DataFrame) -> None:
        if "comparison_id" not in frame.columns:
            return
        ordered_ids = [item.comparison_id for item in self.project_service.selected_comparisons_for_display()]
        available = frame["comparison_id"].dropna().astype(str).unique().tolist()
        values = [item for item in ordered_ids if item in available]
        for item in available:
            if item not in values:
                values.append(item)
        current = self.comparison_filter.currentText()
        self.comparison_filter.blockSignals(True)
        self.comparison_filter.clear()
        self.comparison_filter.addItems(values)
        if current in values:
            self.comparison_filter.setCurrentText(current)
        elif values:
            self.comparison_filter.setCurrentIndex(0)
        self.comparison_filter.blockSignals(False)

    def _refresh_comparison_tabs(self, frame: pd.DataFrame) -> None:
        if "comparison_id" not in frame.columns:
            self.comparison_tabs.clear()
            return
        ordered_pairs = [(item.comparison_id, item.display_resolved_name) for item in self.project_service.selected_comparisons_for_display()]
        available = frame["comparison_id"].dropna().astype(str).unique().tolist()
        values = [item for item in ordered_pairs if item[0] in available]
        for comparison_id in available:
            if comparison_id not in {item[0] for item in values}:
                values.append((comparison_id, comparison_id))
        current = self._selected_comparison()
        self.comparison_tabs.blockSignals(True)
        self.comparison_tabs.clear()
        ids = [item[0] for item in values]
        for index, (comparison_id, label) in enumerate(values):
            self.comparison_tabs.addTab(QWidget(), label)
            self.comparison_tabs.setTabToolTip(index, comparison_id)
        target = ids.index(current) if current in ids else 0
        self.comparison_tabs.setCurrentIndex(target)
        self.comparison_tabs.blockSignals(False)

    def _selected_comparison(self) -> str:
        if self.comparison_tabs.count() > 0:
            tooltip = self.comparison_tabs.tabToolTip(self.comparison_tabs.currentIndex())
            if tooltip:
                return tooltip
        return self.comparison_filter.currentText() or ""

    def export_current_view(self) -> None:
        if self._current_frame.empty:
            return
        output_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export Current Candidate View",
            "candidate_genes.tsv",
            "TSV Files (*.tsv)",
        )
        if not output_path:
            return
        self.project_service.export_frame(self._current_frame, output_path)

    def open_output_folder(self) -> None:
        source = self._current_source_path
        if source is None:
            return
        folder = source.parent
        if folder.exists():
            os.startfile(str(folder))

    def download_current_file(self) -> None:
        source = self._current_source_path
        output_path, _ = QFileDialog.getSaveFileName(
            self,
            "Download Current Candidate File",
            "candidate_source.tsv",
            "TSV Files (*.tsv);;All Files (*.*)",
        )
        if not output_path:
            return
        if source is not None and source.exists():
            shutil.copy2(source, output_path)
            return
        if not self._current_frame.empty:
            self.project_service.export_frame(self._current_frame, output_path)

    def _update_context_label(self) -> None:
        selected = self._selected_comparison()
        if not selected:
            self.context_label.setText("No comparison selected.")
            return
        current = self._current_frame.copy()
        unique_gene_count = 0
        if not current.empty:
            gene_keys = current.get("candidate_unique_key")
            if gene_keys is None:
                gene_keys = current.get("gene_symbol", pd.Series("", index=current.index)).fillna("").astype(str).str.strip()
                if "gene_id" in current.columns:
                    gene_keys = gene_keys.where(gene_keys.ne(""), current["gene_id"].fillna("").astype(str).str.strip())
            unique_gene_count = int(gene_keys.fillna("").astype(str).str.strip().replace("", pd.NA).dropna().nunique())
        dex_count = int(current["DEXSeq_significant"].fillna(False).astype(bool).sum()) if "DEXSeq_significant" in current.columns else 0
        dtu_count = int(current["DTU_significant"].fillna(False).astype(bool).sum()) if "DTU_significant" in current.columns else 0
        self.context_label.setText(
            f"{self.project_service.comparison_context_text(selected)} | "
            f"Unique genes: {unique_gene_count} | "
            f"Filtered rows: {len(current)} | "
            f"DEXSeq-significant rows: {dex_count} | "
            f"DTU-significant rows: {dtu_count}"
        )

    @staticmethod
    def _is_numeric_column(frame: pd.DataFrame, column: str) -> bool:
        return column in frame.columns and is_numeric_dtype(frame[column])

    @staticmethod
    def _is_boolean_column(frame: pd.DataFrame, column: str) -> bool:
        return column in frame.columns and is_bool_dtype(frame[column])
