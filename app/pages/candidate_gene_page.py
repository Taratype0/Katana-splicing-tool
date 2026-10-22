from __future__ import annotations

import os
import shutil

import pandas as pd
from pandas.api.types import is_bool_dtype, is_numeric_dtype
from PySide6.QtCore import QObject, QThread, QTimer, Signal, Qt
from PySide6.QtWidgets import (
    QDialog,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QProgressBar,
    QPushButton,
    QComboBox,
    QDoubleSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from app.pages.module_details import module_detail
from src.services.project_service import ProjectService


class CandidateGenerationWorker(QObject):
    finished = Signal(object)
    failed = Signal(str)

    def __init__(self, project_service: ProjectService) -> None:
        super().__init__()
        self.project_service = project_service

    def run(self) -> None:
        try:
            frame = self.project_service.preview_candidate_gene_screening(
                allow_generate=True,
                force_rebuild=False,
            )
        except Exception as exc:
            self.failed.emit(str(exc))
            return
        self.finished.emit(frame)


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
        self.search.textChanged.connect(self._schedule_refresh)
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
        self.pattern_filter = QComboBox()
        self.pattern_filter.addItems(
            [
                "All",
                "shared_all",
                "shared_subset",
                "comparison_specific",
                "gained",
                "lost",
                "direction_reversed",
                "tier_changed",
                "evidence_class_changed",
                "DEG_only_in_one_condition",
                "splicing_only_in_one_condition",
            ]
        )
        self.pattern_filter.currentTextChanged.connect(self.refresh)
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
        self.detail_button = QPushButton("Detail")
        self.detail_button.clicked.connect(self._show_detail_dialog)
        self.open_button = QPushButton("Open output folder")
        self.open_button.clicked.connect(self.open_output_folder)
        self.download_button = QPushButton("Download file")
        self.download_button.clicked.connect(self.download_current_file)
        self.export_button = QPushButton("Export Current View")
        self.export_button.clicked.connect(self.export_current_view)
        self.progress = QProgressBar()
        self.progress.setRange(0, 0)
        self.progress.setVisible(False)
        self.progress.setMaximumWidth(170)
        self.table = QTableWidget(0, 0)
        self.table.setSortingEnabled(True)
        self._current_frame = pd.DataFrame()
        self._current_source_path = None
        self._candidate_thread: QThread | None = None
        self._candidate_worker: CandidateGenerationWorker | None = None
        self._max_display_rows = 1000
        self._search_timer = QTimer(self)
        self._search_timer.setSingleShot(True)
        self._search_timer.setInterval(300)
        self._search_timer.timeout.connect(self.refresh)

        top = QHBoxLayout()
        top.addWidget(self.label)
        top.addWidget(self.search)
        top.addWidget(self.comparison_filter)
        top.addWidget(self.class_filter)
        top.addWidget(self.evidence_filter)
        top.addWidget(self.direction_filter)
        top.addWidget(self.pattern_filter)
        top.addWidget(self.min_abs_log2fc)
        top.addWidget(self.min_abs_dpsi)
        top.addWidget(self.generate_button)
        top.addWidget(self.detail_button)
        top.addWidget(self.progress)
        top.addWidget(self.open_button)
        top.addWidget(self.download_button)
        top.addWidget(self.export_button)

        layout = QVBoxLayout(self)
        layout.addLayout(top)
        layout.addWidget(self.context_label)
        layout.addWidget(self.comparison_tabs)
        layout.addWidget(self.path_label)
        layout.addWidget(self.table)

    def _show_detail_dialog(self) -> None:
        text = module_detail("5.3.ranking")
        if self.mode == "tiers":
            text += "\n\nTier view note:\nThis page uses the same unique gene-level ranking source, filtered by candidate_tier."
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
                "cross_candidate_pattern",
                "cross_pattern_reason",
                "cross_candidate_comparisons",
                "cross_n_comparisons_candidate",
                "cross_best_overall_rank",
                "cross_comparison_effect_summary",
                "shared_candidate",
                "specific_candidate",
                "gained_candidate",
                "lost_candidate",
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
                "ranking_score",
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
                "cross_candidate_pattern",
                "cross_pattern_reason",
                "cross_candidate_comparisons",
                "cross_n_comparisons_candidate",
                "cross_best_overall_rank",
                "cross_comparison_effect_summary",
                "shared_candidate",
                "specific_candidate",
                "gained_candidate",
                "lost_candidate",
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
                "ranking_score",
            ]
        available_columns = [column for column in columns if column in filtered.columns]
        self.table.setSortingEnabled(False)
        self.table.setColumnCount(len(available_columns))
        self.table.setHorizontalHeaderLabels(available_columns)
        display_frame = filtered.head(self._max_display_rows).copy()
        self.table.setRowCount(len(display_frame))
        for row_idx, (_, row) in enumerate(display_frame[available_columns].iterrows()):
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
        shown = len(display_frame)
        suffix = f"showing first {shown} of {len(filtered)} rows" if len(filtered) > shown else f"{len(filtered)} rows"
        self.label.setText(f"{label_prefix} ({suffix})")

    def _schedule_refresh(self) -> None:
        self._search_timer.start()

    def _generate_candidate_results(self) -> None:
        if self._candidate_thread is not None:
            self.context_label.setText("Candidate ranking is already running. Please wait for this run to finish.")
            return
        self.context_label.setText("Candidate ranking is running... Existing cache will be reused when valid.")
        self.progress.setVisible(True)
        self.generate_button.setEnabled(False)
        self._candidate_thread = QThread(self)
        self._candidate_worker = CandidateGenerationWorker(self.project_service)
        self._candidate_worker.moveToThread(self._candidate_thread)
        self._candidate_thread.started.connect(self._candidate_worker.run)
        self._candidate_worker.finished.connect(self._on_candidate_generation_finished)
        self._candidate_worker.failed.connect(self._on_candidate_generation_failed)
        self._candidate_worker.finished.connect(self._candidate_thread.quit)
        self._candidate_worker.failed.connect(self._candidate_thread.quit)
        self._candidate_thread.finished.connect(self._cleanup_candidate_generation)
        self._candidate_thread.start()

    def _on_candidate_generation_finished(self, frame: pd.DataFrame) -> None:
        if frame is None or frame.empty:
            if self.project_service.run_state.results.tx_splicing_gene_table.empty:
                self.context_label.setText(
                    "Missing input: transcript + splicing integration not run for this project. Please run 5.2.1 first."
                )
            else:
                self.context_label.setText(
                    "Candidate ranking generation finished, but no candidate rows were produced for the current project."
                )
        else:
            self.context_label.setText(f"Candidate ranking ready: {len(frame)} unique comparison-gene rows loaded.")
        self.refresh()

    def _on_candidate_generation_failed(self, message: str) -> None:
        self.context_label.setText(f"Candidate ranking failed: {message}")

    def _cleanup_candidate_generation(self) -> None:
        self.progress.setVisible(False)
        self.generate_button.setEnabled(True)
        if self._candidate_worker is not None:
            self._candidate_worker.deleteLater()
        if self._candidate_thread is not None:
            self._candidate_thread.deleteLater()
        self._candidate_worker = None
        self._candidate_thread = None

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

        selected_pattern = self.pattern_filter.currentText()
        if selected_pattern != "All" and "cross_candidate_pattern" in filtered.columns:
            filtered = filtered.loc[filtered["cross_candidate_pattern"].astype(str) == selected_pattern].copy()

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
        for column in [
            "comparison_id",
            "gene_symbol",
            "gene_id",
            "candidate_tier",
            "evidence_class",
            "direction_class",
            "cross_candidate_pattern",
            "cross_candidate_comparisons",
            "cross_comparison_effect_summary",
            "candidate_reason",
        ]:
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
        dex_count = int(current["DEXSeq_significant"].map(self._safe_bool).sum()) if "DEXSeq_significant" in current.columns else 0
        dtu_count = int(current["DTU_significant"].map(self._safe_bool).sum()) if "DTU_significant" in current.columns else 0
        self.context_label.setText(
            f"{self.project_service.comparison_context_text(selected)} | "
            f"Unique genes: {unique_gene_count} | "
            f"Filtered rows: {len(current)} | "
            f"DEXSeq-significant rows: {dex_count} | "
            f"DTU-significant rows: {dtu_count}"
        )

    @staticmethod
    def _safe_bool(value) -> bool:
        if value is None:
            return False
        try:
            if pd.isna(value):
                return False
        except (TypeError, ValueError):
            pass
        if isinstance(value, str):
            return value.strip().lower() in {"1", "true", "yes", "y", "t"}
        return bool(value)

    @staticmethod
    def _is_numeric_column(frame: pd.DataFrame, column: str) -> bool:
        return column in frame.columns and is_numeric_dtype(frame[column])

    @staticmethod
    def _is_boolean_column(frame: pd.DataFrame, column: str) -> bool:
        return column in frame.columns and is_bool_dtype(frame[column])
