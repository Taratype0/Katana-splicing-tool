from __future__ import annotations

import os
import shutil
from pathlib import Path

import pandas as pd
from PySide6.QtCore import Qt, QUrl
from PySide6.QtWidgets import (
    QFileDialog,
    QComboBox,
    QDoubleSpinBox,
    QFrame,
    QGridLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QSpinBox,
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

from app.pages.module_details import module_detail
from app.widgets.artifact_image_preview import ArtifactImagePreview
from app.widgets.flow_layout import FlowLayout
from app.widgets.preview_dialog import TextPreviewDialog
from src.services.project_service import ProjectService


class JutilsPage(QWidget):
    def __init__(self, project_service: ProjectService, view_mode: str = "browser") -> None:
        super().__init__()
        self.project_service = project_service
        self.view_mode = view_mode
        self.label = QLabel("Jutils output browser")
        self.label.setObjectName("pageTitle")
        self.subtitle = QLabel(
            "Sample-level PSI patterns, replicate QC and comparison similarity with explicit denominators."
        )
        self.subtitle.setObjectName("pageSubtitle")
        self.subtitle.setWordWrap(True)
        self.status_label = QLabel("")
        self.status_label.setWordWrap(True)
        self.summary_frame = QFrame()
        self.summary_frame.setObjectName("analysisSummary")
        summary_layout = QGridLayout(self.summary_frame)
        summary_layout.setContentsMargins(12, 10, 12, 10)
        summary_layout.setHorizontalSpacing(18)
        self.summary_values: dict[str, QLabel] = {}
        for column, (key, title) in enumerate(
            (
                ("comparison_count", "Comparisons"),
                ("event_count", "Events"),
                ("sample_count", "Samples"),
                ("pattern_cluster_count", "Pattern clusters"),
                ("pair_count", "Comparison pairs"),
                ("median_direction_similarity", "Median direction agreement"),
            )
        ):
            title_label = QLabel(title)
            title_label.setObjectName("summaryLabel")
            value_label = QLabel("—")
            value_label.setObjectName("summaryValue")
            summary_layout.addWidget(title_label, 0, column)
            summary_layout.addWidget(value_label, 1, column)
            self.summary_values[key] = value_label
        self.comparison_tabs = QTabWidget()
        self.comparison_tabs.setDocumentMode(True)
        self.comparison_tabs.setUsesScrollButtons(True)
        self.comparison_tabs.setMaximumHeight(46)
        self.comparison_tabs.currentChanged.connect(self.refresh)
        self.command_label = QTextEdit()
        self.command_label.setReadOnly(True)
        self.command_label.setMaximumHeight(90)
        self.output_dir_label = QLabel("")
        self.output_dir_label.setWordWrap(True)

        self.fdr_cutoff = QDoubleSpinBox()
        self.fdr_cutoff.setRange(0.0, 1.0)
        self.fdr_cutoff.setDecimals(4)
        self.fdr_cutoff.setSingleStep(0.01)
        self.fdr_cutoff.setValue(0.05)
        self.dpsi_cutoff = QDoubleSpinBox()
        self.dpsi_cutoff.setRange(0.0, 1.0)
        self.dpsi_cutoff.setDecimals(4)
        self.dpsi_cutoff.setSingleStep(0.01)
        self.dpsi_cutoff.setValue(0.05)
        self.top_events = QSpinBox()
        self.top_events.setRange(10, 5000)
        self.top_events.setValue(100)
        self.selection_mode = QComboBox()
        self.selection_mode.addItem("Supervised significant AS", "supervised")
        self.selection_mode.addItem("Unsupervised top variable PSI", "unsupervised")
        self.gene_input = QLineEdit()
        self.gene_input.setPlaceholderText("Candidate genes: comma, tab or space separated; empty = 5.3.2 selected genes")

        self.run_button = QPushButton(self._run_button_text())
        self.run_button.clicked.connect(self._run_jutils)
        self.preview_button = QPushButton("Preview Selected")
        self.preview_button.clicked.connect(self._preview_selected)
        self.detail_button = QPushButton("Detail")
        self.detail_button.clicked.connect(self._show_details)
        self.open_output_button = QPushButton("Open Output Folder")
        self.open_output_button.clicked.connect(self._open_output_folder)
        self.open_selected_button = QPushButton("Open Selected File")
        self.open_selected_button.clicked.connect(self._open_selected_file)
        self.download_button = QPushButton("Download Selected")
        self.download_button.clicked.connect(self._download_selected_file)
        for button in (
            self.run_button,
            self.preview_button,
            self.detail_button,
            self.open_output_button,
            self.open_selected_button,
            self.download_button,
        ):
            button.setMinimumHeight(36)
            button.setMinimumWidth(136)

        self.table = QTableWidget(0, 0)
        self.table.setSortingEnabled(True)
        self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.table.itemSelectionChanged.connect(self._preview_selected)
        self._frame = pd.DataFrame()

        self.preview_stack = QStackedWidget()
        self.image_scroll = ArtifactImagePreview("No preview selected.")
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

        top_buttons = FlowLayout(spacing=8)
        top_buttons.addWidget(self.run_button)
        top_buttons.addWidget(self.preview_button)
        top_buttons.addWidget(self.detail_button)
        top_buttons.addWidget(self.open_output_button)
        top_buttons.addWidget(self.open_selected_button)
        top_buttons.addWidget(self.download_button)

        controls = FlowLayout(spacing=8)
        controls.addWidget(QLabel("FDR/q <="))
        controls.addWidget(self.fdr_cutoff)
        controls.addWidget(QLabel("|dPSI| >="))
        controls.addWidget(self.dpsi_cutoff)
        controls.addWidget(QLabel("Mode"))
        controls.addWidget(self.selection_mode)
        controls.addWidget(QLabel("Top events"))
        controls.addWidget(self.top_events)
        if self.view_mode == "candidate_heatmap":
            controls.addWidget(QLabel("Genes"))
            controls.addWidget(self.gene_input)

        meta_layout = QVBoxLayout()
        meta_layout.addWidget(self.label)
        meta_layout.addWidget(self.subtitle)
        if self._is_native_mode() and self._native_report_mode() == "pattern":
            meta_layout.addWidget(self.summary_frame)
        meta_layout.addLayout(top_buttons)
        if self.view_mode in {"native_pattern", "candidate_heatmap"}:
            meta_layout.addLayout(controls)
        meta_layout.addWidget(self.comparison_tabs)
        meta_layout.addWidget(self.status_label)
        meta_layout.addWidget(self.output_dir_label)
        meta_layout.addWidget(QLabel("Command log"))
        meta_layout.addWidget(self.command_label)

        self.splitter = QSplitter(Qt.Orientation.Vertical)
        self.splitter.addWidget(self.table)
        self.splitter.addWidget(self.preview_stack)
        self.splitter.setStretchFactor(0, 2)
        self.splitter.setStretchFactor(1, 3)
        self.splitter.setSizes([240, 620])

        layout = QVBoxLayout(self)
        layout.addLayout(meta_layout)
        layout.addWidget(self.splitter, 1)

        cache_only_view = self.view_mode in {
            "native_heatmap",
            "native_pca",
            "native_similarity",
            "native_browser",
        }
        if cache_only_view:
            self.run_button.hide()
            self.fdr_cutoff.hide()
            self.dpsi_cutoff.hide()
            self.top_events.hide()
            self.selection_mode.hide()

    def _is_native_mode(self) -> bool:
        return self.view_mode in {
            "native_pattern",
            "native_heatmap",
            "native_pca",
            "native_similarity",
            "native_browser",
            "candidate_heatmap",
        }

    def _native_report_mode(self) -> str:
        return "candidate" if self.view_mode == "candidate_heatmap" else "pattern"

    def _run_button_text(self) -> str:
        if self.view_mode == "native_pattern":
            return "Run Native AS\nPattern Explorer"
        if self.view_mode == "candidate_heatmap":
            return "Run Candidate Gene\nAS Heatmap"
        if self._is_native_mode():
            return "Run Native Splicing\nPattern Analysis"
        return "Run Legacy\nJutils"

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

        if self._is_native_mode():
            native_mode = self._native_report_mode()
            report = self.project_service.native_jutils_report(native_mode)
            browser = self.project_service.build_native_jutils_output_browser(native_mode)
            if native_mode == "pattern":
                self._update_summary(self.project_service.native_pattern_summary(native_mode))
            if self.view_mode == "candidate_heatmap":
                self.label.setText("Candidate Gene AS Heatmap")
            elif self.view_mode == "native_heatmap":
                self.label.setText("Splicing Pattern Explorer — Sample PSI Heatmaps")
            elif self.view_mode == "native_pca":
                self.label.setText("Splicing Pattern Explorer — PCA & Replicate QC")
            elif self.view_mode == "native_similarity":
                self.label.setText("Splicing Pattern Explorer — Similarity & Overlap")
            elif self.view_mode == "native_browser":
                self.label.setText("Splicing Pattern Explorer — Reproducible Outputs")
            else:
                self.label.setText("Splicing Pattern Explorer")
        else:
            report = self.project_service.jutils_report()
            browser = self.project_service.build_jutils_output_browser()
        if not self._is_native_mode():
            if self.view_mode == "heatmap" and not browser.empty and "relative_path" in browser.columns:
                browser = browser.loc[
                    browser["relative_path"].astype(str).str.contains("heatmap|clustermap", case=False, na=False)
                ].copy()
                self.label.setText("Jutils heatmap outputs")
            elif self.view_mode == "pca" and not browser.empty and "relative_path" in browser.columns:
                browser = browser.loc[browser["relative_path"].astype(str).str.contains("pca", case=False, na=False)].copy()
                self.label.setText("Jutils PCA outputs")
            else:
                self.label.setText("Jutils output browser")
        browser = self._filter_native_browser(browser)
        self._populate_comparison_tabs(browser)
        selected_comparison = self._selected_comparison_id()
        if selected_comparison not in (None, "__all__") and not browser.empty and "comparison_id" in browser.columns:
            browser = browser.loc[browser["comparison_id"].astype(str) == str(selected_comparison)].copy()
        self._frame = browser
        self._populate_file_table(browser)
        if browser.empty:
            self.splitter.setSizes([90, 640])
        else:
            self.splitter.setSizes([240, 620])
        status = str(report.get("status", "idle"))
        output_dir = str(report.get("output_directory", ""))
        native_label = "Native Splicing Pattern Explorer"
        self.status_label.setText(f"{native_label if self._is_native_mode() else 'Jutils'} status: {status}")
        self.output_dir_label.setText(f"Output directory: {output_dir or 'not available'}")
        commands = report.get("command_lines", []) or []
        if self._is_native_mode():
            self.command_label.setPlainText(
                "\n".join(
                    [
                        "Katana-native implementation; no external jutils.py command is used.",
                        f"Active cutoff: FDR/q <= {self.fdr_cutoff.value():g}, |dPSI| >= {self.dpsi_cutoff.value():g}",
                        f"Selection mode: {self.selection_mode.currentText()}",
                        "Top events is applied per comparison in Katana native mode so the heatmap stays readable.",
                    ]
                )
            )
        else:
            self.command_label.setPlainText("\n".join(str(line) for line in commands) or "No command has been run yet.")

        if status in {"finished", "cache_hit"}:
            generated_files = report.get("generated_files", []) or []
            lines = [
                (
                    "Loaded cached Native Splicing Pattern Analysis."
                    if self._is_native_mode() and status == "cache_hit"
                    else "Native Splicing Pattern Analysis finished."
                    if self._is_native_mode()
                    else "Jutils finished."
                ),
                f"Output saved to: {output_dir or 'not available'}",
                "",
                self._native_explanation_text() if self._is_native_mode() else "",
                "",
                "Generated files:",
            ]
            lines.extend(f"- {row.get('relative_path', '')}" for row in generated_files[:200])
            if browser.empty and self.view_mode in {"heatmap", "pca"}:
                lines.append("")
                lines.append(f"No {self.view_mode.upper()} previewable files matched the current browser filter. Use 'Open Output Folder' to inspect raw outputs.")
            self.details.setPlainText("\n".join(lines))
        elif status == "failed":
            self.details.setPlainText(
                "\n".join(
                    [
                        "Native Splicing Pattern Explorer failed." if self._is_native_mode() else "Jutils failed.",
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
            if self._is_native_mode():
                self.details.setPlainText(
                    "\n".join(
                        [
                            "Native Splicing Pattern Explorer has not been run for this mode yet.",
                            self._native_explanation_text(),
                            "Opening this page only reads cached TSV/PNG files; it does not rebuild analysis.",
                        ]
                    )
                )
            else:
                self.details.setPlainText(
                    module_detail("5.6.jutils")
                    + "\n\nNo cached legacy Jutils output is loaded. Run Legacy Jutils only if you need compatibility with external jutils.py outputs."
                )
        if browser.empty and self.view_mode in {"heatmap", "pca"}:
            self.details.setPlainText(
                "\n".join(
                    [
                        f"No {self.view_mode.upper()} outputs are currently available for the selected comparison filter.",
                        f"Output directory: {output_dir or 'not available'}",
                        "Use 'Open Output Folder' to inspect raw files directly.",
                    ]
                )
            )
        self._preview_selected()

    def _update_summary(self, summary: dict[str, object]) -> None:
        for key, label in self.summary_values.items():
            value = summary.get(key)
            if key == "median_direction_similarity":
                try:
                    numeric = float(value)
                    label.setText(f"{numeric:.1%}" if pd.notna(numeric) else "—")
                except (TypeError, ValueError):
                    label.setText("—")
            else:
                try:
                    label.setText(f"{int(value):,}")
                except (TypeError, ValueError):
                    label.setText("—")

    def _native_explanation_text(self) -> str:
        if not self._is_native_mode():
            return ""
        if self.view_mode == "candidate_heatmap":
            scope = (
                "Purpose:\n"
                "Candidate-gene AS heatmap shows AS patterns for selected candidate genes. This uses genes typed in the gene box, or if empty, genes saved in 5.3.2 candidate gene selection. "
                "It is not a global AS overview. The first output is gene-level: every selected gene is retained, and genes without significant AS under the active cutoff are shown as no significant AS instead of disappearing."
            )
        else:
            scope = (
                "Purpose:\n"
                "Selected-comparison AS pattern overview is intended to be checked after basic AS landscape and before DEG+DAS integration. "
                "It is not candidate ranking and not causal modeling."
            )
        return "\n".join(
            [
                scope,
                "",
                "Inputs:",
                "- rMATS event tables from the selected comparisons.",
                "- Default Jutils-like filter: FDR/q <= 0.05 and |standardized dPSI| >= 0.05. The cutoffs are editable but optional.",
                "- dPSI follows Katana canonical AS direction inherited from Pairing.",
                "- Selection mode follows Jutils terminology: supervised = filter by p/q/dPSI; unsupervised = top most variable PSI events.",
                "",
                "Outputs:",
                "- Native AS Pattern Explorer writes per-comparison sample-level PSI heatmaps first, then cross-comparison similarity summaries.",
                "- Candidate Gene AS Heatmap writes a candidate_gene_as_summary table and summary heatmap first, then event-level PSI/dPSI heatmaps for the selected genes that have matching AS events.",
                "- The file browser is split into per-comparison tabs plus a cross_comparison tab, so individual comparison heatmap/PCA outputs do not disappear into one global output.",
                "- Per-comparison PSI heatmap/PCA are the Jutils-like outputs. Cross-comparison dPSI PCA is an extra QC view, not a replacement for legacy Jutils sample PCA.",
                "",
                "Start here for similarity:",
                "- plots/per_comparison/*_sample_psi_heatmap.png answers: within one comparison, what is the PSI pattern in the original sample/replicate order?",
                "- plots/per_comparison/*_sample_psi_heatmap_clustered.png answers: do samples/replicates cluster by PSI pattern when columns are also clustered?",
                "- plots/per_comparison/*_sample_psi_pca.png answers: within one comparison, do sample PSI profiles separate by the expected groups, matching the main purpose of legacy Jutils PCA?",
                "- native_as_similarity_heatmap.png answers: do selected comparisons have similar signed dPSI patterns?",
                "- native_as_direction_similarity_heatmap.png answers: among shared AS events, are the dPSI directions consistent?",
                "- pairwise_as_similarity.tsv is the numeric table behind those plots: shared event count, Jaccard overlap, direction similarity, opposite-direction fraction, Spearman dPSI correlation, weighted cosine similarity, and specificity.",
                "",
                "Interpretation:",
                "- Heatmaps show AS pattern structure; they are not candidate ranking by themselves.",
                "- PCA and similarity matrices are QC/pattern views for selected comparisons or samples.",
                "- PSI-pattern clusters group events with similar sample/replicate PSI trajectories. They are descriptive modules, not proof of a causal splicing program.",
                "",
                "Heatmap parameters:",
                "- Per-comparison sample PSI heatmap: rows = AS events from one comparison, columns = rMATS sample/replicate PSI values parsed from IncLevel1/IncLevel2.",
                "- Cross-comparison dPSI heatmap: rows = AS events, columns = selected comparisons, values = standardized dPSI.",
                "- Supervised mode: events pass FDR/q and |dPSI| filters; Katana then caps display to Top N per comparison for readability.",
                "- Unsupervised mode: events are ranked by PSI variability across samples within each comparison; Top N controls how many are displayed.",
                "- Row order uses Jutils-like hierarchical clustering (weighted/cityblock) when SciPy is available; otherwise it falls back to deterministic signal ordering.",
                "- The primary heatmap preserves selected comparison/sample order so direction and experiment meaning stay readable.",
                "- The *_clustered.png copy also clusters sample columns to check whether replicates group together by PSI pattern.",
                "",
                "PCA:",
                "- per_comparison_sample_PSI_PCA: Jutils-like PCA inside each comparison, using samples/replicates as points.",
                "- combined_sample_PSI_PCA: PCA across all selected comparison sample columns; this is a broader QC view.",
                "- PCA uses SVD after mean imputation and z-scoring. It is a QC/pattern visualization, not a statistical test.",
                "",
                "Similarity quantification:",
                "- Jaccard = overlap of significant AS events.",
                "- Signed Jaccard = overlap after including dPSI direction.",
                "- Direction similarity = fraction of shared events with the same dPSI sign.",
                "- Opposite direction fraction = fraction of shared events with reversed dPSI sign.",
                "- Spearman dPSI correlation = rank similarity of dPSI effects among shared events.",
                "- Weighted cosine similarity = overall signed dPSI pattern similarity across the event union.",
                "- Specificity = fraction of significant AS events unique to one comparison.",
                "",
                "Pattern clusters:",
                "- psi_pattern_clusters groups events by row-z-scored sample/replicate PSI trajectories using weighted-linkage hierarchical clustering.",
                "- The number of clusters is deliberately capped at eight; Katana does not cluster underdetermined two-value comparison dPSI vectors.",
                "- sample_psi_similarity.tsv reports Spearman PSI similarity together with shared_measured_events and a confidence label.",
                "- replicate_group_separation.tsv compares within-group replicate similarity with between-group similarity for descriptive QC.",
                "",
                "Run and cache behavior:",
                "- Opening the page reads cached TSV/PNG files only.",
                "- Run Native Splicing Pattern Analysis rebuilds only this native module cache.",
            ]
        )

    def _filter_native_browser(self, browser: pd.DataFrame) -> pd.DataFrame:
        if browser.empty or not self._is_native_mode() or self.view_mode in {"native_pattern", "native_browser", "candidate_heatmap"}:
            return browser
        if "result_group" not in browser.columns:
            return browser
        prefixes = {
            "native_heatmap": ("1.", "2.", "3."),
            "native_pca": ("4.",),
            "native_similarity": ("5.", "6."),
        }.get(self.view_mode)
        if not prefixes:
            return browser
        return browser.loc[browser["result_group"].astype(str).str.startswith(prefixes)].copy()

    def _show_details(self) -> None:
        text = self.details.toPlainText().strip() or self._native_explanation_text() or module_detail("5.6.jutils")
        dialog = TextPreviewDialog(self.label.text() + " - Detail", text, self)
        dialog.exec()

    def _populate_comparison_tabs(self, frame: pd.DataFrame) -> None:
        current = self._selected_comparison_id()
        comparisons = []
        if "comparison_id" in frame.columns:
            comparisons = list(dict.fromkeys(frame["comparison_id"].dropna().astype(str).tolist()))
        elif self.project_service.current_project is not None:
            comparisons = [item.comparison_id for item in self.project_service.selected_comparisons_for_display()]
        self.comparison_tabs.blockSignals(True)
        self.comparison_tabs.clear()
        self.comparison_tabs.addTab(QWidget(), "All comparisons")
        self.comparison_tabs.setTabToolTip(0, "__all__")
        for index, comparison_id in enumerate(comparisons, start=1):
            label = comparison_id.replace("per_comparison:", "")
            self.comparison_tabs.addTab(QWidget(), label)
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
            self.table.setSortingEnabled(False)
            self.table.setColumnCount(1)
            self.table.setHorizontalHeaderLabels(["Status"])
            self.table.setRowCount(1)
            message = (
                "No cached Native AS Pattern Explorer output. Click Run to generate PSI heatmap/PCA/similarity outputs."
                if self._is_native_mode()
                else "No cached Jutils output."
            )
            item = QTableWidgetItem(message)
            item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.table.setItem(0, 0, item)
            self.table.resizeColumnsToContents()
            self.table.setSortingEnabled(False)
            return
        columns = [
            column
            for column in [
                "result_group",
                "comparison_id",
                "relative_path",
                "file_type",
                "file_size",
                "last_modified",
                "preview_kind",
            ]
            if column in frame.columns
        ]
        self.table.setSortingEnabled(False)
        self.table.setColumnCount(len(columns))
        self.table.setHorizontalHeaderLabels(columns)
        self.table.setRowCount(len(frame))
        for row_idx, (_, row) in enumerate(frame[columns].iterrows()):
            for col_idx, column in enumerate(columns):
                item = QTableWidgetItem("" if pd.isna(row[column]) else str(row[column]))
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(row_idx, col_idx, item)
        self.table.setSortingEnabled(True)
        self.table.resizeColumnsToContents()
        if self.table.rowCount() > 0:
            self.table.selectRow(0)

    def _run_jutils(self) -> None:
        active_label = "Native Splicing Pattern Explorer" if self._is_native_mode() else "Jutils"
        self.status_label.setText(f"{active_label} status: running")
        self.details.setPlainText(f"{active_label} is running...")
        try:
            if self.view_mode in {"native_pattern", "native_heatmap", "native_pca", "native_similarity", "native_browser"}:
                self.project_service.run_native_jutils_pattern_explorer(
                    fdr_cutoff=float(self.fdr_cutoff.value()),
                    abs_dpsi_cutoff=float(self.dpsi_cutoff.value()),
                    top_events=int(self.top_events.value()),
                    selection_mode=str(self.selection_mode.currentData() or "supervised"),
                )
            elif self.view_mode == "candidate_heatmap":
                self.project_service.run_native_jutils_candidate_gene_heatmap(
                    gene_text=self.gene_input.text(),
                    fdr_cutoff=float(self.fdr_cutoff.value()),
                    abs_dpsi_cutoff=float(self.dpsi_cutoff.value()),
                    top_events=int(self.top_events.value()),
                    selection_mode=str(self.selection_mode.currentData() or "supervised"),
                )
            else:
                self.project_service.run_jutils_pipeline()
        except Exception as exc:
            self.details.setPlainText(f"{active_label} failed:\n{exc}")
            self.status_label.setText(f"{active_label} status: failed")
            return
        self.refresh()

    def _selected_row(self) -> pd.Series | None:
        if self._frame.empty:
            return None
        row_idx = self.table.currentRow()
        if row_idx < 0:
            return None
        relative_col = None
        for col_idx in range(self.table.columnCount()):
            header = self.table.horizontalHeaderItem(col_idx)
            if header and header.text() == "relative_path":
                relative_col = col_idx
                break
        if relative_col is not None:
            item = self.table.item(row_idx, relative_col)
            if item is not None:
                matches = self._frame.loc[self._frame["relative_path"].astype(str) == item.text()]
                if not matches.empty:
                    return matches.iloc[0]
        if row_idx >= len(self._frame):
            return None
        return self._frame.iloc[row_idx]

    def _preview_selected(self) -> None:
        row = self._selected_row()
        if row is None:
            self.preview_stack.setCurrentWidget(self.text_preview)
            if self._is_native_mode() and self._frame.empty:
                self.text_preview.setPlainText(
                    "\n".join(
                        [
                            "No cached Native AS Pattern Explorer output is available.",
                            "",
                            "This page only reads existing cache when opened; it does not auto-run analysis.",
                            "Click Run Native AS Pattern Explorer to generate:",
                            "- per-comparison sample PSI heatmap",
                            "- per-comparison sample PSI PCA",
                            "- cross-comparison AS similarity tables/plots",
                            "",
                            "After results exist, the file list and PNG/TSV preview will appear here.",
                        ]
                    )
                )
            else:
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
                self.preview_stack.setCurrentWidget(self.image_scroll)
                self.image_scroll.set_image(path)
        elif preview_kind == "html":
            self.html_browser.setSource(QUrl.fromLocalFile(str(path)))
            self.preview_stack.setCurrentWidget(self.html_browser)
        elif preview_kind in {"table", "excel"}:
            self._load_table_preview(path)
            self.preview_stack.setCurrentWidget(self.table_preview)
        elif preview_kind == "pdf":
            self.preview_stack.setCurrentWidget(self.text_preview)
            self.text_preview.setPlainText(
                "\n".join(
                    [
                        "PDF preview is not embedded in Katana.",
                        "Use 'Open Selected File' or 'Download Selected' to inspect it directly.",
                        "",
                        f"Path: {path}",
                    ]
                )
            )
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
        report = self.project_service.native_jutils_report(self._native_report_mode()) if self._is_native_mode() else self.project_service.jutils_report()
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
            "Download Selected Native Output" if self._is_native_mode() else "Download Selected Jutils File",
            path.name,
            "All Files (*.*)",
        )
        if not output_path:
            return
        shutil.copy2(path, output_path)
