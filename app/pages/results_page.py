from __future__ import annotations

import os
from pathlib import Path
import re
import shutil
import textwrap

import numpy as np
import pandas as pd
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from matplotlib.patches import Circle, Ellipse, Patch
from PySide6.QtCore import QEvent, Qt, QTimer
from PySide6.QtGui import QGuiApplication, QKeySequence, QShortcut
from PySide6.QtWidgets import (
    QAbstractItemView,
    QAbstractScrollArea,
    QComboBox,
    QCheckBox,
    QDialog,
    QDoubleSpinBox,
    QFileDialog,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QScrollArea,
    QSplitter,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from app.pages.module_details import module_detail
from app.pages.cross_comparison_controller import CrossComparisonController
from app.widgets.flow_layout import FlowLayout
from src.services.project_service import ProjectService


class AnalysisHomePage(QWidget):
    def __init__(self, on_open_branch=None) -> None:
        super().__init__()
        self.on_open_branch = on_open_branch
        title = QLabel("Analysis navigation")
        title.setStyleSheet("font-size: 18px; font-weight: 700;")
        description = QLabel(
            "Open only the analysis section you want to run. Each branch is independent: run Landscape for single-comparison "
            "splicing summaries, Group Comparison for comparison-vs-comparison differences, Transcript + Splicing for DEG-vs-dPSI integration, "
            "Mechanism Support for secondary support layers, and Candidate Hook for ranked candidate follow-up."
        )
        description.setWordWrap(True)

        grid = QGridLayout()
        grid.setHorizontalSpacing(14)
        grid.setVerticalSpacing(14)

        sections = [
            ("Landscape", "Single-comparison splicing views and upstream stable-AS filtering.", "landscape"),
            ("Group Comparison", "Comparison-vs-comparison results: Excel preview, program heatmap and direction/strength split.", "group"),
            ("Transcript + Splicing", "Interactive DEG-vs-splicing scatter per comparison with gene drilldown.", "tx"),
            ("Mechanism Support", "DEXSeq / DTU / SUPPA / transcript-level support summary with readable gene symbols.", "mechanism"),
            ("Candidate Gene Screening", "Comparison-first candidate review: heatmaps, cards and isoform follow-up all driven from the integrated candidate table.", "candidate"),
        ]
        for index, (label, body, route) in enumerate(sections):
            card = QWidget()
            card_layout = QVBoxLayout(card)
            card_layout.setContentsMargins(12, 12, 12, 12)
            heading = QLabel(label)
            heading.setStyleSheet("font-size: 15px; font-weight: 700;")
            text = QLabel(body)
            text.setWordWrap(True)
            button = QPushButton(f"Open {label}")
            button.clicked.connect(lambda _checked=False, branch=route: self._open_branch(branch))
            card_layout.addWidget(heading)
            card_layout.addWidget(text)
            card_layout.addStretch(1)
            card_layout.addWidget(button)
            card.setStyleSheet("QWidget { border: 1px solid #3a3a3a; border-radius: 8px; }")
            grid.addWidget(card, index // 2, index % 2)

        layout = QVBoxLayout(self)
        layout.addWidget(title)
        layout.addWidget(description)
        layout.addLayout(grid)
        layout.addStretch(1)

    def _open_branch(self, branch: str) -> None:
        if callable(self.on_open_branch):
            self.on_open_branch(branch)


class ResultsPage(QWidget):
    def __init__(
        self,
        project_service: ProjectService,
        on_run_modules=None,
        fixed_branch: str | None = None,
        fixed_subbranch: str | None = None,
    ) -> None:
        super().__init__()
        self.project_service = project_service
        self.on_run_modules = on_run_modules
        self.fixed_branch = fixed_branch
        self.fixed_subbranch = fixed_subbranch

        self.label = QLabel("Analysis workbench" if fixed_branch is None else "Analysis")
        self.context_header = QLabel("")
        self.context_header.setWordWrap(True)
        self.explanation = QLabel("")
        self.explanation.setWordWrap(True)

        self.tabs = QTabWidget()
        self.landscape_compare_tabs = QTabWidget()
        self.landscape_tabs = QTabWidget()
        self.landscape_tabs.addTab(QWidget(), "Landscape")
        self.landscape_tabs.addTab(QWidget(), "Bar")
        self.landscape_tabs.addTab(QWidget(), "Pie")
        self.landscape_container = QWidget()
        landscape_layout = QVBoxLayout(self.landscape_container)
        landscape_layout.setContentsMargins(0, 0, 0, 0)
        landscape_layout.addWidget(self.landscape_compare_tabs)
        landscape_layout.addWidget(self.landscape_tabs)
        self.group_tabs = QTabWidget()
        self.group_mode_tabs = QTabWidget()
        self.group_mode_tabs.addTab(QWidget(), "Excel Preview")
        self.group_mode_tabs.addTab(QWidget(), "Program Heatmap")
        self.group_mode_tabs.addTab(QWidget(), "Direction vs Strength")
        self.tx_tabs = QTabWidget()
        self.mechanism_compare_tabs = QTabWidget()
        self.mechanism_tabs = QTabWidget()
        self.mechanism_tabs.addTab(QWidget(), "Support Summary")
        self.mechanism_container = QWidget()
        mechanism_layout = QVBoxLayout(self.mechanism_container)
        mechanism_layout.setContentsMargins(0, 0, 0, 0)
        mechanism_layout.addWidget(self.mechanism_compare_tabs)
        mechanism_layout.addWidget(self.mechanism_tabs)
        self.candidate_tabs = QTabWidget()
        self.candidate_mode_tabs = QTabWidget()
        self.candidate_mode_tabs.addTab(QWidget(), "Table")
        self.candidate_mode_tabs.addTab(QWidget(), "Evidence Heatmap")
        self.candidate_mode_tabs.addTab(QWidget(), "Program Heatmap")
        self.candidate_mode_tabs.addTab(QWidget(), "Candidate Cards + Event Follow-up")
        self.candidate_mode_tabs.addTab(QWidget(), "Isoform Follow-up")

        self.group_container = QWidget()
        group_layout = QVBoxLayout(self.group_container)
        group_layout.setContentsMargins(0, 0, 0, 0)
        group_layout.addWidget(self.group_tabs)
        group_layout.addWidget(self.group_mode_tabs)

        self.candidate_container = QWidget()
        candidate_layout = QVBoxLayout(self.candidate_container)
        candidate_layout.setContentsMargins(0, 0, 0, 0)
        candidate_layout.addWidget(self.candidate_tabs)
        candidate_layout.addWidget(self.candidate_mode_tabs)

        self.tabs.addTab(self.landscape_container, "Landscape")
        self.tabs.addTab(self.group_container, "Group Comparison")
        self.tabs.addTab(self.tx_tabs, "Transcript + Splicing")
        self.tabs.addTab(self.mechanism_container, "Mechanism Support")
        self.tabs.addTab(self.candidate_container, "Candidate Screening")

        self.candidate_source_combo = QComboBox()
        self.candidate_source_combo.addItem("Candidate selection", "selection")
        self.candidate_source_combo.addItem("Custom gene list", "custom")
        self.candidate_source_combo.addItem("Manual shortlist genes", "shortlist")
        self.candidate_source_combo.currentIndexChanged.connect(self.refresh)
        self.candidate_top_n = QSpinBox()
        self.candidate_top_n.setRange(1, 500)
        self.candidate_top_n.setValue(20)
        self.candidate_top_n.valueChanged.connect(self.refresh)
        self.candidate_custom_genes = QLineEdit()
        self.candidate_custom_genes.setPlaceholderText("Custom genes: comma, tab or space separated")
        self.candidate_custom_genes.textChanged.connect(self.refresh)
        self.candidate_controls = QWidget()
        candidate_controls_layout = QHBoxLayout(self.candidate_controls)
        candidate_controls_layout.setContentsMargins(0, 0, 0, 0)
        candidate_controls_layout.addWidget(QLabel("Gene source"))
        candidate_controls_layout.addWidget(self.candidate_source_combo)
        candidate_controls_layout.addWidget(QLabel("Top N"))
        candidate_controls_layout.addWidget(self.candidate_top_n)
        candidate_controls_layout.addWidget(QLabel("Custom genes"))
        candidate_controls_layout.addWidget(self.candidate_custom_genes, 1)
        self.candidate_controls.setVisible(False)
        self.card_gene_combo = QComboBox()
        self.card_gene_combo.setMinimumWidth(320)
        self.card_gene_combo.currentIndexChanged.connect(self._on_candidate_card_gene_changed)
        self.card_open_sashimi_button = QPushButton("Open selected event in Sashimi")
        self.card_open_sashimi_button.clicked.connect(self._open_candidate_card_event_in_sashimi)
        self.candidate_card_controls = QWidget()
        candidate_card_layout = QHBoxLayout(self.candidate_card_controls)
        candidate_card_layout.setContentsMargins(0, 0, 0, 0)
        candidate_card_layout.addWidget(QLabel("Card gene"))
        candidate_card_layout.addWidget(self.card_gene_combo, 1)
        candidate_card_layout.addWidget(self.card_open_sashimi_button)
        self.candidate_card_controls.setVisible(False)

        thresholds = self.project_service.current_thresholds
        self.landscape_stable_min_reads = QSpinBox()
        self.landscape_stable_min_reads.setRange(0, 100000)
        self.landscape_stable_min_reads.setValue(int(getattr(thresholds, "stable_min_junction_reads", 20)))
        self.landscape_stable_min_reads.setMinimumWidth(88)
        self.landscape_stable_min_reads.setToolTip("Minimum mean IJC+SJC reads per group for a stable AS event")
        self.landscape_stable_min_reps = QSpinBox()
        self.landscape_stable_min_reps.setRange(1, 20)
        self.landscape_stable_min_reps.setValue(int(getattr(thresholds, "stable_min_replicates_per_group", 2)))
        self.landscape_stable_min_reps.setMinimumWidth(88)
        self.landscape_stable_min_reps.setToolTip("Minimum PSI/read replicates required in each group")
        self.landscape_stable_consistency = QDoubleSpinBox()
        self.landscape_stable_consistency.setRange(0.0, 1.0)
        self.landscape_stable_consistency.setSingleStep(0.01)
        self.landscape_stable_consistency.setDecimals(2)
        self.landscape_stable_consistency.setValue(float(getattr(thresholds, "stable_min_pairwise_consistency", 0.67)))
        self.landscape_stable_consistency.setMinimumWidth(88)
        self.landscape_stable_consistency.setToolTip("Minimum pairwise replicate PSI direction consistency")
        self.landscape_stable_controls = QWidget()
        landscape_stable_layout = QHBoxLayout(self.landscape_stable_controls)
        landscape_stable_layout.setContentsMargins(0, 0, 0, 0)
        landscape_stable_layout.addWidget(QLabel("Stable AS filter: mean IJC+SJC reads >="))
        landscape_stable_layout.addWidget(self.landscape_stable_min_reads)
        landscape_stable_layout.addWidget(QLabel("min reps/group"))
        landscape_stable_layout.addWidget(self.landscape_stable_min_reps)
        landscape_stable_layout.addWidget(QLabel("direction consistency >="))
        landscape_stable_layout.addWidget(self.landscape_stable_consistency)
        landscape_stable_layout.addWidget(QLabel("Applied when 5.1 is rerun; downstream 5.2/5.3 inherit stable AS events."))
        landscape_stable_layout.addStretch(1)
        self.landscape_stable_controls.setVisible(False)

        self.group_pattern_sig_cutoff = QDoubleSpinBox()
        self.group_pattern_sig_cutoff.setRange(0.0, 5.0)
        self.group_pattern_sig_cutoff.setSingleStep(0.01)
        self.group_pattern_sig_cutoff.setDecimals(3)
        self.group_pattern_sig_cutoff.setValue(0.10)
        self.group_pattern_sig_cutoff.setMinimumWidth(88)
        self.group_pattern_sig_cutoff.setToolTip("abs(dPSI) cutoff for considering an AS event significant in this cross-comparison section")
        self.group_pattern_large_delta = QDoubleSpinBox()
        self.group_pattern_large_delta.setRange(0.0, 5.0)
        self.group_pattern_large_delta.setSingleStep(0.01)
        self.group_pattern_large_delta.setDecimals(3)
        self.group_pattern_large_delta.setValue(0.10)
        self.group_pattern_large_delta.setMinimumWidth(88)
        self.group_pattern_large_delta.setToolTip("large cross-comparison delta dPSI cutoff")
        self.group_pattern_view = QComboBox()
        self.group_pattern_view.addItem("Table: same-direction large dPSI delta", "same_direction_large_delta")
        self.group_pattern_view.addItem("Table: opposite-direction shared AS events", "opposite_direction")
        self.group_pattern_view.addItem("Table: all shared significant AS events", "all_shared")
        self.group_pattern_view.addItem("Legacy table: A-only / B-only significant AS events", "condition_specific")
        self.group_pattern_view.setMinimumWidth(360)
        self.group_pattern_view.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToContents)
        self.group_pattern_view.currentIndexChanged.connect(self.refresh)
        self.group_pattern_apply = QPushButton("Apply / Recalculate")
        self.group_pattern_apply.clicked.connect(self._apply_group_pattern_recalculate)
        self.group_pattern_controls = QWidget()
        group_pattern_layout = QHBoxLayout(self.group_pattern_controls)
        group_pattern_layout.setContentsMargins(0, 0, 0, 0)
        group_pattern_layout.addWidget(QLabel("abs dPSI cutoff"))
        group_pattern_layout.addWidget(self.group_pattern_sig_cutoff)
        group_pattern_layout.addWidget(QLabel("large delta cutoff"))
        group_pattern_layout.addWidget(self.group_pattern_large_delta)
        group_pattern_layout.addWidget(QLabel("Result table"))
        group_pattern_layout.addWidget(self.group_pattern_view)
        group_pattern_layout.addWidget(self.group_pattern_apply)
        group_pattern_layout.addStretch(1)
        self.group_pattern_controls.setVisible(False)

        self.group_interaction_sig_cutoff = QDoubleSpinBox()
        self.group_interaction_sig_cutoff.setRange(0.0, 5.0)
        self.group_interaction_sig_cutoff.setSingleStep(0.01)
        self.group_interaction_sig_cutoff.setDecimals(3)
        self.group_interaction_sig_cutoff.setValue(0.10)
        self.group_interaction_sig_cutoff.setMinimumWidth(88)
        self.group_interaction_delta_cutoff = QDoubleSpinBox()
        self.group_interaction_delta_cutoff.setRange(0.0, 5.0)
        self.group_interaction_delta_cutoff.setSingleStep(0.01)
        self.group_interaction_delta_cutoff.setDecimals(3)
        self.group_interaction_delta_cutoff.setValue(0.10)
        self.group_interaction_delta_cutoff.setMinimumWidth(88)
        self.group_interaction_min_reads = QSpinBox()
        self.group_interaction_min_reads.setRange(0, 100000)
        self.group_interaction_min_reads.setValue(20)
        self.group_interaction_min_reads.setMinimumWidth(88)
        self.group_interaction_view = QComboBox()
        self.group_interaction_view.addItem("Unified event result table", "unified_event_table")
        self.group_interaction_view.addItem("All-pair summary table", "batch_summary")
        self.group_interaction_view.setMinimumWidth(260)
        self.group_interaction_view.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToContents)
        self.group_interaction_view.currentIndexChanged.connect(self.refresh)
        self.group_interaction_apply = QPushButton("Apply / Recalculate")
        self.group_interaction_apply.clicked.connect(self._apply_group_interaction_recalculate)
        self.group_interaction_controls = QWidget()
        group_interaction_layout = QHBoxLayout(self.group_interaction_controls)
        group_interaction_layout.setContentsMargins(0, 0, 0, 0)
        group_interaction_layout.addWidget(QLabel("abs dPSI cutoff"))
        group_interaction_layout.addWidget(self.group_interaction_sig_cutoff)
        group_interaction_layout.addWidget(QLabel("large delta cutoff"))
        group_interaction_layout.addWidget(self.group_interaction_delta_cutoff)
        group_interaction_layout.addWidget(QLabel("min reads"))
        group_interaction_layout.addWidget(self.group_interaction_min_reads)
        group_interaction_layout.addWidget(QLabel("Result table"))
        group_interaction_layout.addWidget(self.group_interaction_view)
        group_interaction_layout.addWidget(self.group_interaction_apply)
        group_interaction_layout.addStretch(1)
        self.group_interaction_controls.setVisible(False)

        self.group_venn_sig_cutoff = QDoubleSpinBox()
        self.group_venn_sig_cutoff.setRange(0.0, 1.0)
        self.group_venn_sig_cutoff.setSingleStep(0.01)
        self.group_venn_sig_cutoff.setDecimals(3)
        self.group_venn_sig_cutoff.setValue(0.10)
        self.group_venn_sig_cutoff.setMinimumWidth(88)
        self.group_venn_sig_cutoff.valueChanged.connect(self.refresh)
        self.group_venn_min_reads = QSpinBox()
        self.group_venn_min_reads.setRange(0, 100000)
        self.group_venn_min_reads.setValue(20)
        self.group_venn_min_reads.setMinimumWidth(88)
        self.group_venn_min_reads.valueChanged.connect(self.refresh)
        self.group_venn_filter = QComboBox()
        self.group_venn_filter.addItem("All significant AS events (A-only / shared / B-only)", "all_significant")
        self.group_venn_filter.addItem("A-only/B-only: counterpart has >= min reads", "supported_only")
        self.group_venn_filter.addItem("A-only/B-only: counterpart low/absent reads", "low_power_only")
        self.group_venn_filter.addItem("A-only/B-only: counterpart event not found", "counterpart_not_found_only")
        self.group_venn_filter.setMinimumWidth(520)
        self.group_venn_filter.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToContents)
        self.group_venn_filter.currentIndexChanged.connect(self.refresh)
        self.group_venn_controls = QWidget()
        group_venn_layout = QHBoxLayout(self.group_venn_controls)
        group_venn_layout.setContentsMargins(0, 0, 0, 0)
        group_venn_layout.addWidget(QLabel("Significant AS: FDR <= 0.05, |dPSI| >="))
        group_venn_layout.addWidget(self.group_venn_sig_cutoff)
        group_venn_layout.addWidget(QLabel("counterpart min reads"))
        group_venn_layout.addWidget(self.group_venn_min_reads)
        group_venn_layout.addWidget(QLabel("Exclusive-side filter"))
        group_venn_layout.addWidget(self.group_venn_filter)
        group_venn_layout.addStretch(1)
        self.group_venn_controls.setVisible(False)
        self._group_venn_cache: dict[str, object] = {}

        self.group_venn_multi_sig_cutoff = QDoubleSpinBox()
        self.group_venn_multi_sig_cutoff.setRange(0.0, 1.0)
        self.group_venn_multi_sig_cutoff.setSingleStep(0.01)
        self.group_venn_multi_sig_cutoff.setDecimals(3)
        self.group_venn_multi_sig_cutoff.setValue(0.10)
        self.group_venn_multi_sig_cutoff.setMinimumWidth(88)
        self.group_venn_multi_min_reads = QSpinBox()
        self.group_venn_multi_min_reads.setRange(0, 100000)
        self.group_venn_multi_min_reads.setValue(20)
        self.group_venn_multi_min_reads.setMinimumWidth(88)
        self._syncing_group_venn_multi_checks = False
        self.group_venn_multi_checkboxes: dict[str, QCheckBox] = {}
        self.group_venn_multi_comparison_panel = QWidget()
        self.group_venn_multi_comparison_grid = QGridLayout(self.group_venn_multi_comparison_panel)
        self.group_venn_multi_comparison_grid.setContentsMargins(0, 0, 0, 0)
        self.group_venn_multi_comparison_grid.setHorizontalSpacing(12)
        self.group_venn_multi_comparison_grid.setVerticalSpacing(6)
        self.group_venn_multi_select_all = QPushButton("Select all Pairing comparisons")
        self.group_venn_multi_select_all.clicked.connect(lambda: self._set_group_venn_multi_all_checked(True))
        self.group_venn_multi_clear = QPushButton("Clear")
        self.group_venn_multi_clear.clicked.connect(lambda: self._set_group_venn_multi_all_checked(False))
        self.group_venn_multi_apply = QPushButton("Run Direction-aware Venn")
        self.group_venn_multi_apply.clicked.connect(self._apply_group_venn_multi_recalculate)
        self.group_venn_multi_controls = QWidget()
        group_venn_multi_layout = QVBoxLayout(self.group_venn_multi_controls)
        group_venn_multi_layout.setContentsMargins(0, 0, 0, 0)
        group_venn_multi_top = QHBoxLayout()
        group_venn_multi_top.addWidget(QLabel("Significant AS: FDR <= 0.05, |dPSI| >="))
        group_venn_multi_top.addWidget(self.group_venn_multi_sig_cutoff)
        group_venn_multi_top.addWidget(QLabel("min reads"))
        group_venn_multi_top.addWidget(self.group_venn_multi_min_reads)
        group_venn_multi_top.addWidget(self.group_venn_multi_select_all)
        group_venn_multi_top.addWidget(self.group_venn_multi_clear)
        group_venn_multi_top.addWidget(self.group_venn_multi_apply)
        group_venn_multi_top.addStretch(1)
        group_venn_multi_layout.addLayout(group_venn_multi_top)
        group_venn_multi_layout.addWidget(QLabel("Auto-loaded comparisons from Pairing / selected comparisons:"))
        group_venn_multi_layout.addWidget(self.group_venn_multi_comparison_panel)
        self.group_venn_multi_controls.setVisible(False)
        self._group_venn_multi_cache: dict[str, object] = {}

        for widget in (
            self.tabs,
            self.landscape_compare_tabs,
            self.landscape_tabs,
            self.group_tabs,
            self.group_mode_tabs,
            self.tx_tabs,
            self.mechanism_compare_tabs,
            self.mechanism_tabs,
            self.candidate_tabs,
            self.candidate_mode_tabs,
        ):
            widget.currentChanged.connect(self.refresh)
            widget.setDocumentMode(True)
            widget.setUsesScrollButtons(True)
            widget.setStyleSheet("QTabWidget::pane { border: 0; margin: 0; padding: 0; }")
        for selector in (
            self.landscape_compare_tabs,
            self.landscape_tabs,
            self.group_tabs,
            self.group_mode_tabs,
            self.tx_tabs,
            self.mechanism_compare_tabs,
            self.mechanism_tabs,
            self.candidate_tabs,
            self.candidate_mode_tabs,
        ):
            selector.setMaximumHeight(46)

        self.run_current_button = QPushButton("Run Current Branch")
        self.run_current_button.clicked.connect(self._run_current_branch)
        self.detail_button = QPushButton("Detail")
        self.detail_button.clicked.connect(self._show_detail_dialog)
        self.export_png_button = QPushButton("Export PNG")
        self.export_png_button.clicked.connect(lambda: self._export_current_figure("png"))
        self.export_pdf_button = QPushButton("Export PDF")
        self.export_pdf_button.clicked.connect(lambda: self._export_current_figure("pdf"))
        self.plot_scale = QSpinBox()
        self.plot_scale.setRange(30, 220)
        self.plot_scale.setValue(100)
        self.plot_scale.setMinimumWidth(76)
        self.plot_scale.setToolTip("Plot scale percentage")
        self.plot_scale.valueChanged.connect(self._redraw_scaled_figure)
        self.plot_width_scale = QSpinBox()
        self.plot_width_scale.setRange(50, 260)
        self.plot_width_scale.setValue(100)
        self.plot_width_scale.setMinimumWidth(76)
        self.plot_width_scale.setToolTip("Width scale percentage")
        self.plot_width_scale.valueChanged.connect(self._redraw_scaled_figure)
        self.plot_height_scale = QSpinBox()
        self.plot_height_scale.setRange(50, 260)
        self.plot_height_scale.setValue(100)
        self.plot_height_scale.setMinimumWidth(76)
        self.plot_height_scale.setToolTip("Height scale percentage")
        self.plot_height_scale.valueChanged.connect(self._redraw_scaled_figure)

        self.figure = Figure(figsize=(12, 8))
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.setMinimumHeight(480)
        self.canvas.setMinimumWidth(640)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        self.canvas_host = QWidget()
        self.canvas_host_layout = QVBoxLayout(self.canvas_host)
        self.canvas_host_layout.setContentsMargins(0, 0, 0, 0)
        self.canvas_host_layout.setSpacing(0)
        self.canvas_host_layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.canvas_host_layout.addWidget(self.canvas, 0, Qt.AlignmentFlag.AlignCenter)
        self.canvas_scroll = QScrollArea()
        self.canvas_scroll.setWidgetResizable(False)
        self.canvas_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.canvas_scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.canvas_scroll.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.canvas_scroll.setWidget(self.canvas_host)
        self.canvas_scroll.viewport().installEventFilter(self)
        self.preview_path = QLineEdit()
        self.preview_path.setReadOnly(True)
        self.preview_path.setPlaceholderText("Current preview file path")
        self.preview_button = QPushButton("Data")
        self.preview_button.clicked.connect(self._show_preview_dialog)
        self.open_folder_button = QPushButton("Open output\nfolder")
        self.open_folder_button.clicked.connect(self._open_preview_output_folder)
        self.download_button = QPushButton("Download\nfile")
        self.download_button.clicked.connect(self._download_preview_file)
        self.export_filtered_button = QPushButton("Export current\nfiltered table")
        self.export_filtered_button.clicked.connect(self._export_filtered_preview)
        self.preview_search = QLineEdit()
        self.preview_search.setPlaceholderText("Search gene_symbol / gene_id / comparison_id")
        self.preview_search.textChanged.connect(self._apply_preview_filters)
        self.preview_filter_column = QComboBox()
        self.preview_filter_column.currentTextChanged.connect(self._refresh_preview_filter_values)
        self.preview_filter_value = QComboBox()
        self.preview_filter_value.currentTextChanged.connect(self._apply_preview_filters)
        self.preview_info = QLabel("")
        self.preview_info.setWordWrap(True)
        self.preview_info.setStyleSheet("font-size: 10px; color: #bdbdbd;")
        self.preview_table = QTableWidget(0, 0)
        self.preview_table.setSortingEnabled(True)
        self.preview_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectItems)
        self.preview_table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.preview_table.setHorizontalScrollMode(QAbstractItemView.ScrollMode.ScrollPerPixel)
        self.preview_table.setVerticalScrollMode(QAbstractItemView.ScrollMode.ScrollPerPixel)
        self.preview_table.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.preview_table.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.preview_table.setSizeAdjustPolicy(QAbstractScrollArea.SizeAdjustPolicy.AdjustIgnored)
        self.preview_table.setMinimumWidth(0)
        QShortcut(QKeySequence.StandardKey.Copy, self.preview_table, activated=self._copy_preview_selection)
        self.footer_info = QTextEdit()
        self.footer_info.setReadOnly(True)
        self.footer_info.setStyleSheet("font-size: 10px;")
        self.footer_info.setMinimumHeight(40)
        self.footer_info.setMaximumHeight(56)
        self.details = QTextEdit()
        self.details.setReadOnly(True)
        self.details.setStyleSheet("font-size: 10px;")
        self.details.setMinimumHeight(42)
        self.details.setMaximumHeight(60)

        self._hover_annotation = None
        self._scatter_collection = None
        self._scatter_frame = pd.DataFrame()
        self._bar_patches: list[tuple[object, str]] = []
        self._pie_patches: list[tuple[object, str]] = []
        self._heatmap_image = None
        self._heatmap_rows: list[str] = []
        self._heatmap_cols: list[str] = []
        self._heatmap_values: np.ndarray | None = None
        self._preview_source_frame = pd.DataFrame()
        self._preview_filtered_frame = pd.DataFrame()
        self._preview_columns: list[str] = []
        self._preview_file_path: Path | None = None
        self._preview_output_dir: Path | None = None
        self._preview_name = ""
        self._preview_max_rows = 500
        self._detail_dialog_text = ""
        self._layout_options: dict[str, object] = {"rotate_x": False, "bottom": None, "left": 0.12, "right": 0.76, "top": 0.84}
        self._candidate_card_selected_gene_key: str | None = None
        self._candidate_card_payload: dict[str, object] = {}
        self._pan_active = False
        self._pan_dragged = False
        self._pan_axes = None
        self._pan_origin: tuple[float, float] | None = None
        self._pan_xlim: tuple[float, float] | None = None
        self._pan_ylim: tuple[float, float] | None = None
        self._relayout_pending = False
        self._relayout_generation = 0

        button_row = FlowLayout(spacing=8)
        button_row.addWidget(self.run_current_button)
        button_row.addWidget(self.detail_button)
        button_row.addWidget(QLabel("Scale %"))
        button_row.addWidget(self.plot_scale)
        button_row.addWidget(QLabel("Width %"))
        button_row.addWidget(self.plot_width_scale)
        button_row.addWidget(QLabel("Height %"))
        button_row.addWidget(self.plot_height_scale)
        button_row.addWidget(self.export_png_button)
        button_row.addWidget(self.export_pdf_button)
        for widget in (
            self.run_current_button,
            self.detail_button,
            self.export_png_button,
            self.export_pdf_button,
            self.preview_button,
            self.open_folder_button,
            self.download_button,
            self.export_filtered_button,
            self.plot_scale,
            self.plot_width_scale,
            self.plot_height_scale,
            self.card_open_sashimi_button,
        ):
            widget.setMinimumHeight(32)
        for button in (
            self.run_current_button,
            self.detail_button,
            self.export_png_button,
            self.export_pdf_button,
            self.preview_button,
            self.open_folder_button,
            self.download_button,
            self.export_filtered_button,
            self.card_open_sashimi_button,
        ):
            button.setMinimumWidth(128)

        preview_actions = FlowLayout(spacing=8)
        preview_actions.addWidget(self.preview_button)
        preview_actions.addWidget(self.open_folder_button)
        preview_actions.addWidget(self.download_button)
        preview_actions.addWidget(self.export_filtered_button)

        preview_filters = FlowLayout(spacing=8)
        preview_filters.addWidget(QLabel("File"))
        preview_filters.addWidget(self.preview_path)
        preview_filters.addWidget(QLabel("Filter column"))
        preview_filters.addWidget(self.preview_filter_column)
        preview_filters.addWidget(QLabel("Value"))
        preview_filters.addWidget(self.preview_filter_value)
        preview_filters.addWidget(self.preview_search)

        content_splitter = QSplitter(Qt.Orientation.Vertical)
        content_splitter.addWidget(self.canvas_scroll)
        content_splitter.addWidget(self.details)
        content_splitter.setStretchFactor(0, 9)
        content_splitter.setStretchFactor(1, 1)
        content_splitter.setHandleWidth(1)
        content_splitter.setChildrenCollapsible(False)
        content_splitter.setStyleSheet("QSplitter::handle { background-color: #202124; border: 0px; }")
        self.details.setMinimumHeight(40)
        self.details.setMaximumHeight(52)
        content_splitter.setSizes([1800, 44])

        layout = QVBoxLayout(self)
        layout.addWidget(self.label)
        layout.addLayout(button_row)
        layout.addWidget(self.tabs)
        layout.addWidget(self.landscape_stable_controls)
        layout.addWidget(self.candidate_controls)
        layout.addWidget(self.candidate_card_controls)
        layout.addWidget(self.group_pattern_controls)
        layout.addWidget(self.group_interaction_controls)
        layout.addWidget(self.group_venn_controls)
        layout.addWidget(self.group_venn_multi_controls)
        layout.addWidget(self.toolbar)
        layout.addWidget(content_splitter, 1)
        layout.addLayout(preview_actions)
        layout.addWidget(self.footer_info)
        self.footer_info.setVisible(False)
        self.context_header.setVisible(False)
        self.explanation.setVisible(False)

        self.canvas.mpl_connect("motion_notify_event", self._on_motion)
        self.canvas.mpl_connect("button_press_event", self._on_press)
        self.canvas.mpl_connect("button_release_event", self._on_release)
        self.canvas.mpl_connect("scroll_event", self._on_scroll)

        self._group_pairs: list[tuple[str, str]] = []
        self._landscape_comparisons: list[tuple[str, str]] = []
        self._tx_comparisons: list[tuple[str, str]] = []
        self._candidate_comparisons: list[tuple[str, str]] = []
        self._navigation_branch: str | None = None
        self._navigation_subbranch: str | None = None

    def resizeEvent(self, event) -> None:
        super().resizeEvent(event)
        if not self.isVisible():
            return
        self._schedule_post_render_layout()

    def showEvent(self, event) -> None:
        super().showEvent(event)
        self._schedule_post_render_layout()

    def eventFilter(self, watched, event) -> bool:
        if watched is self.canvas_scroll.viewport() and event.type() in {
            QEvent.Type.Resize,
            QEvent.Type.Show,
            QEvent.Type.LayoutRequest,
        }:
            self._schedule_post_render_layout()
        return super().eventFilter(watched, event)

    def navigate_to(self, branch: str, subbranch: str | None = None) -> None:
        if self.fixed_branch is not None:
            self._navigation_branch = self.fixed_branch
            self._navigation_subbranch = self.fixed_subbranch
            self.refresh()
            return
        branch_map = {
            "landscape": 0,
            "group": 1,
            "tx": 2,
            "mechanism": 3,
            "candidate": 4,
        }
        if branch not in branch_map:
            return
        self.tabs.setCurrentIndex(branch_map[branch])
        if branch == "landscape":
            sub_map = {"overview": 0, "bar": 1, "pie": 2}
            if subbranch in sub_map:
                self.landscape_tabs.setCurrentIndex(sub_map[subbranch])
        elif branch == "group":
            sub_map = {"excel": 0, "heatmap": 1, "direction": 2}
            if subbranch in sub_map:
                self.group_mode_tabs.setCurrentIndex(sub_map[subbranch])
        elif branch == "candidate":
            sub_map = {"table": 0, "heatmap": 1, "program": 2, "cards": 3, "isoform": 4}
            if subbranch in sub_map:
                self.candidate_mode_tabs.setCurrentIndex(sub_map[subbranch])
        self._navigation_branch = branch
        self._navigation_subbranch = subbranch
        self.refresh()

    def clear_navigation_mode(self) -> None:
        if self.fixed_branch is not None:
            self._navigation_branch = self.fixed_branch
            self._navigation_subbranch = self.fixed_subbranch
            return
        self._navigation_branch = None
        self._navigation_subbranch = None

    def refresh(self) -> None:
        self._detail_dialog_text = ""
        if self.fixed_branch is not None:
            self._navigation_branch = self.fixed_branch
            self._navigation_subbranch = self.fixed_subbranch
        project = self.project_service.current_project
        self._rebuild_dynamic_tabs()
        self._sync_landscape_stable_controls_from_thresholds()
        if project is None:
            self.run_current_button.setEnabled(False)
            self._apply_navigation_visibility()
            self._set_explanation("Load a project first. After analysis is available, each branch will explain what it shows and how to interpret it.")
            self._draw_placeholder("Load a project first.")
            self.context_header.clear()
            self._sync_footer_info()
            self._reset_plot_view()
            return

        self._update_run_button_text()
        self._apply_navigation_visibility()

        self._autosize_figure()
        self._reset_layout_options()
        branch_key, subbranch_key = self._effective_branch_keys()
        if branch_key == "landscape":
            self._render_landscape_branch()
        elif branch_key == "group":
            self._render_group_comparison_branch()
        elif branch_key == "tx":
            self._render_tx_splicing_branch()
        elif branch_key == "mechanism":
            self._render_mechanism_branch()
        else:
            self._render_candidate_branch()
        self._sync_footer_info()
        self._schedule_post_render_layout()

    def _reset_layout_options(self) -> None:
        self._layout_options = {
            "rotate_x": False,
            "bottom": 0.18,
            "left": 0.12,
            "right": 0.78,
            "top": 0.82,
        }

    def _rebuild_dynamic_tabs(self) -> None:
        project = self.project_service.current_project
        if project is None:
            self._group_pairs = []
            self._landscape_comparisons = []
            self._tx_comparisons = []
            self._mechanism_comparisons = []
            self._candidate_comparisons = []
            for widget in (self.landscape_compare_tabs, self.group_tabs, self.tx_tabs, self.mechanism_compare_tabs, self.candidate_tabs):
                widget.clear()
            return

        enabled_pairs = [(pair.pair_id, pair.resolved_name) for pair in project.comparison_pairs if pair.enabled]
        selected_comparisons = [(comp.comparison_id, comp.resolved_name) for comp in self.project_service.selected_comparisons_for_display()]
        self._group_pairs = enabled_pairs
        self._landscape_comparisons = selected_comparisons
        self._candidate_comparisons = selected_comparisons
        self._tx_comparisons = selected_comparisons
        self._mechanism_comparisons = selected_comparisons

        self._sync_named_tabs(self.landscape_compare_tabs, self._landscape_comparisons, fallback_label="No selected comparisons")
        self._sync_named_tabs(self.group_tabs, self._group_pairs, fallback_label="No enabled comparison sets")
        self._sync_named_tabs(self.tx_tabs, self._tx_comparisons, fallback_label="No selected comparisons")
        self._sync_named_tabs(self.mechanism_compare_tabs, self._mechanism_comparisons, fallback_label="No selected comparisons")
        self._sync_named_tabs(self.candidate_tabs, self._candidate_comparisons, fallback_label="No selected comparisons")
        self._sync_group_venn_multi_comparison_checks()

    def _apply_navigation_visibility(self) -> None:
        branch = self.fixed_branch if self.fixed_branch is not None else self._navigation_branch
        subbranch = self.fixed_subbranch if self.fixed_branch is not None else self._navigation_subbranch

        standalone = self.fixed_branch is not None
        branch_index_map = {
            "landscape": 0,
            "group": 1,
            "tx": 2,
            "mechanism": 3,
            "candidate": 4,
        }
        if standalone:
            self.tabs.setVisible(True)
            self.tabs.tabBar().hide()
            target_index = branch_index_map.get(branch, 0)
            if self.tabs.currentIndex() != target_index:
                self.tabs.blockSignals(True)
                self.tabs.setCurrentIndex(target_index)
                self.tabs.blockSignals(False)
        else:
            self.tabs.tabBar().show()
            self.tabs.setVisible(branch is None)
        landscape_active = (not standalone and branch in (None, "landscape")) or (standalone and branch == "landscape")
        landscape_subbranch = self.fixed_subbranch if standalone and branch == "landscape" else self._current_landscape_subbranch()
        self.landscape_compare_tabs.setVisible(bool(landscape_active and landscape_subbranch == "pie"))
        self.landscape_tabs.setVisible((not standalone and branch in (None, "landscape") and subbranch is None) or (standalone and branch == "landscape" and subbranch is None))
        group_active = (not standalone and branch in (None, "group")) or (standalone and branch == "group")
        group_subbranch_for_visibility = subbranch or self._current_group_subbranch()
        group_needs_pair_tabs = CrossComparisonController.needs_pair_tabs(
            group_subbranch_for_visibility,
            interaction_view=str(self.group_interaction_view.currentData() or ""),
        )
        self.group_tabs.setVisible(bool(group_active and (group_needs_pair_tabs or (not standalone and subbranch is None))))
        self.group_mode_tabs.setVisible((not standalone and branch in (None, "group") and subbranch is None) or (standalone and branch == "group" and subbranch is None))
        self.tx_tabs.setVisible((not standalone and branch in (None, "tx")) or (standalone and branch == "tx"))
        self.mechanism_compare_tabs.setVisible((not standalone and branch in (None, "mechanism")) or (standalone and branch == "mechanism"))
        self.mechanism_tabs.setVisible((not standalone and branch in (None, "mechanism")) or (standalone and branch == "mechanism"))
        self.candidate_tabs.setVisible((not standalone and branch in (None, "candidate")) or (standalone and branch == "candidate"))
        self.candidate_mode_tabs.setVisible((not standalone and branch in (None, "candidate") and subbranch is None) or (standalone and branch == "candidate" and subbranch is None))
        self.landscape_stable_controls.setVisible(branch == "landscape" and (subbranch in {"overview", "bar", "pie"} or subbranch is None))
        self.candidate_controls.setVisible(branch == "candidate")
        self.candidate_card_controls.setVisible(branch == "candidate" and (subbranch == "cards" or (subbranch is None and self._current_candidate_subbranch() == "cards")))
        self.group_pattern_controls.setVisible(branch == "group" and subbranch == "pattern")
        self.group_interaction_controls.setVisible(branch == "group" and subbranch == "interaction")
        self.group_venn_controls.setVisible(branch == "group" and subbranch == "venn")
        self.group_venn_multi_controls.setVisible(branch == "group" and subbranch in {"venn_multi", "venn_multi_similarity"})
        if branch == "group" and subbranch in {"venn_multi", "venn_multi_similarity"}:
            self.group_venn_multi_apply.setText("Run Direction-aware Similarity" if subbranch == "venn_multi_similarity" else "Run Direction-aware Venn")
            self._sync_group_venn_multi_comparison_checks()

        if branch is None and not standalone:
            self.label.setText("Analysis workbench")
        else:
            label_map = {
                "landscape": "Analysis - Landscape",
                "group": "Analysis - Group Comparison",
                "tx": "Analysis - Transcript + Splicing",
                "mechanism": "Analysis - Mechanism Support",
                "candidate": "Analysis - Candidate Gene Screening",
            }
            self.label.setText(label_map.get(branch, "Analysis"))

    def _sync_named_tabs(self, widget: QTabWidget, items: list[tuple[str, str]], fallback_label: str) -> None:
        current_key = widget.tabToolTip(widget.currentIndex()) if widget.count() else None
        widget.blockSignals(True)
        widget.clear()
        if not items:
            placeholder = QWidget()
            widget.addTab(placeholder, fallback_label)
            widget.setTabToolTip(0, "")
        else:
            target_index = 0
            for index, (key, label) in enumerate(items):
                page = QWidget()
                widget.addTab(page, label)
                widget.setTabToolTip(index, key)
                if current_key and current_key == key:
                    target_index = index
            widget.setCurrentIndex(target_index)
        widget.blockSignals(False)

    def _update_run_button_text(self) -> None:
        branch_key, subbranch_key = self._effective_branch_keys()
        if branch_key == "landscape":
            labels = {
                "overview": "Run Landscape Overview",
                "bar": "Run Landscape Bar",
                "pie": "Run Landscape Pie",
            }
            active_subbranch = subbranch_key or self._current_landscape_subbranch()
            self.run_current_button.setText(labels.get(active_subbranch, "Run Landscape"))
        elif branch_key == "group":
            active_subbranch = subbranch_key or self._current_group_subbranch()
            self.run_current_button.setText(CrossComparisonController.run_label(active_subbranch))
        elif branch_key == "tx":
            labels = {
                "se": "Run SE Subanalysis",
                "ri": "Run RI Subanalysis",
                "a3ss": "Run A3SS Subanalysis",
                "a5ss": "Run A5SS Subanalysis",
                "mxe": "Run MXE Subanalysis",
            }
            self.run_current_button.setText(labels.get(subbranch_key or "", "Run Transcript + Splicing"))
        elif branch_key == "mechanism":
            self.run_current_button.setText("Run Mechanism Support")
        else:
            labels = {
                "table": "Refresh Candidate Table",
                "heatmap": "Refresh Evidence Heatmap",
                "program": "Refresh Program Heatmap",
                "cards": "Refresh Candidate Cards + Event Follow-up",
                "isoform": "Refresh Isoform Follow-up",
            }
            self.run_current_button.setText(labels.get(subbranch_key or "", "Refresh Candidate Screening"))
        self.run_current_button.setEnabled(self.project_service.current_project is not None)

    def _run_current_branch(self) -> None:
        if self.project_service.current_project is None:
            return
        branch_key, subbranch_key = self._effective_branch_keys()
        if branch_key == "landscape":
            landscape_subbranch = subbranch_key or self._current_landscape_subbranch()
            if landscape_subbranch in ("overview", "bar", "pie"):
                self._apply_landscape_stable_threshold_controls()
                self.project_service.set_selected_analysis_modules(["splicing_landscape"])
                if callable(self.on_run_modules):
                    self.on_run_modules()
        elif branch_key == "group":
            group_subbranch = subbranch_key or self._current_group_subbranch()
            action = CrossComparisonController.action(group_subbranch)
            if action == "candidate_matrix":
                self.project_service.preview_cross_comparison_candidate_matrix(allow_generate=True)
                self.refresh()
            elif action == "venn":
                self._apply_group_venn_recalculate()
            elif action == "venn_multi":
                self._apply_group_venn_multi_recalculate()
            elif action == "pattern":
                self._apply_group_pattern_recalculate()
            elif action == "interaction":
                self._apply_group_interaction_recalculate()
            else:
                self.project_service.set_selected_analysis_modules(["program_comparison"])
                if callable(self.on_run_modules):
                    self.on_run_modules()
        elif branch_key == "tx":
            self.project_service.set_selected_analysis_modules(["tx_splicing_integration"])
            if callable(self.on_run_modules):
                self.on_run_modules()
        elif branch_key == "mechanism":
            self.project_service.set_selected_analysis_modules(["mechanism_support"])
            if callable(self.on_run_modules):
                self.on_run_modules()
        else:
            self.refresh()

    def _sync_landscape_stable_controls_from_thresholds(self) -> None:
        thresholds = self.project_service.current_thresholds
        self.landscape_stable_min_reads.blockSignals(True)
        self.landscape_stable_min_reps.blockSignals(True)
        self.landscape_stable_consistency.blockSignals(True)
        self.landscape_stable_min_reads.setValue(int(getattr(thresholds, "stable_min_junction_reads", 20)))
        self.landscape_stable_min_reps.setValue(int(getattr(thresholds, "stable_min_replicates_per_group", 2)))
        self.landscape_stable_consistency.setValue(float(getattr(thresholds, "stable_min_pairwise_consistency", 0.67)))
        self.landscape_stable_min_reads.blockSignals(False)
        self.landscape_stable_min_reps.blockSignals(False)
        self.landscape_stable_consistency.blockSignals(False)

    def _apply_landscape_stable_threshold_controls(self) -> None:
        self.project_service.update_thresholds(
            stable_min_junction_reads=float(self.landscape_stable_min_reads.value()),
            stable_min_replicates_per_group=int(self.landscape_stable_min_reps.value()),
            stable_min_pairwise_consistency=float(self.landscape_stable_consistency.value()),
        )

    def _render_landscape_branch(self) -> None:
        landscape = self.project_service.run_state.results.splicing_landscape
        selected_ids = [comparison_id for comparison_id, _label in self._landscape_comparisons]
        filtered = landscape.copy() if landscape is not None else pd.DataFrame()
        if selected_ids and filtered is not None and not filtered.empty and "comparison_id" in filtered.columns:
            filtered = filtered.loc[filtered["comparison_id"].astype(str).isin(selected_ids)].copy()
        self.context_header.clear()
        current = self.fixed_subbranch or self._navigation_subbranch or self._current_landscape_subbranch()
        if current == "overview":
            self.context_header.clear()
            self._draw_landscape_overview(filtered)
        elif current == "bar":
            self.context_header.clear()
            self._draw_landscape_bar(filtered)
        elif current == "pie":
            comparison_id = self._current_landscape_comparison_id()
            self.context_header.setText(self.project_service.comparison_context_text(comparison_id))
            pie_frame = self._filter_frame_by_value(filtered, "comparison_id", comparison_id)
            self._draw_landscape_pie(pie_frame, comparison_id, self._comparison_name_from_id(comparison_id))
        else:
            self.context_header.clear()
            self._draw_placeholder("Landscape subpage not available.")

    def _render_group_comparison_branch(self) -> None:
        pair_id = self._current_tab_key(self.group_tabs)
        pair_name = self._pair_name_from_id(pair_id)
        mode = self.fixed_subbranch or self._navigation_subbranch or self._current_group_subbranch()
        cross_modes = {"matrix", "shared", "reversal"}
        if mode == "interaction":
            abs_cutoff = float(self.group_interaction_sig_cutoff.value())
            delta_cutoff = float(self.group_interaction_delta_cutoff.value())
            min_reads = int(self.group_interaction_min_reads.value())
            selected_key = str(self.group_interaction_view.currentData() or "unified_event_table")
            source_table_key = selected_key
            is_batch_view = selected_key == "batch_summary"
            if is_batch_view:
                batch_tables = self.project_service.preview_cross_comparison_as_interaction_batch(
                    abs_dpsi_cutoff=abs_cutoff,
                    large_delta_dpsi_cutoff=delta_cutoff,
                    min_total_reads=min_reads,
                    allow_generate=False,
                )
                state = self.project_service.module_state("cross_comparison_as_interactions")
                self.context_header.setText(
                    "5.4.3 all-pair batch summary | Runs the AS interaction/read-support model across every selected comparison pair."
                )
                self._set_explanation(
                    "All-pair 5.4.3 batch summary applies the same replicate-aware AS interaction/read-support model to every configured or auto-derived comparison pair. "
                    "It summarizes shared significant AS interactions, opposite-direction events, same-direction large-delta events, and A/B-specific read-support classes."
                )
                if state.get("status") == "not_run":
                    self._draw_placeholder("Not run.\nNo cached all-pair 5.4.3 batch summary is loaded.")
                    self.details.setPlainText(state.get("message", "All-pair 5.4.3 batch summary has not been run yet."))
                    self._sync_footer_info()
                    self._reset_plot_view()
                    return
                selected_table = batch_tables.get("batch_summary", pd.DataFrame())
                self._draw_cross_as_interaction_batch_summary(batch_tables, abs_cutoff, delta_cutoff, min_reads)
                batch_files = self._cross_interaction_batch_file_map()
                preferred_columns = [
                    "batch_pair_order",
                    "pair_id",
                    "pair_name",
                    "comparison_A",
                    "comparison_B",
                    "n_shared_significant",
                    "n_high_confidence_interaction",
                    "n_moderate_interaction",
                    "n_opposite_direction",
                    "n_same_direction_large_delta",
                    "n_exclusive",
                    "n_exclusive_A",
                    "n_exclusive_B",
                    "n_exclusive_supported",
                    "n_exclusive_A_supported",
                    "n_exclusive_B_supported",
                    "n_exclusive_low_power",
                    "n_exclusive_counterpart_not_found",
                    "delta_top20_model_top20_overlap",
                    "delta_vs_model_spearman_rank_correlation",
                    "gene",
                    "event_type",
                    "event_id",
                    "event_matching_key",
                    "dPSI_A",
                    "dPSI_B",
                    "FDR_A",
                    "FDR_B",
                    "abs_delta_dPSI",
                    "interaction_priority_score",
                    "interaction_class",
                    "exclusive_side",
                    "counterpart_event_found",
                    "non_significant_side_mean_total_reads",
                    "non_significant_side_read_support",
                    "top_gene",
                    "top_event_type",
                    "top_abs_delta_dPSI",
                    "top_interaction_priority_score",
                    "simple_delta_score",
                    "simple_delta_priority_score",
                    "delta_dPSI_rank",
                    "interaction_model_rank",
                    "delta_vs_model_rank_shift",
                    "benchmark_class",
                    "benchmark_score_delta",
                    "quasi_binomial_dispersion",
                    "beta_binomial_rho_estimate",
                    "beta_binomial_interaction_log_odds",
                    "beta_binomial_likelihood_ratio",
                    "beta_binomial_pvalue",
                    "beta_binomial_FDR",
                    "beta_binomial_rho",
                    "beta_binomial_model_status",
                    "model_diagnostic_class",
                    "model_warning",
                    "abs_dPSI_significance_cutoff_used",
                    "large_delta_dPSI_cutoff_used",
                    "min_total_reads_used",
                ]
                self._populate_table(
                    selected_table,
                    preferred_columns,
                    file_path=batch_files.get("batch_summary"),
                    output_dir=self.project_service.cross_comparison_interaction_batch_output_dir(),
                    preview_name="5.4.3 all-pair summary table",
                )
                self._set_detail_text(
                    module_detail("5.4.interaction")
                    + "\n\n5.4.3 all-pair batch summary.\n"
                    "Scope: every user-defined comparison pair plus any missing pairwise combinations among selected comparisons.\n"
                    f"Shared-event filter: rMATS_FDR <= 0.05 and abs(dPSI) >= {abs_cutoff:.3f} in both comparisons.\n"
                    f"Large delta cutoff: {delta_cutoff:.3f}; min mean total junction reads: {min_reads}.\n"
                    "Batch summary = one row per pair with counts and the top-ranked interaction event.\n"
                    "All-pair top interactions = top 100 interaction-priority events per pair.\n"
                    "All-pair exclusive support is split explicitly into A-specific and B-specific counts.\n"
                    "- n_exclusive_A: significant only in comparison A, then checked for counterpart read support in comparison B.\n"
                    "- n_exclusive_B: significant only in comparison B, then checked for counterpart read support in comparison A.\n"
                    "- n_exclusive_A_supported / n_exclusive_B_supported: subset where the counterpart event exists and has enough read support but is not significant.\n"
                    "Delta benchmark compares simple abs(delta dPSI) ranking against the replicate-aware interaction model ranking, so you can see which events are only obvious by effect size versus model-supported."
                    " The model-supported view includes quasi-binomial diagnostics and a formal beta-binomial likelihood-ratio interaction test where replicate IJC/SJC counts are available.\n"
                    "A/B-specific read-support definitions:\n"
                    "- supported_non_significant: the same event is found in the other comparison and has enough mean IJC+SJC read support, but does not pass FDR/dPSI significance there.\n"
                    "- low_power_or_absent: the counterpart side has low mean junction read support, so absence of significance should be interpreted cautiously.\n"
                    "- counterpart_event_not_found: the same event key is not found in the other comparison's rMATS event list.\n"
                    "read support = mean total junction reads from replicate IJC+SJC values on the non-significant/counterpart side.",
                    "5.4.3 batch result loaded. Open Detail for A/B-specific support, model and benchmark definitions.",
                )
                return
            interaction_tables = self.project_service.preview_cross_comparison_as_interactions(
                pair_id,
                abs_dpsi_cutoff=abs_cutoff,
                large_delta_dpsi_cutoff=delta_cutoff,
                min_total_reads=min_reads,
                allow_generate=False,
            )
            state = self.project_service.module_state("cross_comparison_as_interactions")
            self.context_header.setText(
                f"Comparison set: {pair_name or 'selected pair'} | 5.4.3 models AS-event differences with replicate-level rMATS counts and read support."
            )
            self._set_explanation(
                "5.4.3 is the main cross-comparison AS event module. It uses rMATS significant events and replicate-level IJC/SJC counts in a binomial GLM: inclusion ~ comparison + group + comparison:group. "
                "The interaction term tests whether the AS effect differs between comparison A and B; Pearson overdispersion is used as a quasi-binomial-style correction. "
                "A-specific and B-specific tables are split explicitly, so events significant only in A are not mixed with events significant only in B."
            )
            if state.get("status") == "not_run":
                self._draw_placeholder("Not run.\nNo cached 5.4.3 interaction/read-support result is loaded.")
                self.details.setPlainText(state.get("message", "5.4.3 has not been run yet."))
                self._sync_footer_info()
                self._reset_plot_view()
                return
            selected_table = interaction_tables.get(source_table_key, pd.DataFrame())
            self._draw_cross_as_interaction_summary(interaction_tables, pair_name, abs_cutoff, delta_cutoff, min_reads)
            interaction_files = self._cross_interaction_file_map(pair_id)
            preferred_columns = [
                "gene",
                "event_type",
                "event_id",
                "event_matching_key",
                "unified_display_class",
                "unified_result_type",
                "support_class",
                "event_class",
                "comparison_A",
                "comparison_B",
                "experiment_A_display_name",
                "experiment_B_display_name",
                "dPSI_A",
                "dPSI_B",
                "FDR_A",
                "FDR_B",
                "delta_dPSI",
                "abs_delta_dPSI",
                "same_direction",
                "opposite_direction",
                "same_direction_strength_delta",
                "interaction_log_odds",
                "interaction_se",
                "interaction_z",
                "interaction_pvalue",
                "interaction_FDR",
                "interaction_model",
                "interaction_model_status",
                "interaction_overdispersion_scale",
                "interaction_model_n_observations",
                "interaction_model_df_resid",
                "interaction_model_pearson_chi2",
                "interaction_model_deviance",
                "interaction_model_total_reads",
                "quasi_binomial_dispersion",
                "quasi_binomial_overdispersed",
                "beta_binomial_rho_estimate",
                "beta_binomial_interaction_log_odds",
                "beta_binomial_likelihood_ratio",
                "beta_binomial_pvalue",
                "beta_binomial_FDR",
                "beta_binomial_rho",
                "beta_binomial_full_loglik",
                "beta_binomial_null_loglik",
                "beta_binomial_model_status",
                "beta_binomial_warning",
                "model_valid",
                "model_diagnostic_class",
                "model_warning",
                "eb_prior_variance",
                "shrinkage_weight",
                "shrunken_interaction_log_odds",
                "shrunken_interaction_z",
                "shrunken_interaction_pvalue",
                "shrunken_interaction_FDR",
                "interaction_priority_score",
                "simple_delta_score",
                "simple_delta_priority_score",
                "delta_dPSI_rank",
                "interaction_model_rank",
                "delta_vs_model_rank_shift",
                "simple_delta_top20",
                "model_top20",
                "benchmark_class",
                "benchmark_score_delta",
                "interaction_class",
                "read_support_class",
                "replicate_consistency_score",
                "mean_total_reads_A",
                "mean_total_reads_B",
                "replicate_count_A",
                "replicate_count_B",
                "exclusive_side",
                "counterpart_event_found",
                "significant_side_mean_total_reads",
                "non_significant_side_mean_total_reads",
                "significant_side_replicate_count",
                "non_significant_side_replicate_count",
                "non_significant_side_dPSI",
                "non_significant_side_FDR",
                "non_significant_side_read_support",
                "abs_dPSI_significance_cutoff_used",
                "large_delta_dPSI_cutoff_used",
                "min_total_reads_used",
                "pair_id",
                "pair_name",
                "n_events_A",
                "n_events_B",
                "n_shared_significant",
                "n_opposite_direction",
                "n_large_delta",
                "n_exclusive_supported",
                "n_exclusive_low_power",
                "n_exclusive_counterpart_not_found",
                "n_exclusive_not_shared_significant",
            ]
            self._populate_table(
                selected_table,
                preferred_columns,
                file_path=interaction_files.get(source_table_key),
                output_dir=self.project_service.cross_comparison_interaction_output_dir(pair_id),
                preview_name="5.4.3 unified AS interaction/read-support event table",
            )
            self._set_detail_text(
                module_detail("5.4.interaction")
                + f"\n\nCurrent selected pair: {pair_name or 'selected pair'}.\n"
                f"Shared-event filter: rMATS_FDR <= 0.05 and abs(dPSI) >= {abs_cutoff:.3f} in both comparisons.\n"
                f"Large delta cutoff: {delta_cutoff:.3f}; min mean total junction reads: {min_reads}.\n"
                "Event matching key = normalized gene + event_type + event_id, so multiple events from the same gene stay separate.\n"
                "Model: binomial GLM on replicate-level IJC/SJC counts, inclusion ~ comparison + group + comparison:group.\n"
                "interaction_log_odds is reported as comparison A group effect minus comparison B group effect, so its sign follows delta_dPSI.\n"
                "Diagnostics include Pearson quasi-binomial dispersion, overdispersion flag and model validity. "
                "Beta-binomial output is a separate likelihood-ratio interaction fit with beta_binomial_pvalue/FDR and beta_binomial_rho; beta_binomial_rho_estimate is kept only as the quick overdispersion-derived diagnostic estimate.\n"
                "Empirical-Bayes shrinkage pulls unstable log-odds toward zero based on sampling SE, and interaction_priority_score ranks events by delta dPSI, shrunken statistical confidence, replicate consistency and read support.\n"
                "Delta benchmark compares simple abs(delta dPSI) ranking with model ranking; benchmark_class marks events that are top-ranked by delta only, model only, or both.\n"
                "A-specific/B-specific tables are not interaction tests by themselves; they flag events significant in only one side of the selected A-vs-B comparison pair. "
                "For each exclusive event, Katana checks the counterpart event in the other comparison and classifies it as supported_non_significant, low_power_or_absent, or counterpart_event_not_found using mean IJC+SJC read support.\n"
                "A/B-specific support definitions:\n"
                "- comparison_A_only: significant in comparison A, then checked for counterpart read support in comparison B.\n"
                "- comparison_B_only: significant in comparison B, then checked for counterpart read support in comparison A.\n"
                "- supported_non_significant: counterpart event exists and has enough read support but fails FDR/dPSI significance.\n"
                "- low_power_or_absent: counterpart event has low read support, so this is not strong absence evidence.\n"
                "- counterpart_event_not_found: same event key was not detected in the other comparison table.\n"
                "replicate_consistency_score = how consistently replicate-level PSI shifts agree with the comparison-level direction; higher is more reproducible.\n"
                "\nLegacy note:\n"
                "- The old summary-only dPSI pattern page has been removed from the main navigation because this read-count-aware module is the stronger default path.\n"
                "- Setting min reads very low can emulate weak read-support filtering, but the model still uses replicate IJC/SJC structure when available.",
                "5.4.3 result loaded. Open Detail for model and A/B-specific read-support definitions.",
            )
            return
        if mode == "pattern":
            abs_cutoff = float(self.group_pattern_sig_cutoff.value())
            delta_cutoff = float(self.group_pattern_large_delta.value())
            pattern_tables = self.project_service.preview_cross_comparison_significant_as_patterns(
                pair_id,
                abs_dpsi_cutoff=abs_cutoff,
                large_delta_dpsi_cutoff=delta_cutoff,
                allow_generate=False,
            )
            state = self.project_service.module_state("cross_comparison_as_patterns")
            self.context_header.setText(
                f"Legacy summary-only dPSI pattern: {pair_name or 'selected pair'} | Hidden from main workflow; use 5.4.3 for read-count-aware AS interaction."
            )
            self._set_explanation(
                "Legacy dPSI-only pattern analysis is kept only for backward cache compatibility. It does not use replicate read counts, so it has been removed from the main 5.4 navigation."
            )
            if state.get("status") == "not_run":
                self._draw_placeholder("Legacy module not run.\nUse 5.4.3 Cross-comparison AS interaction modeling instead.")
                self.details.setPlainText(state.get("message", "Legacy summary-only dPSI pattern has not been run."))
                self._sync_footer_info()
                self._reset_plot_view()
                return
            selected_key = str(self.group_pattern_view.currentData() or "same_direction_large_delta")
            selected_table = pattern_tables.get(selected_key, pd.DataFrame())
            self._draw_cross_as_pattern_summary(pattern_tables, pair_name, abs_cutoff, delta_cutoff)
            cross_files = self._cross_pattern_file_map(pair_id)
            preferred_columns = [
                "gene",
                "event_type",
                "event_id",
                "event_matching_key",
                "comparison_A",
                "comparison_B",
                "experiment_A_display_name",
                "experiment_B_display_name",
                "dPSI_A",
                "dPSI_B",
                "FDR_A",
                "FDR_B",
                "significant_A",
                "significant_B",
                "direction_A",
                "direction_B",
                "delta_dPSI",
                "abs_delta_dPSI",
                "same_direction",
                "opposite_direction",
                "same_direction_strength_delta",
                "event_class",
                "abs_dPSI_significance_cutoff_used",
                "large_delta_dPSI_cutoff_used",
            ]
            self._populate_table(
                selected_table,
                preferred_columns,
                file_path=cross_files.get(selected_key),
                output_dir=self.project_service.cross_comparison_pattern_output_dir(pair_id),
                preview_name=f"5.4.3 {self.group_pattern_view.currentText()}",
            )
            self._set_detail_text(
                module_detail("5.4.dpsi")
                + f"\n\nCurrent selected pair: {pair_name or 'selected pair'}.\n"
                f"abs(dPSI) significance cutoff: {abs_cutoff:.3f}\n"
                f"large cross-comparison delta cutoff: {delta_cutoff:.3f}\n"
                "Event matching key = normalized gene + event_type + event_id (coordinates fallback when event_id is missing).\n"
                "If no rows pass the current filters, the table will show 'No events found'.\n"
                "\nLegacy status:\n"
                "- This summary-only dPSI page was removed from main navigation because it is mostly a weaker version of the read-count-aware 5.4.3 module.\n"
                "- Use 5.4.3 for the main AS-event comparison because it keeps A/B-specific read support and replicate-level IJC/SJC modeling.",
                "Legacy dPSI pattern result loaded. Prefer 5.4.3 for main interpretation.",
            )
            return
        if mode in {"venn", "venn_multi", "venn_multi_similarity"}:
            cross_matrix = pd.DataFrame()
            filtered_summary = pd.DataFrame()
            filtered_events = pd.DataFrame()
        elif mode in cross_modes:
            cross_matrix_all = self.project_service.preview_cross_comparison_candidate_matrix(allow_generate=False)
            cross_matrix = self._pairwise_candidate_matrix(cross_matrix_all, pair_id)
            self.context_header.setText(
                f"Comparison pair: {pair_name or 'selected pair'} | Gene-level candidate views compare comparison A vs comparison B only."
            )
            if cross_matrix_all is None or cross_matrix_all.empty:
                state = self.project_service.module_state("cross_comparison_candidates")
                self._set_explanation(
                    "These cross-comparison candidate pages do not run program comparison. They only read cached 5.3.1 candidate results and build pairwise candidate-pattern summaries for the selected comparison pair."
                )
                self._draw_placeholder("Not run.\nNo cached cross-comparison candidate matrix is loaded.")
                self.details.setPlainText(state.get("message", "Cross-comparison candidate matrix not available yet."))
                self._sync_footer_info()
                self._reset_plot_view()
                return
            if cross_matrix.empty:
                self._set_explanation(
                    "No candidate genes are present in either comparison of the selected pair. Switch pair tabs or run 5.3.1 for those comparisons first."
                )
                self._draw_placeholder(f"No pairwise candidate rows found for {pair_name or 'selected pair'}.")
                self.details.setPlainText(
                    "5.4.1/5.4.2 are now pairwise views.\n"
                    "They need a selected comparison pair with comparison A and comparison B, plus cached 5.3.1 candidate results for at least one side.\n"
                    "They no longer show the all-comparison global matrix as the main page."
                )
                self._sync_footer_info()
                self._reset_plot_view()
                return
            filtered_summary = pd.DataFrame()
            filtered_events = pd.DataFrame()
        else:
            if self.project_service.run_state.results.program_summary.empty and self.project_service.run_state.results.program_events.empty:
                self.context_header.clear()
                self._set_explanation("Group comparison has not been run yet. Open this page, then click 'Run Current Branch' to compute the current comparison set only.")
                self._draw_placeholder("Not run.\nNo cached group-comparison result is loaded.")
                self._sync_footer_info()
                self._reset_plot_view()
                return
            program_events, program_summary = self.project_service.preview_program_comparison(allow_generate=False)
            self.context_header.setText(
                f"Comparison set: {pair_name or 'selected pair'} | This view compares already-screened candidates or significant events between the two selected comparisons."
            )
            filtered_summary = self._filter_frame_by_value(program_summary, "pair_id", pair_id)
            filtered_events = self._filter_frame_by_value(program_events, "pair_id", pair_id)
            cross_matrix = self.project_service.run_state.results.cross_comparison_candidate_matrix
        if mode == "matrix":
            cross_files = self._cross_file_map()
            self._draw_cross_candidate_membership_summary(
                cross_matrix,
                f"Candidate matrix: {pair_name or 'selected pair'}",
                "Pairwise candidate matrix summarises which genes are candidates in comparison A, comparison B, or both.",
            )
            preferred_columns = [
                "pair_name",
                "comparison_A_display_name",
                "comparison_B_display_name",
                "gene_id",
                "gene_symbol",
                "candidate_pattern",
                "pattern_reason",
                "candidate_comparisons",
                "pair_membership_label",
                "pair_membership_count",
                "is_candidate_A",
                "is_candidate_B",
                "tier_A",
                "tier_B",
                "rank_A",
                "rank_B",
                "evidence_class_A",
                "evidence_class_B",
                "direction_class_A",
                "direction_class_B",
                "standardized_log2FC_A",
                "standardized_log2FC_B",
                "dPSI_A",
                "dPSI_B",
                "DE_FDR_A",
                "DE_FDR_B",
                "rMATS_FDR_A",
                "rMATS_FDR_B",
                "event_type_A",
                "event_type_B",
                "pair_gene_level_direction_changed",
                "pair_DEG_direction_reversed",
                "pair_dominant_AS_direction_reversed",
                "pair_tier_changed",
                "pair_evidence_class_changed",
                "best_pair_rank",
            ]
            self._populate_table(
                cross_matrix,
                preferred_columns,
                file_path=cross_files.get("matrix"),
                output_dir=self.project_service.cross_comparison_output_dir(),
                preview_name=f"Pairwise candidate matrix - {pair_name or 'selected pair'}",
            )
            self._set_detail_text(
                module_detail("5.4.matrix")
                + "\n\nPair-specific notes:\n"
                "Level: gene-level candidate membership, not event-level AS interaction modeling.\n"
                f"Scope: selected comparison pair only: {pair_name or 'comparison A vs comparison B'}.\n"
                "Function: this is the base A-vs-B membership matrix. One row = one unique gene that is a candidate in comparison A and/or comparison B.\n"
                "What 5.4.1 shows:\n"
                "- is_candidate_A / is_candidate_B: whether the gene is a candidate in each side of the selected pair\n"
                "- tier_A / tier_B, rank_A / rank_B, evidence_class_A / evidence_class_B and direction_class_A / direction_class_B\n"
                "- pair_membership_label: A only, B only, or shared A+B\n"
                "- pair_membership_count: 1 means A-only/B-only; 2 means shared by both sides of this pair. This replaces the confusing '1 comparison / 2 comparison' chart wording.\n"
                "- pair_gene_level_direction_changed: gene-level DEG and/or dominant AS direction changed between A and B. This is not strict same-event AS reversal.\n"
                "- pair_DEG_direction_reversed: standardized_log2FC changes sign between A and B.\n"
                "- pair_dominant_AS_direction_reversed: dominant_rMATS_standardized_dPSI changes sign between A and B. The dominant event may differ by comparison, so use 5.4.3/5.4.4 for strict same-event AS reversal.\n"
                "This page is not global anymore. Change the pair tab to inspect a different KO-effect / DP-effect style comparison.\n"
                "This page is also not an event-level AS test; use 5.4.3 for pair-specific AS-event differences with replicate read support.\n"
                "Rank meaning:\n"
                "- rank_A / rank_B = the original within-comparison numeric rank from 5.3.1\n"
                "- best_pair_rank = the best (lowest) rank observed across comparison A and B\n"
                "In other words, rank here is inherited from 5.3.1 and is not a new global rank.",
                "5.4.1 pairwise candidate matrix loaded. Open Detail for membership and rank definitions.",
            )
        elif mode == "comparison":
            self._draw_program_summary(filtered_summary, pair_name)
            self._populate_table(
                filtered_events,
                [
                    "pair_name",
                    "comparison_A_name",
                    "comparison_B_name",
                    "gene_symbol",
                    "gene_id",
                    "event_type",
                    "class_label",
                    "dPSI_A",
                    "dPSI_B",
                    "FDR_A",
                    "FDR_B",
                    "abs_delta_between_programs",
                ],
                output_dir=self.project_service.cross_comparison_output_dir(),
                preview_name="Group comparison event table",
            )
            self.details.setPlainText(
                "Group comparison is the candidate-aware event matrix for one selected comparison pair.\n"
                "Use it to inspect which genes/events are shared, pair-specific, or strongly shifted between comparison A and comparison B."
            )
        elif mode == "shared":
            shared = cross_matrix.loc[
                cross_matrix["candidate_pattern"].astype(str).isin(
                    [
                        "shared_pair",
                        "shared_DEG_direction_reversed",
                        "shared_dominant_AS_direction_reversed",
                        "shared_mixed_direction_changed",
                        "shared_tier_changed",
                        "shared_evidence_class_changed",
                        "comparison_A_specific_lost_in_B",
                        "comparison_B_specific_gained_in_B",
                    ]
                )
            ].copy() if cross_matrix is not None and not cross_matrix.empty else pd.DataFrame()
            self._draw_cross_candidate_pattern_breakdown(
                shared,
                f"Shared / Specific / Gained / Lost: {pair_name or 'selected pair'}",
                "Pairwise pattern breakdown derived from the selected A-vs-B candidate matrix.",
            )
            cross_files = self._cross_file_map()
            preferred_columns = [
                "pair_name",
                "comparison_A_display_name",
                "comparison_B_display_name",
                "gene_id",
                "gene_symbol",
                "candidate_pattern",
                "pattern_reason",
                "candidate_comparisons",
                "pair_membership_label",
                "pair_membership_count",
                "is_candidate_A",
                "is_candidate_B",
                "tier_A",
                "tier_B",
                "rank_A",
                "rank_B",
                "evidence_class_A",
                "evidence_class_B",
                "direction_class_A",
                "direction_class_B",
                "standardized_log2FC_A",
                "standardized_log2FC_B",
                "dPSI_A",
                "dPSI_B",
                "DE_FDR_A",
                "DE_FDR_B",
                "rMATS_FDR_A",
                "rMATS_FDR_B",
                "event_type_A",
                "event_type_B",
                "pair_gene_level_direction_changed",
                "pair_DEG_direction_reversed",
                "pair_dominant_AS_direction_reversed",
                "pair_tier_changed",
                "pair_evidence_class_changed",
                "best_pair_rank",
            ]
            self._populate_table(
                shared,
                preferred_columns,
                file_path=cross_files.get("matrix"),
                output_dir=self.project_service.cross_comparison_output_dir(),
                preview_name=f"Pairwise gene pattern view - {pair_name or 'selected pair'}",
            )
            self._set_detail_text(
                module_detail("5.4.pattern")
                + "\n\nPair-specific notes:\n"
                "Level: gene-level candidate-pattern interpretation, not event-level AS interaction modeling.\n"
                f"Scope: selected comparison pair only: {pair_name or 'comparison A vs comparison B'}.\n"
                "Function: 5.4.2 is the pairwise pattern-interpretation page. It asks whether a gene is shared by A/B, A-specific, B-specific, or changed in direction/tier/evidence class between A and B.\n"
                "Definitions:\n"
                "- shared_pair: candidate in both comparison A and comparison B\n"
                "- comparison_A_specific_lost_in_B: candidate in A only; 'lost' is relative to the selected A -> B pair order, not automatically a time-course claim\n"
                "- comparison_B_specific_gained_in_B: candidate in B only; 'gained' is relative to the selected A -> B pair order\n"
                "- shared_DEG_direction_reversed: candidate in both A and B, and standardized_log2FC changes sign between A and B\n"
                "- shared_dominant_AS_direction_reversed: candidate in both A and B, and dominant_rMATS_standardized_dPSI changes sign between A and B\n"
                "- shared_mixed_direction_changed: candidate in both A and B, and both DEG and dominant AS directions change sign. This is gene-level direction change, not strict same-event AS reversal\n"
                "- shared_tier_changed: candidate in both A and B, but candidate tier differs\n"
                "- shared_evidence_class_changed: candidate in both A and B, but evidence class differs\n"
                "Important: 5.4.2 direction classes are gene-level candidate interpretation. For strict same AS-event reversal, use 5.4.3/5.4.4 where event_matching_key is fixed.\n"
                "Why 5.4.1 and 5.4.2 can look similar: both use the same selected-pair source genes, but 5.4.1 is raw A/B membership while 5.4.2 groups genes into interpretable pairwise classes.\n"
                "Rank still means the original within-comparison numeric rank from 5.3.1.\n"
                "best_pair_rank is not a new effect size; it is simply the best (lowest) within-comparison candidate rank across A and B.",
                "5.4.2 pairwise pattern view loaded. Open Detail for shared/specific/gained/lost definitions.",
            )
        elif mode == "venn":
            abs_cutoff = float(self.group_venn_sig_cutoff.value())
            min_reads = int(self.group_venn_min_reads.value())
            signature = self._group_venn_signature()
            cached = self._group_venn_cache.get(signature, {})
            summary = cached.get("summary", pd.DataFrame())
            events = cached.get("events", pd.DataFrame())
            self.context_header.setText(
                f"Comparison pair: {pair_name or 'selected pair'} | Event-level Venn for significant AS events."
            )
            self._set_explanation(
                "5.4.4 shows event-level Venn overlap only: A-specific significant AS events, shared significant AS events, and B-specific significant AS events. "
                "The read filter applies only to the A-only/B-only sides by checking the counterpart comparison's junction read support."
            )
            self._draw_cross_venn_summary(summary, pair_name)
            event_columns = [
                "event_matching_key",
                "gene",
                "gene_symbol",
                "gene_id",
                "event_type",
                "event_id",
                "coordinates",
                "event_class",
                "comparison_A_display_name",
                "comparison_B_display_name",
                "dPSI_A",
                "FDR_A",
                "dPSI_B",
                "FDR_B",
                "direction_A",
                "direction_B",
                "mean_total_reads_A",
                "mean_total_reads_B",
                "counterpart_read_support",
            ]
            summary_columns = [
                "level",
                "comparison_A_display_name",
                "comparison_B_display_name",
                "exclusive_filter_mode",
                "A_only",
                "shared_A_and_B",
                "B_only",
                "total_union",
                "intersection_fraction",
                "abs_dPSI_cutoff_used",
                "min_counterpart_reads_used",
                "source_cache",
                "notes",
            ]
            preview_frame = events if isinstance(events, pd.DataFrame) and not events.empty else summary
            preview_columns = event_columns if isinstance(events, pd.DataFrame) and not events.empty else summary_columns
            self._populate_table(
                preview_frame if isinstance(preview_frame, pd.DataFrame) else pd.DataFrame(),
                [column for column in preview_columns if isinstance(preview_frame, pd.DataFrame) and column in preview_frame.columns],
                output_dir=self.project_service.cross_comparison_output_dir(),
                preview_name=f"5.4.4 Venn / overlap - {pair_name or 'selected pair'}",
            )
            cache_note = ""
            if not cached:
                cache_note = (
                    "\n\nNo cached 5.4.4 event-overlap result is loaded for this pair. "
                    "Click Run 5.4.4 Event Venn to build only this current pair's event-level overlap."
                )
            self._set_detail_text(
                module_detail("5.4.venn")
                + f"\n\nCurrent pair: {pair_name or 'selected pair'}.\n"
                f"Significant AS definition: rMATS_FDR <= 0.05 and abs(standardized dPSI) >= {abs_cutoff:.3f}.\n"
                f"Counterpart read support threshold: mean IJC+SJC reads >= {min_reads}.\n"
                "This page is event-level only. It no longer mixes 5.4.1 gene candidates into the Venn counts.\n"
                "Shared = significant in both comparisons and therefore modeled in 5.4.3.\n"
                "A only = significant only in comparison A.\n"
                "B only = significant only in comparison B.\n"
                f"Exclusive-side filter = {self.group_venn_filter.currentText()}.\n"
                "When the filter is read-supported / low-power / counterpart-not-found, it only subsets the A-only and B-only sides using counterpart read-support classification; the shared-significant center remains unchanged.\n"
                "The Venn counts are descriptive overlap counts. They are not statistical tests and they do not replace 5.4.3 interaction modeling."
                + cache_note,
                "5.4.4 event-level Venn loaded. Open Detail for significant-event and read-support interpretation.",
            )
        elif mode in {"venn_multi", "venn_multi_similarity"}:
            signature = self._group_venn_multi_signature()
            cached = self._group_venn_multi_cache if self._group_venn_multi_cache.get("signature") == signature else {}
            summary = cached.get("summary", pd.DataFrame())
            events = cached.get("events", pd.DataFrame())
            selected_names = cached.get("selected_names", [])
            is_similarity_page = mode == "venn_multi_similarity"
            self.context_header.setText(
                "5.4.6 selected-comparison direction-aware AS-event similarity heatmap."
                if is_similarity_page
                else "5.4.5 selected-comparison direction-aware AS-event Venn."
            )
            self._set_explanation(
                ("5.4.6 compares" if is_similarity_page else "5.4.5 compares")
                + " strict rMATS-significant AS events across user-selected comparisons. "
                "Venn and similarity are direction-aware: shared requires the same event_matching_key and the same standardized dPSI sign. Opposite-direction shared identities are flagged in the table."
            )
            self._draw_group_venn_multi(
                summary if isinstance(summary, pd.DataFrame) else pd.DataFrame(),
                selected_names if isinstance(selected_names, list) else [],
                events if isinstance(events, pd.DataFrame) else pd.DataFrame(),
            )
            event_columns = [
                "event_matching_key",
                "gene",
                "gene_symbol",
                "event_type",
                "event_id",
                "coordinates",
                "n_comparisons_significant",
                "comparison_membership",
                "membership_pattern",
                "direction_pattern",
                "has_opposite_direction",
            ]
            if isinstance(events, pd.DataFrame) and not events.empty:
                for name in selected_names if isinstance(selected_names, list) else []:
                    safe = self._safe_column_label(str(name))
                    event_columns.extend([f"present_{safe}", f"direction_{safe}", f"dPSI_{safe}", f"FDR_{safe}"])
            self._populate_table(
                events if isinstance(events, pd.DataFrame) else pd.DataFrame(),
                [column for column in event_columns if isinstance(events, pd.DataFrame) and column in events.columns],
                output_dir=self.project_service.cross_comparison_output_dir(),
                preview_name="5.4.6 selected-comparison direction-aware AS-event similarity" if is_similarity_page else "5.4.5 selected-comparison direction-aware AS-event Venn",
            )
            if not cached:
                self._set_detail_text(
                    module_detail("5.4.venn_multi")
                    + (
                        "\n\nNo cached 5.4.6 similarity result is loaded. Comparisons are auto-loaded from Pairing / selected comparisons; tick the comparisons you want, then click Run Direction-aware Similarity. Opening this page does not auto-scan rMATS files."
                        if is_similarity_page
                        else "\n\nNo cached 5.4.5 Venn result is loaded. Comparisons are auto-loaded from Pairing / selected comparisons; tick the comparisons you want, then click Run Direction-aware Venn. Opening this page does not auto-scan rMATS files."
                    ),
                    "5.4.6 not run. Click Run Direction-aware Similarity." if is_similarity_page else "5.4.5 not run. Click Run Direction-aware Venn.",
                )
            else:
                self._set_detail_text(
                    module_detail("5.4.venn_multi")
                    + f"\n\nSelected comparisons: {', '.join(str(x) for x in selected_names)}\n"
                    f"Significant AS definition: rMATS_FDR <= 0.05 and abs(standardized dPSI) >= {self.group_venn_multi_sig_cutoff.value():.3f}.\n"
                    f"Minimum mean IJC+SJC reads loaded for event support context: {self.group_venn_multi_min_reads.value()}.\n"
                    "This is direction-aware event-level overlap only. It is descriptive and does not run the 5.4.3 interaction model.",
                    "5.4.6 selected-comparison similarity loaded." if is_similarity_page else "5.4.5 selected-comparison event Venn loaded.",
                )
        elif mode == "reversal":
            reversed_only = cross_matrix.loc[cross_matrix["candidate_pattern"].astype(str) == "direction_reversed"].copy() if cross_matrix is not None and not cross_matrix.empty else pd.DataFrame()
            self._draw_cross_candidate_summary(
                reversed_only,
                f"Direction reversal: {pair_name}",
                "Direction reversal view shows genes whose standardized transcript or dominant splicing direction changes sign across comparisons.",
            )
            cross_files = self._cross_file_map()
            self._populate_table(reversed_only, [column for column in reversed_only.columns if column in [
                "gene_id", "gene_symbol", "candidate_pattern", "pattern_reason", "candidate_comparisons", "best_overall_tier", "best_overall_rank", "has_direction_reversal"
            ] or column.startswith("direction_class_")], file_path=cross_files.get("matrix"), output_dir=self.project_service.cross_comparison_output_dir(), preview_name="Direction reversal candidates")
            self.details.setPlainText(
                "Direction reversal means cross-comparison biological reversal only.\n"
                "Level: gene-level candidate direction summary. It is not a replicate-level interaction test.\n"
                "It does not refer to Pairing/config flip. It means the standardized transcript direction and/or dominant rMATS splicing direction changes sign between comparisons.\n"
                "The direction is inherited from the standardized 5.3 candidate table: DEG uses standardized_log2FC and AS uses standardized/dominant dPSI. "
                "Use 5.4.3 when you need event-level AS interaction modeling with replicate read support."
            )
        elif mode == "heatmap":
            self._draw_program_heatmap(filtered_events, pair_name)
        else:
            self._draw_direction_strength_summary(filtered_events, pair_name)

    def _draw_cross_candidate_summary(self, frame: pd.DataFrame | None, title: str, explanation: str) -> None:
        self._scatter_collection = None
        self._scatter_frame = pd.DataFrame()
        self._pie_patches = []
        self._heatmap_image = None
        self._heatmap_rows = []
        self._heatmap_cols = []
        self._heatmap_values = None
        self._set_explanation(explanation)
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        self._set_axis_title(axis, title)
        if frame is None or frame.empty:
            axis.axis("off")
            axis.text(0.5, 0.5, "No rows available for this cross-comparison view.", ha="center", va="center")
            self.canvas.draw_idle()
            return
        working = frame.copy()
        summary = (
            working.get("candidate_pattern", pd.Series(dtype="object"))
            .fillna("unclassified")
            .astype(str)
            .value_counts()
            .sort_values(ascending=False)
        )
        bars = axis.bar(summary.index.tolist(), summary.values.tolist(), color="#4C78A8")
        axis.set_ylabel("Genes")
        axis.tick_params(axis="x", rotation=22)
        self._bar_patches = [
            (patch, f"{label}\nGenes: {value}")
            for patch, label, value in zip(bars, summary.index.tolist(), summary.values.tolist(), strict=False)
        ]
        self._finalize_axes_layout(rotate_x=True, bottom=0.28)
        self.canvas.draw_idle()

    def _draw_cross_candidate_membership_summary(self, frame: pd.DataFrame | None, title: str, explanation: str) -> None:
        self._scatter_collection = None
        self._scatter_frame = pd.DataFrame()
        self._pie_patches = []
        self._heatmap_image = None
        self._heatmap_rows = []
        self._heatmap_cols = []
        self._heatmap_values = None
        self._set_explanation(explanation)
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        self._set_axis_title(axis, title)
        if frame is None or frame.empty:
            axis.axis("off")
            axis.text(0.5, 0.5, "No rows available for this cross-comparison view.", ha="center", va="center")
            self.canvas.draw_idle()
            return
        working = frame.copy()
        if "pair_membership_label" in working.columns:
            summary = working["pair_membership_label"].fillna("unclassified").astype(str).value_counts()
            order = [label for label in ["A only", "B only", "shared A+B"] if label in summary.index]
            order.extend(label for label in summary.index.tolist() if label not in order)
            summary = summary.loc[order]
            labels = summary.index.tolist()
        else:
            membership = pd.to_numeric(
                working.get("n_comparisons_candidate", pd.Series(dtype="float")),
                errors="coerce",
            ).dropna()
            if membership.empty:
                axis.axis("off")
                axis.text(0.5, 0.5, "No comparison-membership counts available.", ha="center", va="center")
                self.canvas.draw_idle()
                return
            summary = membership.astype(int).value_counts().sort_index()
            labels = [f"{count} comparison{'s' if count != 1 else ''}" for count in summary.index.tolist()]
        if summary.empty:
            axis.axis("off")
            axis.text(0.5, 0.5, "No comparison-membership counts available.", ha="center", va="center")
            self.canvas.draw_idle()
            return
        bars = axis.bar(labels, summary.values.tolist(), color="#4C78A8")
        axis.set_ylabel("Genes")
        axis.set_xlabel("Pairwise candidate membership")
        self._bar_patches = [
            (patch, f"{label}\nGenes: {value}")
            for patch, label, value in zip(bars, labels, summary.values.tolist(), strict=False)
        ]
        self._finalize_axes_layout(rotate_x=False, bottom=0.22)
        self.canvas.draw_idle()

    def _draw_cross_candidate_pattern_breakdown(self, frame: pd.DataFrame | None, title: str, explanation: str) -> None:
        self._scatter_collection = None
        self._scatter_frame = pd.DataFrame()
        self._pie_patches = []
        self._heatmap_image = None
        self._heatmap_rows = []
        self._heatmap_cols = []
        self._heatmap_values = None
        self._set_explanation(explanation)
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        self._set_axis_title(axis, title)
        if frame is None or frame.empty:
            axis.axis("off")
            axis.text(0.5, 0.5, "No rows available for this cross-comparison pattern view.", ha="center", va="center")
            self.canvas.draw_idle()
            return
        working = frame.copy()
        pattern_order = [
            "shared_pair",
            "shared_DEG_direction_reversed",
            "shared_dominant_AS_direction_reversed",
            "shared_mixed_direction_changed",
            "shared_tier_changed",
            "shared_evidence_class_changed",
            "comparison_A_specific_lost_in_B",
            "comparison_B_specific_gained_in_B",
            "shared_all",
            "shared_subset",
            "comparison_specific",
            "gained",
            "lost",
            "tier_changed",
            "evidence_class_changed",
            "DEG_only_in_one_condition",
            "splicing_only_in_one_condition",
        ]
        pattern_colors = {
            "shared_pair": "#4C78A8",
            "shared_DEG_direction_reversed": "#E45756",
            "shared_dominant_AS_direction_reversed": "#B279A2",
            "shared_mixed_direction_changed": "#D62728",
            "shared_tier_changed": "#B279A2",
            "shared_evidence_class_changed": "#FF9DA6",
            "comparison_A_specific_lost_in_B": "#F58518",
            "comparison_B_specific_gained_in_B": "#54A24B",
            "shared_all": "#4C78A8",
            "shared_subset": "#72B7B2",
            "comparison_specific": "#F58518",
            "gained": "#54A24B",
            "lost": "#E45756",
            "tier_changed": "#B279A2",
            "evidence_class_changed": "#FF9DA6",
            "DEG_only_in_one_condition": "#9D755D",
            "splicing_only_in_one_condition": "#BAB0AC",
        }
        summary = (
            working.get("candidate_pattern", pd.Series(dtype="object"))
            .fillna("unclassified")
            .astype(str)
            .value_counts()
        )
        labels = [label for label in pattern_order if label in summary.index]
        if not labels:
            axis.axis("off")
            axis.text(0.5, 0.5, "No candidate-pattern rows available.", ha="center", va="center")
            self.canvas.draw_idle()
            return
        values = [int(summary[label]) for label in labels]
        colors = [pattern_colors.get(label, "#4C78A8") for label in labels]
        bars = axis.bar(labels, values, color=colors)
        axis.set_ylabel("Genes")
        axis.set_xlabel("Candidate pattern")
        self._bar_patches = [
            (patch, f"{label}\nGenes: {value}")
            for patch, label, value in zip(bars, labels, values, strict=False)
        ]
        self._finalize_axes_layout(rotate_x=True, bottom=0.33)
        self.canvas.draw_idle()

    def _draw_cross_as_pattern_summary(
        self,
        tables: dict[str, pd.DataFrame],
        pair_name: str | None,
        abs_cutoff: float,
        delta_cutoff: float,
    ) -> None:
        self._scatter_collection = None
        self._scatter_frame = pd.DataFrame()
        self._pie_patches = []
        self._heatmap_image = None
        self._heatmap_rows = []
        self._heatmap_cols = []
        self._heatmap_values = None
        self._bar_patches = []
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        self._set_axis_title(axis, f"Legacy summary-only dPSI patterns: {pair_name or 'selected pair'}")
        counts = {
            "Same dir + large delta": len(tables.get("same_direction_large_delta", pd.DataFrame())),
            "Opposite direction": len(tables.get("opposite_direction", pd.DataFrame())),
            "All shared significant": len(tables.get("all_shared", pd.DataFrame())),
            "A-only / B-only": len(tables.get("condition_specific", pd.DataFrame())),
        }
        if not any(counts.values()):
            axis.axis("off")
            axis.text(0.5, 0.5, "No events found", ha="center", va="center")
            self.canvas.draw_idle()
            return
        bars = axis.bar(
            list(counts.keys()),
            list(counts.values()),
            color=["#2563EB", "#7C3AED", "#0F766E", "#D97706"],
        )
        axis.set_ylabel("Events")
        axis.tick_params(axis="x", rotation=18)
        self._bar_patches = [
            (patch, f"{label}\nEvents: {value}\nabs(dPSI) cutoff: {abs_cutoff:.3f}\nlarge delta cutoff: {delta_cutoff:.3f}")
            for patch, (label, value) in zip(bars, counts.items(), strict=False)
        ]
        self._finalize_axes_layout(rotate_x=True, bottom=0.28)
        self.canvas.draw_idle()

    def _group_venn_signature(self) -> tuple[str | None, float, int, str]:
        return (
            self._current_tab_key(self.group_tabs),
            round(float(self.group_venn_sig_cutoff.value()), 6),
            int(self.group_venn_min_reads.value()),
            str(self.group_venn_filter.currentData() or "all_significant"),
        )

    def _build_pair_venn_tables(self, pair_id: str | None) -> tuple[pd.DataFrame, pd.DataFrame]:
        pair = self._pair_definition_from_id(pair_id)
        abs_cutoff = float(self.group_venn_sig_cutoff.value())
        min_reads = int(self.group_venn_min_reads.value())
        venn_filter = str(self.group_venn_filter.currentData() or "all_significant")
        if pair is None or not pair.comparison_a or not pair.comparison_b:
            return pd.DataFrame(), pd.DataFrame()

        comp_a = str(pair.comparison_a)
        comp_b = str(pair.comparison_b)
        name_a = self._comparison_name_from_id(comp_a) or comp_a
        name_b = self._comparison_name_from_id(comp_b) or comp_b
        frame_a = self.project_service._cross_as_interaction_source_frame(comp_a, abs_cutoff, min_reads, significant_only=True)
        frame_b = self.project_service._cross_as_interaction_source_frame(comp_b, abs_cutoff, min_reads, significant_only=True)
        keys_a = set(frame_a.get("event_matching_key", pd.Series(dtype="object")).dropna().astype(str).tolist())
        keys_b = set(frame_b.get("event_matching_key", pd.Series(dtype="object")).dropna().astype(str).tolist())
        shared_keys = keys_a & keys_b
        a_only_keys = keys_a - keys_b
        b_only_keys = keys_b - keys_a
        counterpart_b = (
            self.project_service._cross_as_interaction_source_frame(comp_b, abs_cutoff, min_reads, significant_only=False, candidate_keys=a_only_keys)
            if a_only_keys
            else pd.DataFrame()
        )
        counterpart_a = (
            self.project_service._cross_as_interaction_source_frame(comp_a, abs_cutoff, min_reads, significant_only=False, candidate_keys=b_only_keys)
            if b_only_keys
            else pd.DataFrame()
        )

        def by_key(frame: pd.DataFrame) -> dict[str, pd.Series]:
            if frame.empty or "event_matching_key" not in frame.columns:
                return {}
            result: dict[str, pd.Series] = {}
            for _, row in frame.iterrows():
                key = str(row.get("event_matching_key", ""))
                if key and key not in result:
                    result[key] = row
            return result

        lookup_a = by_key(frame_a)
        lookup_b = by_key(frame_b)
        lookup_counterpart_a = by_key(counterpart_a)
        lookup_counterpart_b = by_key(counterpart_b)

        def support_label(row: pd.Series | None) -> str:
            if row is None:
                return "counterpart_event_not_found"
            mean_reads = pd.to_numeric(row.get("mean_total_reads"), errors="coerce")
            replicate_count = pd.to_numeric(row.get("replicate_count"), errors="coerce")
            if pd.notna(mean_reads) and float(mean_reads) >= min_reads and pd.notna(replicate_count) and int(replicate_count) >= 2:
                return "supported_non_significant"
            return "low_power_or_absent"

        def should_keep_exclusive(label: str) -> bool:
            if venn_filter == "supported_only":
                return label == "supported_non_significant"
            if venn_filter == "low_power_only":
                return label == "low_power_or_absent"
            if venn_filter == "counterpart_not_found_only":
                return label == "counterpart_event_not_found"
            return True

        def base_fields(primary: pd.Series | None, counterpart: pd.Series | None, key: str) -> dict[str, object]:
            source = primary if primary is not None else counterpart
            source = source if source is not None else pd.Series(dtype="object")
            return {
                "event_matching_key": key,
                "gene": source.get("gene", ""),
                "gene_symbol": source.get("gene_symbol", ""),
                "gene_id": source.get("gene_id", ""),
                "event_type": source.get("event_type", ""),
                "event_id": source.get("event_id", ""),
                "coordinates": source.get("coordinates", ""),
            }

        rows: list[dict[str, object]] = []

        def append_row(key: str, event_class: str, row_a: pd.Series | None, row_b: pd.Series | None, counterpart_support: str) -> None:
            event_row = base_fields(row_a, row_b, key)
            event_row.update(
                {
                    "event_class": event_class,
                    "comparison_A": comp_a,
                    "comparison_B": comp_b,
                    "comparison_A_display_name": name_a,
                    "comparison_B_display_name": name_b,
                    "dPSI_A": row_a.get("dPSI", pd.NA) if row_a is not None else pd.NA,
                    "FDR_A": row_a.get("FDR", pd.NA) if row_a is not None else pd.NA,
                    "dPSI_B": row_b.get("dPSI", pd.NA) if row_b is not None else pd.NA,
                    "FDR_B": row_b.get("FDR", pd.NA) if row_b is not None else pd.NA,
                    "direction_A": row_a.get("direction", "") if row_a is not None else "",
                    "direction_B": row_b.get("direction", "") if row_b is not None else "",
                    "mean_total_reads_A": row_a.get("mean_total_reads", pd.NA) if row_a is not None else pd.NA,
                    "mean_total_reads_B": row_b.get("mean_total_reads", pd.NA) if row_b is not None else pd.NA,
                    "counterpart_read_support": counterpart_support,
                }
            )
            rows.append(event_row)

        for key in sorted(shared_keys):
            append_row(key, "shared_significant", lookup_a.get(key), lookup_b.get(key), "shared_significant")
        for key in sorted(a_only_keys):
            label = support_label(lookup_counterpart_b.get(key))
            if should_keep_exclusive(label):
                append_row(key, "comparison_A_specific", lookup_a.get(key), lookup_counterpart_b.get(key), label)
        for key in sorted(b_only_keys):
            label = support_label(lookup_counterpart_a.get(key))
            if should_keep_exclusive(label):
                append_row(key, "comparison_B_specific", lookup_counterpart_a.get(key), lookup_b.get(key), label)

        events = pd.DataFrame(rows)
        if not events.empty:
            events = events.sort_values(["event_class", "gene", "event_type", "FDR_A", "FDR_B"], na_position="last").reset_index(drop=True)
        a_only_count = int(events["event_class"].eq("comparison_A_specific").sum()) if "event_class" in events.columns else 0
        b_only_count = int(events["event_class"].eq("comparison_B_specific").sum()) if "event_class" in events.columns else 0
        shared_count = int(events["event_class"].eq("shared_significant").sum()) if "event_class" in events.columns else 0
        total_events = a_only_count + shared_count + b_only_count
        filter_note = {
            "all_significant": "All significant AS events; A-only and B-only are not filtered by counterpart reads",
            "supported_only": f"A-only/B-only filtered to counterpart events with mean IJC+SJC reads >= {min_reads}",
            "low_power_only": f"A-only/B-only filtered to counterpart events with mean IJC+SJC reads < {min_reads} or absent reads",
            "counterpart_not_found_only": "A-only/B-only filtered to events not found in the counterpart comparison",
        }.get(venn_filter, "All significant AS events")
        summary = pd.DataFrame(
            [
                {
                    "level": "Significant AS events",
                    "pair_id": pair_id or "",
                    "pair_name": pair.resolved_name,
                    "comparison_A": comp_a,
                    "comparison_B": comp_b,
                    "comparison_A_display_name": name_a,
                    "comparison_B_display_name": name_b,
                    "exclusive_filter_mode": filter_note,
                    "A_only": a_only_count,
                    "shared_A_and_B": shared_count,
                    "B_only": b_only_count,
                    "total_union": total_events,
                    "intersection_fraction": (shared_count / total_events) if total_events else 0.0,
                    "abs_dPSI_cutoff_used": abs_cutoff,
                    "min_counterpart_reads_used": min_reads,
                    "source_cache": "5.4.4 direct rMATS strict significant AS-event overlap",
                    "notes": "Event-level Venn only. Shared = same event key significant in both comparisons. A only and B only are kept separate.",
                }
            ]
        )
        output_dir = self.project_service.cross_comparison_output_dir()
        if output_dir is not None:
            pair_folder = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(pair_id or "selected_pair")).strip("_") or "selected_pair"
            target_dir = output_dir / "event_venn" / pair_folder
            target_dir.mkdir(parents=True, exist_ok=True)
            summary.to_csv(target_dir / "event_venn_summary.tsv", sep="\t", index=False)
            events.to_csv(target_dir / "event_venn_events.tsv", sep="\t", index=False)
        return summary, events

    def _build_cross_venn_summary(
        self,
        interaction_tables: dict[str, pd.DataFrame] | None,
        pair_name: str | None,
        venn_filter: str,
        min_reads: int,
    ) -> pd.DataFrame:
        pair_id = self._current_tab_key(self.group_tabs)
        pair = self._pair_definition_from_id(pair_id)
        name_a = self._comparison_name_from_id(pair.comparison_a) if pair is not None else None
        name_b = self._comparison_name_from_id(pair.comparison_b) if pair is not None else None
        name_a = name_a or (pair.comparison_a if pair is not None else "comparison A")
        name_b = name_b or (pair.comparison_b if pair is not None else "comparison B")

        rows: list[dict[str, object]] = []

        interaction_tables = interaction_tables or {}
        shared_events = interaction_tables.get("shared_interaction", pd.DataFrame())
        exclusive = interaction_tables.get("exclusive", pd.DataFrame()).copy()
        filter_note = "All significant AS events; A-only and B-only are not filtered by counterpart reads"
        if not exclusive.empty and "non_significant_side_read_support" in exclusive.columns:
            support = exclusive["non_significant_side_read_support"].fillna("").astype(str)
            if venn_filter == "supported_only":
                exclusive = exclusive.loc[support.str.startswith("supported_")].copy()
                filter_note = f"A-only/B-only filtered to counterpart events with mean IJC+SJC reads >= {min_reads}"
            elif venn_filter == "low_power_only":
                exclusive = exclusive.loc[support.str.startswith("low_power")].copy()
                filter_note = f"A-only/B-only filtered to counterpart events with mean IJC+SJC reads < {min_reads} or absent reads"
            elif venn_filter == "counterpart_not_found_only":
                exclusive = exclusive.loc[support.str.startswith("counterpart_event_not_found")].copy()
                filter_note = "A-only/B-only filtered to events not found in the counterpart comparison"
        a_only_events = 0
        b_only_events = 0
        if not exclusive.empty and "exclusive_side" in exclusive.columns:
            side = exclusive["exclusive_side"].fillna("").astype(str)
            a_only_events = int(side.eq("comparison_A_only").sum())
            b_only_events = int(side.eq("comparison_B_only").sum())
        shared_count = len(shared_events)
        total_events = a_only_events + shared_count + b_only_events
        rows.append(
            {
                "level": "Significant AS events",
                "pair_id": pair_id or "",
                "pair_name": pair_name or "",
                "comparison_A_display_name": name_a,
                "comparison_B_display_name": name_b,
                "exclusive_filter_mode": filter_note,
                "A_only": a_only_events,
                "shared_A_and_B": shared_count,
                "B_only": b_only_events,
                "total_union": total_events,
                "intersection_fraction": (shared_count / total_events) if total_events else 0.0,
                "abs_dPSI_cutoff_used": float(self.group_venn_sig_cutoff.value()),
                "min_counterpart_reads_used": min_reads,
                "source_cache": "5.4.3 AS interaction/read-support cache" if total_events else "No events found in current cache/result",
                "notes": "Event-level Venn only. Shared = significant in both comparisons. A only = significant in comparison A only. B only = significant in comparison B only. Optional filter uses counterpart read-support class.",
            }
        )
        return pd.DataFrame(rows)

    def _draw_single_venn_axis(
        self,
        axis,
        title: str,
        name_a: str,
        name_b: str,
        a_only: int,
        shared: int,
        b_only: int,
        note: str,
    ) -> None:
        axis.set_xlim(0, 1)
        axis.set_ylim(0, 1)
        axis.axis("off")
        left_color = "#2563EB"
        right_color = "#D97706"
        axis.add_patch(Circle((0.39, 0.52), 0.29, color=left_color, alpha=0.33, ec=left_color, lw=2))
        axis.add_patch(Circle((0.61, 0.52), 0.29, color=right_color, alpha=0.33, ec=right_color, lw=2))
        axis.text(0.5, 0.92, title, ha="center", va="center", fontsize=13, fontweight="bold")
        axis.text(0.27, 0.79, textwrap.shorten(str(name_a), width=30, placeholder="..."), ha="center", va="center", fontsize=9, color=left_color, fontweight="bold")
        axis.text(0.73, 0.79, textwrap.shorten(str(name_b), width=30, placeholder="..."), ha="center", va="center", fontsize=9, color=right_color, fontweight="bold")
        axis.text(0.24, 0.52, f"A only\n{a_only}", ha="center", va="center", fontsize=13, fontweight="bold")
        axis.text(0.50, 0.52, f"shared\n{shared}", ha="center", va="center", fontsize=13, fontweight="bold")
        axis.text(0.76, 0.52, f"B only\n{b_only}", ha="center", va="center", fontsize=13, fontweight="bold")
        total = a_only + shared + b_only
        fraction = shared / total if total else 0
        axis.text(0.5, 0.16, f"Union = {total} | Intersection fraction = {fraction:.2%}", ha="center", va="center", fontsize=10)
        axis.text(0.5, 0.07, note, ha="center", va="center", fontsize=8, color="#4B5563")

    def _draw_compact_venn_axis(
        self,
        axis,
        name_a: str,
        name_b: str,
        a_only: int,
        shared: int,
        b_only: int,
    ) -> None:
        axis.set_xlim(0, 1)
        axis.set_ylim(0, 1)
        axis.axis("off")
        left_color = "#2563EB"
        right_color = "#D97706"
        axis.add_patch(Circle((0.42, 0.54), 0.30, color=left_color, alpha=0.30, ec=left_color, lw=1.4))
        axis.add_patch(Circle((0.58, 0.54), 0.30, color=right_color, alpha=0.30, ec=right_color, lw=1.4))
        axis.text(0.5, 0.96, f"{name_a}\nvs\n{name_b}", ha="center", va="top", fontsize=7.5, fontweight="bold")
        axis.text(0.30, 0.54, f"A\n{a_only}", ha="center", va="center", fontsize=10, fontweight="bold")
        axis.text(0.50, 0.54, f"shared\n{shared}", ha="center", va="center", fontsize=10, fontweight="bold")
        axis.text(0.70, 0.54, f"B\n{b_only}", ha="center", va="center", fontsize=10, fontweight="bold")
        total = a_only + shared + b_only
        jaccard = shared / total if total else 0.0
        axis.text(0.5, 0.08, f"Jaccard {jaccard:.2f}", ha="center", va="center", fontsize=8, color="#374151")

    def _exact_venn_membership_counts(self, events: pd.DataFrame, selected_names: list[str]) -> dict[tuple[bool, ...], int]:
        if not isinstance(events, pd.DataFrame) or events.empty or not selected_names:
            return {}
        columns = [f"present_{self._safe_column_label(str(name))}" for name in selected_names]
        if any(column not in events.columns for column in columns):
            return {}
        counts: dict[tuple[bool, ...], int] = {}
        for mask_value in range(1, 2 ** len(columns)):
            mask = tuple(bool(mask_value & (1 << index)) for index in range(len(columns)))
            row_mask = pd.Series(True, index=events.index)
            for column, present in zip(columns, mask, strict=False):
                values = events[column].fillna(False).astype(bool)
                row_mask &= values if present else ~values
            counts[mask] = int(row_mask.sum())
        return counts

    def _direction_aware_venn_events(self, events: pd.DataFrame, selected_names: list[str]) -> pd.DataFrame:
        if not isinstance(events, pd.DataFrame) or events.empty or not selected_names:
            return pd.DataFrame()
        rows: list[dict[str, object]] = []
        for _, event_row in events.iterrows():
            key = str(event_row.get("event_matching_key", ""))
            if not key:
                continue
            direction_to_names: dict[str, list[str]] = {"positive": [], "negative": []}
            for name in selected_names:
                safe = self._safe_column_label(str(name))
                if not bool(event_row.get(f"present_{safe}", False)):
                    continue
                direction = str(event_row.get(f"direction_{safe}", self._dpsi_sign_label(event_row.get(f"dPSI_{safe}", pd.NA))))
                if direction in direction_to_names:
                    direction_to_names[direction].append(str(name))
            for direction, names in direction_to_names.items():
                if not names:
                    continue
                signed_row: dict[str, object] = {
                    "signed_event_matching_key": f"{key}::{direction}",
                    "event_matching_key": key,
                    "direction": direction,
                }
                for name in selected_names:
                    safe = self._safe_column_label(str(name))
                    signed_row[f"present_{safe}"] = str(name) in names
                rows.append(signed_row)
        return pd.DataFrame(rows)

    def _draw_three_set_venn_axis(self, axis, names: list[str], counts: dict[tuple[bool, ...], int]) -> None:
        axis.set_xlim(0, 1)
        axis.set_ylim(0, 1)
        axis.axis("off")
        colors = ["#2563EB", "#D97706", "#16A34A"]
        circles = [
            (0.40, 0.56, colors[0]),
            (0.60, 0.56, colors[1]),
            (0.50, 0.38, colors[2]),
        ]
        for x, y, color in circles:
            axis.add_patch(Circle((x, y), 0.28, color=color, alpha=0.28, ec=color, lw=2))
        axis.text(0.5, 0.94, "Direction-aware exact Venn across selected comparisons", ha="center", va="center", fontsize=13, fontweight="bold")
        axis.text(0.25, 0.83, textwrap.shorten(str(names[0]), width=24, placeholder="..."), ha="center", va="center", fontsize=9, color=colors[0], fontweight="bold")
        axis.text(0.75, 0.83, textwrap.shorten(str(names[1]), width=24, placeholder="..."), ha="center", va="center", fontsize=9, color=colors[1], fontweight="bold")
        axis.text(0.50, 0.08, textwrap.shorten(str(names[2]), width=28, placeholder="..."), ha="center", va="center", fontsize=9, color=colors[2], fontweight="bold")
        labels = {
            (True, False, False): ("A only", 0.25, 0.56),
            (False, True, False): ("B only", 0.75, 0.56),
            (False, False, True): ("C only", 0.50, 0.22),
            (True, True, False): ("A+B", 0.50, 0.66),
            (True, False, True): ("A+C", 0.38, 0.42),
            (False, True, True): ("B+C", 0.62, 0.42),
            (True, True, True): ("A+B+C", 0.50, 0.50),
        }
        for mask, (label, x, y) in labels.items():
            axis.text(x, y, f"{label}\n{counts.get(mask, 0)}", ha="center", va="center", fontsize=10.5, fontweight="bold")
        total = sum(counts.values())
        axis.text(0.5, 0.015, f"Union = {total}. Each signed event response is counted once.", ha="center", va="bottom", fontsize=8, color="#4B5563")

    def _draw_four_set_venn_axis(self, axis, names: list[str], counts: dict[tuple[bool, ...], int]) -> None:
        axis.set_xlim(0, 1)
        axis.set_ylim(0, 1)
        axis.axis("off")
        colors = ["#2563EB", "#D97706", "#16A34A", "#9333EA"]
        ellipses = [
            (0.42, 0.58, 0.62, 0.36, 28, colors[0]),
            (0.58, 0.58, 0.62, 0.36, -28, colors[1]),
            (0.42, 0.42, 0.62, 0.36, -28, colors[2]),
            (0.58, 0.42, 0.62, 0.36, 28, colors[3]),
        ]
        for x, y, width, height, angle, color in ellipses:
            axis.add_patch(Ellipse((x, y), width, height, angle=angle, color=color, alpha=0.22, ec=color, lw=2))
        axis.text(0.5, 0.97, "Direction-aware exact Venn across selected comparisons", ha="center", va="top", fontsize=13, fontweight="bold")
        legend_y = [0.90, 0.86, 0.82, 0.78]
        for index, name in enumerate(names[:4]):
            axis.text(
                0.04,
                legend_y[index],
                f"{chr(65 + index)} = {textwrap.shorten(str(name), width=32, placeholder='...')}",
                ha="left",
                va="center",
                fontsize=8.2,
                color=colors[index],
                fontweight="bold",
            )
        positions = {
            (True, False, False, False): ("A", 0.18, 0.65),
            (False, True, False, False): ("B", 0.82, 0.65),
            (False, False, True, False): ("C", 0.18, 0.35),
            (False, False, False, True): ("D", 0.82, 0.35),
            (True, True, False, False): ("AB", 0.50, 0.72),
            (True, False, True, False): ("AC", 0.30, 0.50),
            (True, False, False, True): ("AD", 0.43, 0.61),
            (False, True, True, False): ("BC", 0.57, 0.61),
            (False, True, False, True): ("BD", 0.70, 0.50),
            (False, False, True, True): ("CD", 0.50, 0.28),
            (True, True, True, False): ("ABC", 0.40, 0.54),
            (True, True, False, True): ("ABD", 0.60, 0.54),
            (True, False, True, True): ("ACD", 0.40, 0.46),
            (False, True, True, True): ("BCD", 0.60, 0.46),
            (True, True, True, True): ("ABCD", 0.50, 0.50),
        }
        for mask, (label, x, y) in positions.items():
            axis.text(x, y, f"{label}\n{counts.get(mask, 0)}", ha="center", va="center", fontsize=7.8, fontweight="bold")
        total = sum(counts.values())
        axis.text(0.5, 0.015, f"Union = {total}. Exact 4-set signed-response membership counts.", ha="center", va="bottom", fontsize=8, color="#4B5563")

    def _draw_cross_venn_summary(self, summary: pd.DataFrame, pair_name: str | None) -> None:
        self._scatter_collection = None
        self._scatter_frame = pd.DataFrame()
        self._pie_patches = []
        self._heatmap_image = None
        self._heatmap_rows = []
        self._heatmap_cols = []
        self._heatmap_values = None
        self._bar_patches = []
        self.figure.clear()
        if summary.empty:
            axis = self.figure.add_subplot(111)
            axis.axis("off")
            axis.text(
                0.5,
                0.5,
                "No Venn / overlap cache available.\nClick Run 5.4.4 Event Venn to build this pair's AS-event overlap cache.",
                ha="center",
                va="center",
            )
            self.canvas.draw_idle()
            return

        plot_rows = summary.loc[summary["level"].isin(["Gene candidates", "Significant AS events"])].copy()
        if plot_rows.empty:
            plot_rows = summary.head(1).copy()
        axes = self.figure.subplots(1, len(plot_rows), squeeze=False)[0]
        for axis, (_, row) in zip(axes, plot_rows.iterrows(), strict=False):
            def _count_value(column: str) -> int:
                value = pd.to_numeric(row.get(column, 0), errors="coerce")
                return 0 if pd.isna(value) else int(value)

            self._draw_single_venn_axis(
                axis,
                str(row.get("level", "Overlap")),
                str(row.get("comparison_A_display_name", "comparison A")),
                str(row.get("comparison_B_display_name", "comparison B")),
                _count_value("A_only"),
                _count_value("shared_A_and_B"),
                _count_value("B_only"),
                str(row.get("source_cache", "")),
            )
        self.figure.suptitle(f"5.4.4 Venn / overlap: {pair_name or 'selected pair'}", fontsize=16, fontweight="bold")
        self.figure.subplots_adjust(left=0.05, right=0.97, top=0.86, bottom=0.08, wspace=0.18)
        self.canvas.draw_idle()

    def _sync_group_venn_multi_comparison_checks(self) -> None:
        project = self.project_service.current_project
        comparisons = list(self.project_service.selected_comparisons_for_display()) if project is not None else []
        checked_ids: set[str] = set()
        existing_ids = set(self.group_venn_multi_checkboxes.keys())
        for comparison_id, checkbox in self.group_venn_multi_checkboxes.items():
            if checkbox.isChecked():
                checked_ids.add(comparison_id)

        self._syncing_group_venn_multi_checks = True
        while self.group_venn_multi_comparison_grid.count():
            item = self.group_venn_multi_comparison_grid.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.setParent(None)
        self.group_venn_multi_checkboxes = {}
        max_columns = 3
        for index, comparison in enumerate(comparisons):
            checkbox = QCheckBox(comparison.display_resolved_name)
            checkbox.setMinimumWidth(240)
            checkbox.setToolTip(
                f"{comparison.display_resolved_name}\n"
                f"comparison_id: {comparison.comparison_id}\n"
                f"resolved name: {comparison.resolved_name}"
            )
            should_check = True if not existing_ids else comparison.comparison_id in checked_ids
            checkbox.setChecked(should_check)
            checkbox.stateChanged.connect(self._on_group_venn_multi_selection_changed)
            self.group_venn_multi_checkboxes[comparison.comparison_id] = checkbox
            self.group_venn_multi_comparison_grid.addWidget(checkbox, index // max_columns, index % max_columns)
        self._syncing_group_venn_multi_checks = False

    def _checked_group_venn_multi_comparison_ids(self) -> list[str]:
        ids: list[str] = []
        for comparison_id, checkbox in self.group_venn_multi_checkboxes.items():
            if checkbox.isChecked():
                ids.append(comparison_id)
        return ids

    def _set_group_venn_multi_all_checked(self, checked: bool) -> None:
        self._sync_group_venn_multi_comparison_checks()
        self._syncing_group_venn_multi_checks = True
        for checkbox in self.group_venn_multi_checkboxes.values():
            checkbox.setChecked(checked)
        self._syncing_group_venn_multi_checks = False
        self._group_venn_multi_cache = {}
        self.refresh()

    def _on_group_venn_multi_selection_changed(self, _state: int | None = None) -> None:
        if self._syncing_group_venn_multi_checks:
            return
        self._group_venn_multi_cache = {}
        self.refresh()

    def _group_venn_multi_signature(self) -> tuple[float, int, str]:
        checked_ids = ",".join(self._checked_group_venn_multi_comparison_ids())
        return (
            round(float(self.group_venn_multi_sig_cutoff.value()), 6),
            int(self.group_venn_multi_min_reads.value()),
            checked_ids,
        )

    def _sync_group_venn_multi_pair_selector(self, summary: pd.DataFrame) -> None:
        previous_key = str(self.group_venn_multi_pair_selector.currentData() or "")
        self.group_venn_multi_pair_selector.blockSignals(True)
        self.group_venn_multi_pair_selector.clear()
        if isinstance(summary, pd.DataFrame) and not summary.empty:
            for _, row in summary.iterrows():
                comparison_a = str(row.get("comparison_A", ""))
                comparison_b = str(row.get("comparison_B", ""))
                name_a = str(row.get("comparison_A_display_name", comparison_a))
                name_b = str(row.get("comparison_B_display_name", comparison_b))
                pair_key = f"{comparison_a}||{comparison_b}"
                self.group_venn_multi_pair_selector.addItem(f"{name_a} vs {name_b}", pair_key)
            if previous_key:
                index = self.group_venn_multi_pair_selector.findData(previous_key)
                if index >= 0:
                    self.group_venn_multi_pair_selector.setCurrentIndex(index)
        self.group_venn_multi_pair_selector.setEnabled(self.group_venn_multi_pair_selector.count() > 1)
        self.group_venn_multi_pair_selector.blockSignals(False)

    def _selected_group_venn_multi_pair_row(self, summary: pd.DataFrame) -> pd.Series | None:
        if not isinstance(summary, pd.DataFrame) or summary.empty:
            return None
        pair_key = str(self.group_venn_multi_pair_selector.currentData() or "")
        if pair_key:
            keys = summary["comparison_A"].astype(str) + "||" + summary["comparison_B"].astype(str)
            matched = summary.loc[keys == pair_key]
            if not matched.empty:
                return matched.iloc[0]
        return summary.iloc[0]

    @staticmethod
    def _safe_column_label(value: str) -> str:
        label = re.sub(r"[^A-Za-z0-9]+", "_", str(value).strip()).strip("_")
        return label or "comparison"

    @staticmethod
    def _dpsi_sign_label(value: object) -> str:
        number = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
        if pd.isna(number):
            return "not_significant"
        if float(number) > 0:
            return "positive"
        if float(number) < 0:
            return "negative"
        return "zero"

    def _selected_venn_multi_comparisons(self) -> list:
        project = self.project_service.current_project
        if project is None:
            return []
        self._sync_group_venn_multi_comparison_checks()
        selected = list(self.project_service.selected_comparisons_for_display())
        checked_ids = set(self._checked_group_venn_multi_comparison_ids())
        if not checked_ids:
            return []
        chosen = []
        for comparison in selected:
            if comparison.comparison_id in checked_ids:
                chosen.append(comparison)
        return chosen

    def _build_group_venn_multi_tables(self) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
        comparisons = self._selected_venn_multi_comparisons()
        abs_cutoff = float(self.group_venn_multi_sig_cutoff.value())
        min_reads = int(self.group_venn_multi_min_reads.value())
        if not comparisons:
            return pd.DataFrame(), pd.DataFrame(), []
        event_frames: dict[str, pd.DataFrame] = {}
        display_names: list[str] = []
        for comparison in comparisons:
            display_name = comparison.display_resolved_name
            display_names.append(display_name)
            frame = self.project_service._cross_as_interaction_source_frame(
                comparison.comparison_id,
                abs_cutoff,
                min_reads,
                significant_only=True,
            )
            event_frames[comparison.comparison_id] = frame.copy()
        all_keys: list[str] = []
        seen_keys: set[str] = set()
        for frame in event_frames.values():
            if frame.empty or "event_matching_key" not in frame.columns:
                continue
            for key in frame["event_matching_key"].dropna().astype(str).tolist():
                if key not in seen_keys:
                    seen_keys.add(key)
                    all_keys.append(key)
        rows: list[dict[str, object]] = []
        for key in all_keys:
            present_names: list[str] = []
            positive_names: list[str] = []
            negative_names: list[str] = []
            base: dict[str, object] = {"event_matching_key": key}
            for comparison in comparisons:
                frame = event_frames.get(comparison.comparison_id, pd.DataFrame())
                display_name = comparison.display_resolved_name
                safe = self._safe_column_label(display_name)
                match = frame.loc[frame["event_matching_key"].astype(str) == key].head(1) if not frame.empty and "event_matching_key" in frame.columns else pd.DataFrame()
                present = not match.empty
                base[f"present_{safe}"] = bool(present)
                base[f"dPSI_{safe}"] = match.iloc[0].get("dPSI", pd.NA) if present else pd.NA
                base[f"FDR_{safe}"] = match.iloc[0].get("FDR", pd.NA) if present else pd.NA
                direction_label = self._dpsi_sign_label(base[f"dPSI_{safe}"]) if present else "not_significant"
                base[f"direction_{safe}"] = direction_label
                if present:
                    present_names.append(display_name)
                    if direction_label == "positive":
                        positive_names.append(display_name)
                    elif direction_label == "negative":
                        negative_names.append(display_name)
                    for column in ("gene", "gene_symbol", "gene_id", "event_type", "event_id", "coordinates"):
                        if column not in base or pd.isna(base.get(column)):
                            base[column] = match.iloc[0].get(column, pd.NA)
            base["n_comparisons_significant"] = len(present_names)
            base["comparison_membership"] = ";".join(present_names)
            base["membership_pattern"] = "shared_all" if len(present_names) == len(comparisons) else "comparison_specific" if len(present_names) == 1 else "shared_subset"
            direction_groups: list[str] = []
            if positive_names:
                direction_groups.append("positive:" + ";".join(positive_names))
            if negative_names:
                direction_groups.append("negative:" + ";".join(negative_names))
            base["direction_pattern"] = " | ".join(direction_groups) if direction_groups else "not_significant"
            base["has_opposite_direction"] = bool(positive_names and negative_names)
            rows.append(base)
        events = pd.DataFrame(rows)
        if not events.empty:
            events = events.sort_values(["n_comparisons_significant", "gene", "event_type"], ascending=[False, True, True], na_position="last").reset_index(drop=True)

        summary_rows: list[dict[str, object]] = []
        sets = {
            comparison.comparison_id: set(event_frames.get(comparison.comparison_id, pd.DataFrame()).get("event_matching_key", pd.Series(dtype="object")).dropna().astype(str).tolist())
            for comparison in comparisons
        }
        signed_sets: dict[str, set[str]] = {}
        direction_maps: dict[str, dict[str, str]] = {}
        for comparison in comparisons:
            frame = event_frames.get(comparison.comparison_id, pd.DataFrame())
            signed: set[str] = set()
            directions: dict[str, str] = {}
            if not frame.empty and "event_matching_key" in frame.columns:
                for _, row in frame.iterrows():
                    key = str(row.get("event_matching_key", ""))
                    if not key:
                        continue
                    direction = self._dpsi_sign_label(row.get("dPSI", pd.NA))
                    if direction not in {"positive", "negative"}:
                        continue
                    directions[key] = direction
                    signed.add(f"{key}::{direction}")
            signed_sets[comparison.comparison_id] = signed
            direction_maps[comparison.comparison_id] = directions
        for i, comp_a in enumerate(comparisons):
            set_a = sets.get(comp_a.comparison_id, set())
            signed_a = signed_sets.get(comp_a.comparison_id, set())
            direction_a = direction_maps.get(comp_a.comparison_id, {})
            for comp_b in comparisons[i + 1:]:
                set_b = sets.get(comp_b.comparison_id, set())
                signed_b = signed_sets.get(comp_b.comparison_id, set())
                direction_b = direction_maps.get(comp_b.comparison_id, {})
                identity_shared_keys = set_a & set_b
                identity_union = len(set_a | set_b)
                same_direction_shared = len(signed_a & signed_b)
                signed_union = len(signed_a | signed_b)
                opposite_direction_shared = sum(
                    1
                    for key in identity_shared_keys
                    if direction_a.get(key) in {"positive", "negative"}
                    and direction_b.get(key) in {"positive", "negative"}
                    and direction_a.get(key) != direction_b.get(key)
                )
                summary_rows.append(
                    {
                        "comparison_A": comp_a.comparison_id,
                        "comparison_B": comp_b.comparison_id,
                        "comparison_A_display_name": comp_a.display_resolved_name,
                        "comparison_B_display_name": comp_b.display_resolved_name,
                        "A_only": len(signed_a - signed_b),
                        "shared_A_and_B": same_direction_shared,
                        "B_only": len(signed_b - signed_a),
                        "same_direction_shared": same_direction_shared,
                        "opposite_direction_shared": opposite_direction_shared,
                        "identity_shared": len(identity_shared_keys),
                        "total_union": signed_union,
                        "identity_union": identity_union,
                        "jaccard": same_direction_shared / signed_union if signed_union else 0.0,
                        "identity_jaccard": len(identity_shared_keys) / identity_union if identity_union else 0.0,
                        "abs_dPSI_cutoff_used": abs_cutoff,
                        "min_reads_loaded": min_reads,
                    }
                )
        summary = pd.DataFrame(summary_rows)
        output_dir = self.project_service.cross_comparison_output_dir()
        if output_dir is not None:
            target_dir = output_dir / "selected_comparison_event_venn"
            target_dir.mkdir(parents=True, exist_ok=True)
            events.to_csv(target_dir / "selected_comparison_event_venn_events.tsv", sep="\t", index=False)
            summary.to_csv(target_dir / "selected_comparison_event_venn_pairwise_summary.tsv", sep="\t", index=False)
        return summary, events, display_names

    def _group_venn_multi_view_key(self) -> str:
        if self.fixed_subbranch == "venn_multi_similarity":
            return "similarity"
        return "venn"

    def _group_venn_multi_similarity_matrix(self, summary: pd.DataFrame, selected_names: list[str]) -> pd.DataFrame:
        matrix = pd.DataFrame(0.0, index=selected_names, columns=selected_names)
        for name in selected_names:
            matrix.loc[name, name] = 1.0
        if summary.empty:
            return matrix
        for _, row in summary.iterrows():
            name_a = str(row.get("comparison_A_display_name", ""))
            name_b = str(row.get("comparison_B_display_name", ""))
            value = float(pd.to_numeric(pd.Series([row.get("jaccard", 0.0)]), errors="coerce").fillna(0.0).iloc[0])
            if name_a in matrix.index and name_b in matrix.columns:
                matrix.loc[name_a, name_b] = value
                matrix.loc[name_b, name_a] = value
        return matrix

    def _draw_group_venn_multi_similarity(self, summary: pd.DataFrame, selected_names: list[str]) -> None:
        axis = self.figure.add_subplot(111)
        matrix = self._group_venn_multi_similarity_matrix(summary, selected_names)
        if len(selected_names) < 2:
            axis.axis("off")
            axis.text(
                0.5,
                0.5,
                "Need at least two selected comparisons for AS-event similarity.",
                ha="center",
                va="center",
            )
            self.canvas.draw_idle()
            return
        image = axis.imshow(matrix.values, vmin=0, vmax=1, cmap="YlGnBu", aspect="auto")
        axis.set_xticks(range(len(matrix.columns)))
        axis.set_xticklabels([str(x) for x in matrix.columns], rotation=90, ha="center", va="top", fontsize=8)
        axis.set_yticks(range(len(matrix.index)))
        axis.set_yticklabels([str(x) for x in matrix.index], fontsize=8)
        for y in range(matrix.shape[0]):
            for x in range(matrix.shape[1]):
                axis.text(x, y, f"{matrix.iat[y, x]:.2f}", ha="center", va="center", fontsize=8)
        self.figure.colorbar(image, ax=axis, fraction=0.036, pad=0.03, label="Direction-aware Jaccard")
        self._set_axis_title(axis, "5.4.6 AS-event response similarity heatmap")
        axis.set_xlabel("Selected comparisons")
        axis.set_ylabel("Selected comparisons")
        self.figure.suptitle("Similarity of strict significant AS-event responses", fontsize=15, fontweight="bold")
        self.figure.subplots_adjust(left=0.24, right=0.88, top=0.86, bottom=0.40)
        self.canvas.draw_idle()

    def _draw_group_venn_multi(self, summary: pd.DataFrame, selected_names: list[str], events: pd.DataFrame) -> None:
        self._scatter_collection = None
        self._scatter_frame = pd.DataFrame()
        self._pie_patches = []
        self._heatmap_image = None
        self._heatmap_rows = []
        self._heatmap_cols = []
        self._heatmap_values = None
        self._bar_patches = []
        self.figure.clear()
        if summary.empty or not selected_names:
            axis = self.figure.add_subplot(111)
            axis.axis("off")
            axis.text(
                0.5,
                0.5,
                "No cached 5.4.5 selected-comparison Venn result.\nClick Run Selected Venn to build strict AS-event overlap.",
                ha="center",
                va="center",
            )
            self.canvas.draw_idle()
            return
        if self._group_venn_multi_view_key() == "similarity":
            self._draw_group_venn_multi_similarity(summary, selected_names)
            return
        axis_venn = self.figure.add_subplot(111)
        direction_events = self._direction_aware_venn_events(events, selected_names)
        exact_counts = self._exact_venn_membership_counts(direction_events, selected_names)
        if len(selected_names) == 1:
            axis_venn.set_xlim(0, 1)
            axis_venn.set_ylim(0, 1)
            axis_venn.axis("off")
            axis_venn.add_patch(Circle((0.5, 0.5), 0.30, color="#2563EB", alpha=0.30, ec="#2563EB", lw=2))
            axis_venn.text(0.5, 0.5, f"{selected_names[0]}\n{sum(exact_counts.values())}", ha="center", va="center", fontsize=13, fontweight="bold")
            self._set_axis_title(axis_venn, "Direction-aware significant AS responses")
        elif len(selected_names) == 2:
            row = summary.iloc[0]
            self._draw_single_venn_axis(
                axis_venn,
                "Direction-aware exact Venn across selected comparisons",
                str(row.get("comparison_A_display_name", selected_names[0])),
                str(row.get("comparison_B_display_name", selected_names[1])),
                int(exact_counts.get((True, False), row.get("A_only", 0))),
                int(exact_counts.get((True, True), row.get("shared_A_and_B", 0))),
                int(exact_counts.get((False, True), row.get("B_only", 0))),
                "Shared requires same event identity and same dPSI sign.",
            )
        elif len(selected_names) == 3:
            self._draw_three_set_venn_axis(axis_venn, selected_names, exact_counts)
        elif len(selected_names) == 4:
            self._draw_four_set_venn_axis(axis_venn, selected_names, exact_counts)
        else:
            axis_venn.axis("off")
            top_patterns = sorted(exact_counts.items(), key=lambda item: item[1], reverse=True)[:14]
            axis_venn.text(0.5, 0.96, "Direction-aware exact intersections across selected comparisons", ha="center", va="top", fontsize=13, fontweight="bold")
            axis_venn.text(
                0.5,
                0.90,
                "More than 4 selected sets cannot be shown as a readable true Venn; signed event-response membership counts are shown below.",
                ha="center",
                va="top",
                fontsize=9,
                color="#4B5563",
            )
            for index, (mask, count) in enumerate(top_patterns):
                membership = " + ".join(name for name, present in zip(selected_names, mask, strict=False) if present)
                axis_venn.text(0.05, 0.80 - index * 0.052, textwrap.shorten(membership, width=58, placeholder="..."), ha="left", va="center", fontsize=8.5)
                axis_venn.text(0.95, 0.80 - index * 0.052, str(count), ha="right", va="center", fontsize=8.5, fontweight="bold")
        self.figure.suptitle("5.4.5 Direction-aware intersections across selected AS-event sets", fontsize=15, fontweight="bold")
        self.figure.subplots_adjust(left=0.05, right=0.97, top=0.86, bottom=0.12)
        self.canvas.draw_idle()

    def _draw_cross_as_interaction_summary(
        self,
        tables: dict[str, pd.DataFrame],
        pair_name: str | None,
        abs_cutoff: float,
        delta_cutoff: float,
        min_reads: int,
    ) -> None:
        self._scatter_collection = None
        self._scatter_frame = pd.DataFrame()
        self._pie_patches = []
        self._heatmap_image = None
        self._heatmap_rows = []
        self._heatmap_cols = []
        self._heatmap_values = None
        self._bar_patches = []
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        self._set_axis_title(axis, f"5.4.3 AS interaction modeling: {pair_name or 'selected pair'}")
        counts = {
            "Shared tested": len(tables.get("shared_interaction", pd.DataFrame())),
            "Opposite direction": len(tables.get("opposite_direction", pd.DataFrame())),
            "Same dir + large delta": len(tables.get("large_delta", pd.DataFrame())),
        }
        exclusive = tables.get("exclusive", pd.DataFrame()).copy()
        if not exclusive.empty and "exclusive_side" in exclusive.columns:
            side = exclusive["exclusive_side"].fillna("").astype(str)
            counts["A-specific only"] = int(side.eq("comparison_A_only").sum())
            counts["B-specific only"] = int(side.eq("comparison_B_only").sum())
        else:
            counts["A-specific only"] = 0
            counts["B-specific only"] = 0
        if not any(counts.values()):
            axis.axis("off")
            axis.text(0.5, 0.5, "No events found", ha="center", va="center")
            self.canvas.draw_idle()
            return
        bars = axis.bar(
            list(counts.keys()),
            list(counts.values()),
            color=["#0F766E", "#7C3AED", "#2563EB", "#2563EB", "#D97706"],
        )
        axis.set_ylabel("Events")
        supported_a = 0
        supported_b = 0
        if not exclusive.empty:
            support_side = exclusive.get("non_significant_side_read_support", pd.Series(dtype="object")).astype(str).str.startswith("supported_")
            side = exclusive.get("exclusive_side", pd.Series(dtype="object")).astype(str)
            supported_a = int((side.eq("comparison_A_only") & support_side).sum())
            supported_b = int((side.eq("comparison_B_only") & support_side).sum())
        self._bar_patches = [
            (
                patch,
                (
                    f"{label}\nEvents: {value}\n"
                    + (
                        f"Supported non-significant counterparts: {supported_a}\n"
                        if label == "A-specific only"
                        else f"Supported non-significant counterparts: {supported_b}\n"
                        if label == "B-specific only"
                        else ""
                    )
                    + f"abs(dPSI): {abs_cutoff:.3f}\nlarge delta: {delta_cutoff:.3f}\nmin reads: {min_reads}"
                ),
            )
            for patch, (label, value) in zip(bars, counts.items(), strict=False)
        ]
        self._finalize_axes_layout(rotate_x=True, bottom=0.28)
        self.canvas.draw_idle()

    def _draw_cross_as_interaction_batch_summary(
        self,
        tables: dict[str, pd.DataFrame],
        abs_cutoff: float,
        delta_cutoff: float,
        min_reads: int,
    ) -> None:
        self._scatter_collection = None
        self._scatter_frame = pd.DataFrame()
        self._pie_patches = []
        self._heatmap_image = None
        self._heatmap_rows = []
        self._heatmap_cols = []
        self._heatmap_values = None
        self._bar_patches = []
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        self._set_axis_title(axis, "5.4.3 all-pair AS interaction batch summary")
        summary = tables.get("batch_summary", pd.DataFrame()).copy()
        if summary.empty:
            axis.axis("off")
            axis.text(0.5, 0.5, "No all-pair summary available", ha="center", va="center")
            self.canvas.draw_idle()
            return
        labels = summary.get("pair_name", pd.Series(dtype="object")).fillna(summary.get("pair_id", "")).astype(str).tolist()
        x = np.arange(len(summary))
        width = 0.18
        high = pd.to_numeric(summary.get("n_high_confidence_interaction"), errors="coerce").fillna(0).to_numpy()
        opposite = pd.to_numeric(summary.get("n_opposite_direction"), errors="coerce").fillna(0).to_numpy()
        exclusive_a = pd.to_numeric(summary.get("n_exclusive_A_supported"), errors="coerce").fillna(0).to_numpy()
        exclusive_b = pd.to_numeric(summary.get("n_exclusive_B_supported"), errors="coerce").fillna(0).to_numpy()
        bars_high = axis.bar(x - (1.5 * width), high, width, label="High-confidence interaction", color="#0F766E")
        bars_opposite = axis.bar(x - (0.5 * width), opposite, width, label="Opposite direction", color="#7C3AED")
        bars_exclusive_a = axis.bar(x + (0.5 * width), exclusive_a, width, label="A-specific supported", color="#2563EB")
        bars_exclusive_b = axis.bar(x + (1.5 * width), exclusive_b, width, label="B-specific supported", color="#D97706")
        axis.set_ylabel("Events")
        axis.set_xticks(x)
        axis.set_xticklabels(labels)
        axis.legend(loc="upper right", frameon=False)
        self._bar_patches = []
        for patches, metric, values in (
            (bars_high, "High-confidence interaction", high),
            (bars_opposite, "Opposite direction", opposite),
            (bars_exclusive_a, "A-specific supported", exclusive_a),
            (bars_exclusive_b, "B-specific supported", exclusive_b),
        ):
            self._bar_patches.extend(
                (
                    patch,
                    f"{label}\n{metric}: {int(value)}\nabs(dPSI): {abs_cutoff:.3f}\nlarge delta: {delta_cutoff:.3f}\nmin reads: {min_reads}",
                )
                for patch, label, value in zip(patches, labels, values, strict=False)
            )
        self._finalize_axes_layout(rotate_x=True, bottom=0.32)
        self.canvas.draw_idle()

    def _render_tx_splicing_branch(self) -> None:
        comparison_id = self._current_tab_key(self.tx_tabs)
        comparison_name = self._comparison_name_from_id(comparison_id)
        self.context_header.setText(self.project_service.comparison_context_text(comparison_id))
        mode = self.fixed_subbranch or self._navigation_subbranch or ""
        event_type_map = {"se": "SE", "ri": "RI", "a3ss": "A3SS", "a5ss": "A5SS", "mxe": "MXE"}
        if self.project_service.run_state.results.tx_splicing_gene_table.empty:
            self._set_explanation("Transcript + splicing integration has not been run yet. Click 'Run Current Branch' to compute this analysis for the current project.")
            self._draw_placeholder("Not run.\nNo cached transcript + splicing integration result is loaded.")
            self._set_detail_text(
                module_detail("5.2.integration"),
                "Not run. Open Detail for 5.2 purpose, inputs and direction rules.",
            )
            self._sync_footer_info()
            self._reset_plot_view()
            return
        integration, _ = self.project_service.preview_tx_splicing_integration(allow_generate=False)
        filtered = self._filter_frame_by_value(integration, "comparison_id", comparison_id)
        if mode in event_type_map:
            event_type = event_type_map[mode]
            if "representative_event_type" in filtered.columns:
                subset = filtered.loc[
                    filtered["representative_event_type"].astype(str) == event_type
                ].copy()
            else:
                subset = pd.DataFrame(columns=filtered.columns)
            subset = self._integration_plot_frame(subset)
            self._draw_event_type_subanalysis(subset, comparison_id, comparison_name, event_type)
            self._populate_table(
                subset,
                [
                    "comparison_id",
                    "comparison_name",
                    "gene_symbol",
                    "gene_id",
                    "representative_event_type",
                    "representative_event_id",
                    "representative_event_FDR",
                    "standardized_log2FC",
                    "standardized_dPSI",
                    "DE_FDR",
                    "rMATS_FDR",
                    "log2fc_direction_label",
                    "dpsi_direction_label",
                    "transcript_direction_flipped",
                    "splicing_direction_flipped",
                    "plot_class",
                    "direction_class",
                ],
                preview_name=f"{event_type} event-type subanalysis",
            )
            return

        filtered = self._integration_plot_frame(filtered)
        self._draw_integration_scatter(filtered, comparison_name)
        self._populate_table(
            filtered,
            [
                "comparison_id",
                "comparison_name",
                "gene_symbol",
                "gene_id",
                "DE_FDR",
                "standardized_log2FC",
                "rMATS_FDR",
                "standardized_dPSI",
                "log2fc_direction_label",
                "dpsi_direction_label",
                "transcript_direction_flipped",
                "splicing_direction_flipped",
                "plot_class",
                "agreement_class",
                "n_sig_events_all",
            ],
        )
        self._set_detail_text(
            module_detail("5.2.integration"),
            "5.2.1 Transcript + Splicing Integration loaded. Open Detail for DEG/DAS class and direction definitions.",
        )

    def _render_mechanism_branch(self) -> None:
        if self.project_service.run_state.results.mechanism_support.empty:
            self.context_header.clear()
            self._set_explanation("Mechanism support has not been run yet. Click 'Run Current Branch' to compute the support layer.")
            self._draw_placeholder("Not run.\nNo cached mechanism-support result is loaded.")
            self._sync_footer_info()
            self._reset_plot_view()
            return
        mechanism = self.project_service.preview_mechanism_support(allow_generate=False)
        comparison_id = self._current_tab_key(self.mechanism_compare_tabs)
        if comparison_id and mechanism is not None and not mechanism.empty and "comparison_id" in mechanism.columns:
            mechanism = mechanism.loc[mechanism["comparison_id"].astype(str) == comparison_id].copy()
        self.context_header.setText(self.project_service.comparison_header_text(comparison_id))
        self._draw_mechanism_support(mechanism)

    def _render_candidate_branch(self) -> None:
        comparison_id = self._current_tab_key(self.candidate_tabs)
        mode = self.fixed_subbranch or self._navigation_subbranch or self._current_candidate_subbranch()
        source = str(self.candidate_source_combo.currentData() or "selection")
        top_n = int(self.candidate_top_n.value())
        custom_genes = self._candidate_custom_gene_tokens()
        self.context_header.setText(self.project_service.comparison_header_text(comparison_id))
        title = f"Candidate gene screening: {self._comparison_title_label(comparison_id)}"
        if self.project_service.run_state.results.candidate_gene_table.empty:
            self._set_explanation("Candidate screening has not been generated yet. Open 5.3.1 or 5.3.2 and explicitly generate/load candidate results first.")
            self._draw_placeholder("Not run.\nNo cached candidate-screening result is loaded.")
            self._set_detail_text(
                module_detail("5.3.ranking"),
                "Not run. Open Detail for 5.3 candidate-ranking inputs and cache behavior.",
            )
            self._sync_footer_info()
            self._reset_plot_view()
            return

        if mode == "isoform":
            frame, meta = self.project_service.candidate_isoform_followup_frame(
                comparison_id,
                source=source,
                top_n=top_n,
                custom_genes=custom_genes,
                allow_generate=False,
            )
            self._draw_isoform_followup(frame, title, meta)
            return

        frame, meta = self.project_service.candidate_selection_frame(
            comparison_id,
            source=source,
            top_n=top_n,
            custom_genes=custom_genes,
            allow_generate=False,
        )
        if mode == "table":
            self._draw_text_panel(title)
            self._populate_table(
                frame,
                [
                    "comparison_id",
                    "rank",
                    "gene_symbol",
                    "gene_id",
                    "candidate_tier",
                    "evidence_class",
                    "evidence_count",
                    "DEG_significant",
                    "rMATS_significant",
                    "DEXSeq_significant",
                    "DTU_significant",
                    "DE_FDR",
                    "standardized_log2FC",
                    "best_rMATS_FDR",
                    "dominant_rMATS_standardized_dPSI",
                    "direction_class",
                    "shortlist_gene",
                    "blacklist_gene",
                    "candidate_reason",
                ],
                file_path=self._candidate_file_map(comparison_id).get("gene_level_integrated_candidates"),
                output_dir=self.project_service.candidate_screening_output_dir(comparison_id),
                preview_name="Candidate screening table",
            )
            self._set_detail_text(
                module_detail("5.3.ranking"),
                "5.3.1 Candidate ranking loaded. Open Detail for rank, evidence and uniqueness definitions.",
            )
        elif mode == "heatmap":
            self._draw_candidate_evidence_heatmap(frame, title, meta)
        elif mode == "program":
            self._draw_candidate_program_heatmap(frame, title, meta)
        else:
            payload, card_meta = self.project_service.candidate_gene_card_payload(
                comparison_id,
                source=source,
                top_n=top_n,
                custom_genes=custom_genes,
                allow_generate=False,
                selected_gene_key=self._candidate_card_selected_gene_key,
            )
            self._sync_candidate_card_gene_selector(payload)
            self._draw_candidate_cards(payload, title, card_meta)

    def _effective_branch_keys(self) -> tuple[str | None, str | None]:
        if self.fixed_branch is not None:
            return self.fixed_branch, self.fixed_subbranch
        branch = self._navigation_branch
        if branch is not None:
            return branch, self._navigation_subbranch
        branch_map = {
            0: "landscape",
            1: "group",
            2: "tx",
            3: "mechanism",
            4: "candidate",
        }
        return branch_map.get(self.tabs.currentIndex()), None

    def _current_landscape_subbranch(self) -> str:
        return {0: "overview", 1: "bar", 2: "pie"}.get(self.landscape_tabs.currentIndex(), "overview")

    def _current_group_subbranch(self) -> str:
        return {0: "excel", 1: "heatmap", 2: "direction"}.get(self.group_mode_tabs.currentIndex(), "excel")

    def _current_candidate_subbranch(self) -> str:
        return {0: "table", 1: "heatmap", 2: "program", 3: "cards", 4: "isoform"}.get(self.candidate_mode_tabs.currentIndex(), "table")

    def _current_landscape_comparison_id(self) -> str | None:
        return self._current_tab_key(self.landscape_compare_tabs)

    def _draw_landscape_overview(self, splicing_landscape: pd.DataFrame | None) -> None:
        self._scatter_collection = None
        self._scatter_frame = pd.DataFrame()
        self._bar_patches = []
        self._pie_patches = []
        self._heatmap_image = None
        self._heatmap_rows = []
        self._heatmap_cols = []
        self._heatmap_values = None
        self._set_explanation(
            "Landscape overview compares all selected pairings using stable AS events after rMATS significance, read-support and replicate-consistency filters."
        )
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        self._set_axis_title(axis, "Landscape: significant event composition across selected pairings")
        self._hover_annotation = axis.annotate(
            "",
            xy=(0, 0),
            xytext=(12, 12),
            textcoords="offset points",
            bbox={"boxstyle": "round", "fc": "white", "ec": "#666666", "alpha": 0.9},
            arrowprops={"arrowstyle": "->", "color": "#666666"},
        )
        self._hover_annotation.set_visible(False)
        if splicing_landscape is None or splicing_landscape.empty:
            axis.text(0.5, 0.5, "Landscape results are not available yet.", ha="center", va="center")
            self._populate_table(pd.DataFrame(), [])
            self._set_detail_text(
                module_detail("5.1.landscape"),
                "Not run. Open Detail for Landscape purpose, inputs and cache behavior.",
            )
            self._finalize_axes_layout(bottom=0.12)
            self.canvas.draw_idle()
            return

        overview = splicing_landscape.copy()
        count_column = "n_stable" if "n_stable" in overview.columns else "n_sig"
        overview[count_column] = pd.to_numeric(overview.get(count_column), errors="coerce").fillna(0)
        overview = overview.loc[overview[count_column] > 0].copy()
        if overview.empty:
            axis.text(0.5, 0.5, "No stable AS events available for the selected pairings.", ha="center", va="center")
            self._populate_table(pd.DataFrame(), [])
            self._finalize_axes_layout(bottom=0.12)
            self.canvas.draw_idle()
            return
        comparison_order = self._selected_comparison_display_order()
        pivot = (
            overview.pivot_table(
                index="comparison_name",
                columns="event_type",
                values=count_column,
                aggfunc="sum",
                fill_value=0,
            )
            .reindex(index=[label for label in comparison_order if label in overview["comparison_name"].astype(str).tolist()], fill_value=0)
            .reindex(columns=["SE", "RI", "A3SS", "A5SS", "MXE"], fill_value=0)
        )
        colors = self._event_color_map()
        bottom = np.zeros(len(pivot.index))
        for event_type in pivot.columns:
            values = pivot[event_type].to_numpy(dtype=float)
            bars = axis.bar(
                pivot.index.tolist(),
                values,
                bottom=bottom,
                label=event_type,
                color=colors.get(event_type, "#6C8AE4"),
                edgecolor="white",
                linewidth=0.6,
            )
            self._bar_patches.extend(
                (patch, f"{comparison}\n{event_type}: {int(value)} stable AS events")
                for patch, comparison, value in zip(bars, pivot.index.tolist(), values, strict=False)
            )
            bottom = bottom + values
        axis.set_ylabel("Stable AS events", fontweight="bold")
        axis.legend(
            loc="upper left",
            bbox_to_anchor=(1.02, 1.00),
            title="Event type",
            frameon=False,
            borderaxespad=0.0,
        )

        preview = overview[
            [
                column
                for column in [
                    "comparison_id",
                    "comparison_name",
                    "event_type",
                    "n_total",
                    "n_sig",
                    "n_stable",
                    "stable_fraction_of_sig",
                    "stable_min_junction_reads",
                    "stable_min_replicates_per_group",
                    "stable_min_pairwise_consistency",
                ]
                if column in overview.columns
            ]
        ].copy()
        self._populate_table(
            preview,
            list(preview.columns),
            preview_name="Landscape overview summary",
        )
        axis.set_xlabel("")
        self._finalize_axes_layout(rotate_x=True, bottom=0.28, left=0.12, right=0.76, top=0.82)
        self._set_detail_text(
            module_detail("5.1.landscape"),
            "5.1.1 Landscape loaded. Open Detail for purpose, inputs, outputs and interpretation.",
        )
        self.canvas.draw_idle()

    def _draw_landscape_bar(self, splicing_landscape: pd.DataFrame | None) -> None:
        self._scatter_collection = None
        self._scatter_frame = pd.DataFrame()
        self._bar_patches = []
        self._pie_patches = []
        self._heatmap_image = None
        self._heatmap_rows = []
        self._heatmap_cols = []
        self._heatmap_values = None
        self._set_explanation(
            "Landscape bar plot compares stable AS-event counts by event type across all selected pairings."
        )
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        self._set_axis_title(axis, "Landscape bar: significant events by type across selected pairings")
        self._hover_annotation = axis.annotate(
            "",
            xy=(0, 0),
            xytext=(12, 12),
            textcoords="offset points",
            bbox={"boxstyle": "round", "fc": "white", "ec": "#666666", "alpha": 0.9},
            arrowprops={"arrowstyle": "->", "color": "#666666"},
        )
        self._hover_annotation.set_visible(False)
        if splicing_landscape is None or splicing_landscape.empty:
            axis.text(0.5, 0.5, "Landscape results are not available yet.", ha="center", va="center")
            self._populate_table(pd.DataFrame(), [])
            self._set_detail_text(
                module_detail("5.1.bar"),
                "Not run. Open Detail for Bar purpose, inputs and cache behavior.",
            )
            self._finalize_axes_layout(bottom=0.12)
            self.canvas.draw_idle()
            return
        working = splicing_landscape.copy()
        count_column = "n_stable" if "n_stable" in working.columns else "n_sig"
        working[count_column] = pd.to_numeric(working[count_column], errors="coerce").fillna(0)
        working = working.loc[working[count_column] > 0].copy()
        if working.empty:
            axis.text(0.5, 0.5, "No stable AS events available for the selected pairings.", ha="center", va="center")
            self._populate_table(pd.DataFrame(), [])
            self._finalize_axes_layout(bottom=0.12)
            self.canvas.draw_idle()
            return
        comparison_order = self._selected_comparison_display_order()
        pivot = (
            working.pivot_table(
                index="event_type",
                columns="comparison_name",
                values=count_column,
                aggfunc="sum",
                fill_value=0,
            )
            .reindex(index=["SE", "RI", "A3SS", "A5SS", "MXE"], fill_value=0)
        )
        pivot = pivot.reindex(columns=[label for label in comparison_order if label in pivot.columns], fill_value=0)
        x = np.arange(len(pivot.index))
        comparison_labels = list(pivot.columns)
        width = 0.14 if comparison_labels else 0.5
        colors = ["#4C78A8", "#F58518", "#54A24B", "#E45756", "#9D6BCB", "#72B7B2"]
        for index, comparison_label in enumerate(comparison_labels):
            values = pivot[comparison_label].to_numpy(dtype=float)
            offsets = x + (index - (len(comparison_labels) - 1) / 2) * width
            bars = axis.bar(
                offsets,
                values,
                width=width,
                label=comparison_label,
                color=colors[index % len(colors)],
                edgecolor="white",
                linewidth=0.6,
            )
            self._bar_patches.extend(
                (patch, f"{comparison_label}\n{event_type}: {int(value)} stable AS events")
                for patch, event_type, value in zip(bars, pivot.index.tolist(), values, strict=False)
            )
        axis.set_xticks(x)
        axis.set_xticklabels(pivot.index.tolist(), fontweight="bold")
        axis.set_ylabel("Stable AS events", fontweight="bold")
        axis.set_xlabel("Event type", fontweight="bold")
        if comparison_labels:
            axis.legend(
                loc="upper left",
                bbox_to_anchor=(1.02, 1.00),
                title="Comparison",
                frameon=False,
                borderaxespad=0.0,
            )
        self._populate_table(
            working,
            [
                column
                for column in [
                    "comparison_id",
                    "comparison_name",
                    "event_type",
                    "n_stable",
                    "n_sig",
                    "n_total",
                    "stable_fraction_of_sig",
                    "stable_min_junction_reads",
                    "stable_min_replicates_per_group",
                    "stable_min_pairwise_consistency",
                ]
                if column in working.columns
            ],
            preview_name="Landscape bar summary",
        )
        self._finalize_axes_layout(rotate_x=True, bottom=0.24, left=0.12, right=0.74, top=0.82)
        self._set_detail_text(
            module_detail("5.1.bar"),
            "5.1.2 Bar loaded. Open Detail for purpose, inputs, outputs and interpretation.",
        )
        self.canvas.draw_idle()

    def _draw_landscape_pie(
        self,
        splicing_landscape: pd.DataFrame | None,
        comparison_id: str | None,
        comparison_name: str | None,
    ) -> None:
        self._scatter_collection = None
        self._scatter_frame = pd.DataFrame()
        self._bar_patches = []
        self._pie_patches = []
        self._heatmap_image = None
        self._heatmap_rows = []
        self._heatmap_cols = []
        self._heatmap_values = None
        self._set_explanation("Landscape pie chart: event-type composition of stable AS events for the selected comparison.")
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        self._set_axis_title(axis, f"Landscape pie: {self._comparison_title_label(comparison_id)}")
        self._hover_annotation = axis.annotate(
            "",
            xy=(0, 0),
            xytext=(12, 12),
            textcoords="offset points",
            bbox={"boxstyle": "round", "fc": "white", "ec": "#666666", "alpha": 0.9},
            arrowprops={"arrowstyle": "->", "color": "#666666"},
        )
        self._hover_annotation.set_visible(False)
        if splicing_landscape is None or splicing_landscape.empty:
            axis.text(0.5, 0.5, "Landscape results are not available yet.", ha="center", va="center")
            self._populate_table(pd.DataFrame(), [])
            self._set_detail_text(
                module_detail("5.1.pie"),
                "Not run. Open Detail for Pie purpose, inputs and cache behavior.",
            )
            self._finalize_axes_layout(bottom=0.12)
            self.canvas.draw_idle()
            return
        event_order = ["SE", "RI", "A3SS", "A5SS", "MXE"]
        count_column = "n_stable" if "n_stable" in splicing_landscape.columns else "n_sig"
        grouped = (
            splicing_landscape.assign(**{count_column: pd.to_numeric(splicing_landscape[count_column], errors="coerce").fillna(0)})
            .groupby("event_type", dropna=False)[count_column]
            .sum()
            .reindex(event_order, fill_value=0)
        )
        grouped = grouped.loc[grouped > 0]
        if grouped.empty or grouped.sum() == 0:
            axis.text(0.5, 0.5, "No non-zero event counts available.", ha="center", va="center")
            self._populate_table(pd.DataFrame(), [])
            self.details.setPlainText("Pie chart needs non-zero event totals.")
            self._finalize_axes_layout(bottom=0.12)
            self.canvas.draw_idle()
            return
        event_colors = self._event_color_map()
        colors = [event_colors.get(str(label), "#9CA3AF") for label in grouped.index.tolist()]
        wedges, _texts, autotexts = axis.pie(
            grouped.values.tolist(),
            labels=grouped.index.tolist(),
            autopct=lambda pct: f"{pct:.1f}%",
            startangle=90,
            radius=0.72,
            pctdistance=0.68,
            labeldistance=1.02,
            colors=colors,
            textprops={"fontweight": "bold", "fontsize": 11},
        )
        for autotext in autotexts:
            autotext.set_color("white")
            autotext.set_fontweight("bold")
        axis.legend(
            wedges,
            grouped.index.tolist(),
            title="Event type",
            loc="upper left",
            bbox_to_anchor=(1.02, 1.00),
            borderaxespad=0.0,
            frameon=False,
        )
        self._pie_patches = [
            (wedge, f"{label}\nStable AS events: {value}\nShare: {value / max(grouped.sum(), 1) * 100:.1f}%")
            for wedge, label, value in zip(wedges, grouped.index.tolist(), grouped.values.tolist(), strict=False)
        ]
        axis.axis("equal")
        self._populate_table(
            grouped.reset_index().rename(columns={"event_type": "Event Type", count_column: "Stable AS Events"}),
            ["Event Type", "Stable AS Events"],
            preview_name="Landscape pie summary",
        )
        self._finalize_axes_layout(bottom=0.10, left=0.12, right=0.76, top=0.82)
        self._set_detail_text(
            module_detail("5.1.pie"),
            "5.1.3 Pie loaded. Open Detail for purpose, inputs, outputs and interpretation.",
        )
        self.canvas.draw_idle()

    def _draw_program_summary(self, program_summary: pd.DataFrame | None, pair_name: str | None) -> None:
        self._set_explanation(
            "Group comparison summary: this branch compares comparison A against comparison B. Use it to see whether splicing programs are shared, opposite, A-specific, or B-specific between the two comparisons."
        )
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        axis.set_title(f"Group comparison: {pair_name or 'No comparison set'}")
        self._hover_annotation = axis.annotate(
            "",
            xy=(0, 0),
            xytext=(12, 12),
            textcoords="offset points",
            bbox={"boxstyle": "round", "fc": "white", "ec": "#666666", "alpha": 0.9},
            arrowprops={"arrowstyle": "->", "color": "#666666"},
        )
        self._hover_annotation.set_visible(False)
        self._bar_patches = []
        if program_summary is None or program_summary.empty:
            axis.text(0.5, 0.5, "No group-comparison summary available.", ha="center", va="center")
            self.canvas.draw_idle()
            return
        grouped = program_summary.groupby("class_label", dropna=False)["n_events"].sum().sort_values(ascending=False)
        bars = axis.bar(grouped.index.tolist(), grouped.values.tolist(), color=["#264653", "#2A9D8F", "#E9C46A", "#E76F51"])
        axis.tick_params(axis="x", rotation=20)
        axis.set_ylabel("Events")
        self._bar_patches = [
            (patch, f"{label}\nEvents: {value}")
            for patch, label, value in zip(bars, grouped.index.tolist(), grouped.values.tolist(), strict=False)
        ]
        self._finalize_axes_layout(rotate_x=True, bottom=0.30)
        self.canvas.draw_idle()

    def _draw_program_heatmap(self, program_events: pd.DataFrame | None, pair_name: str | None) -> None:
        self._scatter_collection = None
        self._scatter_frame = pd.DataFrame()
        self._bar_patches = []
        self._pie_patches = []
        self._set_explanation(
            "Program heatmap: rows are genes, columns are the two compared conditions. Color encodes dPSI, so opposite colors indicate direction differences and similar colors with different intensity indicate strength differences."
        )
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        axis.set_title(f"Program heatmap: {pair_name or 'No comparison set'}")
        if program_events is None or program_events.empty or "gene_symbol" not in program_events.columns:
            axis.text(0.5, 0.5, "No program events available for heatmap.", ha="center", va="center")
            self._populate_table(pd.DataFrame(), [])
            self.details.setPlainText("Run Group Comparison first to generate program heatmap inputs.")
            self.canvas.draw_idle()
            return
        frame = program_events.copy()
        frame["gene_label"] = frame["gene_symbol"].fillna(frame.get("gene_id")).astype(str)
        frame = frame[frame["gene_label"].ne("")].copy()
        if frame.empty:
            axis.text(0.5, 0.5, "No gene labels available for heatmap.", ha="center", va="center")
            self.canvas.draw_idle()
            return
        grouped = (
            frame.groupby("gene_label", dropna=False)["abs_delta_between_programs"]
            .max()
            .sort_values(ascending=False)
            .head(30)
        )
        top_genes = grouped.index.tolist()
        subset = frame[frame["gene_label"].isin(top_genes)].drop_duplicates("gene_label")
        comparison_a_label = "Comparison A"
        comparison_b_label = "Comparison B"
        if "comparison_A_name" in subset.columns and not subset.empty:
            comparison_a_label = f"A: {subset['comparison_A_name'].iloc[0]}"
        if "comparison_B_name" in subset.columns and not subset.empty:
            comparison_b_label = f"B: {subset['comparison_B_name'].iloc[0]}"
        heatmap = pd.DataFrame(
            {
                comparison_a_label: subset.set_index("gene_label")["dPSI_A"],
                comparison_b_label: subset.set_index("gene_label")["dPSI_B"],
            }
        ).reindex(top_genes)
        image = axis.imshow(heatmap.fillna(0).values, cmap="coolwarm", aspect="auto", vmin=-1, vmax=1)
        self._heatmap_image = image
        self._heatmap_rows = heatmap.index.astype(str).tolist()
        self._heatmap_cols = heatmap.columns.astype(str).tolist()
        self._heatmap_values = heatmap.fillna(0).values
        axis.set_yticks(range(len(heatmap.index)))
        axis.set_yticklabels(heatmap.index)
        axis.set_xticks(range(len(heatmap.columns)))
        axis.set_xticklabels(heatmap.columns, rotation=20)
        self.figure.colorbar(image, ax=axis, fraction=0.03, pad=0.02, label="dPSI")
        preview = subset[["gene_label", "class_label", "dPSI_A", "dPSI_B", "abs_delta_between_programs"]].copy()
        self._populate_table(preview, list(preview.columns), output_dir=self.project_service.cross_comparison_output_dir(), preview_name="Program heatmap source genes")
        self._finalize_axes_layout(rotate_x=True, bottom=0.28)
        self.details.setPlainText(
            "Heatmap of top genes with the largest program-level dPSI differences.\n"
            "Use this to spot direction differences and same-direction large-strength shifts."
        )
        self.canvas.draw_idle()

    def _draw_direction_strength_summary(self, program_events: pd.DataFrame | None, pair_name: str | None) -> None:
        self._scatter_collection = None
        self._scatter_frame = pd.DataFrame()
        self._pie_patches = []
        self._heatmap_image = None
        self._heatmap_rows = []
        self._heatmap_cols = []
        self._heatmap_values = None
        self._set_explanation(
            "Direction vs strength summary: opposite_direction means sign reversal between comparisons; same_direction_large_delta means both are significant in the same direction but differ strongly in magnitude."
        )
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        axis.set_title(f"Direction vs strength summary: {pair_name or 'No comparison set'}")
        if program_events is None or program_events.empty:
            axis.text(0.5, 0.5, "No program events available.", ha="center", va="center")
            self._populate_table(pd.DataFrame(), [])
            self.details.setPlainText("Run Group Comparison first to classify direction and strength differences.")
            self.canvas.draw_idle()
            return
        filtered = program_events.loc[
            program_events["class"].astype(str).isin(
                ["opposite_direction", "same_direction_large_delta", "shared_same_direction", "A_only", "B_only"]
            )
        ].copy()
        if filtered.empty:
            axis.text(0.5, 0.5, "No significant direction/strength differences found.", ha="center", va="center")
            self._populate_table(pd.DataFrame(), [])
            self.details.setPlainText("No significant direction or strength differences were detected for the current comparison pair.")
            self.canvas.draw_idle()
            return
        summary = (
            filtered.groupby("class_label", dropna=False)["event_uid"]
            .count()
            .sort_values(ascending=False)
        )
        color_map = {
            "Opposite direction": "#C1121F",
            "Same direction, large dPSI gap": "#F4A261",
            "Shared same direction": "#2A9D8F",
        }
        bar_colors = [color_map.get(label, "#6C8AE4") for label in summary.index.tolist()]
        bars = axis.bar(summary.index.tolist(), summary.values.tolist(), color=bar_colors)
        axis.tick_params(axis="x", rotation=20)
        axis.set_ylabel("Events")
        self._bar_patches = [
            (patch, f"{label}\nEvents: {value}")
            for patch, label, value in zip(bars, summary.index.tolist(), summary.values.tolist(), strict=False)
        ]
        preview = filtered[
            [
                "gene_symbol",
                "event_type",
                "class_label",
                "dPSI_A",
                "dPSI_B",
                "abs_delta_between_programs",
            ]
        ].copy()
        self._populate_table(preview, list(preview.columns), output_dir=self.project_service.cross_comparison_output_dir(), preview_name="Direction vs strength summary")
        self._finalize_axes_layout(rotate_x=True, bottom=0.30)
        self.details.setPlainText(
            "Summary split into direction differences and strength differences.\n"
            "Opposite direction = condition-dependent splicing reversal.\n"
            "Same direction, large dPSI gap = same sign but strong magnitude shift."
        )
        self.canvas.draw_idle()

    def _integration_plot_frame(self, integration: pd.DataFrame | None) -> pd.DataFrame:
        if integration is None or integration.empty:
            return pd.DataFrame()
        frame = integration.copy()
        frame["standardized_log2FC"] = pd.to_numeric(frame.get("standardized_log2FC"), errors="coerce")
        frame["standardized_dPSI"] = pd.to_numeric(frame.get("standardized_dPSI"), errors="coerce")
        frame["DE_FDR"] = pd.to_numeric(frame.get("DE_FDR"), errors="coerce")
        frame["rMATS_FDR"] = pd.to_numeric(frame.get("rMATS_FDR"), errors="coerce")
        deg_series = frame.get("deg_sig")
        if deg_series is None:
            deg_series = frame.get("DEG_significant")
        if deg_series is None:
            deg_series = pd.Series(False, index=frame.index)
        frame["deg_sig"] = deg_series.map(self._safe_bool)
        splicing_series = frame.get("splicing_sig")
        if splicing_series is None:
            splicing_series = frame.get("rMATS_significant")
        if splicing_series is None:
            splicing_series = pd.Series(False, index=frame.index)
        frame["splicing_sig"] = splicing_series.map(self._safe_bool)
        frame = frame.dropna(subset=["standardized_log2FC", "standardized_dPSI"]).reset_index(drop=True)
        if frame.empty:
            return frame
        frame["plot_class"] = frame.apply(self._integration_plot_class, axis=1)
        return frame.loc[frame["plot_class"].astype(str) != "Not significant"].copy()

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
    def _safe_text(value, default: str = "") -> str:
        if value is None:
            return default
        try:
            if pd.isna(value):
                return default
        except (TypeError, ValueError):
            pass
        text = str(value).strip()
        return text if text else default

    @classmethod
    def _first_text(cls, *values: object, default: str = "") -> str:
        for value in values:
            text = cls._safe_text(value)
            if text:
                return text
        return default

    @staticmethod
    def _integration_plot_class(row: pd.Series) -> str:
        deg_sig = ResultsPage._safe_bool(row.get("deg_sig"))
        splicing_sig = ResultsPage._safe_bool(row.get("splicing_sig"))
        agreement = str(row.get("agreement_class", ""))
        if deg_sig and splicing_sig:
            return "DEG + DAS same direction" if agreement == "concordant" else "DEG + DAS opposite direction"
        if deg_sig:
            return "DEG only"
        if splicing_sig:
            return "DAS only"
        return "Not significant"

    def _draw_integration_scatter(self, integration: pd.DataFrame | None, comparison_name: str | None, *, event_type: str | None = None) -> None:
        self._bar_patches = []
        self._pie_patches = []
        self._heatmap_image = None
        self._heatmap_rows = []
        self._heatmap_cols = []
        self._heatmap_values = None
        self._set_explanation(
            "Transcript + splicing scatter: x = standardized DEG log2FC, y = standardized representative dPSI. Colors distinguish DEG-only, DAS-only, both significant same direction, and both significant opposite direction."
        )
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        subtitle = f" ({event_type})" if event_type else ""
        self._set_axis_title(axis, f"Transcript + splicing: {comparison_name or 'selected comparison'}{subtitle}")
        axis.set_xlabel("standardized log2FC", fontweight="bold")
        axis.set_ylabel("standardized dPSI", fontweight="bold")
        axis.axhline(0, color="#888888", linewidth=1)
        axis.axvline(0, color="#888888", linewidth=1)
        self._hover_annotation = axis.annotate(
            "",
            xy=(0, 0),
            xytext=(12, 12),
            textcoords="offset points",
            bbox={"boxstyle": "round", "fc": "white", "ec": "#666666", "alpha": 0.9},
            arrowprops={"arrowstyle": "->", "color": "#666666"},
        )
        self._hover_annotation.set_visible(False)
        self._bar_patches = []
        frame = self._integration_plot_frame(integration)
        if frame.empty:
            axis.axis("off")
            axis.text(0.5, 0.5, "No transcript/splicing rows available.", ha="center", va="center", transform=axis.transAxes)
            self._scatter_collection = None
            self._scatter_frame = pd.DataFrame()
            self._finalize_axes_layout(bottom=0.10)
            self.canvas.draw_idle()
            return
        color_map = {
            "DEG only": "#577590",
            "DAS only": "#43AA8B",
            "DEG + DAS same direction": "#D1495B",
            "DEG + DAS opposite direction": "#8B5CF6",
        }
        colors = [color_map.get(str(value), "#9CA3AF") for value in frame.get("plot_class", pd.Series(dtype=str))]
        self._scatter_collection = axis.scatter(
            frame["standardized_log2FC"],
            frame["standardized_dPSI"],
            c=colors,
            s=42,
            alpha=0.85,
            edgecolors="white",
            linewidths=0.4,
        )
        legend_handles = []
        legend_labels = []
        for class_name, color in color_map.items():
            if any(str(value) == class_name for value in frame.get("plot_class", pd.Series(dtype=str))):
                legend_handles.append(axis.scatter([], [], c=color, s=36, alpha=0.85))
                legend_labels.append(class_name)
        if legend_handles:
            axis.legend(
                legend_handles,
                legend_labels,
                title="DEG + DAS class",
                loc="upper left",
                bbox_to_anchor=(1.02, 1.00),
                borderaxespad=0.0,
            )
        self._scatter_frame = frame
        self._finalize_axes_layout(bottom=0.18, left=0.12, right=0.76, top=0.82)
        self.canvas.draw_idle()

    def _draw_deg_evidence(self, frame: pd.DataFrame | None, comparison_id: str | None, comparison_name: str | None) -> None:
        self._set_explanation(
            "DEG evidence: top genes ranked by |standardized_log2FC| among DEG-significant rows for the selected comparison."
        )
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        title_name = self._comparison_title_label(comparison_id)
        axis.set_title(f"DEG evidence: {title_name}")
        if frame is None or frame.empty:
            axis.text(0.5, 0.5, "No DEG-significant genes for this comparison.", ha="center", va="center")
            self.canvas.draw_idle()
            return
        working = frame.copy()
        working["gene_label"] = working.get("gene_symbol", pd.Series("", index=working.index)).fillna("").astype(str)
        working.loc[working["gene_label"].eq(""), "gene_label"] = working.get("gene_id", pd.Series("", index=working.index)).fillna("").astype(str)
        working["abs_log2FC"] = pd.to_numeric(working.get("abs_log2FC"), errors="coerce").fillna(0)
        working["standardized_log2FC"] = pd.to_numeric(working.get("standardized_log2FC"), errors="coerce").fillna(0)
        top = working.sort_values(["abs_log2FC", "DE_FDR"], ascending=[False, True], na_position="last").head(20)
        colors = ["#EF4444" if value >= 0 else "#2563EB" for value in top["standardized_log2FC"]]
        axis.bar(top["gene_label"].tolist(), top["abs_log2FC"].tolist(), color=colors)
        axis.set_ylabel("|standardized log2FC|")
        axis.set_xlabel("Gene")
        self._bar_patches = [
            (patch, f"{gene}\nlog2FC: {log2fc:+.3f}\nDE_FDR: {fdr}")
            for patch, gene, log2fc, fdr in zip(
                axis.patches,
                top["gene_label"].tolist(),
                top["standardized_log2FC"].tolist(),
                top["DE_FDR"].tolist(),
                strict=False,
            )
        ]
        self._finalize_axes_layout(rotate_x=True, bottom=0.34)
        self.details.setPlainText(
            "DEG evidence uses standardized log2FC inherited from Pairing/comparison configuration.\n"
            "Red = positive in configured numerator/group1; blue = negative."
        )
        self.canvas.draw_idle()

    def _draw_rmats_evidence(self, frame: pd.DataFrame | None, comparison_id: str | None, comparison_name: str | None) -> None:
        self._set_explanation(
            "rMATS evidence: primary splicing support for the selected comparison, ranked by max |standardized dPSI|."
        )
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        title_name = self._comparison_title_label(comparison_id)
        axis.set_title(f"rMATS evidence: {title_name}")
        if frame is None or frame.empty:
            axis.text(0.5, 0.5, "No rMATS-significant genes for this comparison.", ha="center", va="center")
            self.canvas.draw_idle()
            return
        working = frame.copy()
        working["gene_label"] = working.get("gene_symbol", pd.Series("", index=working.index)).fillna("").astype(str)
        working.loc[working["gene_label"].eq(""), "gene_label"] = working.get("gene_id", pd.Series("", index=working.index)).fillna("").astype(str)
        working["max_abs_dPSI"] = pd.to_numeric(working.get("max_abs_dPSI"), errors="coerce").fillna(0)
        working["dominant_rMATS_standardized_dPSI"] = pd.to_numeric(working.get("dominant_rMATS_standardized_dPSI"), errors="coerce").fillna(0)
        top = working.sort_values(["max_abs_dPSI", "best_rMATS_FDR"], ascending=[False, True], na_position="last").head(20)
        colors = ["#10B981" if value >= 0 else "#8B5CF6" for value in top["dominant_rMATS_standardized_dPSI"]]
        axis.bar(top["gene_label"].tolist(), top["max_abs_dPSI"].tolist(), color=colors)
        axis.set_ylabel("max |standardized dPSI|")
        axis.set_xlabel("Gene")
        self._bar_patches = [
            (patch, f"{gene}\nevent: {event_type}\ndPSI: {dpsi:+.3f}\nrMATS_FDR: {fdr}")
            for patch, gene, event_type, dpsi, fdr in zip(
                axis.patches,
                top["gene_label"].tolist(),
                top["dominant_rMATS_event_type"].tolist(),
                top["dominant_rMATS_standardized_dPSI"].tolist(),
                top["best_rMATS_FDR"].tolist(),
                strict=False,
            )
        ]
        self._finalize_axes_layout(rotate_x=True, bottom=0.34)
        self.details.setPlainText(
            "rMATS is the primary splicing evidence layer. Green = positive standardized dPSI; purple = negative standardized dPSI."
        )
        self.canvas.draw_idle()

    def _draw_dexseq_evidence(self, frame: pd.DataFrame | None, comparison_id: str | None, comparison_name: str | None) -> None:
        self._set_explanation(
            "DEXSeq evidence: exon-usage support genes for the selected comparison, ranked by best DEXSeq q-value."
        )
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        title_name = self._comparison_title_label(comparison_id)
        axis.set_title(f"DEXSeq evidence: {title_name}")
        if frame is None or frame.empty:
            axis.text(0.5, 0.5, "No DEXSeq-support genes for this comparison.", ha="center", va="center")
            self.canvas.draw_idle()
            return
        working = frame.copy()
        working["gene_label"] = working.get("gene_symbol", pd.Series("", index=working.index)).fillna("").astype(str)
        working.loc[working["gene_label"].eq(""), "gene_label"] = working.get("gene_id", pd.Series("", index=working.index)).fillna("").astype(str)
        q = pd.to_numeric(working.get("best_DEXSeq_qvalue"), errors="coerce")
        working["score"] = (-np.log10(q.clip(lower=1e-300)) if not q.empty else 0).fillna(0)
        top = working.sort_values(["best_DEXSeq_qvalue", "n_DEXSeq_significant_exons"], ascending=[True, False], na_position="last").head(20)
        axis.bar(top["gene_label"].tolist(), top["score"].tolist(), color="#F59E0B")
        axis.set_ylabel("-log10(best DEXSeq q)")
        axis.set_xlabel("Gene")
        self._bar_patches = [
            (patch, f"{gene}\nDEXSeq q: {qvalue}\nSignificant exons: {n_exons}")
            for patch, gene, qvalue, n_exons in zip(
                axis.patches,
                top["gene_label"].tolist(),
                top["best_DEXSeq_qvalue"].tolist(),
                top["n_DEXSeq_significant_exons"].tolist(),
                strict=False,
            )
        ]
        self._finalize_axes_layout(rotate_x=True, bottom=0.34)
        self.details.setPlainText("DEXSeq is support evidence only. This page shows genes with strongest exon-usage support.")
        self.canvas.draw_idle()

    def _draw_dtu_evidence(self, frame: pd.DataFrame | None, comparison_id: str | None, comparison_name: str | None) -> None:
        self._set_explanation(
            "DTU / Isoform evidence: transcript-usage support genes for the selected comparison, ranked by best DTU q-value."
        )
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        title_name = self._comparison_title_label(comparison_id)
        axis.set_title(f"DTU / Isoform evidence: {title_name}")
        if frame is None or frame.empty:
            axis.text(0.5, 0.5, "No DTU-support genes for this comparison.", ha="center", va="center")
            self.canvas.draw_idle()
            return
        working = frame.copy()
        working["gene_label"] = working.get("gene_symbol", pd.Series("", index=working.index)).fillna("").astype(str)
        working.loc[working["gene_label"].eq(""), "gene_label"] = working.get("gene_id", pd.Series("", index=working.index)).fillna("").astype(str)
        q = pd.to_numeric(working.get("best_DTU_qvalue"), errors="coerce")
        working["score"] = (-np.log10(q.clip(lower=1e-300)) if not q.empty else 0).fillna(0)
        top = working.sort_values(["best_DTU_qvalue", "n_DTU_significant_isoforms"], ascending=[True, False], na_position="last").head(20)
        axis.bar(top["gene_label"].tolist(), top["score"].tolist(), color="#14B8A6")
        axis.set_ylabel("-log10(best DTU q)")
        axis.set_xlabel("Gene")
        self._bar_patches = [
            (patch, f"{gene}\nDTU q: {qvalue}\nSignificant isoforms: {n_iso}")
            for patch, gene, qvalue, n_iso in zip(
                axis.patches,
                top["gene_label"].tolist(),
                top["best_DTU_qvalue"].tolist(),
                top["n_DTU_significant_isoforms"].tolist(),
                strict=False,
            )
        ]
        self._finalize_axes_layout(rotate_x=True, bottom=0.34)
        self.details.setPlainText("DTU / Isoform support is downstream validation evidence. This page shows genes with the strongest transcript-usage signal.")
        self.canvas.draw_idle()

    def _draw_event_type_subanalysis(self, frame: pd.DataFrame | None, comparison_id: str | None, comparison_name: str | None, event_type: str) -> None:
        self._set_explanation(
            f"{event_type} event-type subanalysis: DEG + DAS integration restricted to representative {event_type} events."
        )
        self._draw_integration_scatter(frame, comparison_name, event_type=event_type)
        self._set_detail_text(
            module_detail("5.2.event"),
            f"{event_type} event-type subanalysis loaded. Open Detail for representative-event and direction rules.",
        )
        self.canvas.draw_idle()

    def _draw_mechanism_support(self, mechanism_support: pd.DataFrame | None) -> None:
        self._scatter_collection = None
        self._scatter_frame = pd.DataFrame()
        self._pie_patches = []
        self._heatmap_image = None
        self._heatmap_rows = []
        self._heatmap_cols = []
        self._heatmap_values = None
        self._set_explanation(
            "Mechanism support summary: this aggregates which genes have supporting evidence from complementary analysis layers, such as DEXSeq, DTU/stageR, SUPPA, and transcript-level evidence."
        )
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        axis.set_title("Mechanism support")
        if mechanism_support is None or mechanism_support.empty:
            axis.text(0.5, 0.5, "No mechanism support available.", ha="center", va="center")
            self._populate_table(pd.DataFrame(), [])
            self.details.setPlainText("Run the Mechanism Support branch to generate support summaries.")
            self.canvas.draw_idle()
            return
        support_columns = [column for column in mechanism_support.columns if column.endswith("_support")]
        counts = {
            column.replace("_support", ""): int(pd.to_numeric(mechanism_support[column], errors="coerce").fillna(0).gt(0).sum())
            for column in support_columns
        }
        axis.bar(list(counts.keys()), list(counts.values()), color=["#577590", "#43AA8B", "#F9C74F", "#F9844A"])
        axis.tick_params(axis="x", rotation=20)
        axis.set_ylabel("Genes with support")
        preferred_columns = [
            "gene_symbol",
            "geneSymbol",
            "gene_id",
            "GeneID",
        ] + [column for column in mechanism_support.columns if column not in {"gene_symbol", "geneSymbol", "gene_id", "GeneID"}]
        self._populate_table(mechanism_support, preferred_columns, preview_name="Mechanism support table")
        self._finalize_axes_layout(rotate_x=True, bottom=0.28)
        self.details.setPlainText("Mechanism-support summary and loaded support table.")
        self.canvas.draw_idle()

    def _draw_candidate_evidence_heatmap(self, frame: pd.DataFrame | None, title: str, meta: dict[str, object]) -> None:
        self._scatter_collection = None
        self._scatter_frame = pd.DataFrame()
        self._bar_patches = []
        self._pie_patches = []
        self._set_explanation(
            "Candidate evidence heatmap: rows are candidate genes from gene_level_integrated_candidates.tsv and columns are evidence axes. Use this to review which genes are supported by DEG, rMATS, DEXSeq and DTU layers in the selected comparison."
        )
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        axis.set_title(title)
        if frame is None or frame.empty:
            axis.text(0.5, 0.5, "No candidate rows available.", ha="center", va="center")
            self._populate_table(pd.DataFrame(), [])
            self._set_detail_text(
                module_detail("5.3.heatmap"),
                "No candidate evidence available. Open Detail for heatmap input rules.",
            )
            self.canvas.draw_idle()
            return
        working = frame.copy()
        working["gene_label"] = working.get("gene_symbol", pd.Series("", index=working.index)).fillna("").astype(str)
        working.loc[working["gene_label"].eq(""), "gene_label"] = working.get("gene_id", pd.Series("", index=working.index)).fillna("").astype(str)
        working = working.loc[working["gene_label"].ne("")].copy()
        if working.empty:
            axis.text(0.5, 0.5, "No gene labels available.", ha="center", va="center")
            self._populate_table(pd.DataFrame(), [])
            self.canvas.draw_idle()
            return
        evidence = pd.DataFrame(index=working["gene_label"].tolist())
        evidence["DEG"] = working["DEG_significant"].fillna(False).astype(int).values if "DEG_significant" in working.columns else 0
        evidence["rMATS"] = working["rMATS_significant"].fillna(False).astype(int).values if "rMATS_significant" in working.columns else 0
        evidence["DEXSeq"] = working["DEXSeq_significant"].fillna(False).astype(int).values if "DEXSeq_significant" in working.columns else 0
        evidence["DTU"] = working["DTU_significant"].fillna(False).astype(int).values if "DTU_significant" in working.columns else 0
        evidence["Evidence count"] = pd.to_numeric(working.get("evidence_count"), errors="coerce").fillna(0).values
        evidence["|log2FC|"] = pd.to_numeric(working.get("abs_log2FC"), errors="coerce").fillna(0).values
        evidence["|dPSI|"] = pd.to_numeric(working.get("max_abs_dPSI"), errors="coerce").fillna(0).values
        for column in ("Evidence count", "|log2FC|", "|dPSI|"):
            max_value = float(evidence[column].max()) if not evidence.empty else 0.0
            if max_value > 0:
                evidence[column] = evidence[column] / max_value
        evidence = evidence.sort_values(["Evidence count", "|dPSI|", "|log2FC|"], ascending=False).head(40)
        image = axis.imshow(evidence.values, cmap="YlOrRd", aspect="auto", vmin=0, vmax=1)
        self._heatmap_image = image
        self._heatmap_rows = evidence.index.astype(str).tolist()
        self._heatmap_cols = evidence.columns.astype(str).tolist()
        self._heatmap_values = evidence.values
        axis.set_yticks(range(len(evidence.index)))
        axis.set_yticklabels(evidence.index)
        axis.set_xticks(range(len(evidence.columns)))
        axis.set_xticklabels(evidence.columns, rotation=25, ha="right")
        self.figure.colorbar(image, ax=axis, fraction=0.03, pad=0.02, label="Normalized evidence strength")
        preview = working[
            [
                column
                for column in [
                    "comparison_id",
                    "rank",
                    "gene_symbol",
                    "candidate_tier",
                    "evidence_class",
                    "evidence_count",
                    "DEG_significant",
                    "rMATS_significant",
                    "DEXSeq_significant",
                    "DTU_significant",
                    "DE_FDR",
                    "best_rMATS_FDR",
                    "best_DEXSeq_qvalue",
                    "best_DTU_qvalue",
                    "candidate_reason",
                ]
                if column in working.columns
            ]
        ].copy()
        output_dir = self.project_service.write_candidate_selection_outputs(
            meta.get("comparison_id"),
            kind="evidence_heatmap",
            selected=working,
            parameters=meta,
        )
        self._populate_table(
            preview,
            list(preview.columns),
            file_path=(output_dir / "heatmap_gene_list.tsv") if output_dir else None,
            output_dir=output_dir,
            preview_name="Evidence heatmap source genes",
        )
        self._finalize_axes_layout(rotate_x=True, bottom=0.24)
        self._set_detail_text(
            "\n".join(
                [
                    module_detail("5.3.heatmap"),
                    "",
                    f"Active gene selection rule: {meta.get('selection_rule', '')}",
                    f"Genes shown: {meta.get('gene_count', 0)}",
                    f"Output folder: {output_dir or 'not available'}",
                ]
            ),
            "5.3.3 Evidence heatmap loaded. Open Detail for gene source, blacklist and output rules.",
        )
        self.canvas.draw_idle()

    def _draw_candidate_program_heatmap(self, frame: pd.DataFrame | None, title: str, meta: dict[str, object]) -> None:
        self._scatter_collection = None
        self._scatter_frame = pd.DataFrame()
        self._bar_patches = []
        self._pie_patches = []
        self._set_explanation(
            "Program heatmap: per-comparison candidate heatmap built from standardized_log2FC and standardized_dPSI-centered evidence. Gene selection is explicit and comparison-specific."
        )
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        axis.set_title(title)
        if frame is None or frame.empty:
            axis.text(0.5, 0.5, "No candidate rows available.", ha="center", va="center")
            self._populate_table(pd.DataFrame(), [])
            self._set_detail_text(
                module_detail("5.3.heatmap"),
                "No candidate program heatmap rows available. Open Detail for heatmap input rules.",
            )
            self.canvas.draw_idle()
            return
        working = frame.copy()
        working["gene_label"] = working.get("gene_symbol", pd.Series("", index=working.index)).fillna("").astype(str)
        working.loc[working["gene_label"].eq(""), "gene_label"] = working.get("gene_id", pd.Series("", index=working.index)).fillna("").astype(str)
        working = working.loc[working["gene_label"].ne("")].copy()
        if working.empty:
            axis.text(0.5, 0.5, "No candidate genes available.", ha="center", va="center")
            self._populate_table(pd.DataFrame(), [])
            self.canvas.draw_idle()
            return
        matrix = pd.DataFrame(index=working["gene_label"].tolist())
        matrix["standardized log2FC"] = pd.to_numeric(working.get("standardized_log2FC"), errors="coerce").fillna(0).values
        matrix["standardized dPSI"] = pd.to_numeric(working.get("dominant_rMATS_standardized_dPSI"), errors="coerce").fillna(0).values
        display_matrix = matrix.copy()
        for column in display_matrix.columns:
            max_abs = float(display_matrix[column].abs().max()) if not display_matrix.empty else 0.0
            if max_abs > 0:
                display_matrix[column] = display_matrix[column] / max_abs
        image = axis.imshow(display_matrix.values, cmap="coolwarm", aspect="auto", vmin=-1, vmax=1)
        self._heatmap_image = image
        self._heatmap_rows = matrix.index.astype(str).tolist()
        self._heatmap_cols = display_matrix.columns.astype(str).tolist()
        self._heatmap_values = matrix.values
        axis.set_yticks(range(len(matrix.index)))
        axis.set_yticklabels(matrix.index)
        axis.set_xticks(range(len(display_matrix.columns)))
        axis.set_xticklabels(display_matrix.columns, rotation=25, ha="right")
        self.figure.colorbar(image, ax=axis, fraction=0.03, pad=0.02, label="Column-normalized effect")
        preview = working[
            [
                column
                for column in [
                    "comparison_id",
                    "gene_symbol",
                    "candidate_tier",
                    "evidence_class",
                    "standardized_log2FC",
                    "dominant_rMATS_standardized_dPSI",
                    "abs_log2FC",
                    "max_abs_dPSI",
                    "direction_class",
                    "candidate_reason",
                ]
                if column in working.columns
            ]
        ].copy()
        output_dir = self.project_service.write_candidate_selection_outputs(
            meta.get("comparison_id"),
            kind="program_heatmap",
            selected=working,
            parameters=meta,
        )
        self._populate_table(
            preview,
            list(preview.columns),
            file_path=(output_dir / "heatmap_gene_list.tsv") if output_dir else None,
            output_dir=output_dir,
            preview_name="Program heatmap source genes",
        )
        self._finalize_axes_layout(rotate_x=True, bottom=0.24)
        self._set_detail_text(
            "\n".join(
                [
                    module_detail("5.3.heatmap"),
                    "",
                    f"Active gene selection rule: {meta.get('selection_rule', '')}",
                    f"Output folder: {output_dir or 'not available'}",
                    "Program heatmap note: log2FC and dPSI are normalized per column for display so their different numeric ranges do not collapse one another.",
                ]
            ),
            "5.3.3 Program heatmap loaded. Open Detail for gene source, blacklist and normalization rules.",
        )
        self.canvas.draw_idle()

    def _draw_candidate_cards(self, payload: dict[str, object] | None, title: str, meta: dict[str, object]) -> None:
        self._scatter_collection = None
        self._scatter_frame = pd.DataFrame()
        self._bar_patches = []
        self._pie_patches = []
        self._heatmap_image = None
        self._heatmap_rows = []
        self._heatmap_cols = []
        self._heatmap_values = None
        self._set_explanation(
            "Gene-centered candidate card: one gene card per selection, with the current comparison tab showing DEG/transcriptome support and the full AS event list for that gene. Select a specific event row if you want to open Sashimi for that event."
        )
        self.figure.clear()
        payload = payload or {}
        self._candidate_card_payload = payload
        event_frame = payload.get("event_frame")
        if not isinstance(event_frame, pd.DataFrame):
            event_frame = pd.DataFrame()
        expression_support = payload.get("expression_support")
        if not isinstance(expression_support, pd.DataFrame):
            expression_support = pd.DataFrame()
        comparison_expression = payload.get("comparison_expression")
        if not isinstance(comparison_expression, pd.DataFrame):
            comparison_expression = pd.DataFrame()
        comparison_dpsi = payload.get("comparison_dpsi")
        if not isinstance(comparison_dpsi, pd.DataFrame):
            comparison_dpsi = pd.DataFrame()
        candidate_row = payload.get("candidate_row", {})
        if not isinstance(candidate_row, dict):
            candidate_row = {}
        deg_row = payload.get("deg_row", {})
        if not isinstance(deg_row, dict):
            deg_row = {}
        dominant_event = payload.get("dominant_event", {})
        if not isinstance(dominant_event, dict):
            dominant_event = {}
        comparison_context = payload.get("comparison_context", {})
        if not isinstance(comparison_context, dict):
            comparison_context = {}
        gene_symbol = self._first_text(payload.get("gene_symbol"), candidate_row.get("gene_symbol"), dominant_event.get("gene_symbol"))
        gene_id = self._first_text(payload.get("gene_id"), candidate_row.get("gene_id"), dominant_event.get("gene_id"))
        if not gene_symbol and not gene_id:
            axis = self.figure.add_subplot(111)
            self._set_axis_title(axis, title)
            axis.axis("off")
            axis.text(0.5, 0.5, "No selected candidate gene is available for card rendering.", ha="center", va="center")
            self._populate_table(pd.DataFrame(), [])
            self.card_open_sashimi_button.setEnabled(False)
            self._set_detail_text(
                module_detail("5.3.cards"),
                "No gene-centered candidate card is available yet. Open Detail for card inputs and behavior.",
            )
            self.canvas.draw_idle()
            return
        comparison_name = self._first_text(comparison_context.get("display_name"), meta.get("comparison_id"))
        dominant_rule = self._safe_text(payload.get("dominant_rule"), "Dominant event rule unavailable.")
        self.card_open_sashimi_button.setEnabled(not event_frame.empty)
        as_direction_short = self._safe_text(comparison_context.get("rmats_canonical_direction"), "NA").replace("rMATS canonical direction: ", "")
        deg_direction_short = self._safe_text(comparison_context.get("deg_final_direction"), "NA").replace("DEG selected/final direction: ", "")

        def _float(value: object) -> float | None:
            series = pd.to_numeric(pd.Series([value]), errors="coerce")
            if series.empty or pd.isna(series.iloc[0]):
                return None
            return float(series.iloc[0])

        def _fmt_num(value: object, digits: int = 3) -> str:
            number = _float(value)
            return f"{number:.{digits}f}" if number is not None else "NA"

        def _fmt_signed(value: object, digits: int = 3) -> str:
            number = _float(value)
            return f"{number:+.{digits}f}" if number is not None else "NA"

        def _fmt_bool(value: object) -> str:
            return "yes" if self._safe_bool(value) else "no"

        def _short_text(value: object, limit: int = 42) -> str:
            text = self._safe_text(value)
            return text if len(text) <= limit else f"{text[: limit - 3]}..."

        def _short_event_label(row: pd.Series) -> str:
            event_type = self._safe_text(row.get("event_type", ""))
            event_id = self._safe_text(row.get("event_id", ""))
            if len(event_id) > 32:
                event_id = f"{event_id[:29]}..."
            return f"{event_type} | {event_id}" if event_type or event_id else "event"

        event_display_frame = event_frame.loc[event_frame["significant_event"].fillna(False)].copy() if "significant_event" in event_frame.columns else event_frame.copy()
        expression_axis, dpsi_axis = self.figure.subplots(
            1,
            2,
            gridspec_kw={"width_ratios": [1.08, 1.0], "wspace": 0.26},
        )
        self.figure.suptitle(f"Gene card: {gene_symbol or gene_id}", fontsize=16, fontweight="bold", y=0.96)
        self._bar_patches = []

        def _plot_label(*parts: object, limit: int = 18) -> str:
            cleaned = []
            for part in parts:
                text = str(part or "").strip()
                if not text:
                    continue
                cleaned.append(text if len(text) <= limit else f"{text[: limit - 3]}...")
            return "\n".join(cleaned) if cleaned else "NA"

        # Left panel: expression bars for the selected gene in each experiment/control group.
        expr_frame = comparison_expression.copy()
        if not expr_frame.empty and {"comparison_display_name", "group", "expr"}.issubset(expr_frame.columns):
            expr_plot = expr_frame.copy()
            expr_plot["expr"] = pd.to_numeric(expr_plot["expr"], errors="coerce")
            expr_plot = expr_plot.dropna(subset=["expr"]).reset_index(drop=True)
            expr_plot["group_key"] = expr_plot["group"].fillna("").astype(str).str.strip()
            expr_plot["_current_priority"] = expr_plot.get("current_comparison", pd.Series(False, index=expr_plot.index)).map(lambda value: 0 if self._safe_bool(value) else 1)
            expr_plot = (
                expr_plot.sort_values(["group_key", "_current_priority", "comparison_display_name"], ascending=[True, True, True], na_position="last")
                .drop_duplicates(subset=["group_key"], keep="first")
                .drop(columns=["_current_priority"], errors="ignore")
                .reset_index(drop=True)
            )
        else:
            expr_plot = pd.DataFrame()
        if expr_plot.empty:
            base_mean = _float(deg_row.get("baseMean"))
            if base_mean is None:
                expression_axis.axis("off")
                expression_axis.text(
                    0.5,
                    0.5,
                    "No transcriptome expression values were found for this gene.",
                    ha="center",
                    va="center",
                    transform=expression_axis.transAxes,
                    wrap=True,
                )
            else:
                bars = expression_axis.bar(["DEG baseMean"], [base_mean], color="#64748B", alpha=0.90)
                self._bar_patches.extend([(bars[0], f"DEG baseMean\nExpression: {_fmt_num(base_mean, digits=2)}")])
                expression_axis.set_ylabel("Expression", fontweight="bold")
                expression_axis.tick_params(axis="x", rotation=18)
                self._set_axis_title(expression_axis, "Expression")
                expression_axis.text(
                    0.02,
                    0.96,
                    "Group means unavailable; showing DEG baseMean fallback.",
                    transform=expression_axis.transAxes,
                    ha="left",
                    va="top",
                    fontsize=8.5,
                    color="#475569",
                )
        else:
            x_values = np.arange(len(expr_plot))
            expr_plot["plot_label"] = expr_plot.apply(
                lambda row: _plot_label(row.get("group"), limit=16),
                axis=1,
            )
            role_colors = {"experiment": "#2563EB", "control": "#94A3B8"}
            bar_colors = [role_colors.get(str(role).lower(), "#64748B") for role in expr_plot.get("group_role", pd.Series("", index=expr_plot.index))]
            edge_colors = ["#F59E0B" if bool(value) else "#E2E8F0" for value in expr_plot.get("current_comparison", pd.Series(False, index=expr_plot.index))]
            line_widths = [2.2 if bool(value) else 0.7 for value in expr_plot.get("current_comparison", pd.Series(False, index=expr_plot.index))]
            bars = expression_axis.bar(
                x_values,
                expr_plot["expr"].tolist(),
                color=bar_colors,
                edgecolor=edge_colors,
                linewidth=line_widths,
                alpha=0.92,
            )
            expression_axis.set_xticks(x_values)
            expression_axis.set_xticklabels(expr_plot["plot_label"].tolist(), rotation=35, ha="right", fontsize=8.5)
            expression_axis.set_ylabel("Mean expression", fontweight="bold")
            self._set_axis_title(expression_axis, "Gene expression across included conditions")
            expression_axis.legend(
                handles=[
                    Patch(facecolor="#2563EB", label="experiment/final DEG side"),
                    Patch(facecolor="#94A3B8", label="control/reference side"),
                    Patch(facecolor="white", edgecolor="#F59E0B", label="current comparison"),
                ],
                loc="upper right",
                fontsize=8,
                frameon=False,
            )
            max_expr = max(float(expr_plot["expr"].max()), 1.0)
            expression_axis.set_ylim(0, max_expr * 1.20)
            for patch, (_, row) in zip(bars, expr_plot.iterrows(), strict=False):
                value = _float(row.get("expr"))
                expression_axis.text(
                    patch.get_x() + patch.get_width() / 2,
                    patch.get_height() + max_expr * 0.02,
                    _fmt_num(value, digits=1),
                    ha="center",
                    va="bottom",
                    fontsize=7.8,
                    rotation=0,
                )
                self._bar_patches.append(
                    (
                        patch,
                        "\n".join(
                            [
                                str(row.get("comparison_display_name", "")),
                                f"group: {row.get('group', '')}",
                                f"mean expression: {_fmt_num(value, digits=3)}",
                                f"samples: {row.get('sample_count', 'NA')}",
                            ]
                        ),
                    )
                )

        # Right panel: dominant significant AS dPSI per comparison.
        dpsi_plot = comparison_dpsi.copy()
        if not dpsi_plot.empty and {"comparison_display_name", "dominant_dPSI"}.issubset(dpsi_plot.columns):
            dpsi_plot["dominant_dPSI"] = pd.to_numeric(dpsi_plot["dominant_dPSI"], errors="coerce")
            dpsi_plot["dominant_FDR"] = pd.to_numeric(dpsi_plot.get("dominant_FDR"), errors="coerce")
            dpsi_plot = dpsi_plot.reset_index(drop=True)
        else:
            dpsi_plot = pd.DataFrame()
        if dpsi_plot.empty:
            dpsi_axis.axis("off")
            dpsi_axis.text(0.5, 0.5, "No significant AS dPSI summary is available for this gene.", ha="center", va="center", transform=dpsi_axis.transAxes, wrap=True)
        else:
            x_values = np.arange(len(dpsi_plot))
            dpsi_values = dpsi_plot["dominant_dPSI"].fillna(0.0).tolist()
            has_event = dpsi_plot["dominant_dPSI"].notna()
            dpsi_colors = [
                "#0F766E" if pd.notna(value) and float(value) >= 0 else "#B91C1C" if pd.notna(value) else "#CBD5E1"
                for value in dpsi_plot["dominant_dPSI"].tolist()
            ]
            edge_colors = ["#F59E0B" if bool(value) else "#E2E8F0" for value in dpsi_plot.get("current_comparison", pd.Series(False, index=dpsi_plot.index))]
            line_widths = [2.2 if bool(value) else 0.7 for value in dpsi_plot.get("current_comparison", pd.Series(False, index=dpsi_plot.index))]
            bars = dpsi_axis.bar(
                x_values,
                dpsi_values,
                color=dpsi_colors,
                edgecolor=edge_colors,
                linewidth=line_widths,
                alpha=0.92,
            )
            dpsi_axis.axhline(0, color="#475569", linewidth=1)
            dpsi_axis.set_xticks(x_values)
            dpsi_axis.set_xticklabels(
                [
                    _plot_label(row.get("comparison_display_name"), row.get("dominant_event_type") or "no sig AS", limit=18)
                    for _, row in dpsi_plot.iterrows()
                ],
                rotation=35,
                ha="right",
                fontsize=8.5,
            )
            dpsi_axis.set_ylabel("standardized dPSI", fontweight="bold")
            self._set_axis_title(dpsi_axis, "Dominant significant AS dPSI")
            finite_abs = pd.to_numeric(dpsi_plot["dominant_dPSI"], errors="coerce").abs().dropna()
            max_abs = max(float(finite_abs.max()) if not finite_abs.empty else 0.10, 0.10)
            dpsi_axis.set_ylim(-max_abs * 1.35, max_abs * 1.35)
            dpsi_axis.legend(
                handles=[
                    Patch(facecolor="#0F766E", label="positive dPSI"),
                    Patch(facecolor="#B91C1C", label="negative dPSI"),
                    Patch(facecolor="#CBD5E1", label="no significant AS"),
                ],
                loc="upper right",
                fontsize=8,
                frameon=False,
            )
            for patch, (_, row), event_exists in zip(bars, dpsi_plot.iterrows(), has_event.tolist(), strict=False):
                value = _float(row.get("dominant_dPSI"))
                label = _fmt_signed(value) if event_exists else "no sig"
                y_offset = max_abs * 0.06
                y_text = (value or 0.0) + (y_offset if (value or 0.0) >= 0 else -y_offset)
                dpsi_axis.text(
                    patch.get_x() + patch.get_width() / 2,
                    y_text,
                    label,
                    ha="center",
                    va="bottom" if (value or 0.0) >= 0 else "top",
                    fontsize=8,
                    fontweight="bold" if event_exists else "normal",
                )
                if event_exists:
                    dpsi_axis.text(
                        patch.get_x() + patch.get_width() / 2,
                        -max_abs * 1.20,
                        f"n={row.get('n_significant_as_events', 0)}",
                        ha="center",
                        va="bottom",
                        fontsize=7.8,
                        color="#475569",
                    )
                self._bar_patches.append(
                    (
                        patch,
                        "\n".join(
                            [
                                str(row.get("comparison_display_name", "")),
                                f"dominant dPSI: {_fmt_signed(value)}",
                                f"FDR: {_fmt_num(row.get('dominant_FDR'))}",
                                f"event: {row.get('dominant_event_type', '')} | {row.get('dominant_event_id', '')}",
                                f"significant AS events: {row.get('n_significant_as_events', 0)}",
                                str(row.get("dominant_inclusion_direction", "")),
                            ]
                        ),
                    )
                )

        preview = event_display_frame[
            [
                column
                for column in [
                    "comparison_display_name",
                    "gene_symbol",
                    "gene_id",
                    "event_type",
                    "event_id",
                    "FDR",
                    "dPSI",
                    "psi_experiment_group",
                    "psi_experiment",
                    "psi_control_group",
                    "psi_control",
                    "direction",
                    "inclusion_direction",
                    "significant_event",
                    "is_dominant",
                    "coordinates",
                ]
                if column in event_frame.columns
            ]
        ].copy()
        if "psi_experiment_group" in preview.columns and "psi_control_group" in preview.columns:
            preview = preview.rename(
                columns={
                    "psi_experiment_group": "experiment_group",
                    "psi_control_group": "control_group",
                }
            )
        self._populate_table(
            preview,
            list(preview.columns),
            preview_name="Gene-centered AS event list",
        )
        self._finalize_axes_layout(left=0.06, right=0.99, top=0.86, bottom=0.28)
        self._set_detail_text(
            "\n".join(
                [
                    module_detail("5.3.cards"),
                    "",
                    f"Current gene: {gene_symbol or gene_id}",
                    f"Current comparison: {comparison_name}",
                    f"Gene selection rule: {meta.get('selection_rule', '')}",
                    f"Dominant event rule: {dominant_rule}",
                    f"AS canonical direction: {as_direction_short}",
                    f"DEG final direction: {deg_direction_short}",
                    f"Candidate tier: {candidate_row.get('candidate_tier', 'NA')}",
                    f"Evidence class: {candidate_row.get('evidence_class', 'NA')}",
                    f"Direction class: {candidate_row.get('direction_class', 'NA')}",
                    f"DEG significant: {_fmt_bool(candidate_row.get('DEG_significant', deg_row.get('DEG_significant')))}",
                    f"DE_FDR: {_fmt_num(candidate_row.get('DE_FDR', deg_row.get('DE_FDR')))}",
                    f"standardized_log2FC: {_fmt_signed(candidate_row.get('standardized_log2FC', deg_row.get('standardized_log2FC')))}",
                    f"DEG baseMean: {_fmt_num(deg_row.get('baseMean'), digits=2)}",
                    "Left plot: unique included experiment/condition mean expression values.",
                    "Right plot: dominant significant AS-event dPSI across selected comparisons.",
                    "The table lists only significant AS events for the current comparison.",
                ]
            ),
            "5.3.4 Gene card loaded. Open Detail for event-list, dominant-event and sashimi-link rules.",
        )
        self.canvas.draw_idle()

    def _draw_isoform_followup(self, frame: pd.DataFrame | None, title: str, meta: dict[str, object]) -> None:
        self._scatter_collection = None
        self._scatter_frame = pd.DataFrame()
        self._pie_patches = []
        self._heatmap_image = None
        self._heatmap_rows = []
        self._heatmap_cols = []
        self._heatmap_values = None
        self._set_explanation(
            "Isoform follow-up: use this after candidate classification. DTU/isoform evidence is a validation layer for selected candidate genes, not the first-pass filter for defining candidates."
        )
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        axis.set_title(title)
        if frame is None or frame.empty:
            axis.text(0.5, 0.5, "No selected candidate genes available for isoform follow-up.", ha="center", va="center")
            self._populate_table(pd.DataFrame(), [])
            self.details.setPlainText(
                "Isoform follow-up needs selected candidate genes or shortlist genes from the current comparison."
            )
            self.canvas.draw_idle()
            return
        working = frame.copy()
        summary = working["isoform_followup_status"].fillna("no DTU evidence").value_counts().sort_values(ascending=False)
        bars = axis.bar(summary.index.astype(str).tolist(), summary.values.tolist(), color=["#43AA8B", "#9CA3AF"][: len(summary)])
        axis.set_ylabel("Candidate genes")
        axis.tick_params(axis="x", rotation=20)
        self._bar_patches = [
            (patch, f"{label}\nGenes: {value}")
            for patch, label, value in zip(bars, summary.index.astype(str).tolist(), summary.values.tolist(), strict=False)
        ]
        preview = working[
            [
                column
                for column in [
                    "comparison_id",
                    "gene_symbol",
                    "gene_id",
                    "candidate_tier",
                    "DTU_significant",
                    "best_DTU_qvalue",
                    "n_DTU_significant_isoforms",
                    "transcript_id",
                    "isoform_id",
                    "DTU_qvalue",
                    "isoform_followup_status",
                ]
                if column in working.columns
            ]
        ].copy()
        output_dir = self.project_service.write_candidate_selection_outputs(
            meta.get("comparison_id"),
            kind="isoform_followup",
            selected=working,
            parameters=meta,
        )
        self._populate_table(
            preview,
            list(preview.columns),
            file_path=(output_dir / "heatmap_gene_list.tsv") if output_dir else None,
            output_dir=output_dir,
            preview_name="Isoform follow-up table",
        )
        self._finalize_axes_layout(rotate_x=True, bottom=0.26)
        self.details.setPlainText(
            f"Isoform follow-up for {meta.get('comparison_id', '')}.\n"
            f"Gene selection rule: {meta.get('selection_rule', '')}\n"
            f"Output folder: {output_dir or 'not available'}\n"
            "If DTU evidence is absent, rows are still shown and marked as 'no DTU evidence' instead of being dropped."
        )
        self.canvas.draw_idle()

    def _draw_text_panel(self, title: str) -> None:
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        axis.axis("off")
        axis.text(0.02, 0.95, title, ha="left", va="top", wrap=True, transform=axis.transAxes)
        self._scatter_collection = None
        self._scatter_frame = pd.DataFrame()
        self._bar_patches = []
        self._pie_patches = []
        self._heatmap_image = None
        self._heatmap_rows = []
        self._heatmap_cols = []
        self._heatmap_values = None
        self._finalize_axes_layout(bottom=0.08)
        self.canvas.draw_idle()

    def _handle_host_resize(self) -> None:
        if not self.isVisible():
            return
        self._post_render_layout()

    def _autosize_figure(self) -> None:
        scale = max(self.plot_scale.value(), 10) / 100.0
        width_scale = max(self.plot_width_scale.value(), 10) / 100.0
        height_scale = max(self.plot_height_scale.value(), 10) / 100.0
        viewport = self.canvas_scroll.viewport()
        available_width = max(int(viewport.width()) - 20, 640)
        available_height = max(int(viewport.height()) - 20, 480)
        width_px = max(int(available_width * scale * width_scale), 640)
        height_px = max(int(available_height * scale * height_scale), 480)
        host_width = width_px
        host_height = height_px
        self.figure.set_size_inches(width_px / self.figure.dpi, height_px / self.figure.dpi, forward=True)
        self.canvas.resize(width_px, height_px)
        self.canvas_host.setMinimumSize(0, 0)
        self.canvas_host.resize(host_width, host_height)

    def _reset_plot_view(self) -> None:
        self._center_plot_view()
        QTimer.singleShot(0, self._center_plot_view)
        QTimer.singleShot(80, self._center_plot_view)

    def _center_plot_view(self) -> None:
        horizontal = self.canvas_scroll.horizontalScrollBar()
        vertical = self.canvas_scroll.verticalScrollBar()
        if self.plot_scale.value() <= 100 and self.plot_width_scale.value() <= 100:
            horizontal.setValue(horizontal.minimum())
        else:
            horizontal.setValue((horizontal.minimum() + horizontal.maximum()) // 2)
        if self.plot_scale.value() <= 100 and self.plot_height_scale.value() <= 100:
            vertical.setValue(vertical.minimum())
        else:
            vertical.setValue((vertical.minimum() + vertical.maximum()) // 2)

    def _redraw_scaled_figure(self) -> None:
        self._post_render_layout()

    def _schedule_post_render_layout(self) -> None:
        self._relayout_generation += 1
        generation = self._relayout_generation
        if not self._relayout_pending:
            self._relayout_pending = True

        def _run() -> None:
            if generation != self._relayout_generation:
                return
            self._post_render_layout()

        def _finish() -> None:
            if generation != self._relayout_generation:
                return
            self._relayout_pending = False

        # The first show/tab-switch often reports an incomplete viewport size on the
        # immediate event loop tick. Run a few lightweight relayout passes after the
        # widget settles so the user does not need to touch width/height controls.
        for delay in (0, 60, 180, 360, 650):
            QTimer.singleShot(delay, _run)
        QTimer.singleShot(760, _finish)

    def _post_render_layout(self) -> None:
        if not self.isVisible():
            return
        self._autosize_figure()
        self._apply_saved_layout()
        self._reset_plot_view()
        self.canvas.draw_idle()

    def _apply_saved_layout(self) -> None:
        self._finalize_axes_layout(
            rotate_x=bool(self._layout_options.get("rotate_x", False)),
            bottom=self._layout_options.get("bottom"),
            left=float(self._layout_options.get("left", 0.12)),
            right=float(self._layout_options.get("right", 0.82)),
            top=float(self._layout_options.get("top", 0.86)),
            remember=False,
        )

    def _finalize_axes_layout(
        self,
        *,
        rotate_x: bool = False,
        bottom: float | None = None,
        left: float = 0.10,
        right: float = 0.80,
        top: float = 0.92,
        remember: bool = True,
    ) -> None:
        axes = list(self.figure.axes)
        if not axes:
            return
        if remember:
            self._layout_options = {
                "rotate_x": rotate_x,
                "bottom": bottom,
                "left": left,
                "right": right,
                "top": top,
            }
        longest_tick = 0
        if rotate_x:
            for axis in axes:
                labels = axis.get_xticklabels()
                for label in labels:
                    label.set_rotation(28)
                    label.set_ha("right")
                    label.set_rotation_mode("anchor")
                    longest_tick = max(longest_tick, len(label.get_text() or ""))
        bottom_margin = bottom if bottom is not None else (0.24 if rotate_x else 0.12)
        if rotate_x and longest_tick >= 18:
            bottom_margin = max(bottom_margin, 0.36)
        elif rotate_x:
            bottom_margin = max(bottom_margin, 0.30)
        self.figure.subplots_adjust(
            left=max(left, 0.10),
            right=min(right, 0.86),
            top=min(top, 0.90),
            bottom=min(bottom_margin, 0.42),
        )

    @staticmethod
    def _event_color_map() -> dict[str, str]:
        return {
            "SE": "#4C78A8",
            "RI": "#F58518",
            "A3SS": "#54A24B",
            "A5SS": "#E45756",
            "MXE": "#9D6BCB",
        }

    def _selected_comparison_display_order(self) -> list[str]:
        return [item.display_resolved_name for item in self.project_service.selected_comparisons_for_display()]

    @staticmethod
    def _wrapped_title(text: str, width: int = 58) -> str:
        clean = str(text or "").strip()
        return textwrap.fill(clean, width=width) if clean else ""

    def _set_axis_title(self, axis, title: str) -> None:
        axis.set_title(self._wrapped_title(title), fontweight="bold", fontsize=16, pad=12, loc="center")

    def _set_explanation(self, text: str) -> None:
        self.explanation.setText(text)

    def _set_detail_text(self, text: str, inline_summary: str | None = None) -> None:
        cleaned = str(text or "").strip()
        self._detail_dialog_text = cleaned
        self.details.setPlainText(inline_summary or ("Open Detail for full definitions and interpretation." if cleaned else ""))

    def _sync_footer_info(self) -> None:
        blocks: list[str] = []
        context = self.context_header.text().strip()
        explanation = self.explanation.text().strip()
        if context:
            blocks.append(context)
        if explanation:
            blocks.append(explanation)
        self.footer_info.setPlainText("\n\n".join(blocks))

    def _show_detail_dialog(self) -> None:
        blocks = []
        footer_text = self.footer_info.toPlainText().strip()
        detail_text = self._detail_dialog_text.strip()
        inline_text = self.details.toPlainText().strip()
        if footer_text:
            blocks.append(footer_text)
        if detail_text:
            blocks.append(detail_text)
        elif inline_text and inline_text != footer_text:
            blocks.append(inline_text)
        text = "\n\n".join(blocks).strip()
        if not text:
            return
        dialog = QDialog(self)
        dialog.setWindowTitle("Comparison details")
        dialog.resize(760, 420)
        layout = QVBoxLayout(dialog)
        viewer = QTextEdit(dialog)
        viewer.setReadOnly(True)
        viewer.setPlainText(text)
        viewer.setStyleSheet("font-size: 11px;")
        close_button = QPushButton("Close", dialog)
        close_button.clicked.connect(dialog.accept)
        layout.addWidget(viewer, 1)
        layout.addWidget(close_button)
        QTimer.singleShot(0, dialog.showMaximized)
        dialog.exec()

    @staticmethod
    def _table_item_for_value(frame: pd.DataFrame, column: str, value: object) -> QTableWidgetItem:
        item = QTableWidgetItem()
        if pd.isna(value):
            item.setText("")
        elif column in frame.columns and pd.api.types.is_bool_dtype(frame[column]):
            item.setData(Qt.ItemDataRole.EditRole, bool(value))
        elif column in frame.columns and pd.api.types.is_numeric_dtype(frame[column]):
            numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
            if pd.notna(numeric):
                item.setData(Qt.ItemDataRole.EditRole, float(numeric) if not float(numeric).is_integer() else int(numeric))
            else:
                item.setText(str(value))
        else:
            numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
            if pd.notna(numeric) and str(value).strip() not in {"", "nan", "None"}:
                item.setData(Qt.ItemDataRole.EditRole, float(numeric) if not float(numeric).is_integer() else int(numeric))
            else:
                item.setText(str(value))
        item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
        return item

    def _show_preview_dialog(self) -> None:
        frame = self._preview_filtered_frame.copy() if not self._preview_filtered_frame.empty else self._preview_source_frame.copy()
        dialog = QDialog(self)
        dialog.setWindowTitle(self._preview_name or "Data preview")
        dialog.resize(1080, 680)
        dialog.setWindowFlags(
            dialog.windowFlags()
            | Qt.WindowType.Window
            | Qt.WindowType.WindowMaximizeButtonHint
            | Qt.WindowType.WindowMinimizeButtonHint
        )
        layout = QVBoxLayout(dialog)

        path_label = QLineEdit(dialog)
        path_label.setReadOnly(True)
        path_label.setText(str(self._preview_file_path) if self._preview_file_path else "In-memory preview (no backing file)")
        path_label.setStyleSheet("font-size: 10px;")
        layout.addWidget(path_label)

        info = QLabel(dialog)
        info.setWordWrap(True)
        info.setStyleSheet("font-size: 10px; color: #bdbdbd;")
        layout.addWidget(info)

        table = QTableWidget(dialog)
        table.setSortingEnabled(True)
        table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectItems)
        table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        table.setHorizontalScrollMode(QAbstractItemView.ScrollMode.ScrollPerPixel)
        table.setVerticalScrollMode(QAbstractItemView.ScrollMode.ScrollPerPixel)
        table.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        table.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        layout.addWidget(table, 1)

        actions = QHBoxLayout()
        open_button = QPushButton("Open output folder", dialog)
        open_button.clicked.connect(self._open_preview_output_folder)
        download_button = QPushButton("Download file", dialog)
        download_button.clicked.connect(self._download_preview_file)
        export_button = QPushButton("Export current filtered table", dialog)
        export_button.clicked.connect(self._export_filtered_preview)
        actions.addWidget(open_button)
        actions.addWidget(download_button)
        actions.addWidget(export_button)
        close_button = QPushButton("Close", dialog)
        close_button.clicked.connect(dialog.accept)
        actions.addWidget(close_button)
        layout.addLayout(actions)

        if frame.empty:
            table.setRowCount(0)
            table.setColumnCount(0)
            info.setText(f"{self._preview_name}: no rows available.")
        else:
            display = frame[self._preview_columns].copy() if self._preview_columns else frame.copy()
            truncated = False
            if len(display) > self._preview_max_rows:
                display = display.head(self._preview_max_rows).copy()
                truncated = True
            table.setColumnCount(len(display.columns))
            table.setHorizontalHeaderLabels(display.columns.tolist())
            table.setRowCount(len(display))
            for row_idx, (_, row) in enumerate(display.iterrows()):
                for col_idx, column in enumerate(display.columns):
                    table.setItem(row_idx, col_idx, self._table_item_for_value(display, column, row[column]))
                info.setText(
                    f"{self._preview_name}: {len(frame)} filtered row(s)"
                    + (f"; preview limited to first {self._preview_max_rows} rows" if truncated else "")
                )

        QTimer.singleShot(0, dialog.showMaximized)
        dialog.exec()

    def _draw_placeholder(self, text: str) -> None:
        self._draw_text_panel(text)
        self._populate_table(pd.DataFrame(), [], preview_name="No preview available")
        self.details.setPlainText(text)

    def _populate_table(
        self,
        frame: pd.DataFrame | None,
        preferred_columns: list[str],
        *,
        file_path: str | Path | None = None,
        output_dir: str | Path | None = None,
        preview_name: str | None = None,
    ) -> None:
        self._preview_source_frame = pd.DataFrame() if frame is None else frame.copy()
        self._preview_filtered_frame = pd.DataFrame()
        self._preview_file_path = Path(file_path) if file_path else None
        self._preview_output_dir = Path(output_dir) if output_dir else (self._preview_file_path.parent if self._preview_file_path else None)
        self._preview_name = preview_name or "table preview"
        if frame is None or frame.empty:
            self.preview_table.setRowCount(0)
            self.preview_table.setColumnCount(0)
            self.preview_info.setText(
                f"{self._preview_name}: no rows available."
                + (f" File: {self._preview_file_path}" if self._preview_file_path else "")
            )
            self.preview_path.setText(str(self._preview_file_path) if self._preview_file_path else "")
            self.preview_filter_column.blockSignals(True)
            self.preview_filter_column.clear()
            self.preview_filter_column.addItem("No filters")
            self.preview_filter_column.blockSignals(False)
            self.preview_filter_value.blockSignals(True)
            self.preview_filter_value.clear()
            self.preview_filter_value.addItem("All values")
            self.preview_filter_value.blockSignals(False)
            return
        self._preview_columns = [column for column in preferred_columns if column in frame.columns] or list(frame.columns)
        self.preview_path.setText(str(self._preview_file_path) if self._preview_file_path else "In-memory preview (no backing file)")
        self.preview_filter_column.blockSignals(True)
        current_column = self.preview_filter_column.currentText()
        self.preview_filter_column.clear()
        self.preview_filter_column.addItem("No filters")
        for column in self._preview_columns:
            self.preview_filter_column.addItem(column)
        if current_column:
            index = self.preview_filter_column.findText(current_column)
            if index >= 0:
                self.preview_filter_column.setCurrentIndex(index)
        self.preview_filter_column.blockSignals(False)
        self._refresh_preview_filter_values()
        self._apply_preview_filters()

    def _apply_preview_filters(self) -> None:
        if self._preview_source_frame.empty:
            return
        filtered = self._preview_source_frame.copy()
        filter_column = self.preview_filter_column.currentText().strip()
        filter_value = self.preview_filter_value.currentText().strip()
        if filter_column and filter_column != "No filters" and filter_value and filter_value != "All values" and filter_column in filtered.columns:
            filtered = filtered.loc[filtered[filter_column].astype(str) == filter_value].copy()
        query = self.preview_search.text().strip().lower()
        if query:
            search_columns = [
                column
                for column in ("gene_symbol", "geneSymbol", "gene_id", "GeneID", "comparison_id", "comparison_name", "candidate_tier", "direction_class")
                if column in filtered.columns
            ]
            if not search_columns:
                search_columns = list(filtered.columns[: min(6, len(filtered.columns))])
            mask = pd.Series(False, index=filtered.index)
            for column in search_columns:
                mask = mask | filtered[column].astype(str).str.lower().str.contains(query, na=False)
            filtered = filtered.loc[mask].copy()
        self._preview_filtered_frame = filtered.copy()
        display = filtered[self._preview_columns].copy() if self._preview_columns else filtered.copy()
        truncated = False
        if len(display) > self._preview_max_rows:
            display = display.head(self._preview_max_rows).copy()
            truncated = True
        self.preview_table.setColumnCount(len(display.columns))
        self.preview_table.setHorizontalHeaderLabels(display.columns.tolist())
        self.preview_table.setRowCount(len(display))
        for row_idx, (_, row) in enumerate(display.iterrows()):
            for col_idx, column in enumerate(display.columns):
                self.preview_table.setItem(row_idx, col_idx, self._table_item_for_value(display, column, row[column]))
        file_status = "in-memory result"
        if self._preview_file_path is not None:
            file_status = f"file: {self._preview_file_path}"
            if not self._preview_file_path.exists():
                file_status += " (missing file)"
        info = f"{self._preview_name}: {len(self._preview_filtered_frame)} filtered row(s)"
        if truncated:
            info += f"; preview limited to first {self._preview_max_rows} rows"
        info += f". Source = {file_status}."
        self.preview_info.setText(info)

    def _refresh_preview_filter_values(self) -> None:
        self.preview_filter_value.blockSignals(True)
        current_value = self.preview_filter_value.currentText()
        self.preview_filter_value.clear()
        self.preview_filter_value.addItem("All values")
        column = self.preview_filter_column.currentText().strip()
        if column and column != "No filters" and not self._preview_source_frame.empty and column in self._preview_source_frame.columns:
            values = sorted(self._preview_source_frame[column].dropna().astype(str).unique().tolist())
            for value in values[:500]:
                self.preview_filter_value.addItem(value)
        if current_value:
            index = self.preview_filter_value.findText(current_value)
            if index >= 0:
                self.preview_filter_value.setCurrentIndex(index)
        self.preview_filter_value.blockSignals(False)
        self._apply_preview_filters()

    def _filter_frame_by_value(self, frame: pd.DataFrame | None, column: str, value: str | None) -> pd.DataFrame:
        if frame is None or frame.empty:
            return pd.DataFrame()
        if not value or column not in frame.columns:
            return frame.copy()
        return frame.loc[frame[column].astype(str) == str(value)].copy()

    def _current_tab_key(self, widget: QTabWidget) -> str | None:
        if widget.count() == 0:
            return None
        value = widget.tabToolTip(widget.currentIndex())
        return value or None

    def _pair_name_from_id(self, pair_id: str | None) -> str | None:
        project = self.project_service.current_project
        if project is None or not pair_id:
            return None
        pair = next((item for item in project.comparison_pairs if item.pair_id == pair_id), None)
        return pair.resolved_name if pair is not None else None

    def _pair_definition_from_id(self, pair_id: str | None):
        project = self.project_service.current_project
        if project is None or not pair_id:
            return None
        return next((item for item in project.comparison_pairs if item.pair_id == pair_id), None)

    def _pairwise_candidate_matrix(self, matrix: pd.DataFrame | None, pair_id: str | None) -> pd.DataFrame:
        pair = self._pair_definition_from_id(pair_id)
        if matrix is None or matrix.empty or pair is None or not pair.comparison_a or not pair.comparison_b:
            return pd.DataFrame()
        comp_a = str(pair.comparison_a)
        comp_b = str(pair.comparison_b)
        name_a = self._comparison_name_from_id(comp_a) or comp_a
        name_b = self._comparison_name_from_id(comp_b) or comp_b

        def value(row: pd.Series, prefix: str, comparison_id: str, default=""):
            column = f"{prefix}_{comparison_id}"
            if column in row.index:
                return row.get(column, default)
            return default

        def truthy(value_obj) -> bool:
            if value_obj is None:
                return False
            try:
                if pd.isna(value_obj):
                    return False
            except (TypeError, ValueError):
                pass
            text = str(value_obj).strip().lower()
            return text in {"1", "true", "yes", "y"} or value_obj is True

        def sign(value_obj) -> int:
            numeric = pd.to_numeric(value_obj, errors="coerce")
            if pd.isna(numeric) or float(numeric) == 0:
                return 0
            return 1 if float(numeric) > 0 else -1

        rows: list[dict[str, object]] = []
        for _, row in matrix.iterrows():
            candidate_a = truthy(value(row, "is_candidate", comp_a, False))
            candidate_b = truthy(value(row, "is_candidate", comp_b, False))
            if not (candidate_a or candidate_b):
                continue
            direction_a = str(value(row, "direction_class", comp_a, "") or "")
            direction_b = str(value(row, "direction_class", comp_b, "") or "")
            tier_a = str(value(row, "tier", comp_a, "") or "")
            tier_b = str(value(row, "tier", comp_b, "") or "")
            evidence_a = str(value(row, "evidence_class", comp_a, "") or "")
            evidence_b = str(value(row, "evidence_class", comp_b, "") or "")
            dpsi_a = pd.to_numeric(value(row, "dominant_rMATS_standardized_dPSI", comp_a, pd.NA), errors="coerce")
            dpsi_b = pd.to_numeric(value(row, "dominant_rMATS_standardized_dPSI", comp_b, pd.NA), errors="coerce")
            logfc_a = pd.to_numeric(value(row, "standardized_log2FC", comp_a, pd.NA), errors="coerce")
            logfc_b = pd.to_numeric(value(row, "standardized_log2FC", comp_b, pd.NA), errors="coerce")
            deg_direction_reversed = bool(
                candidate_a
                and candidate_b
                and sign(logfc_a) != 0
                and sign(logfc_b) != 0
                and sign(logfc_a) != sign(logfc_b)
            )
            dominant_as_direction_reversed = bool(
                candidate_a
                and candidate_b
                and sign(dpsi_a) != 0
                and sign(dpsi_b) != 0
                and sign(dpsi_a) != sign(dpsi_b)
            )
            gene_level_direction_changed = deg_direction_reversed or dominant_as_direction_reversed
            tier_changed = candidate_a and candidate_b and tier_a and tier_b and tier_a != tier_b
            evidence_changed = candidate_a and candidate_b and evidence_a and evidence_b and evidence_a != evidence_b
            if candidate_a and candidate_b and deg_direction_reversed and dominant_as_direction_reversed:
                pattern = "shared_mixed_direction_changed"
                reason = (
                    "Candidate is present in both comparisons, and both gene-level DEG direction "
                    "and dominant AS direction change sign between A and B. This is not strict same-event AS reversal."
                )
            elif candidate_a and candidate_b and deg_direction_reversed:
                pattern = "shared_DEG_direction_reversed"
                reason = (
                    "Candidate is present in both comparisons, and standardized DEG log2FC changes sign between A and B."
                )
            elif candidate_a and candidate_b and dominant_as_direction_reversed:
                pattern = "shared_dominant_AS_direction_reversed"
                reason = (
                    "Candidate is present in both comparisons, and dominant rMATS standardized dPSI changes sign between A and B. "
                    "The dominant AS event may differ by comparison; use 5.4.3/5.4.4 for strict same-event reversal."
                )
            elif candidate_a and candidate_b and tier_changed:
                pattern = "shared_tier_changed"
                reason = "Candidate is present in both comparisons, but candidate tier differs between A and B."
            elif candidate_a and candidate_b and evidence_changed:
                pattern = "shared_evidence_class_changed"
                reason = "Candidate is present in both comparisons, but evidence class differs between A and B."
            elif candidate_a and candidate_b:
                pattern = "shared_pair"
                reason = "Candidate is present in both sides of the selected comparison pair."
            elif candidate_a:
                pattern = "comparison_A_specific_lost_in_B"
                reason = "Candidate is present in comparison A but absent in comparison B for this selected pair."
            else:
                pattern = "comparison_B_specific_gained_in_B"
                reason = "Candidate is absent in comparison A but present in comparison B for this selected pair."
            pair_membership_count = int(candidate_a) + int(candidate_b)
            if candidate_a and candidate_b:
                pair_membership_label = "shared A+B"
            elif candidate_a:
                pair_membership_label = "A only"
            else:
                pair_membership_label = "B only"
            rows.append(
                {
                    "pair_id": pair.pair_id,
                    "pair_name": pair.resolved_name,
                    "comparison_A": comp_a,
                    "comparison_B": comp_b,
                    "comparison_A_display_name": name_a,
                    "comparison_B_display_name": name_b,
                    "gene_id": row.get("gene_id", ""),
                    "gene_symbol": row.get("gene_symbol", ""),
                    "candidate_pattern": pattern,
                    "pattern_reason": reason,
                    "candidate_comparisons": ";".join([name for name, present in [(name_a, candidate_a), (name_b, candidate_b)] if present]),
                    "pair_membership_label": pair_membership_label,
                    "pair_membership_count": pair_membership_count,
                    "n_comparisons_candidate": pair_membership_count,
                    "is_candidate_A": candidate_a,
                    "is_candidate_B": candidate_b,
                    "tier_A": tier_a,
                    "tier_B": tier_b,
                    "rank_A": value(row, "rank", comp_a, pd.NA),
                    "rank_B": value(row, "rank", comp_b, pd.NA),
                    "evidence_class_A": evidence_a,
                    "evidence_class_B": evidence_b,
                    "direction_class_A": direction_a,
                    "direction_class_B": direction_b,
                    "standardized_log2FC_A": logfc_a,
                    "standardized_log2FC_B": logfc_b,
                    "dPSI_A": dpsi_a,
                    "dPSI_B": dpsi_b,
                    "DE_FDR_A": value(row, "DE_FDR", comp_a, pd.NA),
                    "DE_FDR_B": value(row, "DE_FDR", comp_b, pd.NA),
                    "rMATS_FDR_A": value(row, "best_rMATS_FDR", comp_a, pd.NA),
                    "rMATS_FDR_B": value(row, "best_rMATS_FDR", comp_b, pd.NA),
                    "event_type_A": value(row, "dominant_rMATS_event_type", comp_a, ""),
                    "event_type_B": value(row, "dominant_rMATS_event_type", comp_b, ""),
                    "pair_gene_level_direction_changed": gene_level_direction_changed,
                    "pair_DEG_direction_reversed": deg_direction_reversed,
                    "pair_dominant_AS_direction_reversed": dominant_as_direction_reversed,
                    "pair_tier_changed": tier_changed,
                    "pair_evidence_class_changed": evidence_changed,
                    "best_pair_rank": pd.to_numeric(pd.Series([value(row, "rank", comp_a, pd.NA), value(row, "rank", comp_b, pd.NA)]), errors="coerce").min(),
                }
            )
        if not rows:
            return pd.DataFrame()
        result = pd.DataFrame(rows)
        return result.sort_values(["candidate_pattern", "best_pair_rank", "gene_symbol"], na_position="last").reset_index(drop=True)

    def _comparison_name_from_id(self, comparison_id: str | None) -> str | None:
        project = self.project_service.current_project
        if project is None or not comparison_id:
            return None
        comparison = next((item for item in project.available_comparisons if item.comparison_id == comparison_id), None)
        return comparison.resolved_name if comparison is not None else None

    def _comparison_title_label(self, comparison_id: str | None) -> str:
        if not comparison_id:
            return "selected comparison"
        name = self._comparison_name_from_id(comparison_id)
        if name:
            return f"{name} [{comparison_id}]"
        return comparison_id

    def _candidate_file_map(self, comparison_id: str | None) -> dict[str, Path]:
        return self.project_service.candidate_screening_file_map(comparison_id)

    def _cross_file_map(self) -> dict[str, Path]:
        return self.project_service.cross_comparison_file_map()

    def _cross_pattern_file_map(self, pair_id: str | None) -> dict[str, Path]:
        return self.project_service.cross_comparison_pattern_file_map(pair_id)

    def _cross_interaction_file_map(self, pair_id: str | None) -> dict[str, Path]:
        return self.project_service.cross_comparison_interaction_file_map(pair_id)

    def _cross_interaction_batch_file_map(self) -> dict[str, Path]:
        return self.project_service.cross_comparison_interaction_batch_file_map()

    def _apply_group_pattern_recalculate(self) -> None:
        pair_id = self._current_tab_key(self.group_tabs)
        self.project_service.preview_cross_comparison_significant_as_patterns(
            pair_id,
            abs_dpsi_cutoff=float(self.group_pattern_sig_cutoff.value()),
            large_delta_dpsi_cutoff=float(self.group_pattern_large_delta.value()),
            allow_generate=True,
            force_rebuild=True,
        )
        self.refresh()

    def _apply_group_interaction_recalculate(self) -> None:
        pair_id = self._current_tab_key(self.group_tabs)
        selected_key = str(self.group_interaction_view.currentData() or "unified_event_table")
        if selected_key == "batch_summary":
            self.project_service.preview_cross_comparison_as_interaction_batch(
                abs_dpsi_cutoff=float(self.group_interaction_sig_cutoff.value()),
                large_delta_dpsi_cutoff=float(self.group_interaction_delta_cutoff.value()),
                min_total_reads=int(self.group_interaction_min_reads.value()),
                allow_generate=True,
                force_rebuild=True,
            )
        else:
            self.project_service.preview_cross_comparison_as_interactions(
                pair_id,
                abs_dpsi_cutoff=float(self.group_interaction_sig_cutoff.value()),
                large_delta_dpsi_cutoff=float(self.group_interaction_delta_cutoff.value()),
                min_total_reads=int(self.group_interaction_min_reads.value()),
                allow_generate=True,
                force_rebuild=True,
            )
        self.refresh()

    def _apply_group_venn_recalculate(self) -> None:
        pair_id = self._current_tab_key(self.group_tabs)
        summary, events = self._build_pair_venn_tables(pair_id)
        signature = self._group_venn_signature()
        self._group_venn_cache[signature] = {
            "summary": summary,
            "events": events,
        }
        self.refresh()

    def _apply_group_venn_multi_recalculate(self) -> None:
        summary, events, selected_names = self._build_group_venn_multi_tables()
        self._group_venn_multi_cache = {
            "signature": self._group_venn_multi_signature(),
            "summary": summary,
            "events": events,
            "selected_names": selected_names,
        }
        self.refresh()

    def _candidate_custom_gene_tokens(self) -> list[str]:
        text = self.candidate_custom_genes.text().strip()
        if not text:
            return []
        parts = re.split(r"[\s,\t;|]+", text)
        seen: set[str] = set()
        tokens: list[str] = []
        for part in parts:
            value = part.strip()
            if not value:
                continue
            key = value.casefold()
            if key in seen:
                continue
            seen.add(key)
            tokens.append(value)
        return tokens

    def _on_candidate_card_gene_changed(self) -> None:
        self._candidate_card_selected_gene_key = str(self.card_gene_combo.currentData() or "")
        self.refresh()

    def _sync_candidate_card_gene_selector(self, payload: dict[str, object]) -> None:
        options = payload.get("gene_options", []) if isinstance(payload, dict) else []
        selected_key = str(payload.get("selected_gene_key", "") or "")
        self.card_gene_combo.blockSignals(True)
        self.card_gene_combo.clear()
        for item in options:
            if not isinstance(item, dict):
                continue
            self.card_gene_combo.addItem(str(item.get("label", "")), str(item.get("key", "")))
        if self.card_gene_combo.count():
            target_index = max(self.card_gene_combo.findData(selected_key), 0)
            self.card_gene_combo.setCurrentIndex(target_index)
            self._candidate_card_selected_gene_key = str(self.card_gene_combo.currentData() or "")
        else:
            self._candidate_card_selected_gene_key = None
        self.card_gene_combo.blockSignals(False)

    def _open_candidate_card_event_in_sashimi(self) -> None:
        event_id = ""
        comparison_id = ""
        preview = self._preview_filtered_frame.copy()
        if not preview.empty and "event_id" in preview.columns:
            selected_rows = sorted({index.row() for index in self.preview_table.selectedIndexes()})
            if selected_rows:
                row = min(selected_rows[0], len(preview) - 1)
                if row >= 0:
                    event_id = str(preview.iloc[row].get("event_id", "") or "")
                    comparison_id = str(preview.iloc[row].get("comparison_id", "") or "")
        if not event_id and isinstance(self._candidate_card_payload, dict):
            dominant = self._candidate_card_payload.get("dominant_event", {})
            if isinstance(dominant, dict):
                event_id = str(dominant.get("event_id", "") or "")
                comparison_id = str(dominant.get("comparison_id", "") or comparison_id)
        if not event_id:
            self.details.setPlainText("No event is selected for Sashimi. Select an event row from the gene card event list first.")
            return
        opener = getattr(self.window(), "_open_named_page", None)
        if callable(opener):
            opener("5.5.2 Run Sashimi")
        sashimi_page = getattr(self.window(), "sashimi_page", None)
        if sashimi_page is not None:
            tabs = getattr(sashimi_page, "comparison_tabs", None)
            if tabs is not None and comparison_id:
                for index in range(tabs.count()):
                    if tabs.tabToolTip(index) == comparison_id:
                        tabs.setCurrentIndex(index)
                        break
            search = getattr(sashimi_page, "search", None)
            if search is not None:
                search.setText(event_id)
            refresh = getattr(sashimi_page, "refresh", None)
            if callable(refresh):
                refresh()
        self.details.setPlainText(
            f"Open Sashimi for the currently selected gene event.\n"
            f"Selected event_id: {event_id}\n"
            "Gene cards do not auto-run Sashimi; use the Sashimi page to review and generate the plot for this event."
        )

    def _on_motion(self, event) -> None:
        if self._pan_active and self._pan_axes is not None and event.inaxes == self._pan_axes and event.xdata is not None and event.ydata is not None and self._pan_origin is not None and self._pan_xlim is not None and self._pan_ylim is not None:
            dx = float(event.xdata) - self._pan_origin[0]
            dy = float(event.ydata) - self._pan_origin[1]
            self._pan_axes.set_xlim(self._pan_xlim[0] - dx, self._pan_xlim[1] - dx)
            self._pan_axes.set_ylim(self._pan_ylim[0] - dy, self._pan_ylim[1] - dy)
            self._pan_dragged = True
            self.canvas.draw_idle()
            return
        if event.inaxes is None or self._hover_annotation is None:
            return
        for patch, tooltip in self._bar_patches:
            contains, _ = patch.contains(event)
            if contains:
                x = patch.get_x() + patch.get_width() / 2
                y = patch.get_height()
                self._hover_annotation.xy = (x, y)
                self._hover_annotation.set_text(tooltip)
                self._hover_annotation.set_visible(True)
                self.canvas.draw_idle()
                return
        for patch, tooltip in self._pie_patches:
            contains, _ = patch.contains(event)
            if contains:
                center = getattr(patch, "center", (0.0, 0.0))
                self._hover_annotation.xy = center
                self._hover_annotation.set_text(tooltip)
                self._hover_annotation.set_visible(True)
                self.canvas.draw_idle()
                return
        if self._heatmap_image is not None and self._heatmap_values is not None and event.xdata is not None and event.ydata is not None:
            col = int(round(float(event.xdata)))
            row = int(round(float(event.ydata)))
            if 0 <= row < len(self._heatmap_rows) and 0 <= col < len(self._heatmap_cols):
                value = float(self._heatmap_values[row, col])
                self._hover_annotation.xy = (col, row)
                self._hover_annotation.set_text(
                    f"{self._heatmap_rows[row]}\n{self._heatmap_cols[col]}\nvalue={value:.3f}"
                )
                self._hover_annotation.set_visible(True)
                self.canvas.draw_idle()
                return
        if self._scatter_collection is not None and not self._scatter_frame.empty:
            contains, info = self._scatter_collection.contains(event)
            if contains and info.get("ind"):
                index = int(info["ind"][0])
                row = self._scatter_frame.iloc[index]
                self._hover_annotation.xy = (float(row["standardized_log2FC"]), float(row["standardized_dPSI"]))
                self._hover_annotation.set_text(
                    f"{row.get('gene_symbol', row.get('geneSymbol', row.get('gene_id', '')))}\n"
                    f"standardized log2FC={row.get('standardized_log2FC', '')}\n"
                    f"standardized dPSI={row.get('standardized_dPSI', '')}"
                )
                self._hover_annotation.set_visible(True)
                self.canvas.draw_idle()
                return
        if self._hover_annotation.get_visible():
            self._hover_annotation.set_visible(False)
            self.canvas.draw_idle()

    def _on_press(self, event) -> None:
        if event.button == 1 and event.inaxes is not None and event.xdata is not None and event.ydata is not None:
            self._pan_active = True
            self._pan_dragged = False
            self._pan_axes = event.inaxes
            self._pan_origin = (float(event.xdata), float(event.ydata))
            self._pan_xlim = tuple(event.inaxes.get_xlim())
            self._pan_ylim = tuple(event.inaxes.get_ylim())

    def _on_release(self, event) -> None:
        was_active = self._pan_active
        was_dragged = self._pan_dragged
        self._pan_active = False
        self._pan_dragged = False
        self._pan_axes = None
        self._pan_origin = None
        self._pan_xlim = None
        self._pan_ylim = None
        if was_active and not was_dragged:
            self._handle_click(event)

    def _on_scroll(self, event) -> None:
        if event.inaxes is None or event.xdata is None or event.ydata is None:
            return
        axis = event.inaxes
        step = getattr(event, "step", 0) or 0
        zoom_in = step > 0 or getattr(event, "button", "") == "up"
        scale_factor = 0.9 if zoom_in else 1.1
        xlim = axis.get_xlim()
        ylim = axis.get_ylim()
        xdata = float(event.xdata)
        ydata = float(event.ydata)
        new_width = (xlim[1] - xlim[0]) * scale_factor
        new_height = (ylim[1] - ylim[0]) * scale_factor
        xfrac = 0.5 if xlim[1] == xlim[0] else (xdata - xlim[0]) / (xlim[1] - xlim[0])
        yfrac = 0.5 if ylim[1] == ylim[0] else (ydata - ylim[0]) / (ylim[1] - ylim[0])
        axis.set_xlim([xdata - new_width * xfrac, xdata + new_width * (1 - xfrac)])
        axis.set_ylim([ydata - new_height * yfrac, ydata + new_height * (1 - yfrac)])
        self.canvas.draw_idle()

    def _handle_click(self, event) -> None:
        if self._scatter_collection is None or self._scatter_frame.empty:
            return
        contains, info = self._scatter_collection.contains(event)
        if not contains or not info.get("ind"):
            return
        index = int(info["ind"][0])
        row = self._scatter_frame.iloc[index]
        lines = [
            f"Comparison: {row.get('comparison_name', '')}",
            f"Gene: {row.get('gene_symbol', row.get('geneSymbol', row.get('gene_id', '')))}",
            f"Gene ID: {row.get('gene_id', '')}",
            f"standardized log2FC: {row.get('standardized_log2FC', row.get('log2FC', ''))}",
            f"DEG padj: {row.get('DE_FDR', row.get('deg_padj', ''))}",
            f"standardized dPSI: {row.get('standardized_dPSI', row.get('representative_dPSI', ''))}",
            f"Representative event FDR: {row.get('rMATS_FDR', row.get('representative_event_FDR', ''))}",
            f"DEG/DAS class: {row.get('plot_class', row.get('combined_class', ''))}",
            f"Agreement: {row.get('agreement_class', '')}",
            f"Significant splicing events: {row.get('n_sig_events_all', '')}",
            f"log2FC direction: {row.get('log2fc_direction_label', '')}",
            f"dPSI direction: {row.get('dpsi_direction_label', '')}",
        ]
        self.details.setPlainText("\n".join(lines))

    def _export_current_figure(self, fmt: str) -> None:
        if fmt not in {"png", "pdf"}:
            return
        output_path, _ = QFileDialog.getSaveFileName(
            self,
            f"Export Analysis as {fmt.upper()}",
            str(self.project_service.default_export_path(f"analysis_export.{fmt}")),
            f"{fmt.upper()} Files (*.{fmt})",
        )
        if not output_path:
            return
        self.figure.savefig(output_path, dpi=300)

    def _open_preview_output_folder(self) -> None:
        target = self._preview_output_dir or (self._preview_file_path.parent if self._preview_file_path else None)
        if target is None or not target.exists():
            self.details.setPlainText("Preview output folder is not available.")
            return
        os.startfile(str(target))

    def _download_preview_file(self) -> None:
        source_path = self._preview_file_path if self._preview_file_path and self._preview_file_path.exists() else None
        suffix = source_path.suffix if source_path else ".tsv"
        output_path, _ = QFileDialog.getSaveFileName(
            self,
            "Download Preview File",
            str(self.project_service.default_export_path(f"preview_export{suffix}")),
            "All Files (*.*)",
        )
        if not output_path:
            return
        destination = Path(output_path)
        if source_path is not None:
            shutil.copy2(source_path, destination)
            return
        frame = self._preview_source_frame if not self._preview_source_frame.empty else self._preview_filtered_frame
        if frame.empty:
            return
        self.project_service.export_frame(frame, destination)

    def _export_filtered_preview(self) -> None:
        if self._preview_filtered_frame.empty:
            return
        output_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export Current Filtered Table",
            str(self.project_service.default_export_path("filtered_preview.tsv")),
            "TSV Files (*.tsv)",
        )
        if not output_path:
            return
        self.project_service.export_frame(self._preview_filtered_frame, output_path)

    def _copy_preview_selection(self) -> None:
        indexes = self.preview_table.selectedIndexes()
        if not indexes:
            return
        rows = sorted({index.row() for index in indexes})
        cols = sorted({index.column() for index in indexes})
        lines = []
        header = [self.preview_table.horizontalHeaderItem(col).text() if self.preview_table.horizontalHeaderItem(col) else "" for col in cols]
        lines.append("\t".join(header))
        for row in rows:
            values = []
            for col in cols:
                item = self.preview_table.item(row, col)
                values.append(item.text() if item is not None else "")
            lines.append("\t".join(values))
        QGuiApplication.clipboard().setText("\n".join(lines))
