from __future__ import annotations

import sys
from pathlib import Path

from PySide6.QtCore import QObject, QThread, Signal
from PySide6.QtGui import QIcon
from PySide6.QtWidgets import (
    QApplication,
    QFileDialog,
    QFrame,
    QLabel,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QSplitter,
    QStatusBar,
    QToolBar,
    QVBoxLayout,
    QWidget,
    QProgressBar,
    QSizePolicy,
)

from app.pages.comparison_page import ComparisonPage
from app.pages.comparison_sets_page import ComparisonSetsPage
from app.pages.cosplice_page import CoSplicePage
from app.pages.cross_comparison_page import CrossComparisonResultsPage
from app.pages.deg_page import DegPage
from app.pages.enrichment_page import EnrichmentPage
from app.pages.candidate_gene_page import CandidateGenePage
from app.pages.candidate_selection_page import CandidateSelectionPage
from app.pages.cards_page import CardsPage
from app.pages.blacklist_page import BlacklistPage
from app.pages.input_check_page import InputCheckPage
from app.pages.input_paths_dialog import InputPathsDialog
from app.pages.isoform_page import IsoformPage
from app.pages.jutils_page import JutilsPage
from app.pages.project_page import ProjectPage
from app.pages.rbp_page import RbpPage
from app.pages.results_page import ResultsPage
from app.pages.section_navigation_page import SectionNavigationPage
from app.pages.sashimi_page import SashimiPage
from app.pages.shortlist_page import ShortlistPage
from app.pages.samples_page import SamplesPage
from app.pages.settings_page import SettingsPage
from app.pages.thresholds_page import ThresholdsPage
from app.pages.visualization_groups_page import VisualizationGroupsPage
from app.pages.workflow_dashboard_page import WorkflowDashboardPage
from app.theme import KATANA_STYLESHEET
from app.widgets.navigation_tree import WorkflowNavigation
from app.widgets.responsive_stack import ResponsiveStackedWidget
from src.runtime_paths import bundled_assets_root
from src.services.project_service import ProjectService


class ProjectLoadWorker(QObject):
    finished = Signal(object)
    failed = Signal(str)
    progress = Signal(str)

    def __init__(
        self,
        project_service: ProjectService,
        folder: str | Path,
        input_paths: dict[str, str] | None = None,
        force_auto_scan: bool = False,
    ) -> None:
        super().__init__()
        self.project_service = project_service
        self.folder = folder
        self.input_paths = input_paths or {}
        self.force_auto_scan = force_auto_scan

    def run(self) -> None:
        try:
            self.progress.emit("start load project")
            project = self.project_service.load_project(
                self.folder,
                input_paths=self.input_paths,
                force_auto_scan=self.force_auto_scan,
                progress_callback=self.progress.emit,
            )
            self.progress.emit("finish load project")
        except Exception as exc:
            self.failed.emit(str(exc))
            return
        self.finished.emit(project)


class AnalysisRunWorker(QObject):
    finished = Signal(object)
    failed = Signal(str)
    progress = Signal(str)

    def __init__(self, project_service: ProjectService) -> None:
        super().__init__()
        self.project_service = project_service

    def run(self) -> None:
        try:
            run_state = self.project_service.run_main_analysis_with_progress(self.progress.emit)
        except Exception as exc:
            self.failed.emit(str(exc))
            return
        self.finished.emit(run_state)


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("Katana Splicing Tool")
        self.resize(1360, 860)
        icon_path = bundled_assets_root() / "katana_icon.png"
        if icon_path.exists():
            self.setWindowIcon(QIcon(str(icon_path)))

        self.project_service = ProjectService()
        self.pages = ResponsiveStackedWidget()
        self.nav = WorkflowNavigation()
        self.nav.setMinimumWidth(252)
        self.nav.setMaximumWidth(360)
        self.nav.currentRowChanged.connect(self._switch_page)
        self._page_routes: list[tuple[int, QWidget, tuple[str, str | None] | None]] = []
        self._named_nav_rows: dict[str, int] = {}
        self.toolbar_title: QLabel | None = None

        self.project_page = ProjectPage(self.project_service, self._advance_after_project_confirmation)
        self.comparison_page = ComparisonPage(self.project_service, self._advance_after_pairing_confirmation)
        self.comparison_sets_page = ComparisonSetsPage(self.project_service, self._advance_after_comparison_sets_confirmation)
        self.visualization_groups_page = VisualizationGroupsPage(self.project_service, self._advance_after_visual_groups_confirmation)
        self.per_comparison_page = WorkflowDashboardPage(self.project_service, self._open_named_page)
        self.overview_page = SectionNavigationPage(
            "Overview",
            "Overview is the first per-comparison layer. Use Landscape for a comparison-level overview, Bar for significant event counts, and Pie for event-type composition.",
            [
                ("Landscape", "Single-comparison overview with total events, significant events and significant fractions.", "5.1.1 Landscape"),
                ("Bar", "Single-comparison significant event counts by event type.", "5.1.2 Bar"),
                ("Pie", "Event-type composition pie chart per comparison.", "5.1.3 Pie"),
            ],
            self._open_named_page,
        )
        self.tx_hub_page = SectionNavigationPage(
            "Transcript + Splicing Integration",
            "This section strictly inherits Pairing direction. First open the overall transcript + splicing integration for one comparison, then open event-type subanalysis only if you need SE/RI/A3SS/A5SS/MXE-specific views. Thresholds shown in the analysis sidebar are used when these analyses are rerun.",
            [
                ("Transcript + Splicing Integration", "Main interactive DEG + splicing integration view for one comparison at a time.", "5.2.1 Transcript + Splicing Integration"),
                ("Event-type Subanalysis", "Open SE/RI/A3SS/A5SS/MXE-specific subanalysis after double-checking the comparison direction.", "5.2.2 Event-type Subanalysis"),
            ],
            self._open_named_page,
        )
        self.tx_event_hub_page = SectionNavigationPage(
            "Event-type Subanalysis",
            "Open one splicing event type at a time. This subanalysis is separate from the main DEG + splicing page so it does not slow down the main 5.2 view.",
            [
                ("SE", "Skipped-exon DEG + DAS integration view.", "5.2.2.1 SE"),
                ("RI", "Retained-intron DEG + DAS integration view.", "5.2.2.2 RI"),
                ("A3SS", "Alternative 3' splice-site DEG + DAS integration view.", "5.2.2.3 A3SS"),
                ("A5SS", "Alternative 5' splice-site DEG + DAS integration view.", "5.2.2.4 A5SS"),
                ("MXE", "Mutually exclusive exon DEG + DAS integration view.", "5.2.2.5 MXE"),
            ],
            self._open_named_page,
            columns=1,
        )
        self.candidate_screen_page = SectionNavigationPage(
            "Candidate Gene Screening",
            "This is the core per-comparison candidate layer. First review the ranked candidate table, then build the per-comparison candidate gene selection that feeds the evidence heatmap and candidate/event follow-up views.",
            [
                ("Candidate Ranking Table", "Rank per-comparison genes using DEG+rMATS+DEXSeq+DTU evidence.", "5.3.1 Candidate Ranking Table"),
                ("Candidate Gene Selection", "Per-comparison candidate gene selection with default top 20, manual add/remove and per-comparison blacklist.", "5.3.2 Candidate Gene Selection"),
                ("Evidence Heatmap", "Heatmap of multi-layer evidence across candidate genes.", "5.3.3 Evidence Heatmap"),
                ("Candidate Cards + Event Follow-up", "Candidate cards driven from the selected genes; event-level follow-up continues into sashimi.", "5.3.4 Candidate Cards + Event Follow-up"),
                ("Candidate Gene AS Heatmap", "Sample-level PSI and dPSI heatmaps restricted to the genes selected in 5.3.2.", "5.3.5 Candidate Gene AS Heatmap"),
            ],
            self._open_named_page,
        )
        self.cross_comparison_page = SectionNavigationPage(
            "Cross-comparison Candidate Comparison",
            "Cross-comparison is split by level: 5.4.1/5.4.2 are pairwise gene-level candidate A-vs-B views, and 5.4.3/5.4.4/5.4.5/5.4.6 are strict AS event-level views.",
            [
                ("Candidate Matrix", "Pairwise A-vs-B matrix of candidate membership for the selected comparison pair.", "5.4.1 Candidate Matrix"),
                ("Gene Pattern Explanation", "Pairwise shared / A-specific / B-specific / gained / lost explanation for gene-level candidate patterns.", "5.4.2 Shared / Specific / Gained / Lost"),
                ("AS interaction / read support", "Replicate-aware AS interaction follow-up for shared significant events, A/B-specific event support and effect-size patterns.", "5.4.3 Cross-comparison AS interaction modeling"),
                ("Strict-event Venn / overlap", "Strict rMATS event-level A-specific / shared / B-specific overlap.", "5.4.4 Venn / Overlap"),
                ("Selected-comparison Venn", "Choose multiple comparisons and show one large direction-aware AS-event Venn.", "5.4.5 Selected-comparison Event Venn"),
                ("Selected-comparison Similarity Heatmap", "Choose multiple comparisons and show one direction-aware AS-event similarity heatmap.", "5.4.6 Selected-comparison Event Similarity"),
            ],
            self._open_named_page,
        )
        self.sashimi_hub_page = SectionNavigationPage(
            "Sashimi",
            "Sashimi is the candidate-event validation layer after shortlist/blacklist. Start with input checks, run selected rMATS events, then preview generated PNGs.",
            [
                ("Input Check", "Confirm candidate events, BAM lists and manifest inputs.", "5.5.1 Input Check"),
                ("Run Sashimi", "Open the sashimi runner and launch rmats2sashimi.", "5.5.2 Run Sashimi"),
                ("Sashimi Preview", "Preview generated rMATS event PNGs.", "5.5.3 Sashimi Preview"),
            ],
            self._open_named_page,
        )
        self.landscape_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="landscape", fixed_subbranch="overview")
        self.landscape_bar_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="landscape", fixed_subbranch="bar")
        self.landscape_pie_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="landscape", fixed_subbranch="pie")
        self.group_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="group")
        self.group_matrix_page = CrossComparisonResultsPage(self.project_service, self.run_analysis, subbranch="matrix")
        self.group_shared_page = CrossComparisonResultsPage(self.project_service, self.run_analysis, subbranch="shared")
        self.group_reversal_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="group", fixed_subbranch="reversal")
        self.group_comparison_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="group", fixed_subbranch="comparison")
        self.group_direction_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="group", fixed_subbranch="direction")
        self.group_pattern_page = CrossComparisonResultsPage(self.project_service, self.run_analysis, subbranch="pattern")
        self.group_interaction_page = CrossComparisonResultsPage(self.project_service, self.run_analysis, subbranch="interaction")
        self.group_venn_page = CrossComparisonResultsPage(self.project_service, self.run_analysis, subbranch="venn")
        self.group_venn_multi_page = CrossComparisonResultsPage(self.project_service, self.run_analysis, subbranch="venn_multi")
        self.group_venn_similarity_page = CrossComparisonResultsPage(
            self.project_service,
            self.run_analysis,
            subbranch="venn_multi_similarity",
        )
        self.tx_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="tx")
        self.tx_direction_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="tx", fixed_subbranch="direction")
        self.tx_deg_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="tx", fixed_subbranch="deg")
        self.tx_rmats_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="tx", fixed_subbranch="rmats")
        self.tx_dexseq_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="tx", fixed_subbranch="dexseq")
        self.tx_dtu_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="tx", fixed_subbranch="dtu")
        self.tx_se_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="tx", fixed_subbranch="se")
        self.tx_ri_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="tx", fixed_subbranch="ri")
        self.tx_a3ss_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="tx", fixed_subbranch="a3ss")
        self.tx_a5ss_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="tx", fixed_subbranch="a5ss")
        self.tx_mxe_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="tx", fixed_subbranch="mxe")
        self.mechanism_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="mechanism")
        self.candidate_hook_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="candidate")
        self.candidate_hook_table_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="candidate", fixed_subbranch="table")
        self.candidate_hook_heatmap_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="candidate", fixed_subbranch="heatmap")
        self.candidate_hook_program_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="candidate", fixed_subbranch="program")
        self.candidate_hook_cards_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="candidate", fixed_subbranch="cards")
        self.candidate_hook_isoform_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="candidate", fixed_subbranch="isoform")
        self.candidate_ranking_page = CandidateGenePage(self.project_service, mode="ranking")
        self.candidate_selection_page = CandidateSelectionPage(self.project_service)
        self.candidate_tier_page = CandidateGenePage(self.project_service, mode="tiers")
        self.shortlist_page = ShortlistPage(self.project_service)
        self.blacklist_page = BlacklistPage(self.project_service)
        self.cards_page = CardsPage(self.project_service)
        self.sashimi_page = SashimiPage(self.project_service, view_mode="run")
        self.sashimi_preview_page = SashimiPage(self.project_service, view_mode="preview")
        self.samples_page = SamplesPage(self.project_service)
        self.input_check_page = InputCheckPage(self.project_service)
        self.rbp_page = RbpPage(self.project_service)
        self.rbp_per_page = RbpPage(self.project_service, initial_tab=0)
        self.rbp_cross_page = RbpPage(self.project_service, initial_tab=1)
        self.isoform_page = IsoformPage(self.project_service)
        self.isoform_per_page = IsoformPage(self.project_service, initial_tab=1)
        self.isoform_cross_page = IsoformPage(self.project_service, initial_tab=2)
        self.candidate_as_heatmap_page = JutilsPage(self.project_service, view_mode="candidate_heatmap")
        self.native_pattern_overview_page = JutilsPage(self.project_service, view_mode="native_pattern")
        self.native_psi_heatmap_page = JutilsPage(self.project_service, view_mode="native_heatmap")
        self.native_psi_pca_page = JutilsPage(self.project_service, view_mode="native_pca")
        self.native_similarity_page = JutilsPage(self.project_service, view_mode="native_similarity")
        self.native_output_browser_page = JutilsPage(self.project_service, view_mode="native_browser")
        self.thresholds_page = ThresholdsPage(self.project_service)
        self.deg_page = DegPage(self.project_service)
        self.enrichment_page = EnrichmentPage(self.project_service)
        self.cosplice_page = CoSplicePage(self.project_service)
        self.settings_page = SettingsPage(self.project_service, self.reload_current_project)

        self._add_page("1. Project", self.project_page)
        self._add_page("2. Pairing", self.comparison_page)
        self._add_page("3. Comparison Sets", self.comparison_sets_page)
        self._add_page("4. Visual Groups", self.visualization_groups_page)
        self._add_page("5. Per-comparison Analysis", self.per_comparison_page)
        self._add_page("5.0 Analysis Thresholds", self.thresholds_page)
        self._add_page("5.1 Overview", self.overview_page)
        self._add_page("5.1.1 Landscape", self.landscape_page, ("landscape", "overview"))
        self._add_page("5.1.2 Bar", self.landscape_bar_page, ("landscape", "bar"))
        self._add_page("5.1.3 Pie", self.landscape_pie_page, ("landscape", "pie"))
        self._add_page("5.2 Transcript + Splicing Integration", self.tx_hub_page)
        self._add_page("5.2.1 Transcript + Splicing Integration", self.tx_page)
        self._add_page("5.2.2 Event-type Subanalysis", self.tx_event_hub_page)
        self._add_page("5.2.2.1 SE", self.tx_se_page)
        self._add_page("5.2.2.2 RI", self.tx_ri_page)
        self._add_page("5.2.2.3 A3SS", self.tx_a3ss_page)
        self._add_page("5.2.2.4 A5SS", self.tx_a5ss_page)
        self._add_page("5.2.2.5 MXE", self.tx_mxe_page)
        self._add_page("5.3 Candidate Gene Screening", self.candidate_screen_page)
        self._add_page("5.3.1 Candidate Ranking Table", self.candidate_ranking_page)
        self._add_page("5.3.2 Candidate Gene Selection", self.candidate_selection_page)
        self._add_page("5.3.3 Evidence Heatmap", self.candidate_hook_heatmap_page, ("candidate", "heatmap"))
        self._add_page("5.3.4 Candidate Cards + Event Follow-up", self.candidate_hook_cards_page, ("candidate", "cards"))
        self._add_page("5.3.5 Candidate Gene AS Heatmap", self.candidate_as_heatmap_page)
        self._add_page("5.4 Cross-comparison Candidate Comparison", self.cross_comparison_page)
        self._add_page("5.4.1 Candidate Matrix", self.group_matrix_page)
        self._add_page("5.4.2 Shared / Specific / Gained / Lost", self.group_shared_page)
        self._add_page("5.4.3 Cross-comparison AS interaction modeling", self.group_interaction_page)
        self._add_page("5.4.4 Venn / Overlap", self.group_venn_page)
        self._add_page("5.4.5 Selected-comparison Event Venn", self.group_venn_multi_page)
        self._add_page("5.4.6 Selected-comparison Event Similarity", self.group_venn_similarity_page)
        self._add_page("5.5 Sashimi", self.sashimi_hub_page)
        self._add_page("5.5.1 Input Check", self.input_check_page)
        self._add_page("5.5.2 Run Sashimi", self.sashimi_page)
        self._add_page("5.5.3 Sashimi Preview", self.sashimi_preview_page)
        self._add_page("5.6 Native Splicing Pattern Explorer", self.native_pattern_overview_page)
        self._add_page("5.6.1 Sample PSI Heatmaps", self.native_psi_heatmap_page)
        self._add_page("5.6.2 PCA & Replicate QC", self.native_psi_pca_page)
        self._add_page("5.6.3 Similarity & Overlap", self.native_similarity_page)
        self._add_page("5.6.4 Native Output Browser", self.native_output_browser_page)
        self._add_page("5.7 DEG Explorer", self.deg_page)
        self._add_page("5.8 GO / KEGG Enrichment", self.enrichment_page)
        self._add_page("6. Samples", self.samples_page)
        self._add_page("7. RBP", self.rbp_page)
        self._add_page("7.1 Per-comparison RBP support", self.rbp_per_page, ("rbp", "per"))
        self._add_page("7.2 Cross-comparison RBP support", self.rbp_cross_page, ("rbp", "cross"))
        self._add_page("8. Isoform", self.isoform_page)
        self._add_page("8.1 Per-comparison isoform support", self.isoform_per_page, ("isoform", "per"))
        self._add_page("8.2 Cross-comparison isoform support", self.isoform_cross_page, ("isoform", "cross"))
        self._add_page("9. CoSplice Network", self.cosplice_page)
        self._add_page("10. Input Check", self.input_check_page)
        self._add_page("11. Settings", self.settings_page)
        self.nav.expand_primary_workflow()
        self.nav.setCurrentRow(0)

        left_panel = QFrame()
        left_panel.setObjectName("sidebar")
        left_panel.setMinimumWidth(252)
        left_panel.setMaximumWidth(360)
        left_layout = QVBoxLayout(left_panel)
        left_layout.setContentsMargins(16, 18, 12, 14)
        left_layout.setSpacing(12)

        brand = QLabel("KATANA")
        brand.setObjectName("brandTitle")
        brand_subtitle = QLabel("RNA splicing workspace")
        brand_subtitle.setObjectName("brandSubtitle")
        self.project_badge = QLabel("No project loaded")
        self.project_badge.setObjectName("projectBadge")
        self.project_badge.setWordWrap(True)
        workflow_label = QLabel("WORKFLOW")
        workflow_label.setObjectName("sidebarOverline")

        left_layout.addWidget(brand)
        left_layout.addWidget(brand_subtitle)
        left_layout.addWidget(self.project_badge)
        left_layout.addSpacing(6)
        left_layout.addWidget(workflow_label)
        left_layout.addWidget(self.nav, 1)

        content_surface = QFrame()
        content_surface.setObjectName("contentSurface")
        content_layout = QVBoxLayout(content_surface)
        content_layout.setContentsMargins(14, 14, 14, 14)
        content_layout.addWidget(self.pages, 1)

        splitter = QSplitter()
        splitter.setObjectName("workspaceSplitter")
        splitter.addWidget(left_panel)
        splitter.addWidget(content_surface)
        splitter.setStretchFactor(1, 1)
        splitter.setChildrenCollapsible(False)
        splitter.setHandleWidth(1)
        splitter.setSizes([300, 1060])
        self.setCentralWidget(splitter)

        toolbar = QToolBar("Main")
        toolbar.setObjectName("mainToolbar")
        toolbar.setMovable(False)
        self.addToolBar(toolbar)

        self.toolbar_title = QLabel("1. Project")
        self.toolbar_title.setObjectName("toolbarTitle")
        toolbar.addWidget(self.toolbar_title)
        toolbar_spacer = QWidget()
        toolbar_spacer.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        toolbar.addWidget(toolbar_spacer)

        self.open_action = QPushButton("Open Project")
        self.open_action.setProperty("class", "primary")
        self.open_action.clicked.connect(self.open_project)
        toolbar.addWidget(self.open_action)

        self.run_action = QPushButton("Run current page analysis")
        self.run_action.clicked.connect(self.run_analysis)
        self.run_action.setToolTip(
            "Hidden legacy global runner. Use the Run / Apply / Generate button inside each analysis page instead."
        )
        toolbar.addWidget(self.run_action)
        self.run_action.setVisible(False)

        self.setStatusBar(QStatusBar(self))
        self.progress_bar = QProgressBar(self)
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setVisible(False)
        self.progress_bar.setMaximumWidth(220)
        self.statusBar().addPermanentWidget(self.progress_bar)
        self.statusBar().showMessage("Ready")

        self._loader_thread: QThread | None = None
        self._loader_worker: ProjectLoadWorker | None = None
        self._analysis_thread: QThread | None = None
        self._analysis_worker: AnalysisRunWorker | None = None

    def _add_page(self, name: str, widget: QWidget, analysis_route: tuple[str, str | None] | None = None) -> None:
        page_index = self.pages.indexOf(widget)
        if page_index < 0:
            self.pages.addWidget(widget)
            page_index = self.pages.indexOf(widget)
        self._page_routes.append((page_index, widget, analysis_route))
        route_index = len(self._page_routes) - 1
        self.nav.add_route(route_index, name)
        self._named_nav_rows[name] = route_index

    def _switch_page(self, index: int) -> None:
        if index > 0 and self.project_service.current_project is not None and not self.project_service.current_project.confirmed:
            QMessageBox.information(
                self,
                "Confirm Project First",
                "Confirm the detected inputs on the Project page before moving to the next step.",
            )
            self.nav.blockSignals(True)
            self.nav.setCurrentRow(0)
            self.nav.blockSignals(False)
            index = 0
        elif index > 1 and self.project_service.current_project is not None and not self.project_service.current_project.pairing_confirmed:
            QMessageBox.information(
                self,
                "Confirm Pairing First",
                "Confirm the pairing on the Pairing page before moving to the next step.",
            )
            self.nav.blockSignals(True)
            self.nav.setCurrentRow(1)
            self.nav.blockSignals(False)
            index = 1
        elif index > 2 and self.project_service.current_project is not None and not self.project_service.current_project.comparison_sets_confirmed:
            QMessageBox.information(
                self,
                "Confirm Comparison Sets First",
                "Confirm the comparison sets on the Comparison Sets page before moving to Analysis and downstream pages.",
            )
            self.nav.blockSignals(True)
            self.nav.setCurrentRow(2)
            self.nav.blockSignals(False)
            index = 2
        elif index > 3 and self.project_service.current_project is not None and not self.project_service.current_project.visualization_groups_confirmed:
            QMessageBox.information(
                self,
                "Confirm Visual Groups First",
                "Confirm the visual groups before moving to Analysis and downstream pages.",
            )
            self.nav.blockSignals(True)
            self.nav.setCurrentRow(3)
            self.nav.blockSignals(False)
            index = 3
        if 0 <= index < len(self._page_routes):
            if self.toolbar_title is not None:
                current_name = next(
                    (name for name, route_index in self._named_nav_rows.items() if route_index == index),
                    "Workspace",
                )
                self.toolbar_title.setText(current_name)
            page_index, page, analysis_route = self._page_routes[index]
            self.pages.setCurrentIndex(page_index)
            if analysis_route is not None:
                branch, subbranch = analysis_route
                navigate_to = getattr(page, "navigate_to", None)
                if callable(navigate_to):
                    navigate_to(branch, subbranch)
            refresh = getattr(page, "refresh", None)
            if callable(refresh):
                refresh()

    def open_project(self) -> None:
        folder = QFileDialog.getExistingDirectory(self, "Select Project Directory", str(Path.home()))
        if not folder:
            return
        self._open_project_with_dialog(folder)

    def reload_current_project(self) -> None:
        project = self.project_service.current_project
        if project is None:
            QMessageBox.information(self, "No project", "Open a project before reloading it.")
            return
        self._load_project(project.project_root)

    def _open_project_with_dialog(self, folder: str | Path) -> None:
        dialog = InputPathsDialog(Path(folder))
        if dialog.exec() == InputPathsDialog.DialogCode.Rejected:
            return
        self._load_project(folder, dialog.selected_paths(), dialog.force_auto_scan)

    def _load_project(
        self,
        folder: str | Path,
        input_paths: dict[str, str] | None = None,
        force_auto_scan: bool = False,
    ) -> None:
        if self._loader_thread is not None:
            QMessageBox.information(self, "Project loading", "A project scan is already running.")
            return
        self.statusBar().showMessage(f"Scanning project: {folder}")
        self.open_action.setEnabled(False)
        self.run_action.setEnabled(False)

        self._loader_thread = QThread(self)
        self._loader_worker = ProjectLoadWorker(self.project_service, folder, input_paths, force_auto_scan)
        self._loader_worker.moveToThread(self._loader_thread)
        self._loader_thread.started.connect(self._loader_worker.run)
        self._loader_worker.progress.connect(self._on_project_progress)
        self._loader_worker.finished.connect(self._on_project_loaded)
        self._loader_worker.failed.connect(self._on_project_failed)
        self._loader_worker.finished.connect(self._loader_thread.quit)
        self._loader_worker.failed.connect(self._loader_thread.quit)
        self._loader_thread.finished.connect(self._cleanup_loader)
        self._loader_thread.start()

    def _on_project_progress(self, message: str) -> None:
        print(f"[project-load] {message}")
        self.statusBar().showMessage(message)
        self.progress_bar.setVisible(True)

    def _on_project_loaded(self, _project) -> None:
        output_error = self.project_service.output_root_validation_error()
        if output_error:
            selected_output = QFileDialog.getExistingDirectory(
                self,
                "Select Katana Output Directory",
                str(Path(_project.project_root).parent),
            )
            if not selected_output:
                self.project_service.current_project = None
                self.statusBar().showMessage("Project load cancelled: output directory is required")
                QMessageBox.warning(
                    self,
                    "Output Directory Required",
                    "Katana will not write analysis results into the tool directory. "
                    "Select a separate output directory to load this project.",
                )
                return
            try:
                if self.project_service.can_migrate_current_output_root():
                    self.project_service.migrate_output_root(selected_output)
                else:
                    self.project_service.update_output_root(selected_output)
            except RuntimeError as exc:
                self.project_service.current_project = None
                self.statusBar().showMessage("Project load cancelled: invalid output directory")
                QMessageBox.critical(self, "Invalid Output Directory", str(exc))
                return
        self.statusBar().showMessage("Project loaded")
        project_name = getattr(_project, "project_name", None) or "Loaded project"
        self.project_badge.setText(project_name)
        self.per_comparison_page.refresh()
        self.pages.setCurrentIndex(0)
        self.project_page.refresh()
        self.comparison_page.refresh()
        self.comparison_sets_page.refresh()
        self.visualization_groups_page.refresh()
        self.thresholds_page.refresh()
        print("[project-load] lightweight UI refresh only; heavy 5.x pages will lazy-load on first open")
        self.nav.blockSignals(True)
        self.nav.setCurrentRow(0)
        self.nav.blockSignals(False)
        self.run_action.setEnabled(True)

    def _on_project_failed(self, message: str) -> None:
        self.statusBar().showMessage("Project load failed")
        QMessageBox.critical(self, "Failed to load project", message)

    def _cleanup_loader(self) -> None:
        self.open_action.setEnabled(True)
        self.run_action.setEnabled(self.project_service.current_project is not None)
        self.progress_bar.setVisible(False)
        if self._loader_worker is not None:
            self._loader_worker.deleteLater()
        if self._loader_thread is not None:
            self._loader_thread.deleteLater()
        self._loader_worker = None
        self._loader_thread = None

    def run_analysis(self) -> None:
        if self.project_service.current_project is None:
            QMessageBox.information(self, "No project", "Open a project before running analysis.")
            return
        if not self.project_service.current_project.confirmed:
            QMessageBox.information(self, "Confirm Project First", "Confirm the project on the Project page before running analysis.")
            self.nav.setCurrentRow(0)
            return
        if not self.project_service.current_project.pairing_confirmed:
            QMessageBox.information(self, "Confirm Pairing First", "Confirm the pairing on the Pairing page before running analysis.")
            self.nav.setCurrentRow(1)
            return
        if not self.project_service.current_project.comparison_sets_confirmed:
            QMessageBox.information(self, "Confirm Comparison Sets First", "Confirm the comparison sets on the Comparison Sets page before running analysis.")
            self.nav.setCurrentRow(2)
            return
        if not self.project_service.current_project.visualization_groups_confirmed:
            QMessageBox.information(self, "Confirm Visual Groups First", "Confirm the visual groups on the Visual Groups page before running analysis.")
            self.nav.setCurrentRow(3)
            return
        if self._analysis_thread is not None:
            QMessageBox.information(self, "Analysis running", "Analysis is already running.")
            return
        self.statusBar().showMessage("Running current analysis request...")
        self.open_action.setEnabled(False)
        self.run_action.setEnabled(False)
        self.progress_bar.setVisible(True)

        self._analysis_thread = QThread(self)
        self._analysis_worker = AnalysisRunWorker(self.project_service)
        self._analysis_worker.moveToThread(self._analysis_thread)
        self._analysis_thread.started.connect(self._analysis_worker.run)
        self._analysis_worker.progress.connect(self._on_analysis_progress)
        self._analysis_worker.finished.connect(self._on_analysis_finished)
        self._analysis_worker.failed.connect(self._on_analysis_failed)
        self._analysis_worker.finished.connect(self._analysis_thread.quit)
        self._analysis_worker.failed.connect(self._analysis_thread.quit)
        self._analysis_thread.finished.connect(self._cleanup_analysis)
        self._analysis_thread.start()

    def _on_analysis_progress(self, message: str) -> None:
        self.statusBar().showMessage(message)

    def _on_analysis_finished(self, _run_state) -> None:
        self._switch_page(self.nav.currentRow())
        self.statusBar().showMessage("Analysis completed")

    def _on_analysis_failed(self, message: str) -> None:
        self.statusBar().showMessage("Analysis failed")
        QMessageBox.critical(self, "Analysis failed", message)

    def _cleanup_analysis(self) -> None:
        self.open_action.setEnabled(True)
        self.run_action.setEnabled(self.project_service.current_project is not None)
        self.progress_bar.setVisible(False)
        if self._analysis_worker is not None:
            self._analysis_worker.deleteLater()
        if self._analysis_thread is not None:
            self._analysis_thread.deleteLater()
        self._analysis_worker = None
        self._analysis_thread = None

    def _advance_after_project_confirmation(self) -> None:
        self.nav.setCurrentRow(1)

    def _advance_after_pairing_confirmation(self) -> None:
        self.nav.setCurrentRow(2)

    def _advance_after_comparison_sets_confirmation(self) -> None:
        self.nav.setCurrentRow(3)

    def _advance_after_visual_groups_confirmation(self) -> None:
        self.nav.setCurrentRow(self._named_nav_rows.get("5. Per-comparison Analysis", 4))

    def _open_named_page(self, target: str) -> None:
        row = self._named_nav_rows.get(target)
        if row is not None:
            self.nav.setCurrentRow(row)

def launch() -> None:
    app = QApplication(sys.argv)
    app.setStyleSheet(KATANA_STYLESHEET)
    icon_path = bundled_assets_root() / "katana_icon.png"
    if icon_path.exists():
        app.setWindowIcon(QIcon(str(icon_path)))
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
