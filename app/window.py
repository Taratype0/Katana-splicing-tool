from __future__ import annotations

import sys
from pathlib import Path

from PySide6.QtCore import QObject, QThread, Signal, Qt
from PySide6.QtGui import QIcon
from PySide6.QtWidgets import (
    QApplication,
    QFileDialog,
    QFrame,
    QListWidget,
    QListWidgetItem,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QSplitter,
    QStackedWidget,
    QStatusBar,
    QToolBar,
    QVBoxLayout,
    QWidget,
    QProgressBar,
)

from app.pages.comparison_page import ComparisonPage
from app.pages.comparison_sets_page import ComparisonSetsPage
from app.pages.candidate_gene_page import CandidateGenePage
from app.pages.candidate_selection_page import CandidateSelectionPage
from app.pages.cards_page import CardsPage
from app.pages.blacklist_page import BlacklistPage
from app.pages.input_check_page import InputCheckPage
from app.pages.input_paths_dialog import InputPathsDialog
from app.pages.jutils_page import JutilsPage
from app.pages.project_page import ProjectPage
from app.pages.results_page import ResultsPage
from app.pages.section_navigation_page import SectionNavigationPage
from app.pages.sashimi_page import SashimiPage
from app.pages.shortlist_page import ShortlistPage
from app.pages.samples_page import SamplesPage
from app.pages.settings_page import SettingsPage
from app.pages.thresholds_page import ThresholdsPage
from app.pages.visualization_groups_page import VisualizationGroupsPage
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
        self.pages = QStackedWidget()
        self.nav = QListWidget()
        self.nav.setMinimumWidth(320)
        self.nav.setMaximumWidth(320)
        self.nav.currentRowChanged.connect(self._switch_page)
        self._page_routes: list[tuple[int, QWidget, tuple[str, str | None] | None]] = []
        self._named_nav_rows: dict[str, int] = {}

        self.project_page = ProjectPage(self.project_service, self._advance_after_project_confirmation)
        self.comparison_page = ComparisonPage(self.project_service, self._advance_after_pairing_confirmation)
        self.comparison_sets_page = ComparisonSetsPage(self.project_service, self._advance_after_comparison_sets_confirmation)
        self.visualization_groups_page = VisualizationGroupsPage(self.project_service, self._advance_after_visual_groups_confirmation)
        self.per_comparison_page = SectionNavigationPage(
            "Per-comparison Analysis",
            "Run analysis inside each comparison first. Start with overview and transcript+splicing integration, then move into candidate screening, cross-comparison candidate comparison, sashimi and Jutils follow-up.",
            [
                ("Overview", "Comparison-level overview plots before candidate calling.", "5.1 Overview"),
                ("Transcript + Splicing Integration", "Direction-aware evidence layers inside each comparison.", "5.2 Transcript + Splicing Integration"),
                ("Candidate Gene Screening", "Generate ranked gene-level candidates, then shortlist/blacklist them.", "5.3 Candidate Gene Screening"),
                ("Cross-comparison Candidate Comparison", "Compare candidate genes only after per-comparison screening.", "5.4 Cross-comparison Candidate Comparison"),
                ("Sashimi", "Validate shortlisted candidate events with sashimi plots.", "5.5 Sashimi"),
                ("Jutils", "Heatmap/PCA/output browser for jutils follow-up.", "5.6 Jutils"),
            ],
            self._open_named_page,
        )
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
            ],
            self._open_named_page,
        )
        self.cross_comparison_page = SectionNavigationPage(
            "Cross-comparison Candidate Comparison",
            "Compare candidate genes only after per-comparison screening. Use this layer to inspect shared/specific/gained/lost candidates and direction reversal across comparisons.",
            [
                ("Candidate Matrix", "Matrix of candidate membership across comparisons.", "5.4.1 Candidate Matrix"),
                ("Shared / Specific / Gained / Lost", "Pattern summary across candidate sets.", "5.4.2 Shared / Specific / Gained / Lost"),
                ("Direction Reversal", "Inspect genes whose standardized directions reverse across comparisons.", "5.4.3 Direction Reversal"),
                ("Group Comparison", "Program-level comparison between selected candidate groups.", "5.4.4 Group Comparison"),
                ("Direction vs Strength", "Summarize reversal vs effect-size shifts across comparisons.", "5.4.5 Direction vs Strength"),
                ("Cross-comparison significant AS event dPSI pattern analysis", "Compare shared significant AS events across comparisons using standardized dPSI differences only.", "5.4.6 Cross-comparison significant AS event dPSI pattern analysis"),
            ],
            self._open_named_page,
        )
        self.sashimi_hub_page = SectionNavigationPage(
            "Sashimi",
            "Sashimi is the candidate-event validation layer after shortlist/blacklist. Start with input checks, then run and review preview or failed jobs.",
            [
                ("Input Check", "Confirm candidate events, BAM lists and manifest inputs.", "5.5.1 Input Check"),
                ("Run Sashimi", "Open the sashimi runner and launch rmats2sashimi.", "5.5.2 Run Sashimi"),
                ("Sashimi Preview", "Preview shortlisted candidate events before/after running.", "5.5.3 Sashimi Preview"),
                ("Failed Sashimi Jobs", "Review failed sashimi runs and diagnostics.", "5.5.4 Failed Sashimi Jobs"),
            ],
            self._open_named_page,
        )
        self.jutils_hub_page = SectionNavigationPage(
            "Jutils",
            "Jutils follow-up is grouped in one place. Open heatmap, PCA or the output browser from here.",
            [
                ("Jutils Heatmap", "Preview or run heatmap-oriented jutils output.", "5.6.1 Jutils Heatmap"),
                ("Jutils PCA", "Preview or run PCA-oriented jutils output.", "5.6.2 Jutils PCA"),
                ("Jutils Output Browser", "Inspect generated files, paths and logs.", "5.6.3 Jutils Output Browser"),
            ],
            self._open_named_page,
        )
        self.landscape_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="landscape", fixed_subbranch="overview")
        self.landscape_bar_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="landscape", fixed_subbranch="bar")
        self.landscape_pie_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="landscape", fixed_subbranch="pie")
        self.landscape_heatmap_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="landscape", fixed_subbranch="heatmap")
        self.landscape_pca_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="landscape", fixed_subbranch="pca")
        self.group_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="group")
        self.group_matrix_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="group", fixed_subbranch="matrix")
        self.group_shared_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="group", fixed_subbranch="shared")
        self.group_reversal_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="group", fixed_subbranch="reversal")
        self.group_comparison_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="group", fixed_subbranch="comparison")
        self.group_direction_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="group", fixed_subbranch="direction")
        self.group_pattern_page = ResultsPage(self.project_service, self.run_analysis, fixed_branch="group", fixed_subbranch="pattern")
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
        self.sashimi_failed_page = SashimiPage(self.project_service, view_mode="failed")
        self.samples_page = SamplesPage(self.project_service)
        self.input_check_page = InputCheckPage(self.project_service)
        self.jutils_heatmap_page = JutilsPage(self.project_service, view_mode="heatmap")
        self.jutils_pca_page = JutilsPage(self.project_service, view_mode="pca")
        self.jutils_page = JutilsPage(self.project_service, view_mode="browser")
        self.thresholds_page = ThresholdsPage(self.project_service)
        self.settings_page = SettingsPage(self.project_service, self.reload_current_project)

        self._add_page("1. Project", self.project_page)
        self._add_page("2. Pairing", self.comparison_page)
        self._add_page("3. Comparison Sets", self.comparison_sets_page)
        self._add_page("4. Visual Groups", self.visualization_groups_page)
        self._add_page("5. Per-comparison Analysis", self.per_comparison_page)
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
        self._add_page("5.4 Cross-comparison Candidate Comparison", self.cross_comparison_page)
        self._add_page("5.4.1 Candidate Matrix", self.group_matrix_page)
        self._add_page("5.4.2 Shared / Specific / Gained / Lost", self.group_shared_page)
        self._add_page("5.4.3 Direction Reversal", self.group_reversal_page)
        self._add_page("5.4.4 Group Comparison", self.group_comparison_page)
        self._add_page("5.4.5 Direction vs Strength", self.group_direction_page)
        self._add_page("5.4.6 Cross-comparison significant AS event dPSI pattern analysis", self.group_pattern_page)
        self._add_page("5.5 Sashimi", self.sashimi_hub_page)
        self._add_page("5.5.1 Input Check", self.input_check_page)
        self._add_page("5.5.2 Run Sashimi", self.sashimi_page)
        self._add_page("5.5.3 Sashimi Preview", self.sashimi_preview_page)
        self._add_page("5.5.4 Failed Sashimi Jobs", self.sashimi_failed_page)
        self._add_page("5.6 Jutils", self.jutils_hub_page)
        self._add_page("5.6.1 Jutils Heatmap", self.jutils_heatmap_page)
        self._add_page("5.6.2 Jutils PCA", self.jutils_pca_page)
        self._add_page("5.6.3 Jutils Output Browser", self.jutils_page)
        self._add_page("6. Samples", self.samples_page)
        self._add_page("7. Input Check", self.input_check_page)
        self._add_page("8. Settings", self.settings_page)
        self.nav.setCurrentRow(0)

        left_panel = QWidget()
        left_panel.setMinimumWidth(320)
        left_panel.setMaximumWidth(340)
        left_layout = QVBoxLayout(left_panel)
        left_layout.setContentsMargins(0, 0, 0, 0)
        left_layout.addWidget(self.nav, 1)
        divider = QFrame()
        divider.setFrameShape(QFrame.Shape.HLine)
        divider.setFrameShadow(QFrame.Shadow.Sunken)
        left_layout.addWidget(divider)
        left_layout.addWidget(self.thresholds_page, 0)

        splitter = QSplitter()
        splitter.addWidget(left_panel)
        splitter.addWidget(self.pages)
        splitter.setStretchFactor(1, 1)
        splitter.setChildrenCollapsible(False)
        splitter.setSizes([360, 1000])
        self.setCentralWidget(splitter)

        toolbar = QToolBar("Main")
        self.addToolBar(toolbar)

        self.open_action = QPushButton("Open Project")
        self.open_action.clicked.connect(self.open_project)
        toolbar.addWidget(self.open_action)

        self.run_action = QPushButton("Run Selected")
        self.run_action.clicked.connect(self.run_analysis)
        toolbar.addWidget(self.run_action)

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
        self._update_threshold_visibility(0)

    def _add_page(self, name: str, widget: QWidget, analysis_route: tuple[str, str | None] | None = None) -> None:
        page_index = self.pages.indexOf(widget)
        if page_index < 0:
            self.pages.addWidget(widget)
            page_index = self.pages.indexOf(widget)
        self._page_routes.append((page_index, widget, analysis_route))
        item = QListWidgetItem(self._nav_display_label(name), self.nav)
        item.setData(Qt.ItemDataRole.UserRole, name)
        self._named_nav_rows[name] = len(self._page_routes) - 1

    def _nav_display_label(self, name: str) -> str:
        prefix = name.split(" ", 1)[0]
        depth = max(prefix.count(".") - 1, 0)
        if prefix == "5.":
            depth = 0
        return f"{'    ' * depth}{name}"

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
            page_index, page, analysis_route = self._page_routes[index]
            self.pages.setCurrentIndex(page_index)
            self._update_threshold_visibility(index)
            if analysis_route is not None:
                branch, subbranch = analysis_route
                navigate_to = getattr(page, "navigate_to", None)
                if callable(navigate_to):
                    navigate_to(branch, subbranch)
            refresh = getattr(page, "refresh", None)
            if callable(refresh):
                refresh()

    def _update_threshold_visibility(self, index: int) -> None:
        if not (0 <= index < len(self._page_routes)):
            self.thresholds_page.setVisible(False)
            return
        item = self.nav.item(index)
        label = item.data(Qt.ItemDataRole.UserRole) if item is not None else ""
        visible = label.startswith("5.")
        self.thresholds_page.setVisible(visible)

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
        self.statusBar().showMessage("Project loaded")
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
        self.statusBar().showMessage("Running selected analysis modules...")
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
    icon_path = bundled_assets_root() / "katana_icon.png"
    if icon_path.exists():
        app.setWindowIcon(QIcon(str(icon_path)))
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
