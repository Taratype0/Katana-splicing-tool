from __future__ import annotations

import os
from dataclasses import dataclass

from PySide6.QtWidgets import (
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)

from app.widgets.flow_layout import FlowLayout
from app.widgets.page_header import PageHeader
from src.services.project_service import ProjectService


@dataclass(frozen=True)
class _StageSpec:
    number: str
    title: str
    outcome: str
    handoff: str
    routes: tuple[tuple[str, str], ...]
    module_keys: tuple[str, ...] = ()


_STAGES = (
    _StageSpec(
        "A",
        "Explore the experiment",
        "Establish AS composition, DEG behavior, sample-level PSI structure, and replicate quality before prioritization.",
        "rMATS + DEG + sample PSI → trusted comparison-level context",
        (
            ("Overview", "5.1 Overview"),
            ("DEG", "5.7 DEG Explorer"),
            ("PSI patterns", "5.6 Native Splicing Pattern Explorer"),
        ),
        ("splicing_landscape", "deg_explorer", "jutils_heatmap", "jutils_pca"),
    ),
    _StageSpec(
        "B",
        "Integrate transcript and splicing",
        "Combine standardized DEG and significant rMATS evidence inside each comparison without re-guessing direction.",
        "trusted comparison context → gene/event integration table",
        (("Open integration", "5.2 Transcript + Splicing Integration"),),
        ("tx_splicing_integration",),
    ),
    _StageSpec(
        "C",
        "Prioritize candidate genes",
        "Rank unique genes, curate a comparison-specific selection, then inspect evidence heatmaps and event-centered cards.",
        "integration table → ranking → selected candidate genes",
        (
            ("Ranking", "5.3.1 Candidate Ranking Table"),
            ("Selection", "5.3.2 Candidate Gene Selection"),
            ("Evidence", "5.3.3 Evidence Heatmap"),
            ("Cards", "5.3.4 Candidate Cards + Event Follow-up"),
        ),
        ("candidate_ranking", "candidate_selection", "candidate_heatmap", "candidate_followup"),
    ),
    _StageSpec(
        "D",
        "Compare biological effects",
        "Use explicit comparison sets to distinguish gene-level candidate membership from strict AS event overlap and interaction evidence.",
        "per-comparison candidates/events → A-versus-B evidence",
        (
            ("Gene comparison", "5.4.1 Candidate Matrix"),
            ("AS interaction", "5.4.3 Cross-comparison AS interaction modeling"),
            ("Event overlap", "5.4.4 Venn / Overlap"),
        ),
        ("cross_comparison_candidates", "cross_comparison_as_patterns", "cross_comparison_as_interactions"),
    ),
    _StageSpec(
        "E",
        "Validate selected events",
        "Inspect read support and generate Native Sashimi only for user-selected events.",
        "selected candidates + AS events + BAM/GTF → event validation plots",
        (
            ("Input check", "5.5.1 Input Check"),
            ("Run Sashimi", "5.5.2 Run Sashimi"),
            ("Preview", "5.5.3 Sashimi Preview"),
        ),
        ("sashimi_manifest", "sashimi_plot"),
    ),
    _StageSpec(
        "F",
        "Interpret mechanism and function",
        "Add functional enrichment, RBP support, and isoform support to selected candidates without presenting association as causality.",
        "selected candidates → functional and mechanistic support",
        (
            ("GO / KEGG", "5.8 GO / KEGG Enrichment"),
            ("RBP support", "7. RBP"),
            ("Isoform support", "8. Isoform"),
        ),
        ("functional_enrichment", "rbp_support", "rbp_cross_support", "isoform_support", "isoform_cross_support"),
    ),
    _StageSpec(
        "G",
        "Discover co-splicing programs",
        "Model coordinated event programs across sufficiently sized cohorts, separately from candidate screening.",
        "sample-level PSI matrix + metadata → modules, preservation, and program reports",
        (("Open CoSplice", "9. CoSplice Network"),),
    ),
)


class WorkflowDashboardPage(QWidget):
    """Read-only workflow map that connects existing independent analysis modules."""

    def __init__(self, project_service: ProjectService, on_open_section=None) -> None:
        super().__init__()
        self.project_service = project_service
        self.on_open_section = on_open_section
        self._stage_status_labels: dict[str, QLabel] = {}

        content = QWidget()
        content_layout = QVBoxLayout(content)
        content_layout.setContentsMargins(2, 2, 10, 18)
        content_layout.setSpacing(16)
        content_layout.addWidget(
            PageHeader(
                "Analysis workspace",
                "Follow one connected evidence path. Every page reads existing upstream results; analyses still run independently only when you click Run.",
                "05 · Discover → prioritize → validate",
            )
        )

        self.context = QLabel()
        self.context.setObjectName("workflowContext")
        self.context.setWordWrap(True)
        content_layout.addWidget(self.context)
        report_actions = FlowLayout(spacing=8)
        build_report = QPushButton("Build cached project report")
        build_report.setProperty("class", "primary")
        build_report.clicked.connect(self._build_project_report)
        open_report = QPushButton("Open latest project report")
        open_report.clicked.connect(self._open_project_report)
        report_actions.addWidget(build_report)
        report_actions.addWidget(open_report)
        content_layout.addLayout(report_actions)

        grid = QGridLayout()
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setHorizontalSpacing(16)
        grid.setVerticalSpacing(16)
        for index, stage in enumerate(_STAGES):
            grid.addWidget(self._build_stage(stage), index // 2, index % 2)
        grid.setColumnStretch(0, 1)
        grid.setColumnStretch(1, 1)
        content_layout.addLayout(grid)
        content_layout.addStretch(1)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        scroll.setWidget(content)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(scroll)
        self.refresh()

    def _build_stage(self, stage: _StageSpec) -> QFrame:
        card = QFrame()
        card.setObjectName("workflowStageCard")
        layout = QVBoxLayout(card)
        layout.setContentsMargins(18, 16, 18, 16)
        layout.setSpacing(8)

        heading_row = QHBoxLayout()
        number = QLabel(stage.number)
        number.setObjectName("workflowStageNumber")
        title = QLabel(stage.title)
        title.setObjectName("cardTitle")
        title.setWordWrap(True)
        status = QLabel("Not run")
        status.setObjectName("statusChip")
        status.setProperty("state", "not_run")
        self._stage_status_labels[stage.number] = status
        heading_row.addWidget(number)
        heading_row.addWidget(title, 1)
        heading_row.addWidget(status)

        outcome = QLabel(stage.outcome)
        outcome.setObjectName("cardBody")
        outcome.setWordWrap(True)
        handoff = QLabel(stage.handoff)
        handoff.setObjectName("handoffLabel")
        handoff.setWordWrap(True)

        buttons = FlowLayout(spacing=7)
        for button_index, (label, route) in enumerate(stage.routes):
            button = QPushButton(label)
            if button_index == 0:
                button.setProperty("class", "primary")
            button.clicked.connect(lambda _checked=False, target=route: self._open(target))
            buttons.addWidget(button)
        layout.addLayout(heading_row)
        layout.addWidget(outcome)
        layout.addWidget(handoff)
        layout.addStretch(1)
        layout.addLayout(buttons)
        return card

    def refresh(self) -> None:
        project = self.project_service.current_project
        if project is None:
            self.context.setText("No project loaded. Start with Project, then confirm Pairing and Comparison Sets.")
        else:
            selected = len(getattr(project, "selected_comparison_ids", []) or [])
            pairs = len([pair for pair in getattr(project, "comparison_pairs", []) if getattr(pair, "enabled", True)])
            selections = sum(
                bool(genes)
                for genes in (getattr(project, "candidate_selection_genes", {}) or {}).values()
            )
            self.context.setText(
                f"Project: {project.project_name}   •   Active comparisons: {selected}   •   "
                f"Comparison sets: {pairs}   •   Candidate selections: {selections}\n"
                "Shared context: canonical directions and thresholds flow forward; downstream pages never rerun upstream analyses automatically."
            )

        for stage in _STAGES:
            statuses = [
                str(self.project_service.module_state(key).get("status", "not_run"))
                for key in stage.module_keys
            ]
            state, label = self._aggregate_status(statuses)
            status_label = self._stage_status_labels[stage.number]
            status_label.setText(label)
            status_label.setProperty("state", state)
            status_label.style().unpolish(status_label)
            status_label.style().polish(status_label)

    @staticmethod
    def _aggregate_status(statuses: list[str]) -> tuple[str, str]:
        if not statuses:
            return "on_demand", "On demand"
        if "running" in statuses:
            return "running", "Running"
        if "failed" in statuses:
            return "failed", "Needs attention"
        completed = sum(status in {"finished", "cache_hit"} for status in statuses)
        if completed == len(statuses):
            return "finished", "Complete"
        if completed:
            return "partial", f"{completed}/{len(statuses)} ready"
        return "not_run", "Not run"

    def _open(self, target: str) -> None:
        if callable(self.on_open_section):
            self.on_open_section(target)

    def _build_project_report(self) -> None:
        try:
            outputs = self.project_service.build_project_report()
        except Exception as exc:
            QMessageBox.critical(self, "Project report failed", str(exc))
            return
        report = outputs.get("report")
        self.context.setText(f"Cached project report generated without running analyses:\n{report}")
        if report and report.exists():
            os.startfile(str(report))

    def _open_project_report(self) -> None:
        project = self.project_service.current_project
        if project is None:
            QMessageBox.information(self, "No project", "Load a project first.")
            return
        try:
            output_root = self.project_service.require_output_root(create=False)
        except RuntimeError as exc:
            QMessageBox.information(self, "Output directory required", str(exc))
            return
        report = output_root / "00_project_report" / "KATANA_PROJECT_REPORT.html"
        if not report.exists():
            QMessageBox.information(
                self,
                "No cached report",
                "No cached project report exists. Click 'Build cached project report' first.",
            )
            return
        os.startfile(str(report))
