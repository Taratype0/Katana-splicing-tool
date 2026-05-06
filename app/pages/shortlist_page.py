from __future__ import annotations

import pandas as pd
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QHBoxLayout,
    QComboBox,
    QLabel,
    QLineEdit,
    QPushButton,
    QSpinBox,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from src.services.project_service import ProjectService


class ShortlistPage(QWidget):
    def __init__(self, project_service: ProjectService) -> None:
        super().__init__()
        self.project_service = project_service
        self.label = QLabel("Shortlist")
        self.explanation = QLabel(
            "Shortlist is the candidate-selection layer after Candidates. Genes kept here move forward into Cards and Sashimi unless they are later excluded in Blacklist."
        )
        self.explanation.setWordWrap(True)
        self.top_n_dp = QSpinBox()
        self.top_n_dp.setRange(1, 100)
        self.top_n_dp.setValue(12)
        self.set_tabs = QTabWidget()
        self.set_tabs.currentChanged.connect(self.refresh)
        self.set_filter = QComboBox()
        self.set_filter.currentIndexChanged.connect(self.refresh)
        self.generate_button = QPushButton("Generate Shortlist")
        self.generate_button.clicked.connect(self._generate_shortlist)
        self.custom_input = QLineEdit()
        self.custom_input.setPlaceholderText("Add one shortlist gene")
        self.add_custom_button = QPushButton("Add Gene")
        self.add_custom_button.clicked.connect(self._add_custom_gene)
        self.bulk_input = QTextEdit()
        self.bulk_input.setPlaceholderText("Paste one gene per line or TSV/CSV cell block")
        self.bulk_add_button = QPushButton("Add Pasted Genes")
        self.bulk_add_button.clicked.connect(self._add_bulk_genes)
        self.remove_button = QPushButton("Remove Selected Shortlist Genes")
        self.remove_button.clicked.connect(self._remove_selected_shortlist)
        self.table = QTableWidget(0, 0)
        self.table.setSortingEnabled(True)
        self.shortlist_table = QTableWidget(0, 1)
        self.shortlist_table.setHorizontalHeaderLabels(["Manual shortlist genes"])
        self._frame = pd.DataFrame()

        layout = QVBoxLayout(self)
        layout.addWidget(self.label)
        layout.addWidget(self.explanation)
        layout.addWidget(QLabel("Top N per comparison set"))
        layout.addWidget(self.top_n_dp)
        layout.addWidget(self.set_tabs)
        layout.addWidget(QLabel("Comparison set"))
        layout.addWidget(self.set_filter)
        layout.addWidget(self.generate_button)
        manual_row = QHBoxLayout()
        manual_row.addWidget(self.custom_input)
        manual_row.addWidget(self.add_custom_button)
        layout.addLayout(manual_row)
        layout.addWidget(self.bulk_input)
        layout.addWidget(self.bulk_add_button)
        layout.addWidget(self.remove_button)
        layout.addWidget(self.shortlist_table)
        layout.addWidget(self.table)

    def refresh(self) -> None:
        run_state = self.project_service.run_state
        frame = run_state.results.cards_shortlist if run_state.status == "completed" else pd.DataFrame()
        self._frame = frame
        if frame.empty:
            self.table.setRowCount(0)
            self.table.setColumnCount(0)
            self.set_filter.blockSignals(True)
            self.set_filter.clear()
            self.set_filter.addItem("All sets")
            self.set_filter.blockSignals(False)
            self.set_tabs.clear()
            self.label.setText("Shortlist (run analysis first)")
            self._refresh_manual_shortlist_table()
            return
        self._populate_set_filter(frame)
        self._populate_set_tabs(frame)
        selected_set = self._selected_set()
        if selected_set:
            if "pair_name" in frame.columns:
                frame = frame.loc[frame["pair_name"].astype(str) == str(selected_set)].copy()
            elif "display_group" in frame.columns:
                frame = frame.loc[frame["display_group"].astype(str) == str(selected_set)].copy()
        columns = [
            "pair_name",
            "display_group",
            "geneSymbol",
            "GeneID",
            "event_type",
            "class_interpretation",
            "comparison_A",
            "comparison_B",
            "rank_score",
        ]
        available = [column for column in columns if column in frame.columns]
        self.table.setColumnCount(len(available))
        self.table.setHorizontalHeaderLabels(available)
        self.table.setRowCount(len(frame))
        for row_idx, (_, row) in enumerate(frame[available].iterrows()):
            for col_idx, column in enumerate(available):
                item = QTableWidgetItem("" if pd.isna(row[column]) else str(row[column]))
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(row_idx, col_idx, item)
        self.label.setText(f"Shortlist ({len(frame)} rows)")
        self._refresh_manual_shortlist_table()

    def _generate_shortlist(self) -> None:
        run_state = self.project_service.run_state
        ranked = run_state.results.ranked_candidates
        if ranked.empty:
            self.label.setText("Shortlist (no ranked candidates; run analysis first)")
            return
        shortlist_dp, shortlist_ko = self.project_service.shortlist_service.split_shortlists(ranked, top_n_dp=self.top_n_dp.value(), top_n_ko=self.top_n_dp.value())
        run_state.results.shortlist_dp = shortlist_dp
        run_state.results.shortlist_ko = shortlist_ko
        run_state.results.cards_shortlist = self.project_service.shortlist_service.build_all_shortlists(
            ranked,
            rank_limit=self.top_n_dp.value(),
        )
        if self.project_service.current_project and self.project_service.current_project.output_root:
            self.project_service._write_step_output(
                "04b_shortlist",
                self.project_service._shortlist_output_tables(),
            )
            self.project_service._save_shortlist_plot()
        self.refresh()

    def _add_custom_gene(self) -> None:
        value = self.custom_input.text().strip()
        if not value:
            return
        self.project_service.add_shortlist_gene(value)
        self.custom_input.clear()
        self.refresh()

    def _add_bulk_genes(self) -> None:
        raw = self.bulk_input.toPlainText().strip()
        if not raw:
            return
        normalized = raw.replace(",", "\n").replace("\t", "\n")
        for gene in normalized.splitlines():
            value = gene.strip()
            if value:
                self.project_service.add_shortlist_gene(value)
        self.bulk_input.clear()
        self.refresh()

    def _remove_selected_shortlist(self) -> None:
        for row in sorted({item.row() for item in self.shortlist_table.selectedItems()}, reverse=True):
            gene_item = self.shortlist_table.item(row, 0)
            if gene_item is not None:
                self.project_service.remove_shortlist_gene(gene_item.text().strip())
        self.refresh()

    def _populate_set_filter(self, frame: pd.DataFrame) -> None:
        if "pair_name" in frame.columns:
            values = frame["pair_name"].dropna().astype(str).drop_duplicates().tolist()
        elif "display_group" in frame.columns:
            values = frame["display_group"].dropna().astype(str).drop_duplicates().tolist()
        else:
            values = []
        current = self.set_filter.currentData()
        self.set_filter.blockSignals(True)
        self.set_filter.clear()
        self.set_filter.addItem("All sets", None)
        for value in values:
            self.set_filter.addItem(value, value)
        if current is not None:
            index = self.set_filter.findData(current)
            if index >= 0:
                self.set_filter.setCurrentIndex(index)
        self.set_filter.blockSignals(False)

    def _populate_set_tabs(self, frame: pd.DataFrame) -> None:
        if "pair_name" in frame.columns:
            values = frame["pair_name"].dropna().astype(str).drop_duplicates().tolist()
        elif "display_group" in frame.columns:
            values = frame["display_group"].dropna().astype(str).drop_duplicates().tolist()
        else:
            values = []
        current = self._selected_set()
        self.set_tabs.blockSignals(True)
        self.set_tabs.clear()
        self.set_tabs.addTab(QWidget(), "All sets")
        self.set_tabs.setTabToolTip(0, "__all__")
        for index, value in enumerate(values, start=1):
            self.set_tabs.addTab(QWidget(), value)
            self.set_tabs.setTabToolTip(index, value)
        if current:
            for index in range(self.set_tabs.count()):
                if self.set_tabs.tabToolTip(index) == current:
                    self.set_tabs.setCurrentIndex(index)
                    break
        self.set_tabs.blockSignals(False)

    def _selected_set(self):
        if self.set_tabs.count() > 0:
            tooltip = self.set_tabs.tabToolTip(self.set_tabs.currentIndex())
            if tooltip and tooltip != "__all__":
                return tooltip
        return self.set_filter.currentData()

    def _refresh_manual_shortlist_table(self) -> None:
        genes = self.project_service.shortlist_genes_for_display()
        self.shortlist_table.setRowCount(len(genes))
        for row, gene in enumerate(genes):
            self.shortlist_table.setItem(row, 0, QTableWidgetItem(gene))
