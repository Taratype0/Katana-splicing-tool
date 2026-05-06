from __future__ import annotations

import pandas as pd
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QComboBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from src.services.project_service import ProjectService


class BlacklistPage(QWidget):
    def __init__(self, project_service: ProjectService) -> None:
        super().__init__()
        self.project_service = project_service
        self.label = QLabel("Blacklist")
        self.explanation = QLabel(
            "Blacklist is the exclusion layer after Shortlist. Genes placed here will be excluded from downstream Cards and Sashimi."
        )
        self.explanation.setWordWrap(True)
        self.comparison_tabs = QTabWidget()
        self.comparison_tabs.currentChanged.connect(self.refresh)
        self.source = QComboBox()
        self.source.addItem("Shortlist genes", "shortlist")
        self.source.addItem("All candidate genes", "all")
        self.source.currentIndexChanged.connect(self.refresh)
        self.search = QLineEdit()
        self.search.setPlaceholderText("Search genes...")
        self.search.textChanged.connect(self.refresh)
        self.custom_input = QLineEdit()
        self.custom_input.setPlaceholderText("Add custom gene symbol or gene ID")
        self.add_custom_button = QPushButton("Add")
        self.add_custom_button.clicked.connect(self._add_custom)
        self.add_selected_button = QPushButton("Blacklist Selected")
        self.add_selected_button.clicked.connect(self._add_selected)
        self.remove_selected_button = QPushButton("Remove Selected")
        self.remove_selected_button.clicked.connect(self._remove_selected)

        self.available_list = QListWidget()
        self.available_list.setSelectionMode(QListWidget.SelectionMode.ExtendedSelection)
        self.available_list.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.available_list.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.blacklist_list = QListWidget()
        self.blacklist_list.setSelectionMode(QListWidget.SelectionMode.ExtendedSelection)
        self.blacklist_list.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.blacklist_list.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)

        top = QHBoxLayout()
        top.addWidget(self.source)
        top.addWidget(self.search)
        top.addWidget(self.add_selected_button)
        top.addWidget(self.remove_selected_button)

        custom = QHBoxLayout()
        custom.addWidget(self.custom_input)
        custom.addWidget(self.add_custom_button)

        lists = QHBoxLayout()
        left = QVBoxLayout()
        left.addWidget(QLabel("Available genes"))
        left.addWidget(self.available_list)
        right = QVBoxLayout()
        right.addWidget(QLabel("Blacklisted genes"))
        right.addWidget(self.blacklist_list)
        lists.addLayout(left, 1)
        lists.addLayout(right, 1)

        layout = QVBoxLayout(self)
        layout.addWidget(self.label)
        layout.addWidget(self.explanation)
        layout.addWidget(self.comparison_tabs)
        layout.addLayout(top)
        layout.addLayout(custom)
        layout.addLayout(lists)

    def refresh(self) -> None:
        project = self.project_service.current_project
        if project is None:
            self.comparison_tabs.clear()
            self.available_list.clear()
            self.blacklist_list.clear()
            return

        available = self._candidate_catalog()
        self._populate_comparison_tabs(available)
        blacklist = set(self.project_service.blacklist_genes_for_display())
        selected_comparison = self._selected_comparison()
        if selected_comparison and selected_comparison != "__all__" and "comparison_id" in available.columns:
            available = available.loc[available["comparison_id"].astype(str) == str(selected_comparison)].copy()
        query = self.search.text().strip().lower()
        if query and not available.empty:
            available = available[
                available["gene_key"].astype(str).str.lower().str.contains(query, na=False)
                | available["gene_label"].astype(str).str.lower().str.contains(query, na=False)
            ].copy()

        self.available_list.clear()
        for _, row in available.iterrows():
            gene_key = str(row["gene_key"])
            if gene_key in blacklist:
                continue
            item = QListWidgetItem(f"{row['gene_label']} ({gene_key})")
            item.setData(Qt.ItemDataRole.UserRole, gene_key)
            self.available_list.addItem(item)

        self.blacklist_list.clear()
        for gene in self.project_service.blacklist_genes_for_display():
            item = QListWidgetItem(gene)
            item.setData(Qt.ItemDataRole.UserRole, gene)
            self.blacklist_list.addItem(item)

        self.label.setText(f"Blacklist ({self.blacklist_list.count()} gene(s))")

    def _candidate_catalog(self) -> pd.DataFrame:
        run_state = self.project_service.run_state
        source = self.source.currentData()
        if source == "shortlist" and run_state.status == "completed" and not run_state.results.cards_shortlist.empty:
            frame = run_state.results.cards_shortlist.copy()
            frame["gene_key"] = frame["GeneID"].fillna(frame.get("gene_id")).astype(str)
            frame["gene_label"] = frame["geneSymbol"].fillna(frame["gene_key"]).astype(str)
            if "comparison_id" not in frame.columns:
                if "comparison_A" in frame.columns:
                    frame["comparison_id"] = frame["comparison_A"].fillna("").astype(str)
                else:
                    frame["comparison_id"] = ""
            return frame[["comparison_id", "gene_key", "gene_label"]].drop_duplicates().reset_index(drop=True)
        if run_state.status == "completed" and not run_state.results.tx_splicing_gene_table.empty:
            frame = run_state.results.tx_splicing_gene_table.copy()
            frame["gene_key"] = frame["gene_id"].astype(str)
            frame["gene_label"] = frame["geneSymbol"].fillna(frame["gene_key"]).astype(str)
            if "comparison_id" not in frame.columns:
                if "comparison_name" in frame.columns:
                    frame["comparison_id"] = frame["comparison_name"].fillna("").astype(str)
                else:
                    frame["comparison_id"] = ""
            return frame[["comparison_id", "gene_key", "gene_label"]].drop_duplicates().reset_index(drop=True)
        return pd.DataFrame(columns=["comparison_id", "gene_key", "gene_label"])

    def _populate_comparison_tabs(self, frame: pd.DataFrame) -> None:
        current = self._selected_comparison()
        values = []
        if "comparison_id" in frame.columns:
            values = sorted(frame["comparison_id"].dropna().astype(str).unique().tolist())
        self.comparison_tabs.blockSignals(True)
        self.comparison_tabs.clear()
        self.comparison_tabs.addTab(QWidget(), "All comparisons")
        self.comparison_tabs.setTabToolTip(0, "__all__")
        for index, value in enumerate(values, start=1):
            self.comparison_tabs.addTab(QWidget(), value)
            self.comparison_tabs.setTabToolTip(index, value)
        if current:
            for index in range(self.comparison_tabs.count()):
                if self.comparison_tabs.tabToolTip(index) == current:
                    self.comparison_tabs.setCurrentIndex(index)
                    break
        self.comparison_tabs.blockSignals(False)

    def _selected_comparison(self):
        if self.comparison_tabs.count() > 0:
            tooltip = self.comparison_tabs.tabToolTip(self.comparison_tabs.currentIndex())
            if tooltip:
                return tooltip
        return "__all__"

    def _add_selected(self) -> None:
        for item in self.available_list.selectedItems():
            self.project_service.add_blacklist_gene(str(item.data(Qt.ItemDataRole.UserRole)))
        self.refresh()

    def _remove_selected(self) -> None:
        for item in self.blacklist_list.selectedItems():
            self.project_service.remove_blacklist_gene(str(item.data(Qt.ItemDataRole.UserRole)))
        self.refresh()

    def _add_custom(self) -> None:
        value = self.custom_input.text().strip()
        if not value:
            return
        self.project_service.add_blacklist_gene(value)
        self.custom_input.clear()
        self.refresh()
