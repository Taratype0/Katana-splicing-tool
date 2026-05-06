from __future__ import annotations

import pandas as pd
from pandas.api.types import is_bool_dtype, is_numeric_dtype
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from src.services.project_service import ProjectService


class CandidateSelectionPage(QWidget):
    def __init__(self, project_service: ProjectService) -> None:
        super().__init__()
        self.project_service = project_service

        self.title = QLabel("Candidate gene selection")
        self.title.setStyleSheet("font-size: 16px; font-weight: 700;")
        self.context_label = QLabel("")
        self.context_label.setWordWrap(True)
        self.context_label.setStyleSheet("font-size: 11px; color: #C9CDD3;")

        self.comparison_tabs = QTabWidget()
        self.comparison_tabs.currentChanged.connect(self.refresh)
        self.comparison_tabs.setDocumentMode(True)
        self.comparison_tabs.setStyleSheet("QTabWidget::pane { border: 0; margin: 0; padding: 0; }")
        self.comparison_tabs.setUsesScrollButtons(True)
        self.comparison_tabs.setMaximumHeight(36)

        self.search = QLineEdit()
        self.search.setPlaceholderText("Search candidate genes for this comparison")
        self.search.textChanged.connect(self.refresh)

        self.top_n = QSpinBox()
        self.top_n.setRange(1, 500)
        self.top_n.setValue(20)

        self.load_button = QPushButton("Load Candidate Ranking")
        self.load_button.clicked.connect(self._load_candidate_ranking)
        self.reset_button = QPushButton("Reset Top 20")
        self.reset_button.clicked.connect(self._reset_top_n)

        self.gene_input = QLineEdit()
        self.gene_input.setPlaceholderText("Search/add candidate gene by symbol or gene_id")
        self.add_gene_button = QPushButton("Add gene")
        self.add_gene_button.clicked.connect(self._add_manual_gene)
        self.add_from_catalog_button = QPushButton("Add selected rows")
        self.add_from_catalog_button.clicked.connect(self._add_selected_catalog_genes)
        self.remove_gene_button = QPushButton("Remove selected genes")
        self.remove_gene_button.clicked.connect(self._remove_selected_genes)

        self.blacklist_input = QLineEdit()
        self.blacklist_input.setPlaceholderText("Add blacklist gene by symbol or gene_id")
        self.add_blacklist_button = QPushButton("Add blacklist")
        self.add_blacklist_button.clicked.connect(self._add_manual_blacklist_gene)
        self.add_blacklist_from_catalog_button = QPushButton("Blacklist selected rows")
        self.add_blacklist_from_catalog_button.clicked.connect(self._blacklist_selected_catalog_genes)
        self.remove_blacklist_button = QPushButton("Remove blacklist genes")
        self.remove_blacklist_button.clicked.connect(self._remove_selected_blacklist_genes)

        self.catalog_table = QTableWidget(0, 0)
        self.catalog_table.setSortingEnabled(True)
        self.catalog_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.catalog_table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.catalog_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)

        self.selected_table = QTableWidget(0, 0)
        self.selected_table.setSortingEnabled(True)
        self.selected_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.selected_table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.selected_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)

        self.blacklist_table = QTableWidget(0, 0)
        self.blacklist_table.setSortingEnabled(True)
        self.blacklist_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.blacklist_table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.blacklist_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)

        self.selection_status = QLabel("")
        self.selection_status.setWordWrap(True)
        self.selection_status.setStyleSheet("font-size: 11px; color: #C9CDD3;")

        controls = QHBoxLayout()
        controls.addWidget(self.load_button)
        controls.addWidget(QLabel("Default top N"))
        controls.addWidget(self.top_n)
        controls.addWidget(self.reset_button)
        controls.addSpacing(16)
        controls.addWidget(self.search, 1)

        selection_controls = QHBoxLayout()
        selection_controls.addWidget(self.gene_input, 1)
        selection_controls.addWidget(self.add_gene_button)
        selection_controls.addWidget(self.add_from_catalog_button)
        selection_controls.addWidget(self.remove_gene_button)

        blacklist_controls = QHBoxLayout()
        blacklist_controls.addWidget(self.blacklist_input, 1)
        blacklist_controls.addWidget(self.add_blacklist_button)
        blacklist_controls.addWidget(self.add_blacklist_from_catalog_button)
        blacklist_controls.addWidget(self.remove_blacklist_button)

        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)
        right_layout.setContentsMargins(0, 0, 0, 0)
        right_layout.addWidget(QLabel("Selected candidate genes"))
        right_layout.addWidget(self.selected_table, 1)
        right_layout.addLayout(selection_controls)
        right_layout.addSpacing(8)
        right_layout.addWidget(QLabel("Blacklist genes"))
        right_layout.addWidget(self.blacklist_table, 1)
        right_layout.addLayout(blacklist_controls)
        right_layout.addWidget(self.selection_status)

        content = QGridLayout()
        content.setContentsMargins(0, 0, 0, 0)
        content.setColumnStretch(0, 3)
        content.setColumnStretch(1, 2)
        content.addWidget(QLabel("Candidate ranking table"), 0, 0)
        content.addWidget(self.catalog_table, 1, 0)
        content.addWidget(right_panel, 1, 1)

        layout = QVBoxLayout(self)
        layout.addWidget(self.title)
        layout.addWidget(self.context_label)
        layout.addWidget(self.comparison_tabs)
        layout.addLayout(controls)
        layout.addLayout(content, 1)

        self._current_catalog = pd.DataFrame()
        self._current_selected = pd.DataFrame()
        self._tab_keys: list[str] = []

    def refresh(self) -> None:
        project = self.project_service.current_project
        if project is None:
            self.context_label.setText("Load a project first.")
            self._tab_keys = []
            self.comparison_tabs.clear()
            self._populate_table(self.catalog_table, pd.DataFrame())
            self._populate_table(self.selected_table, pd.DataFrame())
            self._populate_table(self.blacklist_table, pd.DataFrame())
            return
        cached = self.project_service.run_state.results.candidate_gene_table
        if cached is None or cached.empty:
            self.context_label.setText("No cached candidate ranking is loaded yet. Click 'Load Candidate Ranking' for the current project/comparison.")
            self._refresh_tabs()
            self._populate_table(self.catalog_table, pd.DataFrame())
            self._populate_table(self.selected_table, pd.DataFrame())
            self._populate_table(self.blacklist_table, pd.DataFrame())
            return

        self._refresh_tabs()
        comparison_id = self._current_comparison_id()
        context = self.project_service.comparison_context(comparison_id)
        self.context_label.setText(
            f"{context['display_name']} [{context['comparison_id']}] | {context['group1']} vs {context['group2']} | "
            "Default selection = top numeric rank genes for this comparison. Blacklist genes are excluded from downstream heatmap/cards/sashimi."
        )

        catalog = self.project_service.candidate_selection_catalog(comparison_id, allow_generate=False)
        catalog = self._apply_catalog_search(catalog)
        self._current_catalog = catalog.copy()
        self._populate_table(self.catalog_table, catalog)

        selected, _meta = self.project_service.candidate_selection_frame(
            comparison_id,
            source="selection",
            top_n=int(self.top_n.value()),
            allow_generate=False,
        )
        selected_columns = [
            column
            for column in [
                "rank",
                "gene_symbol",
                "gene_id",
                "candidate_tier",
                "evidence_class",
                "direction_class",
                "DEG_significant",
                "rMATS_significant",
                "DEXSeq_significant",
                "DTU_significant",
            ]
            if column in selected.columns
        ]
        self._current_selected = selected.copy()
        self._populate_table(self.selected_table, selected[selected_columns].copy() if selected_columns else pd.DataFrame())

        blacklist = self.project_service.candidate_blacklist_genes_for_display(comparison_id)
        blacklist_frame = pd.DataFrame({"gene": blacklist})
        self._populate_table(self.blacklist_table, blacklist_frame)

        if catalog.empty and selected.empty:
            self.selection_status.setText(
                "No cached candidate ranking is loaded for this comparison yet. Click 'Load Candidate Ranking' first. "
                "Selection state and blacklist are lightweight metadata and are not recomputed at project load."
            )
            return

        self.selection_status.setText(
            f"Selected genes: {len(self.project_service.candidate_selection_genes_for_display(comparison_id))} | "
            f"Blacklist genes: {len(blacklist)} | "
            "Downstream evidence heatmap / candidate cards use this selection and exclude blacklist genes."
        )

    def _load_candidate_ranking(self) -> None:
        comparison_id = self._current_comparison_id()
        if not comparison_id:
            return
        frame = self.project_service.preview_candidate_gene_screening(allow_generate=True)
        if frame is not None and not frame.empty:
            existing = self.project_service.candidate_selection_genes_for_display(comparison_id)
            if not existing:
                self.project_service.reset_candidate_selection(comparison_id, top_n=int(self.top_n.value()))
        if frame is None or frame.empty:
            if self.project_service.run_state.results.tx_splicing_gene_table.empty:
                self.selection_status.setText(
                    "Missing input: transcript + splicing integration not run for this project. Please run 5.2.1 first."
                )
            else:
                self.selection_status.setText(
                    "Candidate ranking is still unavailable for this project/comparison. Check cached inputs and try again."
                )
        self.refresh()

    def _refresh_tabs(self) -> None:
        comparisons = self.project_service.selected_comparisons_for_display()
        current_key = self._current_comparison_id()
        self._tab_keys = [item.comparison_id for item in comparisons]
        self.comparison_tabs.blockSignals(True)
        self.comparison_tabs.clear()
        if not comparisons:
            self.comparison_tabs.addTab(QWidget(), "No selected comparisons")
            self._tab_keys = []
        else:
            current_index = 0
            for index, comparison in enumerate(comparisons):
                self.comparison_tabs.addTab(QWidget(), comparison.display_resolved_name)
                if current_key and comparison.comparison_id == current_key:
                    current_index = index
            self.comparison_tabs.setCurrentIndex(current_index)
        self.comparison_tabs.blockSignals(False)

    def _current_comparison_id(self) -> str | None:
        index = self.comparison_tabs.currentIndex()
        if index < 0 or index >= len(self._tab_keys):
            return None
        return self._tab_keys[index]

    def _apply_catalog_search(self, frame: pd.DataFrame) -> pd.DataFrame:
        query = self.search.text().strip().lower()
        if frame.empty or not query:
            return frame
        mask = pd.Series(False, index=frame.index)
        for column in ("gene_symbol", "gene_id", "candidate_tier", "evidence_class", "direction_class"):
            if column in frame.columns:
                mask |= frame[column].fillna("").astype(str).str.lower().str.contains(query, na=False)
        return frame.loc[mask].copy()

    def _selected_catalog_genes(self) -> list[str]:
        return self._selected_genes_from_table(self.catalog_table)

    def _selected_genes_from_selection_table(self) -> list[str]:
        return self._selected_genes_from_table(self.selected_table)

    def _selected_blacklist_genes(self) -> list[str]:
        return self._selected_genes_from_table(self.blacklist_table)

    def _reset_top_n(self) -> None:
        comparison_id = self._current_comparison_id()
        if not comparison_id:
            return
        self.project_service.reset_candidate_selection(comparison_id, top_n=int(self.top_n.value()))
        self.refresh()

    def _add_manual_gene(self) -> None:
        comparison_id = self._current_comparison_id()
        value = self.gene_input.text().strip()
        if not comparison_id or not value:
            return
        self.project_service.add_candidate_selection_gene(comparison_id, value)
        self.gene_input.clear()
        self.refresh()

    def _add_selected_catalog_gene(self) -> None:
        comparison_id = self._current_comparison_id()
        genes = self._selected_catalog_genes()
        if not comparison_id or not genes:
            return
        for gene in genes:
            self.project_service.add_candidate_selection_gene(comparison_id, gene)
        self.refresh()

    def _remove_selected_genes(self) -> None:
        comparison_id = self._current_comparison_id()
        genes = self._selected_genes_from_selection_table()
        if not comparison_id or not genes:
            return
        for gene in genes:
            self.project_service.remove_candidate_selection_gene(comparison_id, gene)
        self.refresh()

    def _add_selected_catalog_genes(self) -> None:
        self._add_selected_catalog_gene()

    def _add_manual_blacklist_gene(self) -> None:
        comparison_id = self._current_comparison_id()
        value = self.blacklist_input.text().strip()
        if not comparison_id or not value:
            return
        self.project_service.add_candidate_blacklist_gene(comparison_id, value)
        self.blacklist_input.clear()
        self.refresh()

    def _blacklist_selected_catalog_gene(self) -> None:
        comparison_id = self._current_comparison_id()
        genes = self._selected_catalog_genes()
        if not comparison_id or not genes:
            return
        for gene in genes:
            self.project_service.add_candidate_blacklist_gene(comparison_id, gene)
        self.refresh()

    def _blacklist_selected_catalog_genes(self) -> None:
        self._blacklist_selected_catalog_gene()

    def _remove_selected_blacklist_genes(self) -> None:
        comparison_id = self._current_comparison_id()
        genes = self._selected_blacklist_genes()
        if not comparison_id or not genes:
            return
        for gene in genes:
            self.project_service.remove_candidate_blacklist_gene(comparison_id, gene)
        self.refresh()

    @staticmethod
    def _selected_genes_from_table(table: QTableWidget) -> list[str]:
        rows = sorted({index.row() for index in table.selectionModel().selectedRows()})
        if not rows:
            row = table.currentRow()
            if row >= 0:
                rows = [row]
        if not rows:
            return []
        headers = [
            str(table.horizontalHeaderItem(column).text()) if table.horizontalHeaderItem(column) is not None else ""
            for column in range(table.columnCount())
        ]
        gene_symbol_col = headers.index("gene_symbol") if "gene_symbol" in headers else -1
        gene_id_col = headers.index("gene_id") if "gene_id" in headers else -1
        generic_gene_col = headers.index("gene") if "gene" in headers else -1
        values: list[str] = []
        seen: set[str] = set()
        for row in rows:
            candidates = []
            for column_index in (gene_symbol_col, gene_id_col, generic_gene_col):
                if column_index >= 0:
                    item = table.item(row, column_index)
                    if item is not None:
                        candidates.append(item.text().strip())
            gene = next((value for value in candidates if value), "")
            if gene and gene not in seen:
                seen.add(gene)
                values.append(gene)
        return values

    def _populate_table(self, table: QTableWidget, frame: pd.DataFrame) -> None:
        table.setSortingEnabled(False)
        if frame is None or frame.empty:
            table.setRowCount(0)
            table.setColumnCount(0)
            table.setSortingEnabled(True)
            return
        working = frame.reset_index(drop=True).copy()
        table.setColumnCount(len(working.columns))
        table.setHorizontalHeaderLabels([str(column) for column in working.columns])
        table.setRowCount(len(working))
        for row_idx, (_, row) in enumerate(working.iterrows()):
            for col_idx, column in enumerate(working.columns):
                value = row[column]
                item = QTableWidgetItem()
                if pd.isna(value):
                    item.setText("")
                elif self._is_numeric_column(working, column):
                    numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
                    if pd.notna(numeric):
                        item.setData(Qt.ItemDataRole.EditRole, float(numeric) if not float(numeric).is_integer() else int(numeric))
                    else:
                        item.setText(str(value))
                elif self._is_boolean_column(working, column):
                    item.setData(Qt.ItemDataRole.EditRole, bool(value))
                else:
                    item.setText(str(value))
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                table.setItem(row_idx, col_idx, item)
        table.resizeColumnsToContents()
        table.setSortingEnabled(True)

    @staticmethod
    def _is_numeric_column(frame: pd.DataFrame, column: str) -> bool:
        if column not in frame.columns:
            return False
        return is_numeric_dtype(frame[column]) and not is_bool_dtype(frame[column])

    @staticmethod
    def _is_boolean_column(frame: pd.DataFrame, column: str) -> bool:
        if column not in frame.columns:
            return False
        return is_bool_dtype(frame[column])
