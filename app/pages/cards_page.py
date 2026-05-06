from __future__ import annotations

from pathlib import Path

import pandas as pd
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QComboBox,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMessageBox,
    QPushButton,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from src.services.project_service import ProjectService


class CardsPage(QWidget):
    MAX_COMPARISONS = 8

    def __init__(self, project_service: ProjectService) -> None:
        super().__init__()
        self.project_service = project_service
        self.label = QLabel("Card workbench")
        self.explanation = QLabel(
            "Cards summarize one gene across the selected comparisons. "
            "Top = DEG effect, middle = representative splicing effect, bottom = expression support. "
            "Choose one or more genes and up to 8 comparisons, then export one card per gene. "
            "Cards use shortlist-derived genes for downstream validation and respect the Blacklist page."
        )
        self.explanation.setWordWrap(True)
        self.gene_source = QComboBox()
        self.gene_source.addItem("Shortlist genes", "shortlist")
        self.gene_source.addItem("All detected genes", "all")
        self.gene_source.currentIndexChanged.connect(self._refresh_gene_list)
        self.visual_group = QComboBox()
        self.visual_group.currentIndexChanged.connect(self._populate_comparison_list)
        self.visual_group.currentIndexChanged.connect(self._refresh_gene_list)
        self.gene_search = QLineEdit()
        self.gene_search.setPlaceholderText("Search genes...")
        self.gene_search.textChanged.connect(self._refresh_gene_list)
        self.custom_gene_input = QLineEdit()
        self.custom_gene_input.setPlaceholderText("Add custom gene symbol or gene ID")
        self.add_custom_gene_button = QPushButton("Add Gene")
        self.add_custom_gene_button.clicked.connect(self._add_custom_gene)
        self.generate_button = QPushButton("Export Selected Cards")
        self.generate_button.clicked.connect(self._export_selected_cards)
        self.export_png_button = QPushButton("Export Current PNG")
        self.export_png_button.clicked.connect(lambda: self._export_current_card("png"))
        self.export_pdf_button = QPushButton("Export Current PDF")
        self.export_pdf_button.clicked.connect(lambda: self._export_current_card("pdf"))

        self.gene_list = QListWidget()
        self.gene_list.setSelectionMode(QListWidget.SelectionMode.ExtendedSelection)
        self.gene_list.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.gene_list.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.gene_list.currentItemChanged.connect(self._render_preview)
        self.gene_list.itemSelectionChanged.connect(self._render_preview)

        self.comparison_list = QListWidget()
        self.comparison_list.setSelectionMode(QListWidget.SelectionMode.ExtendedSelection)
        self.comparison_list.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.comparison_list.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.comparison_list.itemSelectionChanged.connect(self._enforce_comparison_limit)
        self.comparison_list.itemSelectionChanged.connect(self._render_preview)

        self.figure = Figure(figsize=(10, 8), tight_layout=True)
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        self.details = QTextEdit()
        self.details.setReadOnly(True)

        left = QVBoxLayout()
        left.addWidget(QLabel("Genes"))
        left.addWidget(QLabel("Visual group"))
        left.addWidget(self.visual_group)
        left.addWidget(self.gene_source)
        left.addWidget(self.gene_search)
        custom_row = QHBoxLayout()
        custom_row.addWidget(self.custom_gene_input)
        custom_row.addWidget(self.add_custom_gene_button)
        left.addLayout(custom_row)
        left.addWidget(self.gene_list, 3)
        left.addWidget(QLabel("Comparisons (max 8)"))
        left.addWidget(self.comparison_list, 2)

        right = QVBoxLayout()
        button_row = QHBoxLayout()
        button_row.addWidget(self.generate_button)
        button_row.addWidget(self.export_png_button)
        button_row.addWidget(self.export_pdf_button)
        right.addLayout(button_row)
        right.addWidget(self.toolbar)
        right.addWidget(self.canvas, 3)
        right.addWidget(self.details, 1)

        top = QHBoxLayout()
        top.addLayout(left, 1)
        top.addLayout(right, 3)

        layout = QVBoxLayout(self)
        layout.addWidget(self.label)
        layout.addWidget(self.explanation)
        layout.addLayout(top)

        self._tx_frame = pd.DataFrame()
        self._expr_frame = pd.DataFrame()
        self._gene_catalog = pd.DataFrame()
        self._shortlist_gene_catalog = pd.DataFrame()
        self._custom_gene_catalog = pd.DataFrame(columns=["gene_key", "gene_label"])
        self._selected_gene_key: str | None = None
        self._comparison_lookup: dict[str, str] = {}

    def refresh(self) -> None:
        project = self.project_service.current_project
        if project is None:
            self._tx_frame = pd.DataFrame()
            self._expr_frame = pd.DataFrame()
            self._gene_catalog = pd.DataFrame()
            self._shortlist_gene_catalog = pd.DataFrame()
            self.gene_list.clear()
            self.comparison_list.clear()
            self.details.setPlainText("")
            self._draw_placeholder("Load a project first.")
            return

        run_state = self.project_service.run_state
        if run_state.status != "completed":
            self._tx_frame = pd.DataFrame()
            self._expr_frame = pd.DataFrame()
            self._gene_catalog = pd.DataFrame()
            self.gene_list.clear()
            self.comparison_list.clear()
            self.details.setPlainText("Run analysis first to populate gene-level DEG/splicing integration.")
            self._draw_placeholder("Run analysis first.")
            return

        self._tx_frame = self.project_service.apply_gene_blacklist(run_state.results.tx_splicing_gene_table.copy())
        self._expr_frame = self.project_service.apply_gene_blacklist(run_state.results.cards_expression_support.copy())
        self._comparison_lookup = {
            comparison.comparison_id: comparison.resolved_name
            for comparison in self.project_service.selected_comparisons_for_display()
        }
        self._populate_visual_groups()
        self._gene_catalog = self._build_gene_catalog(self._tx_frame)
        self._shortlist_gene_catalog = self._build_shortlist_gene_catalog(
            self.project_service.apply_gene_blacklist(run_state.results.cards_shortlist.copy())
        )
        self._populate_comparison_list()
        self._refresh_gene_list()
        self._render_preview()

    def _build_gene_catalog(self, frame: pd.DataFrame) -> pd.DataFrame:
        if frame.empty:
            return pd.DataFrame(columns=["gene_key", "gene_label"])
        working = frame.copy()
        working["gene_key"] = working["gene_id"].astype(str)
        working["gene_label"] = working["geneSymbol"].fillna(working["gene_id"]).astype(str)
        catalog = (
            working[["gene_key", "gene_label"]]
            .drop_duplicates()
            .sort_values(["gene_label", "gene_key"], ascending=[True, True])
            .reset_index(drop=True)
        )
        return catalog

    def _build_shortlist_gene_catalog(self, frame: pd.DataFrame) -> pd.DataFrame:
        if frame.empty:
            return pd.DataFrame(columns=["gene_key", "gene_label"])
        working = frame.copy()
        working["gene_key"] = working["GeneID"].fillna(working.get("gene_id")).astype(str)
        working["gene_label"] = working["geneSymbol"].fillna(working["gene_key"]).astype(str)
        if "pair_name" in working.columns:
            working["gene_label"] = working["gene_label"] + " | " + working["pair_name"].fillna("").astype(str)
        catalog = (
            working[["gene_key", "gene_label"]]
            .drop_duplicates()
            .sort_values(["gene_label", "gene_key"], ascending=[True, True])
            .reset_index(drop=True)
        )
        return catalog

    def _populate_comparison_list(self) -> None:
        selected = {item.data(Qt.ItemDataRole.UserRole) for item in self.comparison_list.selectedItems()}
        values: list[str] = []
        allowed_names = self._selected_visual_group_comparisons()
        if allowed_names:
            values = allowed_names
        elif not self._tx_frame.empty and "comparison_name" in self._tx_frame.columns:
            values = self._tx_frame["comparison_name"].dropna().astype(str).drop_duplicates().tolist()
        self.comparison_list.blockSignals(True)
        self.comparison_list.clear()
        for value in values:
            item = QListWidgetItem(value)
            item.setData(Qt.ItemDataRole.UserRole, value)
            item.setSelected(value in selected)
            self.comparison_list.addItem(item)
        if not selected and self.comparison_list.count() > 0:
            for index in range(self.comparison_list.count()):
                self.comparison_list.item(index).setSelected(True)
        self.comparison_list.blockSignals(False)

    def _populate_visual_groups(self) -> None:
        current = self.visual_group.currentData()
        groups = self.project_service.active_visualization_groups_for_display()
        self.visual_group.blockSignals(True)
        self.visual_group.clear()
        self.visual_group.addItem("All selected comparisons", None)
        for group in groups:
            self.visual_group.addItem(group.resolved_name, group.group_id)
        if current is not None:
            index = self.visual_group.findData(current)
            if index >= 0:
                self.visual_group.setCurrentIndex(index)
        self.visual_group.blockSignals(False)

    def _selected_visual_group_comparisons(self) -> list[str]:
        group_id = self.visual_group.currentData()
        if not group_id:
            return []
        group = next(
            (item for item in self.project_service.active_visualization_groups_for_display() if item.group_id == group_id),
            None,
        )
        if group is None:
            return []
        names: list[str] = []
        for comparison_id in group.comparison_ids:
            label = self._comparison_lookup.get(comparison_id)
            if label:
                names.append(label)
        return names

    def _refresh_gene_list(self) -> None:
        query = self.gene_search.text().strip().lower()
        selected = {item.data(Qt.ItemDataRole.UserRole) for item in self.gene_list.selectedItems()}
        current = self.gene_list.currentItem().data(Qt.ItemDataRole.UserRole) if self.gene_list.currentItem() else None
        source = self.gene_source.currentData()
        catalog = self._shortlist_gene_catalog if source == "shortlist" and not self._shortlist_gene_catalog.empty else self._gene_catalog
        if not self._custom_gene_catalog.empty:
            catalog = pd.concat([catalog, self._custom_gene_catalog], axis=0, ignore_index=True).drop_duplicates("gene_key", keep="first")
        if query and not catalog.empty:
            mask = catalog["gene_label"].astype(str).str.lower().str.contains(query, na=False) | catalog["gene_key"].astype(str).str.lower().str.contains(query, na=False)
            catalog = catalog.loc[mask].copy()
        self.gene_list.blockSignals(True)
        self.gene_list.clear()
        for _, row in catalog.iterrows():
            label = f"{row['gene_label']} ({row['gene_key']})"
            item = QListWidgetItem(label)
            item.setData(Qt.ItemDataRole.UserRole, row["gene_key"])
            item.setSelected(row["gene_key"] in selected)
            self.gene_list.addItem(item)
        if self.gene_list.count() > 0:
            target_key = current or self._selected_gene_key
            target_row = 0
            if target_key is not None:
                for index in range(self.gene_list.count()):
                    if self.gene_list.item(index).data(Qt.ItemDataRole.UserRole) == target_key:
                        target_row = index
                        break
            self.gene_list.setCurrentRow(target_row)
        self.gene_list.blockSignals(False)

    def _selected_comparisons(self) -> list[str]:
        return [str(item.data(Qt.ItemDataRole.UserRole)) for item in self.comparison_list.selectedItems()]

    def _selected_gene_keys(self) -> list[str]:
        return [str(item.data(Qt.ItemDataRole.UserRole)) for item in self.gene_list.selectedItems()]

    def _active_gene_key(self) -> str | None:
        current = self.gene_list.currentItem()
        if current is not None:
            return str(current.data(Qt.ItemDataRole.UserRole))
        selected = self._selected_gene_keys()
        return selected[0] if selected else None

    def _enforce_comparison_limit(self) -> None:
        selected = self.comparison_list.selectedItems()
        if len(selected) <= self.MAX_COMPARISONS:
            return
        for item in selected[self.MAX_COMPARISONS :]:
            item.setSelected(False)
        QMessageBox.information(self, "Comparison limit", f"Select up to {self.MAX_COMPARISONS} comparisons for card preview/export.")

    def _render_preview(self) -> None:
        if self._tx_frame.empty:
            self._draw_placeholder("No card inputs available yet.")
            return
        gene_key = self._active_gene_key()
        selected_comparisons = self._selected_comparisons()
        if gene_key is None:
            self._draw_placeholder("Select at least one gene.")
            return
        if not selected_comparisons:
            self._draw_placeholder("Select at least one comparison.")
            return
        self._selected_gene_key = gene_key
        gene_frame = self._tx_frame.loc[
            (self._tx_frame["gene_id"].astype(str) == gene_key)
            & (self._tx_frame["comparison_name"].astype(str).isin(selected_comparisons))
        ].copy()
        expr = self._expr_frame.loc[self._expr_frame["gene_id"].astype(str) == gene_key].copy() if not self._expr_frame.empty and "gene_id" in self._expr_frame.columns else pd.DataFrame()
        self._draw_gene_card(gene_frame, expr, gene_key, selected_comparisons)

    def _draw_placeholder(self, text: str) -> None:
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        axis.axis("off")
        axis.text(0.5, 0.5, text, ha="center", va="center", wrap=True)
        self.details.setPlainText(text)
        self.canvas.draw_idle()

    def _add_custom_gene(self) -> None:
        value = self.custom_gene_input.text().strip()
        if not value:
            return
        gene_key = value
        gene_label = value
        if not self._gene_catalog.empty:
            mask = (
                self._gene_catalog["gene_key"].astype(str).str.lower() == value.lower()
            ) | (
                self._gene_catalog["gene_label"].astype(str).str.lower().str.split(" \\| ").str[0] == value.lower()
            )
            matches = self._gene_catalog.loc[mask].copy()
            if not matches.empty:
                gene_key = str(matches.iloc[0]["gene_key"])
                gene_label = str(matches.iloc[0]["gene_label"])
        new_row = pd.DataFrame([{"gene_key": gene_key, "gene_label": gene_label}])
        self._custom_gene_catalog = pd.concat([self._custom_gene_catalog, new_row], axis=0, ignore_index=True).drop_duplicates("gene_key", keep="first")
        self.custom_gene_input.clear()
        self._refresh_gene_list()
        for index in range(self.gene_list.count()):
            item = self.gene_list.item(index)
            if item.data(Qt.ItemDataRole.UserRole) == gene_key:
                item.setSelected(True)
                self.gene_list.setCurrentRow(index)
                break

    def _draw_gene_card(self, frame: pd.DataFrame, expr: pd.DataFrame, gene_key: str, selected_comparisons: list[str]) -> None:
        self.figure.clear()
        top = self.figure.add_subplot(311)
        middle = self.figure.add_subplot(312)
        bottom = self.figure.add_subplot(313)

        gene_label = frame["geneSymbol"].dropna().astype(str).iloc[0] if not frame.empty and "geneSymbol" in frame.columns else gene_key
        self.label.setText(f"Card workbench | {gene_label} | {len(selected_comparisons)} comparison(s)")

        if frame.empty:
            top.axis("off")
            middle.axis("off")
            bottom.axis("off")
            top.text(0.5, 0.5, "No DEG/splicing rows for this gene and comparison selection.", ha="center", va="center")
            self.details.setPlainText(f"No rows found for {gene_label} ({gene_key}) in the selected comparisons.")
            self.canvas.draw_idle()
            return

        frame = frame.copy()
        frame["log2FC"] = pd.to_numeric(frame["log2FC"], errors="coerce")
        frame["representative_dPSI"] = pd.to_numeric(frame["representative_dPSI"], errors="coerce")
        labels = frame["comparison_name"].astype(str).tolist()

        top.axhline(0, color="#777777", linewidth=1)
        top.bar(labels, frame["log2FC"].fillna(0).tolist(), color="#00798C")
        top.set_ylabel("log2FC")
        top.set_title("DEG effect by comparison")
        top.tick_params(axis="x", rotation=20)

        middle.axhline(0, color="#777777", linewidth=1)
        middle.bar(labels, frame["representative_dPSI"].fillna(0).tolist(), color="#D1495B")
        middle.set_ylabel("dPSI")
        middle.set_title("Representative splicing effect by comparison")
        middle.tick_params(axis="x", rotation=20)

        if expr.empty:
            bottom.axis("off")
            bottom.text(0.5, 0.5, "No expression support rows found for this gene.", ha="center", va="center")
        else:
            groups = expr["group"].astype(str).tolist()
            values = pd.to_numeric(expr["expr"], errors="coerce").fillna(0).tolist()
            bottom.bar(groups, values, color="#6D597A")
            bottom.set_ylabel("Mean expression")
            bottom.set_title("Expression support")
            bottom.tick_params(axis="x", rotation=20)

        detail_lines = [f"Gene: {gene_label}", f"Gene ID: {gene_key}", ""]
        for _, row in frame.iterrows():
            detail_lines.extend(
                [
                    f"Comparison: {row.get('comparison_name', '')}",
                    f"  log2FC: {row.get('log2FC', '')}",
                    f"  DEG padj: {row.get('deg_padj', '')}",
                    f"  dPSI: {row.get('representative_dPSI', '')}",
                    f"  Event FDR: {row.get('representative_event_FDR', '')}",
                    f"  Combined class: {row.get('combined_class', '')}",
                    f"  Agreement: {row.get('agreement_class', '')}",
                    "",
                ]
            )
        self.details.setPlainText("\n".join(detail_lines).rstrip())
        self.canvas.draw_idle()

    def _export_current_card(self, fmt: str) -> None:
        gene_key = self._selected_gene_key or "card"
        default_name = f"card_{gene_key}.{fmt}"
        output_path, _ = QFileDialog.getSaveFileName(
            self,
            f"Export Current Card as {fmt.upper()}",
            str(Path.cwd() / default_name),
            f"{fmt.upper()} Files (*.{fmt})",
        )
        if not output_path:
            return
        self.figure.savefig(output_path, dpi=300, bbox_inches="tight")

    def _export_selected_cards(self) -> None:
        project = self.project_service.current_project
        if project is None or project.output_root is None:
            QMessageBox.information(self, "No output directory", "Configure the project output directory first.")
            return
        selected_genes = self._selected_gene_keys()
        selected_comparisons = self._selected_comparisons()
        if not selected_genes:
            QMessageBox.information(self, "No genes selected", "Select at least one gene to export cards.")
            return
        if not selected_comparisons:
            QMessageBox.information(self, "No comparisons selected", "Select at least one comparison.")
            return

        output_dir = project.output_root / "05_cards" / "manual_cards"
        output_dir.mkdir(parents=True, exist_ok=True)
        generated = []
        current_gene = self._selected_gene_key
        for gene_key in selected_genes:
            gene_frame = self._tx_frame.loc[
                (self._tx_frame["gene_id"].astype(str) == gene_key)
                & (self._tx_frame["comparison_name"].astype(str).isin(selected_comparisons))
            ].copy()
            expr = self._expr_frame.loc[self._expr_frame["gene_id"].astype(str) == gene_key].copy() if not self._expr_frame.empty and "gene_id" in self._expr_frame.columns else pd.DataFrame()
            self._draw_gene_card(gene_frame, expr, gene_key, selected_comparisons)
            out = output_dir / f"card_{gene_key}.png"
            self.figure.savefig(out, dpi=300, bbox_inches="tight")
            generated.append(str(out))
        if current_gene is not None:
            self._selected_gene_key = current_gene
            self._render_preview()
        self.details.setPlainText("Exported cards:\n" + "\n".join(generated))
