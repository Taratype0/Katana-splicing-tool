from __future__ import annotations

import json
import os
import shutil
from pathlib import Path

from PySide6.QtCore import QObject, QThread, Signal, Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QFileDialog,
    QDialog,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from app.pages.module_details import module_detail
from app.widgets.artifact_image_preview import ArtifactImagePreview
from katana_cosplice.config import load_config
from katana_cosplice.discovery import discover_config
from katana_cosplice.pipeline import CoSplicePipeline
from katana_cosplice.stage_registry import stage_spec, stage_specs
from src.services.project_service import ProjectService


class CoSpliceWorker(QObject):
    progress = Signal(str)
    finished = Signal(str)
    failed = Signal(str)

    def __init__(self, config_path: Path, stage: str) -> None:
        super().__init__()
        self.config_path = config_path
        self.stage = stage

    def run(self) -> None:
        try:
            pipeline = CoSplicePipeline(load_config(self.config_path), progress=self.progress.emit)
            pipeline.run_stage(self.stage)
        except Exception as exc:
            self.failed.emit(str(exc))
            return
        self.finished.emit(self.stage)


class CoSplicePage(QWidget):
    """Standalone CoSplice runner embedded in Katana without implicit execution."""

    def __init__(self, project_service: ProjectService) -> None:
        super().__init__()
        self.project_service = project_service
        self._thread: QThread | None = None
        self._worker: CoSpliceWorker | None = None
        title = QLabel("9. CoSplice Network")
        title.setStyleSheet("font-size: 20px; font-weight: 750;")
        subtitle = QLabel(
            "Discover replicate-supported co-splicing programs from sample-level rMATS PSI. "
            "Each stage runs independently; opening this page reads cache/state only."
        )
        subtitle.setWordWrap(True)

        self.config_path = QLineEdit()
        self.config_path.setPlaceholderText("Path to cosplice_config.yaml")
        choose = QPushButton("Choose config")
        choose.clicked.connect(self._choose_config)
        initialize = QPushButton("Initialize from current project")
        initialize.clicked.connect(self._initialize_config)
        detail = QPushButton("Detail")
        detail.clicked.connect(self._show_detail)
        open_output = QPushButton("Open output folder")
        open_output.clicked.connect(self._open_output)

        config_row = QHBoxLayout()
        config_row.addWidget(QLabel("Config"))
        config_row.addWidget(self.config_path, 1)
        config_row.addWidget(choose)
        config_row.addWidget(initialize)
        config_row.addWidget(detail)
        config_row.addWidget(open_output)

        self.matrix_source = QComboBox()
        self.matrix_source.addItem("rMATS event matrix", "rmats_events")
        self.matrix_source.addItem("BAM competing-junction matrix", "bam_junctions")
        save_matrix_source = QPushButton("Save matrix source")
        save_matrix_source.clicked.connect(self._save_matrix_source)
        self.matrix_source_status = QLabel(
            "Matrix source is stored in the config. Changing it invalidates Stage 2 and downstream caches."
        )
        self.matrix_source_status.setWordWrap(True)
        matrix_row = QHBoxLayout()
        matrix_row.addWidget(QLabel("Matrix source"))
        matrix_row.addWidget(self.matrix_source)
        matrix_row.addWidget(save_matrix_source)
        matrix_row.addWidget(self.matrix_source_status, 1)

        run_grid = QGridLayout()
        self.run_buttons: list[QPushButton] = []
        for index, spec in enumerate(stage_specs(visible_only=True)):
            button = QPushButton(spec.title)
            button.setToolTip(spec.description)
            button.setMinimumHeight(34)
            button.clicked.connect(lambda _checked=False, value=spec.key: self._run_stage(value))
            run_grid.addWidget(button, index // 4, index % 4)
            self.run_buttons.append(button)

        self.status_table = QTableWidget(0, 5)
        self.status_table.setHorizontalHeaderLabels(["Stage", "Status", "Updated", "Message", "Outputs"])
        self.status_table.setSortingEnabled(True)
        self.status_table.horizontalHeader().setStretchLastSection(True)

        self.plot_table = QTableWidget(0, 2)
        self.plot_table.setHorizontalHeaderLabels(["Plot", "Path"])
        self.plot_table.setSortingEnabled(True)
        self.plot_table.currentCellChanged.connect(self._preview_plot)
        self.preview_scroll = ArtifactImagePreview("No cached CoSplice plot result.")
        self.preview_scroll.setMinimumHeight(420)
        export_plot = QPushButton("Export selected plot")
        export_plot.clicked.connect(self._export_plot)
        plot_tab = QWidget()
        plot_layout = QVBoxLayout(plot_tab)
        plot_layout.addWidget(export_plot, 0, Qt.AlignmentFlag.AlignLeft)
        plot_layout.addWidget(self.plot_table)
        plot_layout.addWidget(self.preview_scroll, 1)

        self.module_table = QTableWidget(0, 0)
        self.module_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.module_table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.module_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.module_table.setSortingEnabled(True)
        publish_modules = QPushButton("Publish selected modules to downstream")
        publish_modules.clicked.connect(self._publish_selected_modules)
        self.module_handoff_status = QLabel(
            "Select cached modules and publish them as a project artifact. This does not rerun CoSplice."
        )
        self.module_handoff_status.setWordWrap(True)
        module_tab = QWidget()
        module_layout = QVBoxLayout(module_tab)
        module_layout.addWidget(publish_modules, 0, Qt.AlignmentFlag.AlignLeft)
        module_layout.addWidget(self.module_handoff_status)
        module_layout.addWidget(self.module_table, 1)

        self.log = QTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumHeight(150)
        self.log.setPlaceholderText("Stage log appears here. Cache refresh does not run analysis.")

        tabs = QTabWidget()
        status_tab = QWidget()
        status_layout = QVBoxLayout(status_tab)
        status_layout.addWidget(self.status_table)
        tabs.addTab(status_tab, "Stage status")
        tabs.addTab(plot_tab, "Plots / Preview")
        tabs.addTab(module_tab, "Module handoff")

        layout = QVBoxLayout(self)
        layout.addWidget(title)
        layout.addWidget(subtitle)
        layout.addLayout(config_row)
        layout.addLayout(matrix_row)
        layout.addLayout(run_grid)
        layout.addWidget(tabs, 1)
        layout.addWidget(self.log)

    def _resolved_config_path(self) -> Path | None:
        raw = self.config_path.text().strip()
        return Path(raw) if raw else None

    def _choose_config(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Choose CoSplice config", "", "YAML (*.yaml *.yml)")
        if path:
            self.config_path.setText(path)
            self.refresh()

    def _initialize_config(self) -> None:
        project = self.project_service.current_project
        if project is None:
            QMessageBox.information(self, "No project", "Open a Katana project first.")
            return
        try:
            output_root = self.project_service.require_output_root()
        except RuntimeError as exc:
            QMessageBox.information(self, "Output directory required", str(exc))
            return
        output_dir = output_root / "09_cosplice_network"
        try:
            config = discover_config(project.project_root, output_dir=output_dir)
            config.matrix_input.input_mode = str(self.matrix_source.currentData())
            self._validate_matrix_source(config)
            path = config.write(output_dir / "cosplice_config.yaml")
        except Exception as exc:
            QMessageBox.critical(self, "Config discovery failed", str(exc))
            return
        self.config_path.setText(str(path))
        self.log.setPlainText(f"Config initialized: {path}\nNo analysis was run.")
        self.refresh()

    @staticmethod
    def _validate_matrix_source(config) -> None:
        if config.matrix_input.input_mode != "bam_junctions":
            return
        missing = []
        if config.bam_root is None or not config.bam_root.exists():
            missing.append("BAM root")
        if config.gtf is None or not config.gtf.exists():
            missing.append("GTF")
        if missing:
            raise RuntimeError(
                "BAM competing-junction matrix requires: "
                + ", ".join(missing)
                + ". Keep rMATS event matrix selected or configure these inputs first."
            )

    def _save_matrix_source(self) -> None:
        path = self._resolved_config_path()
        if path is None or not path.exists():
            QMessageBox.information(self, "Missing config", "Choose or initialize cosplice_config.yaml first.")
            return
        try:
            config = load_config(path)
            config.matrix_input.input_mode = str(self.matrix_source.currentData())
            self._validate_matrix_source(config)
            config.validate()
            config.write(path)
        except Exception as exc:
            QMessageBox.critical(self, "Matrix source not saved", str(exc))
            return
        self.log.setPlainText(
            f"Saved matrix source: {config.matrix_input.input_mode}\n"
            "No analysis was run. Run Stage 2 explicitly to rebuild the matrix."
        )
        self.refresh()

    def _run_stage(self, stage: str) -> None:
        if self._thread is not None:
            QMessageBox.information(self, "CoSplice running", "Wait for the current CoSplice stage to finish.")
            return
        path = self._resolved_config_path()
        if path is None or not path.exists():
            QMessageBox.information(self, "Missing config", "Choose or initialize cosplice_config.yaml first.")
            return
        for button in self.run_buttons:
            button.setEnabled(False)
        spec = stage_spec(stage)
        self.log.setPlainText(f"Running only stage: {spec.title} [{stage}]\nNo upstream or downstream stage will be started automatically.")
        self._thread = QThread(self)
        self._worker = CoSpliceWorker(path, stage)
        self._worker.moveToThread(self._thread)
        self._thread.started.connect(self._worker.run)
        self._worker.progress.connect(self._append_log)
        self._worker.finished.connect(self._on_finished)
        self._worker.failed.connect(self._on_failed)
        self._worker.finished.connect(self._thread.quit)
        self._worker.failed.connect(self._thread.quit)
        self._thread.finished.connect(self._cleanup_worker)
        self._thread.start()

    def _append_log(self, message: str) -> None:
        self.log.append(message)

    def _on_finished(self, stage: str) -> None:
        self.log.append(f"Finished: {stage}")
        self.refresh()

    def _on_failed(self, message: str) -> None:
        self.log.append(f"Failed: {message}")
        QMessageBox.critical(self, "CoSplice stage failed", message)
        self.refresh()

    def _cleanup_worker(self) -> None:
        for button in self.run_buttons:
            button.setEnabled(True)
        if self._worker is not None:
            self._worker.deleteLater()
        if self._thread is not None:
            self._thread.deleteLater()
        self._worker = None
        self._thread = None

    def refresh(self) -> None:
        path = self._resolved_config_path()
        if path is None and self.project_service.current_project is not None:
            try:
                output_root = self.project_service.require_output_root(create=False)
            except RuntimeError:
                output_root = None
            if output_root is None:
                self.preview_scroll.set_message("Select a project output directory before using CoSplice.")
                return
            candidate = output_root / "09_cosplice_network" / "cosplice_config.yaml"
            if candidate.exists():
                self.config_path.setText(str(candidate))
                path = candidate
        if path is None or not path.exists():
            self.status_table.setRowCount(0)
            self.plot_table.setRowCount(0)
            self.module_table.setRowCount(0)
            self.preview_scroll.set_message("No CoSplice config selected. Page refresh did not run analysis.")
            return
        try:
            config = load_config(path)
        except Exception as exc:
            self.log.setPlainText(f"Config error: {exc}")
            return
        source_index = self.matrix_source.findData(config.matrix_input.input_mode)
        if source_index >= 0:
            self.matrix_source.setCurrentIndex(source_index)
        if config.matrix_input.input_mode == "bam_junctions":
            self.matrix_source_status.setText(
                "BAM competing-junction usage: junction reads / reads in the donor-or-acceptor competitor cluster. "
                f"Backend: {config.matrix_input.backend}; BAM root: {config.bam_root or 'missing'}"
            )
        else:
            self.matrix_source_status.setText(
                "rMATS event usage: sample-level IncLevel/IJC/SJC from the selected rMATS result files."
            )
        state_path = config.output_dir / "state.json"
        state = json.loads(state_path.read_text(encoding="utf-8")) if state_path.exists() else {"stages": {}}
        cached_stages = state.get("stages", {})
        expected_stages = [spec.key for spec in stage_specs(visible_only=True)]
        stages = {
            stage: cached_stages.get(stage, {"status": "not_run", "updated_at": "", "message": "Not run", "outputs": []})
            for stage in expected_stages
        }
        self.status_table.setSortingEnabled(False)
        self.status_table.setRowCount(len(stages))
        for row, (name, details) in enumerate(stages.items()):
            values = [stage_spec(name).title, details.get("status", ""), details.get("updated_at", ""), details.get("message", ""), "\n".join(details.get("outputs", []))]
            for column, value in enumerate(values):
                self.status_table.setItem(row, column, QTableWidgetItem(str(value)))
        self.status_table.setSortingEnabled(True)
        plots = sorted(config.output_dir.glob("**/*.png"))
        self.plot_table.setSortingEnabled(False)
        self.plot_table.setRowCount(len(plots))
        for row, plot in enumerate(plots):
            self.plot_table.setItem(row, 0, QTableWidgetItem(plot.name))
            self.plot_table.setItem(row, 1, QTableWidgetItem(str(plot)))
        self.plot_table.setSortingEnabled(True)
        if plots and self.plot_table.currentRow() < 0:
            self.plot_table.selectRow(0)
            self._show_pixmap(plots[0])
        elif not plots:
            self.preview_scroll.set_message(
                "No cached CoSplice plot result. Run an individual network/association stage first."
            )
        self._populate_module_catalog(config.output_dir)

    def _populate_module_catalog(self, output_dir: Path) -> None:
        frame = self.project_service.cosplice_module_catalog(output_dir)
        if frame.empty:
            self.module_table.setColumnCount(1)
            self.module_table.setHorizontalHeaderLabels(["Status"])
            self.module_table.setRowCount(1)
            self.module_table.setItem(0, 0, QTableWidgetItem("No cached CoSplice modules are available."))
            return
        columns = list(frame.columns)
        self.module_table.setSortingEnabled(False)
        self.module_table.setColumnCount(len(columns))
        self.module_table.setHorizontalHeaderLabels(columns)
        self.module_table.setRowCount(len(frame))
        for row_index, (_, row) in enumerate(frame.iterrows()):
            for column_index, column in enumerate(columns):
                value = row[column]
                self.module_table.setItem(
                    row_index,
                    column_index,
                    QTableWidgetItem("" if value is None else str(value)),
                )
        self.module_table.setSortingEnabled(True)
        self.module_table.resizeColumnsToContents()

    def _publish_selected_modules(self) -> None:
        path = self._resolved_config_path()
        if path is None or not path.exists():
            QMessageBox.information(self, "Missing config", "Choose a CoSplice config first.")
            return
        selected_rows = sorted({index.row() for index in self.module_table.selectionModel().selectedRows()})
        module_column = next(
            (
                column
                for column in range(self.module_table.columnCount())
                if self.module_table.horizontalHeaderItem(column)
                and self.module_table.horizontalHeaderItem(column).text() == "module"
            ),
            None,
        )
        module_ids = [
            self.module_table.item(row, module_column).text()
            for row in selected_rows
            if module_column is not None and self.module_table.item(row, module_column) is not None
        ]
        if not module_ids:
            QMessageBox.information(self, "No modules selected", "Select one or more module rows first.")
            return
        try:
            config = load_config(path)
            selected = self.project_service.publish_cosplice_module_selection(config.output_dir, module_ids)
        except Exception as exc:
            QMessageBox.critical(self, "Module handoff failed", str(exc))
            return
        self.module_handoff_status.setText(
            f"Published {len(selected)} module(s) to AnalysisArtifactRegistry for downstream RBP, isoform and reporting."
        )

    def _preview_plot(self, row: int, _column: int, _previous_row: int, _previous_column: int) -> None:
        item = self.plot_table.item(row, 1) if row >= 0 else None
        if item:
            self._show_pixmap(Path(item.text()))

    def _show_pixmap(self, path: Path) -> None:
        self.preview_scroll.set_image(path)

    def _selected_plot(self) -> Path | None:
        row = self.plot_table.currentRow()
        item = self.plot_table.item(row, 1) if row >= 0 else None
        return Path(item.text()) if item else None

    def _export_plot(self) -> None:
        source = self._selected_plot()
        if source is None or not source.exists():
            QMessageBox.information(self, "No plot", "Select a generated PNG first.")
            return
        target, _ = QFileDialog.getSaveFileName(self, "Export plot", source.name, "PNG (*.png)")
        if target:
            shutil.copy2(source, target)

    def _open_output(self) -> None:
        path = self._resolved_config_path()
        if path is None or not path.exists():
            return
        output = load_config(path).output_dir
        output.mkdir(parents=True, exist_ok=True)
        os.startfile(str(output))

    def _show_detail(self) -> None:
        dialog = QDialog(self)
        dialog.setWindowTitle("CoSplice Network detail")
        dialog.setWindowFlag(Qt.WindowType.WindowMaximizeButtonHint, True)
        layout = QVBoxLayout(dialog)
        text = QTextEdit()
        text.setReadOnly(True)
        text.setPlainText(module_detail("9.cosplice"))
        close = QPushButton("Close")
        close.clicked.connect(dialog.accept)
        layout.addWidget(text, 1)
        layout.addWidget(close, 0, Qt.AlignmentFlag.AlignRight)
        dialog.resize(980, 720)
        dialog.exec()
