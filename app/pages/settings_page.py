from __future__ import annotations

from collections.abc import Callable

from PySide6.QtWidgets import (
    QComboBox,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPushButton,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from src.services.project_service import ProjectService


class SettingsPage(QWidget):
    def __init__(self, project_service: ProjectService, reload_callback: Callable[[], None] | None = None) -> None:
        super().__init__()
        self.project_service = project_service
        self.reload_callback = reload_callback
        self.label = QLabel("Tool configuration")
        self.rmats_mode = QComboBox()
        self.rmats_mode.addItems(["JC", "JCEC"])
        self.rmats_mode.currentTextChanged.connect(self._update_rmats_mode)
        self.rescan_button = QPushButton("Reload Using Current Input Paths")
        self.rescan_button.clicked.connect(self._reload_project)
        self.output_label = QLabel("Output directory: not set")
        self.output_button = QPushButton("Choose Output Directory")
        self.output_button.clicked.connect(self._choose_output_directory)
        self.input_rmats_label = QLabel("Input rMATS root: auto")
        self.input_rmats_button = QPushButton("Choose rMATS Root")
        self.input_rmats_button.clicked.connect(lambda: self._choose_input_directory("rmats_root", "Select rMATS Root"))
        self.input_deg_label = QLabel("Input DEG root: auto")
        self.input_deg_button = QPushButton("Choose DEG Root")
        self.input_deg_button.clicked.connect(lambda: self._choose_input_directory("deg_root", "Select DEG Root"))
        self.input_suppa_label = QLabel("Input SUPPA root: auto")
        self.input_suppa_button = QPushButton("Choose SUPPA Root")
        self.input_suppa_button.clicked.connect(lambda: self._choose_input_directory("suppa_root", "Select SUPPA Root"))
        self.input_dexseq_label = QLabel("Input DEXSeq root: auto")
        self.input_dexseq_button = QPushButton("Choose DEXSeq Root")
        self.input_dexseq_button.clicked.connect(lambda: self._choose_input_directory("dexseq_root", "Select DEXSeq Root"))
        self.input_dtu_label = QLabel("Input DTU root: auto")
        self.input_dtu_button = QPushButton("Choose DTU Root")
        self.input_dtu_button.clicked.connect(lambda: self._choose_input_directory("dtu_root", "Select DTU Root"))
        self.input_quant_label = QLabel("Input quant root: auto")
        self.input_quant_button = QPushButton("Choose quant Root")
        self.input_quant_button.clicked.connect(lambda: self._choose_input_directory("quant_root", "Select quant Root"))
        self.input_counts_label = QLabel("Input counts file: auto")
        self.input_counts_button = QPushButton("Choose Counts File")
        self.input_counts_button.clicked.connect(lambda: self._choose_input_file("counts_path", "Select counts file", "TSV Files (*.tsv *.txt *.csv);;All Files (*)"))
        self.input_contrast_label = QLabel("Input contrastsheet: auto")
        self.input_contrast_button = QPushButton("Choose contrastsheet")
        self.input_contrast_button.clicked.connect(lambda: self._choose_input_file("contrastsheet_path", "Select contrastsheet", "CSV Files (*.csv *.tsv);;All Files (*)"))
        self.rscript_label = QLabel("Rscript: not set")
        self.rscript_button = QPushButton("Choose Rscript")
        self.rscript_button.clicked.connect(self._choose_rscript)
        self.jutils_label = QLabel("Jutils: not set")
        self.jutils_button = QPushButton("Choose Jutils Directory")
        self.jutils_button.clicked.connect(self._choose_jutils)
        self.sashimi_label = QLabel("rmats2sashimiplot: not set")
        self.sashimi_button = QPushButton("Choose rmats2sashimiplot Directory")
        self.sashimi_button.clicked.connect(self._choose_sashimi)
        self.suppa_label = QLabel("SUPPA: not set")
        self.suppa_button = QPushButton("Choose SUPPA Directory")
        self.suppa_button.clicked.connect(self._choose_suppa)
        self.bam_label = QLabel("BAM root: not set")
        self.bam_button = QPushButton("Choose BAM Root")
        self.bam_button.clicked.connect(self._choose_bam_root)
        self.gtf_label = QLabel("GTF: not set")
        self.gtf_button = QPushButton("Choose GTF")
        self.gtf_button.clicked.connect(self._choose_gtf)
        self.fasta_label = QLabel("FASTA: not set")
        self.fasta_button = QPushButton("Choose FASTA")
        self.fasta_button.clicked.connect(self._choose_fasta)
        self.details = QTextEdit()
        self.details.setReadOnly(True)
        layout = QVBoxLayout(self)
        layout.addWidget(self.label)
        layout.addWidget(self.rmats_mode)
        layout.addWidget(self.rescan_button)
        row = QHBoxLayout()
        row.addWidget(self.output_label)
        row.addWidget(self.output_button)
        layout.addLayout(row)
        self._add_row(layout, self.input_rmats_label, self.input_rmats_button)
        self._add_row(layout, self.input_deg_label, self.input_deg_button)
        self._add_row(layout, self.input_suppa_label, self.input_suppa_button)
        self._add_row(layout, self.input_dexseq_label, self.input_dexseq_button)
        self._add_row(layout, self.input_dtu_label, self.input_dtu_button)
        self._add_row(layout, self.input_quant_label, self.input_quant_button)
        self._add_row(layout, self.input_counts_label, self.input_counts_button)
        self._add_row(layout, self.input_contrast_label, self.input_contrast_button)
        row2 = QHBoxLayout()
        row2.addWidget(self.rscript_label)
        row2.addWidget(self.rscript_button)
        layout.addLayout(row2)
        row_tools_1 = QHBoxLayout()
        row_tools_1.addWidget(self.jutils_label)
        row_tools_1.addWidget(self.jutils_button)
        layout.addLayout(row_tools_1)
        row_tools_2 = QHBoxLayout()
        row_tools_2.addWidget(self.sashimi_label)
        row_tools_2.addWidget(self.sashimi_button)
        layout.addLayout(row_tools_2)
        row_tools_3 = QHBoxLayout()
        row_tools_3.addWidget(self.suppa_label)
        row_tools_3.addWidget(self.suppa_button)
        layout.addLayout(row_tools_3)
        row_tools_4 = QHBoxLayout()
        row_tools_4.addWidget(self.bam_label)
        row_tools_4.addWidget(self.bam_button)
        layout.addLayout(row_tools_4)
        row3 = QHBoxLayout()
        row3.addWidget(self.gtf_label)
        row3.addWidget(self.gtf_button)
        layout.addLayout(row3)
        row4 = QHBoxLayout()
        row4.addWidget(self.fasta_label)
        row4.addWidget(self.fasta_button)
        layout.addLayout(row4)
        layout.addWidget(self.details)

    def refresh(self) -> None:
        project = self.project_service.current_project
        if project is None:
            self.details.setPlainText(
                "No project loaded.\n\n"
                "Use 'Open Project' or 'Use Current Folder' in the main toolbar first.\n"
                "After that, you can configure Rscript, GTF, BAM root, and external tools here."
            )
            return
        self.rmats_mode.setCurrentText(project.selected_rmats_mode)
        self.output_label.setText(f"Output directory: {project.output_root}")
        self.input_rmats_label.setText(f"Input rMATS root: {project.input_paths.get('rmats_root', 'auto')}")
        self.input_deg_label.setText(f"Input DEG root: {project.input_paths.get('deg_root', 'auto')}")
        self.input_suppa_label.setText(f"Input SUPPA root: {project.input_paths.get('suppa_root', 'auto')}")
        self.input_dexseq_label.setText(f"Input DEXSeq root: {project.input_paths.get('dexseq_root', 'auto')}")
        self.input_dtu_label.setText(f"Input DTU root: {project.input_paths.get('dtu_root', 'auto')}")
        self.input_quant_label.setText(f"Input quant root: {project.input_paths.get('quant_root', 'auto')}")
        self.input_counts_label.setText(f"Input counts file: {project.input_paths.get('counts_path', 'auto')}")
        self.input_contrast_label.setText(f"Input contrastsheet: {project.input_paths.get('contrastsheet_path', 'auto')}")
        self.rscript_label.setText(f"Rscript: {project.tool_paths.get('rscript', 'not set')}")
        self.jutils_label.setText(f"Jutils: {project.tool_paths.get('jutils', 'not set')}")
        self.sashimi_label.setText(f"rmats2sashimiplot: {project.tool_paths.get('rmats2sashimiplot', 'not set')}")
        self.suppa_label.setText(f"SUPPA: {project.tool_paths.get('suppa', 'not set')}")
        self.bam_label.setText(f"BAM root: {project.tool_paths.get('bam_root', 'not set')}")
        self.gtf_label.setText(f"GTF: {project.tool_paths.get('gtf', 'not set')}")
        self.fasta_label.setText(f"FASTA: {project.tool_paths.get('fasta', 'not set')}")
        self.details.setPlainText(
            "\n".join(
                [
                    "Recommended manual input targets:",
                    "- rMATS root: folder that directly contains comparison subfolders with rmats_post/",
                    "- DEG root: folder that contains *.deseq2.results.tsv",
                    "- SUPPA root: folder that contains *_local_diffsplice.dpsi",
                    "- DEXSeq root: folder that contains perGeneQValue.*.csv",
                    "- DTU root: folder that contains DEXSeqResults.*.tsv and getAdjustedPValues.*.tsv",
                    "- quant root: folder above per-sample quant.sf directories",
                    "- counts file: all.normalised_counts.tsv",
                    "- contrastsheet: contrastsheet.valid.csv",
                    "",
                    f"Jutils path: {project.tool_paths.get('jutils', 'not set')}",
                    f"rmats2sashimiplot path: {project.tool_paths.get('rmats2sashimiplot', 'not set')}",
                    f"SUPPA path: {project.tool_paths.get('suppa', 'not set')}",
                    f"BAM root: {project.tool_paths.get('bam_root', 'not set')}",
                    f"Output directory: {project.output_root}",
                    "",
                    "Manual input paths:",
                    f"rMATS root: {project.input_paths.get('rmats_root', 'auto')}",
                    f"DEG root: {project.input_paths.get('deg_root', 'auto')}",
                    f"SUPPA root: {project.input_paths.get('suppa_root', 'auto')}",
                    f"DEXSeq root: {project.input_paths.get('dexseq_root', 'auto')}",
                    f"DTU root: {project.input_paths.get('dtu_root', 'auto')}",
                    f"quant root: {project.input_paths.get('quant_root', 'auto')}",
                    f"counts file: {project.input_paths.get('counts_path', 'auto')}",
                    f"contrastsheet: {project.input_paths.get('contrastsheet_path', 'auto')}",
                    "",
                    "Runtime status:",
                    *self.project_service.runtime_status_lines(),
                    "",
                    "Scan messages:",
                    *project.scan_messages,
                ]
            )
        )

    def _update_rmats_mode(self, value: str) -> None:
        if self.project_service.current_project is None:
            return
        self.project_service.update_rmats_mode(value)

    def _choose_output_directory(self) -> None:
        project = self.project_service.current_project
        if project is None:
            return
        folder = QFileDialog.getExistingDirectory(
            self,
            "Select Output Directory",
            str(project.output_root or project.project_root),
        )
        if not folder:
            return
        self.project_service.update_output_root(folder)
        self.refresh()

    def _choose_rscript(self) -> None:
        project = self.project_service.current_project
        if project is None:
            return
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Rscript executable",
            str(project.project_root),
            "Executables (*.exe);;All Files (*)",
        )
        if not file_path:
            return
        project.tool_paths["rscript"] = file_path
        self.project_service.save_project_config()
        self.project_service._refresh_tool_adapters()
        self.refresh()

    def _choose_jutils(self) -> None:
        self._choose_tool_directory("jutils", "Select Jutils Directory")

    def _choose_sashimi(self) -> None:
        self._choose_tool_directory("rmats2sashimiplot", "Select rmats2sashimiplot Directory")

    def _choose_suppa(self) -> None:
        self._choose_tool_directory("suppa", "Select SUPPA Directory")

    def _choose_bam_root(self) -> None:
        self._choose_tool_directory("bam_root", "Select BAM Root Directory")

    def _choose_gtf(self) -> None:
        project = self.project_service.current_project
        if project is None:
            return
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select GTF file",
            str(project.project_root),
            "GTF Files (*.gtf *.gtf.gz);;All Files (*)",
        )
        if not file_path:
            return
        project.tool_paths["gtf"] = file_path
        self.project_service.save_project_config()
        self.refresh()

    def _choose_fasta(self) -> None:
        project = self.project_service.current_project
        if project is None:
            return
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select FASTA file",
            str(project.project_root),
            "FASTA Files (*.fa *.fasta *.fa.gz *.fasta.gz);;All Files (*)",
        )
        if not file_path:
            return
        project.tool_paths["fasta"] = file_path
        self.project_service.save_project_config()
        self.refresh()

    def _choose_tool_directory(self, key: str, title: str) -> None:
        project = self.project_service.current_project
        if project is None:
            QMessageBox.information(self, "Open project first", "Load a project before configuring tool paths.")
            return
        folder = QFileDialog.getExistingDirectory(
            self,
            title,
            str(project.project_root),
        )
        if not folder:
            return
        project.tool_paths[key] = folder
        self.project_service.save_project_config()
        self.project_service._refresh_tool_adapters()
        self.refresh()

    def _choose_input_directory(self, key: str, title: str) -> None:
        project = self.project_service.current_project
        if project is None:
            QMessageBox.information(self, "Open project first", "Load a project before configuring input paths.")
            return
        folder = QFileDialog.getExistingDirectory(
            self,
            title,
            str(project.project_root),
        )
        if not folder:
            return
        self.project_service.update_input_path(key, folder)
        self.refresh()

    def _choose_input_file(self, key: str, title: str, filter_text: str) -> None:
        project = self.project_service.current_project
        if project is None:
            QMessageBox.information(self, "Open project first", "Load a project before configuring input paths.")
            return
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            title,
            str(project.project_root),
            filter_text,
        )
        if not file_path:
            return
        self.project_service.update_input_path(key, file_path)
        self.refresh()

    def _reload_project(self) -> None:
        if self.reload_callback is None:
            return
        self.reload_callback()

    def _add_row(self, layout: QVBoxLayout, label: QLabel, button: QPushButton) -> None:
        row = QHBoxLayout()
        row.addWidget(label)
        row.addWidget(button)
        layout.addLayout(row)
