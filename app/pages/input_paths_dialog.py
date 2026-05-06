from __future__ import annotations

from pathlib import Path

import yaml
from PySide6.QtWidgets import (
    QDialog,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QVBoxLayout,
)


class InputPathsDialog(QDialog):
    FIELD_SPECS = [
        ("rmats_root", "rMATS root", "dir"),
        ("deg_root", "DEG root", "dir"),
        ("suppa_root", "SUPPA root", "dir"),
        ("dexseq_root", "DEXSeq root", "dir"),
        ("dtu_root", "DTU root", "dir"),
        ("quant_root", "quant root", "dir"),
        ("counts_path", "counts file", "file"),
        ("contrastsheet_path", "contrastsheet", "file"),
    ]

    def __init__(self, project_root: Path) -> None:
        super().__init__()
        self.project_root = Path(project_root)
        self.setWindowTitle("Choose Input Locations")
        self.resize(900, 420)
        self.inputs: dict[str, QLineEdit] = {}
        self.force_auto_scan = False

        intro = QLabel(
            "Select the subfolders/files to use for this project. "
            "Leave any field blank to let Katana auto-detect that input."
        )
        intro.setWordWrap(True)

        form = QFormLayout()
        saved = self._load_saved_paths()
        for key, label, kind in self.FIELD_SPECS:
            line = QLineEdit()
            line.setPlaceholderText("auto-detect")
            if saved.get(key):
                line.setText(saved[key])
            button = QPushButton("Browse")
            if kind == "dir":
                button.clicked.connect(lambda _=False, k=key: self._choose_directory(k))
            else:
                button.clicked.connect(lambda _=False, k=key: self._choose_file(k))
            row = QHBoxLayout()
            row.addWidget(line, 1)
            row.addWidget(button)
            wrapper = QVBoxLayout()
            wrapper.setContentsMargins(0, 0, 0, 0)
            wrapper.addLayout(row)
            hint = QLabel(self._hint_for(key))
            hint.setStyleSheet("color: #888; font-size: 11px;")
            hint.setWordWrap(True)
            wrapper.addWidget(hint)
            form.addRow(label, wrapper)
            self.inputs[key] = line

        use_auto = QPushButton("Use Auto Scan")
        use_auto.clicked.connect(self._use_auto_scan)
        cancel = QPushButton("Cancel")
        cancel.clicked.connect(self.reject)
        ok = QPushButton("Load Project")
        ok.clicked.connect(self.accept)

        buttons = QHBoxLayout()
        buttons.addWidget(use_auto)
        buttons.addStretch(1)
        buttons.addWidget(cancel)
        buttons.addWidget(ok)

        layout = QVBoxLayout(self)
        layout.addWidget(intro)
        layout.addLayout(form)
        layout.addLayout(buttons)

    def selected_paths(self) -> dict[str, str]:
        output: dict[str, str] = {}
        for key, line in self.inputs.items():
            value = line.text().strip()
            if value:
                output[key] = value
        return output

    def _choose_directory(self, key: str) -> None:
        folder = QFileDialog.getExistingDirectory(self, f"Select {key}", str(self.project_root))
        if folder:
            self.inputs[key].setText(folder)

    def _choose_file(self, key: str) -> None:
        file_path, _ = QFileDialog.getOpenFileName(self, f"Select {key}", str(self.project_root), "All Files (*)")
        if file_path:
            self.inputs[key].setText(file_path)

    def _clear_all(self) -> None:
        for line in self.inputs.values():
            line.clear()

    def _use_auto_scan(self) -> None:
        self.force_auto_scan = True
        self._clear_all()
        self.accept()

    def _load_saved_paths(self) -> dict[str, str]:
        config_path = self.project_root / ".katana_project.yaml"
        if not config_path.exists():
            return {}
        try:
            with config_path.open("r", encoding="utf-8") as handle:
                data = yaml.safe_load(handle) or {}
        except Exception:
            return {}
        input_paths = data.get("input_paths") or {}
        return {key: str(value) for key, value in input_paths.items() if value}

    def _hint_for(self, key: str) -> str:
        hints = {
            "rmats_root": "Folder that directly contains comparison subfolders with rmats_post/.",
            "deg_root": "Folder that contains *.deseq2.results.tsv files.",
            "suppa_root": "Folder that contains *_local_diffsplice.dpsi outputs.",
            "dexseq_root": "Folder that contains perGeneQValue.*.csv outputs.",
            "dtu_root": "Folder that contains DEXSeqResults.*.tsv and getAdjustedPValues.*.tsv.",
            "quant_root": "Folder above per-sample quant.sf directories.",
            "counts_path": "Usually all.normalised_counts.tsv.",
            "contrastsheet_path": "Usually contrastsheet.valid.csv.",
        }
        return hints.get(key, "")
