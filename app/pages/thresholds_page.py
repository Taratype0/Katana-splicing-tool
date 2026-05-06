from __future__ import annotations

from PySide6.QtWidgets import (
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QLabel,
    QVBoxLayout,
    QWidget,
)

from src.services.project_service import ProjectService


class ThresholdsPage(QWidget):
    def __init__(self, project_service: ProjectService) -> None:
        super().__init__()
        self.project_service = project_service
        self._is_refreshing = False
        self.label = QLabel("Analysis thresholds")
        self.preset_labels = {
            "default": "Literature standard",
            "relaxed": "Exploratory",
            "strict": "Stringent",
        }
        self.preset_keys_by_index: list[str] = []
        self.preset = QComboBox()
        for key, label in self.preset_labels.items():
            self.preset.addItem(label, key)
            self.preset_keys_by_index.append(key)
        self.preset.currentTextChanged.connect(self._apply_preset)
        self.note = QLabel(
            "Literature standard preset: rMATS-style AS uses FDR 0.05 and |dPSI| 0.10; "
            "cross-comparison program-shift uses |ΔdPSI| 0.10 as the default strength cutoff; "
            "DESeq2 keeps adjusted p 0.05 with a biologically meaningful |log2FC| 1.0 default; "
            "DEXSeq is shown at q 0.10; stageR/DTU remains at 0.05."
        )
        self.note.setWordWrap(True)
        self.usage_note = QLabel(
            "These thresholds are used when you run analysis modules. Changing them affects future reruns and candidate-screening outputs; "
            "it does not silently rewrite existing cached result files until that module is rerun."
        )
        self.usage_note.setWordWrap(True)

        self.splicing_fdr = self._spin(0.0, 1.0, 0.01)
        self.splicing_dpsi = self._spin(0.0, 5.0, 0.01)
        self.program_delta_dpsi = self._spin(0.0, 5.0, 0.01)
        self.deg_padj = self._spin(0.0, 1.0, 0.01)
        self.deg_log2fc = self._spin(0.0, 10.0, 0.1)
        self.dexseq_qvalue = self._spin(0.0, 1.0, 0.01)
        self.dtu_qvalue = self._spin(0.0, 1.0, 0.01)
        for spin in (
            self.splicing_fdr,
            self.splicing_dpsi,
            self.program_delta_dpsi,
            self.deg_padj,
            self.deg_log2fc,
            self.dexseq_qvalue,
            self.dtu_qvalue,
        ):
            spin.valueChanged.connect(self._apply_custom_values)

        form = QFormLayout()
        form.addRow("Preset", self.preset)
        form.addRow("Splicing FDR", self.splicing_fdr)
        form.addRow("|dPSI|", self.splicing_dpsi)
        form.addRow("Program |ΔdPSI|", self.program_delta_dpsi)
        form.addRow("DEG padj", self.deg_padj)
        form.addRow("|log2FC|", self.deg_log2fc)
        form.addRow("DEXSeq qvalue", self.dexseq_qvalue)
        form.addRow("DTU qvalue", self.dtu_qvalue)

        group = QGroupBox("Thresholds")
        group.setLayout(form)

        layout = QVBoxLayout(self)
        layout.addWidget(self.label)
        layout.addWidget(group)
        layout.addWidget(self.note)
        layout.addWidget(self.usage_note)
        layout.addStretch(1)

    def refresh(self) -> None:
        thresholds = self.project_service.current_thresholds
        self._is_refreshing = True
        self.splicing_fdr.setValue(thresholds.splicing_fdr)
        self.splicing_dpsi.setValue(thresholds.splicing_dpsi)
        self.program_delta_dpsi.setValue(thresholds.program_delta_dpsi)
        self.deg_padj.setValue(thresholds.deg_padj)
        self.deg_log2fc.setValue(thresholds.deg_log2fc)
        self.dexseq_qvalue.setValue(thresholds.dexseq_qvalue)
        self.dtu_qvalue.setValue(thresholds.dtu_qvalue)
        self._update_note()
        self._is_refreshing = False

    def _apply_preset(self, _preset_label: str) -> None:
        if self.project_service.current_project is None:
            return
        preset_name = self.preset.currentData()
        self.project_service.apply_threshold_preset(preset_name)
        self.refresh()

    def _apply_custom_values(self) -> None:
        if self._is_refreshing or self.project_service.current_project is None:
            return
        self.project_service.update_thresholds(
            splicing_fdr=self.splicing_fdr.value(),
            splicing_dpsi=self.splicing_dpsi.value(),
            program_delta_dpsi=self.program_delta_dpsi.value(),
            deg_padj=self.deg_padj.value(),
            deg_log2fc=self.deg_log2fc.value(),
            dexseq_qvalue=self.dexseq_qvalue.value(),
            dtu_qvalue=self.dtu_qvalue.value(),
        )
        self._update_note()

    def _spin(self, low: float, high: float, step: float) -> QDoubleSpinBox:
        spin = QDoubleSpinBox()
        spin.setRange(low, high)
        spin.setSingleStep(step)
        spin.setDecimals(4)
        spin.setKeyboardTracking(False)
        return spin

    def _update_note(self) -> None:
        key = self.preset.currentData() or "default"
        notes = {
            "default": (
                "Literature standard: AS FDR 0.05 and |dPSI| 0.10; "
                "program-strength cutoff |ΔdPSI| 0.10; "
                "DEG adjusted p 0.05 with |log2FC| 1.0; "
                "DEXSeq q 0.10; DTU/stageR q 0.05."
            ),
            "relaxed": (
                "Exploratory: looser cutoffs for hypothesis generation "
                "(AS/DEG q up to 0.10, smaller effect-size thresholds)."
            ),
            "strict": (
                "Stringent: tighter significance and larger effect-size cutoffs "
                "for higher-confidence shortlist generation."
            ),
        }
        self.note.setText(notes.get(key, notes["default"]))
