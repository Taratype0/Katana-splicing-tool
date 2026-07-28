from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QEvent, QSize, QTimer, Qt
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import QLabel, QScrollArea, QWidget


class ArtifactImagePreview(QScrollArea):
    """Reusable cache-only image preview with automatic fit and relayout."""

    def __init__(self, message: str = "No preview selected.", parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.label = QLabel(message)
        self.label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.label.setScaledContents(False)
        self.setWidgetResizable(True)
        self.setWidget(self.label)
        self._pixmap = QPixmap()
        self.viewport().installEventFilter(self)

    def set_image(self, path: str | Path) -> bool:
        pixmap = QPixmap(str(path))
        if pixmap.isNull():
            self.set_message(f"Preview unavailable: {path}")
            return False
        self._pixmap = pixmap
        self.label.setText("")
        QTimer.singleShot(0, self.fit_to_viewport)
        return True

    def set_message(self, message: str) -> None:
        self._pixmap = QPixmap()
        self.label.clear()
        self.label.setText(message)
        self.label.setAlignment(Qt.AlignmentFlag.AlignCenter)

    def fit_to_viewport(self) -> None:
        if self._pixmap.isNull():
            return
        viewport_size = self.viewport().size()
        target = QSize(max(100, viewport_size.width() - 24), max(100, viewport_size.height() - 24))
        scaled = self._pixmap.scaled(
            target,
            Qt.AspectRatioMode.KeepAspectRatio,
            Qt.TransformationMode.SmoothTransformation,
        )
        self.label.setPixmap(scaled)
        self.label.resize(scaled.size())
        self.label.setAlignment(Qt.AlignmentFlag.AlignCenter)

    def eventFilter(self, watched, event) -> bool:
        if watched is self.viewport() and event.type() in {QEvent.Type.Resize, QEvent.Type.Show}:
            QTimer.singleShot(0, self.fit_to_viewport)
        return super().eventFilter(watched, event)
