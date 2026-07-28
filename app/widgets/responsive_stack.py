from __future__ import annotations

from PySide6.QtCore import QSize
from PySide6.QtWidgets import QSizePolicy, QStackedWidget


class ResponsiveStackedWidget(QStackedWidget):
    """Prevent hidden wide pages from forcing the entire application off-screen."""

    def __init__(self) -> None:
        super().__init__()
        self.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Ignored)

    def sizeHint(self) -> QSize:  # noqa: N802 - Qt API
        return QSize(1040, 720)

    def minimumSizeHint(self) -> QSize:  # noqa: N802 - Qt API
        return QSize(720, 480)
