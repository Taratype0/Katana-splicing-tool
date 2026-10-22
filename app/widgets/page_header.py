from __future__ import annotations

from PySide6.QtWidgets import QFrame, QLabel, QVBoxLayout


class PageHeader(QFrame):
    """Shared compact header for workflow pages."""

    def __init__(self, title: str, description: str, overline: str) -> None:
        super().__init__()
        self.setObjectName("pageHeader")
        layout = QVBoxLayout(self)
        layout.setContentsMargins(20, 16, 20, 16)
        layout.setSpacing(5)

        stage = QLabel(overline.upper())
        stage.setObjectName("pageOverline")
        heading = QLabel(title)
        heading.setObjectName("pageTitle")
        body = QLabel(description)
        body.setObjectName("pageDescription")
        body.setWordWrap(True)

        layout.addWidget(stage)
        layout.addWidget(heading)
        layout.addWidget(body)
