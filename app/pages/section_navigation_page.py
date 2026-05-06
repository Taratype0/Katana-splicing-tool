from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QGridLayout, QLabel, QPushButton, QScrollArea, QVBoxLayout, QWidget


class SectionNavigationPage(QWidget):
    def __init__(
        self,
        title: str,
        description: str,
        sections: list[tuple[str, str, str]],
        on_open_section=None,
        columns: int = 2,
    ) -> None:
        super().__init__()
        self.on_open_section = on_open_section
        self.columns = max(1, int(columns))

        title_label = QLabel(title)
        title_label.setStyleSheet("font-size: 18px; font-weight: 700;")

        description_label = QLabel(description)
        description_label.setWordWrap(True)

        grid = QGridLayout()
        grid.setHorizontalSpacing(14)
        grid.setVerticalSpacing(14)

        for index, (label, body, target) in enumerate(sections):
            card = QWidget()
            card_layout = QVBoxLayout(card)
            card_layout.setContentsMargins(12, 12, 12, 12)

            heading = QLabel(label)
            heading.setStyleSheet("font-size: 15px; font-weight: 700;")

            text = QLabel(body)
            text.setWordWrap(True)

            button = QPushButton(f"Open {label}")
            button.clicked.connect(lambda _checked=False, route=target: self._open(route))

            card_layout.addWidget(heading)
            card_layout.addWidget(text)
            card_layout.addStretch(1)
            card_layout.addWidget(button)
            card.setStyleSheet("QWidget { border: 1px solid #3a3a3a; border-radius: 8px; }")
            if self.columns == 1:
                grid.addWidget(card, index, 0)
            else:
                grid.addWidget(card, index // self.columns, index % self.columns)
        for column in range(self.columns):
            grid.setColumnStretch(column, 1)

        content = QWidget()
        content_layout = QVBoxLayout(content)
        content_layout.setContentsMargins(0, 0, 0, 0)
        content_layout.addWidget(title_label)
        content_layout.addWidget(description_label)
        content_layout.addLayout(grid)
        content_layout.addStretch(1)

        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setWidget(content)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(scroll)

    def _open(self, target: str) -> None:
        if callable(self.on_open_section):
            self.on_open_section(target)
