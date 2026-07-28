from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QFrame,
    QGridLayout,
    QLabel,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)


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

        header = QFrame()
        header.setObjectName("pageHeader")
        header_layout = QVBoxLayout(header)
        header_layout.setContentsMargins(22, 20, 22, 20)
        header_layout.setSpacing(6)

        overline = QLabel("WORKFLOW STAGE")
        overline.setObjectName("pageOverline")
        title_label = QLabel(title)
        title_label.setObjectName("pageTitle")
        description_label = QLabel(description)
        description_label.setObjectName("pageDescription")
        description_label.setWordWrap(True)
        header_layout.addWidget(overline)
        header_layout.addWidget(title_label)
        header_layout.addWidget(description_label)

        grid = QGridLayout()
        grid.setContentsMargins(0, 6, 0, 0)
        grid.setHorizontalSpacing(16)
        grid.setVerticalSpacing(16)

        for index, (label, body, target) in enumerate(sections):
            card = QFrame()
            card.setObjectName("navigationCard")
            card.setMinimumHeight(178)
            card.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            card_layout = QVBoxLayout(card)
            card_layout.setContentsMargins(18, 16, 18, 16)
            card_layout.setSpacing(9)

            route_number = target.split(" ", 1)[0]
            route = QLabel(route_number)
            route.setObjectName("routePill")

            heading = QLabel(label)
            heading.setObjectName("cardTitle")
            heading.setWordWrap(True)

            text = QLabel(body)
            text.setObjectName("cardBody")
            text.setWordWrap(True)

            button = QPushButton("Open workspace")
            button.setProperty("class", "primary")
            button.clicked.connect(lambda _checked=False, route=target: self._open(route))

            card_layout.addWidget(route, 0, Qt.AlignmentFlag.AlignLeft)
            card_layout.addWidget(heading)
            card_layout.addWidget(text)
            card_layout.addStretch(1)
            card_layout.addWidget(button)
            if self.columns == 1:
                grid.addWidget(card, index, 0)
            else:
                grid.addWidget(card, index // self.columns, index % self.columns)
        for column in range(self.columns):
            grid.setColumnStretch(column, 1)

        content = QWidget()
        content_layout = QVBoxLayout(content)
        content_layout.setContentsMargins(2, 2, 10, 18)
        content_layout.setSpacing(16)
        content_layout.addWidget(header)
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
