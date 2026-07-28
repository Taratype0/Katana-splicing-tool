from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import QSize
from PySide6.QtWidgets import QApplication, QWidget

from app.widgets.responsive_stack import ResponsiveStackedWidget


class _VeryWidePage(QWidget):
    def minimumSizeHint(self) -> QSize:  # noqa: N802 - Qt API
        return QSize(4000, 900)


def test_hidden_wide_page_does_not_force_workspace_width() -> None:
    _app = QApplication.instance() or QApplication([])
    stack = ResponsiveStackedWidget()
    stack.addWidget(_VeryWidePage())

    assert stack.minimumSizeHint().width() == 720
    assert stack.sizeHint().width() == 1040
