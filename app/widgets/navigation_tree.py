from __future__ import annotations

import re

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QFont
from PySide6.QtWidgets import QTreeWidget, QTreeWidgetItem


_ROUTE_PATTERN = re.compile(r"^(?P<number>\d+(?:\.\d+)*)(?:\.)?\s+(?P<title>.+)$")


class WorkflowNavigation(QTreeWidget):
    """Hierarchical navigation that keeps route indexes compatible with the legacy page router."""

    currentRowChanged = Signal(int)

    def __init__(self) -> None:
        super().__init__()
        self.setObjectName("workflowNavigation")
        self.setHeaderHidden(True)
        self.setIndentation(16)
        self.setAnimated(True)
        self.setUniformRowHeights(True)
        self.setRootIsDecorated(True)
        self.setExpandsOnDoubleClick(False)
        self._route_items: dict[int, QTreeWidgetItem] = {}
        self._path_items: dict[str, QTreeWidgetItem] = {}
        self.currentItemChanged.connect(self._emit_route_changed)

    def add_route(self, route_index: int, name: str) -> None:
        match = _ROUTE_PATTERN.match(name)
        number = match.group("number") if match else str(route_index + 1)
        title = match.group("title") if match else name
        parent_path = number.rsplit(".", 1)[0] if "." in number else ""
        parent = self._path_items.get(parent_path)

        label = f"{number}  {title}"
        item = QTreeWidgetItem([label])
        item.setData(0, Qt.ItemDataRole.UserRole, route_index)
        item.setData(0, Qt.ItemDataRole.UserRole + 1, name)
        item.setToolTip(0, name)
        if parent is None:
            self.addTopLevelItem(item)
            font = QFont(item.font(0))
            font.setBold(True)
            item.setFont(0, font)
        else:
            parent.addChild(item)
        self._route_items[route_index] = item
        self._path_items[number] = item

    def setCurrentRow(self, route_index: int) -> None:  # noqa: N802 - QListWidget compatibility
        item = self._route_items.get(route_index)
        if item is None:
            return
        parent = item.parent()
        while parent is not None:
            parent.setExpanded(True)
            parent = parent.parent()
        self.setCurrentItem(item)
        self.scrollToItem(item)

    def currentRow(self) -> int:  # noqa: N802 - QListWidget compatibility
        item = self.currentItem()
        if item is None:
            return -1
        value = item.data(0, Qt.ItemDataRole.UserRole)
        return int(value) if value is not None else -1

    def expand_primary_workflow(self) -> None:
        for path in ("5", "7", "8"):
            item = self._path_items.get(path)
            if item is not None:
                item.setExpanded(path == "5")

    def _emit_route_changed(self, current: QTreeWidgetItem | None, _previous: QTreeWidgetItem | None) -> None:
        if current is None:
            return
        value = current.data(0, Qt.ItemDataRole.UserRole)
        if value is not None:
            self.currentRowChanged.emit(int(value))
