from __future__ import annotations


KATANA_STYLESHEET = """
QWidget {
    background: #0b141d;
    color: #eaf2f5;
    font-family: "Segoe UI";
    font-size: 10pt;
}
QMainWindow, QDialog {
    background: #071018;
}
QFrame#sidebar {
    background: #08121b;
    border-right: 1px solid #1d303c;
}
QFrame#contentSurface {
    background: #0b141d;
    border: 0;
}
QLabel#brandTitle {
    color: #f5fbfc;
    font-size: 22pt;
    font-weight: 800;
    letter-spacing: 3px;
    background: transparent;
}
QLabel#brandSubtitle {
    color: #7f98a5;
    font-size: 9.5pt;
    background: transparent;
}
QLabel#projectBadge {
    color: #bfeee8;
    background: #10242c;
    border: 1px solid #1d514f;
    border-radius: 9px;
    padding: 8px 10px;
    font-size: 9pt;
}
QLabel#sidebarOverline, QLabel#pageOverline {
    color: #52c9bd;
    font-size: 8pt;
    font-weight: 700;
    letter-spacing: 1.4px;
    background: transparent;
}
QFrame#pageHeader {
    background: #101e29;
    border: 1px solid #203746;
    border-radius: 15px;
}
QLabel#pageTitle {
    color: #f5fbfc;
    font-size: 21pt;
    font-weight: 760;
    padding: 0;
    background: transparent;
}
QLabel#pageDescription {
    color: #a8bac3;
    font-size: 10.5pt;
    line-height: 1.4;
    background: transparent;
}
QLabel#pageSubtitle {
    color: #9eb1bb;
    font-size: 10pt;
    background: transparent;
}
QFrame#analysisSummary {
    background: #0e1d27;
    border: 1px solid #24414e;
    border-radius: 12px;
}
QLabel#summaryLabel {
    color: #8099a5;
    font-size: 8.5pt;
    font-weight: 600;
    background: transparent;
}
QLabel#summaryValue {
    color: #eaf9f7;
    font-size: 16pt;
    font-weight: 760;
    background: transparent;
}
QLabel#contextPanel {
    color: #a9bbc4;
    background: #0d1b24;
    border: 1px solid #203946;
    border-radius: 10px;
    padding: 10px 12px;
}
QLabel#workflowContext {
    color: #b6c7cf;
    background: #0d2028;
    border: 1px solid #24504f;
    border-radius: 11px;
    padding: 11px 14px;
}
QFrame#workflowStageCard {
    background: #101d27;
    border: 1px solid #203745;
    border-radius: 14px;
    min-height: 190px;
}
QFrame#workflowStageCard:hover {
    background: #12242f;
    border-color: #337d78;
}
QLabel#workflowStageNumber {
    color: #071518;
    background: #69d8cd;
    border-radius: 11px;
    min-width: 22px;
    min-height: 22px;
    max-width: 22px;
    max-height: 22px;
    qproperty-alignment: AlignCenter;
    font-size: 9pt;
    font-weight: 800;
}
QLabel#handoffLabel {
    color: #78bdb8;
    background: #0c1820;
    border-left: 3px solid #278f88;
    padding: 7px 9px;
    font-size: 9pt;
}
QLabel#statusChip {
    color: #9fb2bc;
    background: #182832;
    border: 1px solid #2a424e;
    border-radius: 9px;
    padding: 3px 8px;
    font-size: 8.5pt;
    font-weight: 650;
}
QLabel#statusChip[state="finished"] {
    color: #baf5db;
    background: #12352c;
    border-color: #286b56;
}
QLabel#statusChip[state="partial"] {
    color: #f5ddb0;
    background: #352b19;
    border-color: #74603a;
}
QLabel#statusChip[state="running"] {
    color: #c9e2ff;
    background: #172f4a;
    border-color: #35628e;
}
QLabel#statusChip[state="failed"] {
    color: #ffd1d1;
    background: #3d2024;
    border-color: #84444c;
}
QLabel#statusChip[state="on_demand"] {
    color: #d7c9ff;
    background: #292340;
    border-color: #5a4a82;
}
QLabel#cardTitle {
    color: #f4fafb;
    font-size: 13.5pt;
    font-weight: 720;
    background: transparent;
}
QLabel#cardBody {
    color: #9eb1bb;
    background: transparent;
}
QLabel#routePill {
    color: #7de0d6;
    background: #0d302f;
    border: 1px solid #23605c;
    border-radius: 8px;
    padding: 3px 8px;
    font-size: 8.5pt;
    font-weight: 700;
}
QFrame#navigationCard {
    background: #101d27;
    border: 1px solid #203745;
    border-radius: 14px;
}
QFrame#navigationCard:hover {
    background: #132631;
    border-color: #3dbeb3;
}
QTreeWidget#workflowNavigation {
    background: transparent;
    border: 0;
    padding: 2px 0;
    outline: 0;
}
QTreeWidget#workflowNavigation::item {
    min-height: 31px;
    padding: 3px 7px;
    margin: 1px 0;
    border-radius: 8px;
    color: #9fb1ba;
}
QTreeWidget#workflowNavigation::item:selected {
    background: #156e6a;
    color: #ffffff;
}
QTreeWidget#workflowNavigation::item:hover:!selected {
    background: #112733;
    color: #eaf5f5;
}
QTreeWidget#workflowNavigation::branch {
    background: transparent;
}
QPushButton {
    background: #172733;
    border: 1px solid #304957;
    border-radius: 8px;
    min-height: 30px;
    padding: 4px 13px;
    font-weight: 560;
}
QPushButton:hover {
    background: #203946;
    border-color: #43c4b8;
}
QPushButton:pressed {
    background: #105e5b;
}
QPushButton[class="primary"] {
    background: #16867f;
    border-color: #31b8ad;
    color: #ffffff;
    font-weight: 650;
}
QPushButton[class="primary"]:hover {
    background: #1a9b92;
}
QPushButton:disabled {
    color: #647680;
    background: #111c24;
    border-color: #23333d;
}
QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox, QDateEdit {
    background: #0a151d;
    border: 1px solid #304957;
    border-radius: 8px;
    min-height: 30px;
    padding: 2px 9px;
    selection-background-color: #177a75;
}
QLineEdit:focus, QComboBox:focus, QSpinBox:focus, QDoubleSpinBox:focus {
    border: 1px solid #49c9bd;
}
QComboBox QAbstractItemView {
    background: #101d27;
    border: 1px solid #304957;
    selection-background-color: #177a75;
}
QTabWidget::pane {
    border: 1px solid #203744;
    border-radius: 9px;
    top: -1px;
}
QTabBar::tab {
    background: #111f29;
    color: #99adb7;
    min-height: 30px;
    padding: 7px 15px;
    margin-right: 3px;
    border: 1px solid #203744;
    border-radius: 8px 8px 0 0;
}
QTabBar::tab:selected {
    background: #176f6b;
    color: #ffffff;
    border-color: #2a9992;
}
QTabBar::tab:hover:!selected {
    background: #172c37;
}
QGroupBox {
    border: 1px solid #29404d;
    border-radius: 10px;
    margin-top: 13px;
    padding: 14px 9px 9px;
    font-weight: 650;
}
QGroupBox::title {
    subcontrol-origin: margin;
    left: 12px;
    padding: 0 6px;
    color: #72d4cb;
}
QTableView, QTableWidget, QTreeView {
    background: #0b161e;
    alternate-background-color: #10202a;
    gridline-color: #213641;
    border: 1px solid #243b47;
    border-radius: 8px;
    selection-background-color: #176f6b;
    selection-color: #ffffff;
}
QHeaderView::section {
    background: #142632;
    color: #dce8ec;
    padding: 7px;
    border: 0;
    border-right: 1px solid #29414e;
    border-bottom: 1px solid #29414e;
    font-weight: 650;
}
QToolBar#mainToolbar {
    background: #0a141c;
    border-bottom: 1px solid #1e333f;
    spacing: 8px;
    padding: 7px 14px;
}
QLabel#toolbarTitle {
    color: #dcebef;
    font-size: 11pt;
    font-weight: 700;
    background: transparent;
}
QStatusBar {
    background: #081119;
    color: #8fa4af;
    border-top: 1px solid #1e333f;
}
QSplitter#workspaceSplitter::handle {
    background: #1d303b;
}
QScrollArea {
    border: 0;
    background: transparent;
}
QScrollBar:vertical {
    background: transparent;
    width: 11px;
    margin: 2px;
}
QScrollBar::handle:vertical {
    background: #395463;
    min-height: 30px;
    border-radius: 5px;
}
QScrollBar:horizontal {
    background: transparent;
    height: 11px;
    margin: 2px;
}
QScrollBar::handle:horizontal {
    background: #395463;
    min-width: 30px;
    border-radius: 5px;
}
QToolTip {
    background: #142733;
    color: #edf6f7;
    border: 1px solid #3a5a68;
    padding: 5px;
}
"""
