from pathlib import Path


def test_start_page_does_not_disable_document_scrolling():
    styles = (Path(__file__).resolve().parents[2] / "apps/web/src/styles.css").read_text()
    assert "body { overflow: hidden; }" not in styles
    assert "body { overflow: auto; }" in styles


def test_workspace_has_one_scroll_owner_including_configuration_warnings():
    root = Path(__file__).resolve().parents[2]
    styles = (root / "apps/web/src/styles.css").read_text()
    main = next(line for line in styles.splitlines() if line.startswith(".app-main {"))
    workbench = next(line for line in styles.splitlines() if line.startswith(".workbench {"))
    assert "overflow: auto" in main
    assert "overflow: visible" in workbench
    assert 'tabIndex={0} aria-label="Analysis workspace"' in (root / "apps/web/src/App.tsx").read_text()


def test_desktop_entry_preserves_chooser_lifetime_and_small_screens():
    root = Path(__file__).resolve().parents[2]
    launcher = (root / "katana_next/launcher.py").read_text()
    desktop = (root / "katana_next/desktop.py").read_text()
    assert "setQuitOnLastWindowClosed(False)" in launcher
    assert "setQuitOnLastWindowClosed(True)" in desktop
    assert "availableGeometry()" in desktop
    assert "heartbeat.start()" in desktop
