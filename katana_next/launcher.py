"""Open Katana in a browser backed by a local Python analysis service."""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import faulthandler
import json
import os
from pathlib import Path
import socket
import subprocess
import sys
import time
import traceback
import urllib.request
import webbrowser

_chooser_application = None
_windowed_launch = False
_launch_runtime: Path | None = None


def default_runtime_root() -> Path:
    """Keep per-user launcher state outside the installed Python package."""
    home = Path.home()
    if sys.platform == "win32":
        configured = Path(os.environ.get("LOCALAPPDATA", ""))
        base = configured if configured.is_absolute() else home / "AppData" / "Local"
        return base / "Katana" / "runtime"
    if sys.platform == "darwin":
        return home / "Library" / "Application Support" / "Katana" / "runtime"
    configured = Path(os.environ.get("XDG_DATA_HOME", ""))
    base = configured if configured.is_absolute() else home / ".local" / "share"
    return base / "katana" / "runtime"


def prepare_launch_logging(root: Path) -> None:
    """Keep failures observable when a shortcut uses pythonw without a console."""
    global _windowed_launch, _launch_runtime
    root.mkdir(parents=True, exist_ok=True)
    _launch_runtime = root
    _windowed_launch = sys.stderr is None
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S")
    for stream in ("stdout", "stderr"):
        if getattr(sys, stream) is None:
            path = root / f"launch-{stamp}-{os.getpid()}.{stream}.log"
            setattr(sys, stream, path.open("a", encoding="utf-8", buffering=1))
    faulthandler.enable(file=sys.stderr)
    record_launch("started", arguments=sys.argv[1:])


def record_launch(event: str, **details) -> None:
    if _launch_runtime is not None:
        with (_launch_runtime / "launcher-events.jsonl").open("a", encoding="utf-8") as handle:
            handle.write(json.dumps({"at": datetime.now(timezone.utc).isoformat(), "pid": os.getpid(), "event": event, **details}) + "\n")


def service_ready(url: str) -> bool:
    try:
        from katana_next import __version__
        with urllib.request.urlopen(url + "api/health", timeout=1) as response:
            health = json.load(response)
        # Require a Katana inventory, not just another service's generic health endpoint.
        with urllib.request.urlopen(url + "api/workspaces", timeout=2) as response:
            inventory = json.load(response)
        ids = {item["workspace_id"] for item in inventory}
        return health.get("status") == "ok" and health.get("version") == __version__ and {"project", "differential", "programs", "reports"} <= ids
    except (OSError, ValueError, TypeError, KeyError):
        return False


@contextmanager
def startup_lock(root: Path):
    path = root / "startup.lock"
    with path.open("a+b") as handle:
        handle.seek(0, 2)
        if handle.tell() == 0:
            handle.write(b"0")
            handle.flush()
        deadline = time.monotonic() + 45
        while True:
            try:
                handle.seek(0)
                if os.name == "nt":
                    import msvcrt
                    msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
                else:
                    import fcntl
                    fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except OSError:
                if time.monotonic() >= deadline:
                    raise TimeoutError("Another Katana launcher is still starting. Try again shortly.")
                time.sleep(.1)
        try:
            yield
        finally:
            handle.seek(0)
            if os.name == "nt":
                import msvcrt
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                import fcntl
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def ensure_service(port: int, runtime_root: Path) -> str:
    repo = Path(__file__).resolve().parents[1]
    from katana_next.web_assets import frontend_directory
    if frontend_directory() is None:
        raise FileNotFoundError("Katana interface is not built. Install a browser wheel or build apps/web (pnpm build).")
    runtime_root.mkdir(parents=True, exist_ok=True)
    with startup_lock(runtime_root):
        requested_port = port
        url = f"http://127.0.0.1:{port}/"
        if service_ready(url):
            return url
        record = runtime_root / "service.json"
        if record.is_file():
            try:
                saved = json.loads(record.read_text(encoding="utf-8"))
                saved_url = saved.get("url", "")
                if saved.get("requested_port") == port and saved_url.startswith("http://127.0.0.1:") and service_ready(saved_url):
                    return saved_url
            except (OSError, ValueError, TypeError):
                pass
        with socket.socket() as probe:
            try:
                probe.bind(("127.0.0.1", port))
            except OSError:
                probe.bind(("127.0.0.1", 0))
                port = probe.getsockname()[1]
                url = f"http://127.0.0.1:{port}/"
        stdout = runtime_root / f"server-{port}.stdout.log"
        stderr = runtime_root / f"server-{port}.stderr.log"
        executable = Path(sys.executable)
        console_python = executable.with_name("python.exe")
        if executable.name.lower() == "pythonw.exe" and console_python.is_file():
            executable = console_python
        command = [str(executable), "-m", "katana_next.launcher", "--serve", "--port", str(port), "--runtime-root", str(runtime_root)]
        options = {"creationflags": subprocess.CREATE_NO_WINDOW} if os.name == "nt" else {"start_new_session": True}
        with stdout.open("ab") as out, stderr.open("ab") as err:
            process = subprocess.Popen(command, cwd=repo, stdin=subprocess.DEVNULL, stdout=out, stderr=err, **options)
        deadline = time.monotonic() + 40
        while time.monotonic() < deadline:
            if service_ready(url):
                record.write_text(json.dumps({"url": url, "requested_port": requested_port, "pid": process.pid, "stdout": str(stdout), "stderr": str(stderr)}, indent=2), encoding="utf-8")
                return url
            if process.poll() is not None:
                raise RuntimeError(f"Katana startup failed. See {stderr}")
            time.sleep(.2)
        raise TimeoutError(f"Katana is not ready. See {stderr}; no second server was started.")


def choose_mode() -> str | None:
    global _chooser_application
    from PySide6.QtWidgets import QApplication, QDialog, QLabel, QPushButton, QVBoxLayout
    from PySide6.QtGui import QIcon
    from PySide6.QtWebEngineWidgets import QWebEngineView  # Initialize before the chooser's QApplication.
    from katana_next.desktop import desktop_icon_path
    app = QApplication.instance() or QApplication([])
    _chooser_application = app
    # Closing the chooser must not enqueue application exit before the main window.
    app.setQuitOnLastWindowClosed(False)
    app.setWindowIcon(QIcon(str(desktop_icon_path())))
    dialog = QDialog()
    dialog.setWindowTitle("Open Katana")
    dialog.setMinimumWidth(460)
    layout = QVBoxLayout(dialog)
    layout.setContentsMargins(28, 24, 28, 24)
    layout.setSpacing(16)
    heading = QLabel("How would you like to open Katana?")
    heading.setStyleSheet("font-size: 19px; font-weight: 600")
    layout.addWidget(heading)
    layout.addWidget(QLabel("Both options provide the same analyses, settings and results."))
    selected = []
    for mode, text in (("browser", "Web browser\nOpen Katana in your default browser"), ("desktop", "Desktop application\nOpen Katana in a standalone window")):
        button = QPushButton(text)
        button.setMinimumHeight(64)
        button.clicked.connect(lambda _checked=False, choice=mode: (selected.append(choice), dialog.accept()))
        layout.addWidget(button)
    dialog.exec()
    return selected[0] if selected else None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["choose", "browser", "desktop"], default="browser",
                        help="Open in a browser by default; choose and desktop are compatibility modes.")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--runtime-root", type=Path, default=default_runtime_root(),
                        help="Writable user state directory; defaults to the platform's user data directory.")
    parser.add_argument("--no-open", action="store_true", help="Start or check the server and print the URL without opening a window.")
    parser.add_argument("--serve", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if not 1 <= args.port <= 65535:
        parser.error("port must be between 1 and 65535")
    args.runtime_root = args.runtime_root.expanduser().resolve()
    prepare_launch_logging(args.runtime_root)
    from katana_next.desktop import desktop_runtime_root
    desktop_runtime_root(args.runtime_root / "desktop")
    if args.serve:
        desktop_runtime_root(args.runtime_root)
        import uvicorn
        uvicorn.run("katana_next.api.app:create_app", factory=True, host="127.0.0.1", port=args.port)
        return
    mode = args.mode if args.mode != "choose" or args.no_open else choose_mode()
    if mode is None:
        record_launch("chooser_cancelled")
        return
    url = ensure_service(args.port, args.runtime_root.resolve())
    if sys.stdout:
        print(f"Katana: {url}", flush=True)
    if args.no_open:
        record_launch("service_ready", url=url, mode=mode)
        return
    if mode == "browser":
        if not webbrowser.open(url, new=2):
            raise RuntimeError(f"Your browser could not be opened automatically. Open {url}")
        record_launch("browser_opened", url=url)
    else:
        from katana_next.desktop import run_desktop
        record_launch("desktop_opening", url=url)
        raise SystemExit(run_desktop(runtime_root=args.runtime_root / "desktop", server_url=url))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        record_launch("failed", error=str(exc))
        if sys.stderr:
            traceback.print_exc(file=sys.stderr)
        if _windowed_launch:
            from PySide6.QtWidgets import QApplication, QMessageBox
            app = QApplication.instance() or QApplication([])
            QMessageBox.critical(None, "Katana could not open", f"{exc}\n\nStartup logs: {_launch_runtime}")
        raise SystemExit(1) from exc
