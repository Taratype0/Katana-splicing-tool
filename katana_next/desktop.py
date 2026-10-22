from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
import socket
import threading
import time
import sys
import urllib.request
from pathlib import Path

import uvicorn

def desktop_icon_path() -> Path:
    path = Path(__file__).resolve().parents[1] / "assets" / "katana_icon.ico"
    if not path.is_file():
        raise FileNotFoundError(f"Katana application icon is missing: {path}")
    return path

def _available_port(preferred: int = 8765) -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe:
        try:
            probe.bind(("127.0.0.1", preferred))
        except OSError:
            probe.bind(("127.0.0.1", 0))
        return int(probe.getsockname()[1])


class LocalApiServer(threading.Thread):
    def __init__(self, port: int) -> None:
        from katana_next.api.app import create_app
        super().__init__(name="katana-next-api", daemon=True)
        self.port = port
        self.server = uvicorn.Server(
            uvicorn.Config(
                create_app(),
                host="127.0.0.1",
                port=port,
                log_level="warning",
                access_log=False,
            )
        )

    @property
    def url(self) -> str:
        return f"http://127.0.0.1:{self.port}/"

    def run(self) -> None:
        self.server.run()

    def wait_until_ready(self, timeout: float = 15.0) -> None:
        deadline = time.monotonic() + timeout
        health = f"{self.url}api/health"
        last_error: Exception | None = None
        while time.monotonic() < deadline:
            if not self.is_alive():
                raise RuntimeError("Katana local API stopped during startup")
            try:
                with urllib.request.urlopen(health, timeout=1.0) as response:
                    if response.status == 200:
                        return
            except Exception as exc:  # pragma: no cover - platform timing
                last_error = exc
            time.sleep(0.1)
        raise TimeoutError(f"Katana local API was not ready: {last_error}")

    def stop(self) -> None:
        self.server.should_exit = True
        self.join(timeout=8.0)


def desktop_runtime_root(supplied: str | Path | None = None) -> Path:
    base = Path(sys.executable).resolve().parent.parent if getattr(sys, "frozen", False) else Path(__file__).resolve().parents[2] / "katana-next-runs"
    root = Path(supplied or os.environ.get("KATANA_DESKTOP_RUNTIME", "") or base / "desktop-runtime").resolve()
    root.mkdir(parents=True, exist_ok=True)
    for child in ("cache", "storage", "downloads", "temp", "matplotlib", "numba"):
        (root / child).mkdir(exist_ok=True)
    os.environ["TEMP"] = os.environ["TMP"] = str(root / "temp")
    os.environ["MPLCONFIGDIR"] = str(root / "matplotlib")
    os.environ["NUMBA_CACHE_DIR"] = str(root / "numba")
    return root


def run_desktop(*, port: int = 8765, check_only: bool = False, runtime_root: str | Path | None = None, smoke_fixture: Path | None = None, server_url: str | None = None) -> int:
    runtime = desktop_runtime_root(runtime_root)
    icon_path = desktop_icon_path()
    frontend = Path(__file__).resolve().parents[1] / "apps" / "web" / "dist" / "index.html"
    if not frontend.is_file():
        raise FileNotFoundError(
            "Katana web assets are missing. Run `pnpm --dir apps/web build` first."
        )
    server = None if server_url else LocalApiServer(_available_port(port))
    if server:
        server.start()
    try:
        if server:
            server.wait_until_ready()
        else:
            from katana_next.launcher import service_ready
            if not server_url.startswith("http://127.0.0.1:") or not service_ready(server_url):
                raise RuntimeError("The shared Katana service is not ready on loopback.")
        url = server_url or server.url
        if check_only:
            from katana_next.workflows.registry import WorkflowRegistry
            evidence = {"embedded_api_ready": True, "frontend_exists": True, "application_icon": str(icon_path),
                        "workflow_ids": [w.definition.workflow_id for w in WorkflowRegistry().all()],
                        "frozen": bool(getattr(sys, "frozen", False)), "runtime_root": str(runtime),
                        "acceptance_scope": "API and packaged static assets only; not physical monitor/DPI or interactive Qt acceptance"}
            if smoke_fixture:
                from katana_next.desktop_smoke import run_packaged_smoke
                evidence["packaged_workflow_smoke"] = run_packaged_smoke(url, smoke_fixture)
            (runtime / "startup-check.json").write_text(json.dumps(evidence, indent=2), encoding="utf-8")
            return 0

        from PySide6.QtCore import QTimer, QUrl, qInstallMessageHandler
        from PySide6.QtGui import QCloseEvent, QIcon
        from PySide6.QtWidgets import QApplication
        from PySide6.QtWebEngineWidgets import QWebEngineView
        from PySide6.QtWebEngineCore import QWebEnginePage, QWebEngineProfile

        state = {"status": "starting", "loaded": False, "url": url, "pid": os.getpid()}

        def record_state(status: str | None = None) -> None:
            if status:
                state["status"] = status
            state.update(at=datetime.now(timezone.utc).isoformat(), visible=window.isVisible())
            # Each run keeps its own history; the latest heartbeat is not a stale load check.
            path = runtime / f"window-{os.getpid()}.json"
            temporary = path.with_suffix(".tmp")
            temporary.write_text(json.dumps(state, indent=2), encoding="utf-8")
            temporary.replace(path)
            (runtime / "window-load.json").write_text(json.dumps(state, indent=2), encoding="utf-8")

        def qt_message(kind, context, message) -> None:
            with (runtime / f"qt-{os.getpid()}.log").open("a", encoding="utf-8") as handle:
                handle.write(f"{datetime.now(timezone.utc).isoformat()} {kind}: {message}\n")

        previous_handler = qInstallMessageHandler(qt_message)

        class KatanaWindow(QWebEngineView):
            def closeEvent(self, event: QCloseEvent) -> None:  # noqa: N802
                record_state("closed")
                if server:
                    server.stop()
                super().closeEvent(event)

        if sys.platform == "win32":
            import ctypes
            ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("Katana.Next")
        application = QApplication.instance() or QApplication([])
        application.setQuitOnLastWindowClosed(True)
        application.setApplicationDisplayName("Katana")
        application.setOrganizationName("Katana")
        icon = QIcon(str(icon_path))
        if icon.isNull():
            raise ValueError(f"Katana application icon cannot be decoded: {icon_path}")
        application.setWindowIcon(icon)
        window = KatanaWindow()
        window.setWindowIcon(icon)
        profile = QWebEngineProfile("Katana", application)
        profile.setCachePath(str(runtime / "cache"))
        profile.setPersistentStoragePath(str(runtime / "storage"))
        profile.setDownloadPath(str(runtime / "downloads"))
        window.setPage(QWebEnginePage(profile, window))
        window.setWindowTitle("Katana")
        geometry = application.primaryScreen().availableGeometry()
        window.setMinimumSize(min(800, geometry.width()), min(480, geometry.height()))
        window.resize(min(1480, geometry.width()), min(940, geometry.height()))
        def record_loaded(ok: bool) -> None:
            state["loaded"] = ok
            record_state("ready" if ok else "load_failed")
        def record_renderer_exit(reason, code) -> None:
            state.update(renderer_exit_reason=str(reason), renderer_exit_code=code, loaded=False)
            record_state("renderer_failed")
        window.loadFinished.connect(record_loaded)
        window.page().renderProcessTerminated.connect(record_renderer_exit)
        window.setUrl(QUrl(url))
        window.show()
        window.raise_()
        window.activateWindow()
        record_state()
        heartbeat = QTimer(window)
        heartbeat.setInterval(5000)
        heartbeat.timeout.connect(record_state)
        heartbeat.start()
        result = int(application.exec())
        heartbeat.stop()
        record_state("closed")
        # Destroy pages before their profile, including during interpreter shutdown.
        from shiboken6 import delete
        delete(window)
        delete(profile)
        qInstallMessageHandler(previous_handler)
        return result
    finally:
        if server and server.is_alive():
            server.stop()


def main() -> None:
    parser = argparse.ArgumentParser(description="Launch the Katana Next desktop app.")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--server-url", help="Connect to an already running local Katana service.")
    parser.add_argument("--runtime-root", type=Path, help="Explicit desktop cache/storage/temp directory, separate from analysis output.")
    parser.add_argument("--smoke-fixture", type=Path, help="With --check, verify real synthetic BAM/cache/PNG workflows in a new external directory.")
    parser.add_argument(
        "--check",
        action="store_true",
        help="Start the embedded API, verify the packaged web UI, then exit.",
    )
    arguments = parser.parse_args()
    if arguments.smoke_fixture and not arguments.check:
        parser.error("--smoke-fixture requires --check")
    raise SystemExit(run_desktop(port=arguments.port, check_only=arguments.check, runtime_root=arguments.runtime_root, smoke_fixture=arguments.smoke_fixture, server_url=arguments.server_url))


if __name__ == "__main__":
    main()
