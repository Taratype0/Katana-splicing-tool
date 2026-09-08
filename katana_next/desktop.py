from __future__ import annotations

import argparse
import socket
import threading
import time
import urllib.request
from pathlib import Path

import uvicorn

from katana_next.api.app import create_app


def _available_port(preferred: int = 8765) -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe:
        try:
            probe.bind(("127.0.0.1", preferred))
        except OSError:
            probe.bind(("127.0.0.1", 0))
        return int(probe.getsockname()[1])


class LocalApiServer(threading.Thread):
    def __init__(self, port: int) -> None:
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


def run_desktop(*, port: int = 8765, check_only: bool = False) -> int:
    frontend = Path(__file__).resolve().parents[1] / "apps" / "web" / "dist" / "index.html"
    if not frontend.is_file():
        raise FileNotFoundError(
            "Katana web assets are missing. Run `pnpm --dir apps/web build` first."
        )
    server = LocalApiServer(_available_port(port))
    server.start()
    try:
        server.wait_until_ready()
        if check_only:
            return 0

        from PySide6.QtCore import QUrl
        from PySide6.QtGui import QCloseEvent
        from PySide6.QtWidgets import QApplication
        from PySide6.QtWebEngineWidgets import QWebEngineView

        class KatanaWindow(QWebEngineView):
            def closeEvent(self, event: QCloseEvent) -> None:  # noqa: N802
                server.stop()
                super().closeEvent(event)

        application = QApplication.instance() or QApplication([])
        application.setApplicationDisplayName("Katana Next")
        application.setOrganizationName("Katana")
        window = KatanaWindow()
        window.setWindowTitle("Katana Next")
        window.setMinimumSize(1040, 680)
        window.resize(1480, 940)
        window.setUrl(QUrl(server.url))
        window.show()
        return int(application.exec())
    finally:
        if server.is_alive():
            server.stop()


def main() -> None:
    parser = argparse.ArgumentParser(description="Launch the Katana Next desktop app.")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument(
        "--check",
        action="store_true",
        help="Start the embedded API, verify the packaged web UI, then exit.",
    )
    arguments = parser.parse_args()
    raise SystemExit(run_desktop(port=arguments.port, check_only=arguments.check))


if __name__ == "__main__":
    main()
