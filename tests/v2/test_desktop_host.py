from __future__ import annotations

from katana_next.desktop import LocalApiServer, _available_port


def test_desktop_api_host_starts_and_stops() -> None:
    server = LocalApiServer(_available_port(0))
    server.start()
    try:
        server.wait_until_ready(timeout=15)
        assert server.is_alive()
    finally:
        server.stop()
    assert not server.is_alive()
