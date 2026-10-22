import io
import json
import subprocess
from pathlib import Path
from unittest.mock import Mock

from katana_next import launcher


def test_health_requires_a_katana_inventory(monkeypatch):
    responses = iter([{"status": "ok"}, [{"workspace_id": name} for name in ("project", "differential", "programs", "reports")]])
    monkeypatch.setattr(launcher.urllib.request, "urlopen", lambda *a, **kw: io.BytesIO(json.dumps(next(responses)).encode()))
    assert launcher.service_ready("http://127.0.0.1:8765/")


def test_foreign_service_not_mistaken_for_katana(monkeypatch):
    responses = iter([{"status": "ok"}, [{"workspace_id": "other"}]])
    monkeypatch.setattr(launcher.urllib.request, "urlopen", lambda *a, **kw: io.BytesIO(json.dumps(next(responses)).encode()))
    assert not launcher.service_ready("http://127.0.0.1:8765/")


def test_existing_service_is_reused_without_process_start(tmp_path, monkeypatch):
    monkeypatch.setattr(launcher, "service_ready", lambda url: True)
    spawn = Mock()
    monkeypatch.setattr(launcher.subprocess, "Popen", spawn)
    assert launcher.ensure_service(8765, tmp_path) == "http://127.0.0.1:8765/"
    spawn.assert_not_called()


def test_alternate_port_is_reused_when_default_is_occupied(tmp_path, monkeypatch):
    saved_url = "http://127.0.0.1:9876/"
    (tmp_path / "service.json").write_text(json.dumps({"requested_port": 8765, "url": saved_url}))
    monkeypatch.setattr(launcher, "service_ready", lambda url: url == saved_url)
    spawn = Mock()
    monkeypatch.setattr(launcher.subprocess, "Popen", spawn)
    assert launcher.ensure_service(8765, tmp_path) == saved_url
    spawn.assert_not_called()


def test_new_server_is_detached_logged_and_waited_for(tmp_path, monkeypatch):
    statuses = iter([False, True])
    monkeypatch.setattr(launcher, "service_ready", lambda url: next(statuses))
    process = Mock(pid=42)
    process.poll.return_value = None
    spawn = Mock(return_value=process)
    monkeypatch.setattr(launcher.subprocess, "Popen", spawn)
    # Bind an ephemeral port to obtain a valid unused port for this unit test.
    with launcher.socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        port = sock.getsockname()[1]
    url = launcher.ensure_service(port, tmp_path)
    assert url == f"http://127.0.0.1:{port}/"
    assert spawn.call_args.args[0][2:4] == ["katana_next.launcher", "--serve"]
    assert json.loads((tmp_path / "service.json").read_text())["pid"] == 42
    assert (tmp_path / f"server-{port}.stderr.log").is_file()


def test_browser_mode_opens_only_after_health_check(monkeypatch, tmp_path):
    url = "http://127.0.0.1:8765/"
    ensure = Mock(return_value=url)
    opener = Mock(return_value=True)
    monkeypatch.setattr(launcher, "ensure_service", ensure)
    monkeypatch.setattr(launcher.webbrowser, "open", opener)
    monkeypatch.setattr(launcher.sys, "argv", ["launcher", "--mode", "browser", "--runtime-root", str(tmp_path)])
    launcher.main()
    ensure.assert_called_once()
    opener.assert_called_once_with(url, new=2)


def test_script_no_open_does_not_create_a_window(monkeypatch, tmp_path):
    monkeypatch.setattr(launcher, "ensure_service", lambda *a: "http://127.0.0.1:8765/")
    opener = Mock()
    monkeypatch.setattr(launcher.webbrowser, "open", opener)
    chooser = Mock()
    monkeypatch.setattr(launcher, "choose_mode", chooser)
    monkeypatch.setattr(launcher.sys, "argv", ["launcher", "--no-open", "--runtime-root", str(tmp_path)])
    launcher.main()
    chooser.assert_not_called()
    opener.assert_not_called()


def test_desktop_connects_to_same_service(monkeypatch, tmp_path):
    import pytest
    from katana_next import desktop
    monkeypatch.setattr(launcher, "ensure_service", lambda *a: "http://127.0.0.1:8765/")
    run = Mock(return_value=0)
    monkeypatch.setattr(desktop, "run_desktop", run)
    monkeypatch.setattr(launcher.sys, "argv", ["launcher", "--mode", "desktop", "--runtime-root", str(tmp_path)])
    with pytest.raises(SystemExit) as exc:
        launcher.main()
    assert exc.value.code == 0
    assert run.call_args.kwargs["server_url"] == "http://127.0.0.1:8765/"


def test_root_launchers_target_current_application():
    root = Path(__file__).resolve().parents[2]
    assert "katana_next.launcher" in (root / "launch_katana.bat").read_text()
    assert "--mode browser" in (root / "Open_Katana_Web.bat").read_text()
    assert "--mode desktop" in (root / "Open_Katana_App.bat").read_text()
    assert "app.main" not in (root / "launch_katana.bat").read_text()


def test_pythonw_logs_streams_before_starting_qt(tmp_path, monkeypatch):
    monkeypatch.setattr(launcher.sys, "stdout", None)
    monkeypatch.setattr(launcher.sys, "stderr", None)
    monkeypatch.setattr(launcher.faulthandler, "enable", Mock())
    launcher.prepare_launch_logging(tmp_path)
    assert launcher.sys.stdout is not None
    assert launcher.sys.stderr is not None
    print("startup failure remains visible", file=launcher.sys.stderr, flush=True)
    assert "startup failure" in next(tmp_path.glob("*.stderr.log")).read_text()
    assert json.loads((tmp_path / "launcher-events.jsonl").read_text())["event"] == "started"


def test_real_pythonw_starts_a_logged_service_from_another_directory(tmp_path):
    import os
    import sys
    import pytest
    if os.name != "nt":
        pytest.skip("Windows shortcut runtime")
    root = Path(__file__).resolve().parents[2]
    with launcher.socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        port = sock.getsockname()[1]
    env = {**os.environ, "PYTHONPATH": str(root)}
    command = [str(Path(sys.executable).with_name("pythonw.exe")), "-m", "katana_next.launcher", "--no-open", "--mode", "desktop", "--port", str(port), "--runtime-root", str(tmp_path / "runtime")]
    process = subprocess.Popen(command, cwd=tmp_path, env=env)
    try:
        assert process.wait(timeout=55) == 0
        service = json.loads((tmp_path / "runtime/service.json").read_text())
        assert launcher.service_ready(service["url"])
        assert list((tmp_path / "runtime").glob("launch-*.stderr.log"))
        records = [json.loads(line) for line in (tmp_path / "runtime/launcher-events.jsonl").read_text().splitlines()]
        assert any(record["event"] == "service_ready" for record in records)
    finally:
        record = tmp_path / "runtime/service.json"
        if record.exists():
            pid = json.loads(record.read_text())["pid"]
            # Only this test's fresh ephemeral service tree is stopped.
            subprocess.run(["taskkill", "/PID", str(pid), "/T", "/F"], capture_output=True, check=False)
        if process.poll() is None:
            process.terminate()
