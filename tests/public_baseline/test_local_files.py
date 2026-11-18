from fastapi import FastAPI
from fastapi.testclient import TestClient

from katana_next.api.local_files import create_local_files_router


def client():
    app = FastAPI()
    app.include_router(create_local_files_router())
    return TestClient(app, base_url="http://127.0.0.1")


def test_browse_and_parent_do_not_modify_files(tmp_path):
    folder = tmp_path / "a report with spaces"
    folder.mkdir()
    manifest = tmp_path / "report-manifest.json"
    manifest.write_text("unchanged")
    with client() as api:
        data = api.post("/api/local-files/browse", json={"path": str(tmp_path)}).json()
        assert [e["name"] for e in data["entries"]] == [folder.name, manifest.name]
        assert data["entries"][0]["directory"] is True
        assert data["entries"][1]["directory"] is False
        nested = api.post("/api/local-files/browse", json={"path": str(folder)}).json()
        assert nested["parent"] == str(tmp_path)
        assert api.post("/api/local-files/browse", json={"path": str(manifest)}).json()["path"] == str(tmp_path)
        assert api.post("/api/local-files/browse", json={"path": str(tmp_path / "missing")}).status_code == 422
    assert manifest.read_text() == "unchanged"
    assert len(list(tmp_path.iterdir())) == 2


def test_remote_origins_and_hosts_cannot_browse(tmp_path):
    with client() as api:
        for origin in ["https://example.com", "null"]:
            assert api.post("/api/local-files/browse", json={"path": str(tmp_path)}, headers={"Origin": origin}).status_code == 403
        assert api.post("/api/local-files/browse", json={"path": str(tmp_path)}, headers={"Host": "example.com"}).status_code == 403


def test_report_import_automatically_allocates_separate_workspace(tmp_path, monkeypatch):
    from test_report_workspace import report_fixture
    from katana_next.api.app import create_app
    monkeypatch.setenv("LOCALAPPDATA", str(tmp_path / "user-data"))
    source, _ = report_fixture(tmp_path)
    with TestClient(create_app()) as api:
        result = api.post("/api/projects/import-report", json={"report_path": str(source)})
        assert result.status_code == 200, result.text
        first = result.json()["snapshot"]["manifest_path"]
        second = api.post("/api/projects/import-report", json={"report_path": str(source)}).json()["snapshot"]["manifest_path"]
    assert first == second  # Reopening keeps the target basket/results in the same workspace.
    assert str(tmp_path / "user-data" / "Katana" / "workspaces") in first
