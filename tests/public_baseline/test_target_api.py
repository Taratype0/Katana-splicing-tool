"""HTTP lifecycle/ownership tests; numerical engines are not executed here."""
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
from uuid import uuid4
import zipfile

from fastapi import FastAPI
from fastapi.testclient import TestClient
import pytest

from katana_next.api import targets as target_api


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture
def target_http(tmp_path, monkeypatch):
    project = SimpleNamespace(manifest=SimpleNamespace(output_root=tmp_path / "output"))
    def get(project_id):
        if project_id != "study":
            raise KeyError(project_id)
        return project
    calls = {"exports": [], "runs": [], "plans": []}
    def export(context, selections, output_dir, *, modules=None, contrast_ids=None):
        root = Path(output_dir)
        root.mkdir(parents=True)
        request_id = uuid4().hex
        request = {"project_id": "study", "request_id": request_id, "selections": selections}
        path = root / "request.json"
        path.write_text(json.dumps(request), encoding="utf-8")
        archive = root.with_suffix(".zip")
        with zipfile.ZipFile(archive, "w") as packed:
            packed.write(path, "request.json")
        calls["exports"].append({"context": context, "selections": selections, "modules": modules, "contrast_ids": contrast_ids})
        return {"request_id": request_id, "request_path": str(path), "archive_path": str(archive),
                "bundle_directory": str(root), "target_count": len(selections), "targets": [],
                "resources": [{"resource_id": "bam:sample1", "kind": "bam", "sample_id": "sample1"}]}
    def inspect(path):
        path = Path(path)
        root = path.parent if path.is_file() else path
        receipt = json.loads((root / "result-manifest.json").read_text())
        for name, proof in receipt["members"].items():
            source = root / name
            if not source.resolve().is_relative_to(root.resolve()) or source.is_symlink() or digest(source) != proof["sha256"]:
                raise ValueError("Tampered target member")
        return {**receipt, "verified": True, "report": str(root / "INDEX.html"),
                "archive": str(root / "REPORT.zip"), "manifest": str(root / "result-manifest.json")}
    def run(request_path, output_dir, resources=None, *, workers=None, target_ids=None, progress=None):
        calls["runs"].append({"outdir": str(output_dir), "resources": resources, "workers": workers, "target_ids": target_ids})
        root = Path(output_dir)
        root.mkdir(parents=True, exist_ok=True)
        request = json.loads(Path(request_path).read_text())
        if progress:
            progress({"target_id": "event1", "status": "running"})
        (root / "targets/event1").mkdir(parents=True, exist_ok=True)
        (root / "INDEX.html").write_text("<h1>Target index</h1>")
        (root / "targets/event1/REPORT.html").write_text("<h1>Exact event</h1>")
        (root / "targets/event1/data.json").write_text('{"event_id":"exact:event1","FDR":0.17}')
        files = ["INDEX.html", "targets/event1/REPORT.html", "targets/event1/data.json"]
        receipt = {"schema": "katana.target_results/v1", "status": "complete", "project_id": "study",
                   "request_id": request["request_id"], "targets": [{"target_id": "event1", "status": "complete", "report": "targets/event1/REPORT.html"}],
                   "members": {name: {"sha256": digest(root / name), "bytes": (root / name).stat().st_size} for name in files}}
        (root / "result-manifest.json").write_text(json.dumps(receipt))
        with zipfile.ZipFile(root / "REPORT.zip", "w") as packed:
            for name in [*files, "result-manifest.json"]:
                packed.write(root / name, name)
        return inspect(root)
    def plan(path, resources=None, *, outdir=None):
        calls["plans"].append((str(path), resources))
        return {"ready": True, "errors": [], "warnings": [], "tasks": [{"target_id": "event1"}]}
    monkeypatch.setattr(target_api, "export_targets", export)
    monkeypatch.setattr(target_api, "run_targets", run)
    monkeypatch.setattr(target_api, "inspect_results", inspect)
    monkeypatch.setattr(target_api, "plan_targets", plan)
    app = FastAPI()
    app.include_router(target_api.create_targets_router(SimpleNamespace(projects=SimpleNamespace(get=get))))
    with TestClient(app) as client:
        yield client, calls, tmp_path


def selection():
    return {"source_artifact_id": "source-exact", "target_kind": "event", "row_ids": ["row:non-significant"],
            "filters": {"source_parameters": {"fdr": .05}, "browse_settings": {"fdr": .2, "event_set": "stable_significant", "contrast_ids": ["B_vs_A"]}}}


def test_export_run_reopen_and_import_keep_snapshot_and_owned_files(target_http):
    client, calls, tmp = target_http
    base = "/api/projects/study/targets"
    response = client.post(base + "/export", json={"selections": [selection()], "modules": ["event_card"]})
    assert response.status_code == 200
    exported = response.json()
    assert calls["exports"][0]["selections"] == [selection()]
    assert client.get(exported["download_url"]).content.startswith(b"PK")
    resources = {"bam:sample1": "/execution-machine/sample.bam"}
    assert client.post(base + "/plan", json={"request_id": exported["request_id"], "resources": resources}).json()["ready"] is True
    assert len(calls["exports"]) == 1 and not calls["runs"]
    response = client.post(base + "/run", json={"request_id": exported["request_id"], "resources": resources, "workers": 7})
    assert response.status_code == 202
    run = client.get(base + "/runs/" + response.json()["run_id"]).json()
    assert run["status"] == "complete" and "no scheduler" in run["execution"]
    assert calls["runs"][0]["resources"] == resources and calls["runs"][0]["workers"] == 7
    assert "Target index" in client.get(run["report_url"]).text
    data_url = base + f"/runs/{run['run_id']}/files/targets/event1/data.json"
    assert client.get(data_url).json() == {"event_id": "exact:event1", "FDR": .17}
    assert client.get(run["archive_url"]).content.startswith(b"PK")
    # Repeated requests reuse the same owned product directory, allowing the
    # real runner's per-target cache to preserve completed files.
    client.post(base + "/run", json={"request_id": exported["request_id"], "resources": resources})
    assert calls["runs"][0]["outdir"] == calls["runs"][1]["outdir"]
    # Reopen a result copied from another location without serving that location.
    outside = tmp / "returned"
    import shutil
    shutil.copytree(run["output_directory"], outside)
    response = client.post(base + "/results/import", json={"manifest_path": str(outside / "result-manifest.json")})
    assert response.status_code == 200 and response.json()["imported"] is True
    assert Path(response.json()["output_directory"]).is_relative_to(tmp / "output/target-workspace/imports")
    assert "Target index" in client.get(response.json()["report_url"]).text


def test_rejects_unknown_project_unsafe_identity_unregistered_files_and_tampering(target_http):
    client, _, _ = target_http
    base = "/api/projects/study/targets"
    assert client.post("/api/projects/other/targets/export", json={"selections": [selection()]}).status_code == 404
    assert client.post(base + "/export", json={"selections": []}).status_code == 422
    exported = client.post(base + "/export", json={"selections": [selection()]}).json()
    queued = client.post(base + "/run", json={"request_id": exported["request_id"]}).json()
    result = client.get(base + "/runs/" + queued["run_id"]).json()
    root = Path(result["output_directory"])
    (root / "secret.txt").write_text("not a registered member")
    assert client.get(base + f"/runs/{queued['run_id']}/files/secret.txt").status_code == 404
    assert client.get(base + f"/runs/{queued['run_id']}/files/%2e%2e%2Fsecret.txt").status_code == 404
    (root / "targets/event1/data.json").write_text("tampered")
    assert client.get(result["report_url"]).status_code == 404
    assert client.post(base + "/results/import", json={"manifest_path": str(root)}).status_code == 422


def test_export_cannot_register_archive_outside_owned_root(target_http, monkeypatch):
    client, _, tmp = target_http
    outsider = tmp / "unowned.zip"
    outsider.write_bytes(b"secret")
    original = target_api.export_targets
    def escaping(*args, **kwargs):
        value = original(*args, **kwargs)
        value["archive_path"] = str(outsider)
        return value
    monkeypatch.setattr(target_api, "export_targets", escaping)
    response = client.post("/api/projects/study/targets/export", json={"selections": [selection()]})
    assert response.status_code == 422
    assert outsider.read_bytes() == b"secret"


def test_export_download_rejects_changed_archive(target_http):
    client, _, _ = target_http
    exported = client.post("/api/projects/study/targets/export", json={"selections": [selection()]}).json()
    Path(exported["archive_path"]).write_bytes(b"changed after export")
    response = client.get(exported["download_url"])
    assert response.status_code == 404
    assert "modified after verification" in response.json()["detail"]


def test_active_request_is_not_submitted_twice_and_bindings_cannot_silently_change(target_http, monkeypatch):
    client, calls, _ = target_http
    pending = []
    monkeypatch.setattr(target_api.BackgroundTasks, "add_task", lambda _self, *args, **kwargs: pending.append((args, kwargs)))
    base = "/api/projects/study/targets"
    exported = client.post(base + "/export", json={"selections": [selection()]}).json()
    payload = {"request_id": exported["request_id"], "resources": {"bam:sample1": "/first.bam"}}
    first = client.post(base + "/run", json=payload)
    second = client.post(base + "/run", json=payload)
    assert first.status_code == second.status_code == 202
    assert first.json()["run_id"] == second.json()["run_id"] and len(pending) == 1
    changed = client.post(base + "/run", json={**payload, "resources": {"bam:sample1": "/different.bam"}})
    assert changed.status_code == 409 and len(pending) == 1 and not calls["runs"]


def test_import_rejects_foreign_project_before_copy(target_http, monkeypatch):
    client, _, tmp = target_http
    monkeypatch.setattr(target_api, "inspect_results", lambda path: {"verified": True, "project_id": "foreign"})
    response = client.post("/api/projects/study/targets/results/import", json={"manifest_path": str(tmp)})
    assert response.status_code == 422
    assert not (tmp / "output/target-workspace/imports").exists()


def test_real_portable_bundle_runner_http_roundtrip_and_cache(tmp_path):
    """Exercise the actual new bundle/runner through HTTP with saved toy data."""
    # Reuse the synthetic completed-evidence fixture; no original study or BAM
    # worker is involved, and gene/event card rendering is real.
    from test_target_bundle import study
    toy = study.__wrapped__(tmp_path)
    original = {path: (digest(path), path.stat().st_mtime_ns)
                for path in toy.source.rglob("*") if path.is_file()}
    projects = SimpleNamespace(get=lambda project_id: toy.context if project_id == "toy" else (_ for _ in ()).throw(KeyError(project_id)))
    app = FastAPI()
    app.include_router(target_api.create_targets_router(SimpleNamespace(projects=projects)))
    base = "/api/projects/toy/targets"
    with TestClient(app) as client:
        exported_response = client.post(base + "/export", json={"selections": [toy.select_gene, toy.select_event], "modules": ["gene_card", "event_card"]})
        assert exported_response.status_code == 200, exported_response.text
        exported = exported_response.json()
        assert exported["target_count"] == 2
        planned = client.post(base + "/plan", json={"request_id": exported["request_id"]})
        assert planned.status_code == 200 and not planned.json()["errors"]
        response = client.post(base + "/run", json={"request_id": exported["request_id"], "workers": 2})
        assert response.status_code == 202, response.text
        status = client.get(base + "/runs/" + response.json()["run_id"])
        assert status.status_code == 200, status.text
        result = status.json()
        assert result["status"] == "complete", result
        assert len(result["result"]["targets"]) == 2
        assert client.get(result["report_url"]).status_code == 200
        assert client.get(result["archive_url"]).content.startswith(b"PK")
        completed_files = {path: (digest(path), path.stat().st_mtime_ns)
                           for path in Path(result["output_directory"]).rglob("*") if path.is_file()}
        rerun = client.post(base + "/run", json={"request_id": exported["request_id"], "workers": 1})
        assert rerun.status_code == 202
        cached = client.get(base + "/runs/" + rerun.json()["run_id"]).json()
        assert cached["status"] in {"cache_hit", "complete"}
        assert all((digest(path), path.stat().st_mtime_ns) == proof for path, proof in completed_files.items())
        imported = client.post(base + "/results/import", json={"manifest_path": result["result"]["manifest"]})
        assert imported.status_code == 200, imported.text
        assert client.get(imported.json()["report_url"]).status_code == 200
    assert all((digest(path), path.stat().st_mtime_ns) == proof for path, proof in original.items())
