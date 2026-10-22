from __future__ import annotations

from pathlib import Path

from fastapi.testclient import TestClient

from katana_next.api.app import PlatformRuntime, create_app
from katana_next.domain.models import JobStatus
from katana_next.fixtures import create_phase1_fixture


def _open(client: TestClient, manifest: Path) -> dict:
    response = client.post("/api/projects/open", json={"manifest_path": str(manifest)})
    assert response.status_code == 200, response.text
    return response.json()


def test_phase1_load_run_filter_select_reopen_and_cache(tmp_path: Path) -> None:
    manifest = create_phase1_fixture(tmp_path / "fixture")
    runtime = PlatformRuntime(application_root=tmp_path / "application")
    client = TestClient(create_app(runtime))

    opened = _open(client, manifest)
    project_id = opened["snapshot"]["project_id"]
    assert opened["artifact_count"] == 0
    assert opened["selection_count"] == 0
    assert client.get("/api/jobs", params={"project_id": project_id}).json() == []

    started = client.post(
        f"/api/projects/{project_id}/workflows/phase1-profile/run",
        json={"parameters": {}, "force": False},
    )
    assert started.status_code == 200, started.text
    job_id = started.json()["job_id"]
    finished = runtime.jobs.wait(job_id)
    assert finished.status == JobStatus.COMPLETE
    assert finished.cache_status == "rebuild"
    assert len(finished.artifact_ids) == 2
    event_stream = client.get(f"/api/jobs/{job_id}/events")
    assert event_stream.status_code == 200
    assert "event: progress" in event_stream.text
    assert "event: complete" in event_stream.text

    artifacts = client.get(f"/api/projects/{project_id}/artifacts").json()
    assert {artifact["kind"] for artifact in artifacts} == {"table", "plot_spec"}
    table = next(artifact for artifact in artifacts if artifact["kind"] == "table")
    plot = next(artifact for artifact in artifacts if artifact["kind"] == "plot_spec")

    rows_response = client.get(
        f"/api/projects/{project_id}/artifacts/{table['artifact_id']}/rows",
        params={
            "search": "Afg2a",
            "limit": 20,
            "sort_by": "abs_effect",
            "descending": True,
        },
    )
    assert rows_response.status_code == 200, rows_response.text
    page = rows_response.json()
    assert page["total"] == 1
    assert page["rows"][0]["gene_symbol"] == "Afg2a"

    plot_response = client.get(
        f"/api/projects/{project_id}/artifacts/{plot['artifact_id']}/content"
    )
    assert plot_response.status_code == 200
    assert plot_response.json()["mark"]["type"] == "bar"

    selection_response = client.post(
        f"/api/projects/{project_id}/selections",
        json={
            "name": "Review events",
            "entity_type": "splice_event",
            "source_artifact_id": table["artifact_id"],
            "entity_ids": [page["rows"][0]["row_id"]],
            "filters": {"search": "Afg2a"},
        },
    )
    assert selection_response.status_code == 200, selection_response.text
    assert selection_response.json()["entity_ids"] == [page["rows"][0]["row_id"]]

    # Restart retains history, but never replays completed work.
    runtime.jobs.shutdown()
    reopened_runtime = PlatformRuntime(application_root=tmp_path / "application")
    reopened_client = TestClient(create_app(reopened_runtime))
    reopened = _open(reopened_client, manifest)
    assert reopened["artifact_count"] == 2
    assert reopened["selection_count"] == 1
    history = reopened_client.get("/api/jobs", params={"project_id": project_id}).json()
    assert [job["job_id"] for job in history] == [job_id]
    assert history[0]["status"] == "complete"
    assert history[0]["completed_at"] == finished.model_dump(mode="json")["completed_at"]

    cached = reopened_client.post(
        f"/api/projects/{project_id}/workflows/phase1-profile/run",
        json={"parameters": {}, "force": False},
    )
    assert cached.status_code == 200, cached.text
    assert cached.json()["status"] == "cache_hit"
    assert cached.json()["cache_status"] == "hit"

    forced = reopened_client.post(
        f"/api/projects/{project_id}/workflows/phase1-profile/run",
        json={"parameters": {}, "force": True},
    )
    forced_job = reopened_runtime.jobs.wait(forced.json()["job_id"])
    assert forced_job.status == JobStatus.COMPLETE
    assert len(reopened_client.get(f"/api/projects/{project_id}/artifacts").json()) == 2
    reopened_runtime.jobs.shutdown()


def test_output_root_inside_source_tree_is_rejected(tmp_path: Path) -> None:
    fixture_root = tmp_path / "application" / "fixture"
    manifest = create_phase1_fixture(fixture_root)
    runtime = PlatformRuntime(application_root=tmp_path / "application")
    client = TestClient(create_app(runtime))

    response = client.post("/api/projects/open", json={"manifest_path": str(manifest)})

    assert response.status_code == 422
    assert "source tree" in response.json()["detail"]
