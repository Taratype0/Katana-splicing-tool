"""MotifTrace navigation and planning share the public catalog and workflow service."""
import hashlib
import json

from fastapi.testclient import TestClient
import pandas as pd

from katana_next.api.app import PlatformRuntime, create_app
from katana_next.full_fixture import create_full_fixture
from katana_next.workflows.base import ArtifactWriter


def test_motiftrace_catalog_planning_and_saved_browser_are_read_only(tmp_path):
    from katana_motiftrace import mode_catalog

    manifest = create_full_fixture(tmp_path / "project")
    inputs = manifest.parent / "inputs"
    before = {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs.rglob("*") if path.is_file()}
    runtime = PlatformRuntime(application_root=tmp_path / "application")
    try:
        with TestClient(create_app(runtime)) as client:
            response = client.get("/api/motiftrace/modes")
            assert response.status_code == 200, response.text
            assert response.json() == mode_catalog()
            catalog = response.json()
            assert catalog["scientific_release"] == "1.0.0-frozen.20261017"
            assert {mode["mode_id"] for mode in catalog["modes"]} == {
                "boundary", "se_impact", "completed_maps", "consequences", "usage", "parallel",
            }
            opened = client.post("/api/projects/open", json={"manifest_path": str(manifest)})
            assert opened.status_code == 200, opened.text
            project_id = opened.json()["snapshot"]["project_id"]
            contrast = opened.json()["manifest"]["contrasts"][0]["contrast_id"]
            for mode in catalog["modes"]:
                plan = client.post(
                    f"/api/projects/{project_id}/workflows/{mode['workflow_id']}/plan",
                    json={"parameters": {**mode["fixed_parameters"], "contrast_ids": [contrast]}, "force": False},
                )
                assert plan.status_code == 200, plan.text
                assert plan.json()["workflow_id"] == mode["workflow_id"]
                assert plan.json()["parameters"]["contrast_ids"] == [contrast]
                if mode["operation"] == "completed_result_import":
                    assert plan.json()["ready"] is False
                    assert plan.json()["errors"]  # Toy project has no completed motif bundle.
            assert runtime.jobs.list(project_id) == []

            mode = next(item for item in catalog["modes"] if item["mode_id"] == "completed_maps")
            context = runtime.projects.get(project_id)
            workflow = runtime.workflows.registry.get(mode["workflow_id"])
            writer = ArtifactWriter(context, mode["workflow_id"], "toy-completed-browser", mode["fixed_parameters"], workflow.version)
            path = writer.directory / "saved_sequence_browser.html"
            payload = b"<!doctype html><html><body><p>Completed motif sequence results.</p></body></html>"
            path.write_bytes(payload)
            browser = writer.existing_file("Completed sequence browser", path, entity_type="sequence_browser", metadata={"mime_type": "text/html"})
            writer.finish()
            records = client.get(f"/api/projects/{project_id}/artifacts", params={"analysis_id": mode["workflow_id"]})
            assert any(item["artifact_id"] == browser.artifact_id for item in records.json())
            result = client.get(f"/api/projects/{project_id}/artifacts/{browser.artifact_id}/file")
            assert result.status_code == 200
            assert result.content == payload
            assert runtime.jobs.list(project_id) == []
    finally:
        runtime.jobs.shutdown()
    assert before == {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in before}


def test_rest_rejects_a_saved_comparison_a_selection_relabelled_as_comparison_b_before_submission(tmp_path):
    from katana_motiftrace import mode_catalog
    from katana_next.science.bulk import READ_SUPPORT_POLICY_VERSION

    manifest = create_full_fixture(tmp_path / "project")
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    # These sentinels satisfy declaration checks only. The comparison guard must
    # reject before a scanner or completed-bundle reader could consume them.
    for kind in ("rbp_motif_database", "fimo_runtime", "rbp_sequence_results"):
        sentinel = manifest.parent / "inputs" / f"guard-only-{kind}.txt"
        sentinel.write_text("Not a scientific input: comparison guard must reject first.\n", encoding="utf-8")
        payload["assets"].append({"asset_id": "guard-" + kind, "kind": kind, "path": str(sentinel), "required": True, "metadata": {}})
    manifest.write_text(json.dumps(payload), encoding="utf-8")
    before = {path: hashlib.sha256(path.read_bytes()).hexdigest()
              for path in (manifest.parent / "inputs").rglob("*") if path.is_file()}
    runtime = PlatformRuntime(application_root=tmp_path / "application")
    try:
        with TestClient(create_app(runtime)) as client:
            opened = client.post("/api/projects/open", json={"manifest_path": str(manifest)})
            assert opened.status_code == 200, opened.text
            project_id = opened.json()["snapshot"]["project_id"]
            comparison_a, comparison_b = [c["contrast_id"] for c in opened.json()["manifest"]["contrasts"][:2]]
            context = runtime.projects.get(project_id)
            writer = ArtifactWriter(context, "differential.as", "rest-real-source-comparisons", {},
                                    runtime.workflows.registry.get("differential.as").version,
                                    common_metadata={"read_support_policy": READ_SUPPORT_POLICY_VERSION})
            frame = pd.DataFrame([
                {"row_id": "row-A", "event_id": "SE:A", "event_type": "SE", "chr": "chr1", "strand": "+", "gene_id": "gene-A", "comparison_id": comparison_a},
                {"row_id": "row-B", "event_id": "SE:B", "event_type": "SE", "chr": "chr1", "strand": "+", "gene_id": "gene-B", "comparison_id": comparison_b},
            ])
            source = writer.table("Existing measured event catalog", "events.parquet", frame, entity_type="splice_event")
            writer.finish()
            original_source = source.path.read_bytes()
            selection = client.post(f"/api/projects/{project_id}/selections", json={
                "name": "Comparison A events", "entity_type": "splice_event",
                "source_artifact_id": source.artifact_id, "entity_ids": ["row-A"], "filters": {},
            })
            assert selection.status_code == 200, selection.text
            selection_id = selection.json()["selection_id"]
            modes = [mode for mode in mode_catalog()["modes"]
                     if mode["mode_id"] in {"boundary", "se_impact", "completed_maps", "consequences"}]
            for mode in modes:
                parameters = {**mode["fixed_parameters"], "selection_id": selection_id, "contrast_ids": [comparison_b]}
                planned = client.post(f"/api/projects/{project_id}/workflows/{mode['workflow_id']}/plan",
                                      json={"parameters": parameters, "force": False})
                assert planned.status_code == 200, planned.text
                assert not planned.json()["ready"]
                assert any("different comparisons" in error for error in planned.json()["errors"])
                attempted = client.post(f"/api/projects/{project_id}/workflows/{mode['workflow_id']}/run",
                                        json={"parameters": parameters, "force": False})
                assert attempted.status_code == 422, attempted.text
                assert "different comparisons" in attempted.json()["detail"]
                assert runtime.jobs.list(project_id) == []
            # Without an explicit override, the shared boundary entry records
            # its actual saved comparison even when later input checks fail.
            boundary = next(mode for mode in modes if mode["mode_id"] == "boundary")
            for invalid_scope in (1, {"unexpected": comparison_a}, [{}]):
                invalid_parameters = {**boundary["fixed_parameters"], "selection_id": selection_id, "contrast_ids": invalid_scope}
                invalid_plan = client.post(f"/api/projects/{project_id}/workflows/{boundary['workflow_id']}/plan",
                                           json={"parameters": invalid_parameters, "force": False})
                assert invalid_plan.status_code == 200, invalid_plan.text
                assert not invalid_plan.json()["ready"]
                assert invalid_plan.json()["errors"]
                invalid_run = client.post(f"/api/projects/{project_id}/workflows/{boundary['workflow_id']}/run",
                                          json={"parameters": invalid_parameters, "force": False})
                assert invalid_run.status_code == 422, invalid_run.text
                assert invalid_run.json()["detail"]
                assert runtime.jobs.list(project_id) == []
            scoped = client.post(f"/api/projects/{project_id}/workflows/{boundary['workflow_id']}/plan",
                                 json={"parameters": {**boundary["fixed_parameters"], "selection_id": selection_id}, "force": False})
            assert scoped.json()["parameters"]["contrast_ids"] == [comparison_a]
            assert runtime.jobs.list(project_id) == []
            assert source.path.read_bytes() == original_source
    finally:
        runtime.jobs.shutdown()
    assert before == {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in before}
