import json
from zipfile import ZipFile

import numpy as np
import pandas as pd
import pytest

from katana_cosplice.motif_scanner import load_motif_database, motif_database_audit
from katana_next.science.motif_significance import matched_enrichment, measurable_event_background
from katana_next.workflows.registry import WorkflowRegistry


def test_native_cisbp_evidence_and_unavailable_matrices(tmp_path):
    path = tmp_path / "native.zip"
    metadata = "RBP_ID\tRBP_Name\tRBP_Species\tRBP_Status\tMotif_ID\tMotif_Type\n"
    metadata += "T1\tRbp1\tMus_musculus\tD\tM1\tRNAcompete\nT1\tRbp1\tMus_musculus\tD\tM2\tJPLE\n"
    metadata += "T2\tRbp2\tMus_musculus\tI\tM1\tRNAcompete\nT3\tRbp3\tMus_musculus\tD\tM3\tSELEX\n"
    metadata += "T4\tRbp4\tMus_musculus\tD\tM4\tSELEX\n"
    matrix = "Pos\tA\tC\tG\tU\n1\t1\t0\t0\t0\n2\t0\t1\t0\t0\n3\t0\t0\t1\t0\n"
    with ZipFile(path, "w") as archive:
        archive.writestr("RBP_Information_all_motifs.txt", metadata)
        archive.writestr("pwms/M1.txt", matrix)
        archive.writestr("pwms/M2.txt", matrix)
        archive.writestr("pwms/M3.txt", "")
    motifs = load_motif_database(path, organism="Mus_musculus")
    assert {m.evidence_class for m in motifs} == {"direct_experiment", "model_predicted", "homology_inferred"}
    assert len(motifs) == 3 and len({m.motif_id for m in motifs}) == 3
    assert all(np.allclose(m.pwm.sum(axis=1), 1) for m in motifs)
    assert set(motif_database_audit(path, organism="Mus_musculus").status) == {"available", "empty_matrix", "missing_matrix"}
    with pytest.raises(ValueError, match="explicit target species"):
        load_motif_database(path)


def test_motif_evidence_is_not_pooled():
    catalog = pd.DataFrame({"motif_token": ["m1", "m2"], "RBP": ["R", "R"],
                            "evidence_class": ["direct_experiment", "homology_inferred"]})
    pairs = pd.DataFrame({"status": ["matched"] * 3, "foreground_event": ["f1", "f2", "f3"], "background_event": ["b1", "b2", "b3"]})
    sites = pd.DataFrame({"motif_token": ["m1"] * 3 + ["m2"] * 3, "event_id": ["f1", "f2", "f3", "b1", "b2", "b3"]})
    result = matched_enrichment(sites, catalog, pairs).set_index("evidence_class")
    assert len(result) == 2
    assert result.loc["direct_experiment", "enrichment_pvalue"] == .125
    assert result.loc["homology_inferred", "enrichment_pvalue"] == 1
    assert result.loc["direct_experiment", "enrichment_fdr"] == .25


def test_event_background_is_measurement_not_significance():
    frame = pd.DataFrame({"row_id": ["a", "b", "c"], "measurement_eligible": [True, True, False], "significant": [True, False, True]})
    measured, audit = measurable_event_background(frame)
    assert list(measured.row_id) == ["a", "b"]
    assert audit.loc[audit.row_id.eq("c"), "background_reason"].iloc[0] == "event_measurement_support_failed"


def test_registry_primary_and_legacy_roles():
    registry = WorkflowRegistry()
    workflow = registry.get("prioritize.rbp_mechanism")
    assert workflow.definition.presentation_role == "primary"
    assert workflow.normalize_parameters({})["site_pvalue"] == .0001
    assert registry.get("prioritize.rbp_support").definition.presentation_role == "auxiliary"
    assert registry.get("prioritize.event_motifs").definition.presentation_role == "advanced"
    assert "add-on" in registry.get("single_cell.cell_embedding").definition.label


def test_bed_hash_in_peak_name_is_not_a_comment(tmp_path):
    from katana_next.science.versioned_evidence import read_binding_peaks
    path = tmp_path / "real.bed"
    path.write_text("# whole line comment\ntrack name=example\nchr1\t10\t20\tread#0/1_2\t4\t+\n", encoding="ascii")
    frame = read_binding_peaks(path)
    assert frame.shape == (1, 6)
    assert frame.iloc[0, 3] == "read#0/1_2"


def test_rbp_case_aliases_are_one_evidence_hypothesis():
    catalog = pd.DataFrame({"motif_token": ["m1", "m2"], "RBP": ["Pcbp1", "PCBP1"], "evidence_class": ["direct_experiment"] * 2})
    pairs = pd.DataFrame({"status": ["matched"] * 2, "foreground_event": ["f1", "f2"], "background_event": ["b1", "b2"]})
    sites = pd.DataFrame({"motif_token": ["m1", "m2"], "event_id": ["f1", "f2"]})
    result = matched_enrichment(sites, catalog, pairs)
    assert len(result) == 1 and result.iloc[0].foreground_hit_genes == 2
    assert result.iloc[0].RBP == "Pcbp1"


def test_unified_workflow_api_and_complete_child_outputs(tmp_path, monkeypatch):
    from test_supervisor_feedback import _api_catalog_fixture
    from katana_next.domain.models import JobStatus
    from katana_next.resources import sha256
    from katana_next.workflows.base import ArtifactWriter
    from katana_next.workflows.motif import EventMotifWorkflow
    from katana_next.infrastructure.manifest_store import ManifestStore
    runtime, client, context, project, selection, comparisons = _api_catalog_fixture(tmp_path, clip=True)
    manifest_path = context.manifest_path
    payload = context.manifest.model_dump(mode="json")
    motif = tmp_path / "motif.tsv"
    motif.write_text("RBP\tmotif_id\tposition\tA\tC\tG\tT\nR\tm1\t1\t1\t0\t0\t0\n", encoding="ascii")
    runtime_path = tmp_path / "runtime.json"
    runtime_path.write_text(json.dumps({"python": "fixture only"}))
    payload["assets"].extend([
        {"asset_id": "test::motif", "kind": "rbp_motif_database", "path": str(motif), "metadata": {"source_url": "https://example.org/test-only", "version": "test", "species": "Mus_musculus", "sha256": sha256(motif)}},
        {"asset_id": "test::fimo", "kind": "fimo_runtime", "path": str(runtime_path), "metadata": {}},
    ])
    # Only the configured test motif is allowed; no older unversioned fixture resources.
    payload["assets"] = [a for a in payload["assets"] if a["kind"] != "rbp_motif_database" or a["asset_id"] == "test::motif"]
    manifest_path.write_text(json.dumps(payload))
    runtime.projects.open_project(manifest_path)
    context = runtime.projects.get(project)
    old = context.catalog.list_selections(project)[0]
    source = context.catalog.get_artifact(old.source_artifact_id)
    writer = ArtifactWriter(context, "test.measured_event_catalog", "updated", {}, "1.0")
    source = writer.table("Updated test catalog", "events.parquet", pd.read_parquet(source.path), entity_type=source.entity_type)
    writer.finish()
    from katana_next.domain.models import SelectionSet
    saved = runtime.projects.save_selection(SelectionSet(project_id=project, name="New exact selection", entity_type=source.entity_type,
        source_artifact_id=source.artifact_id, entity_ids=old.entity_ids))
    def fake_scan(self, context, job, fingerprint, parameters):
        output = ArtifactWriter(context, self.definition.workflow_id, fingerprint, parameters, self.version)
        output.table("Fixture motif results", "enrichment.parquet", pd.DataFrame({"row_id": ["R"], "RBP": ["R"]}), entity_type="rbp")
        return output.finish()
    monkeypatch.setattr(EventMotifWorkflow, "run", fake_scan)
    started = client.post(f"/api/projects/{project}/workflows/prioritize.rbp_mechanism/run", json={"parameters": {"selection_id": saved.selection_id, "binding_mode": "required"}})
    assert started.status_code == 200, started.text
    finished = runtime.jobs.wait(started.json()["job_id"], timeout=90)
    assert finished.status == JobStatus.COMPLETE, finished.error
    outputs = [a for a in context.catalog.list_artifacts(project) if a.artifact_id in finished.artifact_ids]
    stages = pd.read_parquet(next(a.path for a in outputs if a.path.name == "stages.parquet"))
    assert set(stages.status) == {"complete"}
    assert "Matched-event CLIP enrichment" in {a.name for a in outputs}
    assert "Fixture motif results" in {a.name for a in outputs}
    cached = client.post(f"/api/projects/{project}/workflows/prioritize.rbp_mechanism/run", json={"parameters": {"selection_id": saved.selection_id, "binding_mode": "required"}})
    assert cached.json()["status"] == "cache_hit"
