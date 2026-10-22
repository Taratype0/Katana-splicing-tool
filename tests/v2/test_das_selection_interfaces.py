"""Public controls and report selection must preserve saved scientific identity."""
import hashlib
import json
import sys
from types import SimpleNamespace

import pandas as pd
import pytest

from katana_next.client import KatanaProject
from katana_next.domain.models import ArtifactKind
from katana_next.science.bulk import READ_SUPPORT_POLICY_VERSION
from katana_next.science.das_selection import selection_metadata
from katana_next.workflows.reporting import _report_artifact_matches, _das_top10_frames, _das_top10_html


def test_python_explicit_selection_rejects_conflicts_without_mutating_parameters():
    original = {"event_set": "stable_significant", "quality_artifact_id": "saved"}
    assert KatanaProject._with_event_set(original, None) == original
    with pytest.raises(ValueError, match="Conflicting event_set"):
        KatanaProject._with_event_set(original, "significant")
    assert original == {"event_set": "stable_significant", "quality_artifact_id": "saved"}
    with pytest.raises(ValueError, match="Choose stable_significant"):
        KatanaProject._with_event_set({}, "anything")


def test_public_workflow_api_exposes_real_selection_capabilities(tmp_path):
    from fastapi.testclient import TestClient
    from katana_next.api.app import PlatformRuntime, create_app
    runtime = PlatformRuntime(application_root=tmp_path / "api")
    try:
        client = TestClient(create_app(runtime))
        response = client.get("/api/workflows")
        assert response.status_code == 200
        definitions = {item["workflow_id"]: item for item in response.json()}
        for workflow in ("differential.as", "differential.cross_comparison", "prioritize.enrichment", "reports.bundle"):
            item = definitions[workflow]
            assert item["das_selection"]["mode"] == "switchable"
            parameter = next(parameter for parameter in item["parameters"] if parameter["key"] == "event_set")
            assert parameter["default"] == "stable_significant"
            assert parameter["options"] == ["stable_significant", "significant"]
        assert definitions["prioritize.event_motifs"]["das_selection"]["mode"] == "rs_required"
        assert definitions["programs.cosplice"]["das_selection"]["mode"] == "measurement_population"
    finally:
        runtime.jobs.shutdown()


@pytest.mark.parametrize("action", ["plan", "run"])
def test_cli_event_set_reaches_the_same_parameter_payload(monkeypatch, capsys, tmp_path, action):
    from katana_next import cli
    calls = []
    class FakeProject:
        def __init__(self, path):
            assert path.name == "project.json"
        def __enter__(self):
            return self
        def __exit__(self, *args):
            pass
        def plan(self, workflow, parameters):
            calls.append((workflow, parameters))
            return {"ready": True, "parameters": parameters}
        def run(self, workflow, parameters, force=False):
            calls.append((workflow, parameters))
            return SimpleNamespace(model_dump=lambda **_: {"status": "complete", "parameters": parameters})
    monkeypatch.setattr(cli, "KatanaProject", FakeProject)
    monkeypatch.setattr(sys, "argv", ["katana-next", "--project", str(tmp_path / "project.json"), action, "differential.as", "--event-set", "significant"])
    cli.main()
    assert calls == [("differential.as", {"event_set": "significant"})]
    assert json.loads(capsys.readouterr().out)["parameters"]["event_set"] == "significant"


def test_report_scope_does_not_relabel_unknown_or_other_selection():
    def artifact(workflow, policy=None):
        return SimpleNamespace(analysis_id=workflow, metadata={} if policy is None else {"das_selection": policy})
    for workflow in ("differential.as", "prioritize.enrichment", "reports.bundle"):
        assert not _report_artifact_matches(artifact(workflow), "stable_significant")
        assert _report_artifact_matches(artifact(workflow, selection_metadata("significant")), "significant")
        assert not _report_artifact_matches(artifact(workflow, selection_metadata("stable_significant")), "significant")
    assert _report_artifact_matches(artifact("prioritize.event_motifs"), "significant")
    assert _report_artifact_matches(artifact("programs.cosplice"), "significant")


def test_report_top10_switches_actual_membership_preserving_source_flags_statistics_and_bytes(tmp_path):
    frame = pd.DataFrame({
        "event_id": [f"E{i:02}" for i in range(25)], "comparison_id": ["A"] * 13 + ["B"] * 12,
        "gene_symbol": [f"G{i}" for i in range(25)], "gene_id": [f"ID{i}" for i in range(25)],
        "event_type": ["SE"] * 25, "fdr": [.001 + i * .0001 for i in range(25)],
        "dpsi": [.2] * 25, "significant": [True] * 25,
        "stable_significant": [i % 2 == 0 for i in range(25)],
        "measurement_eligible": [i % 2 == 0 for i in range(25)],
        "read_support_policy": [READ_SUPPORT_POLICY_VERSION] * 25,
    })
    path = tmp_path / "full.parquet"
    frame.to_parquet(path)
    before = hashlib.sha256(path.read_bytes()).hexdigest()
    source = SimpleNamespace(name="Differential AS events", kind=ArtifactKind.TABLE, path=path, run_fingerprint="saved-original")
    ordinary = _das_top10_frames([source], "significant")
    supported = _das_top10_frames([source], "stable_significant")
    assert [len(top) for _, top, _ in ordinary] == [10, 10]
    assert [len(top) for _, top, _ in supported] == [7, 6]
    assert ordinary[0][1].event_id.tolist() == [f"E{i:02}" for i in range(10)]
    assert supported[0][1].event_id.tolist() == [f"E{i:02}" for i in range(0, 13, 2)]
    assert "Top 7 RS-DAS events | A" in _das_top10_html([source], "stable_significant")
    assert "Top 10 DAS events | A" in _das_top10_html([source], "significant")
    assert "event fields" in _das_top10_html([source], "significant")
    assert hashlib.sha256(path.read_bytes()).hexdigest() == before
    pd.testing.assert_frame_equal(frame, pd.read_parquet(path))


def test_real_evidence_export_separates_saved_foregrounds_without_reanalysis(tmp_path):
    from katana_next.full_fixture import create_full_fixture
    from katana_next.services.project_service import ProjectService
    from katana_next.workflows.base import ArtifactWriter
    from katana_next.workflows.registry import WorkflowRegistry
    from katana_next.workflows.reporting import ReproducibilityReportWorkflow

    manifest = create_full_fixture(tmp_path / "project")
    context = ProjectService(application_root=tmp_path / "state").open_project(manifest)
    registry = WorkflowRegistry()
    sources = []
    for workflow, selection, label in (
        ("prioritize.enrichment", "stable_significant", "RS_foreground"),
        ("prioritize.enrichment", "significant", "DAS_foreground"),
        ("prioritize.event_motifs", "stable_significant", "RS_motif_population"),
        ("programs.cosplice", None, "measured_PSI_population"),
    ):
        parameters = {"event_set": selection} if selection else {}
        writer = ArtifactWriter(context, workflow, "saved-" + label, parameters, registry.get(workflow).version,
            common_metadata={"question": label, "das_selection": selection_metadata(selection)} if selection else {"question": label})
        sources.append(writer.table(label, label + ".parquet", pd.DataFrame({"row_id": [label], "gene": [label], "fdr": [.02]}), entity_type="toy_result"))
        writer.finish()
    for selection in ("stable_significant", "significant"):
        writer = ArtifactWriter(context, "differential.as", "saved-as-" + selection,
            {"event_set": selection}, registry.get("differential.as").version,
            common_metadata={"question": "Saved full DAS table", "read_support_policy": READ_SUPPORT_POLICY_VERSION,
                             "das_selection": selection_metadata(selection)})
        sources.append(writer.table("Differential AS events", "events.parquet", pd.DataFrame({
            "row_id": ["a1", "a2", "b1", "b2"], "event_id": ["a1", "a2", "b1", "b2"],
            "comparison_id": ["A", "A", "B", "B"], "gene_symbol": ["X", "Y", "Z", "W"],
            "fdr": [.01, .02, .01, .02], "dpsi": [.2, .3, .2, .3], "significant": [True] * 4,
            "stable_significant": [True, False, True, False], "measurement_eligible": [True, False, True, False],
            "read_support_policy": [READ_SUPPORT_POLICY_VERSION] * 4,
        }), entity_type="splice_event"))
        writer.finish()
    before = {item.path: hashlib.sha256(item.path.read_bytes()).hexdigest() for item in sources}
    report = ReproducibilityReportWorkflow()
    for selection, included, excluded in (
        ("stable_significant", "RS_foreground", "DAS_foreground"),
        ("significant", "DAS_foreground", "RS_foreground"),
    ):
        run = "report-" + selection
        report.run(context, SimpleNamespace(report=lambda *_: None), run,
                   {"event_set": selection, "latest_only": True, "include_result_files": True})
        outputs = context.catalog.list_artifacts(context.manifest.project_id, analysis_id="reports.bundle")
        outputs = [item for item in outputs if item.run_fingerprint == run]
        page = next(item.path for item in outputs if item.path.name == "KATANA_EVIDENCE_REPORT.html")
        content = page.read_text(encoding="utf-8")
        assert included in content and excluded not in content
        assert "RS-DAS required for new analysis" in content
        assert "measured PSI / read-support population" in content
        assert "not computed for this view" in content
        assert ("Top 1 RS-DAS events | A" if selection == "stable_significant" else "Top 2 DAS events | A") in content
        provenance = json.loads(next(item.path for item in outputs if item.path.name == "provenance.json").read_text(encoding="utf-8"))
        assert provenance["das_report_selection"]["event_set"] == selection
    assert {item.path: hashlib.sha256(item.path.read_bytes()).hexdigest() for item in sources} == before
