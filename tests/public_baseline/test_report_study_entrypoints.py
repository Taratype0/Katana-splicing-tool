"""Public report CLI, registry and workflow registration contracts."""
from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace
import zipfile

import pytest

from katana_next.domain.models import AnalysisArtifact, ArtifactKind, ProjectManifest
from katana_next.infrastructure.catalog import ArtifactCatalog
from katana_next.infrastructure.result_store import ResultStore
from katana_next.reports.cli import main
from katana_next.workflows.report_study import StudyReportWorkflow


def _project(tmp_path):
    root = tmp_path / "source-results"
    root.mkdir()
    image = root / "sample_pca.svg"
    image.write_text('<svg xmlns="http://www.w3.org/2000/svg" width="240" height="140">'
                     '<text x="10" y="30">Existing PCA display fixture</text></svg>', encoding="utf-8")
    manifest = ProjectManifest(project_id="report-fixture", name="Completed report fixture",
                               input_root=tmp_path / "input", output_root=root,
                               annotation={"genome_build": "test-build"})
    record = AnalysisArtifact(project_id=manifest.project_id, analysis_id="landscape.bulk",
                              name="Sample PCA", kind=ArtifactKind.FILE, entity_type="plot",
                              fingerprint="existing-source-digest", run_fingerprint="fixture-run",
                              path=image, metadata={"project_snapshot_fingerprint": "snapshot"})
    project = tmp_path / "katana.project.json"
    project.write_text(manifest.model_dump_json(indent=2), encoding="utf-8")
    return manifest, project, record


def _descriptor(tmp_path):
    manifest, project, artifact = _project(tmp_path)
    descriptor = tmp_path / "source-descriptor.json"
    descriptor.write_text(json.dumps({"project": manifest.model_dump(mode="json"),
                                      "artifacts": [artifact.model_dump(mode="json")]}), encoding="utf-8")
    return descriptor, artifact


def test_cli_defaults_to_rs_without_opening_project_service(tmp_path, capsys, monkeypatch):
    import katana_next.services.project_service
    def forbidden(*args, **kwargs):
        raise AssertionError("Read-only report CLI must not open/migrate the source project")
    monkeypatch.setattr(katana_next.services.project_service.ProjectService, "open_project", forbidden)
    descriptor, source = _descriptor(tmp_path)
    before = (source.path.stat().st_mtime_ns, source.path.read_bytes())
    output = tmp_path / "portable-report"
    assert main(["--project", str(descriptor), "--outdir", str(output)]) == 0
    result = json.loads(capsys.readouterr().out)
    assert result["status"] == "completed" and result["event_set"] == "stable_significant"
    assert result["no_upstream_rerun"] and result["no_source_statistics_refit"]
    assert Path(result["report"]).is_file() and Path(result["archive"]).is_file()
    assert before == (source.path.stat().st_mtime_ns, source.path.read_bytes())
    page = (output / "REPORT.html").read_text(encoding="utf-8")
    assert all(f'id="chapter-{number}"' in page for number in range(1, 8))
    with zipfile.ZipFile(output / "REPORT.zip") as archive:
        assert archive.testzip() is None
        assert {"REPORT.html", "report-manifest.json", "provenance.json", "assets/report.css",
                "assets/report.js", "outputs/index.html"} <= set(archive.namelist())


@pytest.mark.parametrize("option", ["--catalog", "--asset-catalog"])
def test_catalog_cli_alias_and_ordinary_das(tmp_path, capsys, option):
    manifest, project, artifact = _project(tmp_path)
    catalog = tmp_path / "completed-artifacts.json"
    catalog.write_text(json.dumps([artifact.model_dump(mode="json")]), encoding="utf-8")
    out = tmp_path / "portable"
    assert main(["--project", str(project), "--outdir", str(out), option, str(catalog),
                 "--event-set", "significant"]) == 0
    result = json.loads(capsys.readouterr().out)
    assert result["event_set"] == "significant"
    assert result["figures"] >= 1


def test_cli_does_not_overwrite_completed_report(tmp_path, capsys):
    descriptor, _ = _descriptor(tmp_path)
    out = tmp_path / "report"
    assert main(["--project", str(descriptor), "--outdir", str(out)]) == 0
    capsys.readouterr()
    original = (out / "REPORT.html").read_bytes()
    assert main(["--project", str(descriptor), "--outdir", str(out)]) == 1
    streams = capsys.readouterr()
    assert not streams.out
    error = json.loads(streams.err)
    assert error["status"] == "failed" and "preserved" in error["error"]
    assert (out / "REPORT.html").read_bytes() == original


def test_module_cli_help_does_not_import_app_or_start_service():
    result = subprocess.run([sys.executable, "-m", "katana_next.reports.cli", "--help"],
                            capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    assert all(option in result.stdout for option in ("--project", "--outdir", "--event-set", "--catalog", "--asset-catalog"))


def test_registered_study_workflow_keeps_old_report_and_rs_switch():
    from katana_next.science.das_selection import workflow_das_capability
    from katana_next.workflows.registry import WorkflowRegistry
    registry = WorkflowRegistry()
    workflow = registry.get("reports.study")
    assert isinstance(workflow, StudyReportWorkflow)
    assert workflow.definition.workspace_id == "reports"
    assert workflow.normalize_parameters({})["event_set"] == "stable_significant"
    assert workflow.normalize_parameters({"event_set": "significant"})["event_set"] == "significant"
    assert workflow_das_capability("reports.study")["mode"] == "switchable"
    assert registry.get("reports.bundle").version == "2.3.0+das_selection.20261022"


def test_workflow_registers_finished_report_and_preserves_parent(tmp_path):
    manifest, project, parent = _project(tmp_path)
    catalog = ArtifactCatalog(manifest.output_root)
    catalog.register_artifact(parent)
    context = SimpleNamespace(manifest=manifest, manifest_path=project,
                              snapshot=SimpleNamespace(manifest_fingerprint="snapshot"),
                              catalog=catalog, results=ResultStore(manifest.output_root))
    # ArtifactWriter.finish exports the catalogue; avoid constructing/opening a
    # real application ProjectService just to register this synthetic fixture.
    catalog.export_manifest = lambda project_id: None
    progress = []
    job = SimpleNamespace(check_cancelled=lambda: None, report=lambda value, label: progress.append((value, label)))
    before = (parent.path.stat().st_mtime_ns, parent.path.read_bytes())
    workflow = StudyReportWorkflow()
    parameters = workflow.normalize_parameters({})
    ids = workflow.run(context, job, "new-report-fixture-run", parameters)
    assert len(ids) == 4
    outputs = [catalog.get_artifact(identifier) for identifier in ids]
    assert all(output.analysis_id == "reports.study" for output in outputs)
    assert {output.path.name for output in outputs} == {"REPORT.html", "REPORT.zip", "provenance.json", "report-manifest.json"}
    report = next(output for output in outputs if output.path.name == "REPORT.html")
    assert report.metadata["preserve_relative_report_links"]
    assert "assets/report.css" in report.metadata["exports"]
    assert parent.artifact_id in report.parent_artifact_ids
    marker = json.loads((manifest.output_root / "artifacts/reports.study/new-report-fixtu/parameters.json").read_text())
    assert marker["completed"]
    assert before == (parent.path.stat().st_mtime_ns, parent.path.read_bytes())
    assert progress[-1][0] == 1.0


def test_workflow_rejects_unfinished_and_foreign_snapshot_sources(tmp_path):
    manifest, project, parent = _project(tmp_path)
    in_progress = parent.model_copy(update={"status": "running"})
    stale = parent.model_copy(update={"metadata": {"project_snapshot_fingerprint": "other"}})
    context = SimpleNamespace(manifest=manifest, snapshot=SimpleNamespace(manifest_fingerprint="snapshot"),
                              catalog=SimpleNamespace(list_artifacts=lambda project_id: [in_progress, stale]))
    assert StudyReportWorkflow().parent_artifacts(context, {}) == []
