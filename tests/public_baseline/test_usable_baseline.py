"""Portable acceptance of the complete application, using generated toy inputs.

These checks preserve the usable workflow inventory and exercise the public
application interfaces. They do not claim that external R tools, reference
databases or a real experimental study were executed.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import zipfile

from fastapi.testclient import TestClient
import pandas as pd
import pytest

from katana_next import __version__
from katana_next.api.app import PlatformRuntime, create_app
from katana_next.domain.models import JobStatus
from katana_next.full_fixture import create_full_fixture
from katana_next.workflows.registry import WorkflowRegistry


ROOT = Path(__file__).resolve().parents[2]
INVENTORY = json.loads(Path(__file__).with_name("workflow_inventory_read_support_v1.json").read_text(encoding="utf-8"))
# The earlier inventory remains immutable; this approved network policy changes
# the count-bearing network method version, with its parameters unchanged.
for entry in INVENTORY["workflows"]:
    if entry["workflow_id"] == "programs.cosplice":
        entry["version"] = "3.0.0"


def test_workflows_match_approved_all_sample_policy_and_retained_other_methods():
    actual = [
        {
            "workflow_id": workflow.definition.workflow_id,
            "version": workflow.version,
            "parameter_defaults": {p.key: p.default for p in workflow.definition.parameters},
        }
        for workflow in WorkflowRegistry().all()
    ]
    assert len(actual) == INVENTORY["workflow_count"] == 54
    assert actual == INVENTORY["workflows"]


def test_only_approved_read_support_workflows_changed_from_preserved_baseline():
    previous = json.loads(Path(__file__).with_name("workflow_inventory.json").read_text(encoding="utf-8"))
    old = {entry["workflow_id"]: entry for entry in previous["workflows"]}
    current = {entry["workflow_id"]: entry for entry in INVENTORY["workflows"]}
    assert old.keys() == current.keys()
    assert {key for key in old if old[key] != current[key]} == {
        "quality.stable_as", "landscape.bulk", "differential.as", "differential.cross_comparison", "programs.cosplice",
    }
    assert current["quality.stable_as"]["parameter_defaults"] == {"fdr": 0.05, "min_abs_dpsi": 0.1}
    assert current["differential.cross_comparison"]["parameter_defaults"] == {
        "membership_mode": "exact_event", "event_set": "stable_significant",
    }
    assert current["programs.cosplice"]["parameter_defaults"] == old["programs.cosplice"]["parameter_defaults"]


@pytest.mark.parametrize("workflow_id", [
    "prioritize.rbp_mechanism", "prioritize.rbp_sequence_results",
    "prioritize.motiftrace_consequences", "prioritize.motiftrace_usage",
    "prioritize.motiftrace_parallel", "prioritize.event_motifs",
    "programs.cosplice", "programs.psi_network", "programs.psi_projection",
    "programs.djexpress", "programs.module_functions", "programs.hub_events",
    "programs.explain_expression", "programs.diagnostics",
    "programs.external_validation", "reports.bundle", "reports.interpreted_bundle",
])
def test_key_scientific_and_report_entrypoints_are_available(workflow_id):
    workflow = WorkflowRegistry().get(workflow_id)
    assert workflow.definition.enabled
    assert callable(workflow.run)
    assert workflow.definition.method
    assert workflow.version


@pytest.mark.parametrize("relative", [
    "katana_next/science/rbp_sequence_browser.html",
    "katana_next/science/motiftrace_windows.html",
    "katana_next/science/motiftrace_consequences.html",
    "katana_next/science/motiftrace_parallel.html",
    "katana_next/science/rbp_evidence_layers.json",
    "katana_next/science/dexseq_interaction.R",
    "katana_next/science/expression_response.R",
    "katana_next/science/hub_metrics.R",
    "katana_next/science/djexpress_driver.R",
    "katana_next/vendor/DJExpress/DESCRIPTION",
    "katana_next/vendor/DJExpress/LICENSE",
    "katana_next/vendor/DJExpress/R/JCNA1pass.R",
    "assets/katana_icon.ico",
    "apps/web/dist/index.html",
    "src/analysis/gene_filter.py",
    "src/domain/models/gene_filter.py",
    "src/adapters/native_sashimi_adapter.py",
])
def test_runtime_resources_and_retained_legacy_dependencies_exist(relative):
    assert (ROOT / relative).is_file(), relative
    assert (ROOT / relative).stat().st_size > 0, relative


def test_gene_screening_does_not_become_transcript_confirmation(tmp_path):
    from katana_next.science.bulk import load_dtu_results

    path = tmp_path / "stage.tsv"
    pd.DataFrame([{"geneID": "G", "txID": "T", "gene": 0.001, "transcript": float("nan")}]).to_csv(path, sep="\t", index=False)
    context = SimpleNamespace(manifest=SimpleNamespace(
        metadata={}, assets=[SimpleNamespace(kind="dtu_result", path=path, metadata={"contrast_id": "A_vs_B"})],
    ))
    row = load_dtu_results(context, ["A_vs_B"]).iloc[0]
    assert pd.isna(row.dtu_fdr)
    assert row.gene_screen_fdr == 0.001


def test_api_opens_a_synthetic_project_runs_a_cached_result_and_packages_english_reports(tmp_path):
    manifest = create_full_fixture(tmp_path / "synthetic-project")
    inputs = manifest.parent / "inputs"
    original_inputs = {
        path.relative_to(inputs).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in inputs.rglob("*") if path.is_file()
    }
    runtime = PlatformRuntime(application_root=tmp_path / "application-state")
    try:
        with TestClient(create_app(runtime)) as client:
            assert client.get("/api/health").json() == {"status": "ok", "version": __version__}
            assert len(client.get("/api/workflows").json()) == 54
            assert {w["workspace_id"] for w in client.get("/api/workspaces").json()} == {
                "project", "quality", "landscape", "differential", "prioritize",
                "validate", "programs", "single-cell", "reports",
            }
            opened = client.post("/api/projects/open", json={"manifest_path": str(manifest)})
            assert opened.status_code == 200, opened.text
            project_id = opened.json()["snapshot"]["project_id"]
            parameters = {"fdr": 0.05, "min_abs_log2fc": 1.0}

            def run(workflow_id, supplied):
                response = client.post(
                    f"/api/projects/{project_id}/workflows/{workflow_id}/run",
                    json={"parameters": supplied, "force": False},
                )
                assert response.status_code == 200, response.text
                result = runtime.jobs.wait(response.json()["job_id"], timeout=60)
                assert result.status in {JobStatus.COMPLETE, JobStatus.CACHE_HIT}, result.error
                return result

            first = run("differential.deg", parameters)
            assert first.status == JobStatus.COMPLETE
            repeated = run("differential.deg", parameters)
            assert repeated.status == JobStatus.CACHE_HIT
            assert first.run_fingerprint == repeated.run_fingerprint
            report = run("reports.bundle", {"include_result_files": True, "latest_only": True})
            assert report.status == JobStatus.COMPLETE
            artifacts = runtime.projects.list_artifacts(project_id, analysis_id="reports.bundle")
            assert any(a.path.suffix == ".html" for a in artifacts)
            assert any(a.path.suffix == ".pdf" for a in artifacts)
            package = next(a for a in artifacts if a.path.suffix == ".zip")
            with zipfile.ZipFile(package.path) as archive:
                assert archive.testzip() is None
                names = archive.namelist()
                assert "report/KATANA_EVIDENCE_REPORT.html" in names
                assert "report/KATANA_EVIDENCE_REPORT.pdf" in names
                assert "provenance/provenance.json" in names
                assert not any(name.lower().endswith((".bam", ".fastq", ".fastq.gz")) for name in names)
                assert not any("REPORT_ZH" in name for name in names)
                assert any(name.startswith("artifacts/differential.deg/") for name in names)
            assert {j.analysis_id for j in runtime.jobs.list(project_id)} == {"differential.deg", "reports.bundle"}
            downloaded = client.get(f"/api/projects/{project_id}/artifacts/{package.artifact_id}/file")
            assert downloaded.status_code == 200
            assert downloaded.content == package.path.read_bytes()
    finally:
        runtime.jobs.shutdown()
    assert original_inputs == {
        path.relative_to(inputs).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in inputs.rglob("*") if path.is_file()
    }


def test_one_click_scientific_report_keeps_original_sections_and_english_download(tmp_path):
    """Execute the real report workflow on deliberately partial toy results."""
    from katana_next.science.report_complete import MODULES

    manifest = create_full_fixture(tmp_path / "report-project")
    results = Path(json.loads(manifest.read_text(encoding="utf-8"))["output_root"])
    totals = results / "20_Presentation_report" / "comparison_totals.tsv"
    totals.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{"comparison_id": "perturbation-vs-baseline", "significant_events": 1, "DEG": 1}]).to_csv(totals, sep="\t", index=False)
    intersection = totals.parent / "DEG_DAS_gene_intersection.tsv"
    pd.DataFrame([{"comparison_id": "perturbation-vs-baseline", "gene_id": "toy-gene", "P": 0.031, "FDR": 0.062}]).to_csv(intersection, sep="\t", index=False)
    cards = results / "23_Module_explanation" / "toy-run" / "module_cards.tsv"
    cards.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{"module": "C1", "n_events": 2}]).to_csv(cards, sep="\t", index=False)
    expression = results / "22_Expression_interaction_pathways" / "toy-design" / "toy-run"
    expression.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{"gene_id": "toy-gene", "gene_symbol": "Example", "padj": 0.062}]).to_csv(expression / "interaction.tsv", sep="\t", index=False)
    pd.DataFrame([{"contrast": "interaction", "pathway": "TOY_PATHWAY", "study_family_FDR": 0.2, "NES": -1.2, "leadingEdge": "Example"}]).to_csv(expression / "ranked_pathways.tsv", sep="\t", index=False)
    exon = results / "16_Interaction" / "toy-design" / "toy-run" / "dexseq_interaction.tsv"
    exon.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{"groupID": "toy-gene", "bin_id": "E001", "interaction_FDR": 0.08, "pvalue": 0.04}]).to_csv(exon, sep="\t", index=False)
    config = tmp_path / "report-config.json"
    config.write_text(json.dumps({"study": "Synthetic_example", "result_root": str(results), "manifest_path": str(manifest)}), encoding="utf-8")
    original = {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in (totals, intersection, cards, exon, expression / "interaction.tsv", expression / "ranked_pathways.tsv")}
    runtime = PlatformRuntime(application_root=tmp_path / "report-state")
    try:
        with TestClient(create_app(runtime)) as client:
            opened = client.post("/api/projects/open", json={"manifest_path": str(manifest)})
            assert opened.status_code == 200, opened.text
            project_id = opened.json()["snapshot"]["project_id"]
            response = client.post(
                f"/api/projects/{project_id}/workflows/reports.interpreted_bundle/run",
                json={"parameters": {"config_path": str(config)}, "force": False},
            )
            assert response.status_code == 200, response.text
            job = runtime.jobs.wait(response.json()["job_id"], timeout=60)
            assert job.status == JobStatus.COMPLETE, job.error
            artifacts = runtime.projects.list_artifacts(project_id, analysis_id="reports.interpreted_bundle")
            report = next(artifact.path for artifact in artifacts if artifact.path.name == "REPORT_EN.html")
            text = report.read_text(encoding="utf-8")
            assert len(MODULES) == 15
            for _, _, english_title, _, _ in MODULES:
                from katana_cosplice.branding import display_module_name
                assert display_module_name(english_title) in text
            assert "This is not a zero result." in text
            assert "REPORT_ZH.html" not in text
            assert not (report.parent / "REPORT_ZH.html").exists()
            assert len(json.loads((report.parent / "MODULE_COVERAGE.json").read_text(encoding="utf-8"))) == 15
            archive = next(artifact.path for artifact in artifacts if artifact.path.suffix == ".zip")
            with zipfile.ZipFile(archive) as package:
                assert package.testzip() is None
                assert any(name.endswith("/REPORT_EN.html") for name in package.namelist())
                assert any(name.endswith("/DEG_DAS_gene_intersection.tsv") for name in package.namelist())
                assert not any("REPORT_ZH" in name for name in package.namelist())
                assert not any("/.katana/" in name or name.endswith(".lock") for name in package.namelist())
            assert {j.analysis_id for j in runtime.jobs.list(project_id)} == {"reports.interpreted_bundle"}
    finally:
        runtime.jobs.shutdown()
    assert original == {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in original}
