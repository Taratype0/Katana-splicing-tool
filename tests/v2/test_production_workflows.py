from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
import zipfile

import pandas as pd
from fastapi.testclient import TestClient

from katana_next.api.app import PlatformRuntime, create_app
from katana_next.domain.models import JobStatus
from katana_next.full_fixture import create_full_fixture
from katana_next.science.bulk import canonical_event_id, canonical_event_ids
from katana_next.science.isoform import dominant_isoform_switches
from katana_next.science.single_cell import (
    _annotate_junction_choices,
    _splice_site_components,
)
from katana_next.workflows.registry import WorkflowRegistry


def test_registry_exposes_all_nine_workspaces_and_production_workflows() -> None:
    definitions = WorkflowRegistry().definitions()
    assert len(definitions) == 17
    assert {item["workspace_id"] for item in definitions} == {
        "project",
        "quality",
        "landscape",
        "differential",
        "prioritize",
        "validate",
        "programs",
        "single-cell",
        "reports",
    }
    assert all(item["enabled"] for item in definitions)
    assert all(item["question"] and item["method"] for item in definitions)
    assert {
        item["workflow_id"]
        for item in definitions
        if item["supports_contrast_selection"]
    } == {
        "quality.stable_as",
        "landscape.bulk",
        "differential.as",
        "differential.deg",
        "differential.dtu",
        "differential.cross_comparison",
        "prioritize.evidence",
    }


def test_generic_run_is_isolated_cached_and_downloadable(tmp_path: Path) -> None:
    manifest = create_full_fixture(tmp_path / "fixture")
    runtime = PlatformRuntime(application_root=tmp_path / "application")
    client = TestClient(create_app(runtime))
    opened = client.post("/api/projects/open", json={"manifest_path": str(manifest)})
    assert opened.status_code == 200, opened.text
    project_id = opened.json()["snapshot"]["project_id"]

    started = client.post(
        f"/api/projects/{project_id}/workflows/differential.deg/run",
        json={"parameters": {"fdr": 0.05, "min_abs_log2fc": 1.0}, "force": False},
    )
    assert started.status_code == 200, started.text
    finished = runtime.jobs.wait(started.json()["job_id"], timeout=60)
    assert finished.status == JobStatus.COMPLETE

    artifacts = client.get(f"/api/projects/{project_id}/artifacts").json()
    assert {item["analysis_id"] for item in artifacts} == {"differential.deg"}
    assert len(artifacts) == 3
    assert all(
        item["metadata"].get("project_snapshot_fingerprint")
        == opened.json()["snapshot"]["manifest_fingerprint"]
        for item in artifacts
    )
    table = next(item for item in artifacts if item["kind"] == "table")
    download = client.get(
        f"/api/projects/{project_id}/artifacts/{table['artifact_id']}/file"
    )
    assert download.status_code == 200
    assert download.content

    cached = client.post(
        f"/api/projects/{project_id}/workflows/differential.deg/run",
        json={"parameters": {"fdr": 0.05, "min_abs_log2fc": 1.0}, "force": False},
    )
    assert cached.status_code == 200
    assert cached.json()["status"] == "cache_hit"

    compatible = client.post(
        f"/api/projects/{project_id}/workflows/differential.deg/cache",
        json={
            "parameters": {"fdr": 0.05, "min_abs_log2fc": 1.0},
            "force": False,
        },
    )
    assert compatible.status_code == 200
    assert compatible.json()["cache_hit"] is True
    assert len(compatible.json()["artifacts"]) == 3

    integer_encoded_float = client.post(
        f"/api/projects/{project_id}/workflows/differential.deg/cache",
        json={
            "parameters": {"fdr": 0.05, "min_abs_log2fc": 1},
            "force": False,
        },
    )
    assert integer_encoded_float.status_code == 200
    assert (
        integer_encoded_float.json()["run_fingerprint"]
        == compatible.json()["run_fingerprint"]
    )
    assert integer_encoded_float.json()["cache_hit"] is True

    all_contrasts = [
        contrast["contrast_id"] for contrast in opened.json()["manifest"]["contrasts"]
    ]
    explicit_all = client.post(
        f"/api/projects/{project_id}/workflows/differential.deg/cache",
        json={
            "parameters": {
                "fdr": 0.05,
                "min_abs_log2fc": 1.0,
                "contrast_ids": list(reversed(all_contrasts)),
            },
            "force": False,
        },
    )
    assert explicit_all.status_code == 200
    assert explicit_all.json()["run_fingerprint"] == compatible.json()["run_fingerprint"]
    assert explicit_all.json()["cache_hit"] is True

    empty_scope = client.post(
        f"/api/projects/{project_id}/workflows/differential.deg/cache",
        json={
            "parameters": {
                "fdr": 0.05,
                "min_abs_log2fc": 1.0,
                "contrast_ids": [],
            },
            "force": False,
        },
    )
    assert empty_scope.status_code == 422

    unsupported_scope = client.post(
        f"/api/projects/{project_id}/workflows/prioritize.isoform_support/cache",
        json={"parameters": {"contrast_ids": all_contrasts}, "force": False},
    )
    assert unsupported_scope.status_code == 422

    changed = client.post(
        f"/api/projects/{project_id}/workflows/differential.deg/cache",
        json={
            "parameters": {"fdr": 0.01, "min_abs_log2fc": 1.0},
            "force": False,
        },
    )
    assert changed.status_code == 200
    assert changed.json()["cache_hit"] is False
    assert changed.json()["artifacts"] == []


def test_bulk_and_single_cell_core_workflows_on_complete_fixture(tmp_path: Path) -> None:
    manifest = create_full_fixture(tmp_path / "fixture")
    runtime = PlatformRuntime(application_root=tmp_path / "application")
    context = runtime.projects.open_project(manifest)
    project_id = context.manifest.project_id
    workflow_ids = (
        "project.audit",
        "quality.stable_as",
        "landscape.bulk",
        "differential.as",
        "differential.deg",
        "differential.dtu",
        "differential.cross_comparison",
        "prioritize.evidence",
        "prioritize.isoform_support",
        "single_cell.audit",
        "single_cell.quantify",
        "single_cell.differential",
        "reports.bundle",
    )
    for workflow_id in workflow_ids:
        parameters = (
            {"numerator": "Perturbation", "denominator": "Baseline"}
            if workflow_id == "single_cell.differential"
            else {}
        )
        submitted = runtime.workflows.run(workflow_id, project_id, parameters=parameters)
        finished = runtime.jobs.wait(submitted.job_id, timeout=90)
        assert finished.status in {JobStatus.COMPLETE, JobStatus.CACHE_HIT}, (
            workflow_id,
            finished.error,
        )
    observed = {item.analysis_id for item in context.catalog.list_artifacts(project_id)}
    assert set(workflow_ids).issubset(observed)
    landscape_artifacts = context.catalog.list_artifacts(
        project_id, analysis_id="landscape.bulk"
    )
    assert landscape_artifacts
    assert all(
        item.metadata.get("source_cache_status") == "hit"
        for item in landscape_artifacts
    )
    archive = next(
        item
        for item in context.catalog.list_artifacts(
            project_id, analysis_id="reports.bundle"
        )
        if item.path.suffix.lower() == ".zip"
    )
    with zipfile.ZipFile(archive.path) as bundle:
        names = bundle.namelist()
    assert len(names) == len(set(names))


def test_splice_site_components_join_shared_donors_and_acceptors() -> None:
    choices = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1", "chr1"],
            "start": [100, 100, 150, 900],
            "end": [200, 250, 250, 950],
            "strand": ["+", "+", "+", "+"],
        }
    )
    labels = _splice_site_components(choices)
    assert labels.iloc[0] == labels.iloc[1] == labels.iloc[2]
    assert labels.iloc[3] != labels.iloc[0]


def test_vectorized_event_ids_match_scalar_contract() -> None:
    frame = pd.DataFrame(
        {
            "ID": [7, 8],
            "chr": ["chr1", ""],
            "strand": ["+", ""],
            "exonStart_0base": [100, None],
            "exonEnd": [200, None],
        }
    )
    scalar = frame.apply(lambda row: canonical_event_id(row, "SE"), axis=1)
    vectorized = canonical_event_ids(frame, "SE")
    assert vectorized.tolist() == scalar.tolist()


def test_dominant_isoform_switches_select_arm_means_without_gene_loops() -> None:
    context = SimpleNamespace(
        manifest=SimpleNamespace(
            contrasts=(
                SimpleNamespace(
                    contrast_id="treated_vs_control",
                    numerator=SimpleNamespace(sample_ids=("T1", "T2")),
                    denominator=SimpleNamespace(sample_ids=("C1", "C2")),
                ),
            )
        )
    )
    transcript = pd.DataFrame(
        [
            {"sample_id": sample, "gene_id": "G1", "gene_symbol": "Gene1", "transcript_id": transcript_id, "tpm": tpm}
            for sample, transcript_id, tpm in (
                ("T1", "TX1", 2.0),
                ("T1", "TX2", 8.0),
                ("T2", "TX1", 4.0),
                ("T2", "TX2", 6.0),
                ("C1", "TX1", 9.0),
                ("C1", "TX2", 1.0),
                ("C2", "TX1", 7.0),
                ("C2", "TX2", 3.0),
            )
        ]
    )
    observed = dominant_isoform_switches(context, transcript)
    assert len(observed) == 1
    assert observed.loc[0, "numerator_dominant_transcript"] == "TX2"
    assert observed.loc[0, "denominator_dominant_transcript"] == "TX1"
    assert observed.loc[0, "numerator_dominant_tpm"] == 7.0
    assert observed.loc[0, "denominator_dominant_tpm"] == 8.0


def test_single_cell_junctions_use_exact_and_locus_gtf_annotation(
    tmp_path: Path,
) -> None:
    gtf = tmp_path / "annotation.gtf"
    attributes = (
        'gene_id "GENE1"; transcript_id "TX1"; gene_name "ExampleGene";'
    )
    gtf.write_text(
        "chr1\ttest\texon\t1\t100\t.\t+\t.\t" + attributes + "\n"
        "chr1\ttest\texon\t201\t300\t.\t+\t.\t" + attributes + "\n",
        encoding="utf-8",
    )
    choices = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr2"],
            "start": [100, 110, 100],
            "end": [200, 190, 200],
            "strand": ["+", "+", "+"],
        }
    )
    observed = _annotate_junction_choices(choices, gtf)
    assert observed.loc[0, "gene_symbol"] == "ExampleGene"
    assert bool(observed.loc[0, "annotated"]) is True
    assert observed.loc[1, "gene_id"] == "GENE1"
    assert bool(observed.loc[1, "annotated"]) is False
    assert observed.loc[2, "gene_id"] == ""
