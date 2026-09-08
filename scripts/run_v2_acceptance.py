from __future__ import annotations

import argparse
import json
import time
import zipfile
from collections import Counter
from pathlib import Path

import pandas as pd

from katana_cosplice.config import load_config
from katana_next.api.app import PlatformRuntime
from katana_next.domain.models import JobStatus, SelectionSet
from katana_next.full_fixture import create_full_fixture


WORKFLOW_RUNS: tuple[tuple[str, dict[str, object]], ...] = (
    ("project.audit", {}),
    ("quality.stable_as", {}),
    ("landscape.bulk", {}),
    ("differential.as", {}),
    ("differential.deg", {}),
    ("differential.dtu", {}),
    ("differential.cross_comparison", {}),
    ("prioritize.evidence", {}),
    ("prioritize.enrichment", {"minimum_overlap": 2}),
    ("prioritize.rbp_support", {"minimum_overlap": 1}),
    ("prioritize.isoform_support", {}),
    (
        "programs.cosplice",
        {"scope": "discovery", "bootstrap_iterations": 10},
    ),
    ("single_cell.audit", {"min_cells_per_donor_cell_type": 10}),
    ("single_cell.quantify", {"source": "precomputed"}),
    (
        "single_cell.differential",
        {
            "numerator": "Perturbation",
            "denominator": "Baseline",
            "min_donors_per_condition": 2,
        },
    ),
    ("reports.bundle", {}),
)


def run_acceptance(root: str | Path) -> Path:
    target = Path(root).expanduser().resolve()
    if target.exists() and any(target.iterdir()):
        raise ValueError(
            f"Acceptance root must be absent or empty, preserving prior runs: {target}"
        )
    target.mkdir(parents=True, exist_ok=True)
    manifest_path = create_full_fixture(target / "project")
    runtime = PlatformRuntime(application_root=Path(__file__).resolve().parents[1])
    context = runtime.projects.open_project(manifest_path)
    project_id = context.manifest.project_id
    records: list[dict[str, object]] = []
    parameter_map: dict[str, dict[str, object]] = {}
    selection_id = ""

    for workflow_id, declared_parameters in WORKFLOW_RUNS:
        parameters = dict(declared_parameters)
        if workflow_id in {"prioritize.enrichment", "prioritize.rbp_support"}:
            if not selection_id:
                selection_id = _create_candidate_selection(context)
            parameters["selection_id"] = selection_id
        before = {
            artifact.artifact_id
            for artifact in context.catalog.list_artifacts(project_id)
        }
        started_at = time.perf_counter()
        submitted = runtime.workflows.run(
            workflow_id,
            project_id,
            parameters=parameters,
        )
        finished = runtime.jobs.wait(submitted.job_id, timeout=900)
        elapsed = time.perf_counter() - started_at
        if finished.status not in {JobStatus.COMPLETE, JobStatus.CACHE_HIT}:
            raise RuntimeError(
                f"{workflow_id} failed: {finished.error or finished.message}"
            )
        after_records = context.catalog.list_artifacts(project_id)
        created = [item for item in after_records if item.artifact_id not in before]
        foreign = sorted({item.analysis_id for item in created} - {workflow_id})
        if foreign:
            raise AssertionError(
                f"{workflow_id} implicitly created artifacts for {foreign}"
            )
        if not created and finished.status != JobStatus.CACHE_HIT:
            raise AssertionError(f"{workflow_id} completed without artifacts")
        _assert_run_files(created)
        _assert_snapshot_binding(created, context.snapshot.manifest_fingerprint)
        parameter_map[workflow_id] = parameters
        records.append(
            {
                "workflow_id": workflow_id,
                "status": finished.status.value,
                "cache_status": finished.cache_status,
                "artifact_count": len(created),
                "elapsed_seconds": round(elapsed, 3),
                "implicit_workflows": foreign,
            }
        )

    all_artifacts = context.catalog.list_artifacts(project_id)
    _assert_run_files(all_artifacts)
    _assert_snapshot_binding(all_artifacts, context.snapshot.manifest_fingerprint)
    _assert_report_archives(all_artifacts)
    context.catalog.export_manifest(project_id)
    resolved_cosplice = (
        context.manifest.output_root
        / "engines"
        / "cosplice"
        / "katana_next_cosplice.yaml"
    )
    cosplice_config = load_config(resolved_cosplice)
    if cosplice_config.consensus.clustering_method != "backbone_refined_consensus":
        raise AssertionError("CoSplice production mode is not bootstrap refined")

    reopened = PlatformRuntime(application_root=Path(__file__).resolve().parents[1])
    reopened_context = reopened.projects.open_project(manifest_path)
    if reopened.jobs.list(project_id):
        raise AssertionError("Opening a project started an implicit job")
    persisted = reopened_context.catalog.list_artifacts(project_id)
    if {item.artifact_id for item in persisted} != {
        item.artifact_id for item in all_artifacts
    }:
        raise AssertionError("Artifact registry changed after cache-only reopen")
    if not reopened_context.catalog.list_selections(project_id):
        raise AssertionError("Saved selections did not persist")

    cache_checks: list[dict[str, str]] = []
    for workflow_id, _ in WORKFLOW_RUNS:
        cached = reopened.workflows.run(
            workflow_id,
            project_id,
            parameters=parameter_map[workflow_id],
        )
        if cached.status != JobStatus.CACHE_HIT:
            cached = reopened.jobs.wait(cached.job_id, timeout=900)
        if cached.status != JobStatus.CACHE_HIT:
            raise AssertionError(
                f"{workflow_id} rebuilt after an unchanged reopen: {cached.status}"
            )
        cache_checks.append(
            {"workflow_id": workflow_id, "status": cached.status.value}
        )

    summary = {
        "schema_version": "2.0",
        "purpose": "Katana Next software/scientific-contract acceptance; not a biological result.",
        "manifest": str(manifest_path),
        "output_root": str(context.manifest.output_root),
        "project_id": project_id,
        "workspace_count": len(
            {workflow.definition.workspace_id for workflow in runtime.workflows.registry.all()}
        ),
        "registered_workflow_count": len(runtime.workflows.registry.all()),
        "executed_fixture_workflow_count": len(WORKFLOW_RUNS),
        "native_sashimi_validation": "Requires a real indexed BAM and is accepted separately.",
        "cosplice_production_mode": cosplice_config.consensus.clustering_method,
        "runs": records,
        "cache_reopen_checks": cache_checks,
        "artifact_count": len(all_artifacts),
        "artifact_kinds": dict(Counter(item.kind.value for item in all_artifacts)),
        "selection_id": selection_id,
        "checks": {
            "no_hidden_workflows": True,
            "all_artifact_files_exist": True,
            "all_run_parameter_files_exist": True,
            "all_artifacts_match_current_snapshot": True,
            "report_archive_paths_are_unique": True,
            "open_is_cache_only": True,
            "unchanged_reruns_are_cache_hits": True,
        },
    }
    destination = target / "acceptance_summary.json"
    destination.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    return destination


def _create_candidate_selection(context) -> str:
    source = next(
        (
            item
            for item in context.catalog.list_artifacts(
                context.manifest.project_id,
                analysis_id="prioritize.evidence",
            )
            if item.name == "Candidate evidence matrix"
        ),
        None,
    )
    if source is None:
        raise RuntimeError("Candidate evidence artifact is missing")
    frame = pd.read_parquet(source.path)
    entity_ids = tuple(frame.head(20)["row_id"].astype(str))
    selection = SelectionSet(
        project_id=context.manifest.project_id,
        name="Acceptance top candidate genes",
        entity_type="gene",
        source_artifact_id=source.artifact_id,
        entity_ids=entity_ids,
        filters={"source": "top 20 transparent evidence rows"},
    )
    return context.catalog.save_selection(selection).selection_id


def _assert_run_files(artifacts) -> None:
    checked_directories: set[Path] = set()
    for artifact in artifacts:
        if not artifact.path.is_file():
            raise FileNotFoundError(artifact.path)
        if "artifacts" not in artifact.path.parts:
            continue
        directory = artifact.path.parent
        if directory in checked_directories:
            continue
        checked_directories.add(directory)
        parameters = directory / "parameters.json"
        if not parameters.is_file():
            raise FileNotFoundError(parameters)
        payload = json.loads(parameters.read_text(encoding="utf-8"))
        if payload.get("workflow_id") != artifact.analysis_id:
            raise AssertionError(f"Invalid parameters.json beside {artifact.path}")


def _assert_snapshot_binding(artifacts, snapshot_fingerprint: str) -> None:
    stale = [
        artifact.artifact_id
        for artifact in artifacts
        if artifact.metadata.get("project_snapshot_fingerprint")
        != snapshot_fingerprint
    ]
    if stale:
        raise AssertionError(
            f"Artifacts are not bound to the current project snapshot: {stale[:5]}"
        )


def _assert_report_archives(artifacts) -> None:
    for artifact in artifacts:
        if artifact.path.suffix.lower() != ".zip":
            continue
        with zipfile.ZipFile(artifact.path) as archive:
            members = archive.namelist()
        if len(members) != len(set(members)):
            raise AssertionError(
                f"Report archive contains duplicate member paths: {artifact.path}"
            )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the isolated Katana Next publication acceptance fixture."
    )
    parser.add_argument("root", help="New or empty external acceptance directory")
    arguments = parser.parse_args()
    print(run_acceptance(arguments.root))


if __name__ == "__main__":
    main()
