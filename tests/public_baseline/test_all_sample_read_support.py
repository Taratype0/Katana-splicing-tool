"""Real application workflows for the approved every-sample ten-count policy.

All inputs are generated in pytest's temporary project.  No frozen application,
experimental data, native statistical caller, or external runtime is required.
The native P/FDR values below are imported fixture evidence, never recomputed.
"""
from __future__ import annotations

import hashlib
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from katana_next.api.app import PlatformRuntime
from katana_next.domain.models import JobStatus
from katana_next.full_fixture import create_full_fixture
from katana_next.infrastructure.manifest_store import ManifestStore
from katana_next.science.bulk import READ_SUPPORT_POLICY_VERSION, load_bulk_events


GOOD = "perturbation-vs-baseline"
BAD = "rescue-vs-baseline"
SELECTED = [GOOD, BAD]
NATIVE_COLUMNS = [
    "PValue", "FDR", "IncLevelDifference", "IncLevel1", "IncLevel2",
    "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2",
]


def _input_hashes(project_root: Path, manifest_path: Path) -> dict[str, str]:
    paths = sorted((project_root / "inputs").rglob("*")) + [manifest_path]
    return {
        str(path.relative_to(project_root)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in paths if path.is_file()
    }


@pytest.fixture(scope="module")
def workflow_case(tmp_path_factory):
    root = tmp_path_factory.mktemp("all-sample-workflows")
    project_root = root / "fixture"
    manifest_path = create_full_fixture(project_root)
    _path, manifest = ManifestStore().load(manifest_path)

    # These are genuine rMATS-shaped files consumed by the normal public loader.
    # The first event's shared Baseline sample has 30 reads in GOOD and 9 in BAD.
    # A median-only matrix would conceal that failure as 19.5 reads.
    raw_tables = {}
    for asset in manifest.assets:
        if asset.kind != "rmats_result":
            continue
        frame = pd.read_csv(asset.path, sep="\t")
        frame["PValue"] = np.arange(1, len(frame) + 1, dtype=float) / 100_000
        if asset.metadata["event_type"] == "SE":
            for arm in (1, 2):
                frame.loc[0, f"IJC_SAMPLE_{arm}"] = "15,15,15,15,15"
                frame.loc[0, f"SJC_SAMPLE_{arm}"] = "15,15,15,15,15"
                frame.loc[1, f"IJC_SAMPLE_{arm}"] = "5,5,5,5,5"
                frame.loc[1, f"SJC_SAMPLE_{arm}"] = "5,5,5,5,5"
                frame.loc[2, f"IJC_SAMPLE_{arm}"] = "0,1,5,9,10"
                frame.loc[2, f"SJC_SAMPLE_{arm}"] = "10,9,5,1,0"
                frame.loc[3, f"IJC_SAMPLE_{arm}"] = "5,5,5,5,5"
                frame.loc[3, f"SJC_SAMPLE_{arm}"] = "5,5,5,5,5"
            frame.loc[1, "IncLevel1"] = "0.05,0.05,0.95,0.95,0.05"
            frame.loc[1, "IncLevel2"] = "0.2,0.2,0.2,0.2,0.2"
            frame.loc[1, "IncLevelDifference"] = 0.21
            frame.loc[1, "FDR"] = 0.0123
            frame.loc[1, "PValue"] = 0.004321
            frame.loc[2, "IncLevel1"] = "0.2,0.25,0.3,0.35,0.4"
            frame.loc[2, "IncLevel2"] = "0.1,0.1,0.1,0.1,0.1"
            frame.loc[2, "IncLevelDifference"] = 0.2
            frame.loc[3, "IncLevel1"] = "0.2,0.3,NA,0.4,0.5"
            if asset.metadata["contrast_id"] == BAD:
                frame.loc[0, "IJC_SAMPLE_2"] = "15,15,4,15,15"
                frame.loc[0, "SJC_SAMPLE_2"] = "15,15,5,15,15"
                # This event has good measurements in GOOD but no Rescue values
                # at all.  A full selected-sample matrix must also reject it.
                frame = frame.drop(index=4)
        frame.to_csv(asset.path, sep="\t", index=False)
        raw_tables[str(asset.path)] = pd.read_csv(asset.path, sep="\t")

    before = _input_hashes(project_root, manifest_path)
    runtime = PlatformRuntime(application_root=root / "application")
    try:
        context = runtime.projects.open_project(manifest_path)

        def run(workflow_id, parameters):
            started = runtime.workflows.run(
                workflow_id, context.manifest.project_id, parameters=parameters
            )
            finished = runtime.jobs.wait(started.job_id, timeout=90)
            assert finished.status in {JobStatus.COMPLETE, JobStatus.CACHE_HIT}, (
                workflow_id, finished.status, finished.error
            )
            artifacts = [context.catalog.get_artifact(key) for key in finished.artifact_ids]
            assert all(item is not None and item.path.is_file() for item in artifacts)
            return {item.name: item for item in artifacts if item.kind == "table"}

        quality = run("quality.stable_as", {
            "contrast_ids": SELECTED, "fdr": 0.02, "min_abs_dpsi": 0.1,
        })["Stable-AS event QC"]
        inherited = {"contrast_ids": SELECTED, "quality_artifact_id": quality.artifact_id}
        das = run("differential.as", inherited)["Differential AS events"]
        landscape = run("landscape.bulk", inherited)
        selected_landscape = run("landscape.bulk", {
            "contrast_ids": [GOOD], "quality_artifact_id": quality.artifact_id,
        })
        deg_default = run("differential.deg", {"contrast_ids": [GOOD]})[
            "Differentially expressed genes"
        ]
        deg_custom = run("differential.deg", {
            "contrast_ids": [GOOD], "fdr": 0.01, "min_abs_log2fc": 1.6,
        })["Differentially expressed genes"]
        assert _input_hashes(project_root, manifest_path) == before
        yield SimpleNamespace(
            context=context, runtime=runtime, quality=quality, das=das,
            landscape=landscape, selected_landscape=selected_landscape,
            deg_default=deg_default, deg_custom=deg_custom,
            raw_tables=raw_tables, input_hashes=before, project_root=project_root,
            manifest_path=manifest_path,
        )
    finally:
        runtime.jobs.shutdown()
        assert _input_hashes(project_root, manifest_path) == before


def _se_row(frame: pd.DataFrame, comparison: str, native_id: int) -> pd.Series:
    rows = frame.loc[
        frame["comparison_id"].eq(comparison)
        & frame["event_type"].eq("SE")
        & frame["ID"].astype(int).eq(native_id)
    ]
    assert len(rows) == 1
    return rows.iloc[0]


def _native_id(case, position: int) -> int:
    table = next(frame for path, frame in case.raw_tables.items()
                 if "Perturbation-Baseline" in path and Path(path).name.startswith("SE."))
    return int(table.iloc[position]["ID"])


def test_quality_to_das_inherits_every_sample_policy_without_altering_native_evidence(workflow_case):
    case = workflow_case
    quality = pd.read_parquet(case.quality.path)
    das = pd.read_parquet(case.das.path)
    assert set(quality["comparison_id"]) == set(SELECTED)
    assert case.das.metadata["measurement_policy_source"] == case.quality.artifact_id
    assert case.das.metadata["read_support_policy"] == READ_SUPPORT_POLICY_VERSION
    assert case.das.metadata["read_support_min_counts"] == 10
    assert quality["read_support_policy"].eq(READ_SUPPORT_POLICY_VERSION).all()
    assert quality["criterion_min_reads"].eq(10).all()
    pd.testing.assert_frame_equal(
        quality.set_index("row_id").sort_index(),
        das.set_index("row_id").sort_index().drop(columns="minus_log10_fdr"),
    )
    for source, observed in quality.groupby("source_path"):
        source_rows = observed["source_row"].astype(int).tolist()
        expected = case.raw_tables[source].iloc[source_rows][NATIVE_COLUMNS].reset_index(drop=True)
        pd.testing.assert_frame_equal(
            observed[NATIVE_COLUMNS].reset_index(drop=True), expected,
        )
    assert _input_hashes(case.project_root, case.manifest_path) == case.input_hashes


def test_one_low_third_sample_fails_but_direction_and_se_are_diagnostic_only(workflow_case):
    case = workflow_case
    for artifact in (case.quality, case.das):
        frame = pd.read_parquet(artifact.path)
        good = _se_row(frame, GOOD, _native_id(case, 0))
        low = _se_row(frame, BAD, _native_id(case, 0))
        assert good["measurement_eligible"]
        assert low["denominator_low_read_samples"] == 1
        assert low["denominator_expected_samples"] == 5
        assert low["denominator_read_count_supported_samples"] == 4
        assert low["passes_data_usability"]
        assert not low["passes_read_support"] and not low["measurement_eligible"]
        assert not low["stable_event"]
        diagnostic = _se_row(frame, GOOD, _native_id(case, 1))
        assert diagnostic["replicate_pairwise_direction_consistency"] < 0.75
        assert diagnostic["max_measurement_se"] > 0.12
        assert not diagnostic["measurement_se_gate_enabled"]
        assert diagnostic["measurement_eligible"] and diagnostic["stable_significant"]
        boundary = _se_row(frame, GOOD, _native_id(case, 2))
        assert boundary["measurement_eligible"]
        unusable = _se_row(frame, GOOD, _native_id(case, 3))
        assert unusable["passes_read_counts"]
        assert not unusable["passes_data_usability"]
        assert not unusable["measurement_eligible"]


def test_landscape_requires_all_contrasts_and_complete_selected_sample_matrix(workflow_case):
    case = workflow_case
    quality = pd.read_parquet(case.quality.path)
    event = _se_row(quality, GOOD, _native_id(case, 0))["event_id"]
    absent_arm_event = _se_row(quality, GOOD, _native_id(case, 4))["event_id"]
    # Verify that this fixture defeats the old any-passing-contrast rule even
    # when the deduplicated read matrix masks the bad replicate with its median.
    loaded = load_bulk_events(case.context, {"contrast_ids": SELECTED})
    assert loaded.read_matrix.loc[event, "baseline-3"] == 19.5
    assert loaded.read_matrix.loc[event].ge(10).all()
    rows = quality[quality["event_id"].eq(event)]
    assert rows["measurement_eligible"].tolist() == [True, False]
    matrix_artifact = case.landscape["Sample-level PSI matrix"]
    assert matrix_artifact.metadata["measurement_policy_source"] == case.quality.artifact_id
    matrix = pd.read_parquet(matrix_artifact.path).set_index("event_id")
    assert event not in matrix.index
    assert absent_arm_event not in matrix.index
    selected_samples = {sample.sample_id for sample in case.context.manifest.samples}
    assert set(matrix.columns) - {"row_id"} == selected_samples
    assert np.isfinite(matrix.drop(columns="row_id").to_numpy()).all()
    assert matrix.drop(columns="row_id").ge(0).all().all()
    assert matrix.drop(columns="row_id").le(1).all().all()
    expected = loaded.psi_matrix.drop(index=[event, absent_arm_event], errors="ignore")
    # Other ineligible rows are excluded as well; all retained values remain
    # the loader's actual finite medians, not newly imputed sample values.
    pd.testing.assert_frame_equal(
        matrix.drop(columns="row_id").sort_index(axis=1),
        expected.loc[matrix.index].sort_index(axis=1),
    )
    pca = pd.read_parquet(case.landscape["Sample PSI PCA"].path)
    assert set(pca["sample_id"]) == selected_samples
    assert pca["events_used"].eq(len(matrix)).all()
    similarity = pd.read_parquet(case.landscape["Sample PSI similarity"].path)
    assert similarity["shared_events"].eq(len(matrix)).all()


def test_unselected_failed_contrast_does_not_exclude_the_selected_arms(workflow_case):
    case = workflow_case
    quality = pd.read_parquet(case.quality.path)
    matrix = pd.read_parquet(case.selected_landscape["Sample-level PSI matrix"].path)
    assert _se_row(quality, GOOD, _native_id(case, 0))["event_id"] in set(matrix["event_id"])
    assert _se_row(quality, GOOD, _native_id(case, 4))["event_id"] in set(matrix["event_id"])
    contrast = next(item for item in case.context.manifest.contrasts if item.contrast_id == GOOD)
    selected_samples = set(contrast.numerator.sample_ids + contrast.denominator.sample_ids)
    assert set(matrix.columns) - {"event_id", "row_id"} == selected_samples
    assert not any(column.startswith("rescue-") for column in matrix.columns)
    pca = pd.read_parquet(case.selected_landscape["Sample PSI PCA"].path)
    assert set(pca["sample_id"]) == selected_samples


def test_retired_qc_controls_are_absent_and_deg_remains_user_configurable(workflow_case):
    case = workflow_case
    retired = {
        "min_reads", "stability_profile", "min_supported_replicates",
        "min_direction_consistency", "max_measurement_se", "apply_measurement_se_gate",
    }
    for workflow_id in ("quality.stable_as", "differential.as", "landscape.bulk"):
        definition = case.runtime.workflows.registry.get(workflow_id).definition
        assert not retired.intersection(parameter.key for parameter in definition.parameters)
    deg = case.runtime.workflows.registry.get("differential.deg").definition
    assert {parameter.key for parameter in deg.parameters} == {"fdr", "min_abs_log2fc"}
    default = pd.read_parquet(case.deg_default.path).set_index("row_id").sort_index()
    custom = pd.read_parquet(case.deg_custom.path).set_index("row_id").sort_index()
    pd.testing.assert_series_equal(default["deg_fdr"], custom["deg_fdr"])
    pd.testing.assert_series_equal(default["log2fc"], custom["log2fc"])
    pd.testing.assert_series_equal(
        custom["significant"],
        (custom["deg_fdr"].le(0.01) & custom["log2fc"].abs().ge(1.6)).rename("significant"),
    )
    assert custom["significant"].sum() > 0
    assert custom["significant"].sum() < default["significant"].sum()
    assert case.deg_custom.parameters["fdr"] == 0.01
    assert case.deg_custom.parameters["min_abs_log2fc"] == 1.6
