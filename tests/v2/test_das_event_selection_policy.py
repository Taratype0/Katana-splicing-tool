"""DAS membership is a reproducible downstream choice, not a caller refit."""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from katana_next.science.bulk import READ_SUPPORT_POLICY_VERSION
from katana_next.science.das_selection import (
    legacy_kernel_event_view, normalize_event_set, normalize_workflow_selection,
    policy_matches_artifact, selected_event_mask, selected_event_view,
    selection_metadata, workflow_das_capability,
)
from katana_next.workflows.bulk import CandidateEvidenceWorkflow, CrossComparisonWorkflow, DifferentialAsWorkflow
from katana_next.workflows import sources


def _events():
    return pd.DataFrame([
        dict(row_id=f"A::{i}", comparison_id="A", event_id=f"e{i}", gene_id=f"g{i}", gene_symbol=f"G{i}",
             event_type="SE", significant=sig, stable_significant=rs,
             measurement_eligible=read, read_support_policy=READ_SUPPORT_POLICY_VERSION,
             criterion_min_reads=10, dpsi=.2, abs_dpsi=.2, fdr=.001, pvalue=.0001)
        for i, sig, rs, read in [(1, True, True, True), (2, True, False, False), (3, False, False, True)]
    ])


@pytest.mark.parametrize("value", [None, "stable_significant", "RS-DAS", "rs_das"])
def test_rs_default_and_aliases(value):
    assert normalize_event_set(value) == "stable_significant"


@pytest.mark.parametrize("value", ["significant", "DAS", "ordinary"])
def test_ordinary_aliases(value):
    assert normalize_event_set(value) == "significant"


def test_selection_changes_only_membership_and_keeps_full_source():
    events = _events()
    before = events.copy(deep=True)
    rs = selected_event_view(events)
    ordinary = selected_event_view(events, "significant")
    assert rs.selected_das.tolist() == [True, False, False]
    assert ordinary.selected_das.tolist() == [True, True, False]
    assert len(rs) == len(ordinary) == len(events)
    pd.testing.assert_frame_equal(events, before)
    for column in before:
        pd.testing.assert_series_equal(ordinary[column], before[column])
        pd.testing.assert_series_equal(rs[column], before[column])
    adapted = legacy_kernel_event_view(events, "significant")
    assert adapted.stable_significant.tolist() == [True, True, False]
    assert adapted.source_stable_significant.tolist() == [True, False, False]
    pd.testing.assert_frame_equal(events, before)


@pytest.mark.parametrize("change", ["missing_flag", "missing_audit", "historical", "wrong_counts", "false_support", "not_subset", "string_flag"])
def test_rs_missing_or_inconsistent_evidence_never_falls_back(change):
    events = _events()
    if change == "missing_flag":
        events = events.drop(columns="stable_significant")
    elif change == "missing_audit":
        events = events.drop(columns="measurement_eligible")
    elif change == "historical":
        events["read_support_policy"] = "legacy_consistency"
    elif change == "wrong_counts":
        events["criterion_min_reads"] = 5
    elif change == "false_support":
        events.loc[0, "measurement_eligible"] = False
    elif change == "not_subset":
        events.loc[0, "significant"] = False
    elif change == "string_flag":
        events["stable_significant"] = ["True", "False", "False"]
    with pytest.raises(ValueError):
        selected_event_mask(events)
    if change not in {"not_subset"}:
        assert selected_event_mask(events, "significant").tolist() == [True, True, False]


def test_definitions_and_normalization_have_one_default_and_declared_boundaries():
    for workflow in [DifferentialAsWorkflow(), CrossComparisonWorkflow(), CandidateEvidenceWorkflow()]:
        assert workflow.normalize_parameters({})["event_set"] == "stable_significant"
        definition = workflow.definition.to_dict()
        assert definition["das_selection"]["mode"] == "switchable"
        assert len([p for p in definition["parameters"] if p["key"] == "event_set"]) == 1
        assert workflow.normalize_parameters({"event_set": "DAS"})["event_set"] == "significant"
    assert workflow_das_capability("programs.cosplice")["mode"] == "measurement_population"
    with pytest.raises(ValueError, match="does not apply"):
        normalize_workflow_selection("programs.cosplice", {"event_set": "significant"})
    with pytest.raises(ValueError, match="requires RS-DAS"):
        normalize_workflow_selection("prioritize.rbp_mechanism", {"event_set": "significant"})
    assert normalize_workflow_selection("prioritize.event_motifs", {})["event_set"] == "stable_significant"


def _source(tmp_path, identity, event_set, analysis_id="prioritize.evidence", name="Candidate evidence matrix"):
    path = tmp_path / f"{identity}.parquet"
    path.write_bytes(b"unchanged source identity fixture")
    metadata = {"project_snapshot_fingerprint": "snapshot"}
    if event_set is not None:
        metadata["das_selection"] = selection_metadata(event_set)
    return SimpleNamespace(project_id="p", analysis_id=analysis_id, name=name, artifact_id=identity,
        run_fingerprint=identity, path=path, metadata=metadata, contrast_ids=("A",), parameters={})


def _context(tmp_path, artifacts):
    return SimpleNamespace(
        manifest=SimpleNamespace(project_id="p", contrasts=[SimpleNamespace(contrast_id="A")]),
        snapshot=SimpleNamespace(manifest_fingerprint="snapshot"),
        catalog=SimpleNamespace(list_artifacts=lambda *args, **kwargs: artifacts,
            get_artifact=lambda key: next((a for a in artifacts if a.artifact_id == key), None)),
        results=SimpleNamespace(artifact_directory=lambda *args: tmp_path / "none"))


def test_derived_source_resolution_matches_policy_and_preserves_history(tmp_path):
    ordinary = _source(tmp_path, "ordinary", "significant")
    historical = _source(tmp_path, "unmarked", None)
    rs = _source(tmp_path, "rs", "stable_significant")
    context = _context(tmp_path, [ordinary, historical, rs])
    before = {a.artifact_id: a.path.read_bytes() for a in [ordinary, historical, rs]}
    assert sources.resolve_table(context, "prioritize.evidence", "Candidate evidence matrix", {"event_set": "stable_significant"}) is rs
    assert sources.resolve_table(context, "prioritize.evidence", "Candidate evidence matrix", {"event_set": "significant"}) is ordinary
    assert not policy_matches_artifact(historical, "stable_significant")
    with pytest.raises(ValueError, match="different or unrecorded"):
        sources.resolve_table(context, "prioritize.evidence", "Candidate evidence matrix", {
            "event_set": "stable_significant", "source_artifact_ids": {"prioritize.evidence::Candidate evidence matrix": "ordinary"}})
    assert all(a.path.read_bytes() == before[a.artifact_id] for a in [ordinary, historical, rs])


def test_only_complete_das_event_table_is_dual_policy_not_its_rs_plot(tmp_path):
    raw = _source(tmp_path, "raw", "stable_significant", "differential.as", "Differential AS events")
    plot = _source(tmp_path, "plot", "stable_significant", "differential.as", "Differential AS volcano")
    assert policy_matches_artifact(raw, "significant")
    assert not policy_matches_artifact(plot, "significant")


def test_service_rs_required_selection_checks_selected_rows_without_dropping_them(tmp_path, monkeypatch):
    from katana_next.services.workflow_service import WorkflowService
    source = _source(tmp_path, "raw", "stable_significant", "differential.as", "Differential AS events")
    _events().to_parquet(source.path, index=False)
    before = source.path.read_bytes()
    context = _context(tmp_path, [source])
    parameters = {"event_set": "stable_significant", "selection_id": "explicit-choice"}
    monkeypatch.setattr(sources, "selected_source", lambda *_: source)
    monkeypatch.setattr(sources, "selection_rows", lambda *args: args[-1].iloc[[1]])
    with pytest.raises(ValueError, match="requires an explicit RS-DAS selection"):
        WorkflowService._validate_das_parent_selection(context, "prioritize.event_motifs", parameters, [source])
    monkeypatch.setattr(sources, "selection_rows", lambda *args: args[-1].iloc[[0]])
    WorkflowService._validate_das_parent_selection(context, "prioritize.event_motifs", parameters, [source])
    assert source.path.read_bytes() == before


def test_real_workflow_fixture_selection_propagation_and_unchanged_raw(tmp_path, monkeypatch):
    from katana_next.api.app import PlatformRuntime
    from katana_next.domain.models import JobStatus
    from katana_next.full_fixture import create_full_fixture
    from katana_next.infrastructure.manifest_store import ManifestStore
    from katana_next.infrastructure.fingerprints import fingerprint_file
    import katana_next.workflows.figure_exports as exports

    # Figure rendering is independently tested by report QA; this test exercises
    # native artifact registration, completion, dependency selection and cache.
    monkeypatch.setattr(exports, "export_vega", lambda *args, **kwargs: {})
    manifest_path = create_full_fixture(tmp_path / "fixture")
    _, manifest = ManifestStore().load(manifest_path)
    target = next(a for a in manifest.assets if a.kind == "rmats_result" and a.metadata["event_type"] == "SE")
    raw = pd.read_csv(target.path, sep="\t")
    n1 = len(target.metadata["numerator_sample_ids"])
    n2 = len(target.metadata["denominator_sample_ids"])
    raw.loc[0, ["FDR", "PValue", "IncLevelDifference", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "SJC_SAMPLE_1"]] = [
        .001, .0001, .2, ",".join([".7"]*n1), ",".join([".5"]*n2), ",".join(["1"]*n1), ",".join(["1"]*n1)]
    raw.to_csv(target.path, sep="\t", index=False)
    hashes = {a.asset_id: fingerprint_file(a.path, force_full=True).digest for a in manifest.assets if a.kind == "rmats_result"}
    runtime = PlatformRuntime(application_root=tmp_path / "app")
    try:
        context = runtime.projects.open_project(manifest_path)
        project_id = context.manifest.project_id
        def run(workflow, event_set=None):
            parameters = {} if event_set is None else {"event_set": event_set}
            job = runtime.workflows.run(workflow, project_id, parameters=parameters)
            finished = runtime.jobs.wait(job.job_id, timeout=60)
            assert finished.status == JobStatus.COMPLETE, finished.error
            return context.catalog.list_artifacts(project_id, analysis_id=workflow, kind="table")
        das_artifacts = run("differential.as")
        das = next(a for a in das_artifacts if a.name == "Differential AS events")
        raw_events = pd.read_parquet(das.path)
        source_bytes = das.path.read_bytes()
        assert (raw_events.significant & ~raw_events.stable_significant).any()
        rs_artifact = next(a for a in run("prioritize.evidence") if a.name == "Candidate evidence matrix")
        rs = pd.read_parquet(rs_artifact.path)
        ordinary_artifact = next(a for a in run("prioritize.evidence", "significant") if a.parameters["event_set"] == "significant" and a.name == "Candidate evidence matrix")
        ordinary = pd.read_parquet(ordinary_artifact.path)
        assert ordinary.selected_significant_as_events.sum() > rs.selected_significant_as_events.sum()
        assert ordinary.stable_significant_as_events.sum() == rs.stable_significant_as_events.sum()
        assert das.path.read_bytes() == source_bytes
        assert ordinary_artifact.metadata["das_selection"]["event_set"] == "significant"
        assert rs_artifact.metadata["das_selection"]["event_set"] == "stable_significant"
        enrichment_artifacts = run("prioritize.enrichment", "significant")
        enrichment = next(a for a in enrichment_artifacts if a.name == "Gene-set enrichment")
        assert ordinary_artifact.artifact_id in enrichment.parent_artifact_ids
        assert rs_artifact.artifact_id not in enrichment.parent_artifact_ids
        genes_artifact = next(a for a in enrichment_artifacts if a.name == "Enrichment gene inputs by comparison")
        inputs = pd.read_parquet(genes_artifact.path)
        expected_foreground = set(zip(ordinary.loc[ordinary.as_supported, "comparison_id"], ordinary.loc[ordinary.as_supported, "gene_id"]))
        observed_foreground = set(zip(inputs.loc[inputs.selected_input, "comparison_id"], inputs.loc[inputs.selected_input, "gene_id"]))
        assert observed_foreground == expected_foreground
        assert ordinary.as_eligible.sum() >= ordinary.as_read_supported_eligible.sum()
        assert das.path.read_bytes() == source_bytes
        assert all(fingerprint_file(a.path, force_full=True).digest == hashes[a.asset_id] for a in manifest.assets if a.kind == "rmats_result")
    finally:
        runtime.jobs.shutdown()
