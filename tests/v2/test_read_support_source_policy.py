"""Historical bulk QC cannot acquire the current policy through legacy import."""
from __future__ import annotations

import json
from types import SimpleNamespace

import pytest

from katana_next.science.bulk import READ_SUPPORT_POLICY_VERSION
from katana_next.workflows import sources


POLICY_SOURCES = (
    ("quality.stable_as", "Stable-AS event QC"),
    ("differential.as", "Differential AS events"),
    ("landscape.bulk", "Sample-level PSI matrix"),
)


class Catalog:
    def __init__(self, artifacts):
        self.artifacts = artifacts

    def list_artifacts(self, project_id, *, analysis_id, kind):
        assert project_id == "toy-project" and kind == "table"
        return [artifact for artifact in self.artifacts if artifact.analysis_id == analysis_id]

    def get_artifact(self, artifact_id):
        return next((artifact for artifact in self.artifacts if artifact.artifact_id == artifact_id), None)


def artifact(tmp_path, analysis_id, name, artifact_id, metadata):
    path = tmp_path / f"{artifact_id}.parquet"
    path.write_bytes(b"synthetic source identity fixture; this test never reads a table")
    return SimpleNamespace(
        project_id="toy-project", analysis_id=analysis_id, name=name,
        artifact_id=artifact_id, run_fingerprint=artifact_id, path=path,
        metadata={"project_snapshot_fingerprint": "current-snapshot", **metadata},
        contrast_ids=("treated-vs-control",), parameters={},
    )


def context(tmp_path, artifacts):
    return SimpleNamespace(
        manifest=SimpleNamespace(project_id="toy-project", contrasts=[SimpleNamespace(contrast_id="treated-vs-control")]),
        snapshot=SimpleNamespace(manifest_fingerprint="current-snapshot"),
        catalog=Catalog(artifacts),
        results=SimpleNamespace(artifact_directory=lambda analysis_id, fingerprint: tmp_path / "runs" / analysis_id / fingerprint),
    )


def completed_source(tmp_path, analysis_id, name):
    current = artifact(tmp_path, analysis_id, name, "current", {
        "artifact_contract_version": "3.0-correctness",
        "read_support_policy": READ_SUPPORT_POLICY_VERSION,
        "workflow_version": sources._workflow_versions()[analysis_id],
    })
    target = tmp_path / "runs" / analysis_id / current.run_fingerprint
    target.mkdir(parents=True)
    (target / "parameters.json").write_text(json.dumps({
        "completed": True, "generated_artifacts": [{"artifact_id": current.artifact_id}],
    }), encoding="utf-8")
    return current


@pytest.mark.parametrize("analysis_id,name", POLICY_SOURCES)
@pytest.mark.parametrize("contract", [None, "2.0-legacy", "3.0-correctness"])
def test_automatic_resolution_skips_historical_policy_and_keeps_the_current_result(tmp_path, analysis_id, name, contract):
    old_metadata = {"workflow_version": "historical-version"}
    if contract is not None:
        old_metadata["artifact_contract_version"] = contract
    old = artifact(tmp_path, analysis_id, name, "historical", old_metadata)
    original = old.path.read_bytes()
    current = completed_source(tmp_path, analysis_id, name)
    active = context(tmp_path, [old, current])

    assert not sources.source_complete(active, old)
    assert sources.resolve_table(active, analysis_id, name) is current
    assert old.path.read_bytes() == original
    assert active.catalog.get_artifact(old.artifact_id) is old


@pytest.mark.parametrize("analysis_id,name", POLICY_SOURCES)
def test_explicitly_binding_a_legacy_policy_is_rejected_and_keeps_history(tmp_path, analysis_id, name):
    old = artifact(tmp_path, analysis_id, name, "legacy-bound", {})
    active = context(tmp_path, [old])
    parameters = {"source_artifact_ids": {sources.source_key(analysis_id, name): old.artifact_id}}

    assert sources.resolve_table(active, analysis_id, name) is None
    with pytest.raises(ValueError, match="Bound source.*stale"):
        sources.resolve_table(active, analysis_id, name, parameters)
    with pytest.raises(ValueError, match="historical read-support policy"):
        sources.validate_source(active, old)
    assert old.path.is_file()
    assert len(active.catalog.artifacts) == 1


@pytest.mark.parametrize("analysis_id,name", POLICY_SOURCES)
@pytest.mark.parametrize("failure", ["missing_policy", "old_policy", "missing_version", "old_version"])
def test_legacy_contract_cannot_bypass_either_current_policy_or_method_identity(tmp_path, analysis_id, name, failure):
    metadata = {"read_support_policy": READ_SUPPORT_POLICY_VERSION, "workflow_version": sources._workflow_versions()[analysis_id]}
    if failure == "missing_policy":
        metadata.pop("read_support_policy")
    elif failure == "old_policy":
        metadata["read_support_policy"] = "balanced-replicate-consistency"
    elif failure == "missing_version":
        metadata.pop("workflow_version")
    else:
        metadata["workflow_version"] = "historical-version"
    old = artifact(tmp_path, analysis_id, name, "incomplete-identity", metadata)
    assert not sources.source_complete(context(tmp_path, [old]), old)


@pytest.mark.parametrize("analysis_id,name", POLICY_SOURCES)
def test_current_policy_still_requires_a_completed_correctness_contract(tmp_path, analysis_id, name):
    current = completed_source(tmp_path, analysis_id, name)
    active = context(tmp_path, [current])
    bound = {"source_artifact_ids": {sources.source_key(analysis_id, name): current.artifact_id}}
    assert sources.source_complete(active, current)
    assert sources.resolve_table(active, analysis_id, name, bound) is current

    marker = active.results.artifact_directory(analysis_id, current.run_fingerprint) / "parameters.json"
    marker.write_text(json.dumps({"completed": False, "generated_artifacts": [{"artifact_id": current.artifact_id}]}), encoding="utf-8")
    assert not sources.source_complete(active, current)


@pytest.mark.parametrize("analysis_id", ["differential.deg", "differential.dtu", "programs.cosplice", "prioritize.evidence"])
def test_other_historical_module_reading_remains_compatible(tmp_path, analysis_id):
    old = artifact(tmp_path, analysis_id, "Historical source table", "unrelated-history", {"workflow_version": "historical-version"})
    active = context(tmp_path, [old])
    assert sources.source_complete(active, old)
    assert sources.resolve_table(active, analysis_id, old.name) is old
    sources.validate_source(active, old)
