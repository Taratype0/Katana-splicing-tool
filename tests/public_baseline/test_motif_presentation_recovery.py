"""Generated calculation + explicit simulated presentation failure contracts.

These validate recovery mechanics; native failure lifecycle is separately
established by the full generated Nextflow acceptance, never by this fixture.
"""
import copy
import json
from pathlib import Path
import shutil

import pandas as pd
import pytest

from katana_next.flow.model import canonical_hash, file_sha256, read_json, write_json
from katana_next.reports.motif_native_index import (
    SCIENTIFIC_FILES, OVERSIZE_ERROR, finalize_recorded_consequences, validate_recorded_consequences_inputs,
)
from katana_next.science.motiftrace_execution import analyze_consequences
from test_motiftrace_execution_budget import _generated


@pytest.fixture(scope="module")
def recorded(tmp_path_factory):
    root = tmp_path_factory.mktemp("generated-render-failure")
    complete, _ids, models, fasta = _generated(root)
    complete = complete.iloc[:4].copy()
    ids = set(complete.loc[complete.stable_significant, "event_id"])
    models = {"consensus": models["consensus"]}
    provenance = {"genome_FASTA": str(fasta), "source_artifact_id": "generated-source",
        "implementation_SHA256": {"motiftrace_clusters.py": "frozen-fixture-identity", "motiftrace_execution.py": "old-adapter"}}
    source = root / "native-completed"
    analyze_consequences(complete, ids, models, fasta, source, permutations=99, provenance=provenance)
    partial = root / "simulated-failed-display"
    partial.mkdir()
    for name in SCIENTIFIC_FILES:
        shutil.copyfile(source / name, partial / name)
    failure = {"schema": "katana.flow.task_receipt/v1", "status": "failed", "task_id": "generated-display-failure",
        "task_contract_sha256": "a"*64, "workflow_id": "prioritize.motiftrace_consequences", "error": OVERSIZE_ERROR,
        "synthetic_lifecycle_fixture": True, "selection": {"contrast_ids": ["B_vs_A"]},
        "inputs": [{"path": str(p), "sha256": file_sha256(p)} for p in (fasta, Path(str(fasta)+".fai"))]}
    failure["receipt_sha256"] = canonical_hash(failure)
    failure_path = write_json(root / "failure.receipt.json", failure)
    evidence = {"failure_receipt": str(failure_path), "failure_receipt_sha256": file_sha256(failure_path),
        "task_id": failure["task_id"], "task_contract_sha256": failure["task_contract_sha256"],
        "scientific_plan_sha256": file_sha256(partial / "ANALYSIS_PLAN_BEFORE_SCORING.json"),
        "members": {name: {"bytes": (partial/name).stat().st_size, "sha256": file_sha256(partial/name)} for name in SCIENTIFIC_FILES}}
    return partial, evidence, complete, ids, models, fasta, provenance


def test_real_calculation_only_rebuilds_presentation_and_preserves_failed_source(recorded, tmp_path):
    source, evidence, complete, ids, models, fasta, provenance = recorded
    current = copy.deepcopy(provenance)
    current["implementation_SHA256"].update({"motiftrace_execution.py": "new-adapter", "motif_native_index.py": "new-view"})
    current["flow_execution_adapter_sha256"] = "new-flow-adapter"
    proof = validate_recorded_consequences_inputs(source, complete, ids, models, fasta,
        failure_receipt=evidence["failure_receipt"], provenance=current, permutations=99)
    assert proof["scientific_calculations_rerun"] == 0
    assert proof["execution_interface_changes"]["flow_execution_adapter_sha256"]["after"] == "new-flow-adapter"
    destination = tmp_path / "recovered"
    summary = finalize_recorded_consequences(source, destination, evidence)
    assert summary == read_json(source / "RESULTS.json")
    assert read_json(evidence["failure_receipt"])["status"] == "failed"
    assert not (source / "COMPLETE.json").exists()
    assert read_json(destination / "COMPLETE.json")["presentation_recovery"] is True
    for name, identity in evidence["members"].items():
        assert file_sha256(source / name) == file_sha256(destination / name) == identity["sha256"]
    assert (destination / "MOTIFTRACE_EVENT_CONSEQUENCES.html").stat().st_size < 20000
    with pytest.raises(FileExistsError):
        finalize_recorded_consequences(source, destination, evidence)


@pytest.mark.parametrize("kind", ["table", "model", "selection", "permutation", "budget", "scientific_source", "source_provenance"])
def test_recovery_rejects_changed_scientific_request(recorded, kind):
    source, evidence, complete, ids, models, fasta, provenance = recorded
    complete, ids, models, provenance = complete.copy(), set(ids), copy.deepcopy(models), copy.deepcopy(provenance)
    kwargs = {"permutations": 99}
    if kind == "table": complete.loc[complete.index[0], "dpsi"] = .7
    elif kind == "model": models["consensus"]["pattern"] = "AAA"
    elif kind == "selection": ids.pop()
    elif kind == "permutation": kwargs["permutations"] = 199
    elif kind == "budget": kwargs["max_tests"] = 2_000_000
    elif kind == "scientific_source": provenance["implementation_SHA256"]["motiftrace_clusters.py"] = "changed"
    elif kind == "source_provenance": provenance["source_artifact_id"] = "other"
    with pytest.raises(ValueError):
        validate_recorded_consequences_inputs(source, complete, ids, models, fasta,
            failure_receipt=evidence["failure_receipt"], provenance=provenance, **kwargs)


@pytest.mark.parametrize("kind", ["missing_member", "changed_member", "wrong_failure", "wrong_plan"])
def test_recovery_does_not_certify_partial_or_wrong_failed_attempt(recorded, tmp_path, kind):
    source, original, *_ = recorded
    evidence = copy.deepcopy(original)
    if kind == "missing_member": evidence["members"].pop("MAX_CLUSTER_NULL.npy")
    elif kind == "changed_member": evidence["members"]["NEW_WINDOW_AND_FULL_SEQUENCE_TESTS.tsv.gz"]["sha256"] = "0"*64
    elif kind == "wrong_failure": evidence["task_contract_sha256"] = "f"*64
    elif kind == "wrong_plan": evidence["scientific_plan_sha256"] = "0"*64
    with pytest.raises(ValueError):
        finalize_recorded_consequences(source, tmp_path / "must-not-complete", evidence)
    assert not (tmp_path / "must-not-complete/COMPLETE.json").exists()
