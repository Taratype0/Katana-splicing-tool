import pytest

from katana_next.flow.motif_recovery import prepare_presentation_recovery, RENDER_FAILURE
from katana_next.flow.runtime import canonical_hash, read_json, write_json
from katana_next.reports.motif_native_index import SCIENTIFIC_FILES


def retained(tmp_path):
    task = {"id": "motiftrace-a", "workflow_id": "prioritize.motiftrace_consequences"}
    plan = {"plan_id": "plan-a"}
    root = tmp_path / "failed"
    source = root / "results/artifacts/prioritize.motiftrace_consequences/run/event_consequences"
    source.mkdir(parents=True)
    for name in SCIENTIFIC_FILES:
        (source / name).write_bytes(b"retained-native-data")
    failure = {"schema": "katana.flow.task_receipt/v1", "status": "failed", "error": RENDER_FAILURE,
               "plan_id": plan["plan_id"], "task_id": task["id"], "task_contract_sha256": canonical_hash(task),
               "workflow_id": task["workflow_id"], "source_project_manifest_sha256": "manifest-a",
               "base_request_sha256": "base-a", "finished_at": "2026-10-23T00:00:00Z"}
    failure["receipt_sha256"] = canonical_hash(failure)
    path = write_json(root / "failure.receipt.json", failure)
    args = dict(plan=plan, task=task, base_sha="base-a", manifest_sha="manifest-a")
    return path, source, args


def test_recovery_retains_failed_status_and_binds_every_scientific_file(tmp_path):
    path, source, args = retained(tmp_path)
    before = path.read_bytes()
    recovery = prepare_presentation_recovery(path, **args)
    assert set(recovery.evidence["members"]) == SCIENTIFIC_FILES
    assert len(recovery.inputs) == len(SCIENTIFIC_FILES) + 1
    assert recovery.source == source
    assert path.read_bytes() == before
    assert read_json(path)["status"] == "failed"
    with pytest.raises(ValueError, match="did not execute"):
        recovery.execution_evidence()


@pytest.mark.parametrize("key,value", [
    ("error", "calculation failed"), ("status", "running"), ("plan_id", "other"),
    ("task_id", "other"), ("base_request_sha256", "changed-dependency"),
    ("source_project_manifest_sha256", "other"), ("finished_at", None),
])
def test_recovery_rejects_another_failure_or_request_even_if_signed(tmp_path, key, value):
    path, source, args = retained(tmp_path)
    failure = read_json(path)
    failure[key] = value
    failure["receipt_sha256"] = canonical_hash({k: v for k, v in failure.items() if k != "receipt_sha256"})
    write_json(path, failure)
    with pytest.raises(ValueError, match="exact signed render failure"):
        prepare_presentation_recovery(path, **args)


def test_recovery_rejects_missing_or_completed_family(tmp_path):
    path, source, args = retained(tmp_path)
    (source / "COMPLETE.json").write_text('{}')
    with pytest.raises(ValueError, match="already completed"):
        prepare_presentation_recovery(path, **args)
    # A separate fixture represents an incomplete numerical family.
    path, source, args = retained(tmp_path / "incomplete")
    (source / "MAX_CLUSTER_NULL.npy").unlink()
    with pytest.raises(FileNotFoundError):
        prepare_presentation_recovery(path, **args)


def test_recovery_only_finalizes_after_input_check_and_detects_changed_source(tmp_path, monkeypatch):
    import katana_next.reports.motif_native_index as native
    path, source, args = retained(tmp_path)
    recovery = prepare_presentation_recovery(path, **args)
    called = []
    def check(*a, **kw):
        called.append("validate")
    def finalize(*a, **kw):
        called.append("present")
        return {"native_statistics_modified": False}
    monkeypatch.setattr(native, "validate_recorded_consequences_inputs", check)
    monkeypatch.setattr(native, "finalize_recorded_consequences", finalize)
    assert recovery.run(None, [], {}, "reference", tmp_path / "new") == {"native_statistics_modified": False}
    assert called == ["validate", "present"]
    assert recovery.execution_evidence()["scientific_computation_performed"] is False
    with pytest.raises(ValueError, match="only one"):
        recovery.run(None, [], {}, "reference", tmp_path / "other")
    second = prepare_presentation_recovery(path, **args)
    (source / "RESULTS.json").write_bytes(b"changed-after-verification")
    with pytest.raises(ValueError, match="Verified content changed"):
        second.run(None, [], {}, "reference", tmp_path / "changed")
