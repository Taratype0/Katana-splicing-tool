"""Deletion tests operate only on new generated files in pytest temporary roots."""
import json
import os
from pathlib import Path

import pytest

from katana_next.flow import storage


SUCCESS = {"state": "COMPLETED", "exit_code": 0, "engine": "nextflow"}


def _scope(tmp_path):
    source = tmp_path / "inputs"
    source.mkdir()
    (source / "reads.fastq").write_text("original generated input\n")
    results = tmp_path / "results"
    results.mkdir()
    (results / "REPORT.html").write_text("retained generated result\n")
    scope = tmp_path / "managed"
    storage.create_managed_scope(scope, protected_roots=(source, results), run_id="generated-test-run")
    return scope, source, results


def _checkpoint(work):
    assert work.is_dir()
    return {"status": "completed", "validated": "independent generated checkpoint"}


def test_cleanup_only_registered_stage_keeps_other_stage_inputs_and_final_results(tmp_path):
    scope, source, results = _scope(tmp_path)
    with storage.managed_run_lock(scope):
        first = storage.register_stage(scope, "upstream")
        work = Path(first["work_root"])
        (work / "native").mkdir()
        (work / "native" / "intermediate.bin").write_bytes(b"generated cache" * 100)
        second = Path(storage.register_stage(scope, "downstream")["work_root"])
        (second / "running-task").write_text("do not remove")
        receipt = storage.cleanup_stage(scope, "upstream", native_status=SUCCESS,
                                        checkpoint_validator=_checkpoint)
        assert receipt["status"] == "completed"
        assert receipt["deleted_entries"] == 3
        assert receipt["accounting"]["logical_file_bytes"] == 1500
        assert receipt["accounting"]["quota_reclaimed_bytes"] is None
        assert "resume" in receipt["resume_lost"]
        assert not work.exists()
        assert (second / "running-task").read_text() == "do not remove"
    assert (source / "reads.fastq").read_text() == "original generated input\n"
    assert (results / "REPORT.html").read_text() == "retained generated result\n"
    assert not (scope / storage._LOCK).exists()


@pytest.mark.parametrize("native", [
    {"state": "RUNNING", "exit_code": 0}, {"state": "FAILED", "exit_code": 1},
    {"state": "COMPLETED", "exit_code": 1}, {"state": "COMPLETED"},
    {"state": "COMPLETED", "exit_code": False}, "COMPLETED", {},
])
def test_running_failed_or_unproven_native_stage_is_never_deleted(tmp_path, native):
    scope, _source, _results = _scope(tmp_path)
    with storage.managed_run_lock(scope):
        work = Path(storage.register_stage(scope, "upstream")["work_root"])
        (work / "keep").write_text("unchanged")
        with pytest.raises(ValueError, match="completed native"):
            storage.cleanup_stage(scope, "upstream", native_status=native,
                                  checkpoint_validator=lambda _: pytest.fail("must not call validator"))
        assert (work / "keep").read_text() == "unchanged"


@pytest.mark.parametrize("checkpoint", [None, True, {"status": "failed"}, {"status": "running"}])
def test_missing_or_failed_checkpoint_preserves_cache(tmp_path, checkpoint):
    scope, _source, _results = _scope(tmp_path)
    with storage.managed_run_lock(scope):
        work = Path(storage.register_stage(scope, "upstream")["work_root"])
        (work / "keep").write_text("unchanged")
        with pytest.raises(ValueError, match="checkpoint"):
            storage.cleanup_stage(scope, "upstream", native_status=SUCCESS,
                                  checkpoint_validator=lambda _: checkpoint)
        assert (work / "keep").read_text() == "unchanged"


def test_dry_run_and_same_registration_are_non_mutating(tmp_path):
    scope, _source, _results = _scope(tmp_path)
    with storage.managed_run_lock(scope):
        first = storage.register_stage(scope, "upstream")
        work = Path(first["work_root"])
        (work / "keep").write_bytes(b"cache")
        assert storage.register_stage(scope, "upstream") == first
        planned = storage.cleanup_stage(scope, "upstream", native_status=SUCCESS,
            checkpoint_validator=_checkpoint, dry_run=True)
        assert planned["status"] == "planned" and planned["dry_run"] is True
        assert (work / "keep").read_bytes() == b"cache"
        assert not list(scope.glob("cleanup-*"))


def test_registration_can_resume_under_a_new_controller_lock(tmp_path):
    scope, _source, _results = _scope(tmp_path)
    with storage.managed_run_lock(scope):
        first = storage.register_stage(scope, "upstream")
        work = Path(first["work_root"])
        (work / "existing-cache").write_text("existing generated product")
    with storage.managed_run_lock(scope):
        assert storage.register_stage(scope, "upstream") == first
        assert (work / "existing-cache").read_text() == "existing generated product"
        result = storage.cleanup_stage(scope, "upstream", native_status=SUCCESS,
                                      checkpoint_validator=_checkpoint)
        assert result["status"] == "completed"


def test_exclusive_lifetime_lock_cannot_be_bypassed_or_removed(tmp_path):
    scope, _source, _results = _scope(tmp_path)
    with pytest.raises(RuntimeError, match="held exclusive"):
        storage.register_stage(scope, "upstream")
    with storage.managed_run_lock(scope):
        original = (scope / storage._LOCK).read_bytes()
        with pytest.raises(FileExistsError):
            with storage.managed_run_lock(scope):
                pytest.fail("second controller acquired a live run")
        assert (scope / storage._LOCK).read_bytes() == original
        assert storage.register_stage(scope, "upstream")["stage"] == "upstream"
    with pytest.raises(RuntimeError, match="held exclusive"):
        storage.cleanup_stage(scope, "upstream", native_status=SUCCESS, checkpoint_validator=_checkpoint)


def test_unknown_directory_and_arbitrary_stage_roots_are_rejected(tmp_path):
    scope, source, _results = _scope(tmp_path)
    with storage.managed_run_lock(scope):
        with pytest.raises(ValueError, match="exact scope/work"):
            storage.register_stage(scope, "upstream", source)
        unknown = scope / "work" / "upstream"
        unknown.mkdir(parents=True)
        (unknown / "unknown").write_text("not owned")
        with pytest.raises(FileExistsError, match="never adopted"):
            storage.register_stage(scope, "upstream")
        assert (unknown / "unknown").read_text() == "not owned"
        with pytest.raises(ValueError, match="upstream or downstream"):
            storage.register_stage(scope, "../outside")


def test_protected_and_parent_traversal_roots_are_rejected(tmp_path):
    protected = tmp_path / "protected"
    protected.mkdir()
    for root in (protected, protected / "child", tmp_path):
        with pytest.raises(ValueError, match="protected"):
            storage.create_managed_scope(root, protected_roots=(protected,))
    with pytest.raises(ValueError, match="parent traversal"):
        storage.create_managed_scope(tmp_path / "nested" / ".." / "escape")
    with pytest.raises(ValueError):
        storage.create_managed_scope(Path(tmp_path.anchor))


def test_stage_symlinks_are_unlinked_without_touching_their_targets(tmp_path):
    scope, source, _results = _scope(tmp_path)
    with storage.managed_run_lock(scope):
        work = Path(storage.register_stage(scope, "upstream")["work_root"])
        try:
            (work / "input-directory").symlink_to(source, target_is_directory=True)
            (work / "input-file").symlink_to(source / "reads.fastq")
        except OSError as error:
            pytest.skip(f"Symlink privilege unavailable: {error}")
        result = storage.cleanup_stage(scope, "upstream", native_status=SUCCESS,
                                      checkpoint_validator=_checkpoint)
        assert result["accounting"]["regular_files"] == 0
        assert result["accounting"]["symlinks_or_reparse_points"] == 2
        assert not work.exists()
        assert (source / "reads.fastq").read_text() == "original generated input\n"


def test_symlinked_control_root_is_never_followed(tmp_path):
    scope, _source, _results = _scope(tmp_path)
    alias = tmp_path / "alias"
    try:
        alias.symlink_to(scope, target_is_directory=True)
    except OSError as error:
        pytest.skip(f"Symlink privilege unavailable: {error}")
    with pytest.raises(ValueError, match="traverse links"):
        with storage.managed_run_lock(alias):
            pytest.fail("must not follow an alias to a deletion scope")


@pytest.mark.skipif(os.name != "nt", reason="Windows directory junction boundary")
def test_windows_junction_is_removed_without_deleting_retained_target(tmp_path):
    import _winapi
    scope, source, _results = _scope(tmp_path)
    with storage.managed_run_lock(scope):
        work = Path(storage.register_stage(scope, "upstream")["work_root"])
        junction = work / "retained-source-junction"
        _winapi.CreateJunction(str(source), str(junction))
        result = storage.cleanup_stage(scope, "upstream", native_status=SUCCESS,
                                      checkpoint_validator=_checkpoint)
        assert result["accounting"]["symlinks_or_reparse_points"] == 1
        assert result["accounting"]["regular_files"] == 0
        assert not work.exists()
    assert (source / "reads.fastq").read_text() == "original generated input\n"


def test_hardlink_with_retained_external_name_does_not_claim_released_allocation(tmp_path):
    scope, source, _results = _scope(tmp_path)
    with storage.managed_run_lock(scope):
        work = Path(storage.register_stage(scope, "upstream")["work_root"])
        try:
            os.link(source / "reads.fastq", work / "linked-input")
        except OSError as error:
            pytest.skip(f"Hardlinks unavailable: {error}")
        result = storage.cleanup_stage(scope, "upstream", native_status=SUCCESS,
                                      checkpoint_validator=_checkpoint)
        assert result["accounting"]["reclaimable_allocated_bytes_estimate"] in (0, None)
    assert (source / "reads.fastq").is_file()


def test_completed_cleanup_is_idempotent_but_recreated_work_is_rejected(tmp_path):
    scope, _source, _results = _scope(tmp_path)
    with storage.managed_run_lock(scope):
        registered = storage.register_stage(scope, "upstream")
        work = Path(registered["work_root"])
        result = storage.cleanup_stage(scope, "upstream", native_status=SUCCESS,
                                      checkpoint_validator=_checkpoint)
        assert storage.register_stage(scope, "upstream") == registered
        assert not work.exists()
        again = storage.cleanup_stage(scope, "upstream", native_status=SUCCESS,
            checkpoint_validator=lambda _: pytest.fail("must not re-clean completed data"))
        assert again == result
        work.mkdir()
        (work / "new").write_text("foreign new file")
        with pytest.raises(ValueError, match="unexpected new data"):
            storage.register_stage(scope, "upstream")
        with pytest.raises(ValueError, match="unexpected new data"):
            storage.cleanup_stage(scope, "upstream", native_status=SUCCESS,
                                  checkpoint_validator=_checkpoint)
        assert (work / "new").read_text() == "foreign new file"


def test_missing_or_replaced_registered_stage_is_not_adopted_on_resume(tmp_path):
    scope, _source, _results = _scope(tmp_path)
    with storage.managed_run_lock(scope):
        work = Path(storage.register_stage(scope, "upstream")["work_root"])
        work.rename(scope / "saved-work")
        with pytest.raises(ValueError, match="missing or changed identity"):
            storage.register_stage(scope, "upstream")
        work.mkdir()
        with pytest.raises(ValueError, match="missing or changed identity"):
            storage.register_stage(scope, "upstream")


def test_identity_replacement_during_checkpoint_is_refused(tmp_path):
    scope, _source, _results = _scope(tmp_path)
    with storage.managed_run_lock(scope):
        work = Path(storage.register_stage(scope, "upstream")["work_root"])
        (work / "original").write_text("owned")
        moved = scope / "preserved-old-work"
        def changed(path):
            path.rename(moved)
            path.mkdir()
            (path / "new").write_text("new data")
            return {"status": "completed"}
        with pytest.raises(RuntimeError, match="changed during"):
            storage.cleanup_stage(scope, "upstream", native_status=SUCCESS, checkpoint_validator=changed)
        assert (work / "new").read_text() == "new data"
        assert (moved / "original").read_text() == "owned"


def test_changed_inventory_at_quarantine_stops_before_deletion_and_keeps_failure_claim(tmp_path, monkeypatch):
    scope, _source, _results = _scope(tmp_path)
    with storage.managed_run_lock(scope):
        work = Path(storage.register_stage(scope, "upstream")["work_root"])
        (work / "original").write_text("owned")
        verify = storage._verify_inventory
        calls = 0
        def alter_after_plan(path, rows):
            nonlocal calls
            calls += 1
            verify(path, rows)
            if calls == 1:
                (path / "new").write_text("changed after preflight")
        monkeypatch.setattr(storage, "_verify_inventory", alter_after_plan)
        with pytest.raises(RuntimeError, match="changed after cleanup planning"):
            storage.cleanup_stage(scope, "upstream", native_status=SUCCESS, checkpoint_validator=_checkpoint)
        pending = scope / "cleanup-upstream.pending.json"
        assert pending.is_file()
        quarantine = Path(json.loads(pending.read_text())["quarantine"])
        assert (quarantine / "original").read_text() == "owned"
        assert (quarantine / "new").read_text() == "changed after preflight"
        with pytest.raises(RuntimeError, match="unresolved outcome"):
            storage.register_stage(scope, "upstream")
        with pytest.raises(RuntimeError, match="unresolved outcome"):
            storage.cleanup_stage(scope, "upstream", native_status=SUCCESS, checkpoint_validator=_checkpoint)


def test_cleanup_record_cannot_be_rebound_to_another_scope(tmp_path):
    scope, _source, _results = _scope(tmp_path)
    foreign = tmp_path / "foreign"
    foreign.mkdir()
    (foreign / storage._MARKER).write_bytes((scope / storage._MARKER).read_bytes())
    with pytest.raises(ValueError, match="exact registered"):
        with storage.managed_run_lock(foreign):
            pytest.fail("foreign root adopted")
