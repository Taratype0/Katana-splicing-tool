"""A failed complete-report gate must not leave a deliverable REPORT.zip.

The renderer and ZIP files are real. Collection admission and gap responses are
controlled here to test archive/receipt handling independently of scientific
acceptance, which has its own contract tests.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import zipfile

import pytest

from katana_next.flow.cli import build_complete_report
from katana_next.flow.model import file_sha256, seal_plan, write_json


def _render_case(tmp_path, monkeypatch, *, gaps=(), missing_chapter=False, missing_manifest=False, gate_error=False, render_error=False):
    import katana_next.flow.acceptance as acceptance
    import katana_next.flow.collect as collection_module
    import katana_next.reports as reports
    from katana_next.full_fixture import create_full_fixture

    source = create_full_fixture(tmp_path / "generated-source")
    collection = tmp_path / "collection"
    collection.mkdir()
    descriptor = {"project_manifest": str(source), "artifacts": []}
    write_json(collection / "report-source.json", descriptor)
    plan = seal_plan({"project_manifest": str(source), "project_manifest_sha256": file_sha256(source),
                      "report_profile": "complete", "report_specification": "seven-chapter-test-contract",
                      "tasks": [{"id": "generated-collection", "workflow_id": "flow.test", "dependencies": []}]})
    write_json(collection / "plan.json", plan)
    monkeypatch.setattr(collection_module, "validate_collection", lambda _: {"status": "completed"})
    monkeypatch.setattr(collection_module, "materialize_report_descriptor", lambda _: descriptor)
    actual_build = reports.build_report
    observed = {}

    def render(*args, **kwargs):
        manifest = actual_build(*args, **kwargs)
        out = Path(args[1])
        archive = out / "REPORT.zip"
        assert kwargs["create_archive"] is False and not archive.exists()
        observed["rendered_files"] = {
            str(path.relative_to(out)): file_sha256(path)
            for path in out.rglob("*") if path.is_file() and path != archive
        }
        # The gate's chapter/artifact branches are distinct from the detailed
        # scientific-gap contract. Keep the actual rendered bytes untouched
        # except in the deliberate missing-manifest case.
        manifest["panels"] = [{"chapter": number, "image": f"chapter-{number}.svg"}
                              for number in range(1, 7 if missing_chapter else 8)]
        if missing_manifest:
            (out / "report-manifest.json").unlink()
            observed["rendered_files"].pop("report-manifest.json")
        if render_error:
            raise ValueError("Late renderer failure")
        return manifest

    def gate(*args, **kwargs):
        if gate_error:
            raise ValueError("Malformed acceptance evidence")
        return list(gaps)

    from katana_next.reports import refresh
    original_archive = refresh.write_archive
    def record_archive(*args, **kwargs):
        result = original_archive(*args, **kwargs)
        observed["zip_sha256"] = file_sha256(Path(args[0]) / "REPORT.zip")
        return result
    monkeypatch.setattr(refresh, "write_archive", record_archive)
    monkeypatch.setattr(reports, "build_report", render)
    monkeypatch.setattr(acceptance, "required_report_gaps", gate)
    return collection, observed


@pytest.mark.parametrize("mode,message", [
    ("gap", "Required fixed-report analyses"),
    ("chapter", "missing rendered figures"),
    ("manifest", "Incomplete final report artifact"),
    ("gate_error", "Malformed acceptance evidence"),
    ("render_error", "Late renderer failure"),
])
def test_real_render_is_preserved_without_compressing_failed_report(tmp_path, monkeypatch, mode, message):
    collection, observed = _render_case(
        tmp_path, monkeypatch,
        gaps=[{"component": "network", "reason": "Missing required comparison result"}] if mode == "gap" else [],
        missing_chapter=mode == "chapter", missing_manifest=mode == "manifest", gate_error=mode == "gate_error",
        render_error=mode == "render_error")
    output = tmp_path / "new-report"
    with pytest.raises(ValueError, match=message):
        build_complete_report(collection, output)
    assert not (output / "REPORT.zip").exists()
    assert not (output / "flow.acceptance.json").exists()
    diagnostic = output / "REPORT.diagnostic.zip"
    assert not diagnostic.exists() and "zip_sha256" not in observed
    for relative, expected in observed["rendered_files"].items():
        assert file_sha256(output / relative) == expected
    receipt = json.loads((output / "flow.acceptance.failed.json").read_text())
    assert receipt["status"] == "failed" and message in receipt["error"]
    assert receipt["diagnostic_archive"] is None
    before = {p.name: file_sha256(p) for p in output.iterdir() if p.is_file()}
    with pytest.raises(FileExistsError, match="existing report and acceptance evidence"):
        build_complete_report(collection, output)
    assert before == {p.name: file_sha256(p) for p in output.iterdir() if p.is_file()}


def test_success_keeps_actual_archive_and_writes_completed_receipt(tmp_path, monkeypatch):
    collection, observed = _render_case(tmp_path, monkeypatch)
    output = tmp_path / "new-report"
    result = build_complete_report(collection, output)
    assert result["status"] == "completed"
    assert file_sha256(output / "REPORT.zip") == result["archive_sha256"] == observed["zip_sha256"]
    assert json.loads((output / "flow.acceptance.json").read_text()) == result
    assert not (output / "REPORT.diagnostic.zip").exists()
    assert not (output / "flow.acceptance.failed.json").exists()


@pytest.mark.parametrize("phase", ["archive_hash", "receipt_write", "receipt_publish"])
def test_final_io_failure_cannot_leave_deliverable_archive(tmp_path, monkeypatch, phase):
    from contextlib import contextmanager
    import katana_next.flow.model as model
    collection, observed = _render_case(tmp_path, monkeypatch)
    output = tmp_path / "new-report"
    original_hash, original_open, original_rename = model.file_sha256, Path.open, Path.rename

    def checked_hash(path):
        if Path(path) == output.resolve() / "REPORT.zip":
            raise OSError("Injected archive read failure")
        return original_hash(path)

    @contextmanager
    def interrupted_write(path, *args, **kwargs):
        with original_open(path, *args, **kwargs) as handle:
            handle.write('{"status": "comp')
            handle.flush()
            raise OSError("Injected disk full while writing success receipt")
            yield handle

    def checked_open(path, *args, **kwargs):
        if path.name == ".flow.acceptance.pending.json":
            return interrupted_write(path, *args, **kwargs)
        return original_open(path, *args, **kwargs)

    def checked_rename(path, target):
        if path.name == ".flow.acceptance.pending.json":
            raise OSError("Injected receipt publish failure")
        return original_rename(path, target)

    if phase == "archive_hash":
        monkeypatch.setattr(model, "file_sha256", checked_hash)
    elif phase == "receipt_write":
        monkeypatch.setattr(Path, "open", checked_open)
    else:
        monkeypatch.setattr(Path, "rename", checked_rename)
    with pytest.raises(OSError, match="Injected"):
        build_complete_report(collection, output)
    assert not (output / "REPORT.zip").exists()
    assert not (output / "flow.acceptance.json").exists()
    assert file_sha256(output / "REPORT.diagnostic.zip") == observed["zip_sha256"]
    assert json.loads((output / "flow.acceptance.failed.json").read_text())["status"] == "failed"


@pytest.mark.parametrize("name", ["flow.acceptance.json", "flow.acceptance.failed.json", "REPORT.zip", "REPORT.diagnostic.zip"])
def test_existing_success_or_failure_evidence_is_never_overwritten(tmp_path, name):
    output = tmp_path / "old-report"
    output.mkdir()
    target = output / name
    original = b"Preserved earlier evidence\n"
    target.write_bytes(original)
    digest = hashlib.sha256(original).hexdigest()
    with pytest.raises(FileExistsError, match="existing report and acceptance evidence"):
        build_complete_report(tmp_path / "not-read", output)
    assert file_sha256(target) == digest
    assert [p.name for p in output.iterdir()] == [name]
