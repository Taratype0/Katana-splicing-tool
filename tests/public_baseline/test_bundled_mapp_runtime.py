from copy import deepcopy
import json
from pathlib import Path
import shutil

import pytest

from katana_next.flow import bundled_runtime as runtime


def test_bundled_binding_is_portable_and_resolves_inside_each_install(tmp_path, monkeypatch):
    root, _ = runtime.verified_mapp_sources()
    inputs = runtime.mapp_input_defaults({})
    assert "native_root" not in inputs
    spec = {"adapter": "mapp.prepare", "inputs": inputs}
    before = deepcopy(spec)
    relocated = tmp_path / "different-install" / "MAPP"
    shutil.copytree(root, relocated)
    monkeypatch.setattr(runtime, "_source_root", lambda: relocated)
    resolved = runtime.resolve_bundled_stage(spec)
    assert Path(resolved["inputs"]["native_root"]) == relocated
    assert len(resolved["inputs"]["native_source_files"]) == 10
    assert all(Path(path).is_file() for path in resolved["inputs"]["native_source_files"])
    assert spec == before


def test_explicit_old_source_does_not_silently_change():
    inputs = {"native_root": "/declared/author-checkout", "native_commit": "old-revision"}
    assert runtime.mapp_input_defaults(inputs) == inputs
    task = {"adapter": "mapp.maei", "inputs": inputs}
    assert runtime.resolve_bundled_stage(task) is task


def test_changed_binding_and_changed_installed_source_fail_before_output(tmp_path, monkeypatch):
    root, _ = runtime.verified_mapp_sources()
    inputs = runtime.mapp_input_defaults({})
    bad = deepcopy(inputs)
    bad["counter_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="content binding"):
        runtime.resolve_bundled_stage({"adapter": "mapp.prepare", "inputs": bad})
    relocated = tmp_path / "broken-install"
    shutil.copytree(root, relocated)
    (relocated / runtime.COUNTER).write_bytes(b"modified source")
    monkeypatch.setattr(runtime, "_source_root", lambda: relocated)
    from katana_next.flow.stages import execute_stage
    output = tmp_path / "no-output"
    with pytest.raises(ValueError, match="source changed"):
        execute_stage({"adapter": "mapp.prepare", "inputs": inputs}, output)
    assert not output.exists()


def test_all_author_files_have_unchanged_recorded_bytes():
    root, descriptor = runtime.verified_mapp_sources()
    assert len(descriptor["source_sha256"]) == 7
    manifest = json.loads((root / "SOURCE.json").read_text(encoding="utf-8"))
    assert manifest["modified"] is False
    assert manifest["commit"] == runtime.MAPP_COMMIT
