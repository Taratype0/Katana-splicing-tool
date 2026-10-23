"""No Docker/network needed: verify fail-closed source preparation."""
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest


def builder(tmp_path, monkeypatch):
    repo = Path(__file__).resolve().parents[2]
    spec = importlib.util.spec_from_file_location("dexseq_build", repo / "scripts/build_dexseq_container.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    root = tmp_path / "source"
    location = root / "containers/dexseq-tools"
    location.mkdir(parents=True)
    base = "example@sha256:" + "a" * 64
    content = b"test tool source\n"
    lock = {"schema": "katana.dexseq_sources/v1", "base_image": base,
            "source_url_prefix": "https://example.invalid/fixed/", "files": [{"name": "tool.py",
            "url": "https://example.invalid/fixed/tool.py", "sha256": hashlib.sha256(content).hexdigest()}]}
    (location / "sources.lock.json").write_text(json.dumps(lock))
    for name in ("LICENSE.GPL-3.0", "THIRD_PARTY.md", "check_runtime.py"):
        (location / name).write_text("fixture\n")
    (location / "Dockerfile").write_text("FROM " + base + "\n")
    monkeypatch.setattr(module, "ROOT", root)
    return module, content


def test_pinned_sources_preserved_and_existing_context_rejected(tmp_path, monkeypatch):
    module, data = builder(tmp_path, monkeypatch)
    destination = tmp_path / "build"
    module.prepare_context(destination, download=lambda _: data)
    assert (destination / "gpl-scripts/tool.py").read_bytes() == data
    with pytest.raises(ValueError, match="new build context"):
        module.prepare_context(destination, download=lambda _: data)
    with pytest.raises(ValueError, match="new build context"):
        module.prepare_context(module.ROOT / "inside", download=lambda _: data)


def test_changed_download_is_rejected_before_build(tmp_path, monkeypatch):
    module, _ = builder(tmp_path, monkeypatch)
    with pytest.raises(ValueError, match="checksum mismatch"):
        module.prepare_context(tmp_path / "bad", download=lambda _: b"changed source")
    assert not (tmp_path / "bad/gpl-scripts/tool.py").exists()
