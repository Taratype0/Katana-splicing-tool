"""An unchanged edited design must not invalidate science on a label toggle."""
from pathlib import Path
from types import SimpleNamespace
import hashlib
import os

from katana_cosplice.config import CoSpliceConfig
from katana_cosplice.pipeline import CoSplicePipeline
from katana_next.services.design_editor import sync_cosplice_design


def test_unchanged_edited_metadata_preserves_mtime_and_real_stage_fingerprint(tmp_path):
    original = tmp_path / "input.tsv"
    original.write_text("sample_id\tcondition\nS1\told\nS2\told\n", encoding="utf-8")
    original_hash = hashlib.sha256(original.read_bytes()).hexdigest()
    samples = [SimpleNamespace(sample_id=f"S{i}", condition="A", biological_sample_id=f"D{i}",
                               batch=None, attributes={}) for i in (1, 2)]
    manifest = SimpleNamespace(metadata={"design_revision": "edited-design"}, samples=samples)
    config = CoSpliceConfig("design-cache", None, original, tmp_path / "engine")
    sync_cosplice_design(config, manifest, config.output_dir)
    resolved = Path(config.sample_metadata)
    os.utime(resolved, ns=(1_700_000_000_000_000_000, 1_700_000_000_000_000_000))
    before = resolved.read_bytes(), resolved.stat().st_mtime_ns
    engine = CoSplicePipeline(config)
    key = engine._fingerprint("matrix", [resolved], {"design": "fixed"})

    config.consensus.apply_threshold_labels = True
    sync_cosplice_design(config, manifest, config.output_dir)

    assert (resolved.read_bytes(), resolved.stat().st_mtime_ns) == before
    assert engine._fingerprint("matrix", [resolved], {"design": "fixed"}) == key
    assert hashlib.sha256(original.read_bytes()).hexdigest() == original_hash
    samples[0].condition = "B"
    sync_cosplice_design(config, manifest, config.output_dir)
    assert resolved.read_bytes() != before[0]
    assert engine._fingerprint("matrix", [resolved], {"design": "fixed"}) != key
    assert hashlib.sha256(original.read_bytes()).hexdigest() == original_hash
