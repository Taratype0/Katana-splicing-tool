import json

import pytest

from katana_next.flow.importer import import_bundle
from katana_next.flow.model import file_sha256, write_json


def bundle(tmp_path):
    root = tmp_path / "upstream"
    root.mkdir()
    native = root / "SE.MATS.JCEC.txt"
    native.write_text("ID\tFDR\tIncLevelDifference\n1\t0.01\t-0.2\n", encoding="utf-8")
    manifest = {"schema": "katana-pipeline-bundle/v1", "status": "complete", "reference": {"genome_build": "test-build"},
                "samples": [{"sample_id": "a", "condition": "A", "subject": "donor-a"},
                            {"sample_id": "b", "condition": "B", "subject": "donor-b"}],
                "comparisons": [{"comparison_id": "B_vs_A", "treatment": "B", "control": "A",
                                 "rmats_group1": "A", "rmats_group2": "B", "deg_target": "B", "deg_reference": "A"}],
                "rmats_sample_order": {"A": ["a"], "B": ["b"]},
                "artifacts": [{"kind": "rmats", "relative_path": native.name, "comparison_id": "B_vs_A",
                               "event_type": "SE", "mode": "JCEC", "status": "available",
                               "bytes": native.stat().st_size, "sha256": file_sha256(native)}]}
    write_json(root / "manifest.json", manifest)
    template = write_json(tmp_path / "template.json", {"metadata": {"species": "test"}})
    return root, manifest, template


def test_import_retains_native_direction_and_independent_units(tmp_path):
    root, manifest, template = bundle(tmp_path)
    original = file_sha256(root / "SE.MATS.JCEC.txt")
    project = import_bundle(root, template, tmp_path / "workspace")
    value = json.loads(project.read_text("utf-8"))
    asset = value["assets"][0]
    assert asset["metadata"]["reverse_direction"] is True
    assert asset["metadata"]["numerator_sample_ids"] == ["b"]
    assert value["samples"][0]["biological_sample_id"] == "donor-a"
    assert value["input_root"] == "bundle"
    assert file_sha256(project.parent / "bundle" / asset["path"]) == original
    assert file_sha256(root / "SE.MATS.JCEC.txt") == original
    with pytest.raises(ValueError, match="exists"):
        import_bundle(root, template, project.parent)


def test_import_rejects_changed_native_results_and_bad_order(tmp_path):
    root, manifest, template = bundle(tmp_path)
    (root / "SE.MATS.JCEC.txt").write_text("changed", encoding="utf-8")
    with pytest.raises(ValueError, match="SHA256"):
        import_bundle(root, template, tmp_path / "bad")
    assert not (tmp_path / "bad").exists()


def test_import_rejects_ambiguous_rmats_sample_vectors(tmp_path):
    root, manifest, template = bundle(tmp_path)
    manifest["rmats_sample_order"]["A"] = ["a", "a"]
    write_json(root / "manifest.json", manifest)
    with pytest.raises(ValueError, match="vector order"):
        import_bundle(root, template, tmp_path / "bad")
