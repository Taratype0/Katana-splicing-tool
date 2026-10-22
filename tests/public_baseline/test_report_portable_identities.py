"""Repeated contrast names and external Top 10 sources stay correct offline."""
from pathlib import Path
import copy
import hashlib
import json

import pandas as pd
import pytest

from katana_next.reports.engine import Builder
from katana_next.reports.sources import load_source


def test_same_named_tables_keep_distinct_complete_browsers_per_comparison(tmp_path):
    output = tmp_path / "report"
    output.mkdir()
    builder = Builder({"project": {}}, output, "stable_significant")
    first, second = tmp_path / "A.parquet", tmp_path / "B.parquet"
    pd.DataFrame({"gene_id": ["A_ONLY"], "fdr": [.01]}).to_parquet(first)
    pd.DataFrame({"gene_id": ["B_ONLY"], "fdr": [.02]}).to_parquet(second)
    a = builder.table({"name": "DEG-full", "path": str(first), "contrast_id": "A"})
    b = builder.table({"name": "DEG-full", "path": str(second), "contrast_id": "B"})
    assert a["view"] != b["view"]
    assert "A_ONLY" in (output / a["view"]).read_text(encoding="utf-8")
    assert "B_ONLY" not in (output / a["view"]).read_text(encoding="utf-8")
    assert "B_ONLY" in (output / b["view"]).read_text(encoding="utf-8")


def test_external_top10_complete_result_is_copied_and_not_mutated(tmp_path):
    output = tmp_path / "report"
    output.mkdir()
    source = tmp_path / "external-results.tsv"
    source.write_text("gene_id\tpvalue\ng1\t5e-310\n", encoding="utf-8")
    entry = {"id": "top10-A-DEG", "kind": "DEG", "full_path": str(source), "full_table": str(source)}
    builder = Builder({"project": {}}, output, "stable_significant")
    builder.add_views({"top10": [entry]})
    assert entry["full_path"] == str(source)
    copied = builder.top10[0]
    assert not Path(copied["full_path"]).is_absolute()
    assert (output / copied["full_path"]).read_bytes() == source.read_bytes()
    assert copied["full_table"] == copied["full_path"]


def _bound_descriptor():
    project = {"project_id": "study-A", "name": "Study A", "samples": [], "contrasts": []}
    digest = hashlib.sha256(json.dumps(project, sort_keys=True, separators=(',', ':'), default=str).encode()).hexdigest()
    return {"project": project, "artifacts": []}, digest


@pytest.mark.parametrize("slot", ["extra_views", "extra_assets"])
@pytest.mark.parametrize("bad", ["other_project", "other_hash", "running"])
def test_extra_result_cannot_borrow_another_cohort_or_incomplete_run(slot, bad):
    descriptor, digest = _bound_descriptor()
    record = {"project_id": "study-A", "status": "complete", "source_project_manifest_sha256": digest}
    if bad == "other_project": record["project_id"] = "study-B"
    if bad == "other_hash": record["source_project_manifest_sha256"] = "different manifest"
    if bad == "running": record["status"] = "running"
    if slot == "extra_assets":
        record["metadata"] = {"source_project_manifest_sha256": record.pop("source_project_manifest_sha256")}
    descriptor[slot] = [record]
    original = copy.deepcopy(descriptor)
    with pytest.raises(ValueError):
        load_source(descriptor)
    assert descriptor == original


def test_loading_extra_views_does_not_mutate_supplied_nested_metadata():
    descriptor, digest = _bound_descriptor()
    descriptor["extra_views"] = [{"project_id": "study-A", "status": "complete", "source_project_manifest_sha256": digest,
                                  "panels": [{"id": "kept", "metadata": {"policy": ["original"]}}]}]
    original = copy.deepcopy(descriptor)
    loaded = load_source(descriptor)
    loaded["extra_views"][0]["panels"][0]["metadata"]["policy"].append("report copy only")
    assert descriptor == original


def test_unrecorded_rbp_foreground_cannot_substitute_for_current_rs(tmp_path):
    image = tmp_path / "historical_rbp.svg"
    image.write_text('<svg xmlns="http://www.w3.org/2000/svg" width="200" height="100"></svg>', encoding="utf-8")
    artifact = {"artifact_id": "old-rbp", "analysis_id": "prioritize.rbp_support", "name": "Historical RBP identity",
                "path": str(image), "status": "complete", "metadata": {}, "parameters": {}, "contrast_ids": ["B_vs_A"]}
    original = copy.deepcopy(artifact)
    output = tmp_path / "report"
    output.mkdir()
    builder = Builder({"project": {}}, output, "stable_significant")
    builder.collect_rich({"artifacts": [artifact]})
    assert len(builder.panels) == 1 and builder.panels[0]["output_only"]
    assert "Historical or unrecorded foreground" in builder.panels[0]["caption"]
    assert builder.panels[0]["selection_scope"]["recorded_event_set"] is None
    assert artifact == original
