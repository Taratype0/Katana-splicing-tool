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


def _write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2), encoding="utf-8")
    return path


def test_descriptor_and_native_project_use_their_own_relative_path_roots(tmp_path, monkeypatch):
    study = tmp_path / "study"
    exported = tmp_path / "different-export"
    cwd = tmp_path / "caller"
    cwd.mkdir()
    monkeypatch.chdir(cwd)
    project = {"project_id": "study-A", "name": "Study A", "samples": [], "contrasts": [], "output_root": "results"}
    manifest = _write_json(study / "katana.project.json", project)
    parent_sha = hashlib.sha256(manifest.read_bytes()).hexdigest()
    native_table = study / "tables/native.tsv"
    native_table.parent.mkdir()
    native_table.write_text("gene_id\tfdr\nNativeOnly\t0.01\n", encoding="utf-8")
    native_artifact = {"project_id": "study-A", "analysis_id": "differential.deg", "name": "Native DEG",
                       "status": "complete", "path": "tables/native.tsv"}
    _write_json(exported / "catalog/native.json", {"artifacts": [native_artifact]})
    asset = exported / "assets/extra.svg"
    asset.parent.mkdir()
    asset.write_text('<svg xmlns="http://www.w3.org/2000/svg"/>', encoding="utf-8")
    table = exported / "tables/details.tsv"
    table.parent.mkdir()
    table.write_text("event_id\nfive\n", encoding="utf-8")
    binding = {"project_id": "study-A", "status": "complete", "source_project_manifest_sha256": parent_sha}
    file_views = {**binding, "panels": [{"id": "file-panel", "path": "../assets/extra.svg", "svg": "../assets/extra.svg"}],
                  "tables": [{"path": "../tables/details.tsv"}],
                  "top10": [{"full_path": "../tables/details.tsv", "full_table": "../tables/details.tsv"}]}
    view_file = _write_json(exported / "views/completed.json", file_views)
    inline_views = {**binding, "panels": [{"id": "inline-panel", "image": "assets/extra.svg"}],
                    "tables": [{"source": "tables/details.tsv"}], "top10": [{"full_path": "tables/details.tsv"}]}
    child_project = {**project, "project_id": "supplemental-A", "metadata": {
        "source_project_id": "study-A", "source_project_manifest_SHA256": parent_sha}}
    child_descriptor = _write_json(exported / "supplemental/descriptor.json", {"project": child_project, "artifacts": []})
    child_manifest = _write_json(study / "supplemental.project.json", child_project)
    descriptor = {"project_manifest": "../study/katana.project.json", "artifact_catalog": "catalog/native.json",
                  "extra_assets": [{"project_id": "study-A", "status": "complete", "path": "assets/extra.svg",
                                    "metadata": {"source_project_manifest_sha256": parent_sha, "exports": {"png": "assets/extra.svg"}}}],
                  "extra_views": ["views/completed.json", inline_views],
                  "supplemental_sources": ["supplemental/descriptor.json",
                                           {"project_manifest": "../study/supplemental.project.json", "artifacts": []}]}
    descriptor_path = _write_json(exported / "descriptor.json", descriptor)
    originals = {path: path.read_bytes() for path in (manifest, view_file, child_descriptor, child_manifest, descriptor_path, asset, table, native_table)}
    loaded = load_source(descriptor_path)
    assert loaded["artifacts"][0]["path"] == str(native_table.resolve())
    assert loaded["extra_assets"][0]["path"] == str(asset.resolve())
    assert loaded["extra_assets"][0]["metadata"]["exports"]["png"] == str(asset.resolve())
    assert loaded["extra_views"][0]["panels"][0]["path"] == str(asset.resolve())
    assert loaded["extra_views"][0]["panels"][0]["svg"] == str(asset.resolve())
    assert loaded["extra_views"][0]["tables"][0]["path"] == str(table.resolve())
    assert loaded["extra_views"][0]["top10"][0]["full_path"] == str(table.resolve())
    assert loaded["extra_views"][0]["top10"][0]["full_table"] == str(table.resolve())
    assert loaded["extra_views"][1]["panels"][0]["image"] == str(asset.resolve())
    assert loaded["extra_views"][1]["tables"][0]["source"] == str(table.resolve())
    assert loaded["extra_views"][1]["top10"][0]["full_path"] == str(table.resolve())
    assert len(loaded["supplemental"]) == 2
    assert child_descriptor.resolve() in loaded["supplemental"][0]["source_files"]
    assert child_manifest.resolve() in loaded["supplemental"][1]["source_files"]
    assert all(path.read_bytes() == content for path, content in originals.items())


@pytest.mark.parametrize("catalog_key", ["artifact_catalog", "catalog_export", "native_artifact_catalog", "artifacts"])
def test_explicit_catalog_is_relative_to_caller_and_descriptor_catalog_to_descriptor(tmp_path, monkeypatch, catalog_key):
    project = {"project_id": "study-A", "samples": [], "contrasts": []}
    manifest = _write_json(tmp_path / "study/project.json", project)
    descriptor = _write_json(tmp_path / "export/descriptor.json", {"project_manifest": "../study/project.json", catalog_key: "catalog.json"})
    asset = manifest.parent / "native.tsv"
    asset.write_text("gene_id\nNative\n", encoding="utf-8")
    _write_json(descriptor.parent / "catalog.json", [{"name": "descriptor catalog", "path": "native.tsv"}])
    cwd = tmp_path / "caller"
    cwd.mkdir()
    _write_json(cwd / "catalog.json", [{"name": "explicit caller catalog", "path": "native.tsv"}])
    monkeypatch.chdir(cwd)
    from_descriptor = load_source(descriptor)
    explicit = load_source(descriptor, catalog="catalog.json")
    assert from_descriptor["artifacts"][0]["name"] == "descriptor catalog"
    assert explicit["artifacts"][0]["name"] == "explicit caller catalog"
    assert from_descriptor["artifacts"][0]["path"] == explicit["artifacts"][0]["path"] == str(asset.resolve())


def test_inline_descriptor_paths_are_relative_to_caller_without_mutating_input(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    descriptor, digest = _bound_descriptor()
    descriptor["extra_views"] = [{"project_id": "study-A", "status": "complete", "source_project_manifest_sha256": digest,
                                  "panels": [{"path": "figures/native.svg"}], "top10": [{"full_path": "tables/full.tsv"}]}]
    original = copy.deepcopy(descriptor)
    loaded = load_source(descriptor)
    assert loaded["extra_views"][0]["panels"][0]["path"] == str((tmp_path / "figures/native.svg").resolve())
    assert loaded["extra_views"][0]["top10"][0]["full_path"] == str((tmp_path / "tables/full.tsv").resolve())
    assert descriptor == original


def test_native_manifest_paths_bind_original_hash_and_protect_actual_result_root(tmp_path, monkeypatch):
    from katana_next.reports.engine import build_report

    study = tmp_path / "study"
    caller = tmp_path / "caller"
    caller.mkdir()
    monkeypatch.chdir(caller)
    project = {"project_id": "study-A", "name": "Study A", "input_root": "raw", "output_root": "results",
               "annotation": {"genome_build": "mm39", "species": "mus musculus", "gtf": "reference/genes.gtf",
                              "genome_fasta": "reference/genome.fa", "transcript_fasta": "reference/transcripts.fa"},
               "assets": [{"asset_id": "counts", "kind": "counts", "path": "tables/counts.tsv", "metadata": {"label": "run/v1"}}],
               "samples": [{"sample_id": "A", "biological_sample_id": "mouse-A", "condition": "condition/A",
                            "assets": {"quant_sf": "quant/A/quant.sf", "unknown_asset": "leave/unknown"}}], "contrasts": []}
    manifest = _write_json(study / "native.json", project)
    original_manifest = manifest.read_bytes()
    parent_hash = hashlib.sha256(original_manifest).hexdigest()
    descriptor = _write_json(tmp_path / "export/descriptor.json", {
        "project_manifest": "../study/native.json", "artifacts": [],
        "extra_assets": [{"project_id": "study-A", "path": "figures/extra.svg",
                          "metadata": {"source_project_manifest_sha256": parent_hash}}]})
    original_descriptor = descriptor.read_bytes()
    loaded = load_source(descriptor)
    native = loaded["project"]
    assert native["output_root"] == str((study / "results").resolve())
    assert native["input_root"] == str((study / "raw").resolve())
    for key, relative in (("gtf", "genes.gtf"), ("genome_fasta", "genome.fa"), ("transcript_fasta", "transcripts.fa")):
        assert native["annotation"][key] == str((study / "reference" / relative).resolve())
    assert native["annotation"]["genome_build"] == "mm39"
    assert native["annotation"]["species"] == "mus musculus"
    assert native["assets"][0]["path"] == str((study / "tables/counts.tsv").resolve())
    assert native["assets"][0]["metadata"]["label"] == "run/v1"
    assert native["samples"][0]["condition"] == "condition/A"
    assert native["samples"][0]["assets"]["quant_sf"] == str((study / "quant/A/quant.sf").resolve())
    assert native["samples"][0]["assets"]["unknown_asset"] == "leave/unknown"
    with pytest.raises(ValueError, match="separate from the source result"):
        build_report(descriptor, study / "results/report-export")
    assert not (study / "results/report-export").exists()
    assert manifest.read_bytes() == original_manifest
    assert descriptor.read_bytes() == original_descriptor


def test_inline_native_normalization_keeps_original_input_and_parent_identity(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    project = {"project_id": "study-A", "input_root": "raw", "output_root": "results", "samples": [], "contrasts": [],
               "annotation": {"genome_build": "unchanged/build", "gtf": "ref/genes.gtf"},
               "assets": [{"path": "matrix/counts.tsv"}]}
    parent_hash = hashlib.sha256(json.dumps(project, sort_keys=True, separators=(',', ':'), default=str).encode()).hexdigest()
    descriptor = {"project": project, "artifacts": [], "extra_views": [{
        "project_id": "study-A", "source_project_manifest_sha256": parent_hash,
        "panels": [{"path": "figures/stored.svg"}]}]}
    original = copy.deepcopy(descriptor)
    loaded = load_source(descriptor)
    assert loaded["project"]["output_root"] == str((tmp_path / "results").resolve())
    assert loaded["project"]["annotation"]["gtf"] == str((tmp_path / "ref/genes.gtf").resolve())
    assert loaded["project"]["assets"][0]["path"] == str((tmp_path / "matrix/counts.tsv").resolve())
    assert loaded["project"]["annotation"]["genome_build"] == "unchanged/build"
    assert descriptor == original


def test_remote_native_uris_are_retained_without_claiming_remote_io(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    project = {"project_id": "study-A", "input_root": "s3://reads/project", "output_root": "s3://results/project",
               "annotation": {"genome_build": "mm39", "genome_fasta": "s3://reference/mouse.fa"},
               "assets": [{"path": "s3://tables/counts.tsv"}],
               "samples": [{"assets": {"quant_sf": "s3://quant/sample/quant.sf"}}]}
    descriptor = {"project": project, "artifacts": []}
    loaded = load_source(descriptor)
    assert loaded["project"] == project
    with pytest.raises(ValueError, match="Remote result URIs are not read"):
        load_source({"project": project})


def test_default_completed_catalog_uses_native_relative_output_root(tmp_path, monkeypatch):
    study = tmp_path / "study"
    manifest = _write_json(study / "native.json", {"project_id": "study-A", "output_root": "results"})
    native_export = _write_json(study / "results/.katana/artifact_manifest.json", {"artifacts": []})
    caller = tmp_path / "caller"
    caller.mkdir()
    monkeypatch.chdir(caller)
    loaded = load_source(manifest)
    assert native_export.resolve() in loaded["source_files"]
    assert loaded["project"]["output_root"] == str((study / "results").resolve())
