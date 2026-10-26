"""Public flat options reach the real importer/planner without study JSON."""
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("katana_configure", ROOT / "pipeline/bin/katana_configure.py")
configure = importlib.util.module_from_spec(spec)
spec.loader.exec_module(configure)


def disabled():
    return {"skip_" + key: True for key in configure.ANALYSES}


def test_default_resources_and_runtime_are_generated(tmp_path):
    names = ["gene_sets", "rbp_gene_sets", "motif_database", "clip_peaks"]
    paths = []
    for name in names:
        path = tmp_path / (name + " with spaces.txt")
        path.write_text("fixture resource\n")
        paths.append(path)
    original = {p: p.read_bytes() for p in paths}
    template, config = configure.generate({"species": "Mus musculus", "resource_names": names}, paths, tmp_path / "generated")
    assert config["network"] == {"scope": "contrast", "bootstrap": {"enabled": False},
                                 "parameters": {"hub_metric": "abs_kME", "min_module_size": 20, "merge_cut_height": .25}}
    assert config["analyses"]["coverage"] is False
    assert all(config["analyses"][key] for key in configure.ANALYSES)
    assert config["event_set"] == "stable_significant"
    assert len(template["assets"]) == 5
    assert json.loads((tmp_path / "generated/fimo-runtime.json").read_text()) == {"python": "python3"}
    assert template["assets"][0]["metadata"]["version"].startswith("sha256:")
    assert {p: p.read_bytes() for p in paths} == original


@pytest.mark.parametrize("updates,match", [
    ({"deg_fdr": 2}, "deg_fdr"), ({"deg_fdr": float("nan")}, "finite"),
    ({"deg_log2fc": -1}, "deg_log2fc"), ({"das_dpsi": True}, "finite"),
    ({"skip_mapp": "false"}, "true or false"), ({"rmats_mode": "guess"}, "rmats_mode"),
    ({"event_set": "all"}, "event_set"), ({"contrast_ids": "A,A"}, "duplicates"),
    ({"with_network_bootstrap": "false"}, "true or false"),
    ({"network_min_module_size": 2.5}, "integer"),
    ({"network_merge_cut_height": 1}, "strictly"),
])
def test_invalid_options_fail_before_science(tmp_path, updates, match):
    with pytest.raises(ValueError, match=match):
        configure.generate({**disabled(), **updates}, [], tmp_path / "output")
    assert not (tmp_path / "output").exists()


def test_missing_resources_explain_flags(tmp_path):
    with pytest.raises(ValueError, match="--gene_sets"):
        configure.generate({"species": "Mus musculus"}, [], tmp_path / "output")
    assert not (tmp_path / "output").exists()


def test_resource_metadata_conflicts_are_not_relabelled(tmp_path):
    path = tmp_path / "motifs.json"
    path.write_text("{}")
    with pytest.raises(ValueError, match="species conflicts"):
        configure.generate({**disabled(), "species": "Mus musculus", "resource_names": ["motif_database"],
                            "resource_metadata": {"motif_database": {"species": "Homo sapiens"}}},
                           [path], tmp_path / "output")


def test_flat_thresholds_reach_real_project_plan(tmp_path):
    from katana_next.full_fixture import create_full_fixture
    from katana_next.flow.planner import create_plan
    project = create_full_fixture(tmp_path / "study")
    # Fixture API returns the project's manifest file.
    _, config = configure.generate({**disabled(), "mode": "downstream-only", "deg_fdr": .01,
                                    "deg_log2fc": .5, "das_fdr": .02, "das_dpsi": .2,
                                    "dtu_fdr": .03, "event_set": "significant"}, [], tmp_path / "configuration")
    plan = create_plan(project, config)
    tasks = {t["workflow_id"]: t for t in plan["tasks"]}
    assert tasks["differential.deg"]["parameters"]["fdr"] == .01
    assert tasks["differential.deg"]["parameters"]["min_abs_log2fc"] == .5
    assert tasks["quality.stable_as"]["parameters"]["min_abs_dpsi"] == .2
    assert tasks["quality.stable_as"]["parameters"]["fdr"] == .02
    assert tasks["differential.dtu"]["parameters"]["fdr"] == .03
    assert tasks["prioritize.evidence"]["parameters"]["deg_fdr"] == .01
    assert tasks["prioritize.evidence"]["parameters"]["dtu_fdr"] == .03
    assert not any(t["workflow_id"].startswith(("flow.network", "flow.mapp", "flow.qc.genebody")) for t in plan["tasks"])


def test_generated_template_imports_true_bundle_without_statistics_changes(tmp_path):
    from katana_next.flow.importer import import_bundle
    fixture_spec = importlib.util.spec_from_file_location("flow_import_fixture", Path(__file__).with_name("test_flow_import.py"))
    fixture = importlib.util.module_from_spec(fixture_spec)
    fixture_spec.loader.exec_module(fixture)
    root, manifest, _ = fixture.bundle(tmp_path)
    original = (root / "SE.MATS.JCEC.txt").read_bytes()
    _, config = configure.generate({**disabled(), "study_name": "Example", "species": "test"}, [], tmp_path / "configuration")
    imported = import_bundle(root, tmp_path / "configuration/project-template.json", tmp_path / "imported",
                             analysis_configuration=config)
    data = json.loads(imported.read_text())
    assert data["name"] == "Example"
    assert data["contrasts"][0]["contrast_id"] == "B_vs_A"
    assert data["assets"][0]["metadata"]["reverse_direction"] is True
    assert (root / "SE.MATS.JCEC.txt").read_bytes() == original


def test_generated_runtime_and_resource_paths_survive_bundle_import(tmp_path):
    from katana_next.flow.importer import import_bundle
    from katana_next.infrastructure.manifest_store import ManifestStore
    fixture_spec = importlib.util.spec_from_file_location("flow_import_resources_fixture", Path(__file__).with_name("test_flow_import.py"))
    fixture = importlib.util.module_from_spec(fixture_spec)
    fixture_spec.loader.exec_module(fixture)
    root, _, _ = fixture.bundle(tmp_path)
    names = ["gene_sets", "rbp_gene_sets", "motif_database", "clip_peaks"]
    paths = []
    for name in names:
        path = tmp_path / (name + ".fixture")
        path.write_text(name)
        paths.append(path)
    _, config = configure.generate({"species": "test", "resource_names": names}, paths, tmp_path / "configuration")
    imported = import_bundle(root, tmp_path / "configuration/project-template.json", tmp_path / "imported",
                             analysis_configuration=config)
    _, manifest = ManifestStore().load(imported)
    assets = {a.asset_id: a for a in manifest.assets}
    assert all(assets[name].path.read_text() == name for name in names)
    assert json.loads(assets["fimo_runtime"].path.read_text()) == {"python": "python3"}


def test_existing_project_resources_not_replaced(tmp_path):
    path = tmp_path / "gene_sets.gmt"
    path.write_text("set\tdescription\tg1\n")
    with pytest.raises(ValueError, match="uses the resources"):
        configure.generate({"mode": "downstream-only", "resource_names": ["gene_sets"]}, [path], tmp_path / "config")


def test_network_overrides_reach_the_engine_config(tmp_path):
    from katana_cosplice.config import CoSpliceConfig
    from katana_next.workflows.programs import configure_production_modules
    _, settings = configure.generate({**disabled(), "network_min_module_size": 30,
                                       "network_merge_cut_height": .15}, [], tmp_path / "config")
    config = CoSpliceConfig(project_name="fixture", rmats_root=tmp_path,
                           sample_metadata=tmp_path / "samples.tsv", output_dir=tmp_path / "network")
    configure_production_modules(config, settings["network"]["parameters"])
    assert config.network.min_module_size == 30
    assert config.network.wgcna_merge_cut_height == .15


def test_complete_yaml_covers_every_public_native_option():
    import re
    import yaml
    configuration = (ROOT / "pipeline/nextflow.config").read_text().split("params {", 1)[1].split("\n}\n", 1)[0]
    names = set(re.findall(r"^\s*(\w+)\s*=", configuration, flags=re.MULTILINE))
    internal = {"project_template", "flow_config", "flow_plan", "flow_executable", "flow_report_revision", "managed_preflight", "help"}
    template = yaml.safe_load((ROOT / "examples/nextflow/params.complete.yaml").read_text())
    assert set(template) == names - internal
    assert template["with_genebody_coverage"] is False
    assert template["with_network_bootstrap"] is False
    assert template["network_min_module_size"] == 20
    assert template["network_merge_cut_height"] == .25
    assert "study.yml" not in template


def test_deg_threshold_does_not_change_dtu_evidence():
    import pandas as pd
    from katana_next.science.bulk import build_candidate_evidence
    deg = pd.DataFrame([{"comparison_id": "B_vs_A", "gene_id": "g1", "gene_symbol": "G1", "log2fc": 2., "deg_fdr": .04}])
    dtu = pd.DataFrame([{"comparison_id": "B_vs_A", "gene_id": "g1", "gene_symbol": "G1", "transcript_id": "t1",
                         "source_kind": "dtu_result", "dtu_fdr": .025}])
    result = build_candidate_evidence(pd.DataFrame(), deg, dtu, {"deg_fdr": .01, "dtu_fdr": .03})
    assert not bool(result.iloc[0].deg_supported)
    assert bool(result.iloc[0].dtu_supported)
    legacy = build_candidate_evidence(pd.DataFrame(), deg, dtu, {"fdr": .02})
    assert not bool(legacy.iloc[0].deg_supported)
    assert not bool(legacy.iloc[0].dtu_supported)
