"""Approved module controls: inheritance, cache separation and actual provenance.

The workflow adapter is exercised with an explicitly fake engine stage. These
checks do not fit WGCNA or claim a new scientific validation experiment.
"""
from __future__ import annotations

import hashlib
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
import yaml

from katana_cosplice.config import CoSpliceConfig, load_config
from katana_next.full_fixture import create_full_fixture
from katana_next.infrastructure.fingerprints import fingerprint_file, fingerprint_value
from katana_next.services.project_service import ProjectService
from katana_next.workflows import programs


def _config(tmp_path):
    return CoSpliceConfig("module-control-fixture", None, tmp_path / "samples.tsv", tmp_path / "output")


def test_default_and_custom_yaml_module_values_survive_blank_runtime_controls(tmp_path):
    config = _config(tmp_path)
    expected = {"min_module_size": 20, "merge_cut_height": .25}
    assert programs.network_module_parameters(config) == expected
    config.network.min_module_size = 37
    config.network.wgcna_merge_cut_height = .12
    path = tmp_path / "custom.yaml"
    programs.write_production_configuration(config, path)
    reread = load_config(path)
    actual, sources, changes = programs.configure_production_modules(
        reread, programs.CoSpliceWorkflow().normalize_parameters({})
    )
    assert actual == {"min_module_size": 37, "merge_cut_height": .12}
    assert changes == {}
    assert set(sources.values()) == {"project_configuration_or_engine_default"}
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert payload["network_module_parameters"] == actual


@pytest.mark.parametrize("key,value", [
    ("min_module_size", 1), ("min_module_size", 20.5), ("min_module_size", True),
    ("min_module_size", float("nan")), ("merge_cut_height", 0),
    ("merge_cut_height", 1), ("merge_cut_height", False),
    ("merge_cut_height", float("inf")), ("merge_cut_height", ""),
])
def test_invalid_module_controls_fail_before_any_stage(tmp_path, key, value):
    with pytest.raises(ValueError, match=key):
        programs.CoSpliceWorkflow().normalize_parameters({key: value})
    with pytest.raises(ValueError, match=key):
        programs.configure_production_modules(_config(tmp_path), {key: value})
    assert not (tmp_path / "output").exists()


def test_explicit_override_and_equal_value_cache_semantics(tmp_path):
    config = _config(tmp_path)
    values = programs.CoSpliceWorkflow().normalize_parameters({"min_module_size": "31", "merge_cut_height": ".2"})
    assert type(values["min_module_size"]) is int
    actual, sources, changes = programs.configure_production_modules(config, values)
    assert actual == changes == {"min_module_size": 31, "merge_cut_height": .2}
    assert set(sources.values()) == {"runtime_override"}
    _actual, _sources, equal_changes = programs.configure_production_modules(config, values)
    assert equal_changes == {}
    programs.configure_production_membership(config, values)
    assert config.consensus.min_module_size == 31


def test_workflow_records_actual_values_and_only_changes_engine_identity_for_changed_settings(tmp_path, monkeypatch):
    manifest = create_full_fixture(tmp_path / "project")
    context = ProjectService(application_root=tmp_path / "state").open_project(manifest)
    source = programs.CoSpliceWorkflow._config_asset(context).path
    payload = yaml.safe_load(source.read_text(encoding="utf-8"))
    payload["network"]["min_module_size"] = 37
    payload["network"]["wgcna_merge_cut_height"] = .12
    source.write_text(yaml.safe_dump(payload), encoding="utf-8")
    before = hashlib.sha256(source.read_bytes()).hexdigest()
    calls = []

    class FakeEngine:
        def __init__(self, config):
            self.config = config
            calls.append(config)

        def run_stage(self, stage, *, force):
            assert stage == "audit" and not force
            out = self.config.output_dir / "01_audit" / "input_audit.tsv"
            out.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame([{"check": "synthetic engine fixture", "status": "ok"}]).to_csv(out, sep="\t", index=False)

    monkeypatch.setattr(programs, "CoSplicePipeline", FakeEngine)
    monkeypatch.setattr(programs, "INTERPRETATION_STAGES", ("audit",))
    workflow = programs.CoSpliceWorkflow()
    job = SimpleNamespace(check_cancelled=lambda: None, report=lambda *_args: None)
    original_payload = {
        "project": context.snapshot.manifest_fingerprint,
        "inputs": [fingerprint_file(path) for path in workflow.input_paths(context, {})],
        "max_events": 0, "bootstrap_iterations": 50,
        "uncertainty_role": "diagnostic_only", "membership_policy": "backbone_confidence_only",
        "production_contract": "backbone-annotated-v9-all-sample-input",
        "network_input_policy": programs.network_input_policy(),
    }
    legacy_compatible_root = context.manifest.output_root / "engines" / "cosplice" / fingerprint_value(original_payload)[:24]

    for run_id, supplied in (
        ("inherit-config-run", {}),
        ("equal-values-run", {"min_module_size": 37, "merge_cut_height": .12}),
        ("changed-values-run", {"min_module_size": 31, "merge_cut_height": .2}),
    ):
        parameters = workflow.normalize_parameters({"scope": "interpretation", **supplied})
        workflow.run(context, job, run_id, parameters)

    assert calls[0].output_dir == calls[1].output_dir == legacy_compatible_root
    assert calls[2].output_dir != legacy_compatible_root
    artifacts = context.catalog.list_artifacts(context.manifest.project_id, analysis_id="programs.cosplice")
    for artifact in artifacts:
        expected = {"min_module_size": 31, "merge_cut_height": .2} if artifact.run_fingerprint == "changed-values-run" else {"min_module_size": 37, "merge_cut_height": .12}
        assert artifact.metadata["network_module_parameters"] == expected
        if artifact.path.name == "resolved_cosplice.yaml":
            assert yaml.safe_load(artifact.path.read_text(encoding="utf-8"))["network_module_parameters"] == expected
    assert len(artifacts) == 6
    assert hashlib.sha256(source.read_bytes()).hexdigest() == before
