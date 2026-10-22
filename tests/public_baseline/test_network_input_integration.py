"""Exercise real network input-QC orchestration without fitting any R network.

Only synthetic matrices in pytest's temporary project are read.  The assertions
cover scientific eligibility, diagnostic-only uncertainty, cache provenance and
protection of historical output; no experimental data or frozen app is touched.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from katana_cosplice.config import CoSpliceConfig, StableASConfig
from katana_cosplice.pipeline import CoSplicePipeline


POLICY_ID = "network-all-samples-counts-v1"
ELIGIBLE = ["tiny_change", "within_condition_change"]


def _hashes(paths: list[Path]) -> dict[str, str]:
    return {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}


def _project(tmp_path: Path, *, all_excluded: bool = False):
    samples = [f"S{i}" for i in range(6)]
    psi = pd.DataFrame(
        [
            [.4999, .5, .5001, .4999, .5, .5001],
            [.2, .4, .6, .2, .4, .6],
            [.1, .2, .3, .6, .7, .8],
            [.1, .2, np.nan, .6, .7, .8],
            [.1, .2, .3, .6, .7, .8],
            [.1, .2, .3, .6, .7, .8],
            [.1, .2, 1.01, .6, .7, .8],
            [.5] * 6,
        ],
        index=pd.Index(
            ELIGIBLE + [
                "one_low_sample", "missing_psi", "missing_inclusion",
                "nonfinite_count", "invalid_psi", "constant",
            ],
            name="event_key",
        ),
        columns=samples,
    )
    ijc = pd.DataFrame(5., index=psi.index, columns=psi.columns)
    sjc = pd.DataFrame(5., index=psi.index, columns=psi.columns)
    sjc.loc["one_low_sample", "S5"] = 4.
    ijc.loc["missing_inclusion", "S1"] = np.nan
    ijc.loc["nonfinite_count", "S2"] = np.inf
    if all_excluded:
        ijc.loc[:, :] = 0.
        sjc.loc[:, :] = 0.
    metadata = pd.DataFrame({"sample_id": samples, "condition": ["A"] * 3 + ["B"] * 3})
    output = tmp_path / "output"
    matrix_root = output / "02_matrix"
    matrix_root.mkdir(parents=True)
    paths = []
    for frame, name in (
        (psi, "sample_level_psi_matrix.tsv"),
        (ijc, "sample_level_ijc_matrix.tsv"),
        (sjc, "sample_level_sjc_matrix.tsv"),
    ):
        path = matrix_root / name
        frame.to_csv(path, sep="\t")
        paths.append(path)
    metadata_path = matrix_root / "sample_metadata_resolved.tsv"
    metadata.to_csv(metadata_path, sep="\t", index=False)
    paths.append(metadata_path)
    config = CoSpliceConfig("read-support-fixture", None, metadata_path, output)
    # Deliberately hostile historical settings must not reinstate the obsolete
    # preselection; even an exact ten-count event with no condition response passes.
    config.stable_as = StableASConfig(
        preset="strict", min_total_junction_reads=1000.,
        min_psi_variance=1., min_repeatability=.99,
        min_condition_mean_psi_range=1., min_sample_psi_range=1.,
    )
    config.uncertainty.enabled = True
    config.uncertainty.filter_low_reliability = True
    config.uncertainty.min_network_reliability = 1.
    config.uncertainty.network_mode = "diagnostic_only"
    return CoSplicePipeline(config), psi, paths


def _read_outputs(outputs):
    metrics = pd.read_csv(outputs["metrics"], sep="\t", index_col=0)
    matrix = pd.read_csv(outputs["matrix"], sep="\t", index_col=0)
    params = json.loads(outputs["parameters"].read_text(encoding="utf-8"))
    return metrics, matrix, params


def test_qc_uses_every_sample_and_uncertainty_cannot_remove_eligible_events(tmp_path):
    pipeline, original, paths = _project(tmp_path)
    before = _hashes(paths)
    outputs = pipeline.stable_qc()
    metrics, matrix, params = _read_outputs(outputs)

    assert matrix.index.tolist() == ELIGIBLE
    pd.testing.assert_frame_equal(matrix, original.loc[ELIGIBLE])
    assert set(metrics.index[metrics["stable_as"]]) == set(ELIGIBLE)
    assert not metrics.loc[metrics.index.difference(ELIGIBLE), "stable_as"].any()
    assert metrics.loc["tiny_change", "psi_variance"] < .0025
    assert metrics.loc["within_condition_change", "replicate_repeatability"] < .30
    assert params["network_input_policy"]["policy_id"] == POLICY_ID
    assert params["n_input_events"] == len(original)
    assert params["n_stable_events"] == len(ELIGIBLE)
    assert params["ignored_legacy_qc_settings"]["min_total_junction_reads"] == 1000.
    assert params["ignored_legacy_qc_settings"]["min_repeatability"] == .99

    diagnostics = pd.read_csv(outputs["uncertainty_events"], sep="\t", index_col=0)
    sample_diagnostics = pd.read_csv(outputs["uncertainty_samples"], sep="\t")
    assert diagnostics.index.tolist() == ELIGIBLE
    assert set(sample_diagnostics["event_key"]) == set(ELIGIBLE)
    assert len(sample_diagnostics) == 6 * len(ELIGIBLE)
    assert diagnostics["network_reliability"].lt(1.).all()
    assert metrics.loc[metrics.index.difference(ELIGIBLE), "network_reliability"].isna().all()
    assert _hashes(paths) == before


def test_all_excluded_writes_complete_empty_diagnostics_instead_of_failure(tmp_path):
    pipeline, original, paths = _project(tmp_path, all_excluded=True)
    before = _hashes(paths)
    outputs = pipeline.stable_qc()
    metrics, matrix, params = _read_outputs(outputs)

    assert all(path.is_file() and path.stat().st_size for path in outputs.values())
    assert matrix.empty and matrix.columns.tolist() == original.columns.tolist()
    assert len(metrics) == len(original) and not metrics["stable_as"].any()
    assert params["n_stable_events"] == 0
    for key in ("posterior_matrix", "posterior_alpha", "posterior_beta", "uncertainty_events"):
        assert pd.read_csv(outputs[key], sep="\t", index_col=0).empty
    sample_diagnostics = pd.read_csv(outputs["uncertainty_samples"], sep="\t")
    assert sample_diagnostics.empty
    assert {"event_key", "sample_id", "posterior_psi"}.issubset(sample_diagnostics.columns)
    assert _hashes(paths) == before


def test_current_policy_cache_reuses_artifacts_without_rewriting_matrices(tmp_path):
    pipeline, _original, inputs = _project(tmp_path)
    first = pipeline.stable_qc()
    before = _hashes(inputs + list(first.values()))
    times = {key: value.stat().st_mtime_ns for key, value in first.items()}

    second = pipeline.stable_qc()

    assert second == first
    assert pipeline.status()["stages"]["stable_qc"]["status"] == "cache_hit"
    assert {key: value.stat().st_mtime_ns for key, value in second.items()} == times
    assert _hashes(inputs + list(second.values())) == before


@pytest.mark.parametrize("force", [False, True])
def test_legacy_policy_output_is_never_overwritten_even_when_forced(tmp_path, force):
    pipeline, _original, inputs = _project(tmp_path)
    out_dir = pipeline.output / "03_stable_as_qc"
    out_dir.mkdir()
    params = out_dir / "parameters.json"
    params.write_text('{"thresholds":{"min_total_junction_reads":10},"historical":true}', encoding="utf-8")
    old_table = out_dir / "stable_psi_matrix.tsv"
    old_table.write_text("event_key\told_sample\nhistorical\t0.25\n", encoding="utf-8")
    before = _hashes(inputs + [params, old_table])

    with pytest.raises(ValueError, match="(?i)(policy|historical|legacy|new output)"):
        pipeline.stable_qc(force=force)

    assert _hashes(inputs + [params, old_table]) == before
    assert set(path.name for path in out_dir.iterdir()) == {params.name, old_table.name}


def test_old_parameters_cannot_reuse_a_matching_current_state_fingerprint(tmp_path):
    pipeline, _original, inputs = _project(tmp_path)
    outputs = pipeline.stable_qc()
    _metrics, _matrix, params = _read_outputs(outputs)
    del params["network_input_policy"]
    outputs["parameters"].write_text(json.dumps(params), encoding="utf-8")
    before = _hashes(inputs + list(outputs.values()))

    with pytest.raises(ValueError, match="(?i)(policy|historical|legacy|new output)"):
        pipeline.stable_qc()

    assert _hashes(inputs + list(outputs.values())) == before


def test_production_config_declares_policy_without_changing_module_defaults(tmp_path):
    import yaml
    from katana_cosplice.config import NetworkConfig, load_config
    from katana_next.workflows.programs import write_production_configuration

    assert NetworkConfig().min_module_size == 20
    assert NetworkConfig().wgcna_merge_cut_height == .25
    pipeline, _, _ = _project(tmp_path)
    pipeline.config.network.min_module_size = 37
    pipeline.config.network.wgcna_merge_cut_height = .12
    target = write_production_configuration(pipeline.config, tmp_path / 'resolved.yaml')
    payload = yaml.safe_load(target.read_text(encoding='utf-8'))
    assert payload['network_input_policy']['policy_id'] == POLICY_ID
    assert 'ignored historical' in payload['legacy_stable_as_settings_role']
    restored = load_config(target)
    assert restored.network.min_module_size == 37
    assert restored.network.wgcna_merge_cut_height == .12


@pytest.mark.parametrize('stage', ['audit', 'matrix', 'network_observed'])
def test_historical_qc_cannot_be_overwritten_before_another_stage_starts(tmp_path, stage):
    pipeline, _, paths = _project(tmp_path)
    qc_root = pipeline.output / '03_stable_as_qc'
    qc_root.mkdir()
    old = qc_root / 'parameters.json'
    old.write_text('{"thresholds":{"min_read_supported_fraction":0.7}}', encoding='utf-8')
    protected = [*paths, old]
    before = _hashes(protected)
    with pytest.raises(ValueError, match='Historical network QC output'):
        pipeline._mark(stage, 'running', 'must not update historical state')
    assert _hashes(protected) == before
    assert not pipeline.state_path.exists()
