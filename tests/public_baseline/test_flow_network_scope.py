"""Contrast networks use exact existing vectors; these tests never fit a network."""
import json
from pathlib import Path

import pandas as pd
import pytest

from katana_cosplice.config import NetworkConfig
from katana_next.flow.model import file_sha256, validate_plan
from katana_next.flow.network_inputs import (
    network_scope_identity, prepare_project_network, resolve_network_rmats_mode, validate_external_network_config,
    verify_network_engine_scope,
)
from katana_next.flow.planner import create_plan
from katana_next.infrastructure.manifest_store import ManifestStore


@pytest.fixture
def study(tmp_path):
    from katana_next.full_fixture import create_full_fixture
    return create_full_fixture(tmp_path / "study")


def test_two_comparisons_have_independent_fit_and_followups(study):
    _, manifest = ManifestStore().load(study)
    ids = [c.contrast_id for c in manifest.contrasts[:2]]
    plan = create_plan(study, {"contrast_ids": ids})
    validate_plan(plan)
    fits = [t for t in plan["tasks"] if t.get("adapter") == "network.fit"]
    assert len(fits) == 2 and fits[0]["id"] != fits[1]["id"]
    for fit, cid in zip(fits, ids):
        contrast = next(c for c in manifest.contrasts if c.contrast_id == cid)
        expected = set(contrast.numerator.sample_ids + contrast.denominator.sample_ids)
        assert fit["scope"]["contrast_ids"] == [cid]
        assert set(fit["parameters"]["sample_ids"]) == expected
        assert fit["parameters"]["n_samples"] == 10
        assert fit["parameters"]["sample_size_status"] == "minimum_met"
        children = [t for t in plan["tasks"] if t.get("adapter", "").startswith("network.") and t["dependencies"] == [fit["id"]]]
        assert {t["adapter"] for t in children} == {"network.bootstrap", "network.associate", "network.annotate"}
        assert all(t["parameters"] == fit["parameters"] for t in children)
    cohort = create_plan(study, {"contrast_ids": ids, "network": {"scope": "cohort"}})
    cohort_fit = next(t for t in cohort["tasks"] if t.get("adapter") == "network.fit")
    assert cohort_fit["id"] == "network-fit"
    assert cohort_fit["parameters"]["n_samples"] == 15
    assert cohort_fit["parameters"]["minimum_network_samples"] == 15
    assert cohort["plan_id"] != plan["plan_id"]
    assert NetworkConfig().minimum_network_samples == 15


def test_native_projection_contains_only_its_contrast_and_preserves_sources(study, tmp_path):
    from katana_cosplice.matrix import build_sample_matrices
    _, manifest = ManifestStore().load(study)
    originals = {a.path: file_sha256(a.path) for a in manifest.assets if a.kind == "rmats_result"}
    matrices = []
    for index, contrast in enumerate(manifest.contrasts[:2]):
        identity = network_scope_identity(manifest, {"network_scope": "contrast", "contrast_ids": [contrast.contrast_id]})
        config, proof, _ = prepare_project_network(study, tmp_path / f"projection-{index}",
            {**identity, "rmats_mode": "JC"}, file_sha256)
        mapping = json.loads(Path(proof["projection_mapping"]).read_text())
        assert {r["contrast_id"] for r in mapping["projection"]} == {contrast.contrast_id}
        assert len(mapping["projection"]) == 5
        expected = set(contrast.numerator.sample_ids + contrast.denominator.sample_ids)
        metadata = pd.read_csv(config.sample_metadata, sep="\t")
        assert set(metadata.sample_id) == expected
        matrix = build_sample_matrices(config)
        assert set(matrix.psi.columns) == expected
        assert config.network.minimum_network_samples == 6
        # Transport must reproduce each declared side's values in its source.
        source = next(a for a in manifest.assets if a.kind == "rmats_result"
                      and a.metadata["contrast_id"] == contrast.contrast_id and a.metadata["event_type"] == "SE")
        native = pd.read_csv(source.path, sep="\t")
        from katana_cosplice.matrix import _event_frame
        native_keys = _event_frame(native, "SE").event_key
        for side, arm in enumerate(("numerator", "denominator"), 1):
            observed = native[f"IncLevel{side}"].str.split(",", expand=True).astype(float)
            observed.columns = source.metadata[arm + "_sample_ids"]
            observed.index = native_keys
            pd.testing.assert_frame_equal(matrix.psi.loc[list(native_keys), list(observed.columns)],
                                          observed, check_names=False)
        source_rows = matrix.sources
        assert set(source_rows.source_file) == {r["projected"] for r in mapping["projection"]}
        matrices.append(matrix)
        assert len(native) > 0
    assert set(matrices[0].psi.columns) != set(matrices[1].psi.columns)
    assert all(file_sha256(path) == digest for path, digest in originals.items())


@pytest.mark.parametrize("mode", ["JC", "JCEC"])
def test_single_declared_mode_is_resolved_and_recorded_for_every_network_task(study, tmp_path, mode):
    _, manifest = ManifestStore().load(study)
    assets = tuple(a.model_copy(update={"metadata": {**a.metadata, "mode": mode}})
                   if a.kind == "rmats_result" else a for a in manifest.assets)
    manifest = manifest.model_copy(update={"assets": assets})
    path = ManifestStore().save(manifest, tmp_path / f"{mode}.project.json")
    cid = manifest.contrasts[0].contrast_id
    plan = create_plan(path, {"contrast_ids": [cid]})
    tasks = [t for t in plan["tasks"] if t.get("adapter", "").startswith("network.")]
    assert len(tasks) == 4
    assert all(t["parameters"]["rmats_mode"] == mode for t in tasks)
    config, _, _ = prepare_project_network(path, tmp_path / f"{mode}-projection",
        {"network_scope": "contrast", "contrast_ids": [cid]}, file_sha256)
    assert config.rmats_mode == mode
    assert len(list(config.rmats_root.glob(f"*/rmats_post/*.MATS.{mode}.txt"))) == 5


def test_ambiguous_mode_requires_explicit_choice_with_parameter_precedence(study, tmp_path):
    _, manifest = ManifestStore().load(study)
    duplicate = tuple(a.model_copy(update={"asset_id": a.asset_id + "::JCEC",
                                         "metadata": {**a.metadata, "mode": "JCEC"}})
                      for a in manifest.assets if a.kind == "rmats_result")
    manifest = manifest.model_copy(update={"assets": manifest.assets + duplicate})
    cid = manifest.contrasts[0].contrast_id
    scope = {"network_scope": "contrast", "contrast_ids": [cid]}
    with pytest.raises(ValueError, match="both JC and JCEC"):
        resolve_network_rmats_mode(manifest, scope)
    path = ManifestStore().save(manifest, tmp_path / "both.project.json")
    with pytest.raises(ValueError, match="both JC and JCEC"):
        create_plan(path, {"contrast_ids": [cid]})
    chosen = manifest.model_copy(update={"metadata": {**manifest.metadata, "rmats_mode": "JCEC"}})
    assert resolve_network_rmats_mode(chosen, scope) == "JCEC"
    assert resolve_network_rmats_mode(chosen, {**scope, "rmats_mode": "JC"}) == "JC"


def test_unavailable_or_missing_mode_fails_during_planning_and_readiness(study, tmp_path, monkeypatch):
    from katana_next.flow import stages
    _, manifest = ManifestStore().load(study)
    cid = manifest.contrasts[0].contrast_id
    scope = {"network_scope": "contrast", "contrast_ids": [cid], "rmats_mode": "JCEC"}
    with pytest.raises(ValueError, match="has no declared scoped"):
        create_plan(study, {"contrast_ids": [cid], "network": {"parameters": {"rmats_mode": "JCEC"}}})
    def forbidden(*args, **kwargs):
        raise AssertionError("A missing count mode must fail before native runtime probing")
    monkeypatch.setattr(stages, "_require_r_packages", forbidden)
    errors = stages.readiness({"adapter": "network.fit", "inputs": {"project": str(study)}, "parameters": scope})
    assert any("has no declared scoped" in error for error in errors)
    no_mode = manifest.model_copy(update={"assets": tuple(a.model_copy(update={"metadata": {
        key: value for key, value in a.metadata.items() if key != "mode"}})
        if a.kind == "rmats_result" else a for a in manifest.assets)})
    with pytest.raises(ValueError, match="must declare an explicit JC/JCEC"):
        resolve_network_rmats_mode(no_mode, {**scope, "rmats_mode": "JC"})


def test_six_samples_formal_and_five_rejected_without_changing_legacy(study):
    _, manifest = ManifestStore().load(study)
    contrast = manifest.contrasts[0]
    arms = {name: getattr(contrast, name).model_copy(update={"sample_ids": getattr(contrast, name).sample_ids[:3]})
            for name in ("numerator", "denominator")}
    six = contrast.model_copy(update=arms)
    scoped = manifest.model_copy(update={"contrasts": (six,)})
    parameters = {"network_scope": "contrast", "contrast_ids": [six.contrast_id]}
    identity = network_scope_identity(scoped, parameters)
    assert identity["n_samples"] == 6
    assert identity["sample_size_status"] == "minimum_met"
    five = six.model_copy(update={"numerator": six.numerator.model_copy(update={"sample_ids": six.numerator.sample_ids[:2]})})
    with pytest.raises(ValueError, match="at least 6"):
        network_scope_identity(scoped.model_copy(update={"contrasts": (five,)}), parameters)
    with pytest.raises(ValueError, match="sample_ids disagrees"):
        network_scope_identity(scoped, {**identity, "sample_ids": identity["sample_ids"][:-1]})
    assert NetworkConfig().minimum_network_samples == 15


def test_recorded_old_identity_is_readable_but_new_plans_have_neutral_status(study):
    _, manifest = ManifestStore().load(study)
    cid = manifest.contrasts[0].contrast_id
    identity = network_scope_identity(manifest, {"network_scope": "contrast", "contrast_ids": [cid]})
    legacy = {**identity, "sample_size_status": "small_sample_exploratory"}
    assert network_scope_identity(manifest, legacy) == legacy
    assert legacy["sample_size_status"] == "small_sample_exploratory"  # Original record unchanged.
    with pytest.raises(ValueError, match="sample_size_status disagrees"):
        network_scope_identity(manifest, {**legacy, "sample_size_status": "validated_stable"})
    plan = create_plan(study, {"contrast_ids": [cid], "network": {"parameters": legacy}})
    assert all(t["parameters"]["sample_size_status"] == "minimum_met"
               for t in plan["tasks"] if t.get("adapter", "").startswith("network."))


def test_external_config_cannot_hide_foreign_samples_or_source_vectors(study, tmp_path):
    _, manifest = ManifestStore().load(study)
    contrast = manifest.contrasts[0]
    identity = network_scope_identity(manifest, {"network_scope": "contrast", "contrast_ids": [contrast.contrast_id]})
    own, _, _ = prepare_project_network(study, tmp_path / "own", {**identity, "rmats_mode": "JC"}, file_sha256)
    validate_external_network_config(own, manifest, identity, file_sha256)
    cohort, _, _ = prepare_project_network(study, tmp_path / "all", {"rmats_mode": "JC"}, file_sha256)
    with pytest.raises(ValueError, match="exactly match"):
        validate_external_network_config(cohort, manifest, identity, file_sha256)
    native = next(own.rmats_root.glob("*/rmats_post/SE.MATS.JC.txt"))
    native.write_text(native.read_text() + "\n")
    with pytest.raises(ValueError, match="outside this exact contrast"):
        validate_external_network_config(own, manifest, identity, file_sha256)


def test_followup_rejects_other_fit_even_with_same_sample_count(study, tmp_path):
    _, manifest = ManifestStore().load(study)
    identities = [network_scope_identity(manifest, {"network_scope": "contrast", "contrast_ids": [c.contrast_id]})
                  for c in manifest.contrasts[:2]]
    engine = tmp_path / "fit"; (engine / "02_matrix").mkdir(parents=True)
    (engine / "network_scope.json").write_text(json.dumps(identities[0]))
    pd.DataFrame({"sample_id": identities[0]["sample_ids"], "biological_sample_id": identities[0]["biological_sample_ids"]}).to_csv(engine / "02_matrix/sample_metadata_resolved.tsv", sep="\t", index=False)
    verify_network_engine_scope(engine, identities[0])
    with pytest.raises(ValueError, match="scope does not match"):
        verify_network_engine_scope(engine, identities[1])
    pd.DataFrame({"sample_id": identities[1]["sample_ids"]}).to_csv(engine / "02_matrix/sample_metadata_resolved.tsv", sep="\t", index=False)
    with pytest.raises(ValueError, match="foreign samples"):
        verify_network_engine_scope(engine, identities[0])


def test_stage_wrong_fit_rejected_before_any_network_computation(study, tmp_path, monkeypatch):
    from katana_cosplice.pipeline import CoSplicePipeline
    from katana_next.flow import stages
    from katana_next.workflows.programs import write_production_configuration
    _, manifest = ManifestStore().load(study)
    identities = [network_scope_identity(manifest, {"network_scope": "contrast", "contrast_ids": [c.contrast_id]})
                  for c in manifest.contrasts[:2]]
    config, _, _ = prepare_project_network(study, tmp_path / "projected", {**identities[0], "rmats_mode": "JC"}, file_sha256)
    yaml = write_production_configuration(config, tmp_path / "resolved.yaml")
    source = tmp_path / "other-fit"; source.mkdir()
    (source / "network_scope.json").write_text(json.dumps(identities[1]))
    def forbidden(*args, **kwargs):
        raise AssertionError("An incompatible parent cannot enter numerical computation")
    monkeypatch.setattr(CoSplicePipeline, "__init__", forbidden)
    with pytest.raises(ValueError, match="scope does not match"):
        stages._network({"adapter": "network.bootstrap", "parameters": identities[0],
                         "scope": identities[0], "inputs": {"project": str(study), "config": str(yaml),
                         "engine_root": str(source)}}, tmp_path / "bootstrap")
    assert not (tmp_path / "bootstrap").exists()


def test_preflight_rejects_external_cohort_config_for_contrast(study, tmp_path, monkeypatch):
    from katana_next.flow import stages
    from katana_next.workflows.programs import write_production_configuration
    _, manifest = ManifestStore().load(study)
    identity = network_scope_identity(manifest, {"network_scope": "contrast", "contrast_ids": [manifest.contrasts[0].contrast_id]})
    config, _, _ = prepare_project_network(study, tmp_path / "cohort", {"rmats_mode": "JC"}, file_sha256)
    path = write_production_configuration(config, tmp_path / "cohort.yaml")
    monkeypatch.setattr(stages, "_require_r_packages", lambda *a, **k: None)
    errors = stages.readiness({"adapter": "network.fit", "inputs": {"project": str(study), "config": str(path)},
                              "parameters": identity, "scope": identity})
    assert any("exactly match" in error for error in errors)


def test_technical_copies_and_changed_biological_units_are_not_independent_samples(study, tmp_path):
    _, manifest = ManifestStore().load(study)
    parameters = {"network_scope": "contrast", "contrast_ids": [manifest.contrasts[0].contrast_id]}
    identity = network_scope_identity(manifest, parameters)
    samples = list(manifest.samples)
    selected = [i for i, sample in enumerate(samples) if sample.sample_id in identity["sample_ids"]]
    samples[selected[1]] = samples[selected[1]].model_copy(update={"biological_sample_id": samples[selected[0]].biological_sample_id})
    with pytest.raises(ValueError, match="independent biological samples"):
        network_scope_identity(manifest.model_copy(update={"samples": tuple(samples)}), parameters)
    config, _, _ = prepare_project_network(study, tmp_path / "own", {**identity, "rmats_mode": "JC"}, file_sha256)
    metadata = pd.read_csv(config.sample_metadata, sep="\t")
    metadata.loc[0, "biological_sample_id"] = "unrelated-animal"
    metadata.to_csv(config.sample_metadata, sep="\t", index=False)
    with pytest.raises(ValueError, match="biological units differ"):
        validate_external_network_config(config, manifest, identity, file_sha256)


def test_reversed_native_arms_and_sample_order_preserve_the_exact_scoped_matrix(study, tmp_path):
    """Native arm reversal and vector permutations are transport, not new fits."""
    from katana_cosplice.matrix import build_sample_matrices
    _, manifest = ManifestStore().load(study)
    cid = manifest.contrasts[0].contrast_id
    identity = network_scope_identity(manifest, {"network_scope": "contrast", "contrast_ids": [cid]})
    expected_config, _, _ = prepare_project_network(study, tmp_path / "expected",
        {**identity, "rmats_mode": "JC"}, file_sha256)
    expected = build_sample_matrices(expected_config)
    assets = []
    untouched = {}
    for asset in manifest.assets:
        if asset.kind != "rmats_result" or asset.metadata.get("contrast_id") != cid:
            assets.append(asset)
            continue
        untouched[asset.path] = file_sha256(asset.path)
        frame = pd.read_csv(asset.path, sep="\t")
        for prefix in ("IncLevel", "IJC_SAMPLE_", "SJC_SAMPLE_"):
            first, second = frame[prefix + "1"].copy(), frame[prefix + "2"].copy()
            frame[prefix + "1"] = second.map(lambda value: ",".join(str(value).split(",")[::-1]))
            frame[prefix + "2"] = first.map(lambda value: ",".join(str(value).split(",")[::-1]))
        frame["IncLevelDifference"] = -frame["IncLevelDifference"]
        path = tmp_path / (asset.metadata["event_type"] + ".reversed.JC.tsv")
        frame.to_csv(path, sep="\t", index=False)
        metadata = {**asset.metadata, "reverse_direction": True}
        for arm in ("numerator", "denominator"):
            metadata[arm + "_sample_ids"] = list(metadata[arm + "_sample_ids"])[::-1]
        assets.append(asset.model_copy(update={"path": path, "metadata": metadata}))
    reversed_project = ManifestStore().save(manifest.model_copy(update={"assets": tuple(assets)}),
                                          tmp_path / "reversed.project.json")
    config, proof, _ = prepare_project_network(reversed_project, tmp_path / "observed",
        {**identity, "rmats_mode": "JC"}, file_sha256)
    observed = build_sample_matrices(config)
    for field in ("psi", "ijc", "sjc"):
        pd.testing.assert_frame_equal(getattr(observed, field).sort_index(),
                                      getattr(expected, field).sort_index())
    mapping = json.loads(Path(proof["projection_mapping"]).read_text())
    assert len(mapping["projection"]) == 5
    assert all(row["vector_order_transport_changed"] and not row["statistics_refitted"]
               for row in mapping["projection"])
    assert untouched == {path: file_sha256(path) for path in untouched}


@pytest.mark.parametrize("scope", ["contrast", "cohort"])
def test_missing_biological_units_fail_before_fit_or_followup(study, tmp_path, scope):
    _, manifest = ManifestStore().load(study)
    identity = network_scope_identity(manifest, {"network_scope": scope, "contrast_ids": [manifest.contrasts[0].contrast_id]})
    config, _, _ = prepare_project_network(study, tmp_path / "prepared", {**identity, "rmats_mode": "JC"}, file_sha256)
    metadata = pd.read_csv(config.sample_metadata, sep="\t").drop(columns="biological_sample_id")
    metadata.to_csv(config.sample_metadata, sep="\t", index=False)
    with pytest.raises(ValueError, match="requires biological_sample_id"):
        validate_external_network_config(config, manifest, identity, file_sha256)
    engine = tmp_path / "fit"; (engine / "02_matrix").mkdir(parents=True)
    (engine / "network_scope.json").write_text(json.dumps(identity))
    metadata.to_csv(engine / "02_matrix/sample_metadata_resolved.tsv", sep="\t", index=False)
    with pytest.raises(ValueError, match="requires biological_sample_id"):
        verify_network_engine_scope(engine, identity)


@pytest.mark.parametrize("scope", ["contrast", "cohort"])
def test_native_plot_metadata_uses_actual_comparison_scope(study, tmp_path, monkeypatch, scope):
    from katana_next.flow import stages
    _, manifest = ManifestStore().load(study)
    identity = network_scope_identity(manifest, {"network_scope": scope, "contrast_ids": [manifest.contrasts[0].contrast_id]})
    def fixture_outputs(spec, out):
        target = out / "engine/04_network/observed"
        (target / "plots").mkdir(parents=True)
        # This is a ledger fixture; no renderer or scientific fitting is invoked.
        data = target / "event_linkage.npy"; data.write_bytes(b"ledger fixture")
        figure = target / "plots/observed_event_dendrogram.png"; figure.write_bytes(b"ledger fixture")
        return {"network_identity": identity}, [data, figure]
    monkeypatch.setattr(stages, "readiness", lambda spec: [])
    monkeypatch.setattr(stages, "_network", fixture_outputs)
    receipt = stages.execute_stage({"adapter": "network.fit", "parameters": identity}, tmp_path / "output")
    figure = next(row for row in receipt["outputs"] if row["path"].endswith(".png"))
    assert figure["metadata"]["network_identity"] == identity
    assert figure["metadata"]["interpretation_status"] == identity["sample_size_status"]
    assert figure["metadata"]["data_links"][0]["contrast_id"] == (identity["contrast_ids"][0] if scope == "contrast" else "all")
