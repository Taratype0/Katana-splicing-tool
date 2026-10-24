"""Observed memberships are usable without fabricated resampling evidence."""
import copy
import json

import numpy as np
import pandas as pd
import pytest

from katana_cosplice.config import ModuleInterpretationConfig
from katana_cosplice.module_interpretation import module_hub_tables
from katana_next.flow.planner import create_plan
from katana_next.flow.stages import _network_hubs


@pytest.fixture
def study(tmp_path):
    from katana_next.full_fixture import create_full_fixture
    return create_full_fixture(tmp_path / "study")


def test_new_plans_explicitly_disable_bootstrap_and_record_membership(study):
    settings = {"network": {"parameters": {"bootstrap_iterations": 50}}}
    unchanged = copy.deepcopy(settings)
    plan = create_plan(study, settings)
    assert settings == unchanged
    assert plan["configuration"]["network"]["bootstrap"] == {"enabled": False}
    assert plan["configuration"]["network"]["parameters"]["hub_metric"] == "abs_kME"
    tasks = [t for t in plan["tasks"] if t.get("adapter", "").startswith("network.")]
    assert {t["adapter"] for t in tasks} == {"network.fit", "network.associate", "network.annotate"}
    assert all(t["parameters"]["bootstrap_enabled"] is False for t in tasks)
    assert all(t["parameters"]["hub_metric"] == "abs_kME" for t in tasks)
    assert all("bootstrap" not in dependency for t in tasks for dependency in t["dependencies"])


def test_explicit_bootstrap_opt_in_has_a_real_task_and_distinct_plan(study):
    disabled = create_plan(study, {})
    enabled = create_plan(study, {"network": {"bootstrap": {"enabled": True},
        "parameters": {"bootstrap_iterations": 3}}})
    tasks = [t for t in enabled["tasks"] if t.get("adapter", "").startswith("network.")]
    assert any(t["adapter"] == "network.bootstrap" for t in tasks)
    assert all(t["parameters"]["bootstrap_enabled"] and t["parameters"]["bootstrap_iterations"] == 3 for t in tasks)
    assert enabled["plan_id"] != disabled["plan_id"]


@pytest.mark.parametrize("network", [
    {"bootstrap": {"enabled": "false"}}, {"bootstrap": {"enabled": 0}},
    {"bootstrap": "disabled"}, {"bootstrap": {"mode": "none"}},
    {"parameters": {"hub_metric": "confidence"}},
    {"parameters": {"bootstrap_enabled": True}},
    {"parameters": {"bootstrap_enabled": 0}},
    {"bootstrap": {"enabled": True}, "parameters": {"bootstrap_iterations": 1}},
    {"bootstrap": {"enabled": True}, "parameters": {"bootstrap_iterations": "50"}},
])
def test_unknown_or_ambiguous_policy_is_rejected(study, network):
    with pytest.raises(ValueError):
        create_plan(study, {"network": network})


def hub_fixture():
    samples = [f"s{i}" for i in range(6)]
    eigengene = np.arange(6, dtype=float)
    psi = pd.DataFrame([eigengene, -eigengene, [0, 1, 0, 1, 0, 1]],
        index=["z", "a", "weak"], columns=samples)
    assignments = pd.DataFrame({"event_key": ["z", "a", "weak"], "module": ["M01"] * 3,
        "consensus_membership_confidence": [0.1, 0.1, 1], "network_connectivity": [0.01, 0.01, 100]})
    activity = pd.DataFrame({"M01": eigengene}, index=samples)
    annotation = pd.DataFrame({"event_key": psi.index, "gene_id": ["GA", "GZ", "GW"],
        "gene_symbol": ["A", "Z", "Weak"], "event_type": ["SE"] * 3})
    return assignments, psi, activity, annotation


def test_absolute_kme_is_primary_signed_values_and_deterministic_ties_are_retained():
    assignments, psi, activity, annotation = hub_fixture()
    events, genes = module_hub_tables(assignments, psi, activity, annotation)
    assert list(events.event_key) == ["a", "z", "weak"]
    assert list(events.hub_event_rank) == [1, 2, 3]
    assert events.set_index("event_key").loc["a", "module_membership_kME"] < 0
    np.testing.assert_allclose(events.hub_event_score, events.module_membership_kME.abs())
    assert list(genes.gene_id) == ["GA", "GZ", "GW"]
    assert (events.hub_metric == "abs_kME").all()
    assert (genes.hub_metric == "abs_kME").all()
    shuffled = module_hub_tables(assignments.sample(frac=1, random_state=5), psi, activity, annotation)[0]
    pd.testing.assert_frame_equal(events.reset_index(drop=True), shuffled.reset_index(drop=True))
    # Frozen/standalone default is retained; the new flow opts into the new rule.
    assert ModuleInterpretationConfig().hub_metric == "legacy_composite"
    legacy, _ = module_hub_tables(assignments, psi, activity, annotation, hub_metric="legacy_composite")
    assert legacy.set_index("event_key").loc["weak", "hub_event_rank"] == 1


def test_hub_export_reuses_cached_fit_and_records_actual_source_hashes(tmp_path, monkeypatch):
    import katana_cosplice.network as network
    monkeypatch.setattr(network, "discover_module_labels", lambda *a, **k: pytest.fail("network refit"))
    assignments, psi, activity, annotation = hub_fixture()
    source = tmp_path / "04_network/observed"; source.mkdir(parents=True)
    (tmp_path / "02_matrix").mkdir()
    assignments.to_csv(source / "event_module_assignments.tsv", sep="\t", index=False)
    psi.to_csv(source / "selected_event_matrix.tsv", sep="\t")
    activity.to_csv(source / "module_program_scores.tsv", sep="\t")
    annotation.to_csv(tmp_path / "02_matrix/event_metadata.tsv", sep="\t", index=False)
    paths = _network_hubs({"parameters": {"hub_metric": "abs_kME"}}, tmp_path)
    proof = json.loads(open(paths["hub_parameters"], encoding="utf-8").read())
    assert proof["hub_metric"] == "abs_kME"
    assert proof["network_refitted"] is False and proof["bootstrap_performed"] is False
    assert all(len(row["sha256"]) == 64 for row in proof["sources"].values())
    assert pd.read_csv(paths["hub_events"], sep="\t").event_key.tolist() == ["a", "z", "weak"]


def test_grey_only_output_is_an_empty_typed_table_and_unknown_metric_is_rejected():
    assignments, psi, activity, annotation = hub_fixture()
    assignments["module"] = "grey"
    events, genes = module_hub_tables(assignments, psi, activity, annotation)
    assert events.empty and genes.empty
    assert "module_membership_kME" in events and "hub_gene_rank" in genes
    with pytest.raises(ValueError, match="hub_metric"):
        module_hub_tables(assignments, psi, activity, annotation, hub_metric="guess")


def test_disabled_fit_never_enters_edge_or_graph_resampling(monkeypatch):
    import katana_cosplice.network as network
    from katana_cosplice.config import NetworkConfig
    monkeypatch.setattr(network, "_edge_bootstrap_support_arrays",
                        lambda *a, **k: pytest.fail("unrequested resampling"))
    assignments, psi, activity, annotation = hub_fixture()
    settings = NetworkConfig(bootstrap_iterations=0, consensus_bootstrap_iterations=0)
    edges = pd.DataFrame({"event_a": ["a"], "event_b": ["z"]})
    assert np.isnan(network._edge_bootstrap_support(psi, edges, settings, None)).all()
    adjacency = np.array([0.75])
    unchanged, support = network._bootstrap_weighted_edges(
        psi, np.array([0]), np.array([1]), adjacency, settings, None)
    np.testing.assert_array_equal(unchanged, adjacency)
    assert np.isnan(support).all()


def test_observed_membership_is_never_relabelled_as_bootstrap_confidence():
    assignments, psi, activity, annotation = hub_fixture()
    assignments = assignments.drop(columns="consensus_membership_confidence")
    assignments["module_membership"] = [.99, -.99, .5]
    events, genes = module_hub_tables(assignments, psi, activity, annotation)
    assert events.consensus_confidence.isna().all()
