import numpy as np
import pandas as pd
import pytest
from scipy.stats import pearsonr
from katana_next.science.module_dashboard import build_module_dashboard, pairwise_pearson
from katana_next.science.psi_network import prepare_psi_network


def fixture(n=40):
    rng = np.random.default_rng(31)
    samples = pd.DataFrame({"sample_id": [f"s{i}" for i in range(n)], "biological_sample_id": [f"d{i}" for i in range(n)], "condition": ["A"]*(n//2)+["B"]*(n-n//2)})
    matrix = pd.DataFrame(rng.uniform(0,1,(45,n)), index=[f"event{i}" for i in range(45)], columns=samples.sample_id)
    matrix.iloc[1] = matrix.iloc[0]
    assignments = pd.DataFrame({"event_key": matrix.index, "consensus_module": "M01"})
    events = pd.DataFrame({"event_key": matrix.index, "gene_symbol": ["Shared"]*2+[f"G{i}" for i in range(43)]})
    return matrix, assignments, samples, events


def test_pairwise_missing_values_match_scipy_and_do_not_impute():
    x = np.random.default_rng(4).uniform(size=(4,35))
    x[0,::3] = np.nan
    x[1,::4] = np.nan
    actual, n = pairwise_pearson(x)
    valid = np.isfinite(x[0]) & np.isfinite(x[1])
    assert n[0,1] == valid.sum()
    assert actual[0,1] == pytest.approx(pearsonr(x[0,valid],x[1,valid]).statistic)
    x[2,:30] = np.nan
    assert np.isnan(pairwise_pearson(x)[0][0,2])


def test_event_nodes_remain_distinct_and_same_gene_is_display_only():
    matrix, labels, samples, events = fixture()
    a = build_module_dashboard(matrix,labels,samples,"M01",metadata=events,condition="A")
    b = build_module_dashboard(matrix,labels,samples,"M01",metadata=events,condition="A",hide_same_gene=True)
    assert a["node_unit"] == "splice_event" and a["total_events"] == 45
    assert a["same_gene_edges"] == 1 and not any(e["same_gene"] for e in b["edges"])
    assert a["hub_events"] == b["hub_events"] and a["total_edges"] == b["total_edges"]
    assert len(a["event_options"]) == 45
    assert a["heatmap"]["scale"].startswith("Row z-score")
    assert matrix.index.is_unique


def test_underpowered_condition_shows_measurements_not_invented_edges():
    matrix, labels, samples, events = fixture(24)
    d = build_module_dashboard(matrix,labels,samples,"M01",metadata=events,condition="A")
    assert d["status"] == "insufficient_biological_samples"
    assert d["edges"] == [] and len(d["event_profile"]) == 12
    assert len(d["heatmap"]["events"]) == 45


def test_technical_repeats_do_not_make_an_eligible_network():
    matrix, labels, samples, events = fixture()
    samples["biological_sample_id"] = [f"d{i//2}" for i in range(40)]
    d = build_module_dashboard(matrix,labels,samples,"M01",metadata=events,condition="A")
    assert d["biological_samples"] == 10 and d["edges"] == []


def test_cached_pooled_edges_do_not_turn_into_within_condition_edges():
    matrix, labels, samples, events = fixture(24)
    edges = pd.DataFrame({"event_a":["event0"],"event_b":["event1"],"adjacency":[.99],"correlation":[.999]})
    pooled = build_module_dashboard(matrix,labels,samples,"M01",metadata=events,edges=edges)
    group = build_module_dashboard(matrix,labels,samples,"M01",metadata=events,edges=edges,condition="A")
    assert pooled["total_edges"] == 1 and "Cached discovery" in pooled["method"]
    assert group["total_edges"] == 0


def test_unknown_condition_event_and_duplicate_key_rejected():
    matrix, labels, samples, events = fixture()
    with pytest.raises(ValueError,match="Unknown condition"):
        build_module_dashboard(matrix,labels,samples,"M01",condition="missing")
    with pytest.raises(ValueError,match="outside"):
        build_module_dashboard(matrix,labels,samples,"M01",selected_event="absent")
    labels.loc[1,"event_key"] = labels.event_key.iloc[0]
    with pytest.raises(ValueError,match="unique"):
        build_module_dashboard(matrix,labels,samples,"M01")


def test_psi_only_contract_never_fabricates_read_support():
    matrix, labels, samples, events = fixture()
    matrix.iloc[4,0] = np.nan
    observed, complete, qc, means = prepare_psi_network(matrix,samples)
    assert qc.read_support_status.eq("not_assessed_PSI_only").all()
    assert qc.repeat_support_status.eq("not_assessed_PSI_only").all()
    assert np.isnan(observed.iloc[4,0]) and not complete.isna().any().any()
    assert complete.iloc[4,0] == pytest.approx(means.iloc[4])
