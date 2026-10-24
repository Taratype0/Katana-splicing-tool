"""Network presentation shows actual edges and complete gene identities."""
from hashlib import sha256
import json
import os
import re

import pandas as pd
import pytest

from katana_next.reports.network_views import (
    _module_display_graph, _module_gene_records, _write_browser, refresh_completed_network_browser,
)


def _payload(path):
    return json.loads(re.search(r'<script id="data" type="application/json">(.*?)</script>',
                               path.read_text(encoding="utf-8"), re.S).group(1))


def _browser(root):
    (root / "browser").mkdir(parents=True)
    (root / "tables").mkdir()
    members = pd.DataFrame({
        "event_key": ["e1", "e2", "e3", "e4", "e5"],
        "module": ["M01", "M01", "M02", "M02", "grey"],
        "gene_id": ["ENSG1", "ENSG1", "ENSG2", "ENSG3", ""],
        "gene_symbol": ["G1", "G1", "G2", "G3", ""],
        "kME": [.9, -.8, .7, .6, None],
    })
    scores = pd.DataFrame({"M01": [.2, .8], "M02": [.8, .2]}, index=["A", "B"])
    matrix = pd.DataFrame([[.1, .9], [.2, .8], [.3, .7], [.4, .6], [.5, .5]], columns=["A", "B"])
    samples = pd.DataFrame({"sample_id": ["A", "B"], "condition": ["control", "treated"]})
    edges_path = root / "tables/native_network_edges.tsv"
    pd.DataFrame({"event_a": ["e1"], "event_b": ["e2"],
                  "adjacency": [.8751345612345678]}).to_csv(edges_path, sep="\t", index=False)
    _write_browser(root, members, matrix, scores, samples, pd.DataFrame(), edges_path, bootstrap_enabled=False)
    return root / "browser/ALL_MODULES.html", members


def test_complete_gene_identity_groups_events_without_top_n(tmp_path):
    _, members = _browser(tmp_path / "source")
    extra = pd.DataFrame({"event_key": [f"more{i}" for i in range(125)], "module": "M01",
                          "gene_id": [f"ENSG_extra_{i}" for i in range(125)],
                          "gene_symbol": [f"G_extra_{i}" for i in range(125)]})
    frame = pd.concat([members, extra], ignore_index=True)
    original = frame.copy(deep=True)
    genes = _module_gene_records(frame)
    assert len(genes) == 128
    first = genes.set_index("gene_id").loc["ENSG1"]
    assert first["gene_symbol"] == "G1" and first["event_count"] == 2
    assert first["event_keys"] == "e1;e2"
    assert set(genes.loc[genes.module.eq("M01"), "gene_id"]) == {"ENSG1", *extra.gene_id}
    assert not genes.gene_symbol.eq("").any()  # Unmapped events are not invented genes.
    pd.testing.assert_frame_equal(frame, original)


def test_graph_and_complete_genes_are_visible_before_collapsed_matrix(tmp_path):
    browser, _ = _browser(tmp_path / "source")
    page = browser.read_text(encoding="utf-8")
    assert page.index('id="module-network"') < page.index('id="module-genes"') < page.index('id="psi"')
    assert '<details><summary>PSI, all event members' in page
    assert 'id="genes-download"' in page and 'id="gene-next"' in page
    assert "No replacement links were calculated" in page
    assert "sel.value=assigned[0]" in page  # It no longer opens on the blank 'all' graph.
    payload = _payload(browser)
    assert payload["edges"] == [[0, 1, .8751345612345678]]
    assert len(payload["members"]) == 5
    genes = pd.read_csv(browser.parent.parent / "tables/all_module_genes.tsv", sep="\t")
    assert len(genes) == 3


def test_completed_browser_refresh_breaks_hardlinks_and_retains_native_payload(tmp_path):
    source, _ = _browser(tmp_path / "source")
    before = sha256(source.read_bytes()).hexdigest()
    target = tmp_path / "independent-copy/browser/ALL_MODULES.html"
    target.parent.mkdir(parents=True)
    os.link(source, target)
    assert os.stat(source).st_ino == os.stat(target).st_ino
    receipt = refresh_completed_network_browser(source, target)
    assert receipt["passed"] and receipt["original_payload_unchanged"]
    assert receipt["network_discovery_runs"] == receipt["association_tests"] == receipt["enrichment_tests"] == 0
    assert receipt["events"] == 5 and receipt["native_exported_edges"] == 1
    assert os.stat(source).st_ino != os.stat(target).st_ino
    assert sha256(source.read_bytes()).hexdigest() == before
    assert _payload(source) == _payload(target)
    assert (target.parent.parent / "tables/all_module_genes.tsv").is_file()
    assert "kME is not resampling confidence" in target.read_text(encoding="utf-8")


def test_refresh_refuses_original_path_or_unknown_payload_fields(tmp_path):
    source, _ = _browser(tmp_path / "source")
    with pytest.raises(ValueError, match="new report copy"):
        refresh_completed_network_browser(source, source)
    payload = _payload(source)
    payload["unsupported_evidence"] = {"must_not_be_lost": True}
    source.write_text('<script id="data" type="application/json">' + json.dumps(payload) + '</script>', encoding="utf-8")
    with pytest.raises(ValueError, match="silently discard"):
        refresh_completed_network_browser(source, tmp_path / "copy/browser/ALL_MODULES.html")


def _native_graph_fixture(root):
    import numpy as np
    from katana_cosplice.config import NetworkConfig
    from katana_cosplice.network import correlation_matrix, signed_soft_adjacency
    root.mkdir()
    rng = np.random.default_rng(71)
    members = pd.DataFrame({"event_key": [f"event_{i:03}" for i in range(110)], "module": "M01",
        "gene_id": [f"GENEID_{i}" for i in range(110)], "gene_symbol": [f"G{i}" for i in range(110)],
        "abs_kME": np.linspace(.99, .2, 110), "hub_rank_abs_kME": np.arange(1, 111)})
    matrix = pd.DataFrame(rng.uniform(.1, .9, (110, 6)), columns=list("ABCDEF"))
    native = matrix.copy(); native.index = members.event_key
    native.to_csv(root / "selected_event_matrix.tsv", sep="\t")
    members[["event_key", "module"]].to_csv(root / "event_module_assignments.tsv", sep="\t", index=False)
    settings = {"engine": "tom_hierarchical", "method": "bicor", "signed": True,
        "soft_power": 6, "edge_threshold": .1, "uncertainty": {"network_mode": "diagnostic_only"},
        "gene_prior": {"enabled": False}}
    (root / "parameters.json").write_text(json.dumps({"settings": settings}), encoding="utf-8")
    pd.DataFrame({"uncertainty_weighting": [False], "posterior_ensemble_correlation": [False],
        "gene_mantel_prior_enabled": [False], "engine": ["tom_hierarchical"], "method": ["bicor"],
        "soft_power": [6], "signed_network": [True]}).to_csv(root / "network_diagnostics.tsv", sep="\t", index=False)
    adjacency = signed_soft_adjacency(correlation_matrix(native, "bicor"), NetworkConfig(method="bicor", signed=True, soft_power=6))
    i, j = np.unravel_index(np.argmax(adjacency[:100, :100]), (100, 100))
    if i > j:
        i, j = j, i
    pd.DataFrame({"event_a": [members.iloc[i].event_key], "event_b": [members.iloc[j].event_key],
        "adjacency": [float(adjacency[i, j])]}).to_csv(root / "network_edges.tsv", sep="\t", index=False)
    return members, matrix, adjacency


def test_module_graph_display_uses_exact_native_kernel_without_changing_membership(tmp_path):
    import numpy as np
    root = tmp_path / "native"
    members, matrix, native_adjacency = _native_graph_fixture(root)
    source_bytes = {p.name: p.read_bytes() for p in root.iterdir()}
    graph, edges = _module_display_graph(members, matrix, root / "parameters.json", max_edges=10)
    selected = graph["modules"]["M01"]
    assert selected["total_module_events"] == 110
    assert selected["node_indices"] == list(range(100))
    assert len(selected["edges"]) == 10 and len(edges) == 10
    assert len(members) == 110  # Drawing limits do not filter scientific membership.
    assert graph["methods"]["matched_native_edges"] == 1
    assert graph["methods"]["maximum_native_adjacency_difference"] <= 1e-7
    for i, j, value in selected["edges"]:
        assert value == float(native_adjacency[i, j])
        assert value >= .1
    assert graph["methods"]["full_topology_claimed"] is False
    assert graph["methods"]["TOM_runs"] == graph["methods"]["network_discovery_runs"] == 0
    assert source_bytes == {p.name: p.read_bytes() for p in root.iterdir()}


@pytest.mark.parametrize("mutation", ["parameters", "native_scope", "weighted"])
def test_module_graph_display_refuses_substituted_parameters_or_another_fit(tmp_path, mutation):
    root = tmp_path / "native"
    members, matrix, _ = _native_graph_fixture(root)
    if mutation == "parameters":
        value = json.loads((root / "parameters.json").read_text())
        value["settings"]["soft_power"] = 9
        (root / "parameters.json").write_text(json.dumps(value))
    elif mutation == "native_scope":
        members.loc[0, "module"] = "another_module"
    else:
        frame = pd.read_csv(root / "network_diagnostics.tsv", sep="\t")
        frame["uncertainty_weighting"] = True
        frame.to_csv(root / "network_diagnostics.tsv", sep="\t", index=False)
    with pytest.raises(ValueError):
        _module_display_graph(members, matrix, root / "parameters.json")
