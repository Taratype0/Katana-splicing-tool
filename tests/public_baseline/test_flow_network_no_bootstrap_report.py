"""Explicit no-bootstrap completion preserves observed products and kME semantics."""
import copy
import json
from pathlib import Path
import shutil

import pandas as pd
import pytest

from katana_next.flow.model import seal_plan
from katana_next.flow.network_acceptance import network_report_gaps
from katana_next.flow.presentation import (build_presentation_assets, network_bootstrap_enabled,
                                            validate_network_view_proof)
from katana_next.flow.runtime import canonical_hash, file_hash, write_json
from katana_next.reports.engine import ASSET_ROOT, Builder
from katana_next.reports.sources import load_source
from test_flow_network_acceptance import _relative
from test_flow_presentation import (_two_contrast_network_collection, _fixture_collection,
                                    _seal_fixture_collection)


def _disabled_collection(tmp_path, *, declared=True, retain_bootstrap=False):
    original = tmp_path / "original"
    original.mkdir()
    old, receipts, _, sources = _two_contrast_network_collection(original)
    # Negative membership must rank by magnitude; equal magnitudes use the
    # unchanged event identity, rather than signed value or connectivity.
    for source in sources:
        path = source / "fit-engine/04_network/observed/event_module_assignments.tsv"
        frame = pd.read_csv(path, sep="\t")
        # Discovery uses its recorded bicor membership; the independently
        # completed Pearson Hub family below must never be replaced by it.
        frame.loc[frame.event_key.eq("e2"), "module_membership"] = .2
        frame.loc[frame.event_key.eq("e4"), "module_membership"] = .3
        frame.to_csv(path, sep="\t", index=False)
    tasks = [task for task in old["tasks"] if retain_bootstrap or task["adapter"] != "network.bootstrap"]
    products = {task["id"]: ([Path(row["path"]) for row in receipts[task["id"]][1]["outputs"]],
                             receipts[task["id"]][1]["execution_evidence"]) for task in tasks}
    for number, source in enumerate(sources):
        # The no-bootstrap completed fixture has the same saved drawing inputs
        # as a fresh native collection. Bind their actual bytes to its fit
        # receipt, rather than bypassing the production completion gate.
        from katana_cosplice.config import NetworkConfig
        from katana_cosplice.network import correlation_matrix, signed_soft_adjacency
        native = source / "fit-engine/04_network/observed"
        settings = {"engine": "tom_hierarchical", "method": "bicor", "signed": True,
                    "soft_power": 6, "edge_threshold": .1,
                    "uncertainty": {"network_mode": "diagnostic_only"},
                    "gene_prior": {"enabled": False}}
        parameters = native / "parameters.json"
        diagnostics = native / "network_diagnostics.tsv"
        write_json(parameters, {"settings": settings, "synthetic_lifecycle_fixture": True})
        pd.DataFrame([{"engine": "tom_hierarchical", "method": "bicor", "soft_power": 6,
                       "signed_network": True, "uncertainty_weighting": False,
                       "posterior_ensemble_correlation": False,
                       "gene_mantel_prior_enabled": False}]).to_csv(diagnostics, sep="\t", index=False)
        psi = pd.read_csv(native / "selected_event_matrix.tsv", sep="\t", index_col=0,
                          float_precision="round_trip")
        adjacency = signed_soft_adjacency(correlation_matrix(psi, "bicor"),
                                         NetworkConfig(method="bicor", signed=True, soft_power=6))
        # The legacy bootstrap fixture's arbitrary .8 edge is not a native
        # bicor/power-6 output; this fixture records the actual unchanged kernel.
        edges = [{"event_a": str(psi.index[i]), "event_b": str(psi.index[j]),
                  "adjacency": float(adjacency[i, j])}
                 for i in range(len(psi)) for j in range(i + 1, len(psi))
                 if adjacency[i, j] >= .1]
        pd.DataFrame(edges).to_csv(native / "network_edges.tsv", sep="\t", index=False)
        products[f"fit-{number}"][0].extend([parameters, diagnostics])
        directory = source / "annotation-engine"
        events = pd.DataFrame({"event_key": ["e1", "e2", "e3", "e4"], "module": ["M01", "M01", "M02", "M02"],
            "module_membership_kME": [.9, -.95, .8, -.8], "hub_event_score": [.9, .95, .8, .8],
            "hub_event_rank": [2, 1, 1, 2], "hub_metric": ["abs_kME"] * 4})
        paths = {"hub_events": directory / "hub-events.tsv", "hub_genes": directory / "hub-genes.tsv",
                 "hub_parameters": directory / "hub_parameters.json"}
        events.to_csv(paths["hub_events"], sep="\t", index=False)
        events.to_csv(paths["hub_genes"], sep="\t", index=False)
        engine = source / "fit-engine"
        files = {"assignments": engine / "04_network/observed/event_module_assignments.tsv",
                 "psi": engine / "04_network/observed/selected_event_matrix.tsv",
                 "activity": engine / "04_network/observed/module_program_scores.tsv",
                 "annotation": engine / "02_matrix/event_metadata.tsv"}
        write_json(paths["hub_parameters"], {"schema": "katana.network.hub_membership/v1", "hub_metric": "abs_kME",
            "signed_kME_retained": True, "network_refitted": False,
            "sources": {key: {"path": str(path.resolve()), "sha256": file_hash(path)} for key, path in files.items()}})
        outputs, evidence = products[f"annotate-{number}"]
        outputs.extend(paths.values())
        evidence["precomputed_views"].update({key: str(path) for key, path in paths.items()})
    project = json.loads((original / "collection/project.json").read_text())
    root = tmp_path / "collection"
    plan, receipts, transport = _fixture_collection(root, project, tasks, products)
    if declared:
        plan = seal_plan({**plan, "configuration": {"network": {"bootstrap": {"enabled": False},
                           "parameters": {"hub_metric": "abs_kME"}}}})
        write_json(root / "plan.json", plan)
        for path, value in receipts.values():
            value["plan_id"] = plan["plan_id"]
            value["receipt_sha256"] = canonical_hash({k: v for k, v in value.items() if k != "receipt_sha256"})
            write_json(path, value)
    return root, plan, receipts, transport, sources


def test_no_bootstrap_needs_explicit_sealed_policy_and_cannot_hide_legacy_missing_stage(tmp_path):
    root, plan, receipts, transport, _ = _disabled_collection(tmp_path, declared=False)
    assert network_bootstrap_enabled(plan) is True
    with pytest.raises(ValueError, match="all four completed network stages"):
        build_presentation_assets(root, plan, receipts, transport)
    with pytest.raises(ValueError, match="must be boolean"):
        network_bootstrap_enabled({"configuration": {"network": {"bootstrap": {"enabled": "false"}}}})


def test_disabled_plan_rejects_fake_completed_bootstrap_stage(tmp_path):
    root, plan, receipts, transport, _ = _disabled_collection(tmp_path, retain_bootstrap=True)
    with pytest.raises(ValueError, match="every selected completed network stage"):
        build_presentation_assets(root, plan, receipts, transport)


def test_disabled_report_cannot_substitute_discovery_membership_for_missing_hub_family(tmp_path):
    root, plan, receipts, transport, _ = _disabled_collection(tmp_path)
    for task in plan["tasks"]:
        if task["adapter"] == "network.annotate":
            receipts[task["id"]][1]["execution_evidence"]["precomputed_views"].pop("hub_events")
    with pytest.raises(ValueError, match="complete already computed contrast and ORA families"):
        build_presentation_assets(root, plan, receipts, transport)


@pytest.fixture
def disabled_report(tmp_path, monkeypatch):
    root, plan, receipts, transport, sources = _disabled_collection(tmp_path)
    before = {str(p): file_hash(p) for source in sources for p in source.rglob("*") if p.is_file()}
    def forbidden(*args, **kwargs):
        raise AssertionError("Read-only rendering cannot refit observed networks, enrichment or association")
    monkeypatch.setattr("katana_next.reports.network_views.matched_module_associations", forbidden)
    monkeypatch.setattr("katana_next.workflows.module_functions.enrich_modules", forbidden)
    result = build_presentation_assets(root, plan, receipts, transport)
    native = []
    for path in sorted((root / "presentation").glob("network-*-inputs.json")):
        proof = json.loads(path.read_text())
        cid, identity = proof["network_identity"]["contrast_ids"][0], proof["network_identity"]
        overview = next(a for a in result["extra_assets"] if a["contrast_ids"] == [cid]
                        and a["metadata"].get("network_view_role") == "all_module_sizes")
        image = root / "presentation" / (cid + "-dendrogram.png")
        shutil.copyfile(overview["path"], image)
        native.append({"artifact_id": cid + "-dendrogram", "analysis_id": "programs.cosplice",
            "name": "Recorded tree fixture", "kind": "file", "path": str(image), "status": "complete",
            "contrast_ids": [cid], "metadata": {"completion_verified": True, "network_identity": identity,
                "flow_stage": "network.fit", "flow_artifact_key": "dendrogram", "network_view_role": "dendrogram",
                "completion_scope": "observed_network_baseline_only", "bootstrap_enabled": False,
                "data_links": [{"role": "plot_data", "label": "Recorded linkage fixture",
                    "path": str(root / proof["linkage"]), "contrast_id": cid,
                    "unit": "linkage fixture", "policy": None}]}})
    write_json(root / "artifacts.json", {"artifacts": _relative(native, root)})
    write_json(root / "report-project.json", {"extra_assets": _relative(result["extra_assets"], root)})
    _seal_fixture_collection(root, plan, receipts)
    report = tmp_path / "report"; report.mkdir()
    shutil.copytree(ASSET_ROOT, report / "assets")
    builder = Builder(load_source({"project_manifest": str(root / "project.json"), "artifacts": native,
                      "extra_assets": result["extra_assets"], "report_require_plot_data": True}),
                      report, "stable_significant")
    builder.collect_rich({"artifacts": [*native, *result["extra_assets"]]})
    manifest = builder.write()
    assert before == {str(p): file_hash(p) for source in sources for p in source.rglob("*") if p.is_file()}
    return root, plan, result, manifest, report


def test_disabled_report_complete_with_actual_kme_and_no_fake_confidence(disabled_report):
    root, plan, result, manifest, report = disabled_report
    assert network_report_gaps(root, plan, manifest, report) == []
    trees = [row for row in manifest['panels']
             if row.get('metadata', {}).get('network_view_role') == 'dendrogram']
    assert len(trees) == 2 and all(row.get('output_only') is True for row in trees)
    main = (report/'REPORT.html').read_text(encoding='utf-8')
    outputs = (report/'outputs/index.html').read_text(encoding='utf-8')
    for tree in trees:
        assert f'id="{tree["id"]}"' not in main
        assert '../'+tree['image'] in outputs
    assert all(not row["unavailable"] for row in result["extra_views"])
    assert not any(row.get("metadata", {}).get("component") == "network_bootstrap" for row in manifest["unavailable"])
    for path in (root / "presentation").glob("network-*-inputs.json"):
        proof = json.loads(path.read_text())
        assert len(proof["task_receipts"]) == 3
        assert proof["bootstrap"] == {"enabled": False, "status": "disabled_by_request", "hub_metric": "abs_kME"}
        assert "bootstrap_members" not in proof and "bootstrap_modules" not in proof
        validate_network_view_proof(path, [], proof["contrasts"])
        directory = path.with_name(path.name.removesuffix("-inputs.json"))
        receipt = json.loads((directory / "network-report.receipt.json").read_text())
        assert receipt["bootstrap_accepted"] is False and receipt["bootstrap_status"] == "disabled_by_request"
        members = pd.read_csv(directory / "tables/all_module_members.tsv", sep="\t")
        assigned = members.module.ne("grey")
        assert members.loc[assigned, "abs_kME"].equals(members.loc[assigned, "kME"].abs())
        assert not members.loc[assigned, "kME"].equals(members.loc[assigned, "module_membership"])
        assert members.loc[assigned, "hub_rank_abs_kME"].tolist() == [2., 1., 1., 2.]
        assert members.set_index("event_key").loc["e2", "kME"] == -.95
        browser = (directory / "browser/ALL_MODULES.html").read_text(encoding="utf-8")
        assert "kME is not resampling confidence" in browser
        import re
        graph = json.loads(re.search(r'<script id="module-display-data"[^>]*>(.*?)</script>',
                                     browser, re.S).group(1))
        assert set(graph["modules"]) == {"M01", "M02"}
        assert graph["methods"]["matched_native_edges"] == 2
        assert graph["methods"]["maximum_native_adjacency_difference"] == 0
        assert graph["methods"]["network_discovery_runs"] == graph["methods"]["TOM_runs"] == 0
        assert graph["methods"]["association_tests"] == graph["methods"]["enrichment_tests"] == 0
        assert all(len(module["node_indices"]) == 2 and len(module["edges"]) == 1
                   for module in graph["modules"].values())
        changed = copy.deepcopy(proof); changed["bootstrap"]["status"] = "completed"
        write_json(path, changed)
        with pytest.raises(ValueError, match="policy differs"):
            validate_network_view_proof(path, [], proof["contrasts"])
        write_json(path, proof)


@pytest.mark.parametrize("role", ["dendrogram", "module_browser", "all_module_continuous_metrics",
                                  "module_enrichment_overview", "matched_module_effects"])
def test_disabled_still_requires_observed_tree_browser_members_and_complete_statistics(disabled_report, role):
    root, plan, _, manifest, report = disabled_report
    changed = copy.deepcopy(manifest)
    for key in ("panels", "complete_results"):
        changed[key] = [row for row in changed[key] if row.get("metadata", {}).get("network_view_role") != role]
    changed["unavailable"] = []
    gaps = network_report_gaps(root, plan, changed, report)
    assert gaps and all("absent" in gap["reason"] for gap in gaps)
