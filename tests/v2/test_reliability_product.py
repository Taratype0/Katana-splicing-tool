import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from fastapi.testclient import TestClient

from katana_cosplice.reliability_audit import calibration_grid, fixed_module_audit, thinning_diagnostic, event_junctions
from katana_next.api.app import PlatformRuntime, create_app
from katana_next.client import KatanaProject
from katana_next.full_fixture import create_full_fixture
from katana_next.services.readiness import experimental_design
from katana_next.workflows.base import ArtifactWriter
from katana_next.workflows.registry import WorkflowRegistry
from katana_cosplice.uncertainty import _condition_rho, _condition_rho_rows
from katana_cosplice.config import UncertaintyConfig
from katana_next.workflows.programs import configure_production_uncertainty


def audit_fixture():
    rng = np.random.default_rng(4)
    samples = pd.DataFrame({"sample_id": [f"s{i}" for i in range(24)], "condition": ["A"]*12+["B"]*12})
    effect = np.array([.2]*12+[.8]*12)
    psi = pd.DataFrame(np.clip(effect+rng.normal(0,.02,(4,24)),0,1), index=["e1","e2","e3","e4"], columns=samples.sample_id)
    events = pd.DataFrame([{"event_key": k, "gene_symbol": f"G{i}", "event_type": "SE", "chr": "chr1", "strand": "+",
                           "upstreamEE": 100, "exonStart_0base": 200, "exonEnd": 250, "downstreamES": 400} for i,k in enumerate(psi.index)])
    events.loc[2:, "downstreamES"] = [500,600]
    assignment = pd.DataFrame({"event_key": psi.index,"module":["M1"]*3+["grey"]})
    return psi,samples,events,assignment


def test_fixed_modules_keep_all_events_and_explain_condition():
    psi,samples,events,assignment = audit_fixture()
    result = fixed_module_audit(psi,samples,events,assignment)
    assert len(result["events"]) == 4
    assert result["coverage"].retention_fraction.iloc[0] == .75
    module = result["modules"].iloc[0]
    assert module.mean_observed_r > .98
    assert abs(module.mean_condition_adjusted_r) < .5
    assert module.known_structural_dependency_pairs == 3
    assert "duplicate_structure" in set(result["pairs"].relation)
    assert assignment.module.tolist() == ["M1"]*3+["grey"]


def test_structural_identity_ignores_gene_label_and_uses_strand():
    _,_,events,_ = audit_fixture()
    assert event_junctions(events.iloc[0]) == event_junctions(events.iloc[1])
    other = events.iloc[1].to_dict(); other["strand"] = "-"
    assert not event_junctions(events.iloc[0]) & event_junctions(other)


def test_missing_structure_is_unassessable_not_an_empty_event_list():
    psi,samples,events,assignments = audit_fixture()
    audit = fixed_module_audit(psi,samples,events[["event_key"]],assignments)
    assert len(audit["events"]) == 4
    assert set(audit["pairs"].relation) == {"not_assessable"}


def test_diagnostic_sampling_and_confounded_batch():
    args = audit_fixture()
    a = fixed_module_audit(*args,max_pairs=1,seed=7)
    b = fixed_module_audit(*args,max_pairs=1,seed=7)
    pd.testing.assert_frame_equal(a["pairs"],b["pairs"])
    assert a["modules"].all_possible_pairs.iloc[0] == 3
    assert len(a["pairs"]) == 1
    args[1]["batch"] = args[1].condition
    with pytest.raises(ValueError,match="confounded"):
        fixed_module_audit(*args,include_batch=True)


def test_calibration_is_deterministic_and_not_tuned_to_win():
    a = calibration_grid(events=15,seed=8)
    b = calibration_grid(events=15,seed=8)
    pd.testing.assert_frame_equal(a,b)
    assert len(a) == 24
    assert a.observed_coverage.between(0,1).all()
    assert a.mean_interval_width.gt(0).all()
    assert set(a.method) == {"Jeffreys binomial","CoSplice effective-count approximation"}


def test_chunked_uncertainty_matches_original_scalar_estimator():
    rng = np.random.default_rng(17)
    conditions = np.repeat(["A","B","C"],3)
    n = rng.integers(0,1000,(200,9)).astype(float)
    i = rng.binomial(n.astype(int),rng.uniform(0,1,n.shape)).astype(float)
    i[::7,0] = np.nan
    n[::11,2] = np.nan
    expected = np.array([_condition_rho(a,b,conditions) for a,b in zip(i,n)])
    actual = _condition_rho_rows(i,n,conditions,chunk_size=17)
    np.testing.assert_allclose(actual,expected,atol=1e-12,rtol=1e-12)


def test_production_keeps_measured_psi_and_does_not_silently_filter_by_uncalibrated_rho():
    from types import SimpleNamespace
    config = SimpleNamespace(uncertainty=UncertaintyConfig())
    assert config.uncertainty.network_mode == "posterior"
    configure_production_uncertainty(config,{})
    assert config.uncertainty.network_mode == "diagnostic_only"
    assert not config.uncertainty.filter_low_reliability
    assert not config.uncertainty.posterior_edge_diagnostics
    with pytest.raises(ValueError,match="observed PSI"):
        configure_production_uncertainty(config,{"uncertainty_role":"posterior"})


def test_thinning_preserves_missingness_and_requires_integer_counts():
    samples = pd.DataFrame({"sample_id": ["a","b","c","d"],"condition":["A","A","B","B"]})
    i = pd.DataFrame([[3,10,np.nan,2],[10,12,8,3]],columns=samples.sample_id)
    s = i+5
    result = thinning_diagnostic(i,s,samples,iterations=2)
    assert (result.total_event_sample_cells == 7).all()
    i = i.astype(float)
    i.iloc[0,0] = .5
    with pytest.raises(ValueError,match="integer"):
        thinning_diagnostic(i,s,samples)
    i.iloc[0,0] = -1
    with pytest.raises(ValueError,match="nonnegative"):
        thinning_diagnostic(i,s,samples)


def test_readonly_design_and_preflight_share_cli_web_engine(tmp_path):
    manifest = create_full_fixture(tmp_path / "project")
    with KatanaProject(manifest) as project:
        review = project.design()
        assert any(p["matched_denominator"] for p in review["comparison_pairs"])
        assert len(project.jobs.list()) == 0
        assert not project.plan("programs.diagnostics")["ready"]
        assert project.plan("quality.read_calibration", {"simulated_events": 15})["ready"]
        first = project.run("quality.read_calibration", {"simulated_events":15})
        repeat = project.run("quality.read_calibration", {"simulated_events":15})
        assert first.status == "complete" and repeat.status == "cache_hit"
        assert {j.analysis_id for j in project.jobs.list()} == {"quality.read_calibration"}
    runtime = PlatformRuntime()
    with TestClient(create_app(runtime)) as client:
        opened = client.post("/api/projects/open",json={"manifest_path":str(manifest)}).json()
        prefix = f"/api/projects/{opened['snapshot']['project_id']}"
        assert client.get(prefix+"/design").status_code == 200
        plan = client.post(prefix+"/workflows/quality.read_calibration/plan",json={"parameters": {"simulated_events":15}})
        assert plan.json()["ready"]
        recorded = runtime.jobs.list()
        assert {job.job_id for job in recorded} == {first.job_id, repeat.job_id}
        assert {job.status for job in recorded} == {"complete", "cache_hit"}
    runtime.jobs.shutdown()


def test_fixed_run_workflow_module_browser_and_sashimi_bridge(tmp_path):
    manifest = create_full_fixture(tmp_path / "project")
    psi, samples, events, assignments = audit_fixture()
    # Complete synthetic lineage tests routing, not scientific module discovery.
    for column, value in {"upstreamES": 50, "downstreamEE": 450}.items():
        events[column] = value
    assignments["consensus_membership_confidence"] = [.95, .8, .6, np.nan]
    with KatanaProject(manifest) as project:
        writer = ArtifactWriter(project.context, "programs.cosplice", "synthetic-frozen-discovery", {},
                                WorkflowRegistry().get("programs.cosplice").version)
        matrix = writer.table("2. Build PSI matrix: psi", "matrix__psi.parquet", psi.rename_axis("event_key").reset_index(), entity_type="matrix")
        writer.table("2. Build PSI matrix: events", "matrix__events.parquet", events, entity_type="splice_event")
        writer.table("2. Build PSI matrix: samples", "matrix__samples.parquet", samples, entity_type="sample")
        writer.table("Observed assignments", "network_observed__assignments.parquet", assignments, entity_type="cosplice_module")
        writer.table("Consensus assignments", "consensus_modules__assignments.parquet", assignments.rename(columns={"module": "analysis_module"}), entity_type="cosplice_module")
        writer.finish()
        das = ArtifactWriter(project.context, "differential.as", "synthetic-das", {}, WorkflowRegistry().get("differential.as").version)
        bulk = events.iloc[:1].drop(columns="event_key").assign(gene_symbol="ANOTHER_LABEL", contrast_id="contrast", event_id="das-event", row_id="das-row")
        das.table("Differential AS events", "events.parquet", bulk, entity_type="splice_event")
        das.finish()
        parameters = {"matrix_artifact_id": matrix.artifact_id}
        assert project.plan("programs.diagnostics", parameters)["ready"]
        result = project.run("programs.diagnostics", parameters)
        assert result.status == "complete"
        artifacts = project.artifacts("programs.diagnostics")
        module = next(a for a in artifacts if a.name == "Fixed-module modules")
        evidence = project.module(module.artifact_id, "M1", limit=1)
        assert evidence["discovery_run"] == matrix.run_fingerprint
        assert all(len(t["rows"]) == 1 for t in evidence["tables"])
        assert any(t["row_count"] == 3 for t in evidence["tables"])
        bridge = next(a for a in artifacts if a.name == "Module events ready for Native Sashimi")
        mapped = pd.read_parquet(bridge.path)
        assert set(mapped.event_key) == {"e1", "e2"}
        assert set(mapped.module) == {"M1"}
        assert mapped.row_id.is_unique
        assert project.run("programs.diagnostics", parameters).status == "cache_hit"
        assert {j.analysis_id for j in project.jobs.list()} == {"programs.diagnostics"}
