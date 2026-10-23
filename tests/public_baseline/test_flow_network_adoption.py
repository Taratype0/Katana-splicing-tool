import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

from katana_next.flow.network_adoption import (prepare_plan, load_plan, run_replicate,
    reduce_replicates, draw_matrix, replay_discoveries, verify_replicate, digest, partition_settings)


def fixture(tmp_path):
    rng = np.random.default_rng(31)
    matrix = pd.DataFrame(rng.uniform(.1, .9, (40, 12)),
                          index=[f"e{i}" for i in range(40)], columns=[f"s{i}" for i in range(12)])
    matrix.to_csv(tmp_path / "matrix.tsv", sep="\t")
    pd.DataFrame({"event_key":matrix.index,"module":["M01"]*20+["M02"]*20}).to_csv(tmp_path / "assignments.tsv",sep="\t",index=False)
    pd.DataFrame({"sample_id":matrix.columns,"biological_sample_id":matrix.columns,
                  "condition":["a"]*6+["b"]*6}).to_csv(tmp_path / "samples.tsv",sep="\t",index=False)
    (tmp_path / "config.yaml").write_text(yaml.safe_dump({"project_name":"fixture",
        "network":{"minimum_network_samples":3},"consensus":{"bootstrap_iterations":3},
        "uncertainty":{"network_mode":"diagnostic_only"}}))
    path = prepare_plan(config=tmp_path / "config.yaml",matrix=tmp_path / "matrix.tsv",
        assignments=tmp_path / "assignments.tsv",samples=tmp_path / "samples.tsv",destination=tmp_path / "plan")
    return path


def labels(matrix, settings, **kwargs):
    # A deterministic test double isolates scheduling/replay equivalence from R installation.
    values = np.asarray(["M01"]*20+["M02"]*20,dtype=object)
    if matrix.iloc[0,0] > .5: values[[0,20]]=values[[20,0]]
    return values


def test_parallel_replay_matches_complete_original_consensus(tmp_path,monkeypatch):
    import katana_cosplice.consensus as consensus
    import katana_cosplice.network as network
    path=fixture(tmp_path)
    plan,cfg,matrix,assigned,samples=load_plan(path)
    monkeypatch.setattr(network,"discover_module_labels",labels)
    monkeypatch.setattr(consensus,"discover_module_labels",labels)
    direct=consensus.build_bootstrap_coassignment_consensus(matrix,assigned,samples,cfg.network,cfg.consensus,
        reliability_floor=cfg.uncertainty.reliability_floor,reliability_power=cfg.uncertainty.reliability_power,
        gene_prior_weight=cfg.gene_prior.weight)
    directories=[]
    for iteration in range(1,4):
        target=tmp_path / f"replicate-{iteration}"
        result=run_replicate(path,iteration,target)
        assert result["status"]=="completed" and result["full_data_network_refit"] is False
        directories.append(target)
    reduced,proof=reduce_replicates(path,directories,require_native=False)
    assert proof["full_data_network_refit"] is False
    for name in ("assignments","module_summary","program_scores","bootstrap_runs",
                 "coassignment_edges","diagnostics","parameter_search","selection_audit"):
        pd.testing.assert_frame_equal(getattr(direct,name),getattr(reduced,name),check_exact=True)
    np.testing.assert_array_equal(direct.coassignment_matrix,reduced.coassignment_matrix)
    assert consensus.discover_module_labels is labels


def test_reuse_keeps_existing_bytes_and_never_calls_discovery(tmp_path,monkeypatch):
    import katana_cosplice.network as network
    path=fixture(tmp_path); target=tmp_path / "replicate"
    monkeypatch.setattr(network,"discover_module_labels",labels)
    run_replicate(path,1,target)
    before={p.name:(digest(p),p.stat().st_mtime_ns) for p in target.iterdir()}
    monkeypatch.setattr(network,"discover_module_labels",lambda *a,**k:pytest.fail("completed resample recomputed"))
    assert run_replicate(path,1,target)["reused_without_recomputation"]
    assert before=={p.name:(digest(p),p.stat().st_mtime_ns) for p in target.iterdir()}


def test_changed_source_or_partial_output_is_rejected(tmp_path,monkeypatch):
    import katana_cosplice.network as network
    path=fixture(tmp_path); target=tmp_path / "replicate"
    monkeypatch.setattr(network,"discover_module_labels",labels)
    run_replicate(path,1,target)
    with (target / "labels.tsv").open("a") as out: out.write("bad\tM03\n")
    with pytest.raises(ValueError,match="output changed"):run_replicate(path,1,target)
    with (tmp_path / "matrix.tsv").open("a") as out:out.write("bad\t0\n")
    with pytest.raises(ValueError,match="input changed"):load_plan(path)


def test_native_completion_is_required_by_default(tmp_path,monkeypatch):
    import katana_cosplice.network as network
    path=fixture(tmp_path); target=tmp_path / "replicate"
    monkeypatch.setattr(network,"discover_module_labels",labels)
    monkeypatch.delenv("SLURM_JOB_ID",raising=False)
    monkeypatch.delenv("SLURM_ARRAY_TASK_ID",raising=False)
    run_replicate(path,1,target)
    with pytest.raises(ValueError,match="native execution identity"):
        verify_replicate(path,1,target,require_native=True)


def test_replay_rejects_different_input_and_restores_callback(tmp_path):
    import katana_cosplice.consensus as consensus
    path=fixture(tmp_path);plan,cfg,matrix,assigned,samples=load_plan(path)
    original=consensus.discover_module_labels
    with pytest.raises(AssertionError):
        with replay_discoveries(plan,matrix,cfg,[np.asarray(["M01"]*40)]*3):
            consensus.discover_module_labels(matrix,partition_settings(cfg))
    assert consensus.discover_module_labels is original


def test_no_baseline_substitution_or_partial_reduction(tmp_path):
    path=fixture(tmp_path)
    with pytest.raises(ValueError,match="Incomplete"):
        reduce_replicates(path,[],require_native=False)
    assignments=pd.read_csv(tmp_path / "assignments.tsv",sep="\t").iloc[:-1]
    assignments.to_csv(tmp_path / "assignments.tsv",sep="\t",index=False)
    with pytest.raises(ValueError,match="exactly cover"):
        prepare_plan(config=tmp_path / "config.yaml",matrix=tmp_path / "matrix.tsv",
            assignments=tmp_path / "assignments.tsv",samples=tmp_path / "samples.tsv",destination=tmp_path / "plan2")


def test_source_acceptance_must_prove_all_three_immutable_inputs(tmp_path):
    fixture(tmp_path)
    proof=tmp_path / "accepted.json"
    proof.write_text(json.dumps({"passed":True,"files":[{"sha256":digest(tmp_path / "matrix.tsv")}] }))
    with pytest.raises(ValueError,match="does not prove"):
        prepare_plan(config=tmp_path / "config.yaml",matrix=tmp_path / "matrix.tsv",
            assignments=tmp_path / "assignments.tsv",samples=tmp_path / "samples.tsv",
            destination=tmp_path / "plan2",source_acceptance=proof)


def test_native_array_success_is_checked_using_actual_job_identity(tmp_path,monkeypatch):
    from types import SimpleNamespace
    import katana_cosplice.network as network
    import katana_next.flow.network_adoption as adoption
    path=fixture(tmp_path);target=tmp_path / "replicate"
    monkeypatch.setattr(network,"discover_module_labels",labels)
    monkeypatch.setenv("SLURM_ARRAY_JOB_ID","900")
    monkeypatch.setenv("SLURM_ARRAY_TASK_ID","1")
    run_replicate(path,1,target)
    monkeypatch.setattr(adoption.subprocess,"run",lambda *args,**kwargs:SimpleNamespace(stdout="900_1|RUNNING|0:0\n"))
    with pytest.raises(ValueError,match="not actually completed"):
        verify_replicate(path,1,target,require_native=True)
    monkeypatch.setattr(adoption.subprocess,"run",lambda *args,**kwargs:SimpleNamespace(stdout="900_1|COMPLETED|0:0\n900_1.batch|COMPLETED|0:0\n"))
    proof,_=verify_replicate(path,1,target,require_native=True)
    assert proof["native_job_id"]=="900_1" and "COMPLETED" in proof["verified_sacct"]
