from types import SimpleNamespace

import numpy as np
import pandas as pd

from scripts import rerun_cosplice_production_benchmark as benchmark
from katana_cosplice.resampling import bootstrap_plan


def test_nested_bootstrap_keeps_duplicate_paired_blocks_separate():
    samples = pd.DataFrame({"sample_id": ["a0", "a1", "b0", "b1"],
                            "donor_id": ["a", "a", "b", "b"],
                            "condition": ["WT", "MUT"]*2})
    result = benchmark.nested_sample_frame(samples, samples.sample_id, [0, 1, 0, 1])
    assert result.donor_id.tolist() == ["a__copy_0", "a__copy_0", "a__copy_1", "a__copy_1"]
    plan = bootstrap_plan(result, result.sample_id)
    assert len(plan.blocks) == 2
    assert [list(block) for block in plan.blocks] == [[0, 1], [2, 3]]


def test_enrichment_background_includes_unassigned_events_in_every_arm():
    events = ["G1|ES|1", "G2|ES|2", "G3|ES|3", "G4|ES|4"]
    arms = {"full": np.array(["M1", "M1", "M2", "M2"]),
            "filtered": np.array(["M1", "grey", "M2", "grey"])}
    result = benchmark.enrichment(arms, events, {"term": {"g1", "g2"}}, "test")
    assert result.background_genes.eq(4).all()
    assert set(result.arm) == set(arms)


def test_partition_reuse_requires_identical_membership_and_grey_status():
    assert benchmark.partition_label_mapping(["M1", "M1", "grey"], ["C1", "C1", "grey"]) == {"M1":"C1", "grey":"grey"}
    assert benchmark.partition_label_mapping(["M1", "M1", "grey"], ["C1", "grey", "grey"]) is None
    assert benchmark.partition_label_mapping(["M1", "M1", "grey"], ["grey", "grey", "C1"]) is None
    assert benchmark.partition_label_mapping(["M1", "M2"], ["C1", "C1"]) is None


def test_matching_uses_per_module_counts_not_total_count(tmp_path, monkeypatch):
    matrix = pd.DataFrame(np.zeros((6, 4)), index=list("abcdef"))
    baseline = np.array(["M1"]*3 + ["M2"]*3)
    refined = np.array(["M1", "grey", "grey", "M2", "M2", "grey"])
    membership = pd.DataFrame({"event_key": matrix.index, "consensus_module": refined,
                                "full_data_kME": [.2, .9, .8, .1, .9, .8]})
    monkeypatch.setattr(benchmark, "discover_module_labels", lambda *a: baseline)
    monkeypatch.setattr(benchmark, "build_bootstrap_coassignment_consensus", lambda *a, **kw:
                        SimpleNamespace(assignments=membership, bootstrap_runs=pd.DataFrame(),
                                        module_summary=pd.DataFrame()))
    arms = benchmark.fit_arms(matrix, None, None, None, tmp_path/"fit")
    assert arms["WGCNA_size_matched"].tolist() == ["grey", "M1", "grey", "grey", "M2", "M2"]
