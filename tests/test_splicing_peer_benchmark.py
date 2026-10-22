import json

import numpy as np
import pandas as pd
import pytest

from scripts.benchmark_splicing_peers import (
    canonical_labels, freeze, group_means, microexon_enrichment,
    microexon_length, prepare_busskamp, summary,
)
from katana_cosplice.partition_metrics import partition_recovery
from scripts.summarize_splicing_peers import gene_enrichment, unique_exon_microexons


def test_original_cluster_zero_is_not_grey():
    labels = canonical_labels([0, 0, 1, 1, -1])
    assert labels.tolist() == ["M0", "M0", "M1", "M1", "grey"]
    assert partition_recovery(labels, labels)["reference_assigned_events"] == 4


def test_condition_means_preserve_declared_order():
    matrix = pd.DataFrame([[.1, .4, .3, .6]])
    np.testing.assert_allclose(group_means(matrix, ["B", "A", "B", "A"]), [[.2, .5]])
    with pytest.raises(ValueError):
        group_means(matrix, ["A"])


def test_frozen_config_rejects_changed_inputs(tmp_path):
    freeze(tmp_path, {"seed": 1, "input_hash": "old"})
    freeze(tmp_path, {"seed": 1, "input_hash": "old"})
    with pytest.raises(ValueError):
        freeze(tmp_path, {"seed": 1, "input_hash": "new"})


def test_headerless_author_files_and_missingness(tmp_path):
    ids = ["first;SE:chr1:10-20:30-40:+", "bad;SE:chr1:10-20:30-40:+", "weak;SE:chr1:10-20:30-40:+"]
    pd.DataFrame([[.3, .01]*3, [.3, .01]*3, [.1, .01]*3], index=ids).to_csv(
        tmp_path/"Busskamp_multicondition.SE.dpsi.gz", sep="\t", header=False)
    values = np.full((3, 12), .3)
    values[1, 0] = -1
    pd.DataFrame(values, index=ids).to_csv(tmp_path/"Busskamp_multicondition.SE.psivec.gz", sep="\t", header=False)
    matrix, meta, _, audit = prepare_busskamp(tmp_path)
    assert matrix.index.tolist() == [ids[0]]
    assert audit["tutorial_significant"] == 2
    assert audit["excluded_missing_or_invalid"] == 1
    assert meta.condition.nunique() == 4


def test_microexon_inclusive_coordinates_and_shared_background():
    assert microexon_length("ENSG1;SE:chr1:5-10:36-90:+") == 27
    ids = [f"G{i};SE:chr1:5-10:{36 if i<3 else 100}-200:+" for i in range(6)]
    matrix = pd.DataFrame([[.1, .1, .6, .6]]*6, index=ids)
    result = microexon_enrichment(matrix, {"A":np.array(["M1"]*3+["grey"]*3)}, ["before"]*2+["after"]*2)
    assert result.iloc[0].universe_events == 6
    assert result.iloc[0].microexons == 3
    assert result.iloc[0].mean_last_minus_first_PSI == pytest.approx(.5)


def test_all_grey_not_reported_as_useful_retained_stability():
    labels=np.repeat("grey", 4)
    stability=pd.DataFrame([{"arm":"A", **partition_recovery(labels, labels)}])
    row=summary({"A":labels},stability,{"A":1}).iloc[0]
    assert row.modules == 0
    assert row.assigned_fraction == 0
    assert np.isnan(row.mean_ari_jointly_assigned)


def test_gene_enrichment_counts_unique_genes_and_includes_grey():
    assignments=pd.DataFrame({"A":["M1","M1","grey","M2"]},index=["e1","e2","e3","e4"])
    mapping=pd.Series(["G1","G1","G2","G3"],index=assignments.index)
    result=gene_enrichment(assignments,mapping,{"pathway":{"g1","g2"}},"test")
    row=result[result.module.eq("M1")].iloc[0]
    assert row.background_genes==3
    assert row.module_genes==1
    assert row.term_background_genes==2


def test_unique_exon_sensitivity_collapses_flanking_variants():
    ids=["G1;SE:chr1:5-10:36-200:+","G1;SE:chr1:6-10:36-201:+","G2;SE:chr1:5-50:120-200:+"]
    matrix=pd.DataFrame([[.1,.5]]*3,index=ids)
    assignments=pd.DataFrame({"A":["M1","M1","grey"]},index=ids)
    result,_=unique_exon_microexons(matrix,assignments,["early","late"])
    assert result.iloc[0].unique_skipped_exons==1
    assert result.iloc[0].background_unique_exons==2
