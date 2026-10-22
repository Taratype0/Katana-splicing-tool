from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd
import pytest
from katana_next.science.bulk_interaction import interaction_observations, fit_read_interaction
from katana_next.science.motif_significance import match_background, matched_enrichment, extract_windows


def _contrast(numerator, denominator):
    return SimpleNamespace(numerator=SimpleNamespace(sample_ids=numerator), denominator=SimpleNamespace(sample_ids=denominator))


def _counts(inc1="80,85,75", inc2="20,25,15", skip1="20,15,25", skip2="80,75,85"):
    return pd.Series({"IJC_SAMPLE_1": inc1, "IJC_SAMPLE_2": inc2, "SJC_SAMPLE_1": skip1, "SJC_SAMPLE_2": skip2,
                      "IncFormLen": 200, "SkipFormLen": 100, "reverse_direction": False})


def test_shared_controls_are_counted_once_and_direction_correct():
    a = _contrast(("a1", "a2", "a3"), ("w1", "w2", "w3"))
    b = _contrast(("b1", "b2", "b3"), ("w1", "w2", "w3"))
    obs, vector, offset = interaction_observations(_counts(), _counts("40,45,35", skip1="60,55,65"), a, b, [])
    assert len(obs) == 9
    assert sum(vector == 0) == 1
    fit = fit_read_interaction(obs, vector, offset)
    assert fit["status"] == "fit"
    assert fit["effect_log_odds_A_minus_B"] > 0
    assert fit["pvalue"] < .05
    reverse = fit_read_interaction(obs, -vector, offset)
    assert reverse["pvalue"] == pytest.approx(fit["pvalue"], abs=1e-5)


def test_count_conflicts_and_partial_overlap_fail_closed():
    a, b = _contrast(("a1", "a2", "a3"), ("w1", "w2", "w3")), _contrast(("b1", "b2", "b3"), ("w1", "w2", "w3"))
    with pytest.raises(ValueError, match="conflicting"):
        interaction_observations(_counts(), _counts(inc2="21,25,15"), a, b, [])
    b.denominator.sample_ids = ("w1", "w2", "w4")
    with pytest.raises(ValueError, match="partially_overlapping"):
        interaction_observations(_counts(), _counts(), a, b, [])


def test_missing_positions_not_shifted_and_bad_lengths_rejected():
    a, b = _contrast(("a1", "a2", "a3"), ("w1", "w2", "w3")), _contrast(("b1", "b2", "b3"), ("v1", "v2", "v3"))
    row = _counts(inc1="80,NA,75")
    obs, c, offset = interaction_observations(row, _counts(), a, b, [])
    assert "a2" not in set(obs.sample_id)
    assert obs.set_index("sample_id").loc["a3", "IJC"] == 75
    assert fit_read_interaction(obs, c, offset, min_replicates=3)["status"] == "insufficient_read_supported_replicates"
    row["IncFormLen"] = np.nan
    with pytest.raises(ValueError, match="missing_effective"):
        interaction_observations(row, _counts(), a, b, [])


def test_matched_background_unique_genes_no_foreground_contamination():
    data = pd.DataFrame({"event_id": ["a", "b", "c", "d", "e"], "gene_id": ["G1", "G1", "G2", "G3", "G4"],
                         "event_type": ["SE"]*5, "gc": [.4, .4, .42, .8, .41], "bases": [200]*5, "ambiguous_fraction": [0]*5})
    pairs = match_background(data, {"a", "b"})
    assert pairs.iloc[0].background_gene == "G4"
    assert pairs.iloc[1].status == "same_gene_additional_event_sites_only"
    assert set(pairs.background_gene) & {"G1"} == set()


def test_enrichment_corrects_full_rbp_family_and_no_hits():
    pairs = pd.DataFrame({"foreground_event": [f"a{i}" for i in range(10)], "background_event": [f"b{i}" for i in range(10)], "status": ["matched"]*10})
    sites = pd.DataFrame({"RBP": ["R1"]*10, "event_id": pairs.foreground_event})
    catalog = pd.DataFrame({"RBP": ["R1", "R2"]})
    result = matched_enrichment(sites, catalog, pairs).set_index("RBP")
    assert result.loc["R1", "enrichment_pvalue"] == 1/1024
    assert result.loc["R1", "enrichment_fdr"] == 2/1024
    assert result.loc["R2", "enrichment_pvalue"] == 1
    assert not result.loc["R2"].isna().any()


def test_minus_strand_windows_keep_genomic_intervals(tmp_path):
    fasta = tmp_path/"genome.fa"
    fasta.write_bytes(b">chr1\nAAAACCCCGGGGTTTT\n")
    Path(str(fasta)+".fai").write_text("chr1\t16\t6\t16\t17\n")
    events = pd.DataFrame({"event_id": ["e"], "event_type": ["SE"], "gene_id": ["G"], "chr": ["chr1"], "strand": ["-"], "exonStart_0base": [4], "exonEnd": [8]})
    windows, _ = extract_windows(events, fasta, 1)
    assert list(windows.sequence) == ["GT", "CG"]
    assert list(windows.start_0based) == [3, 7]


def test_batched_moran_exactly_preserves_stratified_permutation_kernel():
    from katana_next.science.single_cell_programs import _knn_graph, _morans_i, _shuffle_within_groups, moran_permutation_scores
    rng = np.random.default_rng(42)
    values = rng.uniform(size=30)
    values[[0, 7, 19]] = np.nan
    groups = np.repeat(["donor1", "donor2", "donor3"], 10)
    graph = _knn_graph(rng.normal(size=(30, 3)), 6)
    scalar_rng = np.random.default_rng(19)
    expected = np.array([_morans_i(_shuffle_within_groups(values, groups, scalar_rng), graph) for _ in range(31)])
    score, fast = moran_permutation_scores(values, graph, groups, np.random.default_rng(19), 31, batch_size=9)
    assert score == _morans_i(values, graph)
    np.testing.assert_allclose(fast, expected, atol=1e-12, rtol=1e-12)


@pytest.mark.parametrize("extra_background_a", [0, 96])
def test_installed_fimo_matches_exhaustive_background_null(tmp_path, extra_background_a):
    import os
    from itertools import product
    from scipy.stats import binom
    from katana_cosplice.motif_scanner import Motif
    from katana_next.science.motif_significance import calibrated_scan
    runtime = os.environ.get("KATANA_TEST_FIMO_RUNTIME")
    if not runtime:
        pytest.skip("Set KATANA_TEST_FIMO_RUNTIME for the installed MEME-core integration test")
    sequences = ["".join(bases) for bases in product("ACGT", repeat=3)]
    windows = pd.DataFrame([{"window_id": f"w{i}", "event_id": f"e{i}", "sequence": sequence,
                             "role": "background", "chr": "chr1", "strand": "+", "start_0based": i*3, "end_exclusive": i*3+3}
                            for i, sequence in enumerate(sequences)])
    if extra_background_a:
        windows = pd.concat([windows, pd.DataFrame([{
            "window_id": "extra_background", "event_id": "extra_background", "sequence": "A"*extra_background_a,
            "role": "background", "chr": "chr1", "strand": "+", "start_0based": 1000,
            "end_exclusive": 1000+extra_background_a,
        }])], ignore_index=True)
    motif = Motif("RBP", "AAA", "analytical test", pwm=np.tile([.97, .01, .01, .01], (3, 1)))
    sites, _, diagnostics = calibrated_scan([motif], windows, tmp_path/"fimo", Path(runtime), threshold=1.)
    assert diagnostics["tested_positions_all_motifs"] == 64 + max(0, extra_background_a-2)
    pvalues = sites.set_index("window_id").site_pvalue
    for i, sequence in enumerate(sequences):
        expected = binom.sf(sequence.count("A")-1, 3, (49+extra_background_a)/(196+extra_background_a))
        assert pvalues.get(f"w{i}", 1.) == pytest.approx(expected, abs=1e-6)


def test_interaction_selection_is_pinned_not_rebound_to_latest(monkeypatch):
    import katana_next.workflows.bulk_interaction as module
    parent = SimpleNamespace(analysis_id="differential.as")
    monkeypatch.setattr(module, "selected_source", lambda context, parameters: parent)
    def unexpected(*args):
        raise AssertionError("A saved selection must not resolve a different latest artifact")
    monkeypatch.setattr(module, "_latest_table", unexpected)
    assert module.BulkInteractionWorkflow().parent_artifacts(None, {"selection_id": "pinned"}) == [parent]
