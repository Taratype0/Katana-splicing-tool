import numpy as np
from scipy.stats import mannwhitneyu

from katana_next.science.rmaps_sliding_comparison import covered_windows, mw_windows, window_design


def test_full_window_coverage_keeps_overlapping_motifs():
    assert covered_windows("A"*100, "AAA", "positive", 50).tolist() == [1.]*50
    assert covered_windows("A"*100, "AAA", "negative", 50).tolist() == [1.]*50


def test_windows_scan_every_base_and_do_not_fill_missing_as_zero():
    seq = "C"*20+"AAAA"+"C"*76
    score = covered_windows(seq, "AAAA", "positive", 50)
    assert np.isclose(score[0], 4/50)
    assert score[24] == 0
    assert np.isnan(covered_windows("A"*20, "AAA", "positive", 50)).all()
    assert np.isnan(covered_windows("N"+"A"*99, "AAA", "positive", 50)[0])


def test_fixed_window_design_and_splice_mask():
    design = window_design()
    assert len(design) == 1000
    assert (design.window_end_RNA_relative-design.window_start_RNA_relative).eq(50).all()
    upstream = design[design.region.eq("UpstreamIntron")]
    assert upstream.window_end_RNA_relative.iloc[-1] == 0
    assert not upstream.canonical_mask_eligible.iloc[-1]
    downstream = design[design.region.eq("DownstreamIntron")]
    assert not downstream.canonical_mask_eligible.iloc[0]
    assert downstream.canonical_mask_eligible.iloc[6]


def test_vector_mw_matches_independent_scalar_with_ties():
    a = np.tile(np.arange(25)%3, (4, 1)).T.astype(float)
    b = np.tile(np.arange(27)%2, (4, 1)).T.astype(float)
    a[0, 1] = np.nan
    p, _, _, n1, _ = mw_windows(a, b)
    for i in range(4):
        expected = mannwhitneyu(a[np.isfinite(a[:, i]), i], b[:, i], alternative="greater", method="asymptotic").pvalue
        assert np.isclose(p[i], expected)
    assert n1[1] == 24


def test_full_cassette_middle_is_included_in_one_nt_tiles():
    sequence = "C"*110+"AAAA"+"C"*100
    scores = covered_windows(sequence, "AAAA", "positive", len(sequence)-49)
    assert len(scores) == len(sequence)-49
    assert np.isfinite(scores).all()
    assert scores[90] == 4/50
    assert scores[0] == 0
    assert scores[150] == 0
