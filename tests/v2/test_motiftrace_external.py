from collections import Counter

import numpy as np
import pandas as pd
import pytest

from katana_next.science.motiftrace_external import PeakIndex, dinucleotide_shuffle, interval_from_rna, overlapping_hit_positions, tpr_at_fpr


def test_half_open_strand_projection():
    assert interval_from_rna((100,150), "+", 2,8) == (102,108)
    assert interval_from_rna((100,150), "-", 2,8) == (142,148)
    with pytest.raises(ValueError):
        interval_from_rna((100,150), "-", -1,8)


def test_exact_strand_interval_overlap():
    index = PeakIndex(pd.DataFrame([dict(chr="chr1",strand="-",start=100,end=110),dict(chr="chr1",strand="-",start=102,end=103)]))
    assert index.overlap("chr1","-",109,111)
    assert not index.overlap("chr1","-",110,120)
    assert not index.overlap("chr1","+",100,110)


@pytest.mark.parametrize("word", ["ACGT"*62+"AC", "T"*250, "AAAACCCCGGGGTTTT"*15+"ACGTACGTAC"])
def test_null_preserves_all_dinucleotides_and_terminals(word):
    shuffled = dinucleotide_shuffle(word,np.random.default_rng(20261015))
    assert Counter(zip(word,word[1:])) == Counter(zip(shuffled,shuffled[1:]))
    assert word[0] == shuffled[0] and word[-1] == shuffled[-1]


def test_equal_FPR_ties_are_not_broken_in_favor_of_a_method():
    assert tpr_at_fpr(np.zeros(100),np.zeros(100)) == .05
    assert tpr_at_fpr(np.ones(100),np.zeros(100)) == 1
    assert tpr_at_fpr(np.zeros(100),np.ones(100)) == 0


def test_overlapping_words_not_silently_deduplicated():
    assert overlapping_hit_positions("TTTTTTT","TTTTTT") == [(0,6),(1,7)]
