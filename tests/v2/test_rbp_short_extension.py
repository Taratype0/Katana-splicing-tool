import numpy as np
from scipy.stats import mannwhitneyu

from katana_next.science.rbp_short_extension import (
    PositionIndex, consensus_scores, consensus_windows, pwm_windows,
    rank_test_matrix, window_validity,
)


def test_short_exon_can_have_full_motif_but_not_full50_window():
    score,hits,masks=consensus_scores(["ACTT","AC","TTTT"],"CTT")
    assert score[0]==.75 and hits[0]==[(1,4)]
    assert np.isnan(score[1]) and score[2]==0
    assert not window_validity(["ACTT"],50,1).any()


def test_tiny_exon_can_support_a_motif_across_both_inclusion_junctions():
    index=PositionIndex(["GGACTT"],3,junctions=[(2,3)])
    _,owners,local=index.positions(3)
    assert local.tolist()==[0,1,2]
    score,hits,_=consensus_scores(["GGACTT"],"GAC",junctions=[(2,3)])
    # Eligible crossing starts 0/1/2 cover RNA bases [0,5), not terminal base 5.
    assert hits==[[(1,4)]] and score[0]==3/5


def test_pwm_small_windows_equal_independent_enumeration():
    sequences=["ACTTACTT","ACNTACTT","AC"]
    matrix=np.array([[.1,.7,.1,.1],[.1,.1,.1,.7],[.1,.1,.1,.7]])
    background=np.full(4,.25)
    index=PositionIndex(sequences,3)
    full,_,site=index.pwm(matrix,background)
    valid=window_validity(sequences,4,5)
    scores=pwm_windows(site,3,3,4,5,valid)
    smooth=(20*matrix+.1*background)/20.1
    for i,sequence in enumerate(sequences):
        for start in range(5):
            window=sequence[start:start+4]
            if len(window)<4 or set(window)-set("ACGT"):
                assert np.isnan(scores[i,start])
            else:
                odds=[np.prod([smooth[k,"ACGT".index(window[a+k])]/.25 for k in range(3)]) for a in range(2)]
                assert np.isclose(scores[i,start],np.log2(1+np.mean(odds)))
    assert np.isnan(full[-1])


def test_short_window_coverage_includes_overlap():
    _,hits,masks=consensus_scores(["AAAAA"],"AAA")
    assert hits==[[(0,3),(1,4),(2,5)]]
    result=consensus_windows(masks,4,2,window_validity(["AAAAA"],4,2),3)
    np.testing.assert_array_equal(result,[[1.,1.]])


def test_rank_matrix_matches_fixed_single_sided_test_and_preserves_na():
    a=np.arange(25.)[:,None]
    b=np.arange(25.)[:,None]/2
    p,*_=rank_test_matrix(a,b)
    assert p[0]==mannwhitneyu(a[:,0],b[:,0],alternative="greater",method="asymptotic").pvalue
    assert np.isnan(rank_test_matrix(a[:3],b)[0][0])
