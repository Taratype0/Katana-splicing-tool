import numpy as np
import pytest
from scipy.stats import binom
from katana_next.science.splice_copula import validate,fit_margins,intervals,rectangle_loss,fit,predict

def test_rectangles_reduce_to_independent_probability():
    lo=np.array([0.,.2,.99]);hi=np.array([.1,.7,1.])
    expected=-np.log(np.maximum((hi-lo)*(hi[::-1]-lo[::-1]),1e-9)).sum()
    assert rectangle_loss(0,lo,hi,lo[::-1],hi[::-1])==pytest.approx(expected,abs=1e-10)

def test_depth_changes_probability_intervals_for_same_usage():
    margins=dict(coefficients=[[0]],dispersion=[0])
    lo,hi,mu,seen=intervals(margins,[[1,100]],[[2,200]])
    assert hi[0,0]-lo[0,0] > 5*(hi[0,1]-lo[0,1])
    assert mu[0,0]==mu[0,1]==.5

@pytest.mark.parametrize('y,n',[([[1.5]],[[3]]),([[4]],[[3]]),([[0]],[[-1]]),([[np.nan]],[[3]]),([[1]],[[0]])])
def test_invalid_evidence_is_rejected(y,n):
    with pytest.raises(ValueError):validate(y,n)

def test_frozen_projection_does_not_use_focal_site_reads():
    rng=np.random.default_rng(7);n=np.full((4,60),30);y=rng.binomial(n,.3)
    model,audit=fit(y,n,list('abcd'),list('ABCD'),exact=False)
    p,_=predict(model,y,n);other=y.copy();other[0]=0
    q,_=predict(model,other,n)
    np.testing.assert_allclose(p[0],q[0],atol=1e-12)
    assert len(audit['pair_audit'])==6

def test_missing_all_neighbours_has_finite_training_fallback():
    rng=np.random.default_rng(8);n=np.full((4,60),30);y=rng.binomial(n,.3)
    model,_=fit(y,n,list('abcd'),list('ABCD'),exact=False)
    p,audit=predict(model,np.full((4,2),np.nan),np.full((4,2),np.nan))
    assert np.isfinite(p).all() and audit['observed_neighbours']==[0,0]
