import numpy as np
import pytest
from katana_next.science.splice_structure import moments,fit_competition_model,project_competition_model

def test_factorial_variance_removes_known_binomial_noise():
    rng=np.random.default_rng(42);n=np.full((3,200000),12);p=np.array([.1,.4,.8])[:,None]
    y=rng.binomial(n,p);m=moments(y,n)
    assert np.max(abs(m['latent_variance']))<.0002
    assert np.allclose(m['mean_sampling_variance'],(p*(1-p)/12).ravel(),atol=.00008)

def test_missing_input_returns_frozen_reference_and_no_information():
    rng=np.random.default_rng(4);n=np.full((20,50),80);z=rng.normal(size=(1,50));p=1/(1+np.exp(-rng.normal(size=(20,1))*z))
    model,_=fit_competition_model(rng.binomial(n,p),n,[str(i) for i in range(20)],[str(i) for i in range(20)],rank=2)
    pred,scores,audit=project_competition_model(model,np.zeros((20,3)),np.zeros((20,3)))
    assert np.allclose(pred[:,0],model['means']) and not audit.uses_splice_information.any()
    assert not scores.any()

def test_requires_real_counts_and_one_independent_competition():
    with pytest.raises(ValueError,match='integer'):fit_competition_model(np.full((3,20),.2),np.ones((3,20)),list('abc'),list('abc'))
    with pytest.raises(ValueError,match='one focal'):fit_competition_model(np.ones((3,20)),np.full((3,20),5),list('abc'),['same']*3)

def test_frozen_model_projection_and_gene_expression_scale_invariance():
    # Invariance of latent definition, not of read precision: proportions are
    # unchanged by depth; actual sampling uncertainty must still change.
    rng=np.random.default_rng(7);y=rng.integers(1,40,(30,50));n=np.full_like(y,50)
    model,_=fit_competition_model(y,n,[str(i) for i in range(30)],[str(i) for i in range(30)],rank=3)
    import json
    before=json.dumps(model,sort_keys=True)
    a,z,qa=project_competition_model(model,y[:,:5],n[:,:5])
    b,_,qb=project_competition_model(model,y[:,:5],n[:,:5])
    assert np.array_equal(a,b) and json.dumps(model,sort_keys=True)==before
    assert (a>0).all() and (a<1).all()
