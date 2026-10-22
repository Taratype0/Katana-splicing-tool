import numpy as np
import pytest
from scripts.calibrate_public_rmats_counts import split_counts


def test_disjoint_splits_conserve_counts_and_missingness():
    i=np.array([[0,10,np.nan],[3,0,9]],float)
    s=np.array([[0,20,np.nan],[5,10,1]],float)
    a,b,c,d=split_counts(i,s,.25,np.random.default_rng(4))
    np.testing.assert_allclose(a+c,i,equal_nan=True)
    np.testing.assert_allclose(b+d,s,equal_nan=True)
    np.testing.assert_allclose(a+b,np.floor((i+s)*.25),equal_nan=True)
    again=split_counts(i,s,.25,np.random.default_rng(4))
    for x,y in zip((a,b,c,d),again):
        np.testing.assert_allclose(x,y,equal_nan=True)


def test_splits_reject_effective_or_invalid_counts():
    with pytest.raises(ValueError):
        split_counts(np.array([[1.5]]),np.array([[2]]),.5,np.random.default_rng(1))
    with pytest.raises(ValueError):
        split_counts(np.array([[np.inf]]),np.array([[2]]),.5,np.random.default_rng(1))
