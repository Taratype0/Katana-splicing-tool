import numpy as np
import pandas as pd

from scripts.public134_external_validation_20261016 import base_occupancy,score_geometry,tokens


def test_real_stranded_coordinates_and_acceptor_grid():
    source=pd.Series({'sequence':'ACGTACGT','blocks_json':'[[100,108]]','chr':'chr1','strand':'-'})
    occupied,invalid=base_occupancy(source,{('chr1','-'):(np.array([101]),np.array([103]))})
    np.testing.assert_array_equal(np.diff(occupied),[0,0,0,0,0,1,1,0])
    np.testing.assert_array_equal(score_geometry(occupied,invalid,2,7,True),[1,1,1,0,0,0,0])


def test_missing_bases_remain_nan_no_padding():
    source=pd.Series({'sequence':'ACNT','blocks_json':'[[0,4]]','chr':'chr1','strand':'+'})
    occupied,invalid=base_occupancy(source,{})
    value=score_geometry(occupied,invalid,2,5,False)
    assert value[0]==0 and np.isnan(value[1:]).all()
    assert np.isnan(score_geometry(occupied,invalid,0,1,False)).all()


def test_target_tokens_do_not_match_unrelated_substrings():
    assert 'QKI' in tokens('QKI;PCBP2')
    assert 'QKI' not in tokens('QKI_like;QK')
    assert 'RBFOX2' not in tokens('RBFOX1;RBFOX3')
