import numpy as np
import pandas as pd
import pytest

from katana_next.science.motiftrace_usage import candidate_positions, canonical_dexseq
from scripts.report_public_usage_mapp_20261016 import overlap


def test_candidate_tuple_is_not_mistaken_for_number_of_matches():
    starts, scores, valid = candidate_positions('AAAAAA', {'kind':'consensus','width':3,'pattern':'AAA'})
    assert starts.tolist() == [0,1,2,3]
    assert len(starts) == 4
    assert valid == 4
    np.testing.assert_equal(scores, np.ones(4))


def test_shared_overlapping_counts_are_not_claimed_native_nonoverlap_counts():
    starts, _, _ = candidate_positions('AAAAAA', {'kind':'consensus','width':3,'pattern':'AAA'})
    assert len(starts) == 4
    assert 'AAAAAA'.count('AAA') == 2


def test_clip_overlap_respects_strand_and_half_open_coordinates():
    d = pd.DataFrame({'chr':['chr1','chr1','chr1'],'start0':[10,20,10],
                      'end0':[20,30,20],'strand':['+','+','-']})
    d.attrs['accession'] = 'unit_clip_1'
    assert overlap(d,'chr1','+',[(20,21)]) == 1
    assert overlap(d,'chr1','+',[(30,31)]) == 0
    assert overlap(d,'chr1','-',[(10,11)]) == 1
    assert overlap(d,'chr1','+',[(11,12),(13,14)]) == 1


def test_missing_dexseq_bins_do_not_become_accepted_dexseq():
    with pytest.raises(ValueError,match='DEXSeq needs native bin'):
        canonical_dexseq(pd.DataFrame({'GeneID':['g1'],'IJC_SAMPLE_1':['3,4'],'FDR':[.001]}),'KD_minus_Ctrl')


def test_short_exon_remains_scannable_and_zero_valid_starts_is_not_hit():
    starts, _, valid = candidate_positions('CTCG',{'kind':'consensus','width':4,'pattern':'CTCG'})
    assert starts.tolist()==[0] and valid==1
    starts, _, valid = candidate_positions('CTC',{'kind':'consensus','width':4,'pattern':'CTCG'})
    assert len(starts)==0 and valid==0
