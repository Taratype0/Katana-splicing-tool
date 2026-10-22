import numpy as np
import pandas as pd

from katana_next.science.frozen_event_transport import (
    ExonInterval, aggregate_junction_runs, chain_blocks, frozen_scores,
    geo_samples, graph_interval, inclusion_proportion, junction_triplet, translate_interval,
)


def test_exact_junction_coordinates_and_negative_strand():
    exons = [ExonInterval('chr1', a, b, '-') for a, b in [(100, 150), (200, 250), (300, 350)]]
    assert junction_triplet(exons[::-1]) == (('chr1',151,200,'-'),('chr1',251,300,'-'),('chr1',151,300,'-'))


def test_contiguous_subexons_not_noncontiguous_composites():
    graph = {('G','2.1'): ExonInterval('chr1', 100, 150, '+'),
             ('G','2.2'): ExonInterval('chr1', 150, 180, '+'),
             ('G','3'): ExonInterval('chr1', 200, 250, '+')}
    assert graph_interval('G', '2.1:2.2', graph) == ExonInterval('chr1',100,180,'+')
    assert graph_interval('G', '2.1:3', graph) is None


def test_unique_whole_interval_liftover_and_ambiguity():
    blocks = chain_blocks(['chain 100 chr1 1000 + 100 200 chr2 2000 + 300 400 1', '100'])
    assert translate_interval(ExonInterval('chr1', 120, 150, '+'), blocks) == ExonInterval('chr2',320,350,'+')
    assert translate_interval(ExonInterval('chr1', 50, 150, '+'), blocks) is None
    blocks['chr1'].append((100, 200, 'chr3', 400, 2000, '+'))
    assert translate_interval(ExonInterval('chr1', 120, 150, '+'), blocks) is None


def test_reverse_liftover_changes_strand_not_length():
    blocks = chain_blocks(['chain 100 chr1 1000 + 100 200 chr2 2000 - 300 400 1', '100'])
    result = translate_interval(ExonInterval('chr1', 120, 150, '+'), blocks)
    assert result == ExonInterval('chr2',1650,1680,'-')


def test_technical_runs_added_before_psi():
    runs = pd.DataFrame({'rail_id':[1,2,3], 'experiment_acc':['X','X','Y']})
    samples = pd.DataFrame({'experiment_acc':['X','Y'], 'sample_id':['a','b']})
    observed = aggregate_junction_runs(np.array([[3,7,4],[2,8,5]]), [1,2,3], runs, samples)
    np.testing.assert_equal(observed, [[10,4],[10,5]])


def test_low_depth_is_missing_not_zero():
    psi, total = inclusion_proportion([0,10,20], [0,10,0], [0,0,10])
    assert np.isnan(psi[0]) and psi[1] == 1 and psi[2] == .5
    np.testing.assert_equal(total, [0,10,20])


def test_projection_does_not_refit_and_reports_loss():
    loads = pd.DataFrame({'event_key':['a','b'], 'center':[.5,.5], 'scale':[.1,.1], 'loading':[.6,.8]})
    data = pd.DataFrame({'x':[.6,.4], 'y':[.5,np.nan]}, index=['a','b'])
    scores, coverage = frozen_scores(data, loads)
    assert np.isclose(scores.x, -.2) and np.isnan(scores.y)
    assert coverage['member_fraction'] == 1
    scores, coverage = frozen_scores(data.loc[['a']], loads)
    assert scores is None and coverage['member_fraction'] == .5


def test_official_normal_title_alias_keeps_patient_and_histology_selection():
    text = '''^SAMPLE = GSM1
!Sample_title = L511T
!Sample_characteristics_ch1 = histology: 2
!Sample_relation = SRA: https://x/SRX1
^SAMPLE = GSM2
!Sample_title = matched sample_L511N
!Sample_relation = SRA: https://x/SRX2
^SAMPLE = GSM3
!Sample_title = L563T
!Sample_characteristics_ch1 = histology: 1
!Sample_relation = SRA: https://x/SRX3
^SAMPLE = GSM4
!Sample_title = matched sample_L563N
!Sample_relation = SRA: https://x/SRX4
'''
    samples = geo_samples(text)
    assert samples.sample_id.tolist() == ['L511T','L511N','L563T','L563N']
    assert samples.included.tolist() == [True,True,False,False]
    assert samples.patient_id.tolist() == ['L511','L511','L563','L563']
    assert samples.loc[1,'provider_title'] == 'matched sample_L511N'
