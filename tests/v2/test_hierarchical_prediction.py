import numpy as np
import pandas as pd
import pytest

from katana_next.science.hierarchical_prediction import hierarchical_draws, paired_moment_rho
from katana_next.science.psi_calibration import length_correct
from katana_next.science.replicate_prediction import heldout_arm


def test_paired_moment_recovers_independent_known_rho():
    rng = np.random.default_rng(108)
    mu = rng.beta(2, 2, 12000)
    rho = .07
    q = rng.beta(mu[:, None]*(1/rho-1), (1-mu[:, None])*(1/rho-1), (len(mu), 2))
    counts = rng.binomial(100, q).astype(float)
    fit = paired_moment_rho(counts, 100-counts)
    assert abs(fit['rho']-rho) < .006
    assert paired_moment_rho(counts[:3], 100-counts[:3]) is None


def test_quadrature_resolution_and_length_transform():
    inclusion = np.array([[0, 0], [1, 3], [50, 90], [5000, 5100], [10000, 10000]], dtype=float)
    depth = np.array([20, 100, 100, 10000, 10000])
    exclusion = depth[:, None]-inclusion
    quantiles = []
    for points in (256, 512):
        draws = hierarchical_draws(inclusion, exclusion, depth, rho=.02,
            rng=np.random.default_rng(19), draws=20000, points=points)
        quantiles.append(np.quantile(draws/depth[:, None], [.025, .975], axis=1))
    assert np.max(np.abs(quantiles[0]-quantiles[1])) < .02
    low, high = quantiles[0]
    observed = np.array([0, .02, .6, .55, 1])
    for ratio in (.25, 1, 4):
        corrected = length_correct(observed, ratio, 1.)
        lower = length_correct(low, ratio, 1.)
        upper = length_correct(high, ratio, 1.)
        np.testing.assert_array_equal((observed >= low)&(observed <= high),
                                      (corrected >= lower)&(corrected <= upper))


def test_hierarchical_holdout_has_no_outcome_leakage():
    rng = np.random.default_rng(19)
    counts = rng.binomial(50, rng.beta(3, 3, (40, 3)))
    frame = pd.DataFrame({'event_id':[f'e{i}' for i in range(40)], 'gene_id':[f'g{i}' for i in range(40)],
        'IJC_SAMPLE_1':[','.join(map(str, row)) for row in counts],
        'SJC_SAMPLE_1':[','.join(map(str, row)) for row in 50-counts],
        'IncFormLen':200., 'SkipFormLen':100.})
    kwargs = dict(sample_ids=['a','b','c'], biological_ids=['a','b','c'], index=1,
                  draws=100, model='hierarchical')
    first, fits = heldout_arm(frame, **kwargs)
    frame['IJC_SAMPLE_1'] = frame.IJC_SAMPLE_1.map(lambda s:'49,'+s.split(',',1)[1])
    frame['SJC_SAMPLE_1'] = frame.SJC_SAMPLE_1.map(lambda s:'1,'+s.split(',',1)[1])
    second, changed_fits = heldout_arm(frame, **kwargs)
    pd.testing.assert_frame_equal(fits[fits.held_sample_id.eq('a')], changed_fits[changed_fits.held_sample_id.eq('a')])
    columns = ['event_id','rho','predictive_lower','predictive_upper']
    pd.testing.assert_frame_equal(first.loc[first.held_sample_id.eq('a'), columns],
                                  second.loc[second.held_sample_id.eq('a'), columns])
    assert set(first.method) == {'training_gene_crossfit_hierarchical'}
    with pytest.raises(ValueError, match='Unknown'):
        heldout_arm(frame, **{**kwargs, 'model':'unknown'})


def test_invalid_and_unpaired_counts_are_not_used():
    kwargs = dict(rho=.02, rng=np.random.default_rng(5), draws=100)
    with pytest.raises(ValueError, match='integers'):
        hierarchical_draws([[1.5,2]], [[1,2]], [10], **kwargs)
    with pytest.raises(ValueError, match='two observed'):
        hierarchical_draws([[1,2]], [[np.nan,2]], [10], **kwargs)
    with pytest.raises(ValueError, match='positive integer'):
        hierarchical_draws([[1,2]], [[1,2]], [2.2], **kwargs)


def test_failed_dispersion_retains_all_evaluable_event_denominators():
    frame = pd.DataFrame({'event_id':['e'], 'gene_id':['g'], 'IJC_SAMPLE_1':['10,10,10'],
        'SJC_SAMPLE_1':['10,10,10'], 'IncFormLen':[200], 'SkipFormLen':[100]})
    detail, _ = heldout_arm(frame, ['a','b','c'], ['a','b','c'], index=1, model='hierarchical', draws=100)
    assert len(detail) == 3
    assert detail.covered.isna().all()
    assert detail.reason.eq('insufficient_training_or_dispersion_failure').all()
