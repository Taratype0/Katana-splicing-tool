import itertools

import numpy as np
from scipy.stats import false_discovery_control, rankdata, wilcoxon

from katana_next.science.motiftrace_exact import ExactSignedRank, full_family_adjust


def test_tied_zero_na_matches_exhaustive_signs():
    engine = ExactSignedRank()
    rng = np.random.default_rng(20261016)
    for n in range(1, 17):
        d = rng.choice([-3., -1., 0., 1., 2., 3.], n)
        nonzero = d[d != 0]
        weights = 2 * rankdata(np.abs(nonzero))
        observed = weights[nonzero > 0].sum()
        signs = np.array(list(itertools.product((0, 1), repeat=len(nonzero))))
        expected = np.mean(signs @ weights >= observed)
        assert abs(engine.one(np.r_[d, np.nan])[0] - expected) < 2e-14


def test_untied_matches_scipy_and_extreme_tail_not_subtraction():
    engine = ExactSignedRank()
    rng = np.random.default_rng(53)
    for n in (8, 20, 40, 50):
        d = np.arange(1., n + 1) * rng.choice([-1, 1], n)
        assert abs(engine.one(d)[0] - wilcoxon(d, alternative='greater', method='exact').pvalue) < 2e-14
    assert engine.one(np.arange(1., 251))[0] == 2.**-250
    assert engine.one(np.ones(250))[0] == 2.**-250


def test_matrix_denominator_zero_ties_and_na():
    engine = ExactSignedRank()
    d = np.tile([1., 0., -1., np.nan], (25, 1))
    d[:10, 2] = 1
    z, n, p, nonzero, minimum_p = engine.evaluate(d)
    assert n.tolist() == [25, 25, 25, 0]
    assert p[0] == 2.**-25 and p[1] == 1 and np.isnan(p[3])
    assert nonzero.tolist() == [25, 0, 25, 0]
    assert minimum_p[1] == 1
    assert p[2] == engine.one(d[:, 2])[0]
    assert np.isnan(z[3])


def test_adjustments_equal_standard_bh_by_and_holm():
    p = np.random.default_rng(3).uniform(size=2000)
    p[:3] = [1e-30, 1e-10, .001]
    bh, by, holm = full_family_adjust(p)
    np.testing.assert_allclose(bh, false_discovery_control(p, method='bh'))
    np.testing.assert_allclose(by, false_discovery_control(p, method='by'))
    order = np.argsort(p)
    expected = np.maximum.accumulate(p[order] * np.arange(len(p), 0, -1)).clip(0, 1)
    np.testing.assert_allclose(holm[order], expected)


def test_cache_eviction_does_not_change_probability():
    small, big = ExactSignedRank(100), ExactSignedRank()
    for n in range(1, 40):
        d = np.arange(1., n + 1)
        assert small.one(d)[0] == big.one(d)[0]


def test_conditional_null_is_superuniform_with_ties():
    engine = ExactSignedRank()
    for d in (np.ones(30), np.arange(1., 41), np.repeat([1., 2., 3.], [8, 13, 29])):
        weights = np.rint(2 * rankdata(d)).astype(int)
        tail, _ = engine.tail(weights)
        pmf = np.r_[tail[:-1] - tail[1:], tail[-1]]
        for alpha in (.05, .01, 1e-5, 1e-10):
            assert pmf[tail <= alpha].sum() <= alpha + 2e-14
