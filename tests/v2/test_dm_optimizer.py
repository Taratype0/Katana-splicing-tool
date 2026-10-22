import numpy as np
from scipy.optimize._numdiff import approx_derivative
from katana_next.science.statistics import _dm_value_gradient, dirichlet_multinomial_two_group


def test_dm_analytic_gradient_matches_finite_difference():
    rng = np.random.default_rng(18)
    counts = rng.multinomial(40, [.1, .3, .6], size=16).astype(float)
    design = np.column_stack([np.ones(16), np.repeat([0, 1], 8)])
    parameters = np.array([-1., .2, .4, -.5, 2.])
    _, analytic = _dm_value_gradient(parameters, design, counts)
    numeric = approx_derivative(lambda x: _dm_value_gradient(x, design, counts)[0], parameters).ravel()
    np.testing.assert_allclose(analytic, numeric, rtol=1e-5, atol=1e-7)


def test_dm_fit_is_nested_symmetric_and_count_validated():
    rng = np.random.default_rng(22)
    counts = np.vstack([rng.multinomial(50, [.7, .2, .1], size=20),
                        rng.multinomial(50, [.2, .3, .5], size=20)])
    groups = np.repeat([0, 1], 20)
    result = dirichlet_multinomial_two_group(counts, groups)
    reverse = dirichlet_multinomial_two_group(counts, 1-groups)
    assert result["status"] == reverse["status"] == "fit"
    assert result["full_log_likelihood"] >= result["null_log_likelihood"] - 1e-6
    np.testing.assert_allclose(result["likelihood_ratio"], reverse["likelihood_ratio"], rtol=1e-5)
    np.testing.assert_allclose(result["delta_proportions"], -reverse["delta_proportions"], atol=1e-5)
    assert dirichlet_multinomial_two_group(counts + .1, groups)["status"] == "invalid_counts"
