import numpy as np
import pytest
from scipy.stats import rankdata

from katana_cosplice import dynamic_tree


def test_ordinal_keeps_exact_order_and_values_as_valid_integer_indices():
    values = np.array([7.0, 2.0, 2.0, 9.0])
    expected = rankdata(values, method="ordinal")
    actual = dynamic_tree._rankdata_for_dynamic_tree_cut(values, method="ordinal")
    np.testing.assert_array_equal(actual, expected)
    assert actual.dtype == np.dtype(np.intp)
    np.testing.assert_array_equal(np.array([10, 20, 30, 40])[actual - 1], [30, 10, 20, 40])
    assert dynamic_tree._dynamic_tree_module.rankdata is dynamic_tree._rankdata_for_dynamic_tree_cut


@pytest.mark.parametrize("method", ["average", "min", "max", "dense"])
def test_nonordinal_methods_return_original_object_unmodified(monkeypatch, method):
    original = np.array([1.5, 1.5, 3.0])
    monkeypatch.setattr(dynamic_tree, "_scipy_rankdata", lambda *args, **kwargs: original)
    assert dynamic_tree._rankdata_for_dynamic_tree_cut([1, 1, 3], method=method) is original


@pytest.mark.parametrize("invalid", [[1.5, 2], [np.nan, 2], [np.inf, 2], [0, 2], [1, 3]])
def test_invalid_ordinal_indices_are_rejected_not_truncated(monkeypatch, invalid):
    monkeypatch.setattr(dynamic_tree, "_scipy_rankdata", lambda *args, **kwargs: np.array(invalid))
    with pytest.raises(ValueError, match="invalid ordinal"):
        dynamic_tree._rankdata_for_dynamic_tree_cut([1, 2], method="ordinal")
