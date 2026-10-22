import pickle

import numpy as np
import pandas as pd
import pytest

from scripts import finish_public5712_io_compat_v5_20261016 as compat


def test_reads_legacy_text_and_retains_numeric_bits(tmp_path):
    path = tmp_path / 'legacy.npz'
    p = np.array([2.**-235, .3, 1.], dtype=np.float64)
    n = np.array([235, 21, 0], dtype=np.int64)
    np.savez_compressed(path, method=np.array(['count', 'coverage', 'PWM'], dtype=object), P_exact=p, n=n)
    frame = compat.load_arrays(path)
    assert frame.method.tolist() == ['count', 'coverage', 'PWM']
    assert frame.P_exact.to_numpy().tobytes() == p.tobytes()
    assert frame.n.to_numpy().tobytes() == n.tobytes()


def test_object_numeric_arrays_are_rejected_not_coerced(tmp_path):
    path = tmp_path / 'numeric.npz'
    np.savez_compressed(path, P_exact=np.array([.1, .5], dtype=object))
    with pytest.raises(ValueError, match='registered text'):
        compat.load_arrays(path)


def test_non_string_objects_in_label_are_rejected(tmp_path):
    path = tmp_path / 'mixed.npz'
    np.savez_compressed(path, method=np.array(['count', 123], dtype=object))
    with pytest.raises(ValueError, match='non-string'):
        compat.load_arrays(path)


def test_unexpected_pickle_global_is_never_executed(tmp_path):
    path = tmp_path / 'unexpected.npz'
    np.savez_compressed(path, method=np.array([complex(1, 2)], dtype=object))
    with pytest.raises(pickle.UnpicklingError, match='reconstruction global'):
        compat.load_arrays(path)


def test_extension_string_columns_are_saved_without_pickle(tmp_path):
    frame = pd.DataFrame({'method': pd.Series(['count', 'PWM'], dtype='string'),
                          'P_exact': [2.**-40, 1.]})
    path = tmp_path / 'portable.npz'
    compat.write_arrays(path, frame)
    with np.load(path, allow_pickle=False) as arrays:
        assert arrays['method'].dtype.kind == 'U'
        assert arrays['method'].tolist() == ['count', 'PWM']
        assert arrays['P_exact'].tobytes() == frame.P_exact.to_numpy().tobytes()
