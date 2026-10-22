import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from katana_cosplice.qc import _validate_inputs, stable_as_qc
from katana_cosplice.config import StableASConfig
from katana_cosplice.matrix import _validate_matrix_values


class RangeValidationTests(unittest.TestCase):
    def inputs(self, values):
        psi = pd.DataFrame([values], index=['event'], columns=[f's{i}' for i in range(len(values))])
        counts = pd.DataFrame(20., index=psi.index, columns=psi.columns)
        samples = pd.DataFrame({'sample_id': psi.columns, 'condition': 'A'})
        return psi, counts, counts.copy(), samples

    def test_boundary_values_and_missing_measurements_are_valid(self):
        _validate_inputs(*self.inputs([0., 1., .5, np.nan]))

    def test_genuine_out_of_range_values_and_infinities_are_rejected(self):
        for value in [np.nextafter(0., -1.), np.nextafter(1., 2.), -np.inf, np.inf]:
            with self.subTest(value=value):
                result = stable_as_qc(*self.inputs([value, .5]), StableASConfig())
                self.assertFalse(result.metrics.loc['event', 'stable_as'])
                self.assertIn('invalid_or_missing_psi', result.metrics.loc['event', 'failed_rules'])

    def test_negative_read_counts_are_rejected(self):
        for position in (1, 2):
            inputs = self.inputs([0.5, .6])
            inputs[position].iloc[0, 0] = -1
            with self.subTest(position=position):
                result = stable_as_qc(*inputs, StableASConfig())
                self.assertFalse(result.metrics.loc['event', 'stable_as'])
                self.assertIn('invalid_or_missing_counts', result.metrics.loc['event', 'failed_rules'])

    def test_validator_does_not_depend_on_stacked_series_comparison(self):
        with patch.object(pd.DataFrame, 'stack', side_effect=AssertionError('unneeded stacked index')):
            _validate_inputs(*self.inputs([0., 1., np.nan]))

    def test_event_and_sample_order_checks_still_apply(self):
        inputs = self.inputs([.2, .8])
        reordered = inputs[1].iloc[:, ::-1]
        with self.assertRaisesRegex(ValueError, 'same ordered sample columns'):
            _validate_inputs(inputs[0], reordered, inputs[2], inputs[3])

    def test_import_validator_uses_the_same_exact_bounds(self):
        _validate_matrix_values(*self.inputs([0., 1., np.nan])[:3])
        for value in [np.nextafter(0., -1.), np.nextafter(1., 2.), -np.inf, np.inf]:
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, 'IncLevel values'):
                _validate_matrix_values(*self.inputs([value])[:3])
        for position in (1, 2):
            inputs = self.inputs([.5])
            inputs[position].iloc[0, 0] = -1.
            with self.subTest(position=position), self.assertRaisesRegex(ValueError, 'cannot be negative'):
                _validate_matrix_values(*inputs[:3])

    def test_import_validator_keeps_duplicate_identity_checks(self):
        inputs = self.inputs([.2, .8])
        duplicated = inputs[0].copy()
        duplicated.columns = ['s0', 's0']
        with self.assertRaisesRegex(ValueError, 'duplicate sample columns'):
            _validate_matrix_values(duplicated, *inputs[1:3])
        duplicated = pd.concat([inputs[0], inputs[0]])
        with self.assertRaisesRegex(ValueError, 'duplicate event_key'):
            _validate_matrix_values(duplicated, *inputs[1:3])

    def test_both_validators_avoid_the_failing_stacked_comparison(self):
        with patch.object(pd.DataFrame, 'stack', side_effect=AssertionError('unneeded stacked index')):
            inputs = self.inputs([0., 1., np.nan])
            _validate_inputs(*inputs)
            _validate_matrix_values(*inputs[:3])


if __name__ == '__main__':
    unittest.main()
