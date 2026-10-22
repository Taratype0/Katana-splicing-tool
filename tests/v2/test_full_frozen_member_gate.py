import numpy as np
import pandas as pd
import pytest

from katana_next.science.frozen_event_transport import full_frozen_member_gate


def inputs():
    samples = pd.DataFrame({'sample_id': [f's{i}' for i in range(18)],
                            'patient_id': [f'p{i}' for i in range(18)]})
    keys = [f'gene|ES|{i}' for i in range(20)]
    matrix = pd.DataFrame(np.tile(np.linspace(.1,.9,18), (20,1)),
                          index=keys, columns=samples.sample_id)
    return matrix, keys, samples, dict.fromkeys(keys, True)


def test_full_measured_module_has_no_implicit_statistical_result():
    result = full_frozen_member_gate(*inputs())
    assert result['full_module_input_eligible']
    assert result['original_members'] == result['measured_members'] == 20
    assert result['complete_case_patients'] == 18
    assert result['preservation_pvalue'] is result['preservation_fdr'] is None
    assert result['imputed_cells'] == result['dropped_original_members'] == 0


def test_missing_original_event_does_not_shrink_the_module():
    matrix, keys, samples, identity = inputs()
    result = full_frozen_member_gate(matrix.iloc[:-1], keys, samples, identity)
    assert not result['full_module_input_eligible']
    assert result['original_members'] == 20
    assert result['absent_original_events'] == [keys[-1]]
    assert result['complete_case_patients'] == 0
    assert result['dropped_original_members'] == 0


def test_finite_but_gene_only_match_is_not_exact_event_identity():
    matrix, keys, samples, identity = inputs()
    identity[keys[0]] = False
    result = full_frozen_member_gate(matrix, keys, samples, identity)
    assert not result['full_module_input_eligible']
    assert result['exact_identity_members'] == 19


def test_same_member_complete_cases_apply_to_every_original_node():
    matrix, keys, samples, identity = inputs()
    matrix.loc[keys[-1], samples.sample_id[:4]] = np.nan
    result = full_frozen_member_gate(matrix, keys, samples, identity)
    assert result['complete_case_patients'] == 14
    assert not result['full_module_input_eligible']
    assert result['measured_members'] == 20
    assert 'insufficient_full_module_complete_patients' in result['reasons']


def test_partial_missingness_can_remove_patients_but_never_members():
    matrix, keys, samples, identity = inputs()
    matrix.iloc[0,0] = np.nan
    result = full_frozen_member_gate(matrix, keys, samples, identity)
    assert result['full_module_input_eligible']
    assert result['complete_case_patients'] == 17
    assert result['original_members'] == 20
    assert result['imputed_cells'] == result['dropped_original_members'] == 0


def test_constant_member_cannot_disappear_from_full_family():
    matrix, keys, samples, identity = inputs()
    matrix.iloc[0,:] = .5
    result = full_frozen_member_gate(matrix, keys, samples, identity)
    assert not result['full_module_input_eligible']
    assert result['variable_members_in_complete_cases'] == 19
    assert result['original_members'] == 20


def test_technical_duplicates_are_not_independent_patients():
    matrix, keys, samples, identity = inputs()
    samples.loc[1,'patient_id'] = samples.loc[0,'patient_id']
    with pytest.raises(ValueError, match='One biological observation'):
        full_frozen_member_gate(matrix, keys, samples, identity)


@pytest.mark.parametrize('value', [-.01,1.01])
def test_out_of_range_measured_psi_is_rejected(value):
    matrix, keys, samples, identity = inputs()
    matrix.iloc[0,0] = value
    with pytest.raises(ValueError, match='outside'):
        full_frozen_member_gate(matrix, keys, samples, identity)


def test_duplicate_frozen_membership_is_rejected():
    matrix, keys, samples, identity = inputs()
    with pytest.raises(ValueError, match='unique'):
        full_frozen_member_gate(matrix, keys+[keys[0]], samples, identity)


def test_numpy_boolean_exact_identity_is_accepted():
    matrix, keys, samples, identity = inputs()
    identity = dict.fromkeys(keys, np.bool_(True))
    assert full_frozen_member_gate(matrix, keys, samples, identity)['full_module_input_eligible']
