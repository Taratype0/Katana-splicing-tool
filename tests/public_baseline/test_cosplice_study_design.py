"""Independent-unit accounting and explicit pairing across public entry points."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from katana_cosplice.study_design import validate_study_design, project_study_design


def samples(paired=True, n=3):
    return pd.DataFrame([{'sample_id': f'{condition}{i}',
                          'biological_sample_id': str(i) if paired else f'{condition}{i}',
                          'condition': condition}
                         for condition in ('A', 'B') for i in range(n)])


@pytest.mark.parametrize('n', [1, 3, 6, 12, 15])
def test_no_minimum_study_size_gate(n):
    policy = validate_study_design(samples(False, n))
    assert policy['biological_units'] == 2*n
    assert policy['minimum_units_gate'] is False
    assert policy['recommended_biological_units'] == 12


def test_repeated_observations_count_once_and_pairing_is_optional():
    with pytest.raises(ValueError, match='explicit paired'):
        validate_study_design(samples())
    policy = validate_study_design(samples(), paired=True)
    assert policy['samples'] == 6
    assert policy['biological_units'] == 3
    assert policy['below_recommended_units'] is True
    assert validate_study_design(samples(False))['paired'] is False


def test_technical_duplicates_do_not_inflate_units():
    frame = samples()
    frame.loc[len(frame)] = ['Aextra', '0', 'A']
    with pytest.raises(ValueError, match='technical replicates'):
        validate_study_design(frame, paired=True)


def project(paired=True, declare=True):
    frame = samples(paired)
    arm = lambda c: SimpleNamespace(sample_ids=tuple(frame.loc[frame.condition.eq(c), 'sample_id']))
    return SimpleNamespace(samples=[SimpleNamespace(**r) for r in frame.to_dict('records')],
        contrasts=[SimpleNamespace(contrast_id='B_vs_A', numerator=arm('B'), denominator=arm('A'))],
        metadata={'comparison_designs': {'B_vs_A': {'paired': paired}}} if declare else {})


def test_project_pairing_retains_direction_and_counts():
    policy, rows = project_study_design(project())
    assert policy['biological_units'] == 3
    assert rows == [{'contrast_id':'B_vs_A','factor':'condition','a':'B','b':'A','paired':True}]
    assert project_study_design(project(False))[1][0]['paired'] is False
    with pytest.raises(ValueError, match='explicit paired'):
        project_study_design(project(declare=False))
    assert project_study_design(project(declare=False), paired=True)[1][0]['paired'] is True


def test_incomplete_pairing_is_rejected_before_computation():
    p = project()
    p.samples[-1].biological_sample_id = 'unmatched'
    with pytest.raises(ValueError, match='one observation per unit'):
        project_study_design(p)


def test_core_matrix_preparation_has_no_hidden_fifteen_unit_gate():
    from katana_cosplice import _study_core
    path = Path(_study_core.__file__).parent / 'science/cosplice_research.py'
    spec = importlib.util.spec_from_file_location('test_cosplice_research', path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    frame = samples()
    matrix = pd.DataFrame(np.random.default_rng(17).uniform(size=(50,6)), columns=frame.sample_id)
    *_, meta = module.prepare(matrix, frame)
    assert meta.biological_sample_id.nunique() == 3
