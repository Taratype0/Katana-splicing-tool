"""Numerical roundoff is accepted; scientifically invalid TPM is rejected."""
import numpy as np
import pandas as pd
import pytest
from katana_next.flow.stages import _maei_tables


def inputs(tmp_path, included, total):
    samples = ['A1', 'A2', 'B1', 'B2']
    table = pd.DataFrame([[included]*4 + [total]*4], index=['event1'],
                         columns=[s+'_included' for s in samples]+[s+'_total' for s in samples])
    table.to_csv(tmp_path/'inclusion.tsv', sep='\t')
    pd.DataFrame([[1]], index=['event1'], columns=['AAA']).to_csv(tmp_path/'sites.tsv', sep='\t')
    pd.DataFrame({'sample':samples}).to_csv(tmp_path/'design.tsv', sep='\t', index=False)
    return dict(inclusion_table=tmp_path/'inclusion.tsv', sitecount_matrix=tmp_path/'sites.tsv',
                design=tmp_path/'design.tsv', quantity_provenance='verified transcript sums')


def test_roundoff_private_normalization_preserves_source(tmp_path):
    spec=inputs(tmp_path, 13.390504000000002, 13.390504)
    before=spec['inclusion_table'].read_bytes()
    table, _, _ = _maei_tables(spec)
    assert np.all(table.iloc[:,:4].to_numpy() <= table.iloc[:,4:].to_numpy())
    assert table.attrs['tpm_roundoff_cells'] == 4
    assert 0 < table.attrs['tpm_roundoff_max'] < 1e-12
    assert spec['inclusion_table'].read_bytes() == before


@pytest.mark.parametrize('included,total', [(1.001,1.), (1.,0.), (1.,-1.)])
def test_material_invalid_quantities_still_fail(tmp_path,included,total):
    with pytest.raises(ValueError):
        _maei_tables(inputs(tmp_path,included,total))


def test_valid_values_unchanged(tmp_path):
    table, _, _ = _maei_tables(inputs(tmp_path,2.,3.))
    assert table.attrs['tpm_roundoff_cells'] == 0
    assert table.iloc[0].tolist() == [2.]*4+[3.]*4
