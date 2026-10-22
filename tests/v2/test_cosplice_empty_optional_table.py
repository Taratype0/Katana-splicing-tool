import pandas as pd
import pytest
from katana_next.workflows.programs import _read_engine_table


def test_empty_optional_preservation_is_not_a_result(tmp_path):
    path = tmp_path / 'observed_residual_preservation_pairs.tsv'
    path.write_text('\n')
    assert _read_engine_table(path, '\t') is None


def test_empty_required_table_still_fails(tmp_path):
    path = tmp_path / 'module_assignments.tsv'
    path.write_text('\n')
    with pytest.raises(pd.errors.EmptyDataError):
        _read_engine_table(path, '\t')


def test_header_only_table_retains_schema(tmp_path):
    path = tmp_path / 'observed_residual_preservation_summary.tsv'
    path.write_text('module\tvalue\n')
    result = _read_engine_table(path, '\t')
    assert result.empty
    assert list(result.columns) == ['module', 'value']
