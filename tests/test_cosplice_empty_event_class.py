import pandas as pd
import pytest

from katana_cosplice.matrix import EVENT_COORDINATES, _event_frame, _vector_frame


@pytest.mark.parametrize('event_type', list(EVENT_COORDINATES))
def test_empty_class_retains_identity_schema_and_zero_measurements(event_type):
    columns = ['GeneID', 'geneSymbol', 'chr', 'strand', *EVENT_COORDINATES[event_type]]
    events = _event_frame(pd.DataFrame(columns=columns), event_type)
    assert events.empty
    assert events.columns.tolist() == [
        'event_key', 'gene_id', 'gene_symbol', 'event_type', 'chr', 'strand',
        'event_coordinates', *EVENT_COORDINATES[event_type],
    ]
    matrix = _vector_frame(pd.DataFrame(columns=['IncLevel1']), events.event_key,
                           ['mouse1', 'mouse2'], 'IncLevel1')
    assert matrix.shape == (0, 2)
    assert matrix.columns.tolist() == ['mouse1', 'mouse2']


def test_empty_table_with_missing_identity_still_rejected():
    with pytest.raises(ValueError, match='missing event identity columns'):
        _event_frame(pd.DataFrame(columns=['GeneID']), 'SE')


def test_nonempty_event_identity_is_unchanged():
    frame = pd.DataFrame([{'GeneID': 'ENSG1.4', 'geneSymbol': 'Gene1', 'chr': 'chr1',
                           'strand': '+', **{name:100+i for i,name in enumerate(EVENT_COORDINATES['SE'])}}])
    result = _event_frame(frame, 'SE')
    assert len(result) == 1
    assert result.event_key.iloc[0] == 'ENSG1|SE|chr1|+|100|101|102|103|104|105'
    assert result.gene_symbol.iloc[0] == 'Gene1'
