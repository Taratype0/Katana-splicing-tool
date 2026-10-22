import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / 'scripts'
sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location('repeat_report', SCRIPTS / 'build_repeat_supported_reports.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_ordinary_das_never_enters_downstream_selection():
    frame = pd.DataFrame({'event_id': ['ordinary', 'supported', 'negative'],
                          'significant': [True, True, False],
                          'stable_significant': [False, True, False]})
    assert module.repeat_events(frame).event_id.tolist() == ['supported']


def test_inconsistent_source_is_rejected():
    frame = pd.DataFrame({'significant': [False], 'stable_significant': [True]})
    with pytest.raises(AssertionError):
        module.repeat_events(frame)


def test_output_is_isolated_from_prior_delivery():
    assert module.DEST != module.SOURCE


def test_network_coordinate_serialization_matches_rmats_integers():
    frame = pd.DataFrame({'event_type':['SE'], 'chr':['chr1'], 'strand':['+'],
                          'exonStart_0base':[100.0], 'exonEnd':[120.0]})
    assert module.network_event_ids(frame).iloc[0] == 'SE|chr=chr1|strand=+|exonStart_0base=100|exonEnd=120'


def test_coordinates_are_not_rounded():
    frame = pd.DataFrame({'event_type':['SE'], 'exonStart_0base':[100.1]})
    with pytest.raises(ValueError):
        module.network_event_ids(frame)
