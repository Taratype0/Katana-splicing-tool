"""Execution settings cannot introduce scientific changes or site assumptions."""
import json
from pathlib import Path

import pytest

from katana_next.flow.resources import render_resources

ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize('options', [{}, {'defaults': {}}, {'labels': {'process_single': {}}, 'processes': {'QC': {}}}])
def test_no_override_is_a_valid_empty_config(options):
    # Nextflow rejects empty process/selector blocks, so emit only the comment.
    assert render_resources(options).startswith('// Generated')
    assert '{' not in render_resources(options)


@pytest.mark.parametrize('options', [
    {'labels': []}, {'labels': {'process_single|.*': {'queue': 'shared'}}},
    {'labels': {'process_single': {'queue': "shared'; params.x=1"}}},
    {'labels': {'process_single': {'executor': 'slurm'}}},
    {'labels': {'process_single': {'cpus': 0}}},
    {'labels': {'process_single': {'fdr': 0.05}}},
])
def test_label_overrides_are_bounded_execution_settings(options):
    with pytest.raises(ValueError):
        render_resources(options)


def test_site_selectors_keep_exact_process_overrides_last():
    text = render_resources({'defaults': {'queue': 'default'},
        'labels': {'katana_downstream': {'queue': 'shared', 'memory': '8 GB'}},
        'processes': {'NETWORK_FIT': {'queue': 'parallel', 'memory': '96 GB'}}})
    assert text.index("queue = 'default'") < text.index('withLabel:') < text.index('withName:')
    assert 'executor' not in text and 'params.' not in text


def test_generic_example_does_not_choose_a_cluster_queue():
    value = json.loads((ROOT / 'examples/nextflow/resources.json').read_text())
    text = render_resources(value)
    assert 'queue =' not in text and 'clusterOptions' not in text


def test_dardel_example_keeps_shared_allocation_unless_explicitly_overridden():
    value = json.loads((ROOT / 'examples/nextflow/resources.dardel.json').read_text())
    text = render_resources(value)
    assert value['defaults']['queue'] == 'shared'
    assert 'main' not in text and 'memory' not in text
    for alias in ('QC', 'RSEQC_INFEREXPERIMENT', 'GFFREAD_STRAND_MODELS'):
        assert value['processes'][alias]['queue'] == 'shared'
