from dataclasses import asdict
from types import SimpleNamespace
import json
import pytest
from katana_cosplice.config import NetworkConfig
from katana_next.workflows.hub_events import source_network

def test_reads_actual_psi_discovery_configuration_without_changing_power_or_minimum(tmp_path):
    path=tmp_path/'resolved_psi_configuration.json'
    path.write_text(json.dumps({'network':asdict(NetworkConfig(soft_power=9,minimum_network_samples=25))}))
    actual=source_network(SimpleNamespace(path=path))
    assert actual['power']==9 and actual['min_network_independent_samples']==25
    assert actual['method']=='bicor' and actual['signed'] is True

def test_legacy_configuration_is_supported_and_conflicting_power_is_rejected(tmp_path):
    path=tmp_path/'resolved_cosplice.yaml'
    path.write_text('network: {method: pearson, power: 6, signed: true}\n')
    assert source_network(SimpleNamespace(path=path))['power']==6
    path.write_text('network: {method: pearson, power: 6, soft_power: 9, signed: true}\n')
    with pytest.raises(ValueError,match='conflicting power'):
        source_network(SimpleNamespace(path=path))
