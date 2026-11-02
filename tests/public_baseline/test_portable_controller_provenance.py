import hashlib
import importlib.util
import json
from pathlib import Path

import pytest


def test_portable_controller_omits_build_path_but_keeps_distribution_identity(tmp_path):
    path = Path(__file__).resolve().parents[2] / 'tools/runtime/build_portable_controller.py'
    spec = importlib.util.spec_from_file_location('portable_controller_build', path)
    builder = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(builder)
    wheel = tmp_path / 'katana_test-1-py3-none-any.whl'
    wheel.write_bytes(b'accepted fixture wheel')
    digest = hashlib.sha256(wheel.read_bytes()).hexdigest()
    remote = {'download_info': {'url': 'https://example.org/dependency.whl',
                               'archive_info': {'hashes': {'sha256': 'public-digest'}}}}
    data = {'install': [{'download_info': {'url': wheel.as_uri(),
              'archive_info': {'hashes': {'sha256': digest}}}}, remote]}
    report = tmp_path / 'installation-report.json'
    report.write_text(json.dumps(data))
    builder.normalize_installation_report(report, wheel)
    actual = json.loads(report.read_text())
    assert actual['install'][0]['download_info'] == {
        'url': 'artifact:' + wheel.name, 'local_path_omitted': True,
        'archive_info': {'hashes': {'sha256': digest}}}
    assert actual['install'][1] == remote
    data['install'][0]['download_info']['archive_info']['hashes']['sha256'] = 'wrong'
    report.write_text(json.dumps(data))
    with pytest.raises(ValueError, match='hash differs'):
        builder.normalize_installation_report(report, wheel)
