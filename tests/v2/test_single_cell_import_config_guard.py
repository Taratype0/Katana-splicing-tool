"""Isolated config/API regressions; no calls to the current API process."""
import json
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from katana_next.api.app import PlatformRuntime, create_app
from katana_next.full_fixture import create_full_fixture
from katana_next.science.single_cell_imports import read_config


@pytest.mark.parametrize('value', ['', '   ', '\t\n', None, 0, False, [], {}])
def test_empty_or_non_path_config_is_validation_error(value):
    with pytest.raises(ValueError, match='config_path is required'):
        read_config(value)


@pytest.mark.parametrize('kind', ['missing', 'directory', 'empty', 'invalid_json', 'bad_encoding', 'array', 'null', 'number', 'boolean', 'bad_path_type', 'whitespace_source'])
def test_missing_nonfile_and_invalid_config(tmp_path, kind):
    path = tmp_path / 'config.json'
    if kind == 'directory':
        path.mkdir()
    elif kind == 'bad_encoding':
        path.write_bytes(b'\xff\xfe\xfa')
    elif kind != 'missing':
        payload = {'empty': '', 'invalid_json': '{not json', 'array': '[]', 'null': 'null',
                   'number': '42', 'boolean': 'true', 'bad_path_type': '{"metadata_path": ["a"]}',
                   'whitespace_source': '{"psi_path": "   "}'}[kind]
        path.write_text(payload, encoding='utf-8')
    with pytest.raises(ValueError):
        read_config(path)


def test_permission_failure_is_actionable_value_error(tmp_path, monkeypatch):
    path = tmp_path / 'config.json'
    path.write_text('{}', encoding='utf-8')
    original = Path.read_text

    def denied(self, *args, **kwargs):
        if self == path:
            raise PermissionError('test-only denied read')
        return original(self, *args, **kwargs)

    monkeypatch.setattr(Path, 'read_text', denied)
    with pytest.raises(ValueError, match='Cannot read import configuration'):
        read_config(path)


def test_bom_relative_paths_and_source_values_preserved(tmp_path):
    path = tmp_path / 'config.json'
    config = {'mode': 'marvel', 'metadata_path': 'meta.tsv', 'results_path': 'results.tsv',
              'psi_path': None, 'method': 'permutation', 'adjustment': 'none'}
    path.write_text(json.dumps(config), encoding='utf-8-sig')
    observed = read_config(path)
    assert observed['metadata_path'] == str((tmp_path / 'meta.tsv').resolve())
    assert observed['results_path'] == str((tmp_path / 'results.tsv').resolve())
    assert observed['psi_path'] is None
    assert observed['method'] == 'permutation' and observed['adjustment'] == 'none'


@pytest.fixture(scope='module')
def isolated_api(tmp_path_factory):
    root = tmp_path_factory.mktemp('isolated-api')
    manifest = create_full_fixture(root / 'project')
    runtime = PlatformRuntime(application_root=root / 'application')
    with TestClient(create_app(runtime)) as client:
        response = client.post('/api/projects/open', json={'manifest_path': str(manifest)})
        assert response.status_code == 200, response.text
        yield client, runtime, response.json()['snapshot']['project_id']
    runtime.jobs._executor.shutdown(wait=True)


def good_config(tmp_path, mode, project_id):
    paths = ['metadata_path', 'results_path'] if mode == 'marvel' else ['metadata_path', 'counts_path', 'annotation_path']
    config = {'mode': mode, 'project_id': project_id}
    for key in paths:
        file = tmp_path / (key + '.tsv')
        file.write_text('source\n', encoding='utf-8')
        config[key] = str(file)
    return config


@pytest.mark.parametrize('mode', ['marvel', 'longread_isoforms'])
@pytest.mark.parametrize('kind', ['blank', 'whitespace', 'omitted', 'none', 'missing', 'directory', 'invalid_json', 'array', 'bad_encoding', 'missing_required', 'input_directory', 'input_missing', 'wrong_mode'])
def test_bad_config_preflight_and_cache_no500_or_fake_cache(tmp_path, isolated_api, mode, kind, monkeypatch):
    client, runtime, project_id = isolated_api
    path = tmp_path / 'config.json'
    config = good_config(tmp_path, mode, project_id)
    if kind == 'directory':
        path.mkdir()
    elif kind == 'invalid_json':
        path.write_text('{invalid', encoding='utf-8')
    elif kind == 'array':
        path.write_text('[]', encoding='utf-8')
    elif kind == 'bad_encoding':
        path.write_bytes(b'\xff\xfe\xfa')
    elif kind != 'missing':
        if kind == 'missing_required':
            config.pop('metadata_path')
        elif kind == 'input_directory':
            config['metadata_path'] = str(tmp_path)
        elif kind == 'input_missing':
            config['metadata_path'] = str(tmp_path / 'not-created.tsv')
        elif kind == 'wrong_mode':
            config['mode'] = 'unsupported'
        path.write_text(json.dumps(config), encoding='utf-8')
    value = {'blank': '', 'whitespace': '  \t', 'none': None}.get(kind, str(path))
    parameters = {} if kind == 'omitted' else {'config_path': value}
    endpoint = f'/api/projects/{project_id}/workflows/single_cell.import_{mode}'
    before = len(runtime.jobs.list(project_id))

    def no_fingerprint(*args, **kwargs):
        raise AssertionError('Invalid config reached fingerprint/cache resolution')

    monkeypatch.setattr(runtime.workflows, '_fingerprint', no_fingerprint)
    request = {'parameters': parameters}
    plan = client.post(endpoint + '/plan', json=request)
    assert plan.status_code == 200, plan.text
    assert plan.json()['ready'] is False and plan.json()['errors']
    assert plan.json()['warnings']
    cache = client.post(endpoint + '/cache', json=request)
    assert cache.status_code == 422, cache.text
    assert cache.json()['detail'] in plan.json()['errors']
    assert 'cache_hit' not in cache.json() and 'run_fingerprint' not in cache.json()
    run = client.post(endpoint + '/run', json=request)
    assert run.status_code == 422, run.text
    assert len(runtime.jobs.list(project_id)) == before


@pytest.mark.parametrize('mode', ['marvel', 'longread_isoforms'])
def test_valid_config_reaches_genuine_fingerprint_and_invalidation(tmp_path, isolated_api, mode, monkeypatch):
    client, runtime, project_id = isolated_api
    config = good_config(tmp_path, mode, project_id)
    path = tmp_path / 'config.json'
    path.write_text(json.dumps(config), encoding='utf-8')
    endpoint = f'/api/projects/{project_id}/workflows/single_cell.import_{mode}'
    request = {'parameters': {'config_path': str(path)}}
    calls = []
    original = runtime.workflows._fingerprint

    def spy(file):
        calls.append(file)
        return original(file)

    monkeypatch.setattr(runtime.workflows, '_fingerprint', spy)
    plan = client.post(endpoint + '/plan', json=request)
    assert plan.status_code == 200 and plan.json()['ready'] is True
    first = client.post(endpoint + '/cache', json=request)
    assert first.status_code == 200 and first.json()['cache_hit'] is False
    assert first.json()['run_fingerprint'] and calls
    Path(config['metadata_path']).write_text('changed-source\n', encoding='utf-8')
    changed = client.post(endpoint + '/cache', json=request)
    assert changed.status_code == 200
    assert changed.json()['run_fingerprint'] != first.json()['run_fingerprint']
