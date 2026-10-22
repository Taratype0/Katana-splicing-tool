import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest


def module(name):
    path = Path(__file__).resolve().parents[2] / 'scripts' / name
    spec = importlib.util.spec_from_file_location(name.removesuffix('.py'), path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def test_bundle_rejects_modified_input(tmp_path):
    worker = module('motiftrace_native_resource_worker_v3_20261016.py')
    source = tmp_path / 'input.tsv'
    source.write_text('original\n')
    manifest = {'files': {'input.tsv': {'bytes': source.stat().st_size, 'SHA256': worker.sha(source)}}}
    (tmp_path / 'MANIFEST.json').write_text(json.dumps(manifest))
    assert worker.verify_bundle(tmp_path) == manifest
    source.write_text('changed!\n')
    with pytest.raises(AssertionError):
        worker.verify_bundle(tmp_path)


def test_capacity_uses_hard_limit_not_zero_soft_quota(monkeypatch):
    worker = module('motiftrace_native_resource_worker_v3_20261016.py')
    text = 'quota header\n  4631321812 0 5242880000 - 980897 0 1000000 -\n'
    monkeypatch.setattr(worker.subprocess, 'run', lambda *a, **k: SimpleNamespace(returncode=0, stdout=text, stderr=''))
    actual = worker.capacity()
    assert actual['free_bytes'] == (5242880000 - 4631321812) * 1024
    assert actual['freefiles'] == 19103


def test_capacity_rejects_inode_budget(monkeypatch):
    worker = module('motiftrace_native_resource_worker_v3_20261016.py')
    text = 'quota header\n  4631321812 0 5242880000 - 990000 0 1000000 -\n'
    monkeypatch.setattr(worker.subprocess, 'run', lambda *a, **k: SimpleNamespace(returncode=0, stdout=text, stderr=''))
    with pytest.raises(AssertionError, match='inode'):
        worker.capacity()


def test_finalizer_rejects_partial_success_mapping(tmp_path, monkeypatch):
    final = module('finalize_native_resource_v3_20261016.py')
    (tmp_path / 'BULK_SUBMITTED.json').write_text(json.dumps({'job_id':'123'}))
    monkeypatch.setattr(final.subprocess, 'run', lambda *a, **k: SimpleNamespace(returncode=0, stdout='123_0|456|COMPLETED|0:0|\n'))
    with pytest.raises(AssertionError, match='All500'):
        final.finalize(tmp_path)
    assert not (tmp_path / 'ALL_OWNED_NATIVE_SHARDS_ACCEPTED.json').exists()


def test_finalizer_rejects_failed_allocation(tmp_path, monkeypatch):
    final = module('finalize_native_resource_v3_20261016.py')
    (tmp_path / 'BULK_SUBMITTED.json').write_text(json.dumps({'job_id':'123'}))
    rows = ''.join('123_%d|%d|%s|%s|\n' % (i, 456+i, 'FAILED' if i==499 else 'COMPLETED', '1:0' if i==499 else '0:0') for i in range(500))
    monkeypatch.setattr(final.subprocess, 'run', lambda *a, **k: SimpleNamespace(returncode=0, stdout=rows))
    with pytest.raises(AssertionError, match='Not every'):
        final.finalize(tmp_path)
