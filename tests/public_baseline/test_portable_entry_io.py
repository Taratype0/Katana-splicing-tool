"""Real I/O and HTTP regressions for portable resources and cached verification."""
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import threading
from http.server import ThreadingHTTPServer, BaseHTTPRequestHandler

import pytest
from katana_next.flow import file_io

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("configure_portable", ROOT / "pipeline/bin/katana_configure.py")
configure = importlib.util.module_from_spec(spec)
spec.loader.exec_module(configure)


def test_copy_is_independent_and_source_is_read_once(tmp_path, monkeypatch):
    source = tmp_path / 'source'; target = tmp_path / 'copy'
    source.write_bytes(b'complete result' * 1000)
    expected = hashlib.sha256(source.read_bytes()).hexdigest()
    opened = []
    original = Path.open
    def tracked(self, *args, **kwargs):
        if self in (source, target) and args and args[0] == 'rb': opened.append(self)
        return original(self, *args, **kwargs)
    monkeypatch.setattr(Path, 'open', tracked)
    assert file_io.copy_verified(source, target, expected=expected) == expected
    # Windows identity queries open handles; POSIX has an exact one-pass bound.
    if os.name == 'posix': assert opened == [source]
    assert not os.path.samefile(source, target)
    target.write_bytes(b'changed')
    assert source.stat().st_size == len(b'complete result' * 1000)


def test_persistent_digest_reused_across_processes_and_changed_input_rejected(tmp_path, monkeypatch):
    if os.name == 'nt': pytest.skip('Windows fresh ChangeTime deliberately rehashes')
    cache = tmp_path / 'cache'; source = tmp_path / 'result'
    source.write_bytes(b'abcdefgh')
    monkeypatch.setenv('KATANA_DIGEST_CACHE', str(cache))
    expected = file_io.digest(source)
    code = """from katana_next.flow import file_io
import sys
file_io.hashlib.file_digest = lambda *a: (_ for _ in ()).throw(AssertionError('unnecessary read'))
assert file_io.digest(sys.argv[1]) == sys.argv[2]
"""
    subprocess.run([sys.executable, '-c', code, str(source), expected], check=True)
    stat = source.stat(); source.write_bytes(b'XXXXXXXX'); os.utime(source, ns=(stat.st_atime_ns, stat.st_mtime_ns))
    assert file_io.digest(source) != expected
    assert [p.name for p in cache.iterdir()] == ['digests.sqlite3']


def test_http_download_gzip_cache_and_corruption(tmp_path, monkeypatch):
    payload = b'reference gene sets\n' * 100
    checksum = hashlib.sha256(payload).hexdigest()
    hits = []
    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            hits.append(self.path); self.send_response(200); self.end_headers(); self.wfile.write(gzip.compress(payload))
        def log_message(self, *args): pass
    server = ThreadingHTTPServer(('127.0.0.1', 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True); thread.start()
    item = dict(filename='genes.gmt', compression='gzip', url=f'http://127.0.0.1:{server.server_port}/genes.gz', metadata={'sha256': checksum})
    try:
        path = configure.fetch_resource(item, tmp_path / 'cache')
        assert path.read_bytes() == payload
        assert configure.fetch_resource(item, tmp_path / 'cache') == path
        assert hits == ['/genes.gz']
        path.write_bytes(b'corruption')
        with pytest.raises(ValueError, match='Cached resource changed'): configure.fetch_resource(item, tmp_path / 'cache')
        assert hits == ['/genes.gz']
    finally:
        server.shutdown(); server.server_close(); thread.join()


def test_multiple_libraries_preserve_separate_ids_and_metadata(tmp_path):
    paths = [tmp_path / str(i) for i in range(4)]
    for p in paths: p.write_text('fixture')
    request = {**{'skip_' + k: True for k in configure.ANALYSES}, 'species': 'Mus_musculus',
               'resource_names': ['gene_sets', 'gene_sets', 'clip_peaks', 'clip_peaks'],
               'resource_entries': [{'asset_id': x, 'metadata': {'source_url': 'https://example.org/' + x}} for x in ('GO', 'KEGG', 'NOVA1', 'PTBP1')]}
    template, _ = configure.generate(request, paths, tmp_path / 'configuration')
    assert [a['asset_id'] for a in template['assets']] == ['GO', 'KEGG', 'NOVA1', 'PTBP1']
    assert len({a['metadata']['source_url'] for a in template['assets']}) == 4


def test_disabled_modules_do_not_download_catalog(tmp_path, monkeypatch):
    # Selection occurs before catalog binding in main.nf; empty selections have no network effects.
    monkeypatch.setattr(configure.urllib.request, 'urlopen', lambda *a, **k: (_ for _ in ()).throw(AssertionError('network')))
    configure.generate({**{'skip_' + k: True for k in configure.ANALYSES}, 'resource_catalog': []}, [], tmp_path / 'cfg')


def test_forced_verification_bypasses_both_digest_caches(tmp_path, monkeypatch):
    source = tmp_path / 'source'; source.write_bytes(b'reference')
    monkeypatch.setenv('KATANA_DIGEST_CACHE', str(tmp_path / 'cache'))
    expected = file_io.digest(source)
    opened = []; original = file_io.hashlib.file_digest
    def read(*args, **kwargs):
        opened.append(True)
        return original(*args, **kwargs)
    monkeypatch.setattr(file_io.hashlib, 'file_digest', read)
    monkeypatch.setenv('KATANA_VERIFY_CONTENT', '1')
    assert file_io.digest(source) == expected
    assert opened == [True]


def test_owned_native_import_survives_work_removal_and_borrowed_input_is_independent(tmp_path):
    import shutil
    from test_flow_import import bundle
    from katana_next.flow.importer import import_bundle
    from katana_next.flow.storage_project import _NativeCopies, _verified_copy
    work = tmp_path / 'work'; work.mkdir()
    root, _, template = bundle(work)
    source = root / 'SE.MATS.JCEC.txt'
    project = import_bundle(root, template, work / 'import', owned_link_root=work)
    imported = project.parent / 'bundle' / source.name
    assert os.path.samefile(source, imported)
    collection = work / 'collection'; collection.mkdir()
    owners = _NativeCopies(work, project, collection)
    retained = tmp_path / 'retained' / source.name
    _verified_copy(imported, retained, owned_copies=owners)
    assert os.path.samefile(imported, retained)
    borrowed = tmp_path / 'external'; borrowed.write_bytes(b'external unchanged')
    copied = tmp_path / 'retained' / 'external'
    _verified_copy(borrowed, copied, owned_copies=owners)
    assert not os.path.samefile(borrowed, copied)
    copied.write_bytes(b'new bytes')
    assert borrowed.read_bytes() == b'external unchanged'
    expected = retained.read_bytes()
    shutil.rmtree(work)
    assert retained.read_bytes() == expected


def test_digest_accepts_metadata_only_transition_after_identical_stable_recheck(tmp_path, monkeypatch):
    source = tmp_path / 'completed.tsv'; source.write_bytes(b'unchanged native result')
    initial = file_io.identity(source); phase = [0]; reads = []
    real = file_io.hashlib.file_digest
    monkeypatch.setattr(file_io, 'identity', lambda path: [*initial[:4], initial[4] + phase[0]])
    monkeypatch.setenv('KATANA_VERIFY_CONTENT', '1')
    def read(handle, algorithm):
        value = real(handle, algorithm); reads.append(True); phase[0] = 1
        return value
    monkeypatch.setattr(file_io.hashlib, 'file_digest', read)
    assert file_io.digest(source) == hashlib.sha256(source.read_bytes()).hexdigest()
    assert len(reads) == 2


def test_digest_rejects_real_same_size_content_change_during_verification(tmp_path, monkeypatch):
    source = tmp_path / 'changed.tsv'; source.write_bytes(b'original')
    initial = file_io.identity(source); phase = [0]; stat = source.stat()
    real = file_io.hashlib.file_digest
    monkeypatch.setattr(file_io, 'identity', lambda path: [*initial[:4], initial[4] + phase[0]])
    monkeypatch.setenv('KATANA_VERIFY_CONTENT', '1')
    def read(handle, algorithm):
        value = real(handle, algorithm)
        if phase[0] == 0:
            source.write_bytes(b'modified')
            os.utime(source, ns=(stat.st_atime_ns, stat.st_mtime_ns))
            phase[0] = 1
        return value
    monkeypatch.setattr(file_io.hashlib, 'file_digest', read)
    with pytest.raises(ValueError, match='File changed during content verification'):
        file_io.digest(source)


def test_digest_rejects_unsettled_identity_without_unbounded_retries(tmp_path, monkeypatch):
    source = tmp_path / 'unstable.tsv'; source.write_bytes(b'constant')
    initial = file_io.identity(source); reads = []
    def identity(path):
        reads.append(True); return [*initial[:4], initial[4] + len(reads)]
    monkeypatch.setattr(file_io, 'identity', identity)
    monkeypatch.setenv('KATANA_VERIFY_CONTENT', '1')
    with pytest.raises(ValueError, match='File changed during content verification'):
        file_io.digest(source)
    assert len(reads) == 3


def test_generated_projection_checksum_uses_written_bytes_and_preserves_existing_file(tmp_path, monkeypatch):
    target = tmp_path / 'projection.tsv'; content = b'gene\tevent\tpvalue\nG1\tE1\t1e-310\n' * 1000
    monkeypatch.setattr(file_io.hashlib, 'file_digest', lambda *args: (_ for _ in ()).throw(AssertionError('projection reread')))
    checksum = file_io.write_verified(target, content)
    assert target.read_bytes() == content
    assert checksum == hashlib.sha256(content).hexdigest()
    with pytest.raises(FileExistsError): file_io.write_verified(target, b'replacement')
    assert target.read_bytes() == content


def test_generated_projection_flush_failure_does_not_return_success(tmp_path, monkeypatch):
    def fail(fd): raise OSError('simulated flush failure')
    monkeypatch.setattr(file_io.os, 'fsync', fail)
    with pytest.raises(OSError, match='simulated flush failure'):
        file_io.write_verified(tmp_path / 'partial.tsv', b'result')
