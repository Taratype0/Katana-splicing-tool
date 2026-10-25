"""Restore verified target products into a new local execution workspace.

This imports retained completions, not a new analysis or a replacement for
their original provenance. The target renderer and method identities are not
changed by the local ownership and request-bundle records created here.
"""
from __future__ import annotations

from contextlib import contextmanager
import json
import os
from pathlib import Path
import shutil
import stat
import tempfile
import uuid
import zipfile

from ._execution import bounded, fingerprint, sha256, utc_now, write_json
from .bundle import _json, _relative, open_bundle
from .runner import RESULT_SCHEMA, inspect_results


RESTORE_SCHEMA = 'katana.target_restore/v1'
REQUEST_ARCHIVE = 'TARGET_REQUEST.zip'
RESTORE_RECORD = '.target-restore.json'
_LOCAL_NAMES = {'.target-owner.json', RESTORE_RECORD, REQUEST_ARCHIVE}


def _contract(document):
    if (not isinstance(document, dict) or document.get('schema') != RESULT_SCHEMA
            or not isinstance(document.get('members'), dict)
            or not isinstance(document.get('targets'), list)
            or any(not isinstance(row, dict) for row in document['targets'])
            or any(not isinstance(proof, dict) for proof in document['members'].values())):
        raise ValueError('Unsupported target results schema or product ledger')
    for name in document['members']:
        _relative(name)
    for field, default in (('archive', 'REPORT.zip'), ('report', 'INDEX.html')):
        _relative(document.get(field, default))
    return document


def _snapshot(path):
    path = Path(path)
    info = path.stat()
    return {'sha256': sha256(path), 'bytes': info.st_size, 'mtime_ns': info.st_mtime_ns}


def _claims_absent(root):
    claims = root / '.claims'
    if claims.is_symlink() or (claims.exists() and (not claims.is_dir() or any(claims.rglob('*')))):
        raise ValueError('Source has active or unknown target writer claims; do not restore it while writers may be present')


@contextmanager
def _source(path):
    if path.is_file() and path.suffix.lower() == '.zip':
        with tempfile.TemporaryDirectory(prefix='katana-target-restore-') as directory:
            root = Path(directory)
            with zipfile.ZipFile(path) as packed:
                infos = packed.infolist()
                names = [info.filename for info in infos]
                if len(names) != len(set(names)) or len(names) != len({name.casefold() for name in names}):
                    raise ValueError('Duplicate or case-aliased target ZIP member')
                for info in infos:
                    _relative(info.filename)
                    if info.is_dir() or stat.S_ISLNK(info.external_attr >> 16):
                        raise ValueError('Target ZIP must contain regular ledger files, not directories or symlinks')
                if 'result-manifest.json' not in names:
                    raise ValueError('A completed target result ZIP with result-manifest.json is required')
                with packed.open('result-manifest.json') as stream:
                    document = json.load(stream)
                _contract(document)
                if set(names) != set(document['members']) | {'result-manifest.json'}:
                    raise ValueError('Target ZIP member set differs from its complete product ledger')
                for info in infos:
                    if info.filename != 'result-manifest.json':
                        expected = document['members'][info.filename]
                        if type(expected.get('bytes')) is not int or info.file_size != expected['bytes']:
                            raise ValueError('Target ZIP member size differs from its ledger')
                    destination = bounded(root, info.filename)
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    with packed.open(info) as source, destination.open('xb') as target:
                        shutil.copyfileobj(source, target)
            yield root, document, 'zip'
    else:
        root = path.parent if path.is_file() else path
        if path.is_file() and path.name != 'result-manifest.json':
            raise ValueError('Expected a target ZIP, result directory or result-manifest.json')
        _claims_absent(root)
        yield root, _contract(_json(root / 'result-manifest.json')), 'directory'


def _zip_members(root, names, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(destination, 'x', zipfile.ZIP_DEFLATED, compresslevel=6,
                         strict_timestamps=False) as packed:
        for name in sorted(names):
            packed.write(bounded(root, name), name)


def _new_local_owner(output, request):
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise FileExistsError('Restore output is no longer empty; preserve the other writer or existing files')
    expected = {'schema': 'katana.target_output_owner/v1', 'request_id': request['request_id'],
                'project_id': request['project']['project_id'], 'request_sha256': fingerprint(request)}
    # Unlike a normal run, restore never adopts even a same-request existing
    # owner. O_EXCL also keeps simultaneous restores from sharing a writer.
    descriptor = os.open(output / '.target-owner.json', os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    with os.fdopen(descriptor, 'w', encoding='utf-8') as stream:
        json.dump(expected, stream, sort_keys=True)


def _request(root, document, destination):
    """Authenticate the exact copied request, not the origin's project DB."""
    request = _json(bounded(root, 'request.json'))
    if not isinstance(request, dict) or not isinstance(request.get('members'), dict):
        raise ValueError('Result package does not contain a portable target request ledger')
    names = {'request.json', *request.get('members', {})}
    if not names <= set(document['members']):
        raise ValueError('Result package does not contain its full portable target request')
    _zip_members(root, names, destination)
    with open_bundle(destination) as bundle:
        request = bundle.request
    project = request['project']
    if (document.get('request_id') != request['request_id']
            or document.get('project_id') != project['project_id']
            or document.get('source_project_manifest_sha256') != project['source_manifest_sha256']):
        raise ValueError('Result and portable request have different project, source manifest or request identities')
    declared = {target['target_id']: target for target in request['targets']}
    rows = document.get('targets', [])
    if len(rows) != len(declared) or {row['target_id'] for row in rows} != set(declared):
        raise ValueError('Result targets differ from the complete portable request')
    for row in rows:
        wanted = declared[row['target_id']]
        candidates = [row]
        if row.get('manifest'):
            candidates.append(_json(bounded(root, row['manifest'])))
        for actual in candidates:
            if any(actual.get(key) != wanted.get(key) for key in ('target_kind', 'gene_id', 'gene_symbol', 'event_id')):
                raise ValueError('Returned target biological identity differs from its portable request')
            if 'comparison_ids' in actual and actual['comparison_ids'] != wanted['comparison_ids']:
                raise ValueError('Returned target comparison scope differs from its portable request')
    owner = root / '.target-owner.json'
    if owner.exists():
        expected = {'schema': 'katana.target_output_owner/v1', 'request_id': request['request_id'],
                    'project_id': project['project_id'], 'request_sha256': fingerprint(request)}
        if owner.is_symlink() or _json(owner) != expected:
            raise ValueError('Source local owner does not match the authenticated target request')
    return request


def restore_results(source, output_dir, *, expected_project_id=None):
    """Verify and copy existing products to a new empty, resumable directory.

    Resources remain external and are rebound for the next explicit plan/run.
    This operation neither invokes analysis workers nor changes their source
    fingerprints, completion statuses, native results or original provenance.
    """
    source = Path(source).expanduser().resolve()
    output = Path(output_dir).expanduser().resolve()
    source_root = source if source.is_dir() else source.parent
    if source.is_dir() or source.name == 'result-manifest.json':
        if output == source_root or output.is_relative_to(source_root) or source_root.is_relative_to(output):
            raise ValueError('Restore into a separate new directory, not inside or around the source')
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise FileExistsError('Restore requires a new empty directory; unknown or completed outputs are never overwritten')
    if not source.exists():
        raise FileNotFoundError(source)
    source_marker = source if source.is_file() else source / 'result-manifest.json'
    source_before = _snapshot(source_marker)
    with _source(source) as (root, document, kind), tempfile.TemporaryDirectory(prefix='katana-target-request-') as scratch:
        names = set(document.get('members', {}))
        if any(name.split('/')[0].casefold() in {value.casefold() for value in _LOCAL_NAMES | {'.claims'}}
               for name in names):
            raise ValueError('Product ledger cannot replace local ownership, request archive or writer claims')
        if len(names | {'result-manifest.json'}) != len({name.casefold() for name in names | {'result-manifest.json'}}):
            raise ValueError('Case-aliased product paths are not cross-platform portable')
        verified = inspect_results(root)
        if not verified['verified']:
            raise ValueError('Target products are not verified')
        archive_name = document.get('archive', 'REPORT.zip')
        if (archive_name.split('/')[0].casefold() in {value.casefold() for value in _LOCAL_NAMES | {'.claims'}}
                or archive_name.casefold() in {name.casefold() for name in names | {'result-manifest.json', *_LOCAL_NAMES}}):
            raise ValueError('Result archive path collides with a retained product or local control file')
        bounded(root, archive_name)
        request_archive = Path(scratch) / REQUEST_ARCHIVE
        request = _request(root, document, request_archive)
        if expected_project_id is not None and request['project']['project_id'] != expected_project_id:
            raise ValueError('Target package belongs to a different expected project')
        source_manifest_sha = sha256(root / 'result-manifest.json')
        source_times = {name: (bounded(root, name).stat().st_size, bounded(root, name).stat().st_mtime_ns)
                        for name in names | {'result-manifest.json'}}
        source_archive_before = _snapshot(verified['archive']) if verified.get('archive') else None
        _new_local_owner(output, request)
        record = {'schema': RESTORE_SCHEMA, 'restore_id': uuid.uuid4().hex, 'created_at': utc_now(),
                  'status': 'copying_verified_products', 'source_kind': kind,
                  'source_result_manifest_sha256': source_manifest_sha,
                  'source_request_sha256': sha256(root / 'request.json'),
                  'source_archive_sha256': source_before['sha256'] if kind == 'zip' else None,
                  'source_completion_status': document['status'], 'request_id': request['request_id'],
                  'project_id': document['project_id'], 'original_completion_history_preserved': True,
                  'new_local_owner_created': True, 'statistical_refit': False,
                  'native_BAM_counts_run': False, 'network_fit': False, 'motif_rescan': False}
        try:
            for name in sorted(names | {'result-manifest.json'}):
                target = bounded(output, name)
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(bounded(root, name), target)
            shutil.copy2(request_archive, output / REQUEST_ARCHIVE)
            target_archive = bounded(output, archive_name)
            target_archive.parent.mkdir(parents=True, exist_ok=True)
            if kind == 'zip':
                shutil.copy2(source, target_archive)
                record['archive_operation'] = 'original_archive_copied'
            elif verified.get('archive'):
                shutil.copy2(verified['archive'], target_archive)
                record['source_archive_sha256'] = sha256(Path(verified['archive']))
                record['archive_operation'] = 'original_archive_copied'
            else:
                _zip_members(output, names | {'result-manifest.json'}, target_archive)
                record['archive_operation'] = 'new_archive_of_unchanged_products'
            result = inspect_results(output)
            if sha256(output / 'result-manifest.json') != source_manifest_sha or _snapshot(source_marker) != source_before:
                raise ValueError('Source result changed during restore')
            if kind == 'directory':
                _claims_absent(root)
                inspect_results(root)
                if any((bounded(root, name).stat().st_size, bounded(root, name).stat().st_mtime_ns) != identity
                       for name, identity in source_times.items()):
                    raise ValueError('Source product metadata changed during restore')
                if source_archive_before and _snapshot(verified['archive']) != source_archive_before:
                    raise ValueError('Source archive changed during restore')
            record.update(status='restored_verified_products', completed_at=utc_now(),
                          source_control_unchanged=True, request_archive=REQUEST_ARCHIVE,
                          request_archive_sha256=sha256(output / REQUEST_ARCHIVE))
            write_json(output / RESTORE_RECORD, record)
            return {**result, 'restored': True, 'request': str(output / REQUEST_ARCHIVE),
                    'restore_provenance': str(output / RESTORE_RECORD), 'restoration': record}
        except Exception as exc:
            record.update(status='restore_failed', error=f'{type(exc).__name__}: {exc}')
            write_json(output / RESTORE_RECORD, record)
            raise
