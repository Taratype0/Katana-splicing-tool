"""Bind complete, already fitted native MAPP families into a portable report.

Original catalogues and numerical files are preserved. The additional catalogue
contains navigational paths relative to its own location; it does not introduce
a new correction family, fit a model or change an activity or P value.
"""
from __future__ import annotations

from copy import deepcopy
from hashlib import sha256
import json
import os
from pathlib import Path, PurePosixPath, PureWindowsPath
import re

from .sources import digest


def _bounded_relative(value, root, *, label):
    if not isinstance(value, str) or not value or '\\' in value:
        raise ValueError(label + " requires a portable relative path")
    relative = PurePosixPath(value)
    if relative.is_absolute() or PureWindowsPath(value).drive or '..' in relative.parts:
        raise ValueError(label + " path escapes its declared root")
    path = (root / value).resolve()
    if not path.is_relative_to(root) or not path.is_file():
        raise ValueError(label + " is missing or escapes its declared root")
    return path


def _sha_identity(value):
    if not isinstance(value, str) or not re.fullmatch(r'[0-9a-f]{64}', value):
        raise ValueError("Native MAPP fit requires a complete SHA256 identity")
    return value


def prepare_portable_mapp_panel(panel, *, report_root, copy_file, collection_root=None):
    """Return ``(panel_copy, audit)`` before Builder binds its data links.

    ``copy_file(path, prefix='results')`` must copy a source without changing it
    and return a report-relative destination. ``Builder.copy`` provides this
    contract. The explicit collection root is supplied by the verified producer
    on its native catalogue link, or by a caller adapting a legacy collection.
    Other panel types return an independent copy and ``None``.
    """
    result = deepcopy(panel)
    metadata = result.get('metadata', {})
    if (metadata.get('component') != 'mapp_native' and
            metadata.get('analysis_id') != 'prioritize.mapp_native'):
        return result, None
    root = Path(report_root).expanduser().resolve()
    if not root.is_dir():
        raise ValueError("Portable MAPP catalogue requires an existing new report directory")
    cid = result.get('contrast_id')
    if not isinstance(cid, str) or not cid or cid == 'all':
        raise ValueError("Native MAPP catalogue requires one exact comparison")
    native_sha = _sha_identity(metadata.get('native_fit_index_sha256'))
    links = result.get('data_links', [])
    native_links = [link for link in links if link.get('role') == 'full_results' and
                    (link.get('source_sha256') == native_sha or
                     (link.get('path') and Path(link['path']).is_file() and digest(Path(link['path'])) == native_sha))]
    if len(native_links) != 1:
        raise ValueError("Native MAPP panel needs one exact original complete-fit catalogue")
    native_link = native_links[0]
    if native_link.get('source_sha256') not in (None, native_sha):
        raise ValueError("Native MAPP catalogue link has a conflicting source identity")
    declared_root = collection_root if collection_root is not None else native_link.get('source_collection_root')
    if not declared_root or not Path(declared_root).is_absolute():
        raise ValueError("Native MAPP catalogue requires its explicit absolute source collection root")
    collection = Path(declared_root).resolve()
    if root == collection or root.is_relative_to(collection) or collection.is_relative_to(root):
        raise ValueError("Portable MAPP report must be separate from the original collection")
    catalog_path = Path(native_link['path']).resolve()
    if not catalog_path.is_relative_to(collection) or not catalog_path.is_file():
        raise ValueError("Native MAPP catalogue is missing or escapes its declared collection")
    catalog_bytes = catalog_path.read_bytes()
    if sha256(catalog_bytes).hexdigest() != native_sha:
        raise ValueError("Original native MAPP catalogue changed after completion")
    source_index = json.loads(catalog_bytes.decode('utf-8-sig'))
    rows = source_index.get('fits')
    if (source_index.get('comparison_id') != cid or source_index.get('calibrated_motiftrace_parallel') is not False
            or metadata.get('calibrated_motiftrace_parallel') is not False or not isinstance(rows, list) or not rows):
        raise ValueError("Native MAPP catalogue has a different comparison or statistical family")
    declared_ids = metadata.get('fit_task_ids')
    if (not isinstance(declared_ids, list) or not declared_ids
            or any(not isinstance(token, str) or not token for token in declared_ids)
            or len(declared_ids) != len(set(declared_ids))):
        raise ValueError("Native MAPP panel requires distinct completed task identities")
    ids, verified = [], []
    for row in rows:
        if not isinstance(row, dict) or not isinstance(row.get('scope'), dict):
            raise ValueError("Native MAPP catalogue requires complete fit scope records")
        token = row.get('task_id')
        if not isinstance(token, str) or not token or token in ids:
            raise ValueError("Native MAPP catalogue contains a missing or duplicate fit identity")
        ids.append(token)
        if row.get('scope', {}).get('contrast_ids') != [cid]:
            raise ValueError("Native MAPP catalogue fit belongs to another comparison")
        expected = _sha_identity(row.get('sha256'))
        path = _bounded_relative(row.get('path'), collection, label='Native MAPP fit')
        before = path.stat()
        if digest(path) != expected:
            raise ValueError("Native MAPP numerical fit bytes changed after completion")
        verified.append((row, path, expected, before.st_size, before.st_mtime_ns))
    if set(ids) != set(declared_ids) or len(rows) != len(declared_ids):
        raise ValueError("Native MAPP catalogue does not retain every declared completed fit")
    destination_dir = root / 'outputs' / 'figure-data'
    destination_dir.mkdir(parents=True, exist_ok=True)
    portable = deepcopy(source_index)
    portable.update(schema='katana.portable_native_MAPP_fit_index/v1', path_base='index_file_parent',
                    native_source_catalog_sha256=native_sha, statistics_or_scope_changed=False)
    copied = []
    for target_row, (original_row, source, expected, size, mtime) in zip(portable['fits'], verified, strict=True):
        relative = copy_file(source, prefix='results')
        target = _bounded_relative(relative, root, label='Portable MAPP fit')
        if digest(target) != expected or target.stat().st_size != size:
            raise ValueError("Portable MAPP copy does not match the completed native fit")
        target_row['native_source_collection_path'] = original_row['path']
        target_row['path'] = Path(os.path.relpath(target, destination_dir)).as_posix()
        copied.append({'task_id': original_row['task_id'], 'contrast_id': cid,
                       'source_collection_path': original_row['path'], 'report_path': relative,
                       'path_from_catalog': target_row['path'], 'sha256': expected, 'bytes': size})
    payload = (json.dumps(portable, ensure_ascii=False, indent=2, allow_nan=False) + '\n').encode('utf-8')
    portable_sha = sha256(payload).hexdigest()
    portable_path = destination_dir / (portable_sha[:16] + '-complete_native_fits_portable.json')
    if portable_path.exists():
        if portable_path.read_bytes() != payload:
            raise FileExistsError("A different portable MAPP catalogue already occupies this output")
    else:
        with portable_path.open('xb') as handle:
            handle.write(payload)
    # Keep the authenticated source catalogue as an unmodified provenance link.
    # Its native collection-relative paths must not masquerade as report URLs.
    old_portable = metadata.get('portable_fit_index_sha256')
    if old_portable == native_sha:
        raise ValueError("Native and portable MAPP catalogue identities must remain separate")
    if old_portable:
        links[:] = [link for link in links if link.get('source_sha256') != old_portable]
    native_link['label'] = 'Original native catalog (provenance; source-collection paths)'
    portable_link = deepcopy(native_link)
    portable_link.update(path=str(portable_path), source_sha256=portable_sha,
                         label='Every completed native fit — portable file index')
    portable_link.pop('source_collection_root', None)
    links.insert(1 if links else 0, portable_link)
    metadata.update(portable_fit_index_sha256=portable_sha, portable_fit_index_path_base='index_file_parent',
                    original_native_fit_index_sha256=native_sha)
    for _, path, expected, size, mtime in verified:
        if (path.stat().st_size, path.stat().st_mtime_ns, digest(path)) != (size, mtime, expected):
            raise RuntimeError("Completed native MAPP fit changed during portable collection")
    if catalog_path.read_bytes() != catalog_bytes:
        raise RuntimeError("Original native MAPP catalogue changed during portable collection")
    return result, {'schema': 'katana.portable_native_MAPP_fit_index_acceptance/v1', 'passed': True,
                    'contrast_id': cid, 'fit_count': len(copied), 'fits': copied,
                    'native_catalog_sha256': native_sha,
                    'portable_catalog': portable_path.relative_to(root).as_posix(),
                    'portable_catalog_sha256': portable_sha, 'native_statistics_unchanged': True,
                    'new_fits': 0, 'new_statistical_tests': 0}


def validate_portable_mapp_panel(panel, *, report_root):
    """Verify the final report's complete MAPP links without source-machine paths.

    Call after data-link binding (and again after compaction/transport). The
    original copied native catalogue is immutable evidence, while its separate
    portable projection resolves every exact completed fit in this report.
    """
    root = Path(report_root).resolve()
    metadata = panel.get('metadata', {})
    cid = panel.get('contrast_id')
    if (metadata.get('component') != 'mapp_native' or metadata.get('calibrated_motiftrace_parallel') is not False
            or not isinstance(cid, str) or not cid or cid == 'all'):
        raise ValueError("Portable native MAPP acceptance requires one exact completed comparison panel")
    native_sha = _sha_identity(metadata.get('native_fit_index_sha256'))
    portable_sha = _sha_identity(metadata.get('portable_fit_index_sha256'))
    if (metadata.get('original_native_fit_index_sha256') != native_sha or
            metadata.get('portable_fit_index_path_base') != 'index_file_parent' or native_sha == portable_sha):
        raise ValueError("Portable native MAPP catalogue must retain its separate original identity")
    indices = []
    for identity in (native_sha, portable_sha):
        links = [link for link in panel.get('data_links', []) if link.get('role') == 'full_results'
                 and ((identity == native_sha and link.get('source_sha256') == identity) or
                      (identity == portable_sha and
                       link.get('portable_sha256', link.get('source_sha256')) == identity))]
        if len(links) != 1 or links[0].get('contrast_id') != cid or links[0].get('policy') is not None:
            raise ValueError("Portable native MAPP panel lacks its unique comparison-scoped full-fit catalogue")
        path = _bounded_relative(links[0].get('path'), root, label='Portable MAPP catalogue')
        expected = links[0].get('portable_sha256', identity)
        if digest(path) != expected:
            raise ValueError("Portable native MAPP catalogue bytes changed after report binding")
        # A derived portable index may be rewritten by the lossless storage
        # pass, whose separately bound current SHA is updated in metadata. The
        # native source catalogue must always retain its original exact bytes.
        if digest(path) != identity:
            raise ValueError("Portable or original native MAPP catalogue SHA does not match its bound identity")
        indices.append((path, json.loads(path.read_text(encoding='utf-8-sig'))))
    (native_path, native), (portable_path, portable) = indices
    if (native.get('comparison_id') != cid or portable.get('comparison_id') != cid
            or native.get('calibrated_motiftrace_parallel') is not False
            or portable.get('calibrated_motiftrace_parallel') is not False
            or portable.get('schema') != 'katana.portable_native_MAPP_fit_index/v1'
            or portable.get('path_base') != 'index_file_parent'
            or portable.get('native_source_catalog_sha256') != native_sha
            or portable.get('statistics_or_scope_changed') is not False):
        raise ValueError("Portable native MAPP index changed its comparison or correction-family identity")
    if any(portable.get(key) != value for key, value in native.items() if key not in {'schema', 'fits'}):
        raise ValueError("Portable native MAPP index changed an original correction or catalogue record")
    original_rows, portable_rows = native.get('fits'), portable.get('fits')
    expected_ids = metadata.get('fit_task_ids')
    if (not isinstance(original_rows, list) or not original_rows or not isinstance(portable_rows, list)
            or len(original_rows) != len(portable_rows) or not isinstance(expected_ids, list)
            or any(not isinstance(token, str) or not token for token in expected_ids)
            or len(expected_ids) != len(set(expected_ids))):
        raise ValueError("Portable native MAPP index is missing part of the complete fit family")
    seen, records = set(), []
    for original, row in zip(original_rows, portable_rows, strict=True):
        if not isinstance(original, dict) or not isinstance(row, dict):
            raise ValueError("Portable native MAPP index requires complete fit records")
        reconstructed = deepcopy(row)
        source_path = reconstructed.pop('native_source_collection_path', None)
        reconstructed['path'] = source_path
        token = original.get('task_id')
        if (reconstructed != original or not isinstance(token, str) or not token or token in seen
                or not isinstance(original.get('scope'), dict) or original['scope'].get('contrast_ids') != [cid]):
            raise ValueError("Portable native MAPP fit changed its native task, scope or statistical record")
        seen.add(token)
        relative = row.get('path')
        if not isinstance(relative, str) or not relative or '\\' in relative:
            raise ValueError("Portable native MAPP fit requires an index-relative URL")
        parsed = PurePosixPath(relative)
        if parsed.is_absolute() or PureWindowsPath(relative).drive:
            raise ValueError("Portable native MAPP fit URL escapes the report")
        target = (portable_path.parent / relative).resolve()
        if not target.is_relative_to(root) or not target.is_file():
            raise ValueError("Portable native MAPP fit URL is missing or escapes the report")
        if digest(target) != _sha_identity(original.get('sha256')):
            raise ValueError("Portable native MAPP fit bytes differ from the unchanged completed source")
        records.append({'task_id': token, 'report_path': target.relative_to(root).as_posix(),
                        'sha256': original['sha256']})
    if seen != set(expected_ids):
        raise ValueError("Portable native MAPP index does not retain every declared completed fit")
    return {'schema': 'katana.portable_native_MAPP_final_report_acceptance/v1', 'passed': True,
            'contrast_id': cid, 'fit_count': len(records), 'fits': records,
            'native_catalog_sha256': native_sha, 'portable_catalog_sha256': portable_sha,
            'native_statistics_unchanged': True, 'source_machine_paths_required': False,
            'new_fits': 0, 'new_statistical_tests': 0}
