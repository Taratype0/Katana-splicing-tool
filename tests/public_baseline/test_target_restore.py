"""Returned target archives restore without a project DB or implicit refits."""
import json
from pathlib import Path
import shutil
import stat
import zipfile

import pytest

from katana_next.targets import restore_results
from katana_next.targets._execution import sha256
from katana_next.targets.runner import inspect_results, plan_targets, run_targets
from test_target_runner import _fixture


@pytest.fixture(scope='module')
def complete_products(tmp_path_factory):
    root = tmp_path_factory.mktemp('completed-portable-targets')
    request, bindings, _ = _fixture(root)
    output = root / 'completed'
    result = run_targets(request, output, bindings, workers=2)
    assert result['status'] == 'complete'
    return output


def _snapshot(root):
    return {p.relative_to(root).as_posix(): (sha256(p), p.stat().st_size, p.stat().st_mtime_ns)
            for p in root.rglob('*') if p.is_file()}


def _source_copy(source, destination, *, archive=True):
    shutil.copytree(source, destination, copy_function=shutil.copy2)
    if not archive:
        (destination / 'REPORT.zip').unlink()
    return destination


def test_returned_zip_restores_complete_request_owner_and_cache_without_old_paths(tmp_path, complete_products):
    archive = tmp_path / 'returned-targets.zip'
    shutil.copy2(complete_products / 'REPORT.zip', archive)
    original_sha = sha256(archive)
    output = tmp_path / 'different machine with spaces' / 'restored'
    result = restore_results(archive, output, expected_project_id='portable-toy')
    assert result['restored'] and result['verified'] and result['status'] == 'complete'
    assert result['restoration']['source_completion_status'] == 'complete'
    assert result['restoration']['status'] == 'restored_verified_products'
    assert result['restoration']['source_result_manifest_sha256'] == sha256(complete_products / 'result-manifest.json')
    assert (output / 'result-manifest.json').read_bytes() == (complete_products / 'result-manifest.json').read_bytes()
    assert sha256(output / 'REPORT.zip') == original_sha
    request = Path(result['request'])
    assert request.name == 'TARGET_REQUEST.zip'
    planned = plan_targets(request, outdir=output)
    assert planned['ready'] and all(task['reusable'] for task in planned['tasks'])
    before = _snapshot(output)
    cached = run_targets(request, output, workers=2)
    assert cached['status'] == 'cache_hit' and cached['cache']['hits'] == 2
    assert _snapshot(output) == before and sha256(archive) == original_sha
    assert all(not result['restoration'][flag] for flag in ('statistical_refit', 'native_BAM_counts_run', 'network_fit', 'motif_rescan'))


@pytest.mark.parametrize('as_manifest', [False, True])
def test_unpacked_products_restore_with_a_new_verified_archive(tmp_path, complete_products, as_manifest):
    source = _source_copy(complete_products, tmp_path / 'unpacked-source', archive=False)
    before = _snapshot(source)
    result = restore_results(source / 'result-manifest.json' if as_manifest else source, tmp_path / 'restored')
    assert result['archive_status'] == 'verified'
    assert result['restoration']['archive_operation'] == 'new_archive_of_unchanged_products'
    assert inspect_results(Path(result['manifest']))['status'] == 'complete'
    assert _snapshot(source) == before


def test_foreign_expected_project_cannot_create_a_local_owner(tmp_path, complete_products):
    output = tmp_path / 'rejected'
    with pytest.raises(ValueError, match='different expected project'):
        restore_results(complete_products / 'REPORT.zip', output, expected_project_id='another-study')
    assert not output.exists()


@pytest.mark.parametrize('changed', ['source_manifest', 'target_gene'])
def test_result_request_biological_identity_is_checked_before_new_owner(tmp_path, complete_products, changed):
    source = _source_copy(complete_products, tmp_path / 'mixed', archive=False)
    path = source / 'result-manifest.json'
    document = json.loads(path.read_text('utf-8'))
    if changed == 'source_manifest':
        document['source_project_manifest_sha256'] = 'd' * 64
    else:
        document['targets'][0]['gene_id'] = 'DIFFERENT-GENE'
    path.write_text(json.dumps(document), encoding='utf-8')
    with pytest.raises(ValueError, match='identit'):
        restore_results(source, tmp_path / 'rejected')
    assert not (tmp_path / 'rejected').exists()


@pytest.mark.parametrize('name,symlink', [('../escape.txt', False), ('index.HTML', False), ('link', True)])
def test_unsafe_or_unledgered_zip_is_rejected_before_output(tmp_path, complete_products, name, symlink):
    archive = tmp_path / 'unsafe.zip'
    shutil.copy2(complete_products / 'REPORT.zip', archive)
    with zipfile.ZipFile(archive, 'a') as packed:
        if symlink:
            entry = zipfile.ZipInfo(name)
            entry.create_system = 3
            entry.external_attr = (stat.S_IFLNK | 0o777) << 16
            packed.writestr(entry, 'elsewhere')
        else:
            packed.writestr(name, 'unlisted')
    with pytest.raises(ValueError):
        restore_results(archive, tmp_path / 'rejected')
    assert not (tmp_path / 'rejected').exists() and not (tmp_path / 'escape.txt').exists()


def test_active_claim_and_foreign_local_owner_are_never_adopted(tmp_path, complete_products):
    source = _source_copy(complete_products, tmp_path / 'source')
    claim = source / '.claims/unknown.json'
    claim.parent.mkdir(exist_ok=True)
    claim.write_text('{"pid":"unknown"}', encoding='utf-8')
    with pytest.raises(ValueError, match='writer claims'):
        restore_results(source, tmp_path / 'claimed-output')
    assert claim.exists() and not (tmp_path / 'claimed-output').exists()
    claim.unlink()
    owner = source / '.target-owner.json'
    owner.write_text('{"project_id":"foreign"}', encoding='utf-8')
    with pytest.raises(ValueError, match='owner does not match'):
        restore_results(source, tmp_path / 'foreign-output')
    assert not (tmp_path / 'foreign-output').exists()


def test_existing_outputs_and_nested_source_are_preserved(tmp_path, complete_products):
    output = tmp_path / 'existing'
    output.mkdir()
    sentinel = output / 'keep.txt'
    sentinel.write_text('keep', encoding='utf-8')
    with pytest.raises(FileExistsError):
        restore_results(complete_products, output)
    with pytest.raises(ValueError, match='separate new directory'):
        restore_results(complete_products, complete_products / 'nested')
    assert sentinel.read_text('utf-8') == 'keep'


def test_source_change_during_copy_never_becomes_a_completed_restore(tmp_path, complete_products, monkeypatch):
    from katana_next.targets import restore
    source = _source_copy(complete_products, tmp_path / 'source', archive=False)
    output = tmp_path / 'restore-failed'
    create = restore._new_local_owner
    def create_then_external_change(root, request):
        create(root, request)
        manifest = source / 'result-manifest.json'
        manifest.write_bytes(manifest.read_bytes() + b'\n')
    monkeypatch.setattr(restore, '_new_local_owner', create_then_external_change)
    with pytest.raises(ValueError, match='changed during restore'):
        restore_results(source, output)
    record = json.loads((output / '.target-restore.json').read_text('utf-8'))
    assert record['status'] == 'restore_failed' and record['source_completion_status'] == 'complete'
    assert not record.get('source_control_unchanged', False)


@pytest.mark.parametrize('malformed', ['members', 'proof', 'archive', 'reserved_archive'])
def test_invalid_result_contract_is_rejected_before_local_ownership(tmp_path, complete_products, malformed):
    source = _source_copy(complete_products, tmp_path / 'invalid', archive=False)
    manifest = source / 'result-manifest.json'
    document = json.loads(manifest.read_text('utf-8'))
    if malformed == 'members':
        document['members'] = []
    elif malformed == 'proof':
        document['members'][next(iter(document['members']))] = None
    elif malformed == 'archive':
        document['archive'] = ['REPORT.zip']
    else:
        document['archive'] = '.claims/archive.zip'
    manifest.write_text(json.dumps(document), encoding='utf-8')
    before = _snapshot(source)
    with pytest.raises(ValueError):
        restore_results(source, tmp_path / 'rejected')
    assert not (tmp_path / 'rejected').exists() and _snapshot(source) == before


def test_partial_completion_remains_partial_after_local_restore(tmp_path):
    request, bindings, _ = _fixture(tmp_path)
    source = tmp_path / 'partial-products'
    partial = run_targets(request, source, bindings, target_ids=['gene-G1'], workers=1)
    assert partial['status'] == 'partial'
    result = restore_results(source / 'REPORT.zip', tmp_path / 'restored')
    assert result['status'] == 'partial'
    assert result['restoration']['source_completion_status'] == 'partial'
    assert {target['status'] for target in result['targets']} == {'complete', 'not_requested'}
    plan = plan_targets(result['request'], outdir=tmp_path / 'restored')
    assert {target['target_id']: target['reusable'] for target in plan['tasks']} == {'gene-G1': True, 'event-G2': False}


def test_simultaneous_restores_cannot_share_or_overwrite_the_same_local_owner(tmp_path, complete_products):
    from concurrent.futures import ThreadPoolExecutor
    output = tmp_path / 'one-owner'
    with ThreadPoolExecutor(max_workers=2) as pool:
        pending = [pool.submit(restore_results, complete_products / 'REPORT.zip', output) for _ in range(2)]
        successes, failures = [], []
        for future in pending:
            try:
                successes.append(future.result())
            except FileExistsError as error:
                failures.append(error)
    assert len(successes) == len(failures) == 1
    assert inspect_results(output)['status'] == 'complete'
    record = json.loads((output / '.target-restore.json').read_text('utf-8'))
    assert record['status'] == 'restored_verified_products'
