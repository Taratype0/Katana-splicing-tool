"""Repeated portable refreshes preserve real same-filesystem source inodes."""
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil

import pytest

from katana_next.reports.engine import ASSET_ROOT, Builder
from katana_next.reports import refresh


def _completed_report(root):
    root.mkdir()
    shutil.copytree(ASSET_ROOT, root / 'assets')
    native = root / 'outputs' / 'native.tsv'
    native.parent.mkdir()
    native.write_text('event_id\tfdr\ne1\t0.25\n', encoding='utf-8')
    project = {'project_id': 'preserved-study', 'samples': [],
               'contrasts': [{'contrast_id': 'B_vs_A', 'numerator': {'label': 'B'},
                              'denominator': {'label': 'A'}}]}
    builder = Builder({'project': project, 'artifacts': [], 'mode': 'bulk'}, root, 'stable_significant')
    builder.tables.append({'name': 'Retained native results', 'source': 'outputs/native.tsv', 'chapter': 3})
    builder.create_archive = False
    builder.write()
    return root


def _snapshot(root):
    return {path.relative_to(root).as_posix():
            (sha256(path.read_bytes()).hexdigest(), path.stat().st_size, path.stat().st_mtime_ns)
            for path in root.rglob('*') if path.is_file()}


@pytest.mark.parametrize('compact', [False, True])
def test_repeated_refresh_atomically_replaces_real_hardlinked_lifecycle_receipt(tmp_path, monkeypatch, compact):
    original = _completed_report(tmp_path / 'original')
    source = tmp_path / 'first-refresh'
    refresh.refresh_report(original, source, create_archive=False, compact=compact)
    source_receipt = source / 'presentation-refresh.json'
    source_bytes = source_receipt.read_bytes()
    before = _snapshot(source)
    destination = tmp_path / 'second-refresh'
    links_observed = []
    native_link = os.link

    def observe_real_link(src, dst, *args, **kwargs):
        native_link(src, dst, *args, **kwargs)
        if Path(src) == source_receipt:
            links_observed.append(os.path.samefile(src, dst))

    monkeypatch.setattr(refresh.os, 'link', observe_real_link)
    refresh.refresh_report(source, destination, create_archive=True, compact=compact)

    assert links_observed == [True], 'Exercise a real shared inode, not a copy fallback.'
    assert _snapshot(source) == before
    assert source_receipt.read_bytes() == source_bytes
    target_receipt = destination / 'presentation-refresh.json'
    assert not os.path.samefile(source_receipt, target_receipt)
    receipt = json.loads(target_receipt.read_text(encoding='utf-8'))
    assert receipt['original_report_unchanged'] is True
    assert receipt['completed_analyses_rerun'] == 0
    assert (destination / 'REPORT.zip').is_file()
    assert (destination / 'outputs/native.tsv').read_bytes() == (source / 'outputs/native.tsv').read_bytes()


def test_source_receipt_guard_runs_after_the_final_archive_write(tmp_path, monkeypatch):
    original = _completed_report(tmp_path / 'original')
    source = tmp_path / 'first-refresh'
    refresh.refresh_report(original, source, create_archive=False)
    source_receipt = source / 'presentation-refresh.json'
    write_archive = refresh.write_archive

    def archive_then_external_mutation(output, archive_path=None):
        result = write_archive(output, archive_path)
        source_receipt.write_bytes(source_receipt.read_bytes() + b'\nexternal change\n')
        return result

    monkeypatch.setattr(refresh, 'write_archive', archive_then_external_mutation)
    with pytest.raises(RuntimeError, match='Original report changed'):
        refresh.refresh_report(source, tmp_path / 'second-refresh')
