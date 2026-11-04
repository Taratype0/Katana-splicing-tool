"""Native CSV provenance is distinct from lossless TSV input identity."""
import csv
import hashlib
import io
import json
import shutil
from pathlib import Path

import pytest

from katana_next.flow.model import file_sha256
from katana_next.flow.projections import dexseq_projection_identity
from katana_next.flow.storage_project import retain_collection_project, validate_retained_project
from katana_next.infrastructure.manifest_store import ManifestStore
from test_flow_import import bundle
from test_flow_storage_project import _collect_fixture


def add_projection(payload, work, *, legacy=True):
    from katana_next.flow.importer import import_bundle
    from katana_next.flow.model import write_json
    upstream, manifest, template = bundle(work)
    native = upstream / 'dexseq_exon/DEXSeqResults.B_vs_A.csv'
    native.parent.mkdir()
    native.write_bytes(b'groupID,featureID,pvalue,padj,log2fold_A_B\r\ng1,E001,1e-310,2e-310,-0.7\r\ng1,E002,NA,NA,0.5\r\n')
    manifest['artifacts'] = [{'kind': 'auxiliary', 'relative_path': 'dexseq_exon/' + native.name,
        'comparison_id': 'B_vs_A', 'status': 'available', 'bytes': native.stat().st_size,
        'sha256': file_sha256(native)}]
    write_json(upstream / 'manifest.json', manifest)
    imported = import_bundle(upstream, template, work / 'imported')
    _, model = ManifestStore().load(imported)
    asset = model.assets[0].model_dump(mode='json')
    # Register actual imported files under the fixture's existing input root.
    inputs = Path(payload['input_root'])
    native_copy = inputs / 'native.csv'
    native_copy.write_bytes(native.read_bytes())
    tsv = inputs / 'native.dexseq-exons.tsv'
    tsv.write_bytes(model.assets[0].path.read_bytes())
    asset['path'] = str(tsv)
    asset['metadata']['native_table_projection']['native_path'] = native_copy.name
    if legacy:
        asset['metadata'].pop('sha256')
    payload['assets'].append(asset)


@pytest.mark.parametrize('legacy', [True, False])
def test_import_collect_retain_relocate_preserves_csv_and_tsv(tmp_path, legacy):
    work, project, collection = _collect_fixture(tmp_path / 'work',
        configure=lambda p,w: add_projection(p,w,legacy=legacy))
    original_project = project.read_bytes()
    before = {str(p): p.read_bytes() for p in work.rglob('*') if p.name in {'native.csv','native.dexseq-exons.tsv'}}
    retained = tmp_path / 'retained'
    retain_collection_project(collection, project, retained)
    assert project.read_bytes() == original_project
    assert all(Path(p).read_bytes() == b for p,b in before.items())
    moved = tmp_path / 'moved'
    retained.rename(moved)
    assert work.resolve().is_relative_to(tmp_path.resolve())
    shutil.rmtree(work)
    context = validate_retained_project(moved)
    asset = next(a for a in context.manifest.assets if a.path.name == 'native.dexseq-exons.tsv')
    native, native_sha, projected_sha = dexseq_projection_identity(asset, context.manifest.input_root)
    assert native.is_relative_to(moved) and native.read_bytes() in before.values()
    assert native_sha == asset.metadata['upstream_sha256']
    assert projected_sha == file_sha256(asset.path) == asset.metadata['sha256']
    assert list(csv.reader(io.StringIO(native.read_text()))) == list(csv.reader(io.StringIO(asset.path.read_text()), delimiter='\t'))


@pytest.mark.parametrize('corruption', ['native', 'projected', 'declared', 'escape'])
def test_projection_identity_does_not_allow_changed_bytes_or_escape(tmp_path, corruption):
    def configure(payload, work):
        add_projection(payload,work,legacy=False)
        asset = payload['assets'][-1]
        if corruption == 'declared': asset['metadata']['sha256'] = '0' * 64
        if corruption == 'escape': asset['metadata']['native_table_projection']['native_path'] = '../outside.csv'
    work, project, collection = _collect_fixture(tmp_path / 'work', configure=configure)
    _, source = ManifestStore().load(project)
    if corruption == 'native': (source.input_root / 'native.csv').write_bytes(b'changed')
    if corruption == 'projected': (source.input_root / 'native.dexseq-exons.tsv').write_bytes(b'changed')
    with pytest.raises(ValueError, match='SHA256|changed|projection source path'):
        retain_collection_project(collection, project, tmp_path / 'rejected')
    assert not (tmp_path / 'rejected').exists()
