"""Legacy completed collections receive the approved display on a new copy."""
from hashlib import sha256
import json
from pathlib import Path

import pytest

from katana_next.reports import network_display_upgrade as upgrade
from katana_next.reports.rich_results import copy_rich_asset


def _browser(path, *, payload=None, current=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    text = '<script id="data" type="application/json">' + json.dumps(payload or {'members': ['e1']}) + '</script>'
    if current:
        text += '<script id="module-display-data" type="application/json">{"modules":{"M01":{"edges":[[0,1,0.5]]}}}</script><div id="gene-body"></div><div id="network-display-size"></div>'
    path.write_text(text, encoding='utf-8')
    return path


def _asset(source):
    return {'id': 'registered-browser', 'path': str(source), 'status': 'available',
            'metadata': {'portable_root': str(source.parent.parent), 'relative_path': 'browser/ALL_MODULES.html',
                         'network_view_role': 'module_browser', 'local_dependencies': []}}


def test_generic_rich_copy_invokes_upgrade_on_independent_destination(tmp_path, monkeypatch):
    source = _browser(tmp_path / 'source/browser/ALL_MODULES.html')
    before = source.read_bytes()
    def refresh(asset, target):
        assert Path(asset['path']) == source and target != source
        assert target.read_bytes() == before
        target.write_text('approved complete genes and graph', encoding='utf-8')
        return {'passed': True, 'network_discovery_runs': 0}
    monkeypatch.setattr(upgrade, 'upgrade_copied_network_browser', refresh)
    result = copy_rich_asset(_asset(source), tmp_path / 'new-report')
    assert source.read_bytes() == before
    assert result['network_presentation_refresh']['passed']
    row = result['files'][0]
    assert row['source_sha256'] == sha256(before).hexdigest()
    assert row['portable_sha256'] == sha256(b'approved complete genes and graph').hexdigest()
    assert row['original_preserved'] is True


def test_existing_current_browser_does_not_repeat_drawing_export(tmp_path, monkeypatch):
    source = _browser(tmp_path / 'source/browser/ALL_MODULES.html', current=True)
    monkeypatch.setattr(upgrade, '_completed_inputs', lambda *_: pytest.fail('Do not repeat drawing'))
    assert upgrade.upgrade_copied_network_browser(_asset(source), tmp_path / 'copy/browser/ALL_MODULES.html') is None


def test_null_display_payload_is_not_mistaken_for_completed_graph(tmp_path, monkeypatch):
    source = _browser(tmp_path / 'source/browser/ALL_MODULES.html')
    source.write_text(source.read_text() + '<script id="module-display-data" type="application/json">null</script><div id="gene-body"></div><div id="network-display-size"></div>')
    asset = _asset(source)
    asset['metadata']['completion_scope'] = 'flow_verified_network_complete'
    monkeypatch.setattr(upgrade, '_completed_inputs', lambda *_: (_ for _ in ()).throw(ValueError('authenticate this exact fit')))
    with pytest.raises(ValueError, match='authenticate this exact fit'):
        upgrade.upgrade_copied_network_browser(asset, tmp_path / 'copy/browser/ALL_MODULES.html')


def test_network_browser_identity_must_match_completed_ledger_and_payload(tmp_path):
    native = _browser(tmp_path / 'collection/presentation/network/browser/ALL_MODULES.html')
    source = _browser(tmp_path / 'portable/browser/ALL_MODULES.html')
    root = tmp_path / 'collection'
    metadata = {'portable_root': str(native.parent.parent), 'relative_path': 'browser/ALL_MODULES.html'}
    row = {'path': native.relative_to(root).as_posix(), 'sha256': sha256(native.read_bytes()).hexdigest(), 'size_bytes': native.stat().st_size}
    receipt = {'outputs': [row]}
    upgrade._verify_registered_browser(source, metadata, root, receipt)
    source.write_text('<script id="data" type="application/json">{"members":["substituted-event"]}</script>')
    with pytest.raises(ValueError, match='numerical payload'):
        upgrade._verify_registered_browser(source, metadata, root, receipt)
    native.write_text(native.read_text() + 'changed source')
    with pytest.raises(ValueError, match='collection ledger'):
        upgrade._verify_registered_browser(source, metadata, root, receipt)


def test_native_table_plural_scope_is_preserved_in_generic_report(tmp_path):
    from katana_next.reports.engine import Builder
    table = tmp_path / 'original.tsv'
    table.write_text('event_key\tmodule\ne1\tM01\n', encoding='utf-8')
    builder = Builder({'project': {'project_id': 'fixture', 'contrasts': []}, 'artifacts': []}, tmp_path / 'report', 'stable_significant')
    record = builder.table({'path': str(table), 'artifact_id': 'native-member-table', 'analysis_id': 'programs.cosplice',
                            'contrast_ids': ['B_vs_A'], 'name': 'An arbitrary native name'})
    assert record['contrast_id'] == 'B_vs_A'
    assert table.read_text() == 'event_key\tmodule\ne1\tM01\n'


def test_legacy_missing_table_scopes_use_verified_record_ids_not_title_inference(tmp_path, monkeypatch):
    root = tmp_path / 'collection'
    root.mkdir()
    completion = root / 'proof.json'
    completion.write_text('{}')
    table = root / 'actual.tsv'
    table.write_text('module\tvalue\nM01\t0.5\n')
    digest = sha256(table.read_bytes()).hexdigest()
    artifacts = [{'artifact_id': 'native-id', 'analysis_id': 'programs.cosplice', 'name': 'Uninformative registered name',
                  'path': str(table), 'contrast_ids': ['B_vs_A']}]
    monkeypatch.setattr(upgrade, '_completed_inputs', lambda *_: (root, {}, completion, {}, None))
    import katana_next.flow.collect
    monkeypatch.setattr(katana_next.flow.collect, 'materialize_report_descriptor', lambda *_: {'artifacts': artifacts})
    manifest = {'network_browsers': [{'metadata': {'completion_scope': 'flow_verified_network_complete', 'completion_manifest': str(completion)}}],
                'source_files': [{'report_path': 'copied.tsv', 'source_sha256': digest}],
                'complete_results': [
                    {'artifact_id': 'native-id', 'name': 'Title contains no contrast', 'source': 'native.tsv', 'analysis_id': 'programs.cosplice'},
                    {'name': 'Uninformative registered name', 'source': 'copied.tsv', 'analysis_id': 'programs.cosplice'},
                    {'name': 'Looks like B_vs_A but has no registered identity', 'source': 'other.tsv', 'analysis_id': 'programs.cosplice'}]}
    result = upgrade.bind_legacy_network_table_scopes(manifest)
    assert len(result['bound_tables']) == 2
    assert [x.get('contrast_id') for x in manifest['complete_results']] == ['B_vs_A', 'B_vs_A', None]
    assert all(not x['metadata']['legacy_scope_binding']['title_parsed_for_contrast'] for x in manifest['complete_results'][:2])


def test_browser_upgrade_retains_compacted_exact_table_destinations(tmp_path):
    from katana_next.reports.network_views import refresh_completed_network_browser
    payload = {'members': [{'event_key': 'e1', 'module': 'M01', 'gene_id': 'g1', 'gene_symbol': 'G1'}],
               'samples': ['A', 'B'], 'psi': [[0.2, 0.7]], 'modules': ['M01'],
               'score_samples': ['A', 'B'], 'score_modules': ['M01'], 'scores': [[0.2], [0.7]],
               'conditions': ['control', 'treated'], 'edges': []}
    source = _browser(tmp_path / 'source/browser/ALL_MODULES.html', payload=payload)
    target = _browser(tmp_path / 'copy/browser/ALL_MODULES.html', payload=payload)
    links = '<a href="../../../results/exact-members.tsv">Complete members and continuous metrics</a><a href="../../../results/exact-background.tsv">Measured background</a><a href="../../../results/exact-input-genes.tsv">Module input genes</a>'
    target.write_text(target.read_text() + links, encoding='utf-8')
    before = source.read_bytes()
    refresh_completed_network_browser(source, target)
    text = target.read_text('utf-8')
    assert 'href="../../../results/exact-members.tsv"' in text
    assert 'href="../../../results/exact-background.tsv"' in text
    assert 'href="../../../results/exact-input-genes.tsv"' in text
    assert 'href="../tables/all_module_genes.tsv"' in text
    assert source.read_bytes() == before
