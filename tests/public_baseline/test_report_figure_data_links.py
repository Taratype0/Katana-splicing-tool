"""Every declared plot-data link has real, portable, correctly scoped evidence."""
from pathlib import Path
import copy
import hashlib
import json
import shutil

import pandas as pd
import pytest

from katana_next.reports.engine import ASSET_ROOT, Builder


def _builder(tmp_path, **flags):
    output = tmp_path / 'report'
    output.mkdir(parents=True)
    project = {'project_id': 'study', 'name': 'Study', 'samples': [], 'annotation': {},
               'contrasts': [{'contrast_id': 'B_vs_A', 'label': 'B versus A',
                              'numerator': {'label': 'B', 'sample_ids': []},
                              'denominator': {'label': 'A', 'sample_ids': []}}]}
    return Builder({'project': project, 'artifacts': [], 'mode': 'bulk', **flags}, output, 'stable_significant')


def _panel():
    return {'id': 'exact-chart', 'title': 'Actual event chart', 'chapter': 3,
            'contrast_id': 'B_vs_A', 'policy': 'stable_significant'}


def _link(path, **updates):
    return {'role': 'plot_data', 'path': str(path), 'label': 'Exact plot data',
            'contrast_id': 'B_vs_A', 'policy': 'stable_significant', 'unit': 'exact_event', **updates}


def _data(path, comparison='B_vs_A', policy='stable_significant'):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text('comparison_id\treport_policy\tevent_id\tfdr\tdpsi\n' +
                    f'{comparison}\t{policy}\tevent-1\t3.4e-310\t-0.125\n', encoding='utf-8')
    return path


def test_external_plot_table_is_copied_byte_exact_with_scope_and_hash_evidence(tmp_path):
    builder = _builder(tmp_path)
    source = _data(tmp_path / 'existing/unrelated-file-name.tsv')
    before = source.read_bytes(), source.stat().st_mtime_ns
    raw = _link(source)
    bound = builder.bind_data_links(_panel(), [raw])[0]
    assert not Path(bound['path']).is_absolute()
    assert (builder.output / bound['path']).read_bytes() == before[0]
    assert bound['source_sha256'] == hashlib.sha256(before[0]).hexdigest() == bound['portable_sha256']
    assert builder.data_link_audit[0]['scope_checks'] == {'comparison_values': ['B_vs_A'], 'policy_values': ['stable_significant']}
    assert (source.read_bytes(), source.stat().st_mtime_ns) == before
    assert raw == _link(source)


def test_existing_report_plot_data_is_not_copied_to_itself(tmp_path):
    builder = _builder(tmp_path)
    source = _data(builder.output / 'outputs/exact.tsv')
    before = source.read_bytes(), source.stat().st_mtime_ns
    linked = builder.bind_data_links(_panel(), [_link(source)])[0]
    assert linked['path'] == 'outputs/exact.tsv'
    assert not builder.audit
    assert (source.read_bytes(), source.stat().st_mtime_ns) == before


@pytest.mark.parametrize('changed', [{'contrast_id': 'Other_vs_A'}, {'policy': 'significant'}, {'project_id': 'another-study'}])
def test_declared_foreign_scope_is_rejected(tmp_path, changed):
    builder = _builder(tmp_path)
    path = _data(tmp_path / 'existing.tsv')
    with pytest.raises(ValueError, match='different'):
        builder.bind_data_links(_panel(), [_link(path, **changed)])


@pytest.mark.parametrize('comparison,policy', [('Other_vs_A', 'stable_significant'), ('all', 'stable_significant'), ('B_vs_A', 'significant')])
def test_matching_metadata_cannot_disguise_foreign_actual_plot_rows(tmp_path, comparison, policy):
    builder = _builder(tmp_path)
    source = _data(tmp_path / 'wrong-actual-rows.tsv', comparison, policy)
    with pytest.raises(ValueError, match='Plot data'):
        builder.bind_data_links(_panel(), [_link(source)])


def test_missing_or_unfingerprinted_placeholder_cannot_be_reported_as_data(tmp_path):
    builder = _builder(tmp_path)
    with pytest.raises(FileNotFoundError):
        builder.bind_data_links(_panel(), [_link(tmp_path / 'absent.tsv')])
    blank = tmp_path / 'blank.tsv'
    blank.write_bytes(b'')
    with pytest.raises(ValueError, match='empty'):
        builder.bind_data_links(_panel(), [_link(blank)])
    actual = _data(tmp_path / 'native.tsv')
    with pytest.raises(ValueError, match='SHA'):
        builder.bind_data_links(_panel(), [_link(actual, source_sha256='fabricated')])


def test_relative_traversal_does_not_silently_collect_another_file(tmp_path):
    builder = _builder(tmp_path)
    _data(tmp_path / 'outside.tsv')
    with pytest.raises(ValueError, match='escapes'):
        builder.bind_data_links(_panel(), [_link('../outside.tsv')])


def test_rendered_chart_has_specific_links_and_manifest_binding_not_total_outputs_substitute(tmp_path):
    builder = _builder(tmp_path, report_require_plot_data=True)
    shutil.copytree(ASSET_ROOT, builder.output / 'assets')
    source = _data(tmp_path / 'actual.tsv')
    image = tmp_path / 'figure.svg'
    image.write_text('<svg xmlns="http://www.w3.org/2000/svg"/>', encoding='utf-8')
    panel = {**_panel(), 'path': str(image), 'data_links': [_link(source)]}
    original = copy.deepcopy(panel)
    builder.add_views({'panels': [panel]})
    manifest = builder.write()
    page = (builder.output / 'REPORT.html').read_text(encoding='utf-8')
    link = manifest['panels'][0]['data_links'][0]
    assert f'data-data-role="plot_data" href="{link["path"]}"' in page
    assert 'Exact plot data' in page
    assert manifest['figure_data_links']['main_figures_with_plot_data'] == 1
    assert not manifest['figure_data_links']['unbound_main_figures']
    assert panel == original


def test_full_master_is_not_accepted_as_exact_plot_data_in_strict_report(tmp_path):
    builder = _builder(tmp_path, report_require_plot_data=True)
    source = _data(tmp_path / 'master.tsv')
    image = tmp_path / 'figure.svg'
    image.write_text('<svg xmlns="http://www.w3.org/2000/svg"/>', encoding='utf-8')
    builder.add_views({'panels': [{**_panel(), 'path': str(image),
                                 'data_links': [_link(source, role='full_results')]}]})
    with pytest.raises(ValueError, match='requires exact plot_data'):
        builder.write()


def test_original_native_figures_can_remain_output_only_without_fake_data(tmp_path):
    builder = _builder(tmp_path, report_native_figures_in_main=False, report_require_plot_data=True)
    image = tmp_path / 'old.svg'
    image.write_text('<svg xmlns="http://www.w3.org/2000/svg"/>', encoding='utf-8')
    record = builder.figure({'artifact_id': 'old', 'name': 'Native QC', 'path': str(image), 'analysis_id': 'quality.reads'}, 1)
    assert record['output_only'] and not record['data_links']
    assert builder.figure_data_coverage()['main_figures'] == 0


def test_unselected_comparison_remains_output_only_instead_of_leaking_into_main(tmp_path):
    builder = _builder(tmp_path)
    image = tmp_path / 'extra.svg'
    image.write_text('<svg xmlns="http://www.w3.org/2000/svg"/>', encoding='utf-8')
    builder.add_views({'panels': [{**_panel(), 'path': str(image), 'contrast_id': 'Excluded_vs_A'}]})
    assert builder.panels[0]['output_only']


def test_explicit_qc_metadata_links_are_carried_to_actual_figure(tmp_path):
    builder = _builder(tmp_path)
    source = _data(tmp_path / 'qc.tsv')
    image = tmp_path / 'qc.svg'
    image.write_text('<svg xmlns="http://www.w3.org/2000/svg"/>', encoding='utf-8')
    artifact = {'artifact_id': 'qc', 'name': 'Bound QC', 'path': str(image), 'contrast_ids': ['B_vs_A'],
                'metadata': {'event_set': 'stable_significant', 'data_links': [_link(source)]}}
    original = copy.deepcopy(artifact)
    panel = builder.figure(artifact, 1)
    assert panel['data_links'][0]['role'] == 'plot_data'
    assert panel['data_links'][0]['path'].startswith('outputs/figure-data/')
    assert artifact == original


def test_external_browser_copies_explicit_local_dependencies_without_source_changes(tmp_path):
    builder = _builder(tmp_path)
    bundle = tmp_path / 'browser-source'
    bundle.mkdir()
    browser = bundle / 'view.html'
    browser.write_text('<html><link rel="stylesheet" href="style.css"><a href="data.tsv">All rows</a></html>', encoding='utf-8')
    (bundle / 'style.css').write_text('body {color:#2863A5}', encoding='utf-8')
    (bundle / 'data.tsv').write_text('event_id\ne1\n', encoding='utf-8')
    originals = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in bundle.iterdir()}
    result = builder.bind_data_links(_panel(), [_link(browser, role='full_results', portable_root=str(bundle))])[0]
    copied = builder.output / result['path']
    assert copied.is_file() and (copied.parent / 'style.css').is_file() and (copied.parent / 'data.tsv').is_file()
    assert all((p.read_bytes(), p.stat().st_mtime_ns) == before for p, before in originals.items())


def test_browser_with_missing_dependency_is_not_falsely_marked_portable(tmp_path):
    builder = _builder(tmp_path)
    browser = tmp_path / 'broken.html'
    browser.write_text('<html><a href="absent.tsv">All rows</a></html>', encoding='utf-8')
    with pytest.raises(ValueError, match='dependencies'):
        builder.bind_data_links(_panel(), [_link(browser, role='full_results')])


def test_identical_browser_html_keeps_distinct_comparison_dependencies_isolated(tmp_path):
    builder = _builder(tmp_path)
    builder.source['project']['contrasts'].append({'contrast_id': 'C_vs_A'})
    markup = '<html><a href="data.tsv">All rows</a><script src="view.js"></script></html>'
    originals, bound = {}, []
    for comparison in ('B_vs_A', 'C_vs_A'):
        bundle = tmp_path / 'completed' / comparison
        bundle.mkdir(parents=True)
        browser = bundle / 'ALL_MODELS.html'
        browser.write_text(markup, encoding='utf-8')
        _data(bundle / 'data.tsv', comparison)
        (bundle / 'view.js').write_text('window.comparison=' + json.dumps(comparison), encoding='utf-8')
        originals.update({p: (p.read_bytes(), p.stat().st_mtime_ns) for p in bundle.iterdir()})
        panel = {**_panel(), 'id': comparison + '-panel', 'contrast_id': comparison}
        link = _link(browser, role='full_results', portable_root=str(bundle), contrast_id=comparison)
        result = builder.bind_data_links(panel, [link])[0]
        # Repeated references reuse this bundle without rewriting a prior copy.
        assert builder.bind_data_links(panel, [link])[0]['path'] == result['path']
        bound.append(result)
    assert bound[0]['source_sha256'] == bound[1]['source_sha256']
    assert bound[0]['path'] != bound[1]['path']
    for comparison, result in zip(('B_vs_A', 'C_vs_A'), bound):
        copied = builder.output / result['path']
        assert copied.read_text() == markup
        frame = pd.read_csv(copied.parent / 'data.tsv', sep='\t')
        assert frame.comparison_id.tolist() == [comparison]
        assert (copied.parent / 'view.js').read_text() == 'window.comparison=' + json.dumps(comparison)
    assert all((p.read_bytes(), p.stat().st_mtime_ns) == before for p, before in originals.items())
    assert len([a for a in builder.derived_audits if a['schema'] == 'katana.rich_asset_copy/v1']) == 2


def test_browser_cache_respects_same_file_with_different_declared_portable_roots(tmp_path):
    builder = _builder(tmp_path)
    root = tmp_path / 'collection'
    bundle = root / 'comparison'
    bundle.mkdir(parents=True)
    browser = bundle / 'view.html'
    browser.write_text('<html><a href="data.tsv">Full table</a></html>', encoding='utf-8')
    _data(bundle / 'data.tsv')
    first = builder.bind_data_links(_panel(), [_link(browser, role='full_results', portable_root=str(bundle))])[0]
    second = builder.bind_data_links(_panel(), [_link(browser, role='full_results', portable_root=str(root))])[0]
    assert first['path'] != second['path']
    assert (builder.output / first['path']).parent.name != 'comparison'
    assert (builder.output / second['path']).parent.name == 'comparison'
    assert (builder.output / first['path']).read_bytes() == (builder.output / second['path']).read_bytes()
    assert ((builder.output / first['path']).parent / 'data.tsv').read_bytes() == ((builder.output / second['path']).parent / 'data.tsv').read_bytes()


def test_grouped_upset_compound_plot_has_actual_patterns_and_all_set_sizes(tmp_path):
    from katana_next.reports.comparisons import _Views, _cross_patterns
    builder = _builder(tmp_path, report_require_plot_data=True)
    axes = ['A', 'B', 'C', 'D']
    builder.source['project']['contrasts'] = [{'contrast_id': cid} for cid in axes]
    rows = [{'comparison_id': cid, 'event_id': event, 'gene_id': event, 'dpsi': .2, 'fdr': .01}
            for cid, event in [('A', 'shared'), ('B', 'shared'), ('C', 'unique-C')]]
    producer = _Views(builder.source['project'], [], builder.output)
    _cross_patterns(producer, pd.DataFrame(rows), 'event_id', axes, 'DAS-events', 'stable_significant',
                    group={'group_id': 'declared', 'label': 'Declared exact comparisons'})
    builder.add_views({'panels': producer.panels, 'tables': producer.tables})
    bound = builder.panels[0]
    exact = [link for link in bound['data_links'] if link['role'] == 'plot_data']
    assert len(exact) == 2
    assert all(link['comparison_ids'] == axes and (builder.output / link['path']).is_file() for link in exact)
    sizes = next(link for link in exact if link['unit'] == 'unique_exact_object')
    counts = pd.read_csv(builder.output / sizes['path'], sep='\t').set_index('comparison_id').object_count
    assert counts.to_dict() == {'A': 1, 'B': 1, 'C': 1, 'D': 0}
    assert not builder.figure_data_coverage()['unbound_main_figures']


@pytest.mark.parametrize('ids', [['Foreign_vs_A'], ['B_vs_A', 'B_vs_A']])
def test_global_panels_cannot_authenticate_foreign_or_duplicate_axes_by_agreeing_with_the_link(tmp_path, ids):
    builder = _builder(tmp_path)
    source = _data(tmp_path / 'actual.tsv')
    panel = {**_panel(), 'contrast_id': 'all', 'comparison_ids': ids}
    with pytest.raises(ValueError, match='comparison'):
        builder.bind_data_links(panel, [_link(source, contrast_id='all', comparison_ids=ids)])


def test_compound_plot_links_keep_the_declared_axis_order(tmp_path):
    builder = _builder(tmp_path)
    builder.source['project']['contrasts'].append({'contrast_id': 'C_vs_A'})
    source = _data(tmp_path / 'actual.tsv')
    panel = {**_panel(), 'contrast_id': 'all', 'comparison_ids': ['B_vs_A', 'C_vs_A']}
    with pytest.raises(ValueError, match='axis order'):
        builder.bind_data_links(panel, [_link(source, contrast_id='all', comparison_ids=['C_vs_A', 'B_vs_A'])])


@pytest.mark.parametrize('workflow,mode', [('prioritize.motiftrace_consequences', 'rs_required'),
                                          ('prioritize.motiftrace_usage', 'fixed_completed_result')])
def test_real_rich_adapter_keeps_rs_measurement_bound_while_display_remains_always_visible(tmp_path, workflow, mode):
    builder = _builder(tmp_path)
    source = _data(tmp_path / 'actual-rs-sequence.tsv')
    image = tmp_path / 'actual-sequence.svg'
    image.write_text('<svg xmlns="http://www.w3.org/2000/svg"/>', encoding='utf-8')
    artifact = {'artifact_id': 'existing-sequence', 'analysis_id': workflow, 'path': str(image),
                'name': 'Actual completed MotifTrace sequence', 'status': 'complete', 'contrast_ids': ['B_vs_A'],
                'parameters': {'event_set': 'stable_significant'},
                'metadata': {'data_links': [_link(source)], 'read_support_policy': 'every selected sample IJC+SJC>=10'}}
    original = copy.deepcopy(artifact)
    builder.collect_rich({'artifacts': [artifact]})
    panel = next(p for p in builder.panels if p['title'].startswith('Actual completed'))
    assert panel['metadata']['das_scope']['mode'] == mode
    assert panel['policy'] is None
    assert panel['data_links'][0]['policy'] == 'stable_significant'
    assert builder.data_link_audit[0]['display_policy'] is None
    assert builder.data_link_audit[0]['measurement_policy'] == 'stable_significant'
    assert artifact == original
    from katana_next.reports.layout import _measurement_note, _scope
    assert 'scope' in _measurement_note(panel)
    assert 'data-policy="common"' in _scope(panel)
    artifact['metadata']['data_links'] = [_link(source, policy='significant')]
    with pytest.raises(ValueError, match='different DAS policy'):
        _builder(tmp_path / 'negative').collect_rich({'artifacts': [artifact]})
