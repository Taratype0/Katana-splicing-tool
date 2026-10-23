"""Report axes are an immutable projection, not new biological inference."""
from copy import deepcopy
from hashlib import sha256
import json

import pandas as pd
import pytest

from katana_next.reports.sources import load_source, report_contrast_selection
from katana_next.reports.comparisons import _Views, _cross_patterns


def _native():
    return {'project_id': 'same-native-study', 'samples': [],
            'contrasts': [{'contrast_id': cid, 'numerator': {'sample_ids': []},
                           'denominator': {'sample_ids': []}} for cid in ['A', 'B', 'C', 'D', 'E']]}


def test_report_selection_keeps_native_hash_and_complete_originals(tmp_path):
    native = _native()
    path = tmp_path / 'project.json'
    path.write_text(json.dumps(native), encoding='utf-8')
    before = path.read_bytes()
    descriptor = {'project_manifest': str(path), 'artifacts': [],
        'report_contrast_ids': ['C', 'A'], 'report_contrast_groups': [
            {'group_id': 'declared-condition', 'label': 'Declared condition', 'contrast_ids': ['C', 'A']}],
        'report_native_figures_in_main': False, 'report_require_plot_data': True}
    original = deepcopy(descriptor)
    loaded = load_source(descriptor)
    assert [c['contrast_id'] for c in loaded['project']['contrasts']] == ['C', 'A']
    assert [c['contrast_id'] for c in loaded['native_project']['contrasts']] == ['A', 'B', 'C', 'D', 'E']
    assert loaded['report_selection']['excluded_contrast_ids'] == ['B', 'D', 'E']
    assert loaded['native_project_manifest_sha256'] == sha256(before).hexdigest()
    assert loaded['report_native_figures_in_main'] is False
    assert loaded['report_require_plot_data'] is True
    loaded['project']['contrasts'][0]['numerator']['sample_ids'].append('report-copy-only')
    assert loaded['native_project'] == native
    assert descriptor == original and path.read_bytes() == before


@pytest.mark.parametrize('changes', [
    {'report_contrast_ids': ['A', 'A']}, {'report_contrast_ids': ['unknown']},
    {'report_contrast_ids': []},
    {'report_contrast_ids': ['A', 'B'], 'report_contrast_groups': [{'group_id': 'g', 'contrast_ids': ['A']}]},
    {'report_contrast_groups': [{'group_id': 'g', 'contrast_ids': ['A', 'B', 'C', 'D', 'E']}]},
    {'report_contrast_ids': ['A'], 'report_contrast_groups': [{'group_id': 'g', 'contrast_ids': ['E']}]},
])
def test_invalid_report_projection_fails_without_changing_native(changes):
    native = _native()
    before = deepcopy(native)
    with pytest.raises(ValueError):
        report_contrast_selection(native, changes)
    assert native == before


def test_duplicate_native_ids_rejected_before_selection():
    native = _native()
    native['contrasts'].append(deepcopy(native['contrasts'][0]))
    with pytest.raises(ValueError, match='Native contrast IDs'):
        report_contrast_selection(native, {'report_contrast_ids': ['A']})


def test_narrowing_report_to_one_axis_does_not_authenticate_unscoped_native_table(tmp_path):
    table = tmp_path / 'unscoped.tsv'
    table.write_text('gene_id\tfdr\ng1\t0.01\n')
    view = _Views({'contrasts': [{'contrast_id': 'A'}], 'report_native_contrast_ids': ['A', 'B']}, [], tmp_path / 'report')
    assert view.frame({'path': str(table)}, 'A').empty
    assert view.audit[-1]['status'] == 'unavailable'


def test_grouped_upset_recomputes_axes_and_shows_every_actual_pattern(tmp_path):
    rows = []
    axes = ['A', 'B', 'C', 'D']
    for mask in range(1, 16):
        for i, cid in enumerate(axes):
            if mask & (1 << i):
                rows.append({'comparison_id': cid, 'event_id': 'event-' + str(mask),
                    'gene_id': 'gene-' + str(mask), 'dpsi': .2, 'fdr': 5e-310})
    rows += [{'comparison_id': cid, 'event_id': 'also-outside-group',
              'gene_id': 'extra', 'dpsi': -.3, 'fdr': .001} for cid in ['A', 'E']]
    source = pd.DataFrame(rows)
    original = source.copy(deep=True)
    view = _Views({}, [], tmp_path / 'report')
    result = _cross_patterns(view, source, 'event_id', axes, 'DAS-events', 'stable_significant',
                            group={'group_id': 'no-treatment', 'label': 'NO treatment'})
    assert result['patterns'] == 15 and result['objects'] == 16
    counts = pd.read_csv(view.root / result['intersection_counts'], sep='\t')
    assert len(counts) == 15 and counts.object_count.sum() == 16
    assert counts.loc[counts.membership_mask.eq(1), 'object_count'].iloc[0] == 2
    members = pd.read_parquet(view.root / result['intersection_members'])
    extra = members.loc[members.event_id.eq('also-outside-group')]
    assert extra.included_comparisons.tolist() == ['A']
    assert 'E' not in set(members.comparison_id)
    assert members.fdr.min() == 5e-310
    sizes = pd.read_csv(view.root / result['set_sizes'], sep='\t').set_index('comparison_id')
    assert sizes.object_count.to_dict() == {'A': 9, 'B': 8, 'C': 8, 'D': 8}
    panel = view.panels[-1]
    assert panel['comparison_ids'] == axes
    assert len(panel['data_links']) == 4
    assert all(row['comparison_ids'] == axes and row['policy'] == 'stable_significant' for row in panel['data_links'])
    assert all((view.root / row['path']).is_file() for row in panel['data_links'])
    assert '20 largest' not in panel['caption']
    pd.testing.assert_frame_equal(source, original)


def test_new_data_links_resolve_relative_to_the_declared_view_file(tmp_path):
    native = _native()
    manifest = tmp_path / 'project.json'
    manifest.write_text(json.dumps(native))
    digest = sha256(manifest.read_bytes()).hexdigest()
    folder = tmp_path / 'views'
    folder.mkdir()
    data = tmp_path / 'tables/values.tsv'
    data.parent.mkdir()
    data.write_text('x\n1\n')
    record = {'project_id': native['project_id'], 'status': 'complete',
        'source_project_manifest_sha256': digest, 'panels': [{'id': 'p', 'data_links': [
            {'role': 'plot_data', 'path': '../tables/values.tsv', 'contrast_id': 'A', 'policy': None, 'unit': 'gene'}]}]}
    view = folder / 'extra.json'
    view.write_text(json.dumps(record))
    loaded = load_source({'project_manifest': str(manifest), 'artifacts': [], 'extra_views': [str(view)]})
    assert loaded['extra_views'][0]['panels'][0]['data_links'][0]['path'] == str(data.resolve())


def test_twelve_selected_axes_form_three_groups_and_leave_sixteen_native_results_unchanged(tmp_path, monkeypatch):
    from katana_next.reports.comparisons import build_comparison_views
    ids = ['comparison-' + str(i) for i in range(16)]
    project = {'contrasts': [{'contrast_id': cid,
        'numerator': {'sample_ids': ['B1', 'B2', 'B3']},
        'denominator': {'sample_ids': ['A1', 'A2', 'A3']}} for cid in ids],
        'report_contrast_ids': ids[:12], 'report_contrast_groups': [
            {'group_id': label, 'label': label, 'contrast_ids': ids[start:start + 4]}
            for label, start in [('NO_DP', 0), ('DP', 4), ('DP_response', 8)]]}
    before_project = deepcopy(project)
    rows = []
    for cid in ids:
        for key in ['shared', 'specific-' + cid]:
            rows.append({'row_id': cid + '-' + key, 'comparison_id': cid, 'event_id': key,
                'event_type': 'SE', 'gene_id': 'gene-' + key, 'gene_symbol': 'G-' + key,
                'fdr': 5e-310, 'dpsi': .2, 'significant': True,
                'IJC_SAMPLE_1': '10,10,10', 'SJC_SAMPLE_1': '0,0,0',
                'IJC_SAMPLE_2': '10,10,10', 'SJC_SAMPLE_2': '0,0,0'})
    source = tmp_path / 'all-sixteen.parquet'
    pd.DataFrame(rows).to_parquet(source)
    before_sha = sha256(source.read_bytes()).hexdigest()
    artifact = {'analysis_id': 'differential.as', 'name': 'Differential AS events',
                'path': str(source), 'parameters': {}}
    def capture(self, fig, identity, chapter, title, comparison, caption, policy=None):
        fig.clear()
        self.panels.append({'id': identity, 'chapter': chapter, 'contrast_id': comparison or 'all',
                            'policy': policy, 'title': title, 'caption': caption})
    monkeypatch.setattr(_Views, 'figure', capture)
    result = build_comparison_views(project, [artifact], tmp_path / 'report')
    assert set(result['policy_data']['significant']) == set(ids[:12])
    assert set(result['cross_comparison']['significant']['groups']) == {'NO_DP', 'DP', 'DP_response'}
    assert all(p['contrast_id'] == 'all' or p['contrast_id'] in ids[:12] for p in result['panels'])
    assert all(p.get('data_links') for p in result['panels'])
    for policy in ['significant', 'stable_significant']:
        for group in project['report_contrast_groups']:
            cross = result['cross_comparison'][policy]['groups'][group['group_id']]
            assert cross['comparison_ids'] == group['contrast_ids']
            assert cross['patterns'] == 5 and cross['objects'] == 5
            members = pd.read_parquet(tmp_path / 'report' / cross['intersection_members'])
            assert set(members.comparison_id) == set(group['contrast_ids'])
            assert members.fdr.min() == 5e-310
    assert sha256(source.read_bytes()).hexdigest() == before_sha and project == before_project
