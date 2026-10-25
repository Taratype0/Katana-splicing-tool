"""Positive network-report acceptance using explicit completed-view fixtures.

Native fitting is covered by generated-fixture integration. These tests cover
the omission/identity boundary with real renderers and content-bound receipts.
"""
import copy
import json
from pathlib import Path
import shutil

import pytest
import pandas as pd

from katana_next.flow.network_acceptance import network_report_gaps
from katana_next.flow.presentation import build_presentation_assets
from katana_next.flow.runtime import write_json
from katana_next.reports.engine import ASSET_ROOT, Builder
from katana_next.reports.sources import load_source
from test_flow_presentation import _two_contrast_network_collection, _fixture_collection, _seal_fixture_collection


def _relative(value, root):
    if isinstance(value, dict):
        return {key: _relative(child, root) for key, child in value.items()}
    if isinstance(value, list):
        return [_relative(child, root) for child in value]
    if isinstance(value, str):
        path = Path(value)
        if path.is_absolute() and path.is_relative_to(root):
            return path.relative_to(root).as_posix()
    return value


@pytest.fixture
def complete_network_report(tmp_path, request):
    empty = getattr(request, 'param', False)
    initial = tmp_path / 'initial' if empty else tmp_path
    initial.mkdir(exist_ok=True)
    plan, receipts, transport, _ = _two_contrast_network_collection(initial)
    root = tmp_path / 'collection'
    if empty:
        for task in plan['tasks']:
            record = receipts[task['id']][1]
            evidence = record['execution_evidence']
            if task['adapter'] == 'network.fit':
                engine = Path(evidence['engine_root'])
                path = engine / '04_network/observed/event_module_assignments.tsv'
                frame = pd.read_csv(path, sep='\t'); frame['module'] = 'grey'
                frame.to_csv(path, sep='\t', index=False)
                path = engine / '04_network/observed/module_program_scores.tsv'
                scores = pd.read_csv(path, sep='\t', index_col=0)
                scores.iloc[:, :0].to_csv(path, sep='\t')
            elif task['adapter'] == 'network.bootstrap':
                directory = Path(evidence['engine_root']) / '04c_consensus_modules'
                path = directory / 'consensus_event_module_assignments.tsv'
                frame = pd.read_csv(path, sep='\t'); frame['module'] = 'grey'
                frame.to_csv(path, sep='\t', index=False)
            elif task['adapter'] in {'network.associate', 'network.annotate'}:
                views = evidence['precomputed_views']
                keys = ('full_tests', 'selected_tests') if task['adapter'] == 'network.associate' else ('tested',)
                for key in keys:
                    path = Path(views[key]); frame = pd.read_csv(path, sep='\t')
                    frame.iloc[:0].to_csv(path, sep='\t', index=False)
                if task['adapter'] == 'network.annotate':
                    path = Path(views['input_genes']); frame = pd.read_csv(path, sep='\t')
                    frame['module'] = 'grey'; frame.drop_duplicates().to_csv(path, sep='\t', index=False)
        products = {task['id']: ([Path(row['path']) for row in receipts[task['id']][1]['outputs']],
                                 receipts[task['id']][1]['execution_evidence']) for task in plan['tasks']}
        project = json.loads((initial / 'collection/project.json').read_text())
        plan, receipts, transport = _fixture_collection(root, project, plan['tasks'], products)
    result = build_presentation_assets(root, plan, receipts, transport)
    native = []
    for proof_path in sorted((root / 'presentation').glob('network-*-inputs.json')):
        proof = json.loads(proof_path.read_text())
        identity = proof['network_identity']
        cid = identity['contrast_ids'][0]
        overview = next(a for a in result['extra_assets']
                        if a['contrast_ids'] == [cid] and a['metadata'].get('network_view_role') == 'all_module_sizes')
        for role, stage, key, data_key in (
                ('dendrogram', 'network.fit', 'dendrogram', 'linkage'),
                ('bootstrap_confidence', 'network.bootstrap', 'confidence', 'bootstrap_members')):
            # Explicit plot/ledger fixture; no claim of a native fit here.
            image = root / 'presentation' / (cid + '-' + role + '.png')
            shutil.copyfile(overview['path'], image)
            native.append({'artifact_id': cid + '-' + role, 'analysis_id': 'programs.cosplice',
                'name': role, 'kind': 'file', 'path': str(image), 'status': 'complete',
                'contrast_ids': [cid], 'metadata': {'completion_verified': True,
                    'network_identity': identity, 'flow_stage': stage, 'flow_artifact_key': key,
                    'completion_scope': 'full_network_with_bootstrap',
                    'data_links': [{'role': 'plot_data', 'label': 'Fixture numerical source', 'path': str(root / proof[data_key]),
                                    'contrast_id': cid, 'unit': 'network fixture', 'policy': None}]}})
    write_json(root / 'artifacts.json', {'artifacts': _relative(native, root)})
    write_json(root / 'report-project.json', {'extra_assets': _relative(result['extra_assets'], root)})
    _seal_fixture_collection(root, plan, receipts)
    report = tmp_path / 'report'; report.mkdir()
    shutil.copytree(ASSET_ROOT, report / 'assets')
    descriptor = {'project_manifest': str(root / 'project.json'), 'artifacts': native,
                  'extra_assets': result['extra_assets'], 'report_require_plot_data': True}
    builder = Builder(load_source(descriptor), report, 'stable_significant')
    builder.collect_rich({'artifacts': [*native, *result['extra_assets']]})
    manifest = builder.write()
    return root, plan, manifest, report


def test_positive_gate_accepts_each_exact_completed_comparison(complete_network_report):
    root, plan, manifest, report = complete_network_report
    assert network_report_gaps(root, plan, manifest, report) == []
    browsers = [row for row in manifest['complete_results']
                if row.get('metadata', {}).get('network_view_role') == 'module_browser']
    assert {row['contrast_id'] for row in browsers} == {'B_vs_A', 'C_vs_A'}
    assert len({row['source_artifact_id'] for row in browsers}) == 2
    trees = [row for row in manifest['panels'] if row.get('metadata', {}).get('flow_artifact_key') == 'dendrogram']
    assert len(trees) == 2 and all(row.get('output_only') is True for row in trees)
    body = (report/'REPORT.html').read_text(encoding='utf-8')
    outputs = (report/'outputs/index.html').read_text(encoding='utf-8')
    for tree in trees:
        assert f'id="{tree["id"]}"' not in body
        assert '../'+tree['image'] in outputs
        assert (report/tree['image']).is_file()


def test_missing_network_tasks_or_catalogue_return_concrete_gaps(tmp_path):
    gaps = network_report_gaps(tmp_path, {'tasks': []}, {}, tmp_path / 'report')
    assert len(gaps) == 1 and gaps[0]['component'] == 'network'
    assert 'no planned network fit' in gaps[0]['reason']
    gaps = network_report_gaps(tmp_path, {'tasks': [{'id': 'fit', 'adapter': 'network.fit'}]}, {}, tmp_path / 'report')
    assert len(gaps) == 1 and 'catalogue is absent' in gaps[0]['reason']


@pytest.mark.parametrize('complete_network_report', [True], indirect=True)
def test_zero_assigned_modules_and_zero_testable_terms_do_not_require_invented_panels(complete_network_report):
    root, plan, manifest, report = complete_network_report
    assert network_report_gaps(root, plan, manifest, report) == []
    roles = {r.get('metadata', {}).get('network_view_role') for r in manifest['panels']}
    assert not roles & {'all_module_activity', 'matched_module_effects', 'module_enrichment_overview'}
    assert {'all_module_sizes', 'all_module_event_types'} <= roles


def test_positive_gate_rejects_silent_loss_of_whole_comparison(complete_network_report):
    root, plan, manifest, report = complete_network_report
    changed = copy.deepcopy(manifest)
    for key in ('panels', 'complete_results'):
        changed[key] = [r for r in changed[key] if r.get('contrast_id') != 'C_vs_A']
    changed['unavailable'] = []  # A missing disclosure must not conceal a loss.
    gaps = network_report_gaps(root, plan, changed, report)
    assert len(gaps) == 1 and gaps[0]['contrast_ids'] == ['C_vs_A']
    assert 'absent' in gaps[0]['reason']


def test_positive_gate_rejects_cross_contrast_browser_and_changed_plot_bytes(complete_network_report):
    root, plan, manifest, report = complete_network_report
    changed = copy.deepcopy(manifest)
    browsers = [r for r in changed['complete_results']
                if r.get('metadata', {}).get('network_view_role') == 'module_browser']
    browsers[1]['metadata']['network_identity'] = browsers[0]['metadata']['network_identity']
    assert any('identity differs' in gap['reason'] for gap in network_report_gaps(root, plan, changed, report))
    target = next(p for p in manifest['panels'] if p.get('metadata', {}).get('network_view_role') == 'all_module_sizes')
    (report / target['image']).write_bytes(b'changed fixture image')
    assert any('image differs' in gap['reason'] for gap in network_report_gaps(root, plan, manifest, report))


@pytest.mark.parametrize('mutation', ['missing', 'duplicate', 'main', 'foreign_comparison', 'wrong_image', 'wrong_linkage_scope'])
def test_output_only_dendrogram_still_requires_exact_complete_source_and_linkage(complete_network_report, mutation):
    root, plan, manifest, report = complete_network_report
    changed = copy.deepcopy(manifest)
    tree = next(row for row in changed['panels'] if row.get('metadata', {}).get('flow_artifact_key') == 'dendrogram')
    if mutation == 'missing':
        changed['panels'].remove(tree)
    elif mutation == 'duplicate':
        changed['panels'].append(copy.deepcopy(tree))
    elif mutation == 'main':
        tree['output_only'] = False
    elif mutation == 'foreign_comparison':
        tree['contrast_id'] = 'C_vs_A'
    elif mutation == 'wrong_image':
        altered = report/'altered-retained-tree.png'
        altered.write_bytes(b'Not the authenticated native source figure')
        tree['image'] = altered.name
    elif mutation == 'wrong_linkage_scope':
        tree['data_links'][0]['contrast_id'] = 'OTHER'
    assert network_report_gaps(root, plan, changed, report)


def test_required_network_overviews_cannot_be_hidden_using_output_only(complete_network_report):
    root, plan, manifest, report = complete_network_report
    changed = copy.deepcopy(manifest)
    panel = next(row for row in changed['panels'] if row.get('metadata', {}).get('network_view_role') == 'all_module_sizes')
    panel['output_only'] = True
    assert any('all_module_sizes' in gap['reason'] for gap in network_report_gaps(root, plan, changed, report))
