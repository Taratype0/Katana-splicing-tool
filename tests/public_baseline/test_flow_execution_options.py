"""Selected calculations, shared MAPP preparation and execution-only recovery."""
import copy
import json
import os
from pathlib import Path

import pytest

from katana_next.flow import file_io, managed, stages
from katana_next.flow.planner import create_plan, preflight_plan
from katana_next.flow.resources import render_resources
from katana_next.flow.runtime import file_record, verify_files
from katana_next.flow.selection import DEFAULTS, analysis_selection
from test_flow_managed_run import lifecycle
from test_flow_mapp_applicability import fixture as mapp_fixture, forbidden
from test_flow_planner import study


def test_disabled_calculations_are_absent_from_dag_and_requirements(study):
    options = {key: False for key in DEFAULTS}
    plan = create_plan(study, {"analyses": options, "design_factors": ["a", "b"],
                               "network": {"scope": "invalid-disabled-configuration"}})
    assert {t['workflow_id'] for t in plan['tasks']} == {
        'project.audit', 'quality.stable_as', 'differential.as', 'differential.deg',
        'differential.dtu', 'prioritize.evidence'}
    errors = preflight_plan(plan)['errors']
    assert not any(word in json.dumps(errors) for word in ('rbp_binding_peaks', 'mapp', 'genebody', 'interaction', 'network'))
    from katana_next.flow.cli import required_figure_chapters
    assert required_figure_chapters(plan) == {3}


@pytest.mark.parametrize('options', [{'netwrok': False}, {'network': 'false'}, {'mapp': 0}, []])
def test_misspelled_or_ambiguous_switches_fail(options):
    with pytest.raises(ValueError):
        analysis_selection({'analyses': options})


def test_default_remains_full_and_one_mapp_preparation_per_comparison(study):
    plan = create_plan(study, {})
    assert all(plan['configuration']['analyses'].values())
    preparations = [t for t in plan['tasks'] if t.get('adapter') == 'mapp.prepare']
    fits = [t for t in plan['tasks'] if t.get('adapter') == 'mapp.maei']
    assert len(preparations) == len(plan['contrast_ids'])
    for prep in preparations:
        children = [t for t in fits if t['dependencies'] == [prep['id']]]
        assert len(children) == 14
        assert all(t['scope']['contrast_ids'] == prep['scope']['contrast_ids'] for t in children)


def test_shared_preparation_preserves_all_window_eligibility(tmp_path, monkeypatch):
    from katana_next.flow import mapp_inputs
    spec = mapp_fixture(tmp_path)
    spec['parameters']['all_windows'] = True
    calls = []
    original = mapp_inputs.prepare_shared
    def counted(*args):
        calls.append(1)
        return original(*args)
    monkeypatch.setattr(mapp_inputs, 'prepare_shared', counted)
    monkeypatch.setattr(stages, '_run', forbidden)
    result = stages.execute_stage(spec, tmp_path / 'shared')
    assert len(calls) == 1 and result['evidence']['window_count'] == 14
    for window in mapp_inputs.native_windows():
        key = window['site'] + ':' + mapp_inputs.window_id(window['window_start'], window['window_end'])
        shared = json.loads(Path(result['evidence']['maei_windows'][key]['applicability_proof']).read_text())
        single = copy.deepcopy(spec)
        single['parameters'].pop('all_windows')
        single['parameters'].update(window)
        legacy = stages.execute_stage(single, tmp_path / ('single-' + key.replace(':', '-')))
        old = json.loads(Path(legacy['evidence']['applicability_proof']).read_text())
        for field in ('loss_counts', 'sample_ids', 'source_event_rows', 'comparison_event_rows',
                      'quantity_provenance', 'site', 'window_id', 'measured_events', 'eligible_events'):
            assert shared[field] == old[field]


def test_verification_reuses_identity_but_detects_same_size_modified_input(tmp_path, monkeypatch):
    source = tmp_path / 'data.tsv'
    source.write_text('abcd')
    record = file_record(source)
    calls = []
    original = file_io.hashlib.file_digest
    def counted(*args):
        calls.append(1)
        return original(*args)
    monkeypatch.setattr(file_io.hashlib, 'file_digest', counted)
    verify_files([record, record])
    assert not calls
    before = source.stat()
    source.write_text('wxyz')
    os.utime(source, ns=(before.st_atime_ns, before.st_mtime_ns))
    with pytest.raises(ValueError, match='changed'):
        verify_files([record])
    assert calls == [1]


def test_owned_product_survives_work_deletion_and_copy_fallback(tmp_path, monkeypatch):
    source = tmp_path / 'owned.tsv'; source.write_text('complete')
    target = tmp_path / 'retained.tsv'
    file_io.retain_owned_file(source, target)
    assert source.stat().st_ino == target.stat().st_ino
    source.unlink()
    assert target.read_text() == 'complete'
    monkeypatch.setattr(file_io.os, 'link', lambda *args: (_ for _ in ()).throw(OSError('different filesystem')))
    copy_path = tmp_path / 'copy.tsv'
    file_io.retain_owned_file(target, copy_path)
    assert copy_path.read_bytes() == target.read_bytes()


@pytest.mark.parametrize('value', [
    {'parameters': {'fdr': 0.1}}, {'defaults': {'queue': "main'; params.foo=1"}},
    {'processes': {'NETWORK_FIT': {'memory': 'infinite'}}}, {'defaults': {'cpus': True}},
    {'processes': {'*': {'time': '1h'}}},
])
def test_execution_resource_channel_cannot_change_science_or_inject_config(value):
    with pytest.raises(ValueError): render_resources(value)


def test_resource_change_resumes_failed_downstream_without_upstream(lifecycle, tmp_path):
    command, _, calls, fail = lifecycle
    resource = tmp_path / 'resources.json'
    resource.write_text(json.dumps({'defaults': {'queue': 'main'}, 'processes': {'NETWORK_FIT': {'memory': '72 GB'}}}))
    command += ['--resource-config', str(resource)]
    fail['downstream'] = True
    with pytest.raises(RuntimeError): managed.run(managed.parser().parse_args(command))
    resource.write_text(json.dumps({'defaults': {'queue': 'main'}, 'processes': {'NETWORK_FIT': {'memory': '96 GB'}}}))
    fail['downstream'] = False
    result = managed.run(managed.parser().parse_args([*command, '--resume']))
    assert result['status'] == 'completed'
    assert calls == [('upstream', False), ('downstream', False), ('downstream', True)]
    assert len(result['execution_history']) == 2
    assert result['execution_history'][0]['resource_id'] != result['execution_history'][1]['resource_id']


def test_selected_report_omits_disabled_sections_and_records_them(tmp_path, monkeypatch):
    from katana_next.reports.engine import ASSET_ROOT, Builder
    from katana_next.reports import landscape, comparisons
    from katana_next.reports.layout import render
    import shutil
    options = {key: False for key in DEFAULTS}
    source = {'project': {'project_id': 'selection-fixture', 'name': 'Selection fixture',
                         'samples': [], 'contrasts': [], 'annotation': {}},
              'artifacts': [], 'mode': 'bulk', 'report_analyses': options}
    monkeypatch.setattr(landscape, 'build_landscape_views', forbidden)
    monkeypatch.setattr(comparisons, 'build_comparison_views', forbidden)
    shutil.copytree(ASSET_ROOT, tmp_path / 'assets')
    builder = Builder(source, tmp_path, 'stable_significant')
    builder.create_archive = False
    builder.collect()
    result = builder.write()
    html = (tmp_path / 'REPORT.html').read_text(encoding='utf-8')
    assert 'id="chapter-1"' in html and 'id="chapter-3"' in html
    assert all(f'id="chapter-{n}"' not in html for n in (2, 4, 5, 6, 7))
    assert all(value == 'not_selected' for value in result['analyses'].values())
    assert not any(row.get('chapter') in (2, 4, 5, 6, 7) for row in result['unavailable'])


def test_shared_mapp_native_counter_matches_single_window(tmp_path):
    """Actual author counter, tiny generated sequence; no biological study run."""
    from katana_next.flow.bundled_runtime import verified_mapp_sources
    import pandas as pd
    spec = mapp_fixture(tmp_path, mapped=True)
    root, descriptor = verified_mapp_sources()
    spec['inputs'].update(native_root=str(root), counter_sha256=descriptor['counter_sha256'])
    spec['parameters']['all_windows'] = True
    shared = stages.execute_stage(spec, tmp_path / 'shared-counter')
    spec['parameters'].pop('all_windows')
    single = stages.execute_stage(spec, tmp_path / 'single-counter')
    from katana_next.flow.mapp_inputs import window_id
    key = spec['parameters']['site'] + ':' + window_id(-50, 0)
    observed = shared['evidence']['maei_windows'][key]
    expected = single['evidence']['maei_inputs']
    for name in ('inclusion_table', 'sitecount_matrix', 'design'):
        pd.testing.assert_frame_equal(pd.read_csv(observed[name], sep='\t'), pd.read_csv(expected[name], sep='\t'))
    assert observed['quantity_provenance'] == expected['quantity_provenance']


def test_core_selected_pipeline_reaches_actual_report_acceptance(study, tmp_path):
    """Execute installed analysis code on generated tables; no public/study rerun."""
    from katana_next.flow.model import write_json
    from katana_next.flow.runtime import execute_task
    from katana_next.flow.collect import collect_results
    from katana_next.flow.cli import build_complete_report
    plan = create_plan(study, {'design_factors': ['condition'], 'analyses': dict.fromkeys(DEFAULTS, False)})
    planned = write_json(tmp_path / 'plan.json', plan)
    receipts = {}
    for task in sorted(plan['tasks'], key=lambda task: task['level']):
        receipts[task['id']] = execute_task(planned, task['id'], tmp_path / 'tasks' / task['id'],
            dependency_receipts=[receipts[name] for name in task['dependencies']])
    collection = tmp_path / 'selected-collection'
    collect_results(planned, list(receipts.values()), collection)
    result = build_complete_report(collection, tmp_path / 'selected-report')
    assert result['status'] == 'completed' and result['verified_calculation_tasks'] == 6
    assert not any(result['analyses'].values())


def test_import_skips_only_resources_of_disabled_analyses(tmp_path):
    from test_flow_import import bundle
    from katana_next.flow.importer import import_bundle
    from katana_next.flow.model import read_json, write_json
    root, _, template = bundle(tmp_path)
    value = read_json(template)
    value['assets'] = [{'asset_id': 'motifs', 'kind': 'rbp_motif_database', 'path': 'missing.meme'}]
    write_json(template, value)
    with pytest.raises(ValueError, match='resource unavailable'):
        import_bundle(root, template, tmp_path / 'default')
    project = import_bundle(root, template, tmp_path / 'selected',
        analysis_configuration={'analyses': {'motiftrace': False, 'motif_binding': False}})
    assert all(asset['asset_id'] != 'motifs' for asset in read_json(project)['assets'])
    assert read_json(project.parent / 'import.receipt.json')['not_selected_resource_ids'] == ['motifs']
    value['assets'].append({'asset_id': 'counts', 'kind': 'normalized_counts', 'path': 'missing.tsv'})
    write_json(template, value)
    with pytest.raises(ValueError, match='resource unavailable: counts'):
        import_bundle(root, template, tmp_path / 'missing-core',
            analysis_configuration={'analyses': dict.fromkeys(DEFAULTS, False)})
