"""All hub measurements survive; heuristic candidate labels require opt-in."""
from dataclasses import replace
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from katana_next.science.hub_events import HubSettings, network_label_policy, summarize_hubs
from katana_next.workflows.hub_events import HubEventsWorkflow


def raw_metrics():
    """A complete fixed module, including weak but estimable members."""
    rows = []
    for bootstrap in range(21):
        for event in range(10):
            rows.append({
                'module': 'C1', 'condition': 'all', 'bootstrap': bootstrap,
                'event_key': f'event-{event}', 'module_events': 10,
                'biological_units': 24, 'kWithin': 10.0-event*.5,
                'kWithin_no_same_gene': 8.0-event*.4,
                'kME': .9 if event == 0 else .2,
                'pair_coverage': 1.0 if event == 0 else .4,
                'cross_gene_pair_coverage': 1.0 if event == 0 else .4,
                'supported_degree': 9 if event == 0 else 0,
                'cross_gene_supported_degree': 9 if event == 0 else 0,
                'gene_known': True, 'gene_annotation_fraction': 1.,
                'cross_gene_possible_pairs': 9,
            })
    return pd.DataFrame(rows)


def test_default_retains_every_event_and_complete_ranks_without_candidate_labels():
    raw = raw_metrics()
    original = raw.copy(deep=True)
    settings = HubSettings(bootstrap_iterations=20)
    metrics, boot = summarize_hubs(raw, settings)
    pd.testing.assert_frame_equal(raw, original)
    assert len(metrics) == 10 and len(boot) == 200
    assert metrics.rank_all_edges.tolist() == list(range(1, 11))
    assert metrics.rank_no_same_gene.tolist() == list(range(1, 11))
    assert metrics.label_status.eq('not_applied').all()
    assert metrics.hub_label_status.eq('not_applied').all()
    assert metrics.estimation_status.eq('estimable').all()
    for column in ('candidate_all_edges', 'candidate_no_same_gene', 'hub_candidate',
                   'pass_all_edges', 'pass_no_same_gene'):
        assert metrics[column].isna().all()
    assert boot.pass_all_edges.isna().all()
    # Weak kME, low pair coverage and no edges above the display threshold do not
    # erase a mathematically estimable member or remove it from the ranking.
    assert metrics.loc[metrics.event_key.eq('event-9'), 'kME'].item() == .2
    assert metrics.bootstrap_valid_fraction_all_edges.eq(1.).all()
    assert metrics.loc[metrics.event_key.eq('event-9'),
                       'bootstrap_preset_valid_fraction_all_edges'].item() == 0.


def test_opt_in_adds_labels_but_changes_no_continuous_metrics_rows_or_ranks():
    settings = HubSettings(bootstrap_iterations=20)
    default, default_boot = summarize_hubs(raw_metrics(), settings)
    labelled, labelled_boot = summarize_hubs(
        raw_metrics(), replace(settings, apply_threshold_labels=True))
    columns = ['module', 'condition', 'event_key', 'kME', 'kWithin',
               'kWithin_no_same_gene', 'rank_all_edges', 'rank_no_same_gene',
               'pair_coverage', 'cross_gene_pair_coverage',
               'preset_selection_frequency_all_edges',
               'bootstrap_valid_fraction_all_edges']
    pd.testing.assert_frame_equal(default[columns], labelled[columns])
    pd.testing.assert_frame_equal(default_boot[['event_key', 'bootstrap', 'rank_all_edges']],
                                  labelled_boot[['event_key', 'bootstrap', 'rank_all_edges']])
    assert labelled.label_status.eq('applied').all()
    assert labelled.hub_candidate.tolist() == [True] + [False]*9
    assert len(labelled) == 10


def test_preset_selection_frequency_is_explicitly_threshold_dependent_sensitivity():
    settings = HubSettings(bootstrap_iterations=20)
    raw = raw_metrics()
    raw['supported_degree'] = raw['cross_gene_supported_degree'] = 9
    strict, _ = summarize_hubs(raw, settings)
    loose, _ = summarize_hubs(raw, replace(
        settings, minimum_kme=.1, minimum_pair_coverage=.1, top_fraction=1.))
    for frame in (strict, loose):
        pd.testing.assert_series_equal(
            frame.bootstrap_frequency_all_edges,
            frame.preset_selection_frequency_all_edges,
            check_names=False)
        assert frame.hub_candidate.isna().all()
    # The full connectivity ranking is independent of which preset is displayed.
    assert strict.rank_all_edges.tolist() == loose.rank_all_edges.tolist()
    assert strict.loc[strict.event_key.eq('event-9'),'preset_selection_frequency_all_edges'].item() == 0.
    assert loose.loc[loose.event_key.eq('event-9'),'preset_selection_frequency_all_edges'].item() == 1.
    policy = network_label_policy(settings)
    assert policy['policy_id'] == 'network-descriptive-labels-opt-in-v1'
    assert policy['apply_threshold_labels'] is False
    assert policy['hard_filter'] is False
    assert 'sensitivity' in policy['preset_frequency_semantics']


@pytest.mark.parametrize('problem', ['sample_count', 'no_variance', 'missing_gene'])
def test_nonestimable_measurements_remain_explicit_without_false_candidate_claims(problem):
    raw = raw_metrics()
    if problem == 'sample_count':
        raw['biological_units'] = 6
        raw.loc[:, ['kWithin', 'kWithin_no_same_gene', 'kME']] = np.nan
    elif problem == 'no_variance':
        raw.loc[:, ['kWithin', 'kWithin_no_same_gene', 'kME']] = np.nan
    else:
        raw['gene_known'] = False
        raw['gene_annotation_fraction'] = 0.
        raw['cross_gene_possible_pairs'] = 0
        raw['cross_gene_pair_coverage'] = np.nan
        raw['kWithin_no_same_gene'] = np.nan
    for enabled in (False, True):
        metrics, boot = summarize_hubs(raw, HubSettings(
            bootstrap_iterations=20, apply_threshold_labels=enabled))
        assert len(metrics) == 10 and len(boot) == 200
        assert metrics.estimation_status.eq('not_estimable').all()
        assert metrics.hub_candidate.isna().all()
        assert metrics.label_status.eq('not_estimable' if enabled else 'not_applied').all()
        assert not metrics.estimation_reason.eq('estimable').any()
        if problem == 'missing_gene':
            assert metrics.all_edge_estimation_status.eq('estimable').all()
            assert metrics.cross_gene_estimation_status.eq('not_estimable').all()


def test_unavailable_bootstraps_are_not_reported_as_measured_zero_frequency():
    raw = raw_metrics().query('bootstrap == 0').copy()
    raw['biological_units'] = 6
    raw.loc[:, ['kWithin', 'kWithin_no_same_gene', 'kME']] = np.nan
    metrics, boot = summarize_hubs(raw, HubSettings(bootstrap_iterations=20))
    assert boot.empty
    assert metrics.preset_selection_frequency_all_edges.isna().all()
    assert metrics.bootstrap_valid_fraction_all_edges.isna().all()
    assert metrics.hub_candidate.isna().all()


def test_invalid_bootstraps_count_in_requested_denominator_and_all_rows_survive():
    raw = raw_metrics()
    raw.loc[raw.bootstrap.gt(4), ['kWithin', 'kWithin_no_same_gene', 'kME']] = np.nan
    metrics, boot = summarize_hubs(raw, HubSettings(
        bootstrap_iterations=20, apply_threshold_labels=True))
    center = metrics.loc[metrics.event_key.eq('event-0')].iloc[0]
    assert center.bootstrap_valid_fraction_all_edges == .2
    assert center.preset_selection_frequency_all_edges == .2
    assert not center.hub_candidate
    assert len(metrics) == 10 and len(boot) == 200


@pytest.mark.parametrize('value,expected', [(False, False), (True, True),
    ('false', False), (' TRUE ', True), (None, False)])
def test_boolean_control_does_not_treat_false_string_as_enabled(value, expected):
    parameters = HubEventsWorkflow().normalize_parameters({'apply_threshold_labels': value})
    assert parameters['apply_threshold_labels'] is expected


@pytest.mark.parametrize('value', ['disabled', '', 1, 0, [], {}])
def test_ambiguous_boolean_controls_are_rejected(value):
    with pytest.raises(ValueError, match='apply_threshold_labels'):
        HubEventsWorkflow().normalize_parameters({'apply_threshold_labels': value})


def test_workflow_records_policy_and_keeps_not_applied_summary_counts_missing(tmp_path, monkeypatch):
    from katana_next.full_fixture import create_full_fixture
    from katana_next.services.project_service import ProjectService
    from katana_next.workflows.base import ArtifactWriter
    from katana_next.workflows.programs import CoSpliceWorkflow
    from katana_next.workflows import hub_events as workflow_module

    context = ProjectService(application_root=tmp_path/'state').open_project(
        create_full_fixture(tmp_path/'project'))
    writer = ArtifactWriter(context, 'programs.cosplice', 'discovery-run', {}, CoSpliceWorkflow.version)
    source = writer.table('Modules', 'consensus_modules__assignments.parquet',
        pd.DataFrame({'event_key': ['event-0'], 'module': ['C1']}), entity_type='module')
    writer.table('PSI', 'matrix__psi.parquet',
        pd.DataFrame({'event_key': ['event-0'], 's1': [.5]}), entity_type='matrix')
    writer.table('Events', 'matrix__events.parquet',
        pd.DataFrame({'event_key': ['event-0'], 'gene_id': ['G1']}), entity_type='event')
    config = writer.directory/'resolved_cosplice.json'
    config.write_text(json.dumps({'network': {'method': 'pearson', 'soft_power': 6, 'signed': True}}))
    writer.existing_file('Configuration', config, entity_type='configuration')
    writer.finish()
    rscript = tmp_path/'Rscript'
    rscript.write_text('fixture only')
    context.manifest.metadata['rscript_path'] = str(rscript)

    scientific_calls=[]

    def fake_run(matrix, assignments, samples, annotation, settings, directory, rscript, **kwargs):
        scientific_calls.append(settings)
        metrics, boot = summarize_hubs(raw_metrics(), settings)
        for name in ('hub_plan.json', 'hub_input_psi.tsv', 'hub_input_events.tsv',
                     'hub_biological_units.tsv', 'technical_repeat_audit.json',
                     'R_session.txt', 'R_stdout.log', 'R_stderr.log', 'hub_raw_metrics.tsv'):
            (Path(directory)/name).write_text('fixture\n')
        pd.DataFrame(columns=['module', 'condition', 'event_a', 'event_b', 'adjacency']).to_csv(
            Path(directory)/'hub_display_edges.tsv', sep='\t', index=False)
        raw_metrics().to_csv(Path(directory)/'hub_raw_metrics.tsv',sep='\t',index=False)
        pd.DataFrame({'event_key':[f'event-{i}' for i in range(10)],
                      'gene':[f'G{i}' for i in range(10)],
                      'host_gene':[f'G{i}' for i in range(10)]}).to_csv(
            Path(directory)/'hub_input_events.tsv',sep='\t',index=False)
        return metrics, boot, {}

    monkeypatch.setattr(workflow_module, 'run_hub_analysis', fake_run)
    workflow = HubEventsWorkflow()
    parameters = workflow.normalize_parameters({'source_artifact_id': source.artifact_id,
                                                'bootstrap_iterations': 20})
    workflow.run(context, SimpleNamespace(report=lambda *args: None), 'hub-run', parameters)
    outputs = context.catalog.list_artifacts(context.manifest.project_id, analysis_id='programs.hub_events')
    for artifact in outputs:
        assert artifact.metadata['network_label_policy']['apply_threshold_labels'] is False
        assert artifact.metadata['network_label_policy']['hard_filter'] is False
    summary = pd.read_parquet(next(a.path for a in outputs if a.path.stem=='hub_events__summary'))
    assert summary.label_status.eq('not_applied').all()
    assert summary[['candidates_all_edges', 'candidates_no_same_gene',
                    'robust_hub_candidates']].isna().all().all()
    assert summary.module_events.tolist() == [10]
    original_bytes={a.path:a.path.read_bytes() for a in outputs}
    enabled=workflow.normalize_parameters({**parameters,'apply_threshold_labels':True})
    workflow.run(context,SimpleNamespace(report=lambda *args:None),'hub-labelled-view',enabled)
    assert len(scientific_calls)==1
    assert all(path.read_bytes()==value for path,value in original_bytes.items())
    labelled_outputs=[a for a in context.catalog.list_artifacts(
        context.manifest.project_id,analysis_id='programs.hub_events')
        if a.run_fingerprint=='hub-labelled-view']
    assert all(a.metadata['scientific_measurements_reused'] is True for a in labelled_outputs)
    labelled=pd.read_parquet(next(a.path for a in labelled_outputs if a.path.stem=='hub_events__metrics'))
    assert labelled.label_status.eq('applied').all()
    assert labelled.hub_candidate.tolist()==[True]+[False]*9
    assert len(labelled)==10

    # Changing only label thresholds likewise reuses the complete raw execution.
    custom=workflow.normalize_parameters({**enabled,'top_fraction':1.,'minimum_kme':.1,
                                         'minimum_pair_coverage':.1})
    workflow.run(context,SimpleNamespace(report=lambda *args:None),'hub-custom-label-view',custom)
    assert len(scientific_calls)==1
    # A completed marker alone is insufficient: altered raw bytes are not reused.
    for artifact in context.catalog.list_artifacts(context.manifest.project_id,
                                                  analysis_id='programs.hub_events'):
        if artifact.path.name=='hub_raw_metrics.tsv':
            artifact.path.write_bytes(artifact.path.read_bytes()+b'\n')
    workflow.run(context,SimpleNamespace(report=lambda *args:None),'hub-after-raw-change',enabled)
    assert len(scientific_calls)==2
    new_outputs=[a for a in context.catalog.list_artifacts(context.manifest.project_id,
                  analysis_id='programs.hub_events') if a.run_fingerprint=='hub-after-raw-change']
    assert new_outputs and all(a.metadata['scientific_measurements_reused'] is False for a in new_outputs)
