"""Annotation choices never filter a partition or silently refit a network."""
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from katana_cosplice.config import CoSpliceConfig, ConsensusConfig, NetworkConfig
from katana_cosplice.pipeline import CoSplicePipeline
from katana_cosplice.utils import file_signature, fingerprint
from katana_next.science.network_labels import member_label_view, network_label_policy
from katana_next.workflows import programs, psi_network


def members():
    return pd.DataFrame({'event_key': ['e1', 'e2', 'e3'], 'consensus_module': ['C1', 'C1', 'grey'],
        'membership_role': ['core', 'uncertain', 'grey'], 'membership_evidence': ['coassignment', 'insufficient_confidence', 'unassigned'],
        'consensus_membership_confidence': [.8, .4, 0.], 'consensus_confidence_margin': [.5, .1, 0.],
        'bootstrap_mean_abs_kME': [.9, .6, 0.], 'bootstrap_q10_abs_kME': [.7, .3, 0.],
        'bootstrap_kME_sign_consistency': [1., .8, 0.],
        'split_half_membership_role_a': ['core', 'uncertain', 'grey'],
        'split_half_membership_role_b': ['core', 'peripheral', 'grey'],
        'split_half_core_confirmed': [True, False, False],
        'continuous_membership_confidence': [.9, .2, 0.],
        'continuous_confidence_rank_within_module': [1., 2., float('nan')],
        'FDR': [.04, .3, float('nan')]})


def test_default_labels_are_not_applied_and_complete_continuous_ranks_are_retained():
    raw = members()
    before = raw.copy(deep=True)
    off, on = member_label_view(raw), member_label_view(raw, True)
    pd.testing.assert_frame_equal(raw, before)
    assert off.event_key.tolist() == on.event_key.tolist() == raw.event_key.tolist()
    assert off.membership_role.eq('not_applied').all()
    assert off.split_half_core_confirmed.isna().all()
    assert on.membership_role.tolist() == ['core', 'uncertain', 'grey']
    for column in ('consensus_module', 'consensus_membership_confidence', 'bootstrap_mean_abs_kME',
                   'consensus_confidence_margin', 'FDR', 'coassignment_rank_within_module',
                   'kme_rank_within_module', 'margin_rank_within_module'):
        pd.testing.assert_series_equal(off[column], on[column])
    assert off.coassignment_rank_within_module.iloc[:2].tolist() == [1., 2.]
    assert pd.isna(off.coassignment_rank_within_module.iloc[2])
    assert 'continuous_membership_confidence' not in off
    assert off.legacy_threshold_dependent_membership_score.tolist() == [.9, .2, 0.]


def test_threshold_dependent_composite_and_aggregates_are_explicitly_legacy():
    from katana_cosplice.membership_validation import add_continuous_membership_confidence
    first = members()
    changed = first.copy()
    changed['split_half_membership_role_a'] = 'uncertain'
    changed['split_half_membership_role_b'] = 'uncertain'
    a = member_label_view(add_continuous_membership_confidence(first))
    b = member_label_view(add_continuous_membership_confidence(changed))
    assert not a.legacy_threshold_dependent_membership_score.equals(b.legacy_threshold_dependent_membership_score)
    pd.testing.assert_series_equal(a.coassignment_rank_within_module, b.coassignment_rank_within_module)
    summary = member_label_view(pd.DataFrame({'n_events': [37], 'n_core_events': [20],
                                            'mean_continuous_confidence': [.7]}))
    assert summary.n_events.tolist() == [37]
    assert summary.n_core_events.isna().all()
    assert summary.legacy_threshold_dependent_mean_confidence.tolist() == [.7]
    residual = pd.DataFrame({'observed_residual_program_class': ['residual_preserved_intrinsic_program'],
                              'FDR': [.02], 'jaccard': [.9], 'recall': [.8]})
    transformed = member_label_view(residual)
    assert transformed.observed_residual_program_class.eq('not_applied').all()
    pd.testing.assert_frame_equal(transformed[['FDR', 'jaccard', 'recall']], residual[['FDR', 'jaccard', 'recall']])


def test_switch_is_boolean_and_refinement_never_becomes_optional_label_mode(tmp_path):
    with pytest.raises(ValueError, match='true or false'):
        network_label_policy('false')
    config = CoSpliceConfig('fixture', None, tmp_path / 'samples.tsv', tmp_path / 'out')
    config.consensus.clustering_method = 'backbone_refined_consensus'
    for enabled in (False, True):
        programs.configure_production_membership(config, {'apply_threshold_labels': enabled})
        assert config.consensus.clustering_method == 'backbone_confidence_only'
        assert config.consensus.apply_threshold_labels is enabled
    with pytest.raises(ValueError, match='true or false'):
        programs.configure_production_membership(config, {'apply_threshold_labels': 'false'})


def raw_outputs(tmp_path):
    root = tmp_path / 'out' / '04c_consensus_modules'
    root.mkdir(parents=True)
    paths = {'assignments': root / 'consensus_event_module_assignments.tsv',
             'modules': root / 'consensus_module_summary.tsv',
             'parameters': root / 'parameters.json'}
    members().to_csv(paths['assignments'], sep='\t', index=False)
    pd.DataFrame({'consensus_module': ['C1'], 'n_events': [2], 'n_core_events': [1]}).to_csv(paths['modules'], sep='\t', index=False)
    paths['parameters'].write_text('{"legacy_raw": true}', encoding='utf-8')
    return paths


def signatures(paths):
    return {str(path): (hashlib.sha256(path.read_bytes()).hexdigest(), path.stat().st_mtime_ns) for path in paths}


def test_cli_actual_entrypoint_exports_two_label_views_without_refitting_or_mutating_raw(tmp_path, monkeypatch):
    from katana_cosplice.cli import main
    from katana_cosplice import __version__

    raw = raw_outputs(tmp_path)
    config = CoSpliceConfig('fixture', None, tmp_path / 'samples.tsv', tmp_path / 'out')
    settings = {'consensus': asdict(config.consensus), 'network': asdict(config.network)}
    old_settings = {**settings, 'consensus': {k: v for k, v in settings['consensus'].items() if k != 'apply_threshold_labels'}}
    old_hash = fingerprint({'schema': 1, 'network_implementation': 'biological-blocks-coverage-audit-20260921',
                           'version': __version__, 'stage': 'consensus_modules', 'settings': old_settings,
                           'inputs': [file_signature(raw['assignments'])]})
    observations = []

    def cached_science(self, force=False):
        assert not force
        actual = self._fingerprint('consensus_modules', [raw['assignments']],
            {'consensus': asdict(self.config.consensus), 'network': asdict(self.config.network)})
        assert actual == old_hash  # A mismatch would invoke a new scientific fit.
        observations.append(actual)
        return raw

    monkeypatch.setattr(CoSplicePipeline, 'consensus_modules', cached_science)
    before = signatures(raw.values())
    config_path = tmp_path / 'config.yaml'
    for enabled in (False, True):
        config.consensus.apply_threshold_labels = enabled
        config.write(config_path)
        assert main(['consensus-modules', '--config', str(config_path)]) == 0
    manifests = list((config.output_dir / 'label_views').rglob('parameters.json'))
    assert len(manifests) == 2 and len(observations) == 2
    views = {}
    for manifest in manifests:
        record = json.loads(manifest.read_text())
        enabled = record['network_label_policy']['apply_threshold_labels']
        assert record['network_label_policy']['hard_filter'] is False
        views[enabled] = pd.read_csv(record['view_files']['assignments']['path'], sep='\t')
    assert views[False].membership_role.eq('not_applied').all()
    assert views[True].membership_role.tolist() == ['core', 'uncertain', 'grey']
    pd.testing.assert_series_equal(views[False].consensus_membership_confidence, views[True].consensus_membership_confidence)
    assert before == signatures(raw.values())


def test_cli_rejects_historical_refined_partition_before_scientific_action(tmp_path, monkeypatch):
    from katana_cosplice.cli import main
    config = CoSpliceConfig('fixture', None, tmp_path / 'samples.tsv', tmp_path / 'out')
    config.consensus.clustering_method = 'backbone_refined_consensus'
    config_path = config.write(tmp_path / 'refined.yaml')
    monkeypatch.setattr(CoSplicePipeline, 'consensus_modules', lambda *_a, **_k: pytest.fail('Scientific action must not start'))
    with pytest.raises(ValueError, match='requires backbone_confidence_only'):
        main(['consensus-modules', '--config', str(config_path)])
    assert not list(config.output_dir.rglob('label_views'))


def test_cli_companion_preserves_raw_report_and_links_actual_isolated_table_views(tmp_path):
    config = CoSpliceConfig('fixture', None, tmp_path / 'samples.tsv', tmp_path / 'out')
    raw = raw_outputs(tmp_path)
    report = config.output_dir / '07_report/katana_cosplice_report.html'
    report.parent.mkdir()
    report.write_text('<html>Old threshold labels are a legacy reference</html>', encoding='utf-8')
    before = signatures([report, *raw.values()])
    result = CoSplicePipeline(config).run_stage('report')
    assert result.suffix == '.html' and result != report
    assert 'label_views' in result.parts
    assert 'Old threshold labels are a legacy reference' not in result.read_text(encoding='utf-8')
    record = json.loads((result.parent / 'NETWORK_LABEL_REPORT_PROVENANCE.json').read_text())
    assert record['network_label_policy']['apply_threshold_labels'] is False
    assert record['network_label_policy']['partition'] == 'source_partition_not_recorded'
    for item in record['current_table_views'].values():
        assert (result.parent / item['file']).is_file()
    assert before == signatures([report, *raw.values()])


def test_report_cannot_relabel_explicitly_refined_raw_cache_as_full_backbone(tmp_path):
    config = CoSpliceConfig('fixture', None, tmp_path / 'samples.tsv', tmp_path / 'out')
    raw = raw_outputs(tmp_path)
    raw['parameters'].write_text(json.dumps({'consensus': {'clustering_method': 'backbone_refined_consensus'}}))
    report = config.output_dir / '07_report/katana_cosplice_report.html'
    report.parent.mkdir()
    report.write_text('<html>Historical reference</html>')
    before = signatures([report, *raw.values()])
    with pytest.raises(ValueError, match='historical refined partition'):
        CoSplicePipeline(config).run_stage('report')
    assert before == signatures([report, *raw.values()])


def test_psi_view_reuses_completed_raw_fit_and_rejects_modified_evidence(tmp_path):
    from katana_next.full_fixture import create_full_fixture
    from katana_next.services.project_service import ProjectService
    from katana_next.workflows.base import ArtifactWriter
    from katana_next.infrastructure.fingerprints import fingerprint_value

    manifest = create_full_fixture(tmp_path / 'project')
    context = ProjectService(application_root=tmp_path / 'application').open_project(manifest)
    workflow = psi_network.PsiNetworkWorkflow()
    parameters = workflow.normalize_parameters({})
    contract = psi_network._psi_fit_contract(context, parameters, [])
    writer = ArtifactWriter(context, 'programs.psi_network', 'raw-psi-fit', parameters, workflow.version,
        common_metadata={'psi_fit_fingerprint': fingerprint_value(contract), 'psi_fit_contract': contract})
    psi_network._raw_psi_table(writer, 'consensus_modules__assignments', 'assignments.parquet', members(), entity_type='splice_event')
    writer.finish()
    raw = writer.directory / 'raw_science_cache/assignments.parquet'
    before = signatures([raw])
    parameters['apply_threshold_labels'] = True
    # Same fit identity under a changed annotation choice.
    assert psi_network._psi_fit_contract(context, parameters, []) == contract
    ids = psi_network._reuse_psi_fit(context, parameters, 'labels-enabled', contract)
    assert ids
    exported = context.catalog.get_artifact(ids[0])
    assert pd.read_parquet(exported.path).membership_role.tolist() == ['core', 'uncertain', 'grey']
    assert exported.metadata['network_label_policy']['apply_threshold_labels'] is True
    assert before == signatures([raw])
    raw.write_bytes(b'changed')
    with pytest.raises(ValueError, match='cache changed'):
        psi_network._reuse_psi_fit(context, parameters, 'must-not-fit', contract)


def test_actual_psi_workflow_switch_exports_a_new_view_without_a_second_fit(tmp_path, monkeypatch):
    from katana_next.full_fixture import create_full_fixture
    from katana_next.services.project_service import ProjectService
    from katana_next.science.frozen_psi import _digest

    manifest = create_full_fixture(tmp_path / 'project')
    content = json.loads(manifest.read_text(encoding='utf-8'))
    sample_ids = [sample['sample_id'] for sample in content['samples']]
    matrix = pd.DataFrame({sample: [.1 + .05*i, .8 - .04*i, .3 + .02*i]
                           for i, sample in enumerate(sample_ids)}, index=['e1', 'e2', 'e3'])
    path = tmp_path / 'psi.tsv'
    matrix.to_csv(path, sep='\t')
    annotations = tmp_path / 'events.tsv'
    pd.DataFrame({'event_key': matrix.index, 'gene_id': ['g1', 'g2', 'g3']}).to_csv(annotations, sep='\t', index=False)
    content['assets'].extend([
        {'asset_id': 'psi', 'kind': 'psi_matrix', 'path': str(path)},
        {'asset_id': 'events', 'kind': 'event_annotation', 'path': str(annotations)}])
    manifest.write_text(json.dumps(content), encoding='utf-8')
    context = ProjectService(application_root=tmp_path / 'application').open_project(manifest)
    calls = []

    def fake_fit(values, samples, **kwargs):
        calls.append('fit')
        scores = pd.DataFrame({'C1': range(len(sample_ids))}, index=pd.Index(sample_ids, name='sample_id'))
        baseline = SimpleNamespace(edges=pd.DataFrame(columns=['source', 'target', 'weight']), program_scores=scores)
        confidence = SimpleNamespace(assignments=members(), module_summary=pd.DataFrame({'consensus_module': ['C1'], 'n_events': [2], 'n_core_events': [1]}),
            bootstrap_runs=pd.DataFrame({'iteration': [1, 2]}), program_scores=scores)
        return values, pd.DataFrame(), values.mean(axis=1), baseline, confidence

    def fake_freeze(*_args, **_kwargs):
        calls.append('freeze')
        model = {'schema': 'katana-frozen-psi-1.0', 'network_status': 'modules_available'}
        model['sha256'] = _digest(model)
        return model

    monkeypatch.setattr(psi_network, 'fit_psi_network', fake_fit)
    monkeypatch.setattr(psi_network, 'freeze_psi_modules', fake_freeze)
    workflow = psi_network.PsiNetworkWorkflow()
    job = SimpleNamespace(report=lambda *_args: None)
    disabled = workflow.normalize_parameters({})
    first = workflow.run(context, job, 'first-fit-off', disabled)
    artifacts = [context.catalog.get_artifact(identity) for identity in first]
    raw = [Path(artifact.metadata['psi_fit_source']['path']) for artifact in artifacts]
    before = signatures(raw)
    enabled = workflow.normalize_parameters({'apply_threshold_labels': True})
    second = workflow.run(context, job, 'same-fit-on', enabled)
    assert calls == ['fit', 'freeze']
    assert signatures(raw) == before
    changed = [context.catalog.get_artifact(identity) for identity in second]
    assignment = next(a for a in changed if a.name == 'consensus_modules__assignments')
    assert pd.read_parquet(assignment.path).membership_role.tolist() == ['core', 'uncertain', 'grey']
    assert all(a.metadata['psi_fit_reuse_source_run'] == 'first-fit-off' for a in changed)


def test_cosplice_workflow_registers_offline_companion_and_attachments(tmp_path, monkeypatch):
    from bs4 import BeautifulSoup
    from katana_next.full_fixture import create_full_fixture
    from katana_next.services.project_service import ProjectService
    from katana_next.science.network_label_report import create_network_label_report

    context = ProjectService(application_root=tmp_path / 'application').open_project(
        create_full_fixture(tmp_path / 'project'))
    raw = tmp_path / 'old-report.html'
    raw.write_text('<html>Frozen reference</html>', encoding='utf-8')
    source = tmp_path / 'current-members.tsv'
    member_label_view(members()).to_csv(source, sep='\t', index=False)
    companion = create_network_label_report(raw, tmp_path / 'companion/report.html',
        policy=network_label_policy(False), table_views={'members': source})

    class CachedEngine:
        def __init__(self, config):
            self.config = config
        def run_stage(self, stage, *, force):
            assert stage == 'report' and not force
            return companion

    monkeypatch.setattr(programs, 'CoSplicePipeline', CachedEngine)
    monkeypatch.setattr(programs, 'INTERPRETATION_STAGES', ('report',))
    workflow = programs.CoSpliceWorkflow()
    ids = workflow.run(context, SimpleNamespace(check_cancelled=lambda: None, report=lambda *_args: None),
        'offline-report', workflow.normalize_parameters({'scope': 'interpretation'}))
    report = next(context.catalog.get_artifact(identity) for identity in ids
                  if context.catalog.get_artifact(identity).kind.value == 'report')
    assert report.path.parent.name == 'report-view'
    assert report.metadata['preserve_relative_report_links'] is True
    assert len(report.metadata['exports']) == 3
    for link in BeautifulSoup(report.path.read_text(encoding='utf-8'), 'html.parser').find_all('a', href=True):
        assert (report.path.parent / link['href']).is_file()
    assert raw.read_text(encoding='utf-8') == '<html>Frozen reference</html>'
