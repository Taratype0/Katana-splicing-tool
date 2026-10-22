from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import pytest

from katana_next.domain.models import AnalysisArtifact, ArtifactKind
from katana_next.science.rbp_evidence_layers import (
    CONTRACT, LAYERS, evidence_layer, layer_inventory, overview_html, report_groups,
)


def item(workflow, scope=None, contrast='PBS', name='Saved result'):
    return AnalysisArtifact(project_id='test', analysis_id=workflow, name=name,
        kind=ArtifactKind.FILE, entity_type='saved_evidence', path=Path('not_read.tsv'),
        fingerprint='same', run_fingerprint='same', contrast_ids=(contrast,),
        parameters={} if scope is None else {'sequence_scope':scope},
        metadata={'original_P':1e-24, 'original_FDR':.03})


@pytest.mark.parametrize('workflow,scope,expected', [
    ('prioritize.rbp_support', None, 'identity'),
    ('prioritize.event_motifs', None, 'regulation'),
    ('prioritize.binding_evidence', None, 'regulation'),
    ('prioritize.binding_enrichment', None, 'regulation'),
    ('prioritize.rbp_mechanism', None, 'regulation'),
    ('prioritize.rbp_mechanism', 'boundary_association', 'regulation'),
    ('prioritize.rbp_mechanism', 'event_consequences', 'consequence'),
    ('prioritize.rbp_mechanism', 'completed_parallel_models', 'consequence'),
    ('prioritize.rbp_sequence_results', None, 'consequence'),
    ('prioritize.motiftrace_usage', None, 'consequence'),
])
def test_layer_matches_workflow_and_actual_mode(workflow, scope, expected):
    assert evidence_layer(workflow, {} if scope is None else {'sequence_scope':scope})['id'] == expected


def test_unknown_ids_and_modes_are_not_guessed():
    assert evidence_layer('programs.cosplice') is None
    assert evidence_layer('prioritize.rbp_mechanism', {'sequence_scope':'unknown_future_mode'}) is None


def test_same_entry_different_modes_and_contrasts_remain_separate():
    rows = [item('prioritize.rbp_mechanism'),
            item('prioritize.rbp_mechanism', 'event_consequences'),
            item('prioritize.rbp_mechanism', 'event_consequences', 'KA')]
    before = deepcopy([a.model_dump() for a in rows])
    groups = report_groups(rows)
    assert len(groups) == 3
    assert [g[1]['id'] for g in groups] == ['regulation', 'consequence', 'consequence']
    assert len(layer_inventory(rows)[0]['artifact_ids']) == 0
    assert [a.model_dump() for a in rows] == before


def test_missing_evidence_is_not_a_negative_and_only_english_is_supported():
    en = overview_html([])
    assert 'not a negative result' in en
    with pytest.raises(ValueError, match='Only English'):
        overview_html([], 'ZH')
    for layer in LAYERS:
        assert layer['title_en'] in en
    assert 'Switching views does not run an analysis' in en


def test_existing_html_report_uses_three_layers_without_reading_source_data():
    from katana_next.workflows.reporting import _html_report
    manifest = SimpleNamespace(metadata={}, name='Fixture <study>', assets=[],
        annotation=SimpleNamespace(genome_build='hg19', gtf=None), samples=[], contrasts=[])
    rows = [item('prioritize.rbp_support'), item('prioritize.rbp_mechanism'),
            item('prioritize.rbp_mechanism','completed_sequence_maps')]
    before = [a.model_dump() for a in rows]
    result = _html_report(SimpleNamespace(manifest=manifest), rows)
    assert result.count('data-rbp-layer="identity"') == 2
    assert result.count('data-rbp-layer="regulation"') == 2
    assert result.count('data-rbp-layer="consequence"') == 2
    assert 'Fixture &lt;study&gt;' in result
    assert 'not a mature-RNA sequence removed by DAS' in result
    assert [a.model_dump() for a in rows] == before


def test_presentation_keeps_original_rbp_classes_ids_defaults_and_versions():
    from katana_next.workflows.registry import WorkflowRegistry
    registry = WorkflowRegistry()
    for layer in LAYERS:
        workflow = registry.get(layer['workflow_id'])
        before = workflow.normalize_parameters({})
        data = workflow.definition.to_dict()
        assert data['workflow_id'] == layer['workflow_id']
        assert data['presentation_role'] == 'primary'
        assert workflow.normalize_parameters({}) == before
    assert registry.get('prioritize.rbp_support').definition.presentation_role == 'auxiliary'
    assert registry.get('prioritize.rbp_mechanism').version == '1.6.1'
    assert registry.get('prioritize.motiftrace_consequences').version == '1.0.1'


def test_frontend_imports_the_same_contract():
    repo = Path(__file__).resolve().parents[2]
    text = (repo/'apps/web/src/rbpEvidence.ts').read_text(encoding='utf-8')
    assert '../../../katana_next/science/rbp_evidence_layers.json' in text
    assert [l['number'] for l in CONTRACT['layers']] == [1, 2, 3]


def test_same_layer_different_modes_are_not_grouped_together():
    artifacts = [item('prioritize.rbp_mechanism', 'completed_sequence_maps'),
                 item('prioritize.rbp_mechanism', 'completed_parallel_models')]
    groups = report_groups(artifacts)
    assert len(groups) == 2
    assert [group[1]['id'] for group in groups] == ['consequence', 'consequence']


def test_latest_only_keeps_original_rbp_modes_and_their_own_P_values(tmp_path):
    from katana_next.workflows.reporting import ReproducibilityReportWorkflow
    rows = [item('prioritize.rbp_mechanism', mode) for mode in
            ('completed_sequence_maps', 'completed_parallel_models')]
    for i, artifact in enumerate(rows):
        path = tmp_path / f'{i}.tsv'
        path.write_text('original_P\n0.001\n')
        rows[i] = artifact.model_copy(update={'path':path, 'metadata':{'project_snapshot_fingerprint':'same'}})
    context = SimpleNamespace(manifest=SimpleNamespace(project_id='test'), snapshot=SimpleNamespace(manifest_fingerprint='same'),
                              catalog=SimpleNamespace(list_artifacts=lambda _:rows))
    from unittest.mock import patch
    with patch('katana_next.workflows.sources.source_complete', return_value=True):
        selected = ReproducibilityReportWorkflow().parent_artifacts(context, {'latest_only':True})
    assert len(selected) == 2
    assert [a.parameters['sequence_scope'] for a in selected] == ['completed_sequence_maps', 'completed_parallel_models']
