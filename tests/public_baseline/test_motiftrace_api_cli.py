"""The public API and CLI route to established scientific workflows only."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from katana_motiftrace import MotifTrace, SCIENTIFIC_RELEASE, mode_catalog
from katana_motiftrace.cli import main
from katana_next.client import KatanaProject
from katana_next.full_fixture import create_full_fixture
from katana_next.science.bulk import READ_SUPPORT_POLICY_VERSION
from katana_next.workflows.base import ArtifactWriter
from katana_next.workflows.registry import WorkflowRegistry


def test_catalog_is_detached_and_maps_only_existing_frozen_entries():
    catalog = mode_catalog()
    assert catalog['scientific_release'] == SCIENTIFIC_RELEASE == '1.0.0-frozen.20261017'
    registry = WorkflowRegistry()
    assert len(catalog['modes']) == 6
    for item in catalog['modes']:
        assert registry.get(item['workflow_id']).definition.enabled
    catalog['modes'][0]['fixed_parameters']['sequence_scope'] = 'altered'
    assert mode_catalog()['modes'][0]['fixed_parameters'] == {'sequence_scope': 'boundary_association'}


def test_python_modes_delegate_without_new_algorithms_and_borrowed_client_stays_open():
    class ExistingClient:
        def __init__(self):
            self.calls = []
            self.closed = False
        def plan(self, workflow, parameters):
            self.calls.append(('plan', workflow, parameters))
            return {'ready': True}
        def run(self, workflow, parameters, **options):
            self.calls.append(('run', workflow, parameters, options))
            return 'existing-job-record'
        def artifacts(self, workflow):
            return [SimpleNamespace(run_fingerprint='saved', analysis_id=workflow,
                                    parameters={'sequence_scope': scope})
                    for scope in ('boundary_association', 'event_consequences', 'unknown')]
        def close(self):
            self.closed = True

    existing = ExistingClient()
    with MotifTrace(existing) as motiftrace:
        for item in mode_catalog()['modes']:
            name = item['mode_id']
            assert motiftrace.plan(name, selection_id='exact-rows')['ready']
            assert motiftrace.run(name, contrast_ids=['N-vs-D']) == 'existing-job-record'
            assert existing.calls[-1][1] == item['workflow_id']
        assert len(motiftrace.artifacts(run_fingerprint='saved')) == 18
        assert len(motiftrace.artifacts('boundary')) == 1
        assert motiftrace.artifacts('boundary')[0].parameters['sequence_scope'] == 'boundary_association'
        with pytest.raises(ValueError, match='fixes'):
            motiftrace.plan('boundary', {'sequence_scope': 'completed_parallel_models'})
        with pytest.raises(ValueError, match='Conflicting'):
            motiftrace.plan('usage', {'contrast_ids': ['A']}, contrast_ids=['B'])
        with pytest.raises(ValueError, match='Unknown'):
            motiftrace.plan('arbitrary')
        with pytest.raises(ValueError, match='dict'):
            motiftrace.plan('usage', [])
    assert not existing.closed
    with pytest.raises(RuntimeError, match='closed'):
        motiftrace.artifacts()


def test_actual_python_selection_and_cli_plan_preserve_inputs_and_never_submit(tmp_path, capsys):
    path = create_full_fixture(tmp_path / 'project')
    before = {p: hashlib.sha256(p.read_bytes()).hexdigest()
              for p in path.parent.rglob('*') if p.is_file()}
    with KatanaProject(path) as project:
        from katana_next.science.bulk import READ_SUPPORT_POLICY_VERSION
        writer = ArtifactWriter(project.context, 'differential.as', 'measured-toy-source', {}, '3.0.0',
                                common_metadata={'read_support_policy': READ_SUPPORT_POLICY_VERSION})
        rows = pd.DataFrame([dict(row_id='row1', event_id='SE:e1', event_type='SE', chr='chr1',
                                 strand='+', gene_id='g1', comparison_id='perturbation-vs-baseline')])
        source = writer.table('Complete measured AS', 'events.parquet', rows, entity_type='splice_event')
        writer.finish()
        with project.motiftrace as motiftrace:
            selection = motiftrace.select(source.artifact_id, ['row1'])
            assert selection.entity_ids == ('row1',)
            assert selection.source_artifact_id == source.artifact_id
            for invalid in ([], ['missing'], ['row1', 'row1'], 'row1'):
                with pytest.raises(ValueError):
                    motiftrace.select(source.artifact_id, invalid)
            from katana_next.services.selection_routes import compatible_workflows
            assert 'prioritize.motiftrace_consequences' in compatible_workflows(source)
            assert len(motiftrace.inventory()['modes']) == 6
            # Missing original model results must yield readiness errors, not a
            # surrogate fit or a successful empty output.
            plan = motiftrace.plan('parallel', contrast_ids=['perturbation-vs-baseline'])
            assert plan['ready'] is False and plan['errors']
        assert project.jobs.list() == []
        assert project.inventory()
    assert main(['modes']) == 0
    assert len(json.loads(capsys.readouterr().out)['modes']) == 6
    assert main(['--project', str(path), 'plan', 'parallel', '--contrast', 'perturbation-vs-baseline']) == 2
    payload = json.loads(capsys.readouterr().out)
    assert payload['motiftrace']['scientific_release'] == SCIENTIFIC_RELEASE
    assert payload['ready'] is False
    assert main(['--project', str(path), 'plan', 'not-a-mode']) == 1
    assert 'Unknown' in json.loads(capsys.readouterr().err)['error']
    assert {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in before} == before


def test_cli_does_not_execute_or_hide_bad_parameters(tmp_path, capsys):
    path = create_full_fixture(tmp_path / 'project')
    invalid = tmp_path / 'bad.json'
    invalid.write_text('[]', encoding='utf-8')
    assert main(['--project', str(path), 'plan', 'usage', '--parameters', str(invalid)]) == 1
    assert 'JSON object' in json.loads(capsys.readouterr().err)['error']
    with pytest.raises(SystemExit) as missing:
        main(['plan', 'usage'])
    assert missing.value.code == 2


def test_public_api_executes_original_sequence_kernel_then_cli_reuses_cache(tmp_path, capsys):
    """Exercise the new dedicated route, not only its delegated legacy wrapper."""
    path = create_full_fixture(tmp_path / 'project')
    fasta = path.parent / 'inputs/fixture.genome.fa'
    manifest = json.loads(path.read_text(encoding='utf-8'))
    manifest['metadata']['species'] = 'Mus_musculus'
    path.write_text(json.dumps(manifest), encoding='utf-8')
    # Synthetic fixture faidx declares LF offsets; keep Windows text I/O exact.
    fasta.write_text(fasta.read_text(), encoding='ascii', newline='\n')
    length = len(fasta.read_text().splitlines()[1])
    Path(str(fasta) + '.fai').write_text(f'chr1\t{length}\t6\t{length}\t{length + 1}\n')
    motif = tmp_path / 'motif.tsv'
    motif.write_text('motif_id\tRBP\tspecies\tevidence\tsource\tmotif\n'
                     'M\tR\tMus_musculus\tuser_defined\ttest only\tAC\n')
    before = {p: hashlib.sha256(p.read_bytes()).hexdigest()
              for p in (path.parent / 'inputs').rglob('*') if p.is_file()}
    with MotifTrace(path) as mt:
        rows = pd.DataFrame([dict(row_id='event-row', event_id='SE:e1', event_type='SE', chr='chr1',
            strand='+', gene_id='g1', gene_symbol='G1', comparison_id='perturbation-vs-baseline',
            dpsi=.3, fdr=.001, significant=True, stable_significant=True,
            measurement_eligible=True, read_support_policy=READ_SUPPORT_POLICY_VERSION,
            criterion_min_reads=10,
            IJC_SAMPLE_1='60,70', SJC_SAMPLE_1='40,30',
            IJC_SAMPLE_2='30,40', SJC_SAMPLE_2='70,60',
            numerator_mean_reads=100, denominator_mean_reads=100,
            IncLevel1='.6,.7', IncLevel2='.3,.4', upstreamES=0, upstreamEE=40,
            downstreamES=300, downstreamEE=340, exonStart_0base=100, exonEnd=135)])
        writer = ArtifactWriter(mt.project.context, 'test.existing_event_catalog',
                                'measured-toy-event', {}, '1.0')
        source = writer.table('Existing synthetic event', 'events.parquet', rows,
                              entity_type='splicing_event')
        writer.finish()
        selection = mt.select(source.artifact_id, ['event-row'])
        for mode in ('boundary', 'se_impact', 'completed_maps', 'consequences'):
            conflict = mt.plan(mode, selection_id=selection.selection_id,
                               contrast_ids=['rescue-vs-baseline'])
            assert not conflict['ready']
            assert any('different comparisons' in e for e in conflict['errors'])
        with pytest.raises(ValueError, match='different comparisons'):
            mt.run('consequences', selection_id=selection.selection_id,
                   contrast_ids=['rescue-vs-baseline'])
        for invalid in (1, {'wrong': 'type'}, [{}]):
            plan = mt.plan('boundary', selection_id=selection.selection_id, contrast_ids=invalid)
            assert not plan['ready'] and plan['errors']
            with pytest.raises(ValueError):
                mt.run('boundary', selection_id=selection.selection_id, contrast_ids=invalid)
        assert mt.project.jobs.list() == []
        params = {'custom_motif_path': str(motif), 'cluster_permutations': 99}
        plan = mt.plan('consequences', params, selection_id=selection.selection_id)
        assert plan['ready'], plan['errors']
        job = mt.run('consequences', params, selection_id=selection.selection_id, timeout=90)
        assert job.status.value == 'complete', job.error
        outputs = mt.artifacts('consequences', run_fingerprint=job.run_fingerprint)
        rna = next(a for a in outputs if a.path.name == 'ACTUAL_RNA_FORMS.parquet')
        assert set(pd.read_parquet(rna.path).event_type) == {'SE'}
        assert any(a.path.name == 'ACTUAL_LOCAL_RNA.fa' for a in outputs)
        seal = next(a for a in outputs if a.path.name == 'COMPLETE.json')
        complete = json.loads(seal.path.read_text())
        assert complete['complete'] is True
        assert all(a.metadata['inference_status'].startswith('not_estimable') for a in outputs)
        assert all(a.metadata['new_inference'] is False for a in outputs)
        pd.testing.assert_frame_equal(pd.read_parquet(source.path), rows)
    params_file = tmp_path / 'params.json'
    params_file.write_text(json.dumps(params))
    assert main(['--project', str(path), 'run', 'consequences', '--selection-id',
                 selection.selection_id, '--parameters', str(params_file)]) == 0
    assert json.loads(capsys.readouterr().out)['status'] == 'cache_hit'
    assert {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in before} == before
