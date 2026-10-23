"""Completed-family plotting must not turn partial outputs into science results."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from katana_next.reports.motif_views import (
    SEQUENCE_FILES, STAT_FILES, build_completed_motif_views, verify_completed_motif_family,
)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, payload):
    path.write_text(json.dumps(payload, indent=2), encoding='utf-8')


@pytest.fixture
def completed(tmp_path):
    root = tmp_path / 'completed'
    root.mkdir()
    models = {
        '001': {'kind': 'PWM', 'matrix': [[0.7, 0.1, 0.1, 0.1], [0.1, 0.1, 0.7, 0.1]], 'RBPs': 'RbpA'},
        'literal': {'sequence': 'AUG', 'RBPs': 'RbpB'},
        'empty': {'matrix': [[0.25, 0.25, 0.25, 0.25]], 'RBPs': 'RbpC'},
    }
    write_json(root / SEQUENCE_FILES[0], {'comparison_id': 'B_vs_A', 'models': models, 'permutations': 5})
    pd.DataFrame({'event_id': ['E1', 'E2', 'E3'], 'event_type': ['SE', 'RI', 'OTHER']}).to_csv(root / SEQUENCE_FILES[1], sep='\t', index=False)
    pd.DataFrame([
        {'event_id': 'E1', 'event_type': 'SE', 'model_id': '001', 'instance_state': 'inclusion_only', 'consequence_direction': 'increased'},
        {'event_id': 'E2', 'event_type': 'RI', 'model_id': 'literal', 'instance_state': 'retained', 'consequence_direction': 'decreased'},
        {'event_id': 'E3', 'event_type': 'OTHER', 'model_id': '001', 'instance_state': 'exclusion_only', 'consequence_direction': 'unchanged'},
    ]).to_csv(root / SEQUENCE_FILES[2], sep='\t', index=False)
    pd.DataFrame([
        {'event_id': 'E1', 'role': 'foreground', 'layer': 'mature_local_form', 'region': 'inclusion', 'sequence_RNA': 'CCAUGG', 'length': 6},
        {'event_id': 'E1', 'role': 'foreground', 'layer': 'mature_local_form', 'region': 'exclusion', 'sequence_RNA': 'CCGG', 'length': 4},
    ]).to_csv(root / SEQUENCE_FILES[3], sep='\t', index=False)
    pd.DataFrame([{'event_id': 'E1', 'model_id': '001', 'form': 'inclusion', 'RNA_start0': 2, 'RNA_end0': 4, 'hit_RNA': 'AU'}]).to_csv(root / SEQUENCE_FILES[4], sep='\t', index=False)
    pd.DataFrame({'score': [1.25]}).to_csv(root / SEQUENCE_FILES[5], sep='\t', index=False)
    tests = []
    for token, qvalues in [('001', [0.125, np.nan, 0.0625]), ('literal', [0.5, 0.25, 0.75]), ('empty', [np.nan] * 3)]:
        for start, q in zip([0, 7, 19], qvalues):
            tests.append({'event_type': 'SE', 'direction': 'inclusion_up', 'model_id': token,
                          'RBP': models[token]['RBPs'], 'region': 'inclusion', 'method': 'native_method',
                          'window_nt': 5, 'window_start0': start, 'n_matched_gene_pairs': 3,
                          'paired_rank_Z': start / 10, 'mean_paired_score_difference': start / 100,
                          'raw_P_new_paired_swaps': np.nan if pd.isna(q) else q / 2,
                          'q_new_window_family': q, 'status': 'estimable' if pd.notna(q) else 'not_estimable'})
    pd.DataFrame(tests).to_csv(root / STAT_FILES[0], sep='\t', index=False, float_format='%.17g')
    pd.DataFrame({'model_id': ['001'], 'FWER': [0.125]}).to_csv(root / STAT_FILES[1], sep='\t', index=False)
    np.save(root / STAT_FILES[2], [0., 1., 1.5, 2., 2.5])
    members = {name: {'sha256': sha(root / name), 'bytes': (root / name).stat().st_size} for name in (*SEQUENCE_FILES, *STAT_FILES)}
    write_json(root / 'COMPLETE.json', {'complete': True, 'statistical_family_split': False, 'members': members})
    proof = tmp_path / 'local.receipt.json'
    write_json(proof, {'passed': True, 'source_unchanged': True, 'upstream_refits': 0,
                       'summary': {'comparison_id': 'B_vs_A'}, 'complete_manifest': {'SHA256': sha(root / 'COMPLETE.json')}})
    return root, proof


def render(root, proof, out, **kwargs):
    return build_completed_motif_views(root, out, comparison_id='B_vs_A', project_id='study',
                                       source_project_manifest_sha256='a' * 64, completion_proof=proof,
                                       proof_kind='verified_local', **kwargs)


def make_native(root, local):
    marker = root / 'COMPLETE.json'
    binding = {'relative_path': 'COMPLETE.json', 'sha256': sha(marker), 'bytes': marker.stat().st_size, 'members_verified': True}
    proof = root / 'native_completion.receipt.json'
    payload = {'schema': 'katana.native_motiftrace_completion/v1', 'passed': True, 'comparison_id': 'B_vs_A',
               'native_job_id': '123', 'native_job': {'job_id': '123', 'state': 'COMPLETED', 'exit_code': '0:0'},
               'sacct': {'state': 'COMPLETED', 'exit_code': '0:0', 'main_row': {'job_id': '123', 'state': 'COMPLETED', 'exit_code': '0:0'},
                         'actual_command': {'argv': ['sacct', '-X', '-nP', '-j', '123'], 'returncode': 0, 'stdout': '123|COMPLETED|0:0|00:01:00|\n', 'stderr': ''}},
               'nextflow_trace': {'status': 'COMPLETED', 'exit': 0, 'native_id': '123',
                                  'actual_row': {'status': 'COMPLETED', 'exit': '0', 'native_id': '123', 'name': 'MOTIFTRACE (B_vs_A)'}},
               'claim_completed': {'passed': True, 'job_id': '123', 'comparison_id': 'B_vs_A'},
               'scope': 'collection_only', 'collector_RNA_rescans': 0, 'upstream_refits': 0,
               'model_count': 3, 'permutations': 5, 'family_split': False, 'complete_manifest': binding,
               'project_id': 'study', 'study_lineage': {'source_project_manifest_SHA256': 'a' * 64}}
    write_json(proof, payload)
    collection = {'schema': 'katana.motiftrace_collection/v1', 'passed': True, 'comparison_id': 'B_vs_A', 'native_job_id': '123',
                  'all_scientific_members_unchanged': True,
                  'complete_manifest': {**binding, 'sha256_before': binding['sha256'], 'sha256_after': binding['sha256']},
                  'members': json.loads(marker.read_text())['members'],
                  'native_completion_receipt': {'relative_path': proof.name, 'sha256': sha(proof), 'bytes': proof.stat().st_size}}
    write_json(root / 'COLLECTION.receipt.json', collection)
    return proof, payload


def test_complete_family_render_preserves_values_all_models_and_gaps(completed, tmp_path):
    root, proof = completed
    before = {p.name: sha(p) for p in root.iterdir()}
    views = render(root, proof, tmp_path / 'report')
    assert len(views['panels']) == 8
    assert views['audit']['models'] == 3
    assert views['audit']['statistical_test_rows'] == 9
    assert views['audit']['null_values'] == 5
    assert views['audit']['new_scans'] == views['audit']['new_statistical_tests'] == 0
    assert views['audit']['originals_unchanged'] is True
    assert before == {p.name: sha(p) for p in root.iterdir()}
    assert all(p['policy'] is None and p['metadata']['das_scope']['recorded_event_set'] == 'stable_significant' for p in views['panels'])
    for panel in views['panels']:
        assert Path(panel['path']).is_file()
        assert {'plot_data', 'full_results', 'methods'} <= {x['role'] for x in panel['data_links']}
        assert all(Path(x['path']).is_file() for x in panel['data_links'])
        assert all(x['contrast_id'] == 'B_vs_A' and x['policy'] == 'stable_significant' for x in panel['data_links'])
    out = tmp_path / 'report'
    counts = pd.read_csv(out / 'all-models-plot-data.tsv', sep='\t', dtype={'model_id': str})
    assert counts.model_id.tolist() == ['001', 'literal', 'empty']
    assert counts.loc[counts.model_id.eq('empty'), 'inclusion_only'].iloc[0] == 0
    states = pd.read_csv(out / 'instance-states-plot-data.tsv', sep='\t')
    assert states.event_type.tolist() == ['SE', 'RI', 'OTHER']
    profile = pd.read_csv(out / 'window-profile-001-plot-data.tsv', sep='\t', float_precision='round_trip')
    assert profile.window_start0.tolist() == [0, 7, 19]
    assert profile.q_new_window_family.iloc[0] == 0.125
    assert np.isnan(profile.q_new_window_family.iloc[1])
    assert profile.q_new_window_family.iloc[2] == 0.0625
    assert profile.raw_P_new_paired_swaps.iloc[2] == 0.03125
    assert np.isnan(profile.plotted_minus_log10_q.iloc[1])
    assert pd.read_csv(out / 'whole-family-cluster-null-plot-data.tsv', sep='\t')['count'].sum() == 5
    assert all(token in (out / 'ALL_MOTIF_MODELS.html').read_text() for token in ['001', 'literal', 'empty'])


def test_no_eligible_example_is_explicit_not_fabricated(completed, tmp_path):
    root, proof = completed
    views = render(root, proof, tmp_path / 'tiny-display', example_max_nt=1)
    assert len(views['panels']) == 7
    assert any('No eligible' in x['reason'] for x in views['unavailable'])
    assert not (tmp_path / 'tiny-display' / 'target-RNA-context-plot-data.tsv').exists()


@pytest.mark.parametrize('mutation', ['partial', 'split', 'sha', 'missing', 'foreign', 'no_proof'])
def test_invalid_complete_family_rejected_before_writing(completed, tmp_path, mutation):
    root, proof = completed
    marker = root / 'COMPLETE.json'
    if mutation in ('partial', 'split'):
        data = json.loads(marker.read_text()); data['complete' if mutation == 'partial' else 'statistical_family_split'] = mutation == 'split'; write_json(marker, data)
    elif mutation == 'sha':
        (root / SEQUENCE_FILES[1]).write_text('corrupted', encoding='utf-8')
    elif mutation == 'missing':
        (root / STAT_FILES[1]).unlink()
    elif mutation == 'foreign':
        data = json.loads(proof.read_text()); data['summary']['comparison_id'] = 'OTHER'; write_json(proof, data)
    else:
        write_json(proof, {'passed': True})
    with pytest.raises((ValueError, FileNotFoundError)):
        render(root, proof, tmp_path / 'rejected')
    assert not (tmp_path / 'rejected').exists()


def test_existing_output_is_preserved(completed, tmp_path):
    root, proof = completed
    out = tmp_path / 'existing'; out.mkdir(); (out / 'retain').write_text('preserve')
    with pytest.raises(FileExistsError):
        render(root, proof, out)
    assert (out / 'retain').read_text() == 'preserve'


def test_native_completion_requires_all_actual_terminal_proofs(completed):
    root, local = completed
    proof, payload = make_native(root, local)
    evidence = verify_completed_motif_family(root, proof, comparison_id='B_vs_A', native_job_id='123')
    assert evidence['proof_kind'] == 'native'
    assert len(evidence['members']) == len(SEQUENCE_FILES) + len(STAT_FILES)
    with pytest.raises(ValueError, match='Exact native'):
        verify_completed_motif_family(root, proof, comparison_id='B_vs_A', native_job_id='124')


def test_native_lineage_mismatch_rejected_before_any_render(completed, tmp_path):
    root, local = completed
    proof, payload = make_native(root, local)
    with pytest.raises(ValueError, match='source-project manifest lineage'):
        build_completed_motif_views(root, tmp_path/'wrong-lineage', comparison_id='B_vs_A', project_id='study',
                                    source_project_manifest_sha256='b' * 64, completion_proof=proof, native_job_id='123')
    assert not (tmp_path/'wrong-lineage').exists()


@pytest.mark.parametrize('declared', [False, True, None])
def test_native_writer_plan_provenance_binds_unsplit_family(completed, declared):
    root, local = completed
    planpath = root / SEQUENCE_FILES[0]
    plan = json.loads(planpath.read_text())
    plan['provenance'] = {} if declared is None else {'statistical_family_split': declared}
    write_json(planpath, plan)
    marker = root / 'COMPLETE.json'
    complete = json.loads(marker.read_text()); complete.pop('statistical_family_split')
    complete['members'][SEQUENCE_FILES[0]] = {'bytes': planpath.stat().st_size, 'sha256': sha(planpath)}
    write_json(marker, complete)
    proof, _ = make_native(root, local)
    if declared is False:
        assert verify_completed_motif_family(root, proof, comparison_id='B_vs_A', native_job_id='123')['plan']['provenance']['statistical_family_split'] is False
    else:
        with pytest.raises(ValueError, match='unsplit family'):
            verify_completed_motif_family(root, proof, comparison_id='B_vs_A', native_job_id='123')


@pytest.mark.parametrize('fault', ['pending', 'fake_stdout', 'duplicate_stdout', 'trace', 'foreign_trace_name', 'claim', 'ledger', 'refit', 'split_models'])
def test_native_partial_or_mismatched_evidence_cannot_pass(completed, fault):
    root, local = completed
    proof, payload = make_native(root, local)
    if fault == 'pending': payload['native_job']['state'] = 'PENDING'
    elif fault == 'fake_stdout': payload['sacct']['actual_command']['stdout'] = '123|FAILED|1:0|\n'
    elif fault == 'duplicate_stdout': payload['sacct']['actual_command']['stdout'] *= 2
    elif fault == 'trace': payload['nextflow_trace']['native_id'] = '999'
    elif fault == 'foreign_trace_name': payload['nextflow_trace']['actual_row']['name'] = 'MOTIFTRACE (OTHER)'
    elif fault == 'claim': payload['claim_completed']['passed'] = False
    elif fault == 'refit': payload['upstream_refits'] = 1
    elif fault == 'split_models': payload['model_count'] = 1
    if fault == 'ledger':
        ledger = root / 'COLLECTION.receipt.json'; record = json.loads(ledger.read_text()); record['members'] = {}; write_json(ledger, record)
    else:
        write_json(proof, payload)
    with pytest.raises(ValueError):
        verify_completed_motif_family(root, proof, comparison_id='B_vs_A', native_job_id='123')
