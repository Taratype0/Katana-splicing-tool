"""Completed QC is imported, not recalculated, in the planner and report."""
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from test_flow_planner import study
from katana_next.flow.model import write_json
from katana_next.flow.planner import create_plan
from katana_next.flow.native_coverage import profile_selection, stage_native_qc, bind_native_sample_metadata
from katana_next.flow.multiqc import completed_multiqc, reuse_multiqc


def add_profiles(study, count=None):
    payload = json.loads(study.read_text())
    samples = payload['samples'][:count]
    for sample in samples:
        sid = sample['sample_id']
        path = study.parent / (sid + '.coverage.tsv')
        path.write_text('#Transcript position\tTranscript coverage profile\n' +
                        ''.join(f'{i}\t{100+i}\n' for i in range(1, 101)))
        payload['assets'].append({'asset_id': 'qualimap-' + sid, 'kind': 'gene_body_coverage',
                                 'path': str(path), 'metadata': {'sample_id': sid}})
    write_json(study, payload)
    return samples


def test_existing_qualimap_prevents_every_rseqc_scan(study):
    add_profiles(study)
    for requested in (False, True):
        plan = create_plan(study, {'analyses': {'coverage': requested}})
        assert not any(t.get('adapter') in {'qc.genebody', 'reference.bed12'} for t in plan['tasks'])
        summary = next(t for t in plan['tasks'] if t['id'] == 'coverage-summary')
        assert summary['dependencies'] == ['audit']
        assert plan['configuration']['analyses']['coverage']


def test_partial_native_coverage_only_scans_missing_when_requested(study):
    samples = add_profiles(study, 2)
    no_scan = create_plan(study, {})
    assert not any(t.get('adapter') == 'qc.genebody' for t in no_scan['tasks'])
    assert no_scan['configuration']['qc']['coverage_reuse']['unavailable_samples']
    scan = create_plan(study, {'analyses': {'coverage': True}})
    ids = {t['scope']['sample_id'] for t in scan['tasks'] if t.get('adapter') == 'qc.genebody'}
    assert not ids & {s['sample_id'] for s in samples}
    assert len(ids) == len(json.loads(study.read_text())['samples']) - 2
    only = create_plan(study, {'analyses': {'coverage': True}, 'qc': {'coverage_source': 'existing_native'}})
    assert not any(t.get('adapter') == 'qc.genebody' for t in only['tasks'])


def test_prefer_qualimap_and_reject_duplicate_identity(study):
    add_profiles(study, 1)
    from katana_next.infrastructure.manifest_store import ManifestStore
    _, manifest = ManifestStore().load(study)
    selected, _ = profile_selection(manifest)
    sid, aid = next(iter(selected.items()))
    assert aid == 'qualimap-' + sid
    payload = json.loads(study.read_text())
    other = dict(payload['assets'][-1], asset_id='duplicate')
    payload['assets'].append(other)
    write_json(study, payload)
    with pytest.raises(ValueError, match='Multiple coverage'):
        create_plan(study, {})


def native_multiqc(tmp_path, samples=('s1', 's2')):
    html = tmp_path / 'multiqc_report.html'
    html.write_text('<html><title>MultiQC</title><body>Native quality plots</body></html>')
    data = tmp_path / 'multiqc_data.json'
    write_json(data, {'report_saved_raw_data': {'multiqc_fastqc': {sid+'_R1': {'total_sequences': 100} for sid in samples}}})
    return SimpleNamespace(samples=[SimpleNamespace(sample_id=s) for s in samples], assets=[
        SimpleNamespace(path=p, kind='sequencing_qc', metadata={'qc_type': 'multiqc'}) for p in (html, data)])


def test_multiqc_reuses_bytes_and_never_reads_bams(tmp_path):
    manifest = native_multiqc(tmp_path)
    html, data, proof = completed_multiqc(manifest)
    before = {p: p.read_bytes() for p in (html.path, data.path)}
    result = reuse_multiqc(manifest, tmp_path/'qc')
    assert result['evidence']['upstream_jobs'] == 0
    assert not result['evidence']['raw_reads_or_bams_read']
    assert (tmp_path/'qc/multiqc_report.html').read_bytes() == before[html.path]
    assert before == {p: p.read_bytes() for p in before}
    assert len(result['outputs']) == 4


def test_multiqc_missing_samples_and_version_only_not_passed(tmp_path):
    manifest = native_multiqc(tmp_path)
    manifest.samples.append(SimpleNamespace(sample_id='s3'))
    with pytest.raises(ValueError, match='missing declared samples'):
        completed_multiqc(manifest)
    write_json(tmp_path/'multiqc_data.json', {'report_saved_raw_data': {'software_versions': {'s1': {'v': 1}}}})
    with pytest.raises(ValueError, match='no measured QC'):
        completed_multiqc(manifest)


def test_native_directory_imports_only_qc_with_exact_sample_mapping(tmp_path):
    source = tmp_path/'native'
    directory = source/'s1/raw_data_qualimapReport'
    directory.mkdir(parents=True)
    (directory/'coverage_profile_along_genes_(total).txt').write_text('coverage')
    (source/'huge.bam').write_text('must not copy')
    assets = stage_native_qc(source, tmp_path/'configuration/native-qc')
    assert len(assets) == 1
    bind_native_sample_metadata(assets[0], {'s1'})
    assert assets[0]['metadata']['sample_id'] == 's1'
    assert not list((tmp_path/'configuration').rglob('*.bam'))
    with pytest.raises(ValueError, match='Unknown coverage sample'):
        bind_native_sample_metadata(assets[0], {'s2'})


def test_installed_coverage_workflow_reuses_native_curves(study, tmp_path):
    add_profiles(study)
    from katana_next.flow.task import execute_task
    plan = create_plan(study, {})
    path = write_json(tmp_path/'plan.json', plan)
    execute_task(path, 'audit', tmp_path/'audit')
    execute_task(path, 'coverage-summary', tmp_path/'coverage', dependency_receipts=[tmp_path/'audit/task.receipt.json'])
    receipt = json.loads((tmp_path/'coverage/task.receipt.json').read_text())
    assert receipt['status'] == 'completed'
    assert any(a['name'] == 'Gene-body coverage curves' for a in receipt['generated_artifacts'])


def test_multiqc_only_curves_prevent_scans_and_execute_coverage(study, tmp_path):
    payload = json.loads(study.read_text())
    curves = {s['sample_id']: {str(i): 100+i for i in range(1,101)} for s in payload['samples']}
    path = write_json(study.parent/'multiqc_data.json', {'report_saved_raw_data': {'qualimap_rnaseq_cov_hist': curves}})
    payload['assets'].append({'asset_id':'multiqc-data','kind':'sequencing_qc','path':str(path),
                              'metadata':{'qc_type':'multiqc'}})
    write_json(study, payload)
    plan = create_plan(study, {'analyses':{'coverage':True}})
    assert not any(t.get('adapter') == 'qc.genebody' for t in plan['tasks'])
    assert set(plan['configuration']['qc']['coverage_reuse']['selected_asset_by_sample'].values()) == {'multiqc-data'}
    from katana_next.flow.task import execute_task
    saved = write_json(tmp_path/'plan.json', plan)
    execute_task(saved, 'audit', tmp_path/'audit')
    execute_task(saved, 'coverage-summary', tmp_path/'coverage', dependency_receipts=[tmp_path/'audit/task.receipt.json'])
    receipt = json.loads((tmp_path/'coverage/task.receipt.json').read_text())
    assert receipt['status'] == 'completed'
    assert any(a['name'] == 'Gene-body coverage curves' for a in receipt['generated_artifacts'])
    add_profiles(study, 1)
    preferred = create_plan(study, {})['configuration']['qc']['coverage_reuse']['selected_asset_by_sample']
    assert sum(v.startswith('qualimap-') for v in preferred.values()) == 1


def test_multiqc_curve_does_not_accept_renamed_samples(study):
    payload = json.loads(study.read_text())
    path = write_json(study.parent/'multiqc_data.json', {'report_saved_raw_data':{
        'qualimap_rnaseq_cov_hist':{'unrelated':{str(i):1 for i in range(1,101)}}}})
    payload['assets'].append({'asset_id':'multiqc-data','kind':'sequencing_qc','path':str(path)})
    write_json(study, payload)
    with pytest.raises(ValueError, match='exact sample identity'):
        create_plan(study, {})


def test_multiqc_frame_uses_fixed_chapter_one(tmp_path):
    from katana_next.reports.engine import Builder
    from katana_next.reports.layout import render
    from katana_next.full_fixture import create_full_fixture
    from katana_next.infrastructure.manifest_store import ManifestStore
    project = create_full_fixture(tmp_path/'fixture')
    payload = json.loads(project.read_text())
    source = {'project': payload, 'artifacts': [], 'report_analyses': {}, 'mode': 'bulk'}
    output = tmp_path/'report'
    output.mkdir()
    b = Builder(source, output, 'stable_significant')
    html = tmp_path/'multiqc_report.html'
    html.write_text('<html>MultiQC native plots</html>')
    b.add_views({'panels': [{'id': 'multiqc', 'chapter': 1, 'title': 'MultiQC sequencing quality',
                           'path': str(html), 'metadata': {'native_multiqc': True}}]})
    # Rendering the report must use an HTML browser, never an img tag for HTML.
    render(b)
    text = (output/'REPORT.html').read_text()
    assert 'title="MultiQC sequencing quality"' in text
    assert '<img loading="lazy" data-zoom="'+b.panels[0]['image'] not in text
    assert 'chapter-number">01' in text
    assert b.panels[0]['metadata']['native_multiqc']
    assert (output/b.panels[0]['image']).read_bytes() == html.read_bytes()
