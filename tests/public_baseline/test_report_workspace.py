import json
from pathlib import Path
from types import SimpleNamespace
import zipfile
import pandas as pd
import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient

from katana_next.services.report_import import import_report, sha
from katana_next.services.project_service import ProjectService
from katana_next.targets.bundle import export_targets, open_bundle
from katana_next.targets.runner import run_targets, inspect_results
from katana_next.targets.nextflow import prepare_task
from katana_next.api.report_workspace import create_report_router
from target_runner_fixture import create_target_runner_fixture


def report_fixture(tmp_path):
    fixture = create_target_runner_fixture(tmp_path / 'origin')
    root = tmp_path / 'moved-report'
    root.mkdir()
    project = fixture['context'].manifest.model_dump(mode='json')
    project['samples'] = [{k:v for k,v in s.items() if k!='assets'} for s in project['samples']]
    project['annotation'] = {'genome_build':'synthetic-v1'}
    project['metadata']['coordinate_system'] = '0-based-half-open'
    entries, hashes = [], []
    for i, artifact in enumerate([fixture['source_artifact'], *fixture['evidence_artifacts']]):
        target = root / f'table-{i}.parquet'
        target.write_bytes(artifact.path.read_bytes())
        entries.append({'name':artifact.name,'analysis_id':artifact.analysis_id,'source':target.name,
                        'contrast_id':'B_vs_A','parameters':artifact.parameters})
        hashes.append({'report_path':target.name,'source_sha256':sha(target)})
    manifest = {'schema':'katana.approved_report/v1','project_id':project['project_id'],
                'display_project':project,'complete_results':entries,'source_files':hashes}
    (root/'report-manifest.json').write_text(json.dumps(manifest),encoding='utf-8')
    return root, fixture


def test_import_browse_export_cards_and_nextflow(tmp_path):
    root, fixture = report_fixture(tmp_path)
    before = {p.name:sha(p) for p in root.iterdir()}
    projects = ProjectService()
    context = import_report(root, tmp_path/'workspace', projects)
    artifacts = context.catalog.list_artifacts(context.manifest.project_id, kind='table')
    source = next(a for a in artifacts if a.name=='All measured splice events')
    app = FastAPI(); app.include_router(create_report_router(SimpleNamespace(projects=projects)))
    client = TestClient(app)
    url = f'/api/projects/{context.manifest.project_id}/report-tables/{source.artifact_id}'
    result = client.post(url,json={'identifiers':['G2'],'fdr_column':'fdr','fdr_max':1}).json()
    assert result['total']==1
    assert result['rows'][0]['fdr']==.9  # Nonsignificant exact target retained.
    assert client.post(url,json={'comparison':'foreign'}).json()['total']==0
    assert client.post(url,json={'fdr_column':'bad','fdr_max':.1}).status_code==422
    selected = [{'source_artifact_id':source.artifact_id,'row_ids':[result['rows'][0]['row_id']],
                 'target_kind':'event','filters':{'fdr_max':1}}]
    exported = export_targets(context,selected,tmp_path/'request',modules=['gene_card','event_card','sashimi'])
    with open_bundle(exported['archive_path']) as bundle:
        assert bundle.request['targets'][0]['event_id']=='SE:e5'
        assert len(bundle.request['resources'])==9
        assert bundle.request['project']['contrasts'][0]['numerator']['sample_ids']==['B2','B1']
    bindings = fixture['bindings']
    returned = run_targets(exported['request_path'],tmp_path/'returned',resources=bindings,workers=1)
    assert returned['status']=='complete', returned
    assert inspect_results(tmp_path/'returned')['verified']
    prepared = prepare_task(exported['archive_path'],tmp_path/'nf',resources=bindings['resources'])
    assert not prepared['submitted']
    with zipfile.ZipFile(prepared['archive']) as archive:
        assert 'bindings.local.json' not in archive.namelist()
        assert 'katana_targets.nf' in archive.namelist()
    assert before == {p.name:sha(p) for p in root.iterdir()}
    # Immutable imported evidence checked again before target export.
    changed = pd.read_parquet(source.path); changed.loc[0,'fdr']=.99; changed.to_parquet(source.path)
    with pytest.raises(ValueError,match='changed'):
        export_targets(context,selected,tmp_path/'tampered')


def test_tampered_and_unsafe_report_rejected_before_workspace(tmp_path):
    root, _ = report_fixture(tmp_path)
    path = root/'table-0.parquet'; path.write_bytes(path.read_bytes()+b'x')
    with pytest.raises(ValueError,match='checksum'):
        import_report(root,tmp_path/'bad')
    assert not (tmp_path/'bad').exists()
    manifest = json.loads((root/'report-manifest.json').read_text())
    manifest['complete_results'][0]['source']='../outside.parquet'
    (root/'report-manifest.json').write_text(json.dumps(manifest))
    with pytest.raises(ValueError,match='bounded'):
        import_report(root,tmp_path/'bad2')


def test_source_report_cannot_be_destination(tmp_path):
    root, _ = report_fixture(tmp_path)
    with pytest.raises(ValueError,match='separate'):
        import_report(root,root/'workspace')

