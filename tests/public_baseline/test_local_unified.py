"""Full-population exploration and local scientific provenance, independent of UI pages."""
import json
from pathlib import Path
from types import SimpleNamespace
import zipfile

import numpy as np
import pandas as pd
import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient
from scipy.stats import hypergeom

from katana_next.api.report_workspace import create_report_router
from katana_next.services.project_service import ProjectService
from katana_next.services.report_import import import_report, sha
from katana_next.services.report_archive import report_source
from test_report_workspace import report_fixture


def workspace(tmp_path, count=321, mode='bulk'):
    root, fixture = report_fixture(tmp_path)
    manifest = json.loads((root / 'report-manifest.json').read_text())
    manifest['mode'] = mode
    frame = pd.DataFrame({'gene_id': [f'G{i}' for i in range(count)],
        'event_id': [f'SE:{i}' for i in range(count)], 'event_type': ['SE'] * count,
        'comparison_id': ['B_vs_A'] * count, 'cell_type': ['T'] * count,
        'log2FoldChange': np.linspace(-2, 2, count), 'pvalue': [.01] * count,
        'padj': [.01 if i < 21 else .9 for i in range(count)]})
    frame.loc[count-1, 'padj'] = np.nan
    frame.to_parquet(root / 'full.parquet', index=False)
    manifest['complete_results'] = [{'source':'full.parquet','name':'Complete DEG / DAS',
        'analysis_id':'differential.deg','contrast_id':'B_vs_A','metadata':{'full_results':True}}]
    manifest['source_files'] = [{'report_path':'full.parquet','source_sha256':sha(root/'full.parquet')}]
    (root / 'report-manifest.json').write_text(json.dumps(manifest),encoding='utf-8')
    projects = ProjectService()
    context = import_report(root, tmp_path / 'workspace', projects)
    artifact = context.catalog.list_artifacts(context.manifest.project_id)[0]
    runtime = SimpleNamespace(projects=projects, jobs=SimpleNamespace(bind_project=lambda *args:None))
    app = FastAPI();app.include_router(create_report_router(runtime))
    client = TestClient(app)
    url = f'/api/projects/{context.manifest.project_id}/report-tables/{artifact.artifact_id}'
    return client, url, context, artifact, root


def test_full_population_filters_targets_export_and_na(tmp_path):
    client, url, context, artifact, root = workspace(tmp_path)
    filters = {'effect_column':'log2FoldChange','fdr_column':'padj'}
    page = client.post(url,json={**filters,'limit':100}).json()
    assert page['total'] == 321 and len(page['rows']) == 100
    plotted = client.post(url+'/plot',json=filters).json()
    assert plotted['plotted'] == 320 and plotted['unplottable'] == 1
    assert len(plotted['points']) == 320  # Never the first page only.
    all_rows = client.post(url+'/selection',json=filters).json()
    assert all_rows['total'] == 321
    filtered = {**filters,'fdr_max':.05,'effect_min':1.5,'direction':'negative'}
    assert client.post(url,json=filtered).json()['total'] == 21
    assert client.post(url+'/plot',json=filtered).json()['plotted'] == 21
    assert client.post(url+'/selection',json=filtered).json()['total'] == 21
    exported = client.post(url+'/export',json=filtered)
    assert exported.status_code == 200, exported.text
    assert len(exported.text.strip().splitlines()) == 22
    assert client.post(url,json={**filters,'scopes':{'cell_type':'B'}}).json()['total'] == 0
    assert client.post(url,json={**filters,'effect_low':1,'effect_high':2}).json()['total'] == 81
    assert sha(artifact.path) == artifact.metadata['imported_table_sha256']


def test_density_accounts_for_every_valid_row(tmp_path):
    client, url, *_ = workspace(tmp_path,count=150002)
    response=client.post(url+'/plot',json={'effect_column':'log2FoldChange','fdr_column':'padj'})
    assert response.status_code==200,response.text
    result=response.json()
    assert result['mode']=='density'
    assert sum(n for x,y,n in result['bins']) == result['plotted'] == 150001
    assert result['unplottable']==1


def test_large_text_table_mixed_empty_numeric_chunks_and_incomplete_import(tmp_path):
    root, _ = report_fixture(tmp_path)
    manifest = json.loads((root/'report-manifest.json').read_text())
    path = root/'mixed.tsv'
    # Cross pandas' low-memory inference boundary: an empty token in one chunk,
    # numeric data in another previously produced an unwriteable object column.
    path.write_text('gene_id\tmean_paired_score_difference\n001\t\n' + '002\t0.25\n'*600000, encoding='utf-8')
    manifest['complete_results']=[{'source':'mixed.tsv','name':'Complete mixed table'}]
    manifest['source_files']=[]
    (root/'report-manifest.json').write_text(json.dumps(manifest),encoding='utf-8')
    projects=ProjectService()
    context=import_report(root,tmp_path/'mixed-workspace',projects)
    artifact=context.catalog.list_artifacts(context.manifest.project_id)[0]
    frame=pd.read_parquet(artifact.path)
    assert len(frame)==600001 and frame.gene_id.iloc[0]=='001'
    assert frame.mean_paired_score_difference.iloc[0]==''
    assert float(frame.mean_paired_score_difference.iloc[-1])==.25
    (context.manifest_path.parent/'report-import.json').unlink()
    app=FastAPI();app.include_router(create_report_router(SimpleNamespace(projects=projects)))
    response=TestClient(app).post('/api/projects/import-report',json={'report_path':str(root)})
    assert response.status_code==422 and 'did not complete' in response.text


def test_enrichment_uses_full_background_and_all_terms(tmp_path):
    client,url,context,artifact,root=workspace(tmp_path,count=321,mode='single-cell')
    library=tmp_path/'pathways.gmt'
    library.write_text('P1\tdesc\t'+'\t'.join(f'G{i}' for i in range(21))+'\nP2\tdesc\t'+'\t'.join(f'G{i}' for i in range(21,321))+'\n',encoding='utf-8')
    request={'filters':{'comparison':'B_vs_A','scopes':{'cell_type':'T'},'fdr_column':'padj','fdr_max':.05},
             'library_path':str(library),'gene_column':'gene_id','test_column':'pvalue','library_species':'other'}
    result=client.post(url+'/enrichment',json=request)
    assert result.status_code==200,result.text
    result=result.json()
    assert result['audit']['assay_mode']=='single-cell'
    assert result['audit']['selected_unique_genes']==21
    assert result['audit']['annotated_background_genes']==321
    assert result['audit']['tested_terms']==2
    first=result['rows'][0]
    assert first['pvalue']==pytest.approx(hypergeom.sf(20,321,21,21),rel=1e-10,abs=1e-100)
    assert first['fdr']==pytest.approx(first['pvalue']*2,rel=1e-10,abs=1e-100)
    assert result['rows'][1]['overlap_genes']==0 and result['rows'][1]['pvalue']==1
    assert client.get(result['files']['background.tsv']).status_code==200
    changed=pd.read_parquet(artifact.path);changed.loc[0,'padj']=1;changed.to_parquet(artifact.path)
    assert client.post(url+'/enrichment',json=request).status_code==422


def test_zip_and_directory_import_are_equivalent_and_unsafe_zip_rejected(tmp_path):
    _,_,_,_,root=workspace(tmp_path)
    original={p.name:sha(p) for p in root.iterdir() if p.is_file()}
    zipped=tmp_path/'REPORT.zip'
    with zipfile.ZipFile(zipped,'w',zipfile.ZIP_DEFLATED) as z:
        for path in root.iterdir():
            if path.is_file():z.write(path,'study/'+path.name)
    extracted=report_source(zipped,tmp_path/'cache')
    assert report_source(extracted,tmp_path/'cache')==extracted
    assert report_source(zipped,tmp_path/'cache')==extracted
    assert original=={p.name:sha(p) for p in extracted.iterdir() if p.is_file()}
    imported=import_report(extracted,tmp_path/'zip-workspace')
    assert imported.catalog.list_artifacts(imported.manifest.project_id)[0].row_count==321
    bad=tmp_path/'bad.zip'
    with zipfile.ZipFile(bad,'w') as z:
        z.writestr('../outside.txt','no');z.writestr('report-manifest.json',(root/'report-manifest.json').read_bytes())
    with pytest.raises(ValueError,match='Unsafe'):
        report_source(bad,tmp_path/'cache')
    assert not (tmp_path/'outside.txt').exists()
    (extracted/'full.parquet').write_bytes(b'changed')
    with pytest.raises(ValueError,match='changed'):
        report_source(zipped,tmp_path/'cache')
    assert original=={p.name:sha(p) for p in root.iterdir() if p.is_file()}


def test_single_cell_and_unknown_metadata_are_not_silently_bulk(tmp_path):
    from katana_next.services.assay import assay_identity
    assert assay_identity({})['mode']=='unknown'
    assert assay_identity({'assay':'snRNA-seq'})['mode']=='single-cell'
    assert assay_identity({'assay':'bulk_rna_seq'})['mode']=='bulk'
    assert assay_identity({'assay_mode':'single_cell'})['mode']=='single-cell'
    assert assay_identity({},[{'kind':'single_cell_metadata'}])['mode']=='single-cell'
    _,_,context,_,_=workspace(tmp_path,mode='single-cell')
    assert context.manifest.metadata['assay_mode']=='single-cell'
    assert 'donor' in context.manifest.metadata['assay_identity']['replicate_unit']


def test_storage_only_copy_preserves_native_bytes_and_nonsignificant_rows(tmp_path):
    from katana_next.services.report_archive import compact_copy
    _,_,_,_,root=workspace(tmp_path)
    original=(root/'full.parquet').read_bytes()
    (root/'duplicate.parquet').write_bytes(original)
    (root/'REPORT.html').write_text('<html><h1>Original presentation</h1><a href="duplicate.parquet">Complete results</a></html>',encoding='utf-8')
    before={p.name:sha(p) for p in root.iterdir() if p.is_file()}
    audit=compact_copy(root,tmp_path/'compact')
    assert len(audit['deleted_duplicate_aliases'])==1
    assert audit['native_statistics_recomputed'] is False
    manifest=json.loads((tmp_path/'compact/report-manifest.json').read_text())
    source=tmp_path/'compact'/manifest['complete_results'][0]['source']
    assert source.read_bytes()==original
    frame=pd.read_parquet(source)
    assert len(frame)==321 and (frame.padj==.9).sum()==299 and frame.padj.isna().sum()==1
    assert before=={p.name:sha(p) for p in root.iterdir() if p.is_file()}
    with zipfile.ZipFile(tmp_path/'compact/REPORT.zip') as z:
        assert z.read(source.relative_to(tmp_path/'compact').as_posix())==original


def test_empty_report_evidence_target_mask_is_boolean():
    import pyarrow as pa
    import pyarrow.compute as pc
    from katana_next.targets.bundle import _target_mask
    mask = _target_mask(pa.table({"summary": pa.array([], type=pa.string())}), {"target_kind":"gene", "gene_id":"G1"}, set())
    assert mask.type == pa.bool_()
    assert pc.any(mask).as_py() is None


def test_pathway_context_is_not_fabricated_event_evidence():
    from katana_next.targets.runner import _event_cards
    target={"target_kind":"event", "event_id":"SE:exact", "event_type":"SE", "gene_id":"G1", "comparison_ids":["B_vs_A"], "event_geometry":{"chr":"chr1","strand":"+"}}
    context=pd.DataFrame({"term":["GO:1"]*1000,"pvalue":[.01]*1000})
    cards=_event_cards({},target,[({"evidence_id":"context","metadata":{}},context)])
    assert len(cards)==1 and cards[0]["event_id"]=="SE:exact"
    assert cards[0]["evidence"]==[]


def test_repeated_string_redaction_preserves_rows_nulls_and_counts():
    import pyarrow as pa
    from katana_next.targets.bundle import _sanitize_table
    values=['C:/private/reads.bam', 'G1,G2,G3', None]*10000
    table=pa.table({'plain':pa.chunked_array([values[:12000],values[12000:]]),
                    'dictionary':pa.array(values).dictionary_encode()})
    result,counts=_sanitize_table(table,{})
    assert counts=={'plain':10000,'dictionary':10000}
    assert result.num_rows==30000 and result['plain'].null_count==10000
    assert result['plain'][1].as_py()=='G1,G2,G3'
    assert result['plain'][0].as_py().startswith('local-reference:')
    assert result['dictionary'].to_pylist()==result['plain'].to_pylist()
    assert not _sanitize_table(result,{})[1]


def test_empty_tsv_geometry_does_not_duplicate_exact_event_cards():
    from katana_next.targets.runner import _event_cards
    target={'target_kind':'event','event_id':'RI:1','event_type':'RI','gene_id':'G1',
            'comparison_ids':['B_vs_A'],'event_geometry':{'chr':'chr1','strand':'+','riExonStart_0base':10,'riExonEnd':30}}
    frame=pd.DataFrame({'event_id':['RI:1','RI:1'],'event_type':['RI','RI'],
                        'exonStart_0base':[None,''],'FDR':[.01,.01]})
    evidence={'evidence_id':'saved','metadata':{},'sha256':'source'}
    cards=_event_cards({},target,[(evidence,frame)])
    assert len(cards)==1 and len(cards[0]['evidence'])==2
    assert 'exonStart_0base' not in cards[0]['event_geometry']


def test_large_card_context_preserves_all_rows_without_duplicate_text(tmp_path):
    from katana_next.targets.runner import _write_evidence_table
    (tmp_path/'tables').mkdir()
    frame=pd.DataFrame({'gene_id':[f'G{i}' for i in range(6001)],'FDR':[.9]*6000+[None]})
    files=_write_evidence_table(tmp_path,'tables/context',{'evidence_id':'context'},frame)
    pd.testing.assert_frame_equal(pd.read_parquet(tmp_path/files['parquet']),frame)
    assert files['preview_rows']==200 and 'tsv' not in files and 'json' not in files
    html=(tmp_path/files['html']).read_text('utf-8')
    assert '200 of 6,001 rows' in html and 'Download all 6,001 rows' in html
    assert 'Search preview' in html and 'complete rows' not in html


def test_report_cache_uses_configured_drive_and_checks_free_space(tmp_path, monkeypatch):
    from katana_next.api.local_files import workspace_home
    from katana_next.services import report_archive
    monkeypatch.setenv('KATANA_WORKSPACE_HOME',str(tmp_path/'larger-drive'))
    assert workspace_home()==(tmp_path/'larger-drive').resolve()
    source,_=report_fixture(tmp_path)
    archive=tmp_path/'REPORT.zip'
    with zipfile.ZipFile(archive,'w') as z:
        for path in source.iterdir():z.write(path,path.name)
    monkeypatch.setattr(report_archive.shutil,'disk_usage',lambda path:SimpleNamespace(free=0))
    with pytest.raises(ValueError,match='Not enough space'):
        report_source(archive,workspace_home())
    assert not list((workspace_home()/'report-archives').glob('extract-*'))
