import numpy as np
import pandas as pd
import pytest
from katana_next.science.single_cell_imports import import_marvel, import_longread, run_import
from katana_next.workflows.registry import WorkflowRegistry


def fixture(tmp_path):
    meta=pd.DataFrame(dict(cell_id=['c1','c2'],donor_id=['d1','d2'],cell_type=['neuron']*2,condition=['WT','MUT']))
    meta.to_csv(tmp_path/'meta.tsv',sep='\t',index=False)
    table=pd.DataFrame({'tran_id':['e1','e2'],'mean.g1':[20,50],'mean.g2':[70,50],'mean.diff':[50,0],'p.val':[.001,.8],'p.val.adj':[.002,.8]})
    table.to_csv(tmp_path/'results.tsv',sep='\t',index=False)
    pd.DataFrame({'tran_id':['e1','e2'],'c1':[.2,np.nan],'c2':[.7,.5]}).to_csv(tmp_path/'psi.tsv',sep='\t',index=False)
    return dict(mode='marvel',metadata_path=str(tmp_path/'meta.tsv'),results_path=str(tmp_path/'results.tsv'),psi_path=str(tmp_path/'psi.tsv'),method='wilcox',tool_version='declared-test-version',reference='WT',target='MUT',cell_type='neuron',adjustment='BH')


def test_marvel_native_units_direction_missing_and_no_retest(tmp_path):
    c=fixture(tmp_path);tables=import_marvel(c);df=tables['marvel_complete_results']
    assert df.delta_psi.tolist()==[.5,0]
    assert df.source_fdr.tolist()==[.002,.8]
    assert tables['marvel_cell_psi'].psi.isna().sum()==1
    assert len(df)==2


def test_permutation_adjusted_name_is_not_fdr(tmp_path):
    c=fixture(tmp_path);c['method']='permutation'
    with pytest.raises(ValueError,match='not assumed BH'):import_marvel(c)
    c['adjustment']='none'
    assert import_marvel(c)['marvel_complete_results'].source_fdr.isna().all()


def test_wrong_units_rejected(tmp_path):
    c=fixture(tmp_path);f=pd.read_csv(c['psi_path'],sep='\t');f.loc[0,'c1']=20;f.to_csv(c['psi_path'],sep='\t',index=False)
    with pytest.raises(ValueError,match='native PSI'):import_marvel(c)


def test_reversed_effect_rejected(tmp_path):
    c=fixture(tmp_path);f=pd.read_csv(c['results_path'],sep='\t');f.loc[0,'mean.diff']=-50;f.to_csv(c['results_path'],sep='\t',index=False)
    with pytest.raises(ValueError,match='group2 minus group1'):import_marvel(c)


def test_unknown_cell_rejected(tmp_path):
    c=fixture(tmp_path);f=pd.read_csv(c['psi_path'],sep='\t').rename(columns={'c1':'unknown'});f.to_csv(c['psi_path'],sep='\t',index=False)
    with pytest.raises(ValueError,match='absent'):import_marvel(c)


def long_fixture(tmp_path):
    c=fixture(tmp_path)
    counts=pd.DataFrame({'cell_id':['c1','c1','c2','c2'],'gene_id':['g']*4,'transcript_id':['t1','t2']*2,'count':[20,5,5,20]})
    counts.to_csv(tmp_path/'counts.tsv',sep='\t',index=False)
    pd.DataFrame({'isoform':['t1','t2'],'associated_gene':['g','g'],'structural_category':['full-splice_match','novel_in_catalog'],'predicted_NMD':['non_NMD','NMD']}).to_csv(tmp_path/'anno.tsv',sep='\t',index=False)
    c.update(mode='longread_isoforms',counts_path=str(tmp_path/'counts.tsv'),annotation_path=str(tmp_path/'anno.tsv'),reference_build='test-build')
    c.pop('results_path');c.pop('psi_path')
    return c


def test_longread_annotation_counts_and_low_donor_gate(tmp_path):
    tables=import_longread(long_fixture(tmp_path))
    assert tables['structural_landscape'].molecules.sum()==50
    assert tables['isoform_annotations'].predicted_NMD.tolist()==['non_NMD','NMD']
    assert tables['gene_DIU'].empty or tables['gene_DIU'].fdr.isna().all()


def test_longread_gene_identity_rejected(tmp_path):
    c=long_fixture(tmp_path);f=pd.read_csv(c['annotation_path'],sep='\t');f.loc[0,'associated_gene']='other';f.to_csv(c['annotation_path'],sep='\t',index=False)
    with pytest.raises(ValueError,match='gene mismatch'):import_longread(c)


def test_english_report_end_to_end(tmp_path):
    from bs4 import BeautifulSoup
    out=tmp_path/'output';run_import(fixture(tmp_path),out)
    for lang in ('EN',):
        soup=BeautifulSoup((out/f'REPORT_{lang}.html').read_text(encoding='utf-8'),'html.parser')
        assert soup.img and all((out/a['href']).is_file() for a in soup.find_all('a'))
    assert (out/'provenance.json').is_file()
    assert not (out/'REPORT_ZH.html').exists()


def test_workflows_registered():
    registry=WorkflowRegistry()
    assert registry.get('single_cell.import_marvel')
    assert registry.get('single_cell.import_longread_isoforms')


@pytest.mark.parametrize('mode', ['marvel', 'longread_isoforms'])
def test_registered_import_run_and_input_cache_invalidation(tmp_path, mode):
    import json
    from fastapi.testclient import TestClient
    from katana_next.api.app import PlatformRuntime, create_app
    from katana_next.domain.models import JobStatus
    from katana_next.full_fixture import create_full_fixture

    manifest = create_full_fixture(tmp_path / 'project')
    runtime = PlatformRuntime(application_root=tmp_path / 'application')
    client = TestClient(create_app(runtime))
    opened = client.post('/api/projects/open', json={'manifest_path': str(manifest)})
    assert opened.status_code == 200, opened.text
    project_id = opened.json()['snapshot']['project_id']
    config = fixture(tmp_path) if mode == 'marvel' else long_fixture(tmp_path)
    config['project_id'] = project_id
    config_path = tmp_path / 'import.json'
    config_path.write_text(json.dumps(config), encoding='utf-8')
    endpoint = f'/api/projects/{project_id}/workflows/single_cell.import_{mode}'
    request = {'parameters': {'config_path': str(config_path)}}
    started = client.post(endpoint + '/run', json=request)
    assert started.status_code == 200, started.text
    done = runtime.jobs.wait(started.json()['job_id'], timeout=120)
    assert done.status == JobStatus.COMPLETE, done
    assert client.get(f'/api/projects/{project_id}/artifacts').json()
    assert client.post(endpoint + '/cache', json=request).json()['cache_hit'] is True
    meta = pd.read_csv(config['metadata_path'], sep='\t')
    meta.loc[0, 'donor_id'] = 'changed_donor'
    meta.to_csv(config['metadata_path'], sep='\t', index=False)
    assert client.post(endpoint + '/cache', json=request).json()['cache_hit'] is False
