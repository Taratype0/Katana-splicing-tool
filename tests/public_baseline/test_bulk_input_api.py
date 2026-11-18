import json
from fastapi.testclient import TestClient
from katana_next.api.app import PlatformRuntime,create_app

def request(tmp_path):
    gtf=tmp_path/'genes.gtf';gtf.write_text('chr1\ttest\tgene\t1\t100\t.\t+\t.\tgene_id "G1"; gene_biotype "protein_coding";\n')
    deg=tmp_path/'deg.tsv';deg.write_text('gene_id\tlog2FoldChange\tpadj\nG1\t2.5\t0.01\n')
    samples=[{'sample_id':s,'biological_sample_id':s,'condition':s} for s in ['A','B']]
    return {'workspace':str(tmp_path/'new'),'manifest':{'name':'input test','input_root':str(tmp_path),'output_root':str(tmp_path/'new'/'results'),'annotation':{'genome_build':'test','gtf':str(gtf)},'samples':samples,'contrasts':[{'contrast_id':'B_A','label':'B versus A','numerator':{'label':'B','sample_ids':['B']},'denominator':{'label':'A','sample_ids':['A']}}],'assets':[{'asset_id':'deg','kind':'deg_result','path':str(deg),'metadata':{'contrast_id':'B_A'}}]}}

def test_create_bulk_preserves_original_and_does_not_run_workflows(tmp_path):
    runtime=PlatformRuntime(application_root=tmp_path/'app');client=TestClient(create_app(runtime));payload=request(tmp_path)
    original=(tmp_path/'deg.tsv').read_bytes()
    response=client.post('/api/projects/create-bulk',json=payload)
    assert response.status_code==200,response.text
    assert response.json()['snapshot']['sample_count']==2
    saved=tmp_path/'new'/'katana.project.json';before=saved.read_bytes()
    assert (tmp_path/'deg.tsv').read_bytes()==original
    assert client.post('/api/projects/create-bulk',json=payload).status_code==422
    assert saved.read_bytes()==before

def test_unknown_contrast_is_rejected_before_writing(tmp_path):
    client=TestClient(create_app(PlatformRuntime(application_root=tmp_path/'app')));payload=request(tmp_path)
    payload['manifest']['assets'][0]['metadata']['contrast_id']='wrong'
    response=client.post('/api/projects/create-bulk',json=payload)
    assert response.status_code==422
    assert not (tmp_path/'new').exists()
