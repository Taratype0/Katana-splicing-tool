import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from katana_next.science.protein_domains import scan_domains


@pytest.mark.parametrize('change', ['none', 'corrupt_result', 'incomplete_result', 'runtime', 'sequence'])
def test_verified_project_scan_reuse(tmp_path, monkeypatch, change):
    database = tmp_path/'db.hmm'
    database.write_bytes(b'test database')
    db = SimpleNamespace(path=database, metadata={'sha256':hashlib.sha256(database.read_bytes()).hexdigest(),
        'version':'test', 'source_url':'https://example.org/db'})
    runtime_path = tmp_path/'runtime.json'
    runtime_path.write_text(json.dumps({'version':'1', 'python':'python'}))
    runtime = SimpleNamespace(path=runtime_path)
    features = pd.DataFrame([{'transcript_id':'tx', 'protein_sequence':'MKK',
        'protein_sha256':hashlib.sha256(b'MKK').hexdigest(), 'protein_status':'validated_annotated_ORF'}])
    messages, scans = [], []
    job = SimpleNamespace(check_cancelled=lambda:None, report=lambda _, text:messages.append(text))
    class Process:
        def __init__(self, command, **kwargs):
            records=json.loads(Path(command[-2]).read_text())
            result={'engine_version':json.loads(runtime_path.read_text())['version'], 'profiles_scanned':1,
                'proteins':[{'transcript_id':r['transcript_id'],
                    'protein_sha256':hashlib.sha256(r['protein_sequence'].encode()).hexdigest()} for r in records], 'domains':[]}
            Path(command[-1]).write_text(json.dumps(result))
            scans.append(command)
        def wait(self, timeout):
            return 0
    monkeypatch.setattr('katana_next.science.protein_domains.subprocess.Popen', Process)
    prior=tmp_path/'old'
    first,key=scan_domains(features,db,runtime,prior,job)
    if change=='corrupt_result':
        (prior/'domains.json').write_text('{}')
    elif change=='incomplete_result':
        first['proteins']=[]
        (prior/'domains.json').write_text(json.dumps(first))
        (prior/'checkpoint.json').write_text(json.dumps({'key':key,
            'result_sha256':hashlib.sha256((prior/'domains.json').read_bytes()).hexdigest()}))
    elif change=='runtime':
        runtime_path.write_text(json.dumps({'version':'2', 'python':'python'}))
    elif change=='sequence':
        features.loc[0,'protein_sequence']='MKL'
        features.loc[0,'protein_sha256']=hashlib.sha256(b'MKL').hexdigest()
    result,second_key=scan_domains(features,db,runtime,tmp_path/'new',job,
        cache_root=tmp_path/'shared', previous_directories=[prior])
    assert len(scans)==(1 if change=='none' else 2)
    assert result['proteins'][0]['protein_sha256']==features.iloc[0].protein_sha256
    assert (tmp_path/'new'/'reuse.json').exists() == (change=='none')
    third,third_key=scan_domains(features,db,runtime,tmp_path/'third',job,cache_root=tmp_path/'shared')
    assert third==result and third_key==second_key
    assert len(scans)==(1 if change=='none' else 2)
    assert 'cache hit' in messages[-1]
