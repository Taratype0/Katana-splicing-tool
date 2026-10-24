"""Saved views authenticate against the actual native project before rebinding."""
from copy import deepcopy
from hashlib import sha256
import json

import pytest

from katana_next.reports.sources import load_source


def _save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2), encoding='utf-8')
    return path


def _sha(path):
    return sha256(path.read_bytes()).hexdigest()


def _inputs(tmp):
    project = {'project_id':'study','name':'Study','annotation':{'genome_build':'fixture'},
               'samples':[{'sample_id':'a','biological_sample_id':'a','condition':'A','attributes':{}},
                          {'sample_id':'b','biological_sample_id':'b','condition':'B','attributes':{}}],
               'contrasts':[{'contrast_id':'B_vs_A','numerator':{'label':'B','sample_ids':['b']},
                             'denominator':{'label':'A','sample_ids':['a']}}]}
    target = _save(tmp/'target.json', project)
    native = _save(tmp/'native'/'project.json', {**deepcopy(project),'name':'Original compute revision'})
    (tmp/'native'/'native.tsv').write_bytes(b'model\tP\nAAA\t0.001\n')
    views = _save(tmp/'native'/'views.json', {'project_id':project['project_id'],
        'source_project_manifest_sha256':_sha(native),'status':'complete',
        'panels':[{'id':'saved-view','contrast_id':'B_vs_A','path':'native.tsv','chapter':7,
                   'metadata':{'component':'saved-result'}}]})
    config = {'schema':'katana.completed_result_set/v1','project_manifest':'target.json',
              'sources':[{'source_id':'saved','kind':'views','path':'native/views.json',
                          'project_manifest':'native/project.json','target_project_manifest_sha256':_sha(target)}]}
    return target, native, views, config


def _load(tmp, config):
    return load_source(_save(tmp/'results.json',config))


def test_native_view_is_authenticated_then_explicitly_rebound_without_source_change(tmp_path):
    target,native,views,config = _inputs(tmp_path)
    before={p:p.read_bytes() for p in (target,native,views)}
    loaded=_load(tmp_path,config)
    view=loaded['extra_views'][0]
    binding=view['panels'][0]['metadata']['report_source_binding']
    assert view['project_id']=='study' and view['source_project_manifest_sha256']==_sha(target)
    assert binding['native_manifest_sha256']==_sha(native)
    assert binding['target_manifest_sha256']==_sha(target)
    assert binding['native_manifest_sha256']!=binding['target_manifest_sha256']
    assert binding['biological_design_identical'] is True
    assert binding['scientific_computation_performed'] is False
    assert native in loaded['source_files'] and views in loaded['source_files']
    assert before=={p:p.read_bytes() for p in before}


@pytest.mark.parametrize('mutation',['omit_native_manifest','wrong_native_manifest','missing_target_sha',
                                   'wrong_target_sha','different_samples','swapped_effect','different_reference'])
def test_native_view_rebinding_never_ignores_native_or_target_identity(tmp_path,mutation):
    target,native,views,config = _inputs(tmp_path)
    row=config['sources'][0]
    expected='Extra report views require matching'
    if mutation=='omit_native_manifest':row.pop('project_manifest')
    elif mutation=='wrong_native_manifest':row['project_manifest']='target.json'
    elif mutation=='missing_target_sha':
        row.pop('target_project_manifest_sha256'); expected='exact target manifest SHA'
    elif mutation=='wrong_target_sha':
        row['target_project_manifest_sha256']='a'*64; expected='exact target manifest SHA'
    else:
        project=json.loads(native.read_text('utf-8'))
        if mutation=='different_samples':project['samples'][0]['biological_sample_id']='another-mouse'
        elif mutation=='swapped_effect':
            c=project['contrasts'][0]; c['numerator'],c['denominator']=c['denominator'],c['numerator']
        elif mutation=='different_reference':project['annotation']['genome_build']='another-build'
        _save(native,project)
        view=json.loads(views.read_text('utf-8')); view['source_project_manifest_sha256']=_sha(native)
        _save(views,view)
        expected='different biological samples, contrast arms or reference'
    with pytest.raises(ValueError,match=expected):_load(tmp_path,config)
