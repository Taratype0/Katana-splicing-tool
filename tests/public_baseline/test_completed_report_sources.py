import copy
import hashlib
import json
from pathlib import Path
import pytest
from katana_next.reports.sources import load_source
from katana_next.reports.path_mapping import PathMappings
from katana_next.reports.storage import copy_payload


def save(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value),encoding='utf-8')
    return path


def setup_sources(tmp):
    project={'project_id':'study','name':'Study','annotation':{'genome_build':'fixture'},
             'samples':[{'sample_id':'a','biological_sample_id':'a','condition':'A','attributes':{}},
                        {'sample_id':'b','biological_sample_id':'b','condition':'B','attributes':{}}],
             'contrasts':[{'contrast_id':'B_vs_A','numerator':{'label':'B','sample_ids':['b']},
                           'denominator':{'label':'A','sample_ids':['a']}}]}
    parent=save(tmp/'project.json',project)
    sources=[]
    for name in ['first','second']:
        table=tmp/(name+'.tsv');table.write_text('gene\tFDR\ng1\t0.003\n',encoding='utf-8')
        desc=save(tmp/(name+'.json'),{'project_manifest':'project.json','artifacts':[
            {'artifact_id':name,'analysis_id':'project.audit','project_id':'study','status':'complete',
             'path':table.name,'kind':'table','contrast_ids':['B_vs_A']}]})
        sources.append({'source_id':name,'kind':'descriptor','path':desc.name})
    config={'schema':'katana.completed_result_set/v1','project_manifest':parent.name,'sources':sources}
    return project,parent,config


def test_multiple_completed_sources_are_read_only(tmp_path):
    _,parent,config=setup_sources(tmp_path)
    before=parent.read_bytes()
    loaded=load_source(save(tmp_path/'sources.json',config))
    assert len(loaded['artifacts'])==2
    assert loaded['assembly_audit']['scientific_jobs_submitted']==0
    assert all(a['metadata']['report_source_binding']['biological_design_identical'] for a in loaded['artifacts'])
    assert parent.read_bytes()==before


def test_swapped_effect_direction_rejected(tmp_path):
    project,parent,config=setup_sources(tmp_path)
    changed=copy.deepcopy(project)
    c=changed['contrasts'][0];c['numerator'],c['denominator']=c['denominator'],c['numerator']
    save(tmp_path/'other.json',changed)
    desc=json.loads((tmp_path/'second.json').read_text('utf-8'));desc['project_manifest']='other.json'
    save(tmp_path/'second.json',desc)
    config['sources'][1]['target_project_manifest_sha256']=hashlib.sha256(parent.read_bytes()).hexdigest()
    with pytest.raises(ValueError,match='contrast arms'):load_source(save(tmp_path/'sources.json',config))


def test_conflicting_artifact_ids_rejected(tmp_path):
    _,_,config=setup_sources(tmp_path)
    desc=json.loads((tmp_path/'second.json').read_text('utf-8'));desc['artifacts'][0]['artifact_id']='first'
    save(tmp_path/'second.json',desc)
    with pytest.raises(ValueError,match='Conflicting completed'):load_source(save(tmp_path/'sources.json',config))


def test_unverified_status_cannot_be_cleared(tmp_path):
    _,_,config=setup_sources(tmp_path)
    config['superseded_statuses']=[{'source_id':'first','title':'Old missing result','replacement_source_ids':['second']}]
    with pytest.raises(ValueError,match='verified completed replacement'):load_source(save(tmp_path/'sources.json',config))


def test_bound_views_supersede_only_exact_old_status(tmp_path):
    _,parent,config=setup_sources(tmp_path)
    native_hash=hashlib.sha256(parent.read_bytes()).hexdigest()
    first=json.loads((tmp_path/'first.json').read_text('utf-8'))
    first['extra_views']=[{'project_id':'study','source_project_manifest_sha256':native_hash,
                          'status':'complete','unavailable':[{'title':'Old missing result','status':'wip','chapter':7}]}]
    save(tmp_path/'first.json',first)
    view=save(tmp_path/'views.json',{'project_id':'study','source_project_manifest_sha256':native_hash,
                                   'status':'complete','panels':[{'contrast_id':'B_vs_A','path':'second.tsv','chapter':7}]})
    evidence=save(tmp_path/'proof.json',{'schema':'external.completed/v1','status':'completed'})
    config['sources'].append({'source_id':'replacement','kind':'views','path':view.name,
                             'completion_evidence':[{'path':evidence.name,'sha256':hashlib.sha256(evidence.read_bytes()).hexdigest()}]})
    config['superseded_statuses']=[{'source_id':'first','title':'Old missing result','replacement_source_ids':['replacement']}]
    loaded=load_source(save(tmp_path/'sources.json',config))
    assert not any(v.get('unavailable') for v in loaded['extra_views'])
    assert loaded['assembly_audit']['superseded_status_history'][0]['original_status']['status']=='wip'
    evidence.write_text('{}',encoding='utf-8')
    with pytest.raises(ValueError,match='evidence SHA'):load_source(tmp_path/'sources.json')


def test_mapping_is_explicit_boundary_checked_and_read_only(tmp_path):
    mapper=PathMappings([{'source_prefix':'E:/study','target_prefix':str(tmp_path)}])
    assert Path(mapper.value(r'e:\study\tables\x.tsv'))==tmp_path/'tables/x.tsv'
    assert mapper.value('E:/study-other/x.tsv')=='E:/study-other/x.tsv'
    with pytest.raises(ValueError,match='escapes'):mapper.value('E:/study/../../outside')


def test_result_store_owns_bytes_and_does_not_link_original(tmp_path):
    source=tmp_path/'original.tsv';source.write_text('FDR\n0.001\n',encoding='utf-8')
    a,b=tmp_path/'a/result.tsv',tmp_path/'b/result.tsv'
    store=tmp_path/'store'
    copy_payload(source,a,store=store);copy_payload(source,b,store=store)
    assert a.read_bytes()==b.read_bytes()==source.read_bytes()
    assert source.stat().st_ino!=a.stat().st_ino
    source.write_text('FDR\n0.999\n',encoding='utf-8')
    assert a.read_text('utf-8')=='FDR\n0.001\n'
    with pytest.raises(ValueError,match='existing report payload'):copy_payload(source,a,store=store)


def test_tampered_result_store_is_rejected(tmp_path):
    source=tmp_path/'source.tsv';source.write_text('safe',encoding='utf-8')
    store=tmp_path/'store';copy_payload(source,tmp_path/'out/a.tsv',store=store)
    next(store.iterdir()).write_text('changed',encoding='utf-8')
    with pytest.raises(ValueError,match='identity differs'):copy_payload(source,tmp_path/'other/a.tsv',store=store)


def test_shared_store_preserves_source_times_without_retiming_another_report(tmp_path):
    import os
    first, second = tmp_path/'first.tsv', tmp_path/'second.tsv'
    first.write_bytes(b'same completed statistics\n')
    second.write_bytes(first.read_bytes())
    os.utime(first, ns=(0, 0))
    os.utime(second, ns=(2_000_000_000, 2_000_000_000))
    a, b = tmp_path/'one/a.tsv', tmp_path/'two/b.tsv'
    copy_payload(first, a, store=tmp_path/'store')
    copy_payload(second, b, store=tmp_path/'store')
    assert a.stat().st_mtime_ns == first.stat().st_mtime_ns == 0
    assert b.stat().st_mtime_ns == second.stat().st_mtime_ns == 2_000_000_000
