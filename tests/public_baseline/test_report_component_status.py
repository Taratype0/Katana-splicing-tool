"""Final scoped completion evidence supersedes only the matching deferred gap."""
from copy import deepcopy
from hashlib import sha256
from types import SimpleNamespace

import pytest

from katana_next.reports.component_status import reconcile_component_statuses


def _state(family='clip', **updates):
    value={'id':'rich-'+family+'-not-computed','status':'not_computed','contrast_id':'all','cell_type':'all',
           'metadata':{'component':family}}
    value.update(updates)
    return value


def _panel(root,cid,**metadata):
    links=[]
    for role in ('plot_data','full_results','foreground','background','methods'):
        path=root/(cid+'-'+role+'.tsv')
        path.write_bytes(b'comparison\tfdr\n'+cid.encode()+b'\t0.9\n')
        links.append({'role':role,'contrast_id':cid,'path':path.name,
                      'source_sha256':sha256(path.read_bytes()).hexdigest()})
    meta={'analysis_id':'prioritize.binding_enrichment','report_module':'clip','completion_verified':True}
    meta.update(metadata)
    return {'id':'clip-'+cid,'contrast_id':cid,'cell_type':'all','metadata':meta,
            'data_links':links,'output_only':True}


def _builder(root,states=None):
    panels=[_panel(root,'B_vs_A'),_panel(root,'C_vs_A')]
    return SimpleNamespace(output=root,panels=panels,tables=[],derived_audits=[],
        source={'project':{'project_id':'study','contrasts':[{'contrast_id':'B_vs_A'},{'contrast_id':'C_vs_A'}]},
                'native_project_manifest_sha256':'a'*64},exclusions=states if states is not None else [_state()])


def test_completed_nonsignificant_clip_resolves_all_selected_scopes_and_preserves_history(tmp_path):
    b=_builder(tmp_path)
    original=deepcopy(b.exclusions)
    before={p:p.read_bytes() for p in tmp_path.iterdir()}
    audit=reconcile_component_statuses(b)
    assert b.exclusions==[] and audit==b.derived_audits[-1]
    decision=audit['decisions'][0]
    assert decision['original_status']==original[0]
    assert decision['action']=='resolved_completed' and decision['status']=='completed'
    assert decision['significance_used_as_completion_gate'] is False
    assert len(decision['verified_evidence'])==2
    assert all(len(e['bound_native_files'])==4 for e in decision['verified_evidence'])
    assert before=={p:p.read_bytes() for p in before}
    count=len(b.derived_audits)
    assert not reconcile_component_statuses(b)['decisions']
    assert len(b.derived_audits)==count


def test_partial_comparison_scope_does_not_erase_global_required_gap(tmp_path):
    b=_builder(tmp_path);b.panels.pop()
    audit=reconcile_component_statuses(b)
    assert len(b.exclusions)==1
    assert audit['decisions'][0]['missing_contrast_ids']==['C_vs_A']


@pytest.mark.parametrize('mutation',['other_cell','aggregate_panel','other_analysis','sequence_only',
                                    'not_verified','explicit_incomplete','wrong_target_project','modified_source'])
def test_foreign_or_incomplete_layer_never_supersedes_clip(tmp_path,mutation):
    b=_builder(tmp_path)
    for panel in b.panels:
        meta=panel['metadata']
        if mutation=='other_cell':panel['cell_type']='cell-A'
        elif mutation=='aggregate_panel':panel['contrast_id']='all'
        elif mutation=='other_analysis':meta['analysis_id']='prioritize.motiftrace_consequences'
        elif mutation=='sequence_only':meta['completion_scope']='sequence_consequences_only'
        elif mutation=='not_verified':meta.pop('completion_verified')
        elif mutation=='explicit_incomplete':meta['completion_verified']=False
        else:
            meta['report_source_binding']={'schema':'katana.report_source_binding/v1',
                'biological_design_identical':True,'source_files_modified':False,'scientific_computation_performed':False,
                'target_project_id':'study','target_manifest_sha256':'a'*64}
            if mutation=='wrong_target_project':meta['report_source_binding']['target_project_id']='other'
            elif mutation=='modified_source':meta['report_source_binding']['source_files_modified']=True
    assert reconcile_component_statuses(b)['decisions'][0]['action']=='retained_unresolved_required_scope'
    assert b.exclusions


@pytest.mark.parametrize('role',['full_results','foreground','background','methods'])
def test_clip_requires_all_native_evidence_roles(tmp_path,role):
    b=_builder(tmp_path)
    b.panels[0]['data_links']=[link for link in b.panels[0]['data_links'] if link['role']!=role]
    assert reconcile_component_statuses(b)['decisions'][0]['missing_contrast_ids']==['B_vs_A']
    assert b.exclusions


@pytest.mark.parametrize('mutation',['missing','hash','comparison','traversal','absolute','cell'])
def test_missing_tampered_or_unbounded_native_file_retains_gap(tmp_path,mutation):
    b=_builder(tmp_path)
    link=next(l for l in b.panels[0]['data_links'] if l['role']=='full_results')
    if mutation=='missing':(tmp_path/link['path']).unlink()
    elif mutation=='hash':link['source_sha256']='b'*64
    elif mutation=='comparison':link['contrast_id']='other'
    elif mutation=='traversal':link['path']='../outside.tsv'
    elif mutation=='absolute':link['path']=str(tmp_path/link['path'])
    elif mutation=='cell':link['cell_type']='other'
    assert reconcile_component_statuses(b)['decisions'][0]['missing_contrast_ids']==['B_vs_A']
    assert b.exclusions


def test_exact_cell_scope_and_selected_comparison_can_resolve_without_other_cells(tmp_path):
    b=_builder(tmp_path,[_state(contrast_id='B_vs_A',cell_type='cell-A')])
    b.panels[0]['cell_type']='cell-A'
    audit=reconcile_component_statuses(b)
    assert not b.exclusions and audit['decisions'][0]['required_contrast_ids']==['B_vs_A']


def test_unsupplied_conditional_parallel_is_audited_as_absent_not_completed(tmp_path):
    state=_state('mapp_parallel')
    state['metadata'].update(requirement='conditional_import',optional_source_supplied=False)
    b=_builder(tmp_path,[state])
    b.panels[0]['metadata'].update(analysis_id='prioritize.mapp_native',calibrated_motiftrace_parallel=False)
    audit=reconcile_component_statuses(b)
    decision=audit['decisions'][0]
    assert not b.exclusions
    assert decision['status']=='not_computed' and decision['optional_source_supplied'] is False
    assert decision['original_status']==state
    assert decision['action']=='optional_absence_retained_in_provenance'


@pytest.mark.parametrize('supplied',['flag','incomplete_panel','source_artifact'])
def test_supplied_but_unproven_optional_parallel_still_retains_gap(tmp_path,supplied):
    state=_state('mapp_parallel')
    state['metadata'].update(requirement='conditional_import',optional_source_supplied=False)
    b=_builder(tmp_path,[state])
    if supplied=='flag':state['metadata']['optional_source_supplied']=True
    elif supplied=='incomplete_panel':b.panels[0]['metadata']['analysis_id']='prioritize.motiftrace_parallel'
    else:b.source['artifacts']=[{'analysis_id':'prioritize.motiftrace_parallel'}]
    assert reconcile_component_statuses(b)['decisions'][0]['action']=='retained_unresolved_required_scope'
    assert b.exclusions


def test_only_actual_calibrated_parallel_family_resolves_the_parallel_status(tmp_path):
    b=_builder(tmp_path,[_state('mapp_parallel')])
    for panel in b.panels:
        panel['metadata'].update(analysis_id='prioritize.motiftrace_parallel',component='mapp_parallel',
                                 report_module='mapp_parallel',calibrated_motiftrace_parallel=True,
                                 completion_scope='completed_calibrated_parallel')
    assert reconcile_component_statuses(b)['decisions'][0]['action']=='resolved_completed'
    assert not b.exclusions


@pytest.mark.parametrize('state',[_state(status='failed'),_state(status='wip'),_state('unregistered_required')])
def test_required_failures_and_unknown_components_are_never_hidden(tmp_path,state):
    b=_builder(tmp_path,[deepcopy(state)])
    assert not reconcile_component_statuses(b)['decisions']
    assert b.exclusions==[state]
