"""Evidence-bound module candidate gate; never refit or tune to a target count."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from study import sha, verify_frozen

RULE = dict(association_fdr_max=.05, minimum_requested_refits=20,
            minimum_success_fraction=.9, median_member_support_min=.5,
            supported_member_fraction_min=.5, median_match_jaccard_min=.5,
            maximum_unit_variance_fraction_lt=.5)


def summarize_modules(base, refits):
    """Deletion losses count as zero; failed refits are missing, never success."""
    base=base.astype(int); assigned=base[base.gt(0)]; rows=[]; event_rows=[]
    success=sum(r is not None for r in refits)
    for module,members in assigned.groupby(assigned):
        keys=members.index; supports=[]; overlaps=[]
        for labels in refits:
            if labels is None:continue
            labels=labels.astype(int); match=labels.reindex(keys).fillna(0).astype(int)
            counts=match[match.gt(0)].value_counts(); sizes=labels[labels.gt(0)].value_counts()
            supports.append(np.asarray([(counts.get(x,0)-1)/(len(keys)-1) if x>0 and len(keys)>1 else 0. for x in match]))
            overlaps.append(max([n/(len(keys)+sizes[k]-n) for k,n in counts.items()] or [0.]))
        mean=np.mean(supports,axis=0) if supports else np.full(len(keys),np.nan)
        rows.append(dict(module='M'+str(module),stability_events=len(keys),requested_refits=len(refits),successful_refits=success,
                         median_member_support=float(np.median(mean)),
                         supported_member_fraction=float(np.mean(mean>=.5)) if success else np.nan,
                         median_match_jaccard=float(np.median(overlaps)) if success else np.nan))
        event_rows.extend(dict(event_key=k,module=int(module),within_module_coassignment=float(v),
                               requested_replicates=len(refits),successful_replicates=success) for k,v in zip(keys,mean))
    return pd.DataFrame(rows),pd.DataFrame(event_rows)


def bind_stability(run, modules, events, provenance):
    run=Path(run); manifest=verify_frozen(run); target=run/'stability';target.mkdir(exist_ok=True)
    modules.to_csv(target/'module_evidence.tsv',sep='\t',index=False)
    events.to_csv(target/'qualified_event_support.tsv',sep='\t',index=False)
    record=dict(network_id=manifest['network_id'],frozen_files=manifest['frozen_files'],
                full_qc_refits=True,provenance=provenance,
                files={n:sha(target/n) for n in ['module_evidence.tsv','qualified_event_support.tsv']})
    (target/'qualification_manifest.json').write_text(json.dumps(record,indent=2),encoding='utf-8')


def bind_native_stability(run):
    """Only full-QC deletion refits qualify; older bootstrap caches do not."""
    run=Path(run);target=run/'stability'
    if not (target/'summary.json').exists():return False
    summary=json.loads((target/'summary.json').read_text())
    if 'full event QC' not in summary.get('scope',''):return False
    planned=json.loads((target/'planned_draws.json').read_text())['training_sample_draws']
    base=pd.read_csv(run/'model/assignments.tsv',sep='\t',index_col=0).module
    refits=[];sources={}
    for i in range(len(planned)):
        folder=target/f'refit_{i:03}';status=json.loads((folder/'status.json').read_text())
        path=folder/'model/assignments.tsv'
        refits.append(pd.read_csv(path,sep='\t',index_col=0).module if status['status']=='completed' else None)
        if status['status']=='completed':sources[str(path.relative_to(run))]=sha(path)
    modules,events=summarize_modules(base,refits)
    bind_stability(run,modules,events,dict(kind='native_full_qc_deletion_refits',sources=sources))
    return True


def qualification(run, associations, diagnostics):
    run=Path(run);frozen=verify_frozen(run);path=run/'stability/qualification_manifest.json';status='not_measured';evidence=pd.DataFrame()
    if path.exists():
        try:
            manifest=json.loads(path.read_text())
            if manifest['network_id']!=frozen['network_id'] or manifest['frozen_files']!=frozen['frozen_files'] or manifest.get('full_qc_refits') is not True:
                raise ValueError('Network or input mismatch')
            for name,h in manifest['files'].items():
                if Path(name).name!=name or sha(path.parent/name)!=h:raise ValueError('Evidence changed')
            evidence=pd.read_csv(path.parent/'module_evidence.tsv',sep='\t')
            expected=pd.read_csv(run/'study/modules.tsv',sep='\t').set_index('module').events
            if not evidence.module.is_unique or set(evidence.module)!=set(expected.index):raise ValueError('Module universe mismatch')
            if not evidence.set_index('module').stability_events.reindex(expected.index).equals(expected):raise ValueError('Membership count mismatch')
            status='matched_full_qc_refits'
        except (ValueError,KeyError,OSError):status='invalid_or_mismatched_evidence'
    tab=associations.copy()
    if tab.empty:return tab,dict(status=status,rule=RULE,qualified=0)
    from donor38 import METHOD as DONOR_METHOD
    donor_valid=False;donor_reason='not_requested'
    if 'test' in tab and tab.test.eq(DONOR_METHOD).any():
        from inference38 import verify_associations
        donor_valid,donor_reason=verify_associations(run,tab)
    if not evidence.empty and status=='matched_full_qc_refits':tab=tab.merge(evidence,on='module',how='left',validate='many_to_one')
    if 'module' in diagnostics:tab=tab.merge(diagnostics[['module','maximum_unit_variance_fraction']],on='module',how='left',validate='many_to_one')
    for col in ['median_member_support','supported_member_fraction','median_match_jaccard','requested_refits','successful_refits','maximum_unit_variance_fraction']:
        if col not in tab:tab[col]=np.nan
    reasons=[]
    for r in tab.to_dict('records'):
        why=[]
        if r['status']!='estimated' or not np.isfinite(r.get('q_value',np.nan)):why.append('association_not_estimable')
        elif r['q_value']>RULE['association_fdr_max']:why.append('association_FDR_above_0.05')
        if r.get('inference_gate_eligible') is not True:why.append('inference_not_calibrated_for_candidate_gate')
        if r.get('test')==DONOR_METHOD:
            if not donor_valid:why.append('donor_calibration_invalid_or_mismatched')
            if not np.isfinite(r.get('p_fwer',np.nan)) or r['p_fwer']>.05:why.append('donor_Holm_FWER_above_0.05')
        if status!='matched_full_qc_refits':why.append('stability_'+status)
        else:
            if not r['requested_refits']>=RULE['minimum_requested_refits'] or not r['successful_refits']/max(r['requested_refits'],1)>=RULE['minimum_success_fraction']:why.append('insufficient_successful_refits')
            if not r['median_member_support']>=RULE['median_member_support_min']:why.append('low_member_coassignment')
            if not r['supported_member_fraction']>=RULE['supported_member_fraction_min']:why.append('too_few_supported_members')
            if not r['median_match_jaccard']>=RULE['median_match_jaccard_min']:why.append('low_module_match_Jaccard')
        if not np.isfinite(r['maximum_unit_variance_fraction']):why.append('unit_influence_not_measured')
        elif r['maximum_unit_variance_fraction']>=RULE['maximum_unit_variance_fraction_lt']:why.append('single_unit_dominated')
        reasons.append(';'.join(why) if why else 'passes_operational_candidate_gate')
    tab['qualification_reason']=reasons;tab['qualified_candidate']=tab.qualification_reason.eq('passes_operational_candidate_gate')
    tab=tab.sort_values(['q_value','p_value','module'],kind='stable',na_position='last',
        key=lambda s:pd.to_numeric(s.astype(str).str.replace('M','',regex=False)) if s.name=='module' else s).reset_index(drop=True)
    tab['significance_rank']=np.arange(1,len(tab)+1)
    return tab,dict(status=status,rule=RULE,qualified=int(tab.qualified_candidate.sum()),
                    donor_calibration_status=donor_reason,
                    scope='Operational stability-qualified association candidates; no external validation or causal claim. NA and mismatched evidence never pass.',
                    correction='All estimable requested tests are adjusted before stability filtering; FDR is not recalculated on selected modules.')
