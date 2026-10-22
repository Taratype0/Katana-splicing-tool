"""Auditable hub selection over all fixed-module events, independent of drawing."""
from dataclasses import asdict, dataclass
from pathlib import Path
import json, os, subprocess
import numpy as np
import pandas as pd
from katana_next.science.network_labels import network_label_policy as _network_label_policy


NETWORK_LABEL_POLICY_ID = 'network-descriptive-labels-opt-in-v1'


@dataclass(frozen=True)
class HubSettings:
    method: str = 'bicor'
    power: float = 6.
    signed: bool = True
    edge_threshold: float = .10
    bootstrap_iterations: int = 100
    seed: int = 20261018
    minimum_samples: int = 15
    top_fraction: float = .10
    minimum_kme: float = .70
    minimum_selection_frequency: float = .70
    minimum_valid_bootstrap_fraction: float = .80
    minimum_pair_coverage: float = .80
    block_size: int = 256
    apply_threshold_labels: bool = False


def network_label_policy(settings):
    return {
        **_network_label_policy(settings.apply_threshold_labels),
        'label_thresholds': {key:getattr(settings,key) for key in (
            'top_fraction','minimum_kme','minimum_selection_frequency',
            'minimum_valid_bootstrap_fraction','minimum_pair_coverage')},
        'ranking': 'Complete all-edge and cross-gene connectivity ranks; no preset eligibility filter.',
        'preset_frequency_semantics':
            'preset_selection_frequency_* and legacy bootstrap_frequency_* count selection '
            'under the recorded top-fraction/kME/coverage/edge preset. These are descriptive '
            'sensitivity metrics, not threshold-free confidence or default ranking criteria.',
        'bootstrap_valid_semantics':
            'bootstrap_valid_* and bootstrap_valid_fraction_* describe mathematical estimability; '
            'bootstrap_preset_valid_* additionally require the recorded pair-coverage preset.',
    }


def prepare_inputs(matrix, assignments, samples, annotation, settings, scope='all_and_conditions'):
    if type(settings.apply_threshold_labels) is not bool:
        raise ValueError('apply_threshold_labels must be a boolean')
    if settings.method not in {'pearson','bicor','spearman'} or settings.power <= 0:
        raise ValueError('Unsupported source correlation method or power')
    if settings.bootstrap_iterations < 20 or settings.minimum_samples < 3:
        raise ValueError('At least 20 bootstraps and 3 biological units are required')
    for value in (settings.top_fraction,settings.minimum_selection_frequency,settings.minimum_valid_bootstrap_fraction,settings.minimum_pair_coverage):
        if not 0 < value <= 1: raise ValueError('Hub fractions must lie in (0,1]')
    if not 0 <= settings.minimum_kme <= 1 or not 0 < settings.edge_threshold <= 1:
        raise ValueError('Invalid kME or descriptive edge threshold')
    if scope not in {'all','conditions','all_and_conditions'}: raise ValueError('Unknown condition scope')
    required={'sample_id','biological_sample_id','condition'}
    if not required <= set(samples): raise ValueError('Declare biological units and conditions')
    if not matrix.index.is_unique or not matrix.columns.is_unique or not samples.sample_id.is_unique or not assignments.event_key.is_unique:
        raise ValueError('Matrix, sample and event identities must be unique')
    if set(matrix.columns) != set(samples.sample_id): raise ValueError('Sample metadata must exactly match PSI columns')
    if samples[list(required)].isna().any().any() or samples[list(required)].astype(str).apply(lambda x:x.str.strip().eq('')).any().any():
        raise ValueError('Missing biological/sample/condition identity')
    if 'all' in set(samples.condition.astype(str)): raise ValueError('Condition name all is reserved for the discovery-wide reference')
    values=matrix.to_numpy(float)
    if np.isinf(values).any() or np.any(np.isfinite(values)&((values<0)|(values>1))): raise ValueError('Require measured PSI in [0,1]')
    label=next((c for c in ('analysis_module','consensus_module','module','observed_module') if c in assignments),None)
    if label is None: raise ValueError('Frozen module assignments are required')
    labels=assignments.set_index('event_key')[label].astype(str)
    labels=labels.loc[~labels.str.lower().isin({'grey','gray','0','0.0','unassigned'})]
    if labels.empty: raise ValueError('No assigned modules in the completed source')
    if not set(labels.index) <= set(matrix.index): raise ValueError('Do not subset a module to available events')
    if not annotation.event_key.is_unique: raise ValueError('Unique event annotation is required')
    gene=next((c for c in ('gene_id','GeneID','gene_symbol','gene_name') if c in annotation),None)
    annotated=annotation.set_index('event_key').reindex(labels.index)
    host=annotated[gene].fillna('').astype(str) if gene else pd.Series('',index=labels.index)
    symbols=next((c for c in ('gene_symbol','gene_name','gene_id','GeneID') if c in annotated),None)
    events=pd.DataFrame({'event_key':labels.index,'host_gene':host.values,'host_gene_source':gene or 'unavailable','gene':annotated[symbols].fillna('').astype(str).values if symbols else '', 'module':labels.values})
    design=samples.set_index('sample_id').loc[matrix.columns].copy()
    # Only technical repeats within the same unit and condition are averaged.
    groups=design.groupby(['biological_sample_id','condition'],sort=True).groups
    columns=[]; rows=[]; audit=[]
    for i,((unit,condition),ids) in enumerate(groups.items()):
        key=f'unit_{i}'
        columns.append(matrix.loc[labels.index,list(ids)].mean(axis=1).rename(key))
        rows.append(dict(sample_id=key,biological_sample_id=str(unit),condition=str(condition)))
        audit.append(dict(collapsed_sample=key,biological_sample_id=str(unit),condition=str(condition),source_samples=list(ids)))
    measured=pd.concat(columns,axis=1); design=pd.DataFrame(rows)
    rng=np.random.default_rng(settings.seed); scopes=[]
    names=(['all'] if scope in {'all','all_and_conditions'} else [])+(sorted(design.condition.unique()) if scope!='all' else [])
    for name in names:
        chosen=design if name=='all' else design.loc[design.condition.eq(name)]
        profile={unit:tuple(sorted(g.condition.tolist())) for unit,g in chosen.groupby('biological_sample_id')}
        strata={}
        for unit,p in profile.items(): strata.setdefault(p,[]).append(unit)
        boots=[]
        if chosen.biological_sample_id.nunique() >= settings.minimum_samples:
            for b in range(settings.bootstrap_iterations):
                ix=[]; units=[]; draw=0
                for p in sorted(strata):
                    ids=sorted(strata[p])
                    for unit in rng.choice(ids,size=len(ids),replace=True):
                        positions=chosen.index[chosen.biological_sample_id.eq(unit)].tolist()
                        ix.extend(positions); units.extend([f'draw_{draw}']*len(positions)); draw+=1
                boots.append(dict(indices=ix,units=units))
        scopes.append(dict(name=str(name),indices=chosen.index.tolist(),units=chosen.biological_sample_id.tolist(),bootstraps=boots))
    modules=[dict(name=str(m),events=events.loc[events.module.eq(m),'event_key'].tolist()) for m in sorted(events.module.unique())]
    return measured,events,design,dict(scopes=scopes,modules=modules),audit


def summarize_hubs(raw, settings):
    if type(settings.apply_threshold_labels) is not bool:
        raise ValueError('apply_threshold_labels must be a boolean')
    raw=raw.copy()
    rank_cols={'all_edges':('kWithin','pair_coverage','supported_degree'), 'no_same_gene':('kWithin_no_same_gene','cross_gene_pair_coverage','cross_gene_supported_degree')}
    keys=['module','condition','bootstrap']
    effective_kme=raw.kME if settings.signed else raw.kME.abs()
    design_estimable=raw.biological_units.ge(settings.minimum_samples)&raw.module_events.ge(3)
    for mode,(kw,coverage,degree) in rank_cols.items():
        raw['rank_'+mode]=raw.groupby(keys)[kw].rank(ascending=False,method='average')
        estimable=design_estimable&np.isfinite(raw[kw])&np.isfinite(raw.kME)&np.isfinite(raw[coverage])&raw[coverage].gt(0)
        if mode=='no_same_gene':
            estimable &= raw.gene_known.fillna(False)&raw.cross_gene_possible_pairs.gt(0)&raw.gene_annotation_fraction.eq(1.)
        raw['estimable_'+mode]=estimable
        cutoff=np.maximum(1,np.ceil(raw.module_events*settings.top_fraction))
        raw['preset_pass_'+mode]=estimable&(raw['rank_'+mode]<=cutoff)&(effective_kme>=settings.minimum_kme)&(raw[kw]>0)&(raw[degree]>0)&(raw[coverage]>=settings.minimum_pair_coverage)
        raw['pass_'+mode]=raw['preset_pass_'+mode].astype('boolean') if settings.apply_threshold_labels else pd.Series(pd.NA,index=raw.index,dtype='boolean')
        if settings.apply_threshold_labels:
            raw.loc[~estimable,'pass_'+mode]=pd.NA
    baseline=raw.loc[raw.bootstrap.eq(0)].copy()
    boot=raw.loc[raw.bootstrap.gt(0)].copy()
    indexes=['module','condition','event_key']
    for mode,(kw,coverage,degree) in rank_cols.items():
        valid=boot['estimable_'+mode]
        boot['valid_'+mode]=valid
        boot['preset_valid_'+mode]=valid&boot[coverage].ge(settings.minimum_pair_coverage)
        grouped=boot.groupby(indexes)
        counts=grouped.agg(**{
            'bootstrap_attempts_'+mode:('bootstrap','size'),
            'bootstrap_valid_'+mode:('valid_'+mode,'sum'),
            'bootstrap_preset_valid_'+mode:('preset_valid_'+mode,'sum'),
            'bootstrap_selected_'+mode:('preset_pass_'+mode,'sum'),
        }).reset_index()
        baseline=baseline.merge(counts,on=indexes,how='left')
        for prefix in ('bootstrap_attempts_','bootstrap_valid_','bootstrap_preset_valid_','bootstrap_selected_'):
            baseline[prefix+mode]=baseline[prefix+mode].fillna(0).astype(int)
        assessed=baseline['bootstrap_attempts_'+mode].gt(0)
        baseline['bootstrap_valid_fraction_'+mode]=(baseline['bootstrap_valid_'+mode]/settings.bootstrap_iterations).where(assessed)
        baseline['bootstrap_preset_valid_fraction_'+mode]=(baseline['bootstrap_preset_valid_'+mode]/settings.bootstrap_iterations).where(assessed)
        baseline['preset_selection_frequency_'+mode]=(baseline['bootstrap_selected_'+mode]/settings.bootstrap_iterations).where(assessed)
        # Legacy column retained only as a clearly recorded preset-sensitivity alias.
        baseline['bootstrap_frequency_'+mode]=baseline['preset_selection_frequency_'+mode]
        candidate=baseline['preset_pass_'+mode]&baseline['bootstrap_preset_valid_fraction_'+mode].ge(settings.minimum_valid_bootstrap_fraction)&baseline['preset_selection_frequency_'+mode].ge(settings.minimum_selection_frequency)
        candidate_estimable=baseline['estimable_'+mode]&assessed&baseline['bootstrap_valid_'+mode].gt(0)
        baseline['candidate_'+mode]=candidate.astype('boolean') if settings.apply_threshold_labels else pd.Series(pd.NA,index=baseline.index,dtype='boolean')
        if settings.apply_threshold_labels:
            baseline.loc[~candidate_estimable,'candidate_'+mode]=pd.NA
        status_name='all_edge_estimation_status' if mode=='all_edges' else 'cross_gene_estimation_status'
        baseline[status_name]=np.where(baseline['estimable_'+mode],'estimable','not_estimable')
    if len(boot):
        stats=boot.groupby(indexes).agg(bootstrap_kME_q10=('kME',lambda x:x.quantile(.1)),bootstrap_kME_median=('kME','median'),
            bootstrap_kWithin_q10=('kWithin',lambda x:x.quantile(.1)),bootstrap_cross_gene_kWithin_q10=('kWithin_no_same_gene',lambda x:x.quantile(.1)),
            bootstrap_rank_q90=('rank_all_edges',lambda x:x.quantile(.9)),bootstrap_cross_gene_rank_q90=('rank_no_same_gene',lambda x:x.quantile(.9))).reset_index()
        baseline=baseline.merge(stats,on=indexes,how='left')
    else:
        for key in ('bootstrap_kME_q10','bootstrap_kME_median','bootstrap_rank_q90','bootstrap_kWithin_q10','bootstrap_cross_gene_kWithin_q10','bootstrap_cross_gene_rank_q90'): baseline[key]=np.nan
    baseline['rank_change_no_same_gene']=baseline.rank_no_same_gene-baseline.rank_all_edges
    baseline['hub_candidate']=baseline.candidate_all_edges&baseline.candidate_no_same_gene
    estimable=baseline.estimable_all_edges&baseline.estimable_no_same_gene
    # An unassessable arm cannot establish either a positive or a negative candidate label.
    label_estimable=estimable&baseline.candidate_all_edges.notna()&baseline.candidate_no_same_gene.notna()
    baseline.loc[~label_estimable,'hub_candidate']=pd.NA
    baseline['estimation_status']=np.where(estimable,'estimable','not_estimable')
    baseline['label_status']=np.where(label_estimable,'applied','not_estimable') if settings.apply_threshold_labels else 'not_applied'
    baseline['hub_label_status']=baseline['label_status']
    baseline['threshold_labels_applied']=settings.apply_threshold_labels
    def reason(row):
        reasons=[]
        if row.biological_units<settings.minimum_samples: reasons.append('insufficient_biological_units')
        if row.module_events<3: reasons.append('fewer_than_three_module_events')
        if not np.isfinite(row.kWithin): reasons.append('no_estimable_connections')
        if not np.isfinite(row.kME): reasons.append('unavailable_kME')
        if not np.isfinite(row.pair_coverage) or row.pair_coverage<=0: reasons.append('no_estimable_pair_coverage')
        if not row.gene_known or row.cross_gene_possible_pairs==0: reasons.append('cross_gene_comparison_unavailable')
        if row.gene_annotation_fraction<1: reasons.append('incomplete_host_gene_annotation')
        if not np.isfinite(row.kWithin_no_same_gene): reasons.append('no_estimable_cross_gene_connections')
        if not np.isfinite(row.cross_gene_pair_coverage) or row.cross_gene_pair_coverage<=0: reasons.append('no_estimable_cross_gene_pair_coverage')
        return ';'.join(reasons) or 'estimable'
    baseline['estimation_reason']=baseline.apply(reason,axis=1)
    def selection_reason(row):
        if not settings.apply_threshold_labels:
            return 'threshold_labels_not_applied'+(';' + row.estimation_reason if row.estimation_reason!='estimable' else '')
        reasons=[] if row.estimation_reason=='estimable' else [row.estimation_reason]
        if row.supported_degree==0: reasons.append('preset_no_supported_edges')
        if np.isfinite(row.kME) and (row.kME if settings.signed else abs(row.kME))<settings.minimum_kme: reasons.append('preset_low_kME')
        if np.isfinite(row.pair_coverage) and row.pair_coverage<settings.minimum_pair_coverage: reasons.append('preset_insufficient_pair_coverage')
        for mode in rank_cols:
            value=row['candidate_'+mode]
            if pd.isna(value): reasons.append(mode+'_candidate_not_estimable')
            elif not value: reasons.append(mode+'_preset_criteria_not_met')
        return ';'.join(reasons) or 'passed_optional_preset'
    baseline['selection_reason']=baseline.apply(selection_reason,axis=1)
    return baseline,boot


def run_hub_analysis(matrix, assignments, samples, annotation, settings, directory, rscript, *, scope='all_and_conditions', progress=None):
    directory=Path(directory); directory.mkdir(parents=True,exist_ok=True)
    measured,events,design,structure,audit=prepare_inputs(matrix,assignments,samples,annotation,settings,scope)
    measured.to_csv(directory/'hub_input_psi.tsv',sep='\t',index_label='event_key')
    events.to_csv(directory/'hub_input_events.tsv',sep='\t',index=False)
    design.to_csv(directory/'hub_biological_units.tsv',sep='\t',index=False)
    plan={**asdict(settings),**structure,'network_label_policy':network_label_policy(settings),'source_units':design.biological_sample_id.tolist(),'matrix':str(directory/'hub_input_psi.tsv'),'events':str(directory/'hub_input_events.tsv'),
          'output':str(directory/'hub_raw_metrics.tsv'),'edges':str(directory/'hub_display_edges.tsv'),'session':str(directory/'R_session.txt')}
    (directory/'hub_plan.json').write_text(json.dumps(plan,indent=2),encoding='utf-8')
    (directory/'technical_repeat_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    if progress: progress('Computing full-module adjacency, PC1/kME and biological-block bootstrap with official R WGCNA')
    environment=os.environ.copy()
    environment['R_LIBS']=str(Path(__file__).resolve().parents[3]/'r_library')+os.pathsep+environment.get('R_LIBS','')
    result=subprocess.run([str(rscript),'--vanilla',str(Path(__file__).with_name('hub_metrics.R')),str(directory/'hub_plan.json')],capture_output=True,text=True,encoding='utf-8',errors='replace',env=environment)
    (directory/'R_stdout.log').write_text(result.stdout,encoding='utf-8'); (directory/'R_stderr.log').write_text(result.stderr,encoding='utf-8')
    if result.returncode: raise RuntimeError(f'Hub analysis failed (exit {result.returncode}); original results preserved: '+result.stderr[-2000:]+result.stdout[-500:])
    raw=pd.read_csv(directory/'hub_raw_metrics.tsv',sep='\t')
    metrics,boot=summarize_hubs(raw,settings)
    metrics=metrics.merge(events[['event_key','gene','host_gene']],on='event_key',how='left')
    metrics.to_csv(directory/'hub_complete_metrics.tsv',sep='\t',index=False)
    return metrics,boot,plan
