"""Study-wide discovery contract and interpretation of frozen sample-level modules.

No network fitting is imported by the interpretation path. Traits and contrasts
are deliberately absent from the discovery identity. Inference is exploratory:
modules are learned on these samples and robust covariance is not validation.
"""
import hashlib, html, itertools, json, sys, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from cosplice_research import write_json

CORE = ['sample_id', 'biological_sample_id', 'condition']

def sha(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,ensure_ascii=False,allow_nan=False).encode()).hexdigest()

def validate_contract(cfg, meta):
    if cfg.get('value_kind','usage')!='usage':raise ValueError('Study discovery requires individual sample PSI/usage; response profiles are auxiliary')
    if not cfg.get('study_id','').strip() or not cfg.get('cohort_scope','').strip() or cfg.get('cohort_complete') is not True:
        raise ValueError('Declare study_id, comparable cohort_scope and cohort_complete=true for all eligible study samples')
    if not set(CORE).issubset(meta) or meta[CORE].isna().any().any() or not meta.sample_id.is_unique:
        raise ValueError('Study metadata requires unique samples and nonmissing condition/unit identities')
    if meta[CORE].astype(str).apply(lambda x:x.str.strip().eq('').any()).any():raise ValueError('Empty study sample identities')
    if meta.duplicated(['biological_sample_id','condition']).any():raise ValueError('Repeated unit/condition: reconcile technical duplicates before discovery')
    from katana_cosplice.study_design import validate_study_design
    sample_policy=validate_study_design(meta,paired=cfg.get('paired',False))
    trajectory=cfg.get('cohort_design')=='differentiation_trajectory'
    if trajectory:
        if 'time_day' not in meta:raise ValueError('Differentiation trajectories require numeric time_day metadata')
        days=pd.to_numeric(meta.time_day,errors='coerce')
        if not np.isfinite(days).all() or meta.assign(time_day=days).duplicated(['biological_sample_id','time_day']).any():
            raise ValueError('Differentiation trajectories require unique numeric time points per donor')
        if meta.assign(time_day=days).groupby('biological_sample_id').time_day.nunique().lt(2).any():
            raise ValueError('Each differentiation donor must have at least two time points')
        if 'cell_type' in meta and (meta.cell_type.isna().any() or meta.assign(time_day=days).groupby('time_day').cell_type.nunique().gt(1).any()):
            raise ValueError('Different cell types at the same time require stratified analysis, not automatic trajectory pooling')
    for field in ['study_id','tissue','cell_type','assay','genome_build','quantification']:
        if field=='cell_type' and trajectory:continue
        if field in meta and (meta[field].isna().any() or meta[field].nunique()!=1):
            raise ValueError(f'Incompatible or missing {field}: define a comparable cohort explicitly before pooling')
    if 'study_id' in meta and str(meta.study_id.iloc[0])!=cfg['study_id']:raise ValueError('Study ID differs from metadata')
    return {'study_id':cfg['study_id'].strip(),'cohort_scope':cfg['cohort_scope'].strip(),
        'cohort_design':cfg.get('cohort_design','comparable_samples'),
        'cell_type_labels':sorted(meta.cell_type.unique().tolist()) if 'cell_type' in meta else [],
        'completeness':'User-declared all eligible comparable study samples; exact matrix/manifest match checked; unavailable samples cannot be inferred',
        'samples':len(meta),'biological_units':int(meta.biological_sample_id.nunique()),
        'conditions':{str(k):{'samples':len(g),'units':int(g.biological_sample_id.nunique())} for k,g in meta.groupby('condition',sort=True)},
        'repeated_observations':len(meta)-int(meta.biological_sample_id.nunique()),
        'discovery_scope':'One pooled sample-level network; no group averaging, no DEG/DAS significance filter',
        'scope_warning':'Pooled co-response is not within-condition regulation or network rewiring',
        'sample_size_note':sample_policy['sample_size_note'],'sample_size_policy':sample_policy}

def discovery_identity(cfg, measured, meta):
    """Same data + discovery inputs => same key, independent of traits/contrasts."""
    x=measured.sort_index().sort_index(axis=1)
    m=meta.set_index('sample_id').loc[x.columns,CORE[1:]].astype(str)
    # CSV canonicalizes NaNs and row/column ordering, preserving individual values.
    values=hashlib.sha256(x.to_csv(sep='\t',na_rep='NA',float_format='%.17g').encode()).hexdigest()
    source=Path(__file__).parent
    engines={p.name:sha(p) for p in [source/n for n in ['cosplice_research.py','engines.R','stable_network.py','unit_robust.py','audit_labels24.py','module_merge.py','merge_modules.R','advanced.py','stacking.py','product.py']] if p.exists()}
    if cfg['method']=='factor_core':engines['algorithm28.py']=sha(source/'algorithm28.py')
    if cfg['method']=='global_local':engines['algorithm31.py']=sha(source/'algorithm31.py')
    if cfg.get('feature_cap') is not None:engines['feature_selection36.py']=sha(source/'feature_selection36.py')
    extra={k:sha(cfg[k]) for k in ['expression','annotation'] if cfg.get(k) and (k=='expression' and cfg['method']=='mofa' or k=='annotation' and cfg['method']=='stacking')}
    payload={'matrix':values,'sample_units_and_conditions':m.reset_index().to_dict('records'),
        'method':cfg['method'],'merge_modules':bool(cfg.get('merge_modules')),'engine_sources':engines,'discovery_inputs':extra,
        'wgcna_parameters':{'min_module_size':cfg.get('min_module_size',20),'merge_cut_height':cfg.get('merge_cut_height',.25)},
        'qc':{'profile':cfg.get('qc_profile','legacy'),'observed_fraction':.8,'variance_gt':1e-12,'event_cap':cfg.get('feature_cap'),'imputation':'event mean over discovery samples'},'schema':2}
    return digest(payload),payload

def score_modules(measured, assignments):
    scores=pd.DataFrame(index=measured.columns);coverage=pd.DataFrame(index=measured.columns);rows=[];basis=[]
    labels=assignments.reindex(measured.index).module
    for module in sorted(labels.dropna().unique()):
        if str(module) in {'0','0.0','grey'}:continue
        block=measured.loc[labels.eq(module)];name='M'+str(int(module))
        means=block.mean(axis=1);scales=block.std(axis=1).clip(lower=1e-8)
        z=((block.T.fillna(means)-means)/scales).to_numpy()
        u,s,v=np.linalg.svd(z,full_matrices=False);load=v[0];raw=z@load
        sign=1 if float(load.sum())>=0 else -1
        scale=max(float(np.std(raw,ddof=1)),1e-12);load=load*sign/scale
        coverage[name]=block.notna().mean(axis=0)
        scores[name]=pd.Series(z@load,index=measured.columns).where(coverage[name].ge(.8))
        rows.append({'module':name,'events':len(block),'pc1_variance_fraction':float(s[0]**2/(s@s)),
            'score_units':'SD of study PC1; sign aligned with mean standardized PSI','evaluable_samples':int(scores[name].notna().sum())})
        basis.extend({'event_key':key,'module':name,'center':float(means[key]),'scale':float(scales[key]),'loading':float(load[j])} for j,key in enumerate(block.index))
    scores.index.name=coverage.index.name='sample_id'
    return scores,coverage,pd.DataFrame(rows,columns=['module','events','pc1_variance_fraction','score_units','evaluable_samples']),pd.DataFrame(basis,columns=['event_key','module','center','scale','loading'])

def initialize(out,cfg,measured,meta,assignments):
    out=Path(out);target=out/'study';target.mkdir(exist_ok=False)
    contract=validate_contract(cfg,meta);key,identity=discovery_identity(cfg,measured,meta)
    scores,cov,modules,basis=score_modules(measured,assignments)
    for name,tab,index in [('scores',scores,True),('coverage',cov,True),('module_basis',basis,False),('modules',modules,False)]:tab.to_csv(target/(name+'.tsv'),sep='\t',index=index)
    assignment_hash=sha(out/'model/assignments.tsv')
    contract.update(discovery_id=key,network_id=digest({'discovery_id':key,'assignments':assignment_hash}),
        assignment_sha256=assignment_hash,discovery_identity=identity,
        frozen_files={name:sha(out/name) for name in ['usage.tsv.gz','samples.tsv','model/assignments.tsv','study/scores.tsv','study/coverage.tsv','study/module_basis.tsv']})
    write_json(target/'manifest.json',contract)
    if (out/'stability/summary.json').exists():
        from module_qualification import bind_native_stability
        bind_native_stability(out)
    return contract

def verify_frozen(run):
    contract=json.loads((run/'study/manifest.json').read_text(encoding='utf-8'))
    for name,expected in contract['frozen_files'].items():
        if sha(run/name)!=expected:raise ValueError(f'Frozen discovery artifact changed: {name}')
    return contract

def design(meta,fields,check_residual_df=True):
    x=pd.DataFrame({'intercept':np.ones(len(meta))},index=meta.index);terms={}
    for name,kind in fields.items():
        if name not in meta:raise ValueError(f'Metadata field absent: {name}')
        if kind=='continuous':
            column=pd.to_numeric(meta[name],errors='raise').astype(float)
            if not np.isfinite(column).all():raise ValueError(f'Missing/nonfinite covariate: {name}')
            part=column.to_frame('numeric:'+name)
        elif kind=='categorical':
            if meta[name].isna().any() or meta[name].astype(str).str.strip().eq('').any():raise ValueError(f'Missing categorical covariate: {name}')
            part=pd.get_dummies(meta[name].astype(str),prefix=name,drop_first=True,dtype=float)
            if part.shape[1]==0:raise ValueError(f'Constant categorical field: {name}')
        else:raise ValueError('Trait/covariate kind must be categorical or continuous')
        terms[name]=list(range(x.shape[1],x.shape[1]+part.shape[1]));x=pd.concat([x,part],axis=1)
    if not np.isfinite(x).all().all() or np.linalg.matrix_rank(x.to_numpy())!=x.shape[1]:raise ValueError('Confounded/rank-deficient design; effects cannot be separated')
    if check_residual_df and len(x)-x.shape[1]<3:raise ValueError('Fewer than three residual degrees of freedom')
    return x,terms

def robust_fit(y,x,units):
    model=sm.OLS(y,x).fit();repeated=units.duplicated().any()
    if repeated:
        fitted=model.get_robustcov_results(cov_type='cluster',groups=units,use_correction=True,df_correction=True,use_t=True)
        name='OLS + unit-cluster robust covariance (small-sample correction)'
    else:
        fitted=model.get_robustcov_results(cov_type='HC3',use_t=True);name='OLS + HC3 covariance'
    return fitted,name

def adjust(rows):
    tab=pd.DataFrame(rows);tab['q_value']=np.nan
    if 'p_value' not in tab:return tab
    ids=tab.index[np.isfinite(tab.p_value.astype(float))];p=tab.loc[ids,'p_value'].to_numpy(float)
    if np.any((p<0)|(p>1)):raise ValueError('P values must be in [0,1]')
    if len(p):
        order=np.argsort(p);q=np.empty(len(p));q[order]=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];tab.loc[ids,'q_value']=np.minimum(q,1)
    # Stable exports and report rows use the same significance ordering. A module
    # ID is an identity, not its rank. Unestimable tests remain explicit NA rows.
    return tab.sort_values(['q_value','p_value','module'],kind='stable',na_position='last',
        key=lambda s:pd.to_numeric(s.astype(str).str.replace('M','',regex=False)) if s.name=='module' else s).reset_index(drop=True)

def welch_oneway(values,meta,factor):
    """Secondary parametric diagnostic; not the small-N primary P value."""
    from statsmodels.stats.oneway import anova_oneway
    groups=[values.loc[g.index].to_numpy() for _,g in meta.groupby(factor,sort=True)]
    if any(len(g)<3 or not np.isfinite(g).all() or np.var(g,ddof=1)<1e-12 for g in groups):
        raise ValueError('Welch requires >=3 observations and nonzero variance in every group')
    test=anova_oneway(groups,use_var='unequal',welch_correction=True)
    return dict(test='Welch one-way ANOVA + Satterthwaite degrees of freedom',
                statistic=float(test.statistic),p_value=float(test.pvalue),df=float(test.df_denom),df_num=float(test.df_num))

def omnibus(values,meta,factor,kind,covariates,label_blind=True):
    """Calibrated label test for independent one-way designs; explicit robust models otherwise."""
    if kind=='categorical':
        counts=meta.groupby(factor).biological_sample_id.nunique()
        if len(counts)<2 or counts.min()<3:
            raise ValueError('Each observed categorical level needs at least three biological units after score coverage exclusion')
    if kind=='categorical' and not covariates and meta.biological_sample_id.is_unique:
        if not label_blind:
            return dict(welch_oneway(values,meta,factor),inference_gate_eligible=False,
                        inference_limit='Discovery label-invariance not audited for this engine; approximate P remains exploratory')
        from permutation_inference import permutation_anova
        result=permutation_anova(values,meta[factor])
        result['inference_gate_eligible']=True
        try:
            diagnostic=welch_oneway(values,meta,factor)
            result.update(welch_diagnostic_p=diagnostic['p_value'],welch_diagnostic_df=diagnostic['df'])
        except ValueError:result.update(welch_diagnostic_p=np.nan,welch_diagnostic_df=np.nan)
        return result
    x,terms=design(meta,{factor:kind,**covariates});fit,name=robust_fit(values,x,meta.biological_sample_id)
    restriction=np.eye(x.shape[1])[terms[factor]];cov=restriction@fit.cov_params()@restriction.T
    if not np.isfinite(cov).all() or np.linalg.matrix_rank(cov)!=len(restriction) or np.max(np.diag(cov))<1e-24:
        raise ValueError('Nonestimable/singular trait covariance')
    test=fit.f_test(restriction)
    result=dict(test=name,statistic=float(test.fvalue),p_value=float(test.pvalue),df=float(test.df_denom),df_num=float(test.df_num),
                inference_gate_eligible=False,inference_limit='Adjusted/cluster inference requires design-specific calibration before candidate qualification')
    if kind=='continuous':result['effect']=float(fit.params[terms[factor][0]])
    return result

def contrast(y,meta,factor,a,b,covariates,paired=False):
    planned=meta.loc[meta[factor].astype(str).isin([a,b])]
    if planned.duplicated(['biological_sample_id',factor]).any():raise ValueError('Multiple observations per unit/contrast level; explicit longitudinal model needed')
    expected_a=set(planned.loc[planned[factor].astype(str).eq(a),'biological_sample_id']);expected_b=set(planned.loc[planned[factor].astype(str).eq(b),'biological_sample_id'])
    if not isinstance(paired,bool):raise ValueError('paired must be true or false')
    if not paired and expected_a & expected_b:raise ValueError('Repeated units require explicit paired mode')
    if paired and not expected_a & expected_b:raise ValueError('Paired contrast has no shared biological units')
    if paired and expected_a!=expected_b:raise ValueError('Partially paired contrast: specify a longitudinal model; independent test disabled')
    sel=meta[factor].astype(str).isin([a,b]) & y.notna();m=meta.loc[sel].copy();v=y.loc[sel]
    units=m.biological_sample_id;levels=m[factor].astype(str)
    n={c:int(units[levels.eq(c)].nunique()) for c in [a,b]}
    result={'samples':len(m),'biological_units':int(units.nunique()),'units_a':n[a],'units_b':n[b]}
    if min(n.values())<3:raise ValueError('Each contrast level needs at least three biological units')
    if m.duplicated(['biological_sample_id',factor]).any():raise ValueError('Multiple observations per unit/contrast level; explicit longitudinal model needed')
    if not paired and not covariates:
        left=v[levels.eq(a)];right=v[levels.eq(b)]
        va=left.var(ddof=1)/len(left);vb=right.var(ddof=1)/len(right)
        effect=float(left.mean()-right.mean());variance=va+vb
        if not np.isfinite(variance) or variance<1e-24:
            return dict(result,effect=effect,statistic=np.nan,p_value=np.nan,ci_low=np.nan,ci_high=np.nan,df=np.nan,
                        paired_units=0,test='Welch independent-sample t test',reason='Zero/near-zero effect variance')
        df=variance**2/(va**2/(len(left)-1)+vb**2/(len(right)-1))
        se=np.sqrt(variance);t=effect/se;half=stats.t.ppf(.975,df)*se
        return dict(result,effect=effect,statistic=float(t),p_value=float(2*stats.t.sf(abs(t),df)),
                    ci_low=float(effect-half),ci_high=float(effect+half),df=float(df),paired_units=0,test='Welch independent-sample t test')
    if paired:
        # Difference within the same unit; covariates enter as within-unit differences.
        d=pd.DataFrame({'y':v,'unit':units,'level':levels}).pivot(index='unit',columns='level',values='y').reindex(columns=[a,b]).dropna()
        if len(d)<3:raise ValueError('Fewer than three complete score pairs after coverage exclusion')
        m=m[m.biological_sample_id.isin(d.index)]
        result.update(samples=2*len(d),biological_units=len(d),units_a=len(d),units_b=len(d),planned_pairs=len(expected_a),excluded_incomplete_pairs=len(expected_a)-len(d))
        delta=d[a]-d[b];x=pd.DataFrame({'intercept':1.},index=delta.index)
        for field,kind in covariates.items():
            if kind!='continuous':
                t=m.pivot(index='biological_sample_id',columns=factor,values=field)
                if not t[a].equals(t[b]):raise ValueError('Changing categorical covariate in paired contrast requires an explicit model')
                continue  # time-invariant attributes cancel in the paired difference.
            t=m.assign(**{field:pd.to_numeric(m[field],errors='raise')}).pivot(index='biological_sample_id',columns=factor,values=field)
            change=t[a]-t[b]
            if np.any(abs(change)>1e-12):x[field]=change
        if not np.isfinite(x).all().all() or np.linalg.matrix_rank(x)!=x.shape[1]:raise ValueError('Confounded within-unit contrast')
        if len(x)-x.shape[1]<2:raise ValueError('Insufficient paired residual degrees of freedom')
        fitted=sm.OLS(delta,x).fit();test=fitted.t_test(np.eye(x.shape[1])[0]);result['test']='paired difference OLS (t test when no changing covariates)'
        result['paired_units']=len(d)
    else:
        m['_contrast_a']=levels.eq(a).astype(float)
        fields={'_contrast_a':'continuous',**covariates};x,terms=design(m,fields)
        fitted,name=robust_fit(v,x,units);test=fitted.t_test(np.eye(x.shape[1])[terms['_contrast_a'][0]])
        result['test']=name;result['paired_units']=0
    if not np.isfinite(np.asarray(test.sd)).all() or float(np.asarray(test.sd).ravel()[0])<1e-12:
        result.update(effect=float(np.asarray(test.effect).ravel()[0]),statistic=np.nan,p_value=np.nan,ci_low=np.nan,ci_high=np.nan,
            df=float(test.df_denom),reason='Zero/near-zero effect variance; inferential statistic and interval are not estimable')
        return result
    ci=np.asarray(test.conf_int()).ravel()
    result.update(effect=float(np.asarray(test.effect).ravel()[0]),statistic=float(np.asarray(test.statistic).ravel()[0]),
        p_value=float(np.asarray(test.pvalue).ravel()[0]),ci_low=float(ci[0]),ci_high=float(ci[1]),df=float(test.df_denom))
    return result

def analyze(run,cfg,out):
    mode=cfg.get('inference_mode','auto')
    if mode not in ['auto','legacy','donor_trajectory']:raise ValueError('Unknown inference mode')
    if mode=='donor_trajectory' or (mode=='auto' and (Path(run)/'inference38/manifest.json').exists()):
        from inference38 import compatible,analyze as donor_analyze
        if mode=='donor_trajectory' or compatible(cfg):return donor_analyze(run,cfg,out)
    run=Path(run);out=Path(out);out.mkdir(parents=True,exist_ok=False)
    contract=verify_frozen(run)
    meta=pd.read_csv(run/'samples.tsv',sep='\t',dtype={k:str for k in CORE}).set_index('sample_id')
    if cfg.get('traits_path'):
        extra=pd.read_csv(cfg['traits_path'],sep='\t',dtype={k:str for k in CORE}).set_index('sample_id')
        if not extra.index.is_unique or set(extra.index)!=set(meta.index):raise ValueError('Trait table must match every frozen study sample exactly once')
        extra=extra.loc[meta.index]
        for key in CORE[1:]:
            if key in extra and not extra[key].equals(meta[key]):raise ValueError(f'Frozen {key} cannot change during interpretation')
        for key in extra:meta[key]=extra[key]
    scores=pd.read_csv(run/'study/scores.tsv',sep='\t',index_col=0,dtype={'sample_id':str}).reindex(meta.index)
    traits=cfg.get('traits',{'condition':'categorical'});covariates=cfg.get('covariates') or {}
    if set(traits)&set(covariates):raise ValueError('Trait and covariate names overlap')
    if any(k in {'sample_id','biological_sample_id','_contrast_a'} for k in [*traits,*covariates]):raise ValueError('Identity/reserved fields are not traits')
    # Invalid global designs fail visibly before any p-values are produced.
    for trait,kind in traits.items():design(meta,{trait:kind,**covariates},check_residual_df=False)
    contrasts=cfg.get('contrasts',[])
    if contrasts is None:
        contrasts=[{'factor':f,'a':a,'b':b} for f,k in traits.items() if k=='categorical' for a,b in itertools.combinations(sorted(meta[f].astype(str).unique()),2)]
    identities=[(c['factor'],c['a'],c['b']) for c in contrasts]
    if len(identities)!=len(set(identities)):raise ValueError('Duplicate contrasts')
    for c in contrasts:
        if traits.get(c['factor'])!='categorical' or c['a']==c['b'] or not {c['a'],c['b']}.issubset(set(meta[c['factor']].astype(str))):raise ValueError('Contrast levels/factor absent or identical')
    restricted={}
    if cfg.get('permutation_strata'):
        if covariates:raise ValueError('Use categorical permutation strata OR covariate regression')
        if contract['discovery_identity']['method'] not in ['jackknife','factor_core','global_local']:raise ValueError('Restricted permutation requires audited label-blind discovery')
        from statistics31 import stratified_omnibus
        for factor,kind in traits.items():
            if kind!='categorical':raise ValueError('Restricted permutation currently supports categorical traits only')
            restricted[factor]=stratified_omnibus(scores,meta,factor,cfg['permutation_strata'])
    rows=[];means=[];diagnostic=[]
    for module in scores:
        y=scores[module]
        for factor,kind in traits.items():
            ok=y.notna();m=meta.loc[ok];values=y[ok]
            row={'module':module,'kind':'omnibus','factor':factor,'a':None,'b':None,'effect':np.nan,'p_value':np.nan,'status':'NA','samples':int(ok.sum()),'biological_units':int(m.biological_sample_id.nunique())}
            try:
                if kind=='categorical' and set(m[factor].astype(str))!=set(meta[factor].astype(str)):
                    raise ValueError('An entire requested condition is missing after score coverage exclusion')
                inference=restricted[factor][module] if factor in restricted else omnibus(values,m,factor,kind,covariates,label_blind=contract['discovery_identity']['method'] in ['jackknife','factor_core','global_local'])
                row.update(inference,status='estimated')
            except (ValueError,np.linalg.LinAlgError) as exc:row['reason']=str(exc)
            rows.append(row)
            if kind=='categorical':
                frame=pd.DataFrame({'score':y,'level':meta[factor],'unit':meta.biological_sample_id})
                for level,g in frame.groupby('level'):
                    v=g.score.dropna();means.append({'module':module,'factor':factor,'level':level,'samples':len(v),'units':g.loc[g.score.notna(),'unit'].nunique(),'mean_score':v.mean(),'sd_score':v.std()})
        # Descriptive between-condition fraction on observed frozen PC1 scores.
        frame=pd.DataFrame({'score':y,'condition':meta.condition}).dropna();total=((frame.score-frame.score.mean())**2).sum()
        between=sum(len(g)*(g.score.mean()-frame.score.mean())**2 for _,g in frame.groupby('condition'))
        contribution=(frame.score-frame.score.mean())**2
        by_unit=contribution.groupby(meta.loc[frame.index,'biological_sample_id']).sum()
        largest_unit=str(by_unit.idxmax()) if len(by_unit) and total>1e-12 else None
        unit_share=float(by_unit.max()/total) if total>1e-12 else np.nan
        diagnostic.append({'module':module,'condition_mean_variance_fraction':between/total if total>1e-12 else np.nan,'samples':len(frame),
            'maximum_unit_variance_fraction':unit_share,'highest_contribution_unit':largest_unit,'single_unit_dominated':bool(unit_share>=.5) if np.isfinite(unit_share) else None,
            'interpretation':'Descriptive pooled co-response; not within-condition coregulation/rewiring'})
        for c in contrasts:
            row={'module':module,'kind':'contrast',**c,'effect':np.nan,'p_value':np.nan,'status':'NA'}
            try:
                selected_meta = meta
                if 'sample_ids' in c:
                    ids = c['sample_ids']
                    if (not isinstance(ids, list) or len(ids) != len(set(ids))
                            or not set(ids) <= set(meta.index)):
                        raise ValueError('Contrast sample_ids must name unique observed samples')
                    selected_meta = meta.loc[ids]
                row.update(contrast(y.reindex(selected_meta.index),selected_meta,c['factor'],c['a'],c['b'],covariates,paired=c.get('paired',False)))
                row['status']='estimated' if np.isfinite(row['p_value']) else 'NA'
            except (ValueError,np.linalg.LinAlgError) as exc:row['reason']=str(exc)
            rows.append(row)
    results=adjust(rows)
    results.to_csv(out/'associations.tsv',sep='\t',index=False)
    pd.DataFrame(means).to_csv(out/'condition_scores.tsv',sep='\t',index=False)
    pd.DataFrame(diagnostic).to_csv(out/'co_response_diagnostic.tsv',sep='\t',index=False)
    scores.to_csv(out/'sample_scores.tsv',sep='\t')
    meta.to_csv(out/'analysis_samples.tsv',sep='\t')
    summary={'network_id':contract['network_id'],'discovery_id':contract['discovery_id'],'study_id':contract['study_id'],
        'samples':len(meta),'biological_units':int(meta.biological_sample_id.nunique()),'conditions':contract['conditions'],
        'modules':len(scores.columns),'contrasts':contrasts,'traits':traits,'covariates':covariates,
        'tests_requested':len(results),'tests_estimable':int(results.p_value.notna().sum()) if 'p_value' in results else 0,
        'q05':int(results.q_value.le(.05).sum()),'network_refitted':False,'frozen_assignment_sha256':contract['assignment_sha256'],
        'single_unit_dominated_modules':sum(d['single_unit_dominated'] is True for d in diagnostic),
        'unit_dominance_rule':'One biological unit contributes >=50% of centered frozen PC1 sum of squares; descriptive leverage flag, not automatic exclusion or a P-value',
        'multiplicity':'BH over all estimable module × omnibus/contrast tests requested together in this interpretation; not across separate interpretations',
        'inference_revision':'2.7: independent categorical one-way designs use 9999 condition-label permutations +1; null requires exchangeable distributions, not an unequal-variance equal-means null. Welch is a secondary diagnostic. Independent unadjusted contrasts use Welch t; categorical groups require >=3 units after coverage exclusion. Adjusted/cluster models retain explicitly exploratory robust OLS.',
        'sort_order':'FDR ascending, P ascending, stable module ID; NA last',
        'score_definition':'Frozen global standardized-event PC1; sign aligned to mean standardized PSI; score SD=1 in full discovery cohort; sample module coverage >=80%',
        'limits':['Exploratory same-study inference; estimated modules/scores treated as fixed, no independent biological validation',
            'Cluster covariance for repeated-unit omnibus is asymptotic; small biological N limits inference',
            'Paired contrasts use within-unit differences; partial pairing is NA, not silently independent',
            'Pooled condition co-response does not establish within-condition coregulation or rewiring',
            'Covariates adjust module associations only; discovery input is not automatically batch corrected'],
        'metadata_sha256':sha(out/'analysis_samples.tsv'),'interpretation_source_sha256':sha(__file__)}
    if restricted:
        summary['permutation_strata']=cfg['permutation_strata']
        summary['inference_revision']='3.1: 9999 condition-label permutations within saturated categorical nuisance strata; partial-F statistic, +1 P values; conditional distribution exchangeability required. Observational, not causal.'
        summary['stratified_inference_source_sha256']=sha(Path(__file__).with_name('statistics31.py'))
    from module_qualification import qualification
    candidates,qualification_summary=qualification(run,results,pd.DataFrame(diagnostic))
    summary['qualified_candidates']=qualification_summary['qualified']
    summary['qualification']=qualification_summary
    verify_frozen(run)
    write_json(out/'summary.json',summary);write_json(out/'config.json',cfg)
    from study_report import build
    build(run,out,out)
    with zipfile.ZipFile(out/'results.zip','w',zipfile.ZIP_DEFLATED) as z:
        for p in out.rglob('*'):
            if p.is_file() and p.name!='results.zip':z.write(p,p.relative_to(out))
    return summary

if __name__=='__main__':
    print(json.dumps(analyze(Path(sys.argv[1]),json.loads(Path(sys.argv[2]).read_text(encoding='utf-8')),Path(sys.argv[3])),ensure_ascii=False))
