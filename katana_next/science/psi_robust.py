"""Explicit shared robustness policies; these do not change the WGCNA algorithm."""
from __future__ import annotations
import hashlib
import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold

SEED = 20261018
ALPHAS = (1000., 100., 10., 1., .1)
_CALIBRATION_MASKS={}


def coverage_variable_audit(matrix, filters, *, max_features=2000):
    from katana_next.science.psi_network import psi_filter_audit
    if not isinstance(max_features, int) or max_features < 40:
        raise ValueError("max_features must be an integer >=40")
    audit = psi_filter_audit(matrix, filters)
    for col in list(audit):
        if col.startswith('failed_') or col in {'passed', 'filter_reasons'}:
            audit['legacy_'+col] = audit[col]
    audit['failed_valid_fraction'] = ~audit.valid_fraction.ge(filters.min_valid_fraction)
    audit['failed_variance'] = ~audit.psi_variance.gt(1e-12)
    audit['failed_range'] = False
    audit['failed_alternative_fraction'] = False
    audit['threshold_variance'] = 1e-12
    audit['threshold_range'] = 0.
    audit['threshold_alternative_fraction'] = 0.
    eligible = ~(audit.failed_valid_fraction | audit.failed_variance)
    ranked = audit.loc[eligible].assign(event_id=audit.index[eligible].astype(str)).sort_values(
        ['psi_variance','event_id'], ascending=[False,True], kind='stable')
    selected = ranked.index[:max_features]
    audit['failed_cap'] = eligible & ~audit.index.isin(selected)
    flags = ['valid_fraction', 'variance', 'cap']
    audit['failure_count'] = audit[['failed_'+x for x in flags]].sum(axis=1)
    audit['filter_reasons'] = [';'.join(x for x in flags if row['failed_'+x]) or 'passed'
                               for _,row in audit.iterrows()]
    audit['passed'] = audit.failure_count.eq(0)
    audit['filter_profile'] = 'coverage_variable'
    audit['max_features'] = max_features
    audit['variance_rule'] = 'strictly_greater_than_numeric_zero_1e-12'
    return audit


def calibrate_guarded_projection(entry, measured):
    """Use only fully measured training vectors; no heldout observations/traits."""
    keys=entry['events'];p=measured.reindex(keys)
    p=p.loc[:,p.notna().all(axis=0)]
    default=dict(weight=0.,status='insufficient_fully_measured_training_vectors',training_vectors=p.shape[1],
                 seed=SEED,weights=[0.,.25,.5,.75,1.],selection='masked_training_score_fidelity_one_SE_prefer_mean')
    if p.shape[1]<15 or len(keys)<2:return default
    a=p.to_numpy(float);load=np.array(entry['loading']);sing=entry['singular_value']
    center=np.array(entry['centers'])[:,None];scale=np.array(entry['scales'])[:,None]
    z=(a-center)/scale;truth=z.T@load/sing;loss=[]
    weights=np.array(default['weights']);all_samples=tuple(sorted(measured.columns.astype(str)))
    positions=pd.Index(all_samples).get_indexer(p.columns.astype(str))
    for rate in [.05,.1,.2]:
        for repeat in range(8):
            rows=[]
            for key in keys:
                cachekey=(str(key),all_samples,rate,repeat)
                if cachekey not in _CALIBRATION_MASKS:
                    seed=int(hashlib.sha256(f'{SEED}|{rate}|{repeat}|{key}'.encode()).hexdigest()[:12],16)
                    _CALIBRATION_MASKS[cachekey]=np.random.default_rng(seed).random(len(all_samples))>=rate
                rows.append(_CALIBRATION_MASKS[cachekey][positions])
            seen=np.array(rows)
            observed=seen.mean(axis=0);info=seen.T@(load**2)/(load@load)
            eligible=(observed>=.8)&(info>=.8-1e-12)
            if eligible.sum()<3:continue
            old=np.where(seen,z,0).T@load/sing
            ls=np.divide(old,info,out=old.copy(),where=info>1e-12)
            choices=old[:,None]+(ls-old)[:,None]*weights
            loss.append(((choices[eligible]-truth[eligible,None])**2).mean(axis=0))
    if len(loss)<5:
        default['status']='insufficient_eligible_training_masks';return default
    errors=np.array(loss);mean=errors.mean(axis=0);best=int(mean.argmin())
    bound=mean[best]+errors[:,best].std(ddof=1)/np.sqrt(len(errors))
    chosen=int(np.flatnonzero(mean<=bound+1e-15)[0])
    return dict(**{**default,'weight':float(weights[chosen]),'status':'training_calibrated'},
                mask_repeats=len(errors),inner_mean_MSE=mean.tolist(),one_SE_bound=float(bound))


def _design(scores, reference=None):
    if not scores.index.is_unique or not scores.columns.is_unique:
        raise ValueError('Scores require unique sample and module IDs')
    columns=[str(c) for c in scores.columns]
    if any(c.startswith('REF::') for c in columns):raise ValueError('REF:: is reserved for references')
    x=scores.copy();x.columns=columns
    if reference is not None:
        if not reference.index.is_unique or not reference.columns.is_unique or set(reference.index)!=set(x.index):
            raise ValueError('Reference scores require exactly matching unique sample IDs')
        ref=reference.loc[x.index].copy();ref.columns=['REF::'+str(c) for c in ref.columns]
        x=pd.concat([x,ref],axis=1)
    return x


def _fit(kind, alpha, frame, y):
    modules=[c for c in frame if not c.startswith('REF::')]
    refs=[c for c in frame if c.startswith('REF::')]
    cols=[]
    if kind=='single' and modules:
        correlations=[]
        for c in sorted(modules):
            v=frame[c].to_numpy(float)
            corr=abs(float(np.corrcoef(v,y)[0,1])) if np.std(v)>1e-12 and np.std(y)>1e-12 else -1.
            correlations.append((corr if np.isfinite(corr) else -1.,c))
        best=sorted(correlations,key=lambda z:(-z[0],z[1]))[0]
        if best[0]>=0:cols=[best[1]]
    elif kind=='modules':cols=modules
    elif kind=='reference':cols=refs
    elif kind=='combined':cols=modules+refs
    x=frame[cols].to_numpy(float);center=x.mean(axis=0);scale=x.std(axis=0)
    scale=np.where(scale<1e-12,1.,scale);coef=np.zeros(len(cols));intercept=float(np.mean(y))
    if cols:
        z=(x-center)/scale
        if kind=='single':
            b=np.linalg.lstsq(np.column_stack([np.ones(len(y)),z]),y,rcond=None)[0]
            intercept=float(b[0]);coef=b[1:]
        else:
            fit=Ridge(alpha=alpha).fit(z,y);coef=fit.coef_;intercept=float(fit.intercept_)
    return dict(kind=kind,alpha=alpha,features=cols,centers=center.tolist(),scales=scale.tolist(),
                coefficients=coef.tolist(),intercept=intercept)


def _predict(model, frame):
    coef=np.array(model['coefficients']);used=np.abs(coef)>1e-15
    if not used.any():return np.full(len(frame),model['intercept']),np.ones(len(frame),dtype=bool)
    x=frame.reindex(columns=model['features']).to_numpy(float)
    z=(x[:,used]-np.array(model['centers'])[used])/np.array(model['scales'])[used]
    valid=np.isfinite(z).all(axis=1)
    values=model['intercept']+z@coef[used]
    return np.where(valid,values,np.nan),valid


def fit_reliable_predictors(scores, traits, reference_scores=None):
    if not traits.index.is_unique or set(scores.index)!=set(traits.index):
        raise ValueError('Traits must exactly match unique training sample IDs')
    frame=_design(scores,reference_scores)
    if not np.isfinite(frame.to_numpy(float)).all():raise ValueError('Training predictors must be finite')
    modules=list(scores);refs=[] if reference_scores is None else list(reference_scores)
    candidates=[('intercept',None)]
    if refs:candidates += [('reference',a) for a in ALPHAS]
    if modules:candidates += [('single',None)]+[('modules',a) for a in ALPHAS]
    if modules and refs:candidates += [('combined',a) for a in ALPHAS]
    entries=[]
    for trait in traits:
        y=pd.to_numeric(traits.loc[frame.index,trait],errors='raise').to_numpy(float)
        if np.isinf(y).any():raise ValueError('Traits cannot contain infinity')
        good=np.isfinite(y);yy=y[good];xx=frame.loc[good]
        if len(yy)<15 or np.std(yy)<1e-12:continue
        errors=np.zeros((5,len(candidates)))
        for fold,(tr,va) in enumerate(KFold(5,shuffle=True,random_state=SEED).split(yy)):
            for col,(kind,alpha) in enumerate(candidates):
                model=_fit(kind,alpha,xx.iloc[tr],yy[tr]);pred,_=_predict(model,xx.iloc[va])
                errors[fold,col]=np.mean((pred-yy[va])**2)
        mean=errors.mean(axis=0);se=errors.std(axis=0,ddof=1)/np.sqrt(5)
        best=int(mean.argmin());bound=mean[best]+se[best]
        selected=int(np.flatnonzero(mean<=bound+1e-15)[0])
        hierarchy=[(selected,'primary')]
        ref_ids=[i for i,c in enumerate(candidates) if c[0] in {'intercept','reference'}]
        rb=min(ref_ids,key=lambda i:mean[i]);rbound=mean[rb]+se[rb]
        reference=next(i for i in ref_ids if mean[i]<=rbound+1e-15)
        if reference!=selected and candidates[reference][0]!='intercept':hierarchy.append((reference,'reference_fallback'))
        mod_ids=[i for i,c in enumerate(candidates) if c[0] in {'single','modules'}]
        if mod_ids:
            mb=min(mod_ids,key=lambda i:mean[i]);mbound=mean[mb]+se[mb]
            chosen=next(i for i in mod_ids if mean[i]<=mbound+1e-15)
            if chosen not in [i for i,_ in hierarchy] and mean[chosen]+se[chosen]<mean[0]:
                hierarchy.append((chosen,'validated_module_fallback'))
        if 0 not in [i for i,_ in hierarchy]:hierarchy.append((0,'training_mean_fallback'))
        models=[dict(route=route,model=_fit(*candidates[i],xx,yy),candidate=i) for i,route in hierarchy]
        entries.append(dict(trait=str(trait),training_mean=float(yy.mean()),training_n=len(yy),
            hierarchy=models,module_columns=[str(c) for c in modules],reference_columns=[str(c) for c in refs],
            selection='five_training_folds_one_SE_prefer_baseline_reference_single_then_ridge_conditional_on_frozen_modules',
            candidates=[dict(kind=k,alpha=a) for k,a in candidates],inner_mean_MSE=mean.tolist(),
            inner_selected_MSE=float(mean[selected]),inner_one_SE_bound=float(bound)))
    return entries


def predict_reliable(entry, scores, reference_scores=None):
    frame=_design(scores,reference_scores)
    result=np.full(len(frame),np.nan);audit=pd.DataFrame(index=frame.index)
    audit['prediction_route']='unavailable';audit['uses_module']=False;audit['uses_reference']=False
    audit['primary_available']=False
    audit['available_module_count']=scores.reindex(columns=entry['module_columns']).notna().sum(axis=1)
    for i,item in enumerate(entry['hierarchy']):
        model=item['model'];values,valid=_predict(model,frame)
        if i==0:audit['primary_available']=valid
        use=~np.isfinite(result)&valid
        result[use]=values[use]
        nonzero=[f for f,c in zip(model['features'],model['coefficients']) if abs(c)>1e-15]
        audit.loc[use,'prediction_route']=item['route']+':'+model['kind']
        audit.loc[use,'uses_module']=any(not c.startswith('REF::') for c in nonzero)
        audit.loc[use,'uses_reference']=any(c.startswith('REF::') for c in nonzero)
    assert np.isfinite(result).all(),'Every fitted hierarchy must include a finite intercept fallback'
    return result,audit
