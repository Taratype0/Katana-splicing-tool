"""Unit influence QC and leave-one-biological-unit-out Spearman graph.

Influence screening is label-free (unit identities, never condition effects).
Each edge must survive deletion of *every* unit, including all repeated visits.
This is a jackknife robustness guard, not evidence of causality or a new theorem.
"""
import time
from pathlib import Path
import numpy as np
import pandas as pd
from stable_network import standardized,candidates,labels_from_edges,DEFAULTS
from cosplice_research import write_json

QC={'observed_fraction':.8,'variance_gt':1e-12,'maximum_unit_fraction':.5,'minimum_cell_reads':10}

def unit_audit(measured,meta,read_support=None):
    meta=meta.set_index('sample_id').loc[measured.columns];x=measured.copy()
    if read_support is not None:
        if not x.index.equals(read_support.index) or not x.columns.equals(read_support.columns):raise ValueError('Read support identities/order differ')
        v=read_support.to_numpy();obs=v[np.isfinite(v)]
        if np.isinf(v).any() or (obs<0).any():raise ValueError('Read support must be finite nonnegative effective reads or NA')
        x=x.where(read_support.ge(QC['minimum_cell_reads']))
    # Every subject has total weight 1, so repeat visits cannot manufacture units.
    weights=1/meta.biological_sample_id.map(meta.biological_sample_id.value_counts()).to_numpy()
    a=x.to_numpy();ok=np.isfinite(a);den=(ok*weights).sum(axis=1)
    mean=np.divide(np.nansum(a*weights,axis=1),den,out=np.full(len(a),np.nan),where=den>0)
    squared=np.nan_to_num((a-mean[:,None])**2,nan=0)*weights
    units=meta.groupby('biological_sample_id',sort=True).indices
    per_unit=np.stack([squared[:,idx].sum(axis=1) for idx in units.values()],axis=1)
    total=per_unit.sum(axis=1);share=np.divide(per_unit.max(axis=1),total,out=np.full(len(a),np.nan),where=total>1e-12)
    audit=pd.DataFrame({'observed_fraction':ok.mean(axis=1),'variance':x.var(axis=1),'maximum_unit_variance_fraction':share,
        'highest_contribution_unit':np.asarray(list(units))[per_unit.argmax(axis=1)],
        'masked_low_support_cells':measured.notna().sum(axis=1)-x.notna().sum(axis=1)},index=x.index)
    audit['eligible']=audit.observed_fraction.ge(.8)&audit.variance.gt(1e-12)&audit.maximum_unit_variance_fraction.lt(.5)
    audit['reason']=[';'.join(r for failed,r in [(o<.8,'missingness_or_read_support'),(not v>1e-12,'constant'),(not s<.5,'unit_dominated_variation')] if failed) or 'passed' for o,v,s in zip(audit.observed_fraction,audit.variance,audit.maximum_unit_variance_fraction)]
    return x,audit

def edge_minimum(values,meta,a,b):
    z=standardized(values);observed=np.einsum('ij,ij->i',z[a],z[b]).clip(-1,1)
    minimum=observed.copy();groups=meta.groupby('biological_sample_id',sort=True).indices
    trace=[]
    for unit,cols in groups.items():
        keep=np.ones(values.shape[1],dtype=bool);keep[cols]=False
        sub=standardized(values[:,keep]);corr=np.einsum('ij,ij->i',sub[a],sub[b]).clip(-1,1)
        minimum=np.minimum(minimum,corr)
        trace.append({'unit':str(unit),'removed_observations':len(cols),'edges_nonpositive_or_unestimable':int((corr<=0).sum())})
    return observed,minimum,trace

def linear_standardized(values):
    x=np.asarray(values,float).copy();mean=np.nanmean(x,axis=1)
    x=np.where(np.isfinite(x),x,mean[:,None]);x-=x.mean(axis=1,keepdims=True)
    x/=np.maximum(np.linalg.norm(x,axis=1,keepdims=True),1e-12)
    return x

def linear_minimum(values,meta,a,b):
    z=linear_standardized(values);minimum=np.einsum('ij,ij->i',z[a],z[b]).clip(-1,1)
    for cols in meta.groupby('biological_sample_id',sort=True).indices.values():
        keep=np.ones(values.shape[1],bool);keep[cols]=False;z=linear_standardized(values[:,keep])
        minimum=np.minimum(minimum,np.einsum('ij,ij->i',z[a],z[b]).clip(-1,1))
    return minimum

def fit_jackknife(measured,complete,meta,out,linear_guard=True):
    out=Path(out);out.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    z=standardized(measured.to_numpy());a,b=candidates(z,DEFAULTS['neighbors'])
    raw,minimum,trace=edge_minimum(measured.to_numpy(),meta,a,b)
    linear=linear_minimum(measured.to_numpy(),meta,a,b) if linear_guard else minimum
    weights=np.maximum(np.minimum(minimum,linear),0)**2
    labels=labels_from_edges(len(measured),a,b,weights,seed=DEFAULTS['seed'],resolution=DEFAULTS['resolution'],minimum=DEFAULTS['min_module_size'],partition='cpm')
    keys=np.asarray(measured.index,str);good=weights>0
    pd.DataFrame({'source':keys[a[good]],'target':keys[b[good]],'weight':weights[good],'sign':1,'full_spearman':raw[good],'minimum_leave_unit_out_spearman':minimum[good],'minimum_leave_unit_out_pearson':linear[good]}).to_csv(out/'edges.tsv',sep='\t',index=False)
    pd.DataFrame({'event_key':keys,'module':labels}).to_csv(out/'assignments.tsv',sep='\t',index=False)
    complete.to_csv(out/'input.tsv',sep='\t',index_label='event_key');meta.to_csv(out/'samples.tsv',sep='\t',index=False)
    write_json(out/'unit_deletion_audit.json',trace)
    status={'status':'completed','method':'jackknife','events':len(keys),'samples':len(meta),'units':int(meta.biological_sample_id.nunique()),
        'assigned_events':int((labels!=0).sum()),'modules':len(set(labels)-{0}),'graph_edges':int(good.sum()),'candidate_edges':len(a),
        'seconds':time.monotonic()-started,'parameters':dict(DEFAULTS,edge_weight='max(0,min(Spearman,Pearson) over full data and all unit deletions)^2' if linear_guard else 'rank-only',linear_guard=linear_guard),
        'implementation':'Native Leiden CPM with worst-unit-deletion rank-correlation edge guard; imputation/ranks recomputed after each deletion',
        'support_scope':'Construction robustness, not independent stability validation'}
    write_json(out/'status.json',status);return status
