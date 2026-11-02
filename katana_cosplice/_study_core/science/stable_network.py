"""Sparse all-event consensus graph with native Leiden partitions.

Inspired by resampling consensus clustering, not an invented causal model.
Internal coassignment builds the model; external refits assess its stability.
"""
from pathlib import Path
import json,time
import numpy as np
import pandas as pd
from scipy.stats import rankdata
import igraph as ig
import leidenalg as la
from cosplice_research import write_json

DEFAULTS={'neighbors':30,'inner_resamples':24,'subsample_fraction':.8,'consensus_threshold':.6,'resolution':.05,'partition':'cpm','min_module_size':20,'seed':22031}

def standardized(values):
    x=np.asarray(values,float).copy();mean=np.nanmean(x,axis=1)
    x=np.where(np.isfinite(x),x,np.nan_to_num(mean,nan=.5)[:,None])
    x=rankdata(x,axis=1).astype(np.float32);x-=x.mean(axis=1,keepdims=True)
    x/=np.maximum(np.linalg.norm(x,axis=1,keepdims=True),1e-12)
    return x

def candidates(z,k):
    p=len(z);k=min(k,p-1);pairs=[]
    for start in range(0,p,256):
        end=min(start+256,p);r=z[start:end]@z.T
        r[np.arange(end-start),np.arange(start,end)]=-np.inf
        ids=np.argpartition(-r,k-1,axis=1)[:,:k]
        a=np.broadcast_to(np.arange(start,end)[:,None],ids.shape);good=np.take_along_axis(r,ids,axis=1)>0
        a=a[good];b=ids[good];pairs.extend((np.minimum(a,b)*p+np.maximum(a,b)).tolist())
    key=np.unique(pairs);return (key//p).astype(int),(key%p).astype(int)

def labels_from_edges(n,a,b,w,*,seed=22031,resolution=1.,minimum=20,partition='rb'):
    keep=np.isfinite(w)&(w>0);g=ig.Graph(n=n,edges=list(zip(a[keep].tolist(),b[keep].tolist())),directed=False)
    labels=np.zeros(n,dtype=int)
    if g.ecount():
        fitted=la.find_partition(g,la.CPMVertexPartition if partition=='cpm' else la.RBConfigurationVertexPartition,weights=w[keep].astype(float).tolist(),resolution_parameter=resolution,seed=seed,n_iterations=4)
        groups=sorted((sorted(x) for x in fitted if len(x)>=minimum),key=lambda x:(-len(x),x[0]))
        for j,group in enumerate(groups,1):labels[group]=j
    return labels

def draw_columns(meta,rng,fraction):
    groups=meta.groupby('biological_sample_id',sort=False).indices;strata={}
    for unit,cols in groups.items():strata.setdefault(tuple(sorted(meta.iloc[cols].condition)),[]).append(unit)
    chosen=[]
    for units in strata.values():chosen.extend(rng.choice(units,max(1,int(np.ceil(len(units)*fraction))),replace=False))
    return np.concatenate([groups[u] for u in chosen]),list(map(str,chosen))

def build(values,meta,*,consensus=True,settings=None):
    cfg=dict(DEFAULTS,**(settings or {}));z=standardized(values);p=len(z)
    a,b=candidates(z,cfg['neighbors']);r=np.einsum('ij,ij->i',z[a],z[b]);base_w=np.maximum(r,0)**2
    initial=labels_from_edges(p,a,b,base_w,seed=cfg['seed'],resolution=cfg['resolution'],minimum=cfg['min_module_size'],partition=cfg['partition'])
    if not consensus:return initial,a,b,base_w,{'parameters':cfg,'internal_support':None,'draws':[]}
    rng=np.random.default_rng(cfg['seed']);support=np.zeros(len(a));weights=np.zeros(len(a));draws=[];partitions=[]
    for i in range(cfg['inner_resamples']):
        cols,units=draw_columns(meta,rng,cfg['subsample_fraction']);sub=standardized(values[:,cols]);rr=np.einsum('ij,ij->i',sub[a],sub[b]);w=np.maximum(rr,0)**2
        lab=labels_from_edges(p,a,b,w,seed=cfg['seed']+i+1,resolution=cfg['resolution'],minimum=cfg['min_module_size'],partition=cfg['partition'])
        support+=(lab[a]==lab[b])&(lab[a]!=0);weights+=w;partitions.append(lab);draws.append(units)
    support/=cfg['inner_resamples'];weights/=cfg['inner_resamples'];weights*=support
    weights[support<cfg['consensus_threshold']]=0
    labels=labels_from_edges(p,a,b,weights,seed=cfg['seed'],resolution=cfg['resolution'],minimum=cfg['min_module_size'],partition=cfg['partition'])
    return labels,a,b,weights,{'parameters':cfg,'internal_support':support,'partitions':np.asarray(partitions),'draws':draws}

def fit_stable(measured,complete,meta,out,consensus=True,settings=None):
    out=Path(out);out.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    write_json(out/'status.json',{'status':'running','method':'consensus' if consensus else 'leiden'})
    labels,a,b,w,info=build(measured.to_numpy(),meta,consensus=consensus,settings=settings)
    keys=np.asarray(complete.index,str);good=w>0
    edges=pd.DataFrame({'source':keys[a[good]],'target':keys[b[good]],'weight':w[good],'sign':1})
    if consensus:edges['construction_coassignment']=info['internal_support'][good]
    edges.to_csv(out/'edges.tsv',sep='\t',index=False)
    pd.DataFrame({'event_key':keys,'module':labels}).to_csv(out/'assignments.tsv',sep='\t',index=False)
    complete.to_csv(out/'input.tsv',sep='\t',index_label='event_key');meta.to_csv(out/'samples.tsv',sep='\t',index=False)
    z=(complete.T-complete.mean(axis=1))/complete.std(axis=1).clip(lower=1e-8)
    activity={f'M{m}':z.loc[:,keys[labels==m]].mean(axis=1) for m in np.unique(labels) if m}
    pd.DataFrame(activity,index=complete.columns).to_csv(out/'scores.tsv',sep='\t',index_label='sample_id')
    write_json(out/'parameters.json',info['parameters']);write_json(out/'internal_draws.json',info['draws'])
    if consensus:np.savez_compressed(out/'internal_partitions.npz',labels=info['partitions'],event_key=keys)
    status={'status':'completed','method':'consensus' if consensus else 'leiden','events':len(keys),'samples':len(meta),
        'assigned_events':int((labels!=0).sum()),'modules':len(activity),'graph_edges':int(good.sum()),'seconds':time.monotonic()-started,
        'implementation':'native Leiden on sparse positive Spearman graph; biological-unit subsample coassignment consensus' if consensus else 'native Leiden on identical full-data sparse positive Spearman graph',
        'support_scope':'Internal coassignment constructs the graph and is NOT independent stability evidence','parameters':info['parameters']}
    write_json(out/'status.json',status);return status
