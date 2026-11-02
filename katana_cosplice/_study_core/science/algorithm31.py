"""Global programs plus local residual topology, with no event-count cap.

The low-rank + sparse decomposition is inspired by POET; this is not the
published POET estimator. Weighted shared-neighbour topology borrows WGCNA's
TOM formula on a sparse graph. Programs never serve as module labels.
"""
from pathlib import Path
from types import SimpleNamespace
import time
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.linalg import eigh
from stable_network import labels_from_edges
from cosplice_research import write_json

DEFAULTS=dict(global_rank=4,neighbors=30,threshold_scale=.5,resolution=.02,
              minimum_size=20,seed=31,search='auto',covariance='global_plus_sparse')

def decompose(complete,rank=4):
    x=complete.to_numpy(dtype=float,copy=True)
    mean=x.mean(axis=1);scale=x.std(axis=1,ddof=1)
    z=(x-mean[:,None])/np.maximum(scale[:,None],1e-8)
    rank=min(rank,z.shape[1]-3,z.shape[0]-1)
    # PCA gives an actual orthogonal complement, unlike posterior FA scores.
    # Solve the small sample Gram matrix exactly. Randomized SVD can leave
    # components of the global score in the purported residual complement.
    gram=z.T@z
    values,vectors=eigh(gram,subset_by_index=(len(gram)-rank,len(gram)-1))
    values=values[::-1];vectors=vectors[:,::-1]
    scores=vectors*np.sqrt(np.maximum(values,1e-12))
    components=(z@scores/np.maximum(values,1e-12)).T
    model=SimpleNamespace(components_=components,mean_=np.zeros(len(z)),
                          explained_variance_ratio_=values/np.trace(gram))
    global_fit=scores@model.components_
    residual=z-global_fit.T
    return z,residual,scores,model,mean,scale

def local_graph(residual,settings=None,global_fit=None):
    cfg=dict(DEFAULTS,**(settings or {}));p,n=residual.shape
    z=residual.copy();z-=z.mean(axis=1,keepdims=True)
    norm=np.linalg.norm(z,axis=1);sd=norm/np.sqrt(n-1)
    z/=np.maximum(norm[:,None],1e-12)
    penalty=cfg['threshold_scale']*np.sqrt(np.log(max(p,2))/max(n-cfg['global_rank']-1,1))
    g=None if global_fit is None else global_fit/np.sqrt(n-1)
    search_z=z if g is None else (residual+global_fit)/np.sqrt(n-1)
    k=min(cfg['neighbors'],p-1);aa=[];bb=[];rr=[];recall=None
    # Exact blockwise search: all events participate; 256 is a memory block,
    # not a feature selection limit. No dense p-by-p matrix is retained.
    approximate=cfg['search']=='approximate' or (cfg['search']=='auto' and p>40000)
    if approximate:
        from pynndescent import NNDescent
        index=NNDescent(search_z.astype(np.float32),n_neighbors=min(2*k+1,p),metric='euclidean',
                        random_state=cfg['seed'],n_jobs=2,low_memory=True)
        neighbours,_=index.neighbor_graph
        for start in range(0,p,1024):
            stop=min(start+1024,p);ids=neighbours[start:stop]
            valid=ids>=0;safe=np.maximum(ids,0)
            r=np.einsum('ij,ikj->ik',z[start:stop],z[safe])
            if g is not None:
                r=np.einsum('ij,ikj->ik',g[start:stop],g[safe])+np.sign(r)*np.maximum(np.abs(r)-penalty,0)*sd[start:stop,None]*sd[safe]
            r[~valid]=-np.inf;r[ids==np.arange(start,stop)[:,None]]=-np.inf
            pick=np.argpartition(-r,k-1,axis=1)[:,:k]
            aa.append(np.repeat(np.arange(start,stop),k));bb.append(np.take_along_axis(safe,pick,axis=1).ravel())
            rr.append(np.take_along_axis(r,pick,axis=1).ravel())
        check=np.random.default_rng(31031).choice(p,min(128,p),replace=False)
        exact=z[check]@z.T
        if g is not None:exact=g[check]@g.T+np.sign(exact)*np.maximum(np.abs(exact)-penalty,0)*sd[check,None]*sd[None,:]
        exact[np.arange(len(check)),check]=-np.inf
        truth=np.argpartition(-exact,k-1,axis=1)[:,:k]
        selected=np.concatenate(bb).reshape(p,k)
        recall=float(np.mean([len(set(t)&set(selected[i]))/k for i,t in zip(check,truth)]))
        if recall<.8:raise RuntimeError(f'Approximate neighbourhood recall {recall:.3f} <0.8; no validated graph produced')
    else:
        for start in range(0,p,256):
            stop=min(start+256,p);r=z[start:stop]@z.T
            if g is not None:r=g[start:stop]@g.T+np.sign(r)*np.maximum(np.abs(r)-penalty,0)*sd[start:stop,None]*sd[None,:]
            r[np.arange(stop-start),np.arange(start,stop)]=-np.inf
            ids=np.argpartition(-r,k-1,axis=1)[:,:k]
            aa.append(np.repeat(np.arange(start,stop),k));bb.append(ids.ravel())
            rr.append(np.take_along_axis(r,ids,axis=1).ravel())
    a=np.concatenate(aa);b=np.concatenate(bb);r=np.concatenate(rr)
    # Soft threshold uses sample count after global-factor degrees of freedom;
    # it is an estimation regularizer, never an edge P value or FDR.
    weight=np.maximum(r-(penalty if g is None else 0),0)**2
    adj=sparse.csr_matrix((weight,(a,b)),shape=(p,p));adj=adj.maximum(adj.T)
    adj.setdiag(0);adj.eliminate_zeros()
    # Evaluate TOM only for observed candidate edges. Sparse multiplication
    # restricted by the candidate pattern avoids materialising two-hop pairs.
    coo=sparse.triu(adj,k=1).tocoo();a=coo.row;b=coo.col;w=coo.data
    degree=np.asarray(adj.sum(axis=1)).ravel();shared=np.empty(len(a))
    for start in range(0,len(a),20000):
        sl=slice(start,start+20000)
        shared[sl]=np.asarray(adj[a[sl]].multiply(adj[b[sl]]).sum(axis=1)).ravel()
    tom=(shared+w)/(np.minimum(degree[a],degree[b])+1-w)
    return a,b,tom,dict(penalty=float(penalty),candidate_edges=len(a),
        graph='low-rank covariance plus softly thresholded residual covariance; sparse weighted topological overlap' if g is not None else 'positive residual correlation; soft threshold; sparse weighted topological overlap',
        edge_scope='estimated shared-variation association; not causal',
        search='NNDescent all-event candidates with exact correlation weights' if approximate else 'exact all-event nearest neighbours, memory-blocked',
        approximate_recall_check=recall)

def fit_global_local(measured,complete,meta,out,settings=None):
    out=Path(out);out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    cfg=dict(DEFAULTS,**(settings or {}))
    if meta.biological_sample_id.duplicated().any():
        raise ValueError('Global/local engine currently requires one observation per biological unit')
    if len(meta)<15:raise ValueError('At least 15 independent biological units required')
    if list(meta.sample_id)!=list(complete.columns):raise ValueError('Sample identities/order differ')
    if not np.isfinite(complete.to_numpy()).all():raise ValueError('Training imputation required')
    z,residual,scores,model,mean,scale=decompose(complete,cfg['global_rank'])
    a,b,w,diagnostic=local_graph(residual,cfg,z-residual if cfg['covariance']=='global_plus_sparse' else None)
    labels=labels_from_edges(len(complete),a,b,w,seed=cfg['seed'],resolution=cfg['resolution'],
                             minimum=cfg['minimum_size'],partition='cpm')
    keys=complete.index.to_numpy(str);programs=[f'P{i+1}' for i in range(scores.shape[1])]
    scorez=scores/np.maximum(scores.std(axis=0,ddof=1),1e-12)
    loading=z@scorez/(z.shape[1]-1)
    pd.DataFrame({'event_key':keys,'module':labels}).to_csv(out/'assignments.tsv',sep='\t',index=False)
    pd.DataFrame({'source':keys[a],'target':keys[b],'weight':w,'sign':1,
                  'edge_type':'global_local_topological_overlap'}).to_csv(out/'edges.tsv',sep='\t',index=False)
    pd.DataFrame(scorez,index=complete.columns,columns=programs).to_csv(out/'program_scores.tsv',sep='\t',index_label='sample_id')
    pd.DataFrame(loading,index=keys,columns=programs).to_csv(out/'program_loadings.tsv',sep='\t',index_label='event_key')
    pd.DataFrame({'program':programs,'explained_variance_ratio':model.explained_variance_ratio_,
                  'interpretation':'global shared variation; biological or technical'}).to_csv(out/'programs.tsv',sep='\t',index=False)
    links=[];details=[]
    for m in sorted(set(labels)-{0}):
        ix=np.flatnonzero(labels==m)
        _,_,v=np.linalg.svd(z[ix],full_matrices=False);s=v[0]
        if (z[ix]@s).sum()<0:s=-s
        for j,p in enumerate(programs):
            links.append(dict(module=f'M{m}',program=p,correlation=float(np.corrcoef(s,scorez[:,j])[0,1])))
        details.append(dict(module=f'M{m}',events=len(ix),
             residual_variance_fraction=float(np.sum(residual[ix]**2)/np.sum(z[ix]**2))))
    pd.DataFrame(links,columns=['module','program','correlation']).to_csv(out/'module_program_links.tsv',sep='\t',index=False)
    pd.DataFrame(details,columns=['module','events','residual_variance_fraction']).to_csv(out/'local_modules.tsv',sep='\t',index=False)
    np.savez_compressed(out/'global_model.npz',components=model.components_,mean=model.mean_,
                        event_mean=mean,event_scale=scale,event_key=keys)
    complete.to_csv(out/'input.tsv',sep='\t',index_label='event_key')
    meta.to_csv(out/'samples.tsv',sep='\t',index=False)
    status=dict(status='completed',method='global_local',events=len(keys),samples=len(meta),
                programs=len(programs),modules=len(details),assigned_events=int((labels>0).sum()),
                seconds=time.monotonic()-start,parameters=cfg,label_blind=True,
                program_module_separation=True,diagnostic=diagnostic,
                stability='Requires independent full-QC outer refits; internal edge weights are not validation')
    write_json(out/'status.json',status);return status

def deletion_draws(meta,iterations=20,seed=31071):
    """Delete 20% (at least one) per stratum; never merely 1 of 462."""
    if meta.biological_sample_id.duplicated().any():raise ValueError('Independent units required')
    rng=np.random.default_rng(seed);draws=[];seen=set();groups=list(meta.groupby('condition',sort=True))
    for attempt in range(iterations*100):
        keep=[]
        for _,g in groups:
            n=len(g);take=max(1,min(n-1,int(np.floor(.8*n)))) if n>1 else 1
            keep.extend(rng.choice(g.sample_id.to_numpy(),take,replace=False).tolist())
        chosen=set(keep);key=tuple(s for s in meta.sample_id if s in chosen)
        if key not in seen:draws.append(list(key));seen.add(key)
        if len(draws)==iterations:break
    if len(draws)<iterations:raise ValueError('Not enough distinct biological deletion splits')
    if min(map(len,draws))<15:raise ValueError('Deletion refits would have fewer than 15 independent units')
    return draws
