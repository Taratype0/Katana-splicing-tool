"""Cross-field network engines. Native source retained; adaptations named explicitly."""
from pathlib import Path
import json,sys,time
import numpy as np
import pandas as pd
from scipy.linalg import solve
import networkx as nx
from cosplice_research import ROOT,write_json

sys.path.insert(0,str(ROOT/'python_packages'))
sys.path.insert(0,str(ROOT/'vendor/arboreto'))

class OASNetwork:
    """Exact OAS precision in sample space; O(np+n²) stored, streamed edge rows."""
    def __init__(self,values):
        self.mean=values.mean(axis=1);self.scale=np.maximum(values.std(axis=1),1e-8)
        self.x=((values-self.mean[:,None])/self.scale[:,None]).T
        n,p=self.x.shape;g=self.x@self.x.T
        mu=np.sum(self.x*self.x)/(n*p)
        alpha=np.sum(g*g)/(n*n*p*p)
        denominator=(n+1)*(alpha-mu*mu/p)
        self.shrinkage=1. if denominator<=0 else min((alpha+mu*mu)/denominator,1.)
        self.d=max(self.shrinkage*mu,1e-10);self.a=(1-self.shrinkage)/n
        self.bx=solve(np.eye(n)+(self.a/self.d)*g,self.x,assume_a='pos')
        self.diag=1/self.d-(self.a/self.d**2)*np.sum(self.x*self.bx,axis=0)
    def rows(self,start,end):
        r=(self.a/self.d**2)*(self.x[:,start:end].T@self.bx)
        r/=np.sqrt(self.diag[start:end,None]*self.diag[None,:])
        r[np.arange(end-start),np.arange(start,end)]=0
        return r
    def dense(self):return self.rows(0,self.x.shape[1])

def nearest_edges(row_source,keys,k=20,genes=None):
    rows=[];p=len(keys)
    for start in range(0,p,128):
        block=row_source(start,min(p,start+128))
        for offset,r in enumerate(block):
            i=start+offset;valid=np.arange(p)!=i
            if genes is not None:valid&=np.asarray(genes)!=genes[i]
            ids=np.flatnonzero(valid&(np.abs(r)>0));ids=ids[np.argsort(-np.abs(r[ids]),kind='stable')[:k]]
            rows.extend((str(keys[i]),str(keys[j]),float(abs(r[j])),float(np.sign(r[j]))) for j in ids)
    return pd.DataFrame(rows,columns=['source','target','weight','sign'])

def graph_outputs(edges,complete,out,mutual=False):
    keys=list(map(str,complete.index));positions={k:i for i,k in enumerate(keys)}
    g=nx.Graph();g.add_nodes_from(keys);pairs=set(zip(edges.source,edges.target))
    for r in edges.itertuples():
        if r.source not in positions or r.target not in positions:continue
        if mutual and (r.target,r.source) not in pairs:continue
        if r.source!=r.target:g.add_edge(r.source,r.target,weight=max(r.weight,g.get_edge_data(r.source,r.target,{}).get('weight',0)))
    labels={k:0 for k in keys}
    if g.number_of_edges():
        groups=nx.community.louvain_communities(g,seed=31,weight='weight')
        for module,group in enumerate(sorted(groups,key=lambda z:(-len(z),min(z))),1):
            if len(group)>=5:
                for key in group:labels[key]=module
    pd.DataFrame({'event_key':keys,'module':[labels[k] for k in keys]}).to_csv(out/'assignments.tsv',sep='\t',index=False)
    # Mean standardized module activity, explicitly distinct from WGCNA eigengene.
    z=(complete.T-complete.mean(axis=1))/complete.std(axis=1).clip(lower=1e-8)
    activity={f'M{m}':z[[k for k in keys if labels[k]==m]].mean(axis=1) for m in sorted(set(labels.values())) if m}
    pd.DataFrame(activity,index=complete.columns).to_csv(out/'scores.tsv',sep='\t',index_label='sample_id')
    return {'modules':len(activity),'graph_edges':g.number_of_edges(),'assigned_events':sum(v!=0 for v in labels.values())}

def fit_advanced(method,measured,complete,meta,out,expression=None):
    out=Path(out);out.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    complete.to_csv(out/'input.tsv',sep='\t',index_label='event_key');meta.to_csv(out/'samples.tsv',sep='\t',index=False)
    status={'method':method,'events':len(complete),'samples':complete.shape[1],'status':'running'};write_json(out/'status.json',status)
    try:
        if method=='oas':
            model=OASNetwork(complete.to_numpy());edges=nearest_edges(model.rows,complete.index)
            np.savez_compressed(out/'oas.npz',mean=model.mean,scale=model.scale,x=model.x,bx=model.bx,diag=model.diag,a=model.a,d=model.d)
            status.update(shrinkage=model.shrinkage,implementation='exact OAS sample-space precision; mutual top20 graph')
            status.update(graph_outputs(edges,complete,out,mutual=True))
        elif method=='grnboost2':
            from arboreto.core import fit_model,to_links_df,SGBM_KWARGS,EARLY_STOP_WINDOW_LENGTH
            values=complete.to_numpy().T;keys=list(map(str,complete.index));parts=[];diags=[]
            # Scheduling differs from Dask; native fitting/importance/early stopping unchanged.
            for i,key in enumerate(keys):
                ids=np.flatnonzero(np.arange(len(keys))!=i)
                model=fit_model('GBM',SGBM_KWARGS,values[:,ids],values[:,i],EARLY_STOP_WINDOW_LENGTH,31)
                tab=to_links_df('GBM',SGBM_KWARGS,model,[keys[j] for j in ids],key)
                # Store all nonzero inferred edges. No target/event exclusion.
                parts.append(tab);diags.append({'target':key,'trees':model.n_estimators_})
                if i%25==0:print(f'GRNBoost2 target {i+1}/{len(keys)}',flush=True)
            links=pd.concat(parts,ignore_index=True);links.to_csv(out/'native_importance.tsv',sep='\t',index=False)
            pd.DataFrame(diags).to_csv(out/'trees.tsv',sep='\t',index=False)
            # Edge rows use source=predictor,target=response, never a causal assertion.
            edges=links.rename(columns={'TF':'source','importance':'weight'});edges['sign']=np.nan
            status.update(graph_outputs(edges.sort_values('weight',ascending=False).groupby('target',sort=False).head(20),complete,out))
            status['implementation']='original Arboreto GRNBoost2 core; serial target scheduling; predictive importance'
        elif method=='mofa':
            from mofapy2.run.entry_point import entry_point
            view=measured.to_numpy().T;names=['usage'];features=[list(map(str,measured.index))];data=[[view]]
            if expression is not None:
                if set(expression.columns)!=set(measured.columns):raise ValueError('Expression sample IDs must exactly match usage')
                e=expression.loc[:,measured.columns]
                e=e.loc[(e.notna().mean(axis=1)>=.8)&e.var(axis=1).gt(1e-12)]
                if len(e)<2:raise ValueError('Expression view has fewer than two informative features')
                data.append([e.to_numpy().T]);names.append('expression');features.append(list(map(str,e.index)))
            model=entry_point();model.set_data_options(scale_views=True)
            model.set_data_matrix(data,likelihoods=['gaussian']*len(data),views_names=names,
                samples_names=[list(measured.columns)],features_names=features)
            model.set_model_options(factors=max(1,min(8,complete.shape[1]//5,len(complete)-1)),spikeslab_weights=True)
            model.set_train_options(iter=1000,seed=31,quiet=True,convergence_mode='medium',weight_views=len(data)>1)
            model.build();model.run();model.save(str(out/'mofa.hdf5'),save_data=True)
            variance=model.model.calculate_variance_explained()[0]
            pd.DataFrame([{'view':name,'factor':k+1,'variance_explained':float(variance[v,k])} for v,name in enumerate(names) for k in range(variance.shape[1])]).to_csv(out/'factor_variance.tsv',sep='\t',index=False)
            expectation=model.model.getExpectations();w=expectation['W'][0]['E'];z=expectation['Z']['E']
            np.savez_compressed(out/'mofa_loadings.npz',loadings=w)
            pd.DataFrame(w,index=measured.index).to_csv(out/'loadings.tsv',sep='\t',index_label='event_key')
            pd.DataFrame(z,index=measured.columns).to_csv(out/'scores.tsv',sep='\t',index_label='sample_id')
            for v,name in enumerate(names[1:],1):pd.DataFrame(expectation['W'][v]['E'],index=features[v]).to_csv(out/f'{name}_loadings.tsv',sep='\t')
            norm=np.maximum(np.linalg.norm(w,axis=1),1e-12);unit=w/norm[:,None]
            edges=nearest_edges(lambda a,b:unit[a:b]@unit.T,complete.index)
            status.update(graph_outputs(edges,complete,out,mutual=True))
            # Restore native latent factors after graph summary writes activity.
            pd.DataFrame(z,index=measured.columns).to_csv(out/'scores.tsv',sep='\t',index_label='sample_id')
            status.update(views=names,view_features=dict(zip(names,[len(x) for x in features])),weight_views=len(data)>1,
                          variance_explained=dict(zip(names,map(float,model.model.calculate_variance_explained(total=True)[0]))),
                          factors=w.shape[1],implementation='native mofapy2; graph scored by normalized shared factor loadings')
        else:raise ValueError('Unknown method')
        edges.to_csv(out/'edges.tsv',sep='\t',index=False);status['status']='completed'
    except Exception as exc:status.update(status='failed',error=f'{type(exc).__name__}: {exc}')
    status['seconds']=time.monotonic()-started;write_json(out/'status.json',status);return status

def advanced_scores(method,out,keys,values):
    out=Path(out);p=len(keys)
    if method=='oas':return np.abs(OASNetwork(values).dense())
    if method=='mofa':
        w=np.load(out/'mofa_loadings.npz')['loadings'];w/=np.maximum(np.linalg.norm(w,axis=1,keepdims=True),1e-12)
        s=np.abs(w@w.T);np.fill_diagonal(s,0);return s
    links=pd.read_csv(out/'native_importance.tsv',sep='\t');pos={str(k):i for i,k in enumerate(keys)};s=np.zeros((p,p))
    for r in links.itertuples():s[pos[str(r.target)],pos[str(r.TF)]]=r.importance
    return s
