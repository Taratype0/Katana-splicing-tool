"""Label-blind latent-program and graph ablations; revision 28 development."""
from pathlib import Path
import time
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA, FactorAnalysis
from cosplice_research import write_json
from stable_network import candidates, labels_from_edges
from advanced import nearest_edges, graph_outputs

CORE = dict(rank=4, minimum_correlation=.6, minimum_margin=.15, minimum_size=20)

def fit_latent(method, measured, complete, meta, out):
    out=Path(out);out.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    if complete.shape[1]<15:raise ValueError('At least 15 independent samples required')
    k=CORE['rank'];z=((complete.T-complete.mean(axis=1))/complete.std(axis=1).clip(lower=1e-8)).to_numpy()
    factor=method.startswith('factor')
    model=FactorAnalysis(n_components=k,rotation='varimax',random_state=31,max_iter=1000) if factor else PCA(n_components=k,svd_solver='full')
    scores=model.fit_transform(z)
    if factor:
        loading=model.components_.T
        loading=loading/np.sqrt(np.sum(loading**2,axis=1)+model.noise_variance_)[:,None]
    else:
        loading=model.components_.T*np.sqrt(model.explained_variance_)
    keys=complete.index.to_numpy(str)
    if method.endswith('_core'):
        ab=np.abs(loading);winner=np.argmax(ab,axis=1);order=np.sort(ab,axis=1)
        eligible=(order[:,-1]>=CORE['minimum_correlation']) & ((order[:,-1]-order[:,-2])>=CORE['minimum_margin'])
        lab=2*winner+1+(loading[np.arange(len(keys)),winner]<0).astype(int)
        lab[~eligible]=0
        for m in set(lab)-{0}:
            if np.sum(lab==m)<CORE['minimum_size']:lab[lab==m]=0
        pd.DataFrame({'event_key':keys,'module':lab}).to_csv(out/'assignments.tsv',sep='\t',index=False)
        pd.DataFrame({'event_key':keys,'strongest_factor':winner+1,'absolute_loading':order[:,-1],
            'loading_margin':order[:,-1]-order[:,-2],'module':lab}).to_csv(out/'membership_evidence.tsv',sep='\t',index=False)
    elif method=='factor_fixed':
        edges=nearest_edges(lambda a,b:loading[a:b]@loading.T,keys)
        graph_outputs(edges,complete,out,mutual=True);edges.to_csv(out/'edges.tsv',sep='\t',index=False)
        lab=pd.read_csv(out/'assignments.tsv',sep='\t').module.to_numpy()
    elif method=='pca_graph':
        # Dot products retain communalities; no normalization to perfect
        # correlations for low-signal events.
        a,b=candidates(loading.astype(np.float32),30)
        w=np.maximum(np.einsum('ij,ij->i',loading[a],loading[b]),0)**2
        lab=labels_from_edges(len(keys),a,b,w,partition='rb',resolution=1,minimum=20)
        pd.DataFrame({'event_key':keys,'module':lab}).to_csv(out/'assignments.tsv',sep='\t',index=False)
        pd.DataFrame({'source':keys[a],'target':keys[b],'weight':w,'sign':1}).to_csv(out/'edges.tsv',sep='\t',index=False)
    else:raise ValueError(method)
    complete.to_csv(out/'input.tsv',sep='\t',index_label='event_key');meta.to_csv(out/'samples.tsv',sep='\t',index=False)
    pd.DataFrame(loading,index=keys,columns=[f'F{i+1}' for i in range(k)]).to_csv(out/'loadings.tsv',sep='\t',index_label='event_key')
    pd.DataFrame(scores,index=complete.columns).to_csv(out/'native_factor_scores.tsv',sep='\t',index_label='sample_id')
    np.savez_compressed(out/'latent_model.npz',components=model.components_,mean=model.mean_,
                        training_event_mean=complete.mean(axis=1).to_numpy(),training_event_scale=complete.std(axis=1).to_numpy())
    status=dict(status='completed',method=method,events=len(keys),samples=len(meta),modules=len(set(lab)-{0}),
                assigned_events=int(np.sum(lab!=0)),seconds=time.monotonic()-started,parameters=CORE,
                label_blind=True,scope='Shared response modules; not causal regulatory interactions')
    write_json(out/'status.json',status);return status
