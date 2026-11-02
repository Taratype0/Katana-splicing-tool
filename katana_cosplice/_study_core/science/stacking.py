from pathlib import Path
import json,time
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from sklearn.model_selection import GroupKFold
from sklearn.linear_model import Ridge
from cosplice_research import ROOT,prepare,fit,write_json
from advanced import fit_advanced,OASNetwork,graph_outputs

METHODS=['wgcna','pca','mofa','grnboost2','oas','correlation']

def scores(method,out,complete):
    x=complete.to_numpy();z=(x-x.mean(axis=1,keepdims=True))/np.maximum(x.std(axis=1,keepdims=True),1e-8)
    if method=='wgcna':
        labels=pd.read_csv(out/'assignments.tsv',sep='\t').set_index('event_key').loc[complete.index,'module'].to_numpy()
        return lambda a,b:np.maximum(z[a:b]@z.T/z.shape[1],0)*((labels[a:b,None]==labels[None,:])&(labels[a:b,None]!=0))
    if method=='correlation':return lambda a,b:np.abs(z[a:b]@z.T/z.shape[1])
    if method=='oas':
        model=OASNetwork(x);return lambda a,b:np.abs(model.rows(a,b))
    if method in ['pca','mofa']:
        if method=='pca':
            f=np.load(out/'factor.npz');w=f['loadings'];sd=np.sqrt(np.sum(w*w,axis=1)+f['noise']);w=w/sd[:,None]
        else:
            w=np.load(out/'mofa_loadings.npz')['loadings'];w=w/np.maximum(np.linalg.norm(w,axis=1,keepdims=True),1e-12)
        return lambda a,b:np.abs(w[a:b]@w.T)
    links=pd.read_csv(out/'native_importance.tsv',sep='\t');pos={k:i for i,k in enumerate(complete.index)}
    mapping={pos[k]:[(pos[r.TF],r.importance) for r in g.itertuples()] for k,g in links.groupby('target')}
    def rows(a,b):
        s=np.zeros((b-a,len(complete)))
        for i in range(a,b):
            for j,v in mapping.get(i,[]):s[i-a,j]=v
        return s
    return rows

def decoder(row_score,complete,genes):
    x=complete.to_numpy();mu=x.mean(axis=1);sd=np.maximum(x.std(axis=1),1e-8);z=(x-mu[:,None])/sd[:,None];p=len(x)
    genes=np.asarray(genes);known=np.asarray([pd.notna(g) and str(g).strip() not in ('','UNKNOWN') for g in genes])
    neighbor=np.full((p,8),-1,dtype=int);coef=np.zeros((p,8));bias=np.zeros(p)
    for start in range(0,p,128):
        block=row_score(start,min(p,start+128))
        for offset,r in enumerate(block):
            i=start+offset;ids=np.flatnonzero((np.arange(p)!=i)&known&known[i]&(genes!=genes[i])&(r>0))
            ids=ids[np.argsort(-r[ids],kind='stable')[:8]]
            if len(ids):
                m=Ridge(alpha=10).fit(z[ids].T,z[i]);neighbor[i,:len(ids)]=ids;coef[i,:len(ids)]=m.coef_;bias[i]=m.intercept_
    return dict(mean=mu,scale=sd,neighbor=neighbor,coef=coef,bias=bias,keys=np.asarray(complete.index,str))

def predict_decoder(model,test):
    y=test.reindex(model['keys']).to_numpy();z=(y-model['mean'][:,None])/model['scale'][:,None];z=np.nan_to_num(z)
    n=model['neighbor'];valid=n>=0;pred=np.sum(z[np.maximum(n,0)]*model['coef'][:,:,None]*valid[:,:,None],axis=1)
    return model['mean'][:,None]+model['scale'][:,None]*(pred+model['bias'][:,None])

def expert(method,measured,complete,meta,genes,out):
    status=json.loads((out/'status.json').read_text()) if (out/'status.json').exists() else (
        fit_advanced(method,measured,complete,meta,out) if method in ['mofa','grnboost2','oas'] else fit(method,complete,meta,out))
    if status['status']!='completed':raise RuntimeError(f'{method} failed: {status}')
    model=decoder(scores(method,out,complete),complete,genes);np.savez_compressed(out/'decoder.npz',**model);return model

def convex_weights(pred,y):
    mask=np.isfinite(y);a=pred[:,mask].T;b=y[mask]
    # Scale objective to avoid premature SLSQP stopping on small usage MSE.
    gram=a.T@a/len(b);cross=a.T@b/len(b)
    objective=lambda w:float(1000*(w@gram@w-2*w@cross))
    jac=lambda w:2000*(gram@w-cross)
    result=minimize(objective,np.ones(len(pred))/len(pred),jac=jac,bounds=[(0,1)]*len(pred),
        constraints={'type':'eq','fun':lambda w:w.sum()-1,'jac':lambda w:np.ones(len(w))},method='SLSQP',options={'ftol':1e-11,'maxiter':2000})
    if not result.success:raise RuntimeError('Stacking optimization failed: '+result.message)
    return result.x

def train_stacking(measured,meta,genes,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    if meta.biological_sample_id.nunique()<23:raise ValueError('Stacking requires at least 23 biological units')
    started=time.monotonic();p,n=measured.shape;pred=np.full((len(METHODS)+1,p,n),np.nan);records=[]
    splitter=GroupKFold(n_splits=3,shuffle=True,random_state=41)
    for fold,(tr,va) in enumerate(splitter.split(meta,groups=meta.biological_sample_id)):
        raw=measured.iloc[:,tr];m,c,a,means,md=prepare(raw,meta.iloc[tr])
        dest=out/f'fold_{fold}';dest.mkdir()
        baseline=raw.mean(axis=1).fillna(.5).to_numpy();pred[0,:,va]=np.broadcast_to(baseline,(len(va),p))
        records.append({'fold':fold,'train_samples':meta.iloc[tr].sample_id.tolist(),'validation_samples':meta.iloc[va].sample_id.tolist(),
                        'train_units':meta.iloc[tr].biological_sample_id.unique().tolist(),'validation_units':meta.iloc[va].biological_sample_id.unique().tolist()})
        for mi,method in enumerate(METHODS,1):
            model=expert(method,m,c,md,genes.reindex(c.index).to_numpy(),dest/method)
            result=pd.DataFrame(predict_decoder(model,measured.iloc[:,va]),index=c.index,columns=measured.columns[va]).reindex(measured.index)
            result=result.T.fillna(pd.Series(baseline,index=measured.index)).T
            pred[mi,:,va]=result.to_numpy().T
        write_json(out/'folds.json',records)
    weights=convex_weights(pred,measured.to_numpy());calibrated=convex_weights(pred[:2],measured.to_numpy())
    error=np.nanmean((pred-measured.to_numpy()[None,:,:])**2,axis=(1,2))
    pd.DataFrame({'expert':['mean']+METHODS,'weight':weights,'internal_oof_mse':error}).to_csv(out/'selection.tsv',sep='\t',index=False)
    np.savez_compressed(out/'oof.npz',predictions=pred,observed=measured.to_numpy(),weights=weights,calibrated_wgcna_weights=calibrated)
    measured,complete,audit,means,meta=prepare(measured,meta)
    edges=[]
    for mi,method in enumerate(METHODS,1):
        model=expert(method,measured,complete,meta,genes.to_numpy(),out/method)
        for i,key in enumerate(complete.index):
            for j,c in zip(model['neighbor'][i],model['coef'][i]):
                if j>=0 and abs(c)*weights[mi]>1e-12:edges.append((str(key),str(complete.index[j]),abs(float(c))*weights[mi],float(np.sign(c))))
    edge=pd.DataFrame(edges,columns=['source','target','weight','sign']);edge=edge.groupby(['source','target'],as_index=False).agg(weight=('weight','sum'),sign=('sign','mean'))
    edge.to_csv(out/'edges.tsv',sep='\t',index=False);graph=graph_outputs(edge,complete,out)
    np.savez_compressed(out/'stacking.npz',weights=weights,calibrated_wgcna_weights=calibrated,mean=complete.mean(axis=1).to_numpy(),keys=complete.index.to_numpy(str))
    status={'status':'completed','method':'stacking','seconds':time.monotonic()-started,'events':len(complete),'samples':len(meta),
        'weights':dict(zip(['mean']+METHODS,weights.tolist())),'calibrated_wgcna_weights':calibrated.tolist(),
        'implementation':'biological-unit-disjoint cross-fitted convex prediction ensemble; not a causal network',**graph}
    write_json(out/'status.json',status);return status

def predict_stacking(out,test):
    out=Path(out);a=np.load(out/'stacking.npz');y=test.reindex(a['keys']);pred=[np.broadcast_to(a['mean'][:,None],y.shape)]
    for method in METHODS:
        model=dict(np.load(out/method/'decoder.npz'));pred.append(predict_decoder(model,y))
    pred=np.asarray(pred)
    return np.einsum('m,mpn->pn',a['weights'],pred),np.einsum('m,mpn->pn',a['calibrated_wgcna_weights'],pred[:2])
