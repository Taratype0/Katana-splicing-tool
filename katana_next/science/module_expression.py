"""Descriptive and condition-adjusted module/expression associations."""
import numpy as np
import pandas as pd
from scipy.stats import rankdata, t
from katana_next.science.statistics import bh_adjust


def correlate_modules(scores, features, conditions):
    if scores.index.duplicated().any() or features.index.duplicated().any():
        raise ValueError('Duplicate sample identities')
    if set(scores.index)!=set(features.index) or set(scores.index)!=set(conditions.index):
        raise ValueError('Module/expression/design samples must match exactly')
    features=features.loc[scores.index];conditions=conditions.loc[scores.index]
    if not np.isfinite(scores.to_numpy()).all() or not np.isfinite(features.to_numpy()).all():
        raise ValueError('Missing values must be resolved explicitly, not mean-imputed')
    x=rankdata(scores.to_numpy(),axis=0);y=rankdata(features.to_numpy(),axis=0)
    design=np.column_stack([np.ones(len(scores)),pd.get_dummies(conditions,drop_first=True,dtype=float)])
    rows=[]
    for mode in ('overall','condition_adjusted'):
        d=np.ones((len(scores),1)) if mode=='overall' else design
        rank=np.linalg.matrix_rank(d);df=len(scores)-rank-1
        xx=x-d@np.linalg.lstsq(d,x,rcond=None)[0];yy=y-d@np.linalg.lstsq(d,y,rcond=None)[0]
        den=np.sqrt((xx**2).sum(0)[:,None]*(yy**2).sum(0)[None,:])
        with np.errstate(invalid='ignore',divide='ignore'):
            rho=np.clip(xx.T@yy/den,-1,1)
            p=2*t.sf(np.abs(rho)*np.sqrt(max(df,0)/np.maximum(1-rho*rho,1e-15)),df) if df>0 else np.full_like(rho,np.nan)
        p[den<1e-12]=np.nan;rho[den<1e-12]=np.nan
        for i,module in enumerate(scores.columns):
            for j,feature in enumerate(features.columns):
                rows.append(dict(module=module,feature=feature,association=mode,rho=rho[i,j],pvalue=p[i,j],residual_df=df,n_samples=len(scores)))
    out=pd.DataFrame(rows)
    out['FDR']=np.nan
    for mode,idx in out.groupby('association').groups.items():
        valid=out.loc[idx,'pvalue'].dropna()
        out.loc[valid.index,'FDR']=bh_adjust(valid.to_numpy())
    return out
