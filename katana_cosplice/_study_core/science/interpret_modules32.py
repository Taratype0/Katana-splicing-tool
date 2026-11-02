"""Descriptive program links and pathway labels from already-computed evidence."""
import numpy as np,pandas as pd

def interpret(tab,pathways=None,scores=None,factors=None):
 tab=tab.copy();links=pd.DataFrame()
 if pathways is not None and len(pathways):
  p=pathways[pathways.q_value.le(.05)&pathways.overlap.ge(3)].sort_values(['q_value','p_value']).copy()
  p['module']=p.module.astype(str).map(lambda s:s if s.startswith('M') else 'M'+s)
  names=p.groupby('module').pathway.agg(lambda s:'; '.join(s.drop_duplicates().head(3).str.replace('REACTOME_','',regex=False).str.replace('_',' ',regex=False)))
  tab['leading_pathways']=tab.module.map(names).fillna('No significant curated pathway support')
  tab['best_pathway_fdr']=tab.module.map(p.groupby('module').q_value.min())
 else:
  tab['leading_pathways']='Not computed: no gene-set annotation supplied' if pathways is None else 'No significant curated pathway support'
  tab['best_pathway_fdr']=np.nan
 if scores is not None and factors is not None:
  if not set(scores.index)<=set(factors.index):raise ValueError('Factor and module sample identities differ')
  factors=factors.reindex(scores.index);Y=scores.to_numpy();F=factors.to_numpy();yc=Y-Y.mean(axis=0);fc=F-F.mean(axis=0)
  den=np.linalg.norm(yc,axis=0)[:,None]*np.linalg.norm(fc,axis=0)[None,:];cor=np.divide(yc.T@fc,den,out=np.zeros_like(den),where=den>1e-12)
  design=np.column_stack([np.ones(len(F)),F]);pred=design@np.linalg.lstsq(design,Y,rcond=None)[0]
  r2=1-((Y-pred)**2).sum(axis=0)/np.maximum((yc*yc).sum(axis=0),1e-30)
  best=np.argmax(np.abs(cor),axis=1)
  description=pd.DataFrame(dict(module=scores.columns,largest_correlated_factor=best+1,factor_score_correlation=cor[np.arange(len(best)),best],variance_explained_by_all_factors=r2))
  tab=tab.merge(description,on='module',validate='one_to_one')
  links=pd.DataFrame([dict(module=m,factor=j+1,score_correlation=float(cor[i,j])) for i,m in enumerate(scores.columns) for j in range(F.shape[1])])
 return tab,links
