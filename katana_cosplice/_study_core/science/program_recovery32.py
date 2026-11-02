"""Descriptive global-program recovery, distinct from local-module validation."""
from pathlib import Path
import numpy as np,pandas as pd
from scipy.optimize import linear_sum_assignment

def program_recovery(full_input,full_model,refit_models,out):
 full_input,full_model,out=map(Path,[full_input,full_model,out]);out.mkdir(parents=True,exist_ok=True)
 ev=pd.read_csv(full_input/'events.tsv',sep='\t');keys=pd.Index(ev.event_key)
 l=pd.read_csv(full_model/'loadings.tsv',sep='\t',index_col=0).reindex(keys).to_numpy()
 f=pd.read_csv(full_model/'factor_scores.tsv',sep='\t',index_col=0)
 geneweights=1/np.sqrt(ev.groupby('ensembl_gene_id').ensembl_gene_id.transform('size').to_numpy())
 a=l*geneweights[:,None];an=np.linalg.norm(a,axis=0);k=a.shape[1];records=[];spaces=[];successful=0
 def orth(x):
  u,s,_=np.linalg.svd(x,full_matrices=False);return u[:,s>max(s[0]*1e-10,1e-12)]
 qa=orth(a)
 for i,path in enumerate(refit_models):
  p=Path(path)
  if not (p/'COMPLETED').exists():continue
  b=pd.read_csv(p/'loadings.tsv',sep='\t',index_col=0).reindex(keys).fillna(0).to_numpy()*geneweights[:,None]
  g=pd.read_csv(p/'factor_scores.tsv',sep='\t',index_col=0);bn=np.linalg.norm(b,axis=0)
  cos=np.divide(a.T@b,an[:,None]*bn[None,:],out=np.zeros((k,b.shape[1])),where=an[:,None]*bn[None,:]>1e-12)
  ii,jj=linear_sum_assignment(-np.abs(cos));matched=dict(zip(ii,jj));successful+=1
  # Loading-based one-to-one assignment avoids relabelling each program by
  # whichever independent score statistic happens to look best.
  shared=f.index.intersection(g.index);fs=f.loc[shared].to_numpy(copy=True);gs=g.loc[shared].to_numpy(copy=True);fs-=fs.mean(axis=0);gs-=gs.mean(axis=0)
  den=np.linalg.norm(fs,axis=0)[:,None]*np.linalg.norm(gs,axis=0)[None,:]
  score=np.divide(fs.T@gs,den,out=np.zeros_like(cos),where=den>1e-12)
  for j in range(k):
   hit=matched.get(j);c=cos[j,hit] if hit is not None else 0
   records.append(dict(refit=i,factor=j+1,matched_factor=hit+1 if hit is not None else None,loading_cosine=abs(c),score_correlation=float(score[j,hit]*np.sign(c)) if hit is not None else 0,matched=hit is not None))
  qb=orth(b);overlap=np.sum((qa.T@qb)**2)
  spaces.append(dict(refit=i,full_rank=qa.shape[1],refit_rank=qb.shape[1],symmetric_subspace_overlap=2*overlap/(qa.shape[1]+qb.shape[1])))
 tab=pd.read_csv(full_model/'factors.tsv',sep='\t');tab['requested_refits']=len(refit_models);tab['successful_refits']=successful
 if records:
  r=pd.DataFrame(records);agg=r.groupby('factor').agg(median_loading_cosine=('loading_cosine','median'),lower10_loading_cosine=('loading_cosine',lambda v:v.quantile(.1)),median_score_correlation=('score_correlation','median'),matched_refit_fraction=('matched','mean')).reset_index();tab=tab.merge(agg,on='factor')
 else:
  r=pd.DataFrame(columns=['refit','factor','matched_factor','loading_cosine','score_correlation','matched'])
  for col in ['median_loading_cosine','lower10_loading_cosine','median_score_correlation','matched_refit_fraction']:tab[col]=np.nan
 r.to_csv(out/'program_refit_matches.tsv',sep='\t',index=False);pd.DataFrame(spaces,columns=['refit','full_rank','refit_rank','symmetric_subspace_overlap']).to_csv(out/'program_subspace_recovery.tsv',sep='\t',index=False)
 tab.to_csv(out/'global_program_summary.tsv',sep='\t',index=False)
 return tab
