"""Conditional gene-set enrichment controlling detectable event multiplicity.

Within each fixed event-count stratum, the module selects n genes without
replacement. The pathway overlap is a sum of independent hypergeometrics;
convolve their PMFs. No Monte Carlo P-value floor or outcome-driven gene set.
"""
import numpy as np,pandas as pd
from scipy.stats import hypergeom
from scipy.fft import rfft,irfft,next_fast_len
from statsmodels.stats.multitest import multipletests
from validation24 import pathway_sets

def conditional_tail(N,K,n,observed):
 N=np.asarray(N,int);K=np.asarray(K,int);n=np.asarray(n,int);observed=np.asarray(observed,int)
 length=next_fast_len(int(n.sum())+1);frequency=np.ones((len(K),length//2+1),complex)
 for j in range(len(N)):
  if N[j]==0 or n[j]==0:continue
  support=np.arange(n[j]+1);pmf=hypergeom.pmf(support[None,:],N[j],K[:,j,None],n[j])
  frequency*=rfft(pmf,n=length,axis=1)
 dist=np.maximum(irfft(frequency,n=length,axis=1),0);dist/=dist.sum(axis=1,keepdims=True)
 tail=np.flip(np.cumsum(np.flip(dist,axis=1),axis=1),axis=1)
 # Conservative numerical floor: FFT cancellation is not evidence for P=0.
 return np.maximum(tail[np.arange(len(K)),observed],1e-12)

def context(events,path):
 count=events.groupby('gene_id').size();background=set(count.index);sets=pathway_sets(path,background)
 strata=np.digitize(count.to_numpy(),[2,3,5,9]);N=np.bincount(strata,minlength=5)
 lookup=pd.Series(strata,index=count.index);names=list(sets)
 K=np.array([np.bincount(lookup.loc[list(sets[name][1])].to_numpy(),minlength=5) for name in names],dtype=int).reshape(-1,5)
 return dict(background=background,sets=sets,names=names,N=N,K=K,lookup=lookup)

def tests(module_genes,module,ctx):
 genes=set(module_genes)&ctx['background'];n=np.bincount(ctx['lookup'].loc[list(genes)].to_numpy(),minlength=5)
 hits=[genes&ctx['sets'][name][1] for name in ctx['names']];observed=np.array([len(h) for h in hits]);sizes=ctx['K'].sum(axis=1)
 return pd.DataFrame(dict(module=module,pathway=ctx['names'],module_genes=len(genes),background_genes=len(ctx['background']),term_genes=sizes,
  overlap=observed,contributing_genes=[';'.join(sorted(h)) for h in hits],
  p_value=conditional_tail(ctx['N'],ctx['K'],n,observed),
  unadjusted_event_multiplicity_hypergeom_p=hypergeom.sf(observed-1,len(ctx['background']),sizes,len(genes)),
  test='Exact event-multiplicity-stratified gene overlap; five strata (1,2,3-4,5-8,9+ QC events/gene)'))

def enrich(events,path,out):
 ctx=context(events,path);blocks=[tests(part.gene_id,int(m),ctx) for m,part in events[events.module>0].groupby('module')]
 result=pd.concat(blocks,ignore_index=True);result['q_value']=multipletests(result.p_value,method='fdr_bh')[1] if len(result) else []
 result['genes']=result.contributing_genes;result['pathway_background_genes']=result.term_genes
 result.to_csv(out,sep='\t',index=False);return result
