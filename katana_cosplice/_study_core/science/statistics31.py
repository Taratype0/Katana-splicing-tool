"""Restricted label permutations for categorical nuisance strata.

Labels are permuted within saturated nuisance strata. The residualised label
design has the same Gram matrix in every permutation. This is an exact group
invariance argument under conditional exchangeability, with Monte Carlo +1
P values; it is not a heteroscedastic equal-means or causal treatment test.
"""
import numpy as np
import pandas as pd

def stratified_omnibus(scores,meta,factor,strata,permutations=9999,seed=31027):
 if not strata or factor in strata:raise ValueError('Distinct nuisance strata required')
 if meta.biological_sample_id.duplicated().any():raise ValueError('Independent biological units required')
 if not scores.index.equals(meta.index):raise ValueError('Score/sample identities differ')
 if any(s not in meta or meta[s].isna().any() for s in strata):raise ValueError('Missing nuisance stratum')
 result={};patterns={}
 for name in scores:
  key=tuple(np.isfinite(scores[name].to_numpy(float)))
  patterns.setdefault(key,[]).append(name)
 for mask,names in patterns.items():
  m=meta.loc[list(mask)];y=scores.loc[list(mask),names].to_numpy(float)
  levels=m[factor].astype(str);counts=levels.value_counts()
  if len(m)<15 or len(counts)<2 or counts.min()<3:raise ValueError('Insufficient independent units per categorical level')
  group=m[strata].astype(str).agg('\x1f'.join,axis=1)
  groups=[np.asarray(i) for i in group.groupby(group,sort=True).indices.values()]
  d=pd.get_dummies(levels,drop_first=True,dtype=float).to_numpy()
  for ix in groups:y[ix]-=y[ix].mean(axis=0);d[ix]-=d[ix].mean(axis=0)
  q=d.shape[1];gram=d.T@d;df=len(m)-len(groups)-q
  if df<3 or np.linalg.matrix_rank(gram)!=q:raise ValueError('Condition is confounded with nuisance strata')
  inv=np.linalg.inv(gram);total=np.sum(y*y,axis=0)
  if np.any(total<1e-12):raise ValueError('Constant score within nuisance strata')
  cross=d.T@y;observed=np.einsum('im,ij,jm->m',cross,inv,cross)
  rng=np.random.default_rng(seed);exceed=np.zeros(len(names),int)
  for _ in range(permutations):
   order=np.arange(len(m))
   for ix in groups:order[ix]=rng.permutation(ix)
   c=d[order].T@y;null=np.einsum('im,ij,jm->m',c,inv,c)
   exceed+=null>=observed-1e-12*np.maximum(total,1)
  for j,name in enumerate(names):
   result[name]=dict(test='Within-nuisance-stratum condition-label permutation partial F',
     statistic=float((observed[j]/q)/max((total[j]-observed[j])/df,np.finfo(float).tiny)),
     p_value=float((exceed[j]+1)/(permutations+1)),df=float(df),df_num=float(q),
     permutations=permutations,permutation_exceedances=int(exceed[j]),permutation_seed=seed,
     p_resolution=1/(permutations+1),inference_gate_eligible=True,
     nuisance_strata=' x '.join(strata),nuisance_strata_count=len(groups),
     null_hypothesis='Condition labels exchangeable within nuisance strata; not arbitrary unequal-variance equal-means or causal effects')
 return result
