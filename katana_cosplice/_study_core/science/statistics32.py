"""Donor-level inference; no observation-wise pseudo-replication."""
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

def linear_trajectory_tests(scores,meta,permutations=49999,seed=32073):
 """Two-sided whole-donor sign randomization of linear trajectory contrasts.

 Null requires independent, centrally symmetric donor trajectory errors about
 zero mean time response. Arbitrary exchangeability of days is NOT assumed.
 All time points of a donor receive the same sign. Discovery from within-donor
 crossproducts is invariant to these signs; no time labels select membership.
 Tests are associations, not randomized treatment effects or causal claims.
 """
 if not scores.index.equals(pd.Index(meta.sample_id)):raise ValueError('Sample identities differ')
 effects=[];energies=[];donors=[]
 for donor,ix in meta.groupby('biological_sample_id',sort=True).indices.items():
  t=meta.time_day.iloc[ix].to_numpy(float);t-=t.mean()
  if len(set(t))<4:raise ValueError('At least four distinct times per donor required')
  y=scores.iloc[ix].to_numpy(float);y-=y.mean(axis=0)
  effects.append(t@y/(t@t));energies.append(np.mean(y*y,axis=0));donors.append(str(donor))
 effects=np.asarray(effects);energies=np.asarray(energies)
 if len(donors)<15 or not np.isfinite(effects).all():raise ValueError('Insufficient complete independent donor trajectories')
 obs=np.abs(effects.mean(axis=0));rng=np.random.default_rng(seed);exceed=np.zeros(scores.shape[1],int)
 for start in range(0,permutations,1000):
  sign=rng.choice([-1.,1.],size=(min(1000,permutations-start),len(donors)))
  null=np.abs(sign@effects/len(donors));exceed+=(null>=obs[None,:]-1e-14).sum(axis=0)
 p=(exceed+1)/(permutations+1);q=multipletests(p,method='fdr_bh')[1]
 dominance=energies.max(axis=0)/np.maximum(energies.sum(axis=0),1e-30)
 out=pd.DataFrame(dict(module=scores.columns,p_value=p,q_value=q,mean_donor_slope=effects.mean(axis=0),
   slope_standard_error=effects.std(axis=0,ddof=1)/np.sqrt(len(donors)),
   maximum_unit_variance_fraction=dominance,independent_donors=len(donors),permutations=permutations,
   inference_gate_eligible=True,test='whole-donor sign randomization of linear time slope',
   null_assumption='Independent donors; centrally symmetric donor trajectories under no mean linear time response'))
 return out,pd.DataFrame(effects,index=donors,columns=scores.columns)

def trajectory_tests(scores,meta,permutations=49999,seed=32073):
 """Prespecified cubic trajectory omnibus, sign-randomized by whole donor.

 The second-moment studentizer is invariant to donor signs. This allows a
 joint linear/quadratic/cubic time-response test without permuting days, and
 captures transient programs. Linear slopes are descriptive effect summaries.
 """
 if not scores.index.equals(pd.Index(meta.sample_id)):raise ValueError('Sample identities differ')
 t=meta.time_day.to_numpy(float);span=np.ptp(t)
 if not span>0:raise ValueError('Time variation required')
 z=2*(t-t.min())/span-1;design=np.polynomial.legendre.legvander(z,3)
 effects=[];slopes=[];energies=[];donors=[]
 for donor,ix in meta.groupby('biological_sample_id',sort=True).indices.items():
  if np.linalg.matrix_rank(design[ix])<4:raise ValueError('At least four distinct times per donor required')
  y=scores.iloc[ix].to_numpy(float);y-=y.mean(axis=0)
  effects.append(np.linalg.lstsq(design[ix],y,rcond=None)[0][1:]);energies.append(np.mean(y*y,axis=0));donors.append(str(donor))
  tt=t[ix]-t[ix].mean();slopes.append(tt@y/(tt@tt))
 B=np.asarray(effects);slopes=np.asarray(slopes);energies=np.asarray(energies);n=len(donors)
 if n<15 or not np.isfinite(B).all():raise ValueError('Insufficient complete independent donor trajectories')
 C=np.einsum('dqm,drm->mqr',B,B)/n;inverse=np.linalg.pinv(C,hermitian=True);mu=B.mean(axis=0).T
 obs=np.einsum('mq,mqr,mr->m',mu,inverse,mu);rng=np.random.default_rng(seed);exceed=np.zeros(scores.shape[1],int)
 for start in range(0,permutations,1000):
  sign=rng.choice([-1.,1.],size=(min(1000,permutations-start),n));means=np.einsum('bd,dqm->bmq',sign,B)/n
  null=np.einsum('bmq,mqr,bmr->bm',means,inverse,means);exceed+=(null>=obs[None,:]-1e-12).sum(axis=0)
 p=(exceed+1)/(permutations+1);q=multipletests(p,method='fdr_bh')[1]
 out=pd.DataFrame(dict(module=scores.columns,p_value=p,q_value=q,mean_donor_slope=slopes.mean(axis=0),
  slope_standard_error=slopes.std(axis=0,ddof=1)/np.sqrt(n),maximum_unit_variance_fraction=energies.max(axis=0)/np.maximum(energies.sum(axis=0),1e-30),
  independent_donors=n,permutations=permutations,trajectory_statistic=obs,time_basis_degrees=3,inference_gate_eligible=True,
  test='whole-donor sign randomization of cubic time-response omnibus',
  null_assumption='Independent donors; centrally symmetric donor trajectories under no mean cubic-basis time response'))
 return out,pd.DataFrame(slopes,index=donors,columns=scores.columns)

def balanced_scores(x,genes,labels):
 """Each gene contributes equal weight, regardless of its event multiplicity."""
 result={}
 for m in sorted(set(labels)-{0}):
  ix=np.flatnonzero(labels==m);g=genes[ix]
  group=pd.DataFrame(x[ix],index=g).groupby(level=0,sort=True).mean()
  result[f'M{m}']=group.mean(axis=0).to_numpy()
 return pd.DataFrame(result)

def gene_recovery(base,refits,annotation):
 genes=annotation.set_index('event_key').ensembl_gene_id
 def sets(lab):return {m:set(genes.reindex(v.index).dropna()) for m,v in lab[lab>0].groupby(lab[lab>0])}
 b=sets(base);rs=[sets(r) if r is not None else None for r in refits];rows=[]
 for m,g in b.items():
  values=[max([len(g&h)/len(g|h) for h in r.values()] or [0.]) for r in rs if r is not None]
  rows.append(dict(module=f'M{m}',median_gene_jaccard=float(np.median(values)) if values else np.nan))
 return pd.DataFrame(rows)
