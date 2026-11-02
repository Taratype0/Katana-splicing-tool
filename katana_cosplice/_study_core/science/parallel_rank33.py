"""Engineering candidate: block-sign parallel analysis, before outer outcomes.

Adapt Horn's spectral noise comparison to repeated donor trajectories and
multi-event genes. A sign is shared by every time point and event of each
gene/donor block. Thus within-donor sample covariance, gene compositional
structure, observation masks and marginal variances are retained. Cross-gene
agreement across donors is removed. This is a structured adaptation, not the
independent-entry null or a claimed gene/module FDR procedure.
"""
from pathlib import Path
import argparse,time
import numpy as np,pandas as pd
from common32 import ROOT,EVID,write,sha

def select(input,out,repetitions=999,seed=33041):
 input,out=map(Path,[input,out]);out.mkdir(parents=True,exist_ok=True)
 ev=pd.read_csv(input/'events.tsv',sep='\t');meta=pd.read_csv(input/'samples.tsv',sep='\t',dtype={'biological_sample_id':str})
 x=np.ascontiguousarray(np.fromfile(input/'input.bin',dtype='<f8').reshape(len(ev),len(meta)).T)
 units,unames=pd.factorize(meta.biological_sample_id,sort=True);genes,gnames=pd.factorize(ev.ensembl_gene_id,sort=True)
 if len(unames)<15:raise ValueError('At least 15 independent units required')
 cap=min(64,len(meta)-len(unames)-1)
 observed=np.linalg.eigvalsh(x@x.T)[::-1][:cap]
 rng=np.random.default_rng(seed);null=np.empty((repetitions,cap));start=time.time()
 for b in range(repetitions):
  signs=rng.choice(np.array([-1.,1.]),size=(len(unames),len(gnames)))
  z=x*signs[units[:,None],genes[None,:]]
  null[b]=np.linalg.eigvalsh(z@z.T)[::-1][:cap]
  if (b+1)%(100 if repetitions>100 else 20)==0:print('Parallel null',b+1,'/',repetitions,flush=True)
 limit=np.quantile(null,.95,axis=0,method='higher');above=observed>limit
 rank=int(np.flatnonzero(~above)[0]) if (~above).any() else cap
 table=pd.DataFrame(dict(rank=np.arange(1,cap+1),observed_eigenvalue=observed,null_mean=null.mean(axis=0),null_q95=limit,above_null=above,retained=np.arange(1,cap+1)<=rank))
 table.to_csv(out/'parallel_rank.tsv',sep='\t',index=False);np.save(out/'parallel_null_eigenvalues.npy',null)
 result=dict(selected_rank=rank,safety_cap=cap,safety_cap_reached=rank==cap,repetitions=repetitions,seed=seed,seconds=time.time()-start,
  input_sha256=sha(input/'input.bin'),gene_mapping_sha256=sha(input/'events.tsv'),samples_sha256=sha(input/'samples.tsv'),
  rule='Consecutive leading eigenvalues above their 95th percentile gene-by-donor block-sign reference; no outer outcome is used.',
  null='One independent Rademacher sign per gene and donor, shared across events and time points; within-unit covariance and same-gene dependence retained.',
  interpretation='Rank-selection heuristic for shared programs, not a calibrated module or gene FDR test.',
  borrowed_method='FlipPA (Hong, Sheng, Dobriban); structured gene-by-donor sign adaptation, not the original independent-entry algorithm',
  source='https://arxiv.org/abs/2012.02985')
 result['selector_sha256']=sha(Path(__file__))
 write(out/'parallel_rank.json',result);print(result,flush=True);return result

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('split',default='full',nargs='?');p.add_argument('--repetitions',type=int,default=999);a=p.parse_args()
 select(EVID/'inputs'/a.split,EVID/'development/parallel_rank33'/a.split,repetitions=a.repetitions)
