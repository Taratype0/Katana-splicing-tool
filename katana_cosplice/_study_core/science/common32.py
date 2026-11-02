"""Identical donor-aware inputs and immutable splits for both engines."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];DATA=ROOT/'data/revision32/recount3/SRP168391';EVID=ROOT/'evidence/revision32'
def write(p,obj):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(obj,indent=2,default=lambda x:x.item()),encoding='utf8')
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def load():
 events=pd.read_csv(DATA/'events.tsv',sep='\t',dtype={'gene_id':str,'ensembl_gene_id':str})
 meta=pd.read_csv(DATA/'samples.tsv',sep='\t',dtype={'biological_sample_id':str,'sample_id':str})
 return events,meta,np.load(DATA/'usage.npy',mmap_mode='r'),np.load(DATA/'counts.npy',mmap_mode='r'),np.load(DATA/'support.npy',mmap_mode='r')
def make_splits():
 path=EVID/'protocol.json'
 if path.exists():return json.loads(path.read_text())
 meta=pd.read_csv(DATA/'samples.tsv',sep='\t',dtype=str);donors=sorted(meta.biological_sample_id.unique());rng=np.random.default_rng(32041)
 draws=[]
 while len(draws)<20:
  d=sorted(rng.choice(donors,15,replace=False).tolist())
  if d not in draws:draws.append(d)
 folds=[sorted(x.tolist()) for x in np.array_split(rng.permutation(donors),5)]
 splits={'full':donors,**{f'refit_{i:02}':d for i,d in enumerate(draws)},**{f'fold_{i:02}':sorted(set(donors)-set(d)) for i,d in enumerate(folds)}}
 obj=dict(dataset='recount3 SRP168391',donor_splits=splits,heldout_donors=folds,seed=32041,
  qc=dict(observed_fraction=.8,donor_fraction_with_80pct_observed=.8,within_donor_usage_sd=.02,maximum_donor_variance_fraction_lt=.5,minimum_cluster_reads=10),
  transformation='Jeffreys posterior mean junction/cluster usage, arcsin sqrt, within-donor centering, zero fill after centering, sqrt inverse-donor-size row weights, event unit variance',
  comparison_scope='identical donor-centered and equally donor-weighted input; common training-only assay QC, no condition significance filter',
  model_configuration_status='shape-only development; freeze before opening outer outcomes',
  source_data_manifest_sha256=sha(DATA/'manifest.json'))
 write(path,obj);return obj
def prepare(donors,out,source=None):
 out=Path(out);out.mkdir(parents=True,exist_ok=True)
 if (out/'manifest.json').exists():
  cached=json.loads((out/'manifest.json').read_text())
  if cached['donor_ids']!=sorted(map(str,donors)) or sha(out/'input.bin')!=cached['input_sha256']:
   raise ValueError('Cached training input identity mismatch')
  return cached
 events,meta,usage,counts,support=source or load();cols=np.flatnonzero(meta.biological_sample_id.isin(donors));m=meta.iloc[cols].reset_index(drop=True)
 if m.biological_sample_id.nunique()<15:raise ValueError('At least 15 independent donor trajectories required')
 u=np.asarray(usage[:,cols],float);c=np.asarray(counts[:,cols],float);s=np.asarray(support[:,cols],float)
 observed=np.isfinite(u);coverage=observed.mean(axis=1);groups=list(m.groupby('biological_sample_id',sort=True).indices.values())
 donor_coverage=np.column_stack([observed[:,ix].mean(axis=1) for ix in groups]);donor_pass=(donor_coverage>=.8).mean(axis=1)
 centered=u.copy();transformed=np.arcsin(np.sqrt((c+.5)/(s+1)));transformed[~observed]=np.nan
 weights=np.ones(len(m));unitvar=[]
 for ix in groups:
  n=observed[:,ix].sum(axis=1);mu=np.divide(np.nansum(u[:,ix],axis=1),n,out=np.zeros(len(u)),where=n>0)
  centered[:,ix]-=mu[:,None]
  mt=np.divide(np.nansum(transformed[:,ix],axis=1),n,out=np.zeros(len(u)),where=n>0)
  transformed[:,ix]-=mt[:,None];weights[ix]=len(m)/len(groups)/len(ix)
  unitvar.append(np.nansum(centered[:,ix]**2,axis=1)/len(ix))
 uv=np.column_stack(unitvar);total=uv.sum(axis=1);sd=np.sqrt(total/len(groups));dom=np.divide(uv.max(axis=1),total,out=np.ones(len(u)),where=total>0)
 keep=(coverage>=.8)&(donor_pass>=.8)&(sd>=.02)&(dom<.5)
 audit=pd.DataFrame(dict(event_key=events.event_key,observed_fraction=coverage,donor_coverage_pass_fraction=donor_pass,
     within_donor_usage_sd=sd,maximum_donor_variance_fraction=dom,eligible=keep));audit.to_csv(out/'event_qc.tsv.gz',sep='\t',index=False)
 x=np.nan_to_num(transformed[keep],nan=0.0)*np.sqrt(weights)[None,:];x-=x.mean(axis=1,keepdims=True)
 scale=np.sqrt(np.mean(x*x,axis=1));x/=scale[:,None]
 if not np.isfinite(x).all() or len(x)<100:raise ValueError('Invalid or insufficient training matrix')
 e=events.loc[keep].reset_index(drop=True);e['source_event_row']=np.flatnonzero(keep)
 e.to_csv(out/'events.tsv',sep='\t',index=False);m.to_csv(out/'samples.tsv',sep='\t',index=False)
 (out/'event_keys.txt').write_text('\n'.join(e.event_key)+'\n');np.asarray(x,dtype='<f8').tofile(out/'input.bin')
 np.save(out/'transform_scale.npy',scale);np.save(out/'observation_weights.npy',weights)
 result=dict(events=len(e),genes=e.gene_id.nunique(),samples=len(m),donors=len(groups),source_event_rows=len(events),
   donor_ids=sorted(map(str,donors)),input_sha256=sha(out/'input.bin'),event_keys_sha256=sha(out/'event_keys.txt'),
   analysis='within-donor splicing trajectories; technical repeats are not biological units')
 write(out/'manifest.json',result);print('prepared',out.name,result,flush=True);return result

if __name__=='__main__':
 p=make_splits();prepare(p['donor_splits']['full'],EVID/'inputs/full')
