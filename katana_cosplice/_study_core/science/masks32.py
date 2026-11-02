"""Deterministic inner donor-by-gene holdouts, independent of observed values."""
import hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
from common32 import DATA,EVID,write,sha
def prepare_masks(input):
 input=Path(input)
 if (input/'mask_manifest.json').exists() and (input/'sampling_sd.bin').exists():return
 events=pd.read_csv(input/'events.tsv',sep='\t');samples=pd.read_csv(input/'samples.tsv',sep='\t',dtype={'biological_sample_id':str})
 original=pd.read_csv(DATA/'samples.tsv',sep='\t',dtype=str);cols=pd.Index(original.sample_id).get_indexer(samples.sample_id)
 observed=np.isfinite(np.load(DATA/'usage.npy',mmap_mode='r')[np.ix_(events.source_event_row,cols)])
 genes=events.ensembl_gene_id.unique();donors=sorted(samples.biological_sample_id.unique())
 mask=np.array([[int.from_bytes(hashlib.blake2b(f'32051|{g}|{d}'.encode(),digest_size=4).digest(),'little')%5==0 for d in donors] for g in genes])
 gi=pd.Index(genes).get_indexer(events.ensembl_gene_id);di=pd.Index(donors).get_indexer(samples.biological_sample_id)
 blockmask=mask[np.ix_(gi,di)]
 blockmask.astype('u1').tofile(input/'inner_mask.bin');observed.astype('u1').tofile(input/'observed_mask.bin')
 support=np.load(DATA/'support.npy',mmap_mode='r')[np.ix_(events.source_event_row,cols)].astype(float)
 # Delta-method variance of arcsin(sqrt(Jeffreys posterior usage)); propagate
 # within-donor mean removal, observation weights and event scaling explicitly.
 v=1/(4*(support+2));v[~observed]=0;noise=np.zeros_like(v)
 for ix in samples.groupby('biological_sample_id',sort=True).indices.values():
  n=observed[:,ix].sum(axis=1).clip(min=1)
  noise[:,ix]=v[:,ix]*(1-2/n[:,None])+v[:,ix].sum(axis=1,keepdims=True)/(n[:,None]**2)
 weights=np.load(input/'observation_weights.npy');scale=np.load(input/'transform_scale.npy')
 noise=np.maximum(noise,0)*weights[None,:]/scale[:,None]**2
 np.sqrt(np.maximum(noise,1e-12)).astype('<f8').tofile(input/'sampling_sd.bin')
 write(input/'mask_manifest.json',dict(kind='entire gene by entire donor trajectory, hash-assigned 20 percent internal holdout',
    inner_mask_sha256=sha(input/'inner_mask.bin'),observed_mask_sha256=sha(input/'observed_mask.bin'),
    heldout_fraction=float(blockmask.mean()),validation_observed_fraction=float((blockmask&observed).mean()),
    training_only='outer test donors are absent; internal masks never use time labels, pathways or outer recovery',
    sampling_sd_sha256=sha(input/'sampling_sd.bin'),median_sampling_variance_fraction=float(np.median((noise*observed).sum(axis=1)/observed.sum(axis=1))),
    uncertainty='binomial delta approximation, propagated through donor centering; within-donor error correlations not modeled by native diagonal-S likelihood; donor-block validation remains required'))
if __name__=='__main__':
 for p in sorted((EVID/'inputs').iterdir()):
  if p.is_dir():prepare_masks(p)
