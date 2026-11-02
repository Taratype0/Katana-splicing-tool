"""Donor-centered latent covariance, with the unchanged common local partition.

Weighted EBMF can leave fitted score offsets even when its input was centered.
Remove those offsets at the same donor level as the common input. This makes
the covariance construction invariant to a whole-donor sign change, including
with unequal trajectory lengths; no label, pathway or outer endpoint is used.
"""
from pathlib import Path
import numpy as np,pandas as pd,json
from common32 import sha,write
from graph32 import neighbours,make_edges as native_edges
from graph32_product import hierarchical_partition


def centered_factor_covariance(scores,units):
 f=np.asarray(scores,float).copy();units=np.asarray(units)
 for unit in np.unique(units):
  ix=units==unit;f[ix]-=f[ix].mean(axis=0)
 # After block centering the global mean is zero. The second moment avoids
 # an unnecessary second centering and is explicitly invariant to block signs.
 return f.T@f/len(f)


def make_edges(method,input,model,out):
 input,model,out=map(Path,[input,model,out]);out.mkdir(parents=True,exist_ok=True)
 if method!='factor_adaptive':return native_edges(method,input,model,out)
 dest=out/'edges.tsv.gz';manifest=out/'edge_manifest.json'
 if dest.exists():
  info=json.loads(manifest.read_text())
  assert info['source_sha256']==sha(Path(__file__)) and info['input_sha256']==sha(input/'input.bin')
  assert info['loading_sha256']==sha(model/'loadings.tsv') and info['score_sha256']==sha(model/'factor_scores.tsv')
  return pd.read_csv(dest,sep='\t')
 events=pd.read_csv(input/'events.tsv',sep='\t');meta=pd.read_csv(input/'samples.tsv',sep='\t',dtype={'biological_sample_id':str})
 L=pd.read_csv(model/'loadings.tsv',sep='\t',index_col=0);F=pd.read_csv(model/'factor_scores.tsv',sep='\t',index_col=0)
 assert np.array_equal(L.index,events.event_key) and np.array_equal(F.index,meta.sample_id)
 sf=centered_factor_covariance(F.to_numpy(),meta.biological_sample_id.to_numpy());ev,U=np.linalg.eigh(sf)
 embedding=L.to_numpy()@(U*np.sqrt(np.maximum(ev,0)))
 embedding/=np.maximum(np.linalg.norm(embedding,axis=1,keepdims=True),1)
 np.save(out/'latent_embedding.npy',embedding.astype(np.float32));rows=[];genes=events.ensembl_gene_id.to_numpy();p=len(events)
 for start in range(0,p,256):
  ix=np.arange(start,min(start+256,p));cor=np.clip(embedding[ix]@embedding.T,-1,1)
  rows.extend(neighbours(np.maximum(cor,0)**2,ix,np.arange(p),genes))
 t=pd.DataFrame(rows,columns=['a','b','weight']);a=np.minimum(t.a,t.b);b=np.maximum(t.a,t.b);t['a']=a;t['b']=b
 t=t.groupby(['a','b'],sort=True).weight.max().reset_index();t.to_csv(dest,sep='\t',index=False)
 write(manifest,dict(method=method,events=p,edges=len(t),neighbours=30,input_sha256=sha(input/'input.bin'),
  source_sha256=sha(Path(__file__)),loading_sha256=sha(model/'loadings.tsv'),score_sha256=sha(model/'factor_scores.tsv'),
  scope='cross-gene edges; one nearest event per other gene; symmetric union',
  construction='posterior factor covariance after within-donor score centering; shared recursive partition unchanged'))
 return t


def run(method,input,model,out):
 edges=make_edges(method,input,model,out)
 return hierarchical_partition(input,edges,Path(out)/'hierarchical')
