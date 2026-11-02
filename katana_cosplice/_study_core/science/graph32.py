"""Sparse latent-profile modules and a controlled native WGCNA-TOM graph arm.

Both graphs use the same cross-gene neighbour rule and Leiden resolution rule.
Native WGCNA partitions are preserved as a separate, unmodified output.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import igraph as ig
import leidenalg as la
from common32 import write,sha

PARAMETERS=dict(neighbours=30,min_events=20,min_genes=5,max_genes=300,
 target_effective_events=[100,200,400],gamma_grid=[float(x) for x in 2**np.linspace(-2,6,25)],seed=32041)

def neighbours(sim,own_global,candidate_global,genes,k=30):
 """At most one neighbour per other gene; no compositional within-gene edge."""
 rows=[];n=len(candidate_global)
 for i,own in enumerate(own_global):
  v=sim[i];v[genes[candidate_global]==genes[own]]=-np.inf
  width=min(n,max(128,k*4));chosen=[]
  while True:
   ix=np.argpartition(v,n-width)[n-width:] if width<n else np.arange(n)
   ix=ix[np.argsort(-v[ix],kind='stable')];seen=set();chosen=[]
   for j in ix:
    if not np.isfinite(v[j]) or v[j]<=0:break
    g=genes[candidate_global[j]]
    if g in seen:continue
    seen.add(g);chosen.append((int(own),int(candidate_global[j]),float(v[j])))
    if len(chosen)==k:break
   if len(chosen)>=k or width==n:break
   width=min(n,width*2)
  rows.extend(chosen)
 return rows

def make_edges(method,input,model,out):
 input=Path(input);model=Path(model);out=Path(out);out.mkdir(parents=True,exist_ok=True);dest=out/'edges.tsv.gz'
 if dest.exists():return pd.read_csv(dest,sep='\t')
 events=pd.read_csv(input/'events.tsv',sep='\t');genes=events.ensembl_gene_id.to_numpy();p=len(events);rows=[]
 if method=='factor_adaptive':
  L=pd.read_csv(model/'loadings.tsv',sep='\t',index_col=0);F=pd.read_csv(model/'factor_scores.tsv',sep='\t',index_col=0)
  if not np.array_equal(L.index,events.event_key):raise ValueError('Factor event identity mismatch')
  sf=np.cov(F.to_numpy(),rowvar=False,ddof=0);sf=np.atleast_2d(sf);ev,U=np.linalg.eigh(sf)
  embedding=L.to_numpy()@(U*np.sqrt(np.maximum(ev,0)))
  # The common input has per-event variance 1; cap small numerical/model excess.
  embedding/=np.maximum(np.linalg.norm(embedding,axis=1,keepdims=True),1)
  np.save(out/'latent_embedding.npy',embedding.astype(np.float32))
  for start in range(0,p,256):
   ix=np.arange(start,min(start+256,p));cor=np.clip(embedding[ix]@embedding.T,-1,1)
   # Positive model-implied correlation; unsupported zero loadings create no edge.
   sim=np.maximum(cor,0)**2;rows.extend(neighbours(sim,ix,np.arange(p),genes))
 elif method=='wgcna_tom':
  for b in json.loads((model/'tom_blocks.json').read_text()):
   ix=np.asarray(b['indices']);a=np.memmap(model/b['file'],mode='r',dtype='<f4',shape=(b['n'],b['n']))
   for start in range(0,len(ix),256):
    sub=ix[start:start+256];sim=np.asarray(a[start:start+256],float).copy();rows.extend(neighbours(sim,sub,ix,genes))
 else:raise ValueError(method)
 t=pd.DataFrame(rows,columns=['a','b','weight'])
 a=np.minimum(t.a,t.b);b=np.maximum(t.a,t.b);t['a']=a;t['b']=b
 t=t.groupby(['a','b'],sort=True).weight.max().reset_index();t.to_csv(dest,sep='\t',index=False)
 write(out/'edge_manifest.json',dict(method=method,events=p,edges=len(t),neighbours=30,
  input_sha256=sha(input/'input.bin'),scope='cross-gene edges; one nearest event per other gene; symmetric union',
  construction='posterior factor-implied positive covariance' if method=='factor_adaptive' else 'actual native WGCNA blockwise signed TOM'))
 return t

def partition(input,edges,out):
 input=Path(input);out=Path(out);events=pd.read_csv(input/'events.tsv',sep='\t');genes=events.ensembl_gene_id.to_numpy();p=len(events)
 graph=ig.Graph(n=p,edges=list(zip(edges.a,edges.b)),directed=False);weights=edges.weight.to_list()
 rows=[];all_labels={}
 for gamma in PARAMETERS['gamma_grid']:
  part=la.find_partition(graph,la.RBConfigurationVertexPartition,weights=weights,resolution_parameter=gamma,seed=PARAMETERS['seed'],n_iterations=4)
  raw=np.asarray(part.membership);labels=np.zeros(p,int);accepted=[];large=0;small=0
  for m in np.unique(raw):
   ix=np.flatnonzero(raw==m);ng=len(set(genes[ix]))
   if len(ix)<20 or ng<5:small+=len(ix);continue
   if ng>300:large+=len(ix);continue
   accepted.append(ix)
  accepted.sort(key=lambda ix:(-len(ix),str(events.event_key.iloc[ix].min())))
  for i,ix in enumerate(accepted,1):labels[ix]=i
  sz=np.array([len(ix) for ix in accepted]);gs=np.array([len(set(genes[ix])) for ix in accepted]);assigned=int(sz.sum())
  row=dict(gamma=gamma,modules=len(sz),assigned_events=assigned,coverage=assigned/p,
   assigned_genes=len(set(genes[labels>0])),gene_coverage=len(set(genes[labels>0]))/len(set(genes)),
   median_events=float(np.median(sz)) if len(sz) else None,p90_events=float(np.quantile(sz,.9)) if len(sz) else None,
   effective_events=float((sz*sz).sum()/sz.sum()) if assigned else None,
   median_genes=float(np.median(gs)) if len(gs) else None,p90_genes=float(np.quantile(gs,.9)) if len(gs) else None,
   max_genes=int(gs.max()) if len(gs) else 0,large_global_only_events=large,small_unassigned_events=small)
  rows.append(row);all_labels[gamma]=labels
 curve=pd.DataFrame(rows);curve.to_csv(out/'resolution_curve.tsv',sep='\t',index=False)
 configs=[]
 for target in PARAMETERS['target_effective_events']:
  eligible=curve[curve.coverage.ge(.7)&curve.effective_events.notna()].copy();fallback=False
  if eligible.empty:eligible=curve[curve.effective_events.notna()].copy();fallback=True
  if eligible.empty:raise ValueError('No eligible cross-gene local module at any resolution')
  eligible['distance']=np.abs(np.log(eligible.effective_events/target))
  if fallback:eligible=eligible.sort_values(['coverage','distance','gamma'],ascending=[False,True,True])
  else:eligible=eligible.sort_values(['distance','coverage','gamma'],ascending=[True,False,True])
  pick=eligible.iloc[0];labels=all_labels[float(pick.gamma)]
  dest=out/f'target{target}';dest.mkdir(exist_ok=True)
  pd.DataFrame({'event_key':events.event_key,'module':labels}).to_csv(dest/'assignments.tsv',sep='\t',index=False)
  rec=dict(target_effective_events=target,**pick.drop('distance').to_dict(),coverage_requirement_met=not fallback)
  write(dest/'selection.json',rec);configs.append(rec)
 write(out/'partition_manifest.json',dict(parameters=PARAMETERS,selection_scope='training shapes only; no outer recovery, time association or pathway outcomes',selected=configs))
 return configs

def run(method,input,model,out):
 edges=make_edges(method,input,model,out);return hierarchical_partition(input,edges,Path(out)/'hierarchical')

def hierarchical_partition(input,edges,out):
 """Refine only oversized communities; preserve already-resolved branches.

 A global resolution increase fragments every community simultaneously. This
 common adapter instead applies Leiden to each oversized induced subgraph.
 It is applied identically to both native-TOM and latent-profile graphs.
 """
 input=Path(input);out=Path(out);out.mkdir(parents=True,exist_ok=True)
 events=pd.read_csv(input/'events.tsv',sep='\t');genes=events.ensembl_gene_id.to_numpy();p=len(events)
 graph=ig.Graph(n=p,edges=list(zip(edges.a,edges.b)),directed=False);graph.es['weight']=edges.weight.to_list();graph.vs['event_index']=list(range(p))
 cache={};history=[]
 def split(ix):
  key=tuple(ix)
  if key in cache:return cache[key]
  sub=graph.induced_subgraph(ix);indices=np.array(sub.vs['event_index']);children=[]
  for gamma in [1.,2.,4.,8.,16.]:
   part=la.find_partition(sub,la.RBConfigurationVertexPartition,weights='weight',resolution_parameter=gamma,seed=32041,n_iterations=4)
   children=[np.sort(indices[np.asarray(m,int)]) for m in part]
   if len(children)>1:break
  if len(children)<=1:children=[]
  cache[key]=children;history.append(dict(parent_size=len(ix),parent_genes=len(set(genes[ix])),gamma=gamma,child_sizes=[len(c) for c in children]))
  return children
 rows=[];all_labels={}
 for bound in [100,150,200,300,400,600,800]:
  pending=[np.arange(p)];accepted=[];small=0;unresolved=0
  while pending:
   ix=pending.pop();ng=len(set(genes[ix]))
   if len(ix)<20 or ng<5:small+=len(ix);continue
   if len(ix)<=bound and ng<=300:accepted.append(ix);continue
   children=split(ix)
   if not children:unresolved+=len(ix)
   else:pending.extend(children)
  accepted.sort(key=lambda ix:(-len(ix),str(events.event_key.iloc[ix].min())))
  labels=np.zeros(p,int)
  for i,ix in enumerate(accepted,1):labels[ix]=i
  sz=np.array([len(ix) for ix in accepted]);gs=np.array([len(set(genes[ix])) for ix in accepted]);assigned=int(sz.sum())
  row=dict(max_leaf_events=bound,modules=len(sz),assigned_events=assigned,coverage=assigned/p,
   assigned_genes=len(set(genes[labels>0])),gene_coverage=len(set(genes[labels>0]))/len(set(genes)),
   median_events=float(np.median(sz)) if len(sz) else None,p90_events=float(np.quantile(sz,.9)) if len(sz) else None,
   effective_events=float((sz*sz).sum()/sz.sum()) if assigned else None,
   median_genes=float(np.median(gs)) if len(gs) else None,p90_genes=float(np.quantile(gs,.9)) if len(gs) else None,
   max_genes=int(gs.max()) if len(gs) else 0,large_global_only_events=unresolved,small_unassigned_events=small)
  rows.append(row);all_labels[bound]=labels
 curve=pd.DataFrame(rows);curve.to_csv(out/'resolution_curve.tsv',sep='\t',index=False);configs=[]
 for target in PARAMETERS['target_effective_events']:
  eligible=curve[curve.coverage.ge(.7)&curve.effective_events.notna()].copy();fallback=False
  if eligible.empty:eligible=curve[curve.effective_events.notna()].copy();fallback=True
  if eligible.empty:raise ValueError('No local cross-gene modules')
  eligible['distance']=np.abs(np.log(eligible.effective_events/target))
  eligible=eligible.sort_values(['coverage','distance'] if fallback else ['distance','coverage'],ascending=[False,True] if fallback else [True,False])
  pick=eligible.iloc[0];labels=all_labels[int(pick.max_leaf_events)];dest=out/f'target{target}';dest.mkdir(exist_ok=True)
  pd.DataFrame({'event_key':events.event_key,'module':labels}).to_csv(dest/'assignments.tsv',sep='\t',index=False)
  rec=dict(target_effective_events=target,**pick.drop('distance').to_dict(),coverage_requirement_met=not fallback)
  write(dest/'selection.json',rec);configs.append(rec)
 write(out/'partition_manifest.json',dict(parameters=PARAMETERS,rule='recursive refinement of oversized induced subgraphs, with same rule for both methods',
    selection_scope='training shapes only; no stability, holdout or pathway endpoint used',selected=configs,hierarchy=history))
 return configs
if __name__=='__main__':
 import argparse
 from common32 import EVID
 ap=argparse.ArgumentParser();ap.add_argument('method',choices=['factor_adaptive','wgcna_tom']);ap.add_argument('split',default='full',nargs='?');a=ap.parse_args()
 native='ebmf_counts' if a.method=='factor_adaptive' else 'wgcna'
 print(run(a.method,EVID/'inputs'/a.split,EVID/'development'/native/a.split,EVID/'development'/a.method/a.split/'graph'),flush=True)
