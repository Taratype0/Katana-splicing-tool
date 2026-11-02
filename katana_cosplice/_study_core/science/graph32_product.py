"""Production small-input correction; frozen large-cohort graph32 stays intact.

Always discover root communities before applying the size bound. For the
recount3 benchmark (>800 events in every fit), this has identical decisions.
"""
from graph32 import Path,json,np,pd,ig,la,PARAMETERS,write,make_edges

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
   if len(ix)<=bound and ng<=300 and len(ix)!=p:accepted.append(ix);continue
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

def run(method,input,model,out):
 edges=make_edges(method,input,model,out)
 return hierarchical_partition(input,edges,Path(out)/'hierarchical')
