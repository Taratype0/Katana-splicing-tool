"""Motif-blind maximum-cardinality/minimum-cost matching; every event retained."""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching, min_weight_full_bipartite_matching

CALIPERS = np.array([.05, .20, .15, 1.0])
FEATURES = ['GC', 'length', 'baseline_PSI', 'depth']

def candidate_graph(foreground, background):
    required = ['event_id', 'gene_id', 'event_type', 'ambiguous_fraction', *FEATURES]
    for frame in [foreground, background]:
        if not set(required) <= set(frame): raise ValueError('Missing matching features')
        if frame.event_id.duplicated().any(): raise ValueError('Duplicate event identity')
        if frame.gene_id.isna().any() or frame.gene_id.astype(str).str.strip().eq('').any():
            raise ValueError('Gene identity required to define dependence')
    if set(foreground.gene_id) & set(background.gene_id): raise ValueError('Foreground/background gene overlap')
    fg = foreground.sort_values('event_id').reset_index(drop=True).copy()
    bg = background.sort_values('event_id').reset_index(drop=True).copy()
    genes = np.array(sorted(bg.gene_id.astype(str).unique()))
    gene_index = {g:i for i,g in enumerate(genes)}
    bgcodes = bg.gene_id.astype(str).map(gene_index).to_numpy()
    values = bg[FEATURES].to_numpy(float)
    kinds = bg.event_type.to_numpy()
    valid = np.isfinite(values).all(axis=1) & (values[:,1]>0) & (values[:,3]>0) & (bg.ambiguous_fraction.to_numpy()<=.05)
    edges=[]; counts=[]
    for i, row in enumerate(fg.itertuples()):
        x=np.array([getattr(row,c) for c in FEATURES],float)
        if not np.isfinite(x).all() or x[1]<=0 or x[3]<=0 or row.ambiguous_fraction>.05:
            counts.append(0); continue
        with np.errstate(divide='ignore',invalid='ignore'):
            distances=np.column_stack([abs(values[:,0]-x[0]),abs(values[:,1]-x[1])/x[1],abs(values[:,2]-x[2]),abs(np.log2(values[:,3]/x[3]))])
        idx=np.flatnonzero(valid & (kinds==row.event_type) & np.isfinite(distances).all(axis=1) & (distances<=CALIPERS).all(axis=1))
        counts.append(len(idx))
        if not len(idx):continue
        cost=1+np.square(distances[idx]/CALIPERS).sum(axis=1)
        # One edge per background gene; original event ID settles exact cost ties.
        order=np.lexsort((idx,cost,bgcodes[idx])); previous=-1
        for k in order:
            j=idx[k]; gene=bgcodes[j]
            if gene==previous:continue
            previous=gene
            edges.append((i,int(gene),int(j),float(cost[k]),*distances[j]))
    columns=['foreground_index','background_gene_index','background_index','cost','GC_distance','length_distance','PSI_distance','depth_distance']
    e=pd.DataFrame(edges,columns=columns)
    if len(e):e[['foreground_index','background_gene_index','background_index']]=e[['foreground_index','background_gene_index','background_index']].astype(int)
    fg['candidate_background_events']=counts
    return fg,bg,genes,e

def assign_graph(fg,bg,genes,edges):
    n,g=len(fg),len(genes)
    if n==0:raise ValueError('Empty foreground')
    ri=edges.foreground_index.to_numpy(int); ci=edges.background_gene_index.to_numpy(int)
    graph=csr_matrix((np.ones(len(edges)),(ri,ci)),shape=(n,g))
    optimum=int((maximum_bipartite_matching(graph,perm_type='column')>=0).sum())
    # Real edge cost is in [1,5]. This dominates every possible sum of real costs.
    penalty=5.0*n+1
    mandatory=fg.get('retain_original_match',pd.Series(False,index=fg.index)).to_numpy(bool)
    dummy_rows=np.flatnonzero(~mandatory)
    costs=csr_matrix((np.r_[edges.cost.to_numpy(float),np.full(len(dummy_rows),penalty)],(np.r_[ri,dummy_rows],np.r_[ci,g+dummy_rows])),shape=(n,g+n))
    rows, cols=min_weight_full_bipartite_matching(costs)
    lookup={(int(r.foreground_index),int(r.background_gene_index)):r for r in edges.itertuples()}
    assignments={}
    for i,j in zip(rows,cols):
        if j<g:assignments[int(i)]=lookup[(int(i),int(j))]
    assert len(assignments)==optimum
    assert set(np.flatnonzero(mandatory))<=set(assignments)
    result=[]
    for i,row in enumerate(fg.itertuples()):
        rec={'foreground_event':row.event_id,'foreground_gene':row.gene_id,'event_type':row.event_type,
             'direction':getattr(row,'direction',''),'status':'unmatched','background_event':'','background_gene':'',
             'candidate_background_events':row.candidate_background_events,'match_cost':np.nan}
        if i in assignments:
            edge=assignments[i];control=bg.iloc[int(edge.background_index)]
            rec.update(status='matched',background_event=control.event_id,background_gene=control.gene_id,match_cost=edge.cost)
            for name in ['GC_distance','length_distance','PSI_distance','depth_distance']:rec[name]=getattr(edge,name)
            rec['reason']='globally_assigned'
        else:rec['reason']='no_compatible_background' if row.candidate_background_events==0 else 'background_gene_capacity_exhausted'
        result.append(rec)
    result=pd.DataFrame(result)
    matched=result[result.status.eq('matched')]
    assert not matched.background_gene.duplicated().any()
    return result,{'maximum_cardinality_certificate':optimum,'candidate_gene_edges':len(edges),'dummy_penalty':penalty,'total_matched_cost':float(matched.match_cost.sum()),'protected_original_events':int(mandatory.sum())}

def global_match(foreground,background):
    fg,bg,genes,edges=candidate_graph(foreground,background)
    result,certificate=assign_graph(fg,bg,genes,edges)
    return result,certificate,edges

def greedy_all_events(fg,bg):
    """Ablation: original lexicographic distances/order, remove only foreground cap."""
    rows=[]
    values=bg[FEATURES].to_numpy(float);genes=bg.gene_id.to_numpy();kinds=bg.event_type.to_numpy();ids=bg.event_id.to_numpy()
    available=np.ones(len(bg),bool)
    for row in fg.sort_values('event_id').itertuples():
        x=np.array([getattr(row,c) for c in FEATURES]);chosen=None
        with np.errstate(divide='ignore',invalid='ignore'):
            dist=np.array([abs(values[:,0]-x[0]),abs(values[:,1]-x[1])/x[1],abs(values[:,2]-x[2]),abs(np.log2(values[:,3]/x[3]))])
        mask=(kinds==row.event_type)&(bg.ambiguous_fraction.to_numpy()<=.05)&np.isfinite(dist).all(axis=0)&(dist<=CALIPERS[:,None]).all(axis=0)&available
        if np.isfinite(x).all() and x[1]>0 and x[3]>0 and row.ambiguous_fraction<=.05:
            idx=np.flatnonzero(mask)
            if len(idx):chosen=idx[np.lexsort((ids[idx],dist[3,idx],dist[2,idx],dist[1,idx],dist[0,idx]))[0]]
        rec=dict(foreground_event=row.event_id,foreground_gene=row.gene_id,event_type=row.event_type,direction=row.direction,status='unmatched',background_event='',background_gene='')
        if chosen is not None:
            available[genes==genes[chosen]]=False;rec.update(status='matched',background_event=ids[chosen],background_gene=genes[chosen])
        rows.append(rec)
    return pd.DataFrame(rows)
