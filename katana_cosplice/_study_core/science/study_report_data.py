"""Descriptive member annotations over existing modules, PSI and graph weights.

No network fitting or significance testing. kME uses pairwise observed PSI and
the frozen module score. Connectivity is a weighted degree, not confidence.
"""
import numpy as np
import pandas as pd


def member_evidence(events, usage, scores, samples, edges):
    events=events.copy()
    if events.empty:
        return events,pd.DataFrame(columns=['module','gene_id','events','max_abs_kME','within_module_weight']),{},{}
    events['gene_id']=events.get('gene_id',pd.Series(index=events.index,dtype=str)).fillna('Unannotated').astype(str)
    inferred=events.event_key.map(lambda x:'|'.join(str(x).split('|')[2:-1]) or 'Unspecified')
    events['event_type']=events.get('event_type',inferred).fillna(inferred)
    sample_ids=samples.sample_id.tolist();values=usage.reindex(index=events.event_key,columns=sample_ids).astype(float)
    scores=scores.set_index('sample_id').reindex(sample_ids)
    events['kME']=np.nan;events['kME_samples']=0
    for module,idx in events.groupby('module').groups.items():
        if module not in scores:continue
        keys=events.loc[idx,'event_key'];x=values.loc[keys].to_numpy();y=scores[module].to_numpy(float)
        valid=np.isfinite(x)&np.isfinite(y)[None,:];n=valid.sum(axis=1)
        av=np.divide(np.where(valid,x,0).sum(axis=1),n,out=np.zeros(len(n)),where=n>0)
        bv=np.divide(np.where(valid,y[None,:],0).sum(axis=1),n,out=np.zeros(len(n)),where=n>0)
        a=np.where(valid,x-av[:,None],0);b=np.where(valid,y[None,:]-bv[:,None],0)
        den=np.sqrt((a*a).sum(axis=1)*(b*b).sum(axis=1))
        cor=np.divide((a*b).sum(axis=1),den,out=np.full(len(n),np.nan),where=(n>=3)&(den>1e-12))
        events.loc[idx,'kME']=cor;events.loc[idx,'kME_samples']=n
    condition_means=[]
    for _,g in samples.groupby('condition',sort=False):condition_means.append(values[g.sample_id].mean(axis=1))
    means=pd.concat(condition_means,axis=1)
    events['condition_mean_psi_range']=(means.max(axis=1)-means.min(axis=1)).to_numpy()
    events['observed_samples']=values.notna().sum(axis=1).to_numpy()
    events['within_module_weight']=np.nan;graphs={}
    if {'source','target','weight'}.issubset(edges):
        labels=events.set_index('event_key').module
        e=edges[edges.source.isin(labels.index)&edges.target.isin(labels.index)].copy()
        e=e[e.source.ne(e.target)&e.source.map(labels).eq(e.target.map(labels))]
        # Scientific graph is undirected; count each observed edge exactly once.
        pair=e.apply(lambda r:tuple(sorted([r.source,r.target])),axis=1)
        e=e.loc[~pair.duplicated()].copy();e['module']=e.source.map(labels)
        degrees=pd.concat([e.set_index('source').weight,e.set_index('target').weight]).groupby(level=0).sum()
        events['within_module_weight']=events.event_key.map(degrees).fillna(0.)
        for module,g in events.groupby('module',sort=False):
            selected=g.sort_values(['within_module_weight','event_key'],ascending=[False,True]).head(40)
            keys=selected.event_key.tolist();positions={k:i for i,k in enumerate(keys)}
            full=e[e.module.eq(module)]
            view=full[full.source.isin(keys)&full.target.isin(keys)].sort_values('weight',ascending=False).head(120)
            graphs[module]={'nodes':keys,'edges':[[positions[r.source],positions[r.target],float(r.weight)] for r in view.itertuples()],
                            'total_nodes':len(g),'total_edges':len(full),'rule':'At most 40 events by within-module weighted degree; at most 120 strongest edges among displayed nodes'}
    genes=events.assign(abs_kME=events.kME.abs()).groupby(['module','gene_id'],sort=False,dropna=False).agg(
        events=('event_key','size'),max_abs_kME=('abs_kME','max'),within_module_weight=('within_module_weight',lambda x:x.sum(min_count=1))).reset_index()
    genes=genes.sort_values(['module','within_module_weight','gene_id'],ascending=[True,False,True],na_position='last')
    return events,genes,dict(zip(events.event_key,values.to_numpy().tolist())),graphs
