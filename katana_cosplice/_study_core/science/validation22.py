"""Frozen partition adapters, exact response retrieval and outer-refit metrics."""
from pathlib import Path
import json,time,traceback
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score
from cosplice_research import fit,write_json,within_module_support
from stable_network import fit_stable
from advanced import nearest_edges,graph_outputs

METHODS=['consensus','leiden','wgcna','wgcna_enhanced','factor','megena']

def fit_partition(method,measured,complete,meta,out):
    out=Path(out)
    if (out/'status.json').exists():
        status=json.loads((out/'status.json').read_text())
        if status['status'] in ['completed','failed']:return status
        raise RuntimeError(f'Incomplete directory requires inspection: {out}')
    try:
        if method in ['consensus','leiden']:
            return fit_stable(measured,complete,meta,out,consensus=method=='consensus')
        result=fit(method,complete,meta,out,timeout_seconds=900 if method=='megena' else 3600)
        if result['status']!='completed':return result
        if method=='factor':
            d=np.load(out/'factor.npz');w=d['loadings'];w=w/np.sqrt(np.sum(w*w,axis=1)+d['noise'])[:,None]
            edges=nearest_edges(lambda a,b:w[a:b]@w.T,complete.index)
            # Keep native factor scores; display module activity is separate.
            (out/'scores.tsv').rename(out/'native_factor_scores.tsv')
            graph_outputs(edges,complete,out,mutual=True);edges.to_csv(out/'edges.tsv',sep='\t',index=False)
            result['partition_adapter']='Existing mutual top20 absolute loading covariance graph; native Louvain seed31, minimum5'
        if method=='megena':
            tab=pd.read_csv(out/'memberships.tsv',sep='\t',dtype=str)
            tab['size']=tab.groupby('module').event_key.transform('size')
            selected=tab.sort_values(['size','module']).drop_duplicates('event_key').set_index('event_key').module.reindex(complete.index)
            unique=sorted(selected.dropna().unique());lookup={v:i+1 for i,v in enumerate(unique)}
            pd.DataFrame({'event_key':complete.index,'module':selected.map(lookup).fillna(0).astype(int)}).to_csv(out/'assignments.tsv',sep='\t',index=False)
            write_json(out/'module_names.json',lookup);result['partition_adapter']='Smallest native hierarchy module per event; ties lexical; full native hierarchy retained'
        write_json(out/'status.json',result);return result
    except Exception as exc:
        out.mkdir(parents=True,exist_ok=True);(out/'failure.txt').write_text(traceback.format_exc())
        result={'status':'failed','method':method,'error':str(exc)};write_json(out/'status.json',result);return result

def labels_at(folder,keys):
    return pd.read_csv(Path(folder)/'assignments.tsv',sep='\t').set_index('event_key').reindex(keys).module.fillna(0).astype(int).to_numpy()

def choose2(x):
    x=np.asarray(x,dtype=np.int64);return int(np.sum(x*(x-1)//2))

def cross_pairs(genes):
    return len(genes)*(len(genes)-1)//2-choose2(pd.Series(genes).value_counts().values)

def response_retrieval(labels,genes,response,threshold=10):
    """Exact two-score AP over all evaluable cross-gene event pairs. Unmarked != negative truth."""
    tab=pd.DataFrame({'m':labels,'g':genes,'v':response}).dropna()
    N=cross_pairs(tab.g);positive=tp=predicted=0
    for _,group in tab[tab.m.ne(0)].groupby('m',sort=False):predicted+=cross_pairs(group.g)
    for selected in [tab[tab.v.ge(threshold)],tab[tab.v.le(-threshold)]]:
        positive+=cross_pairs(selected.g)
        for _,group in selected[selected.m.ne(0)].groupby('m',sort=False):tp+=cross_pairs(group.g)
    prevalence=positive/N if N else np.nan
    recall=tp/positive if positive else np.nan
    precision=tp/predicted if predicted else 0.
    ap=recall*precision+(1-recall)*prevalence if positive else np.nan
    return dict(evaluable_events=len(tab),candidate_pairs=N,known_positive_pairs=positive,module_pairs=predicted,
        recovered_positive_pairs=tp,prevalence=prevalence,precision=precision,recall=recall,average_precision=ap,
        ap_lift=ap/prevalence if prevalence else np.nan)

def singleton_grey(labels):
    z=np.asarray(labels,dtype=np.int64).copy();grey=z==0;z[grey]=-1-np.flatnonzero(grey);return z

def partition_metrics(base,lab):
    support=within_module_support(base,lab);sizes=pd.Series(lab[lab!=0]).value_counts().values
    chance=choose2(sizes)/(len(lab)*(len(lab)-1)/2)
    pair_base=choose2(pd.Series(base[base!=0]).value_counts().values);pair_refit=choose2(sizes)
    joint=pd.DataFrame({'a':base,'b':lab});joint=joint[(joint.a!=0)&(joint.b!=0)]
    common=choose2(joint.groupby(['a','b']).size().values)
    jacc=[];den=pd.Series(lab).value_counts()
    for m in sorted(set(base)-{0}):
        ids=base==m;counts=pd.Series(lab[ids]).value_counts()
        jacc.append(max([v/(ids.sum()+den[k]-v) for k,v in counts.items() if k!=0] or [0]))
    # Singleton grey vertices contribute no coassigned pair. This is algebraically
    # identical to sklearn ARI on singleton_grey(), without treating many classes
    # as a regression target. A baseline with no module has no stability endpoint.
    expected=pair_base*pair_refit/(len(base)*(len(base)-1)/2);maximum=(pair_base+pair_refit)/2
    ari=(common-expected)/(maximum-expected) if maximum!=expected else 1.
    return support,dict(chance_coassignment=chance,mean_chance_adjusted_support=float(np.nanmean((support-chance)/max(1-chance,1e-12))) if np.isfinite(support).any() else None,
        pair_f1=2*common/(pair_base+pair_refit) if pair_base else np.nan,
        singleton_grey_ari=ari if pair_base else np.nan,
        module_match_jaccard=float(np.mean(jacc)) if jacc else np.nan)

def make_draws(meta,iterations=20,seed=22041):
    groups=meta.groupby('biological_sample_id',sort=False).indices;strata={}
    for u,c in groups.items():strata.setdefault(tuple(sorted(meta.iloc[c].condition)),[]).append(u)
    rng=np.random.default_rng(seed)
    return [[str(u) for units in strata.values() for u in rng.choice(units,len(units),replace=True)] for _ in range(iterations)]

def outer_refit(method,measured,meta,units,out):
    out=Path(out);cached=out/'labels.npz'
    if cached.exists():return np.load(cached)['labels'],json.loads((out/'outer_status.json').read_text())
    groups=meta.groupby('biological_sample_id',sort=False).indices;cols=np.concatenate([groups[u] for u in units])
    x=measured.iloc[:,cols].copy();x.columns=[f'b{i}' for i in range(len(cols))]
    md=meta.iloc[cols].reset_index(drop=True).copy();md.sample_id=x.columns
    md.biological_sample_id=[f'unit{j}' for j,u in enumerate(units) for _ in groups[u]]
    estimable=x.mean(axis=1).notna()&x.var(axis=1).gt(1e-12);sub=x.loc[estimable];complete=sub.T.fillna(sub.mean(axis=1)).T
    result=fit_partition(method,sub,complete,md,out/'model')
    lab=labels_at(out/'model',measured.index) if result['status']=='completed' else np.zeros(len(measured),dtype=int)
    result.update(drawn_units=units,unestimable_events=int((~estimable).sum()),eligible_events=len(sub))
    np.savez_compressed(cached,labels=lab,event_key=np.asarray(measured.index,str));write_json(out/'outer_status.json',result)
    return lab,result

def summarize_stability(base,keys,refits,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);supports=[];records=[]
    for i,lab,status in refits:
        s,metrics=partition_metrics(base,lab);supports.append(s);records.append(dict(replicate=i,**status,**metrics))
    a=np.asarray(supports);count=np.isfinite(a).sum(axis=0)
    avg=np.divide(np.nansum(a,axis=0),count,out=np.full(len(base),np.nan),where=count>0)
    success=sum(r['status']=='completed' for r in records)
    pd.DataFrame(dict(event_key=keys,module=base,within_module_coassignment=avg,requested_replicates=len(records),successful_replicates=success)).to_csv(out/'stability.tsv',sep='\t',index=False)
    write_json(out/'replicates.json',records);np.savez_compressed(out/'per_event_replicates.npz',support=a,event_key=np.asarray(keys,str))
    sizes=pd.Series(base[base!=0]).value_counts()
    summary=dict(replicates=len(records),successful=success,events=len(base),assigned_events=int((base!=0).sum()),modules=len(sizes),largest_module=int(sizes.max()) if len(sizes) else 0,
        event_support_median=float(np.nanmedian(avg)) if np.isfinite(avg).any() else None,event_support_mean=float(np.nanmean(avg)) if np.isfinite(avg).any() else None,
        coverage_weighted_support=float(np.nansum(avg)/len(base)))
    for k in ['mean_chance_adjusted_support','pair_f1','singleton_grey_ari','module_match_jaccard']:
        v=[r[k] for r in records if r[k] is not None];summary[k]=float(np.nanmean(v)) if v else None
    write_json(out/'summary.json',summary);return summary
