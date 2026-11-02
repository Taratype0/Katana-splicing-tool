"""Full-QC unit-deletion refits and heldout biological replication endpoints."""
from pathlib import Path
import json,time
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from unit_robust import unit_audit,fit_jackknife
from cosplice_research import write_json
from validation22 import fit_partition,labels_at,partition_metrics
from study import score_modules

METHODS=['consensus','jackknife','wgcna_robust','wgcna_enhanced','factor','megena']

def fit_common(method,raw,meta,support,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    meta=meta.set_index('sample_id').loc[raw.columns].rename_axis('sample_id').reset_index()
    measured,audit=unit_audit(raw,meta,support);audit.to_csv(out/'event_audit.tsv.gz',sep='\t',index_label='event_key')
    x=measured.loc[audit.eligible]
    if method=='global_local':
        from peer30 import impute
        complete=impute(x)
    else:complete=x.T.fillna(x.mean(axis=1)).T
    if (out/'model/status.json').exists():result=json.loads((out/'model/status.json').read_text())
    elif len(x)<40:result={'status':'failed','method':method,'error':'Fewer than 40 common-QC events'}
    elif method=='jackknife':result=fit_jackknife(x,complete,meta,out/'model')
    elif method=='global_local':
        from algorithm31 import fit_global_local
        result=fit_global_local(x,complete,meta,out/'model')
    elif method=='factor_core':
        if meta.biological_sample_id.duplicated().any():raise ValueError('Factor core currently requires independent biological units')
        from algorithm28 import fit_latent
        result=fit_latent(method,x,complete,meta,out/'model')
    else:result=fit_partition(method,x,complete,meta,out/'model')
    result.update(eligible_events=len(x),input_events=len(raw),biological_units=int(meta.biological_sample_id.nunique()))
    write_json(out/'status.json',result)
    return result,x,meta

def deletion_draws(meta,iterations=20,seed=24031):
    groups=meta.groupby('biological_sample_id',sort=True).indices;strata={}
    for u,ix in groups.items():strata.setdefault(tuple(sorted(meta.iloc[ix].condition)),[]).append(u)
    rng=np.random.default_rng(seed);draws=[]
    for i in range(iterations):
        # Delete one actual unit per exchangeable condition profile. Ceil(.8*3)
        # formerly retained all three units and was not a perturbation at all.
        omit=[str(rng.choice(g)) for g in strata.values() if len(g)>1]
        if not omit:raise ValueError('No exchangeable repeated units for deletion validation')
        draws.append(meta.loc[~meta.biological_sample_id.isin(omit),'sample_id'].tolist())
    return draws

def summarize(base,keys,records,out):
    out=Path(out);out.mkdir(exist_ok=True,parents=True);supports=[];rows=[]
    for i,lab,status in records:
        if status['status']=='completed':s,m=partition_metrics(base,lab)
        else:s=np.full(len(base),np.nan);m={k:np.nan for k in ['pair_f1','singleton_grey_ari','module_match_jaccard','mean_chance_adjusted_support']}
        supports.append(s);rows.append(dict(replicate=i,**status,**m))
    values=np.asarray(supports);n=np.isfinite(values).sum(axis=0)
    avg=np.divide(np.nansum(values,axis=0),n,out=np.full(len(base),np.nan),where=n>0)
    success=sum(r['status']=='completed' for r in rows)
    pd.DataFrame({'event_key':keys,'module':base,'within_module_coassignment':avg,'requested_replicates':len(records),'successful_replicates':success}).to_csv(out/'stability.tsv',sep='\t',index=False)
    result={'requested':len(records),'successful':success,'events':len(base),'assigned_events':int((base!=0).sum()),'modules':len(set(base)-{0}),
        'event_support_median':float(np.nanmedian(avg)) if np.isfinite(avg).any() else None,
        'event_support_mean':float(np.nanmean(avg)) if np.isfinite(avg).any() else None,
        'coverage_weighted_support':float(np.nansum(avg)/len(base)) if success else None,
        'scope':'Condition-stratified biological-unit deletion; full event QC, imputation, graph and partition refitted; failures are NA, not zero support'}
    for k in ['pair_f1','singleton_grey_ari','module_match_jaccard','mean_chance_adjusted_support']:
        v=[r[k] for r in rows if np.isfinite(r[k])];result[k]=float(np.mean(v)) if v else None
    write_json(out/'summary.json',result);write_json(out/'replicates.json',rows)
    np.savez_compressed(out/'per_event_replicates.npz',support=values,event_key=np.asarray(keys,str));return result

def product_stability(method,raw,meta,support,assignment,keys,out,iterations):
    seed=28071 if method=='factor_core' else 24031
    out=Path(out);out.mkdir(exist_ok=False)
    if method=='global_local':
        from algorithm31 import deletion_draws as global_local_draws
        seed=31071;draws=global_local_draws(meta,iterations,seed)
    else:draws=deletion_draws(meta,iterations,seed)
    write_json(out/'planned_draws.json',{'training_sample_draws':draws,'seed':seed});records=[]
    for i,cols in enumerate(draws):
        target=out/f'refit_{i:03}';result,_,_=fit_common(method,raw.loc[:,cols],meta,support.loc[:,cols] if support is not None else None,target)
        labels=labels_at(target/'model',keys) if result['status']=='completed' else None;records.append((i,labels,result))
    return summarize(assignment.reindex(keys).module.to_numpy(),keys,records,out)

def pathway_sets(path,background):
    sets={}
    for line in Path(path).read_text(encoding='utf-8').splitlines():
        name,url,*genes=line.split('\t');g=set(genes)&background
        if 10<=len(g)<=500:sets[name]=(url,g)
    return sets

def single_gene(value):
    return isinstance(value,str) and value not in ['nan','NA','UNKNOWN',''] and not any(x in value for x in ';,|')

def enrich(assignment,genes,background,pathways,out):
    from statsmodels.stats.multitest import multipletests
    sets=pathway_sets(pathways,background);rows=[]
    lab=assignment.module.reindex(genes.index);known=genes[genes.map(single_gene)]
    for m in sorted(set(lab.dropna())-{0}):
        g=set(known[lab.eq(m)])&background
        if len(g)<5:continue
        for name,(url,term) in sets.items():
            hits=g&term;rows.append({'module':int(m),'pathway':name,'url':url,'module_genes':len(g),'background_genes':len(background),'pathway_background_genes':len(term),'overlap':len(hits),'genes':';'.join(sorted(hits)),
                'p_value':float(hypergeom.sf(len(hits)-1,len(background),len(term),len(g)))})
    d=pd.DataFrame(rows,columns=['module','pathway','url','module_genes','background_genes','pathway_background_genes','overlap','genes','p_value']);d['q_value']=multipletests(d.p_value,method='fdr_bh')[1] if len(d) else []
    d.to_csv(Path(out)/'pathway_enrichment.tsv',sep='\t',index=False);sig=d[d.q_value.lt(.05)&d.overlap.ge(3)]
    return {'tests':len(d),'tested_modules':int(d.module.nunique()),'significant_modules':int(sig.module.nunique()),'significant_pathways':int(sig.pathway.nunique()),'background_genes':len(background),'gene_set_collection':'MSigDB Reactome 2025.1.Hs','interpretation':'Curated functional association, not causal regulatory truth; same gene counted once per module'}

def pair_correlations(values,a,b):
    left=values[a];right=values[b];ok=np.isfinite(left)&np.isfinite(right);n=ok.sum(axis=1)
    left=np.where(ok,left,0);right=np.where(ok,right,0)
    left=np.where(ok,left-left.sum(axis=1)[:,None]/np.maximum(n[:,None],1),0)
    right=np.where(ok,right-right.sum(axis=1)[:,None]/np.maximum(n[:,None],1),0)
    den=np.linalg.norm(left,axis=1)*np.linalg.norm(right,axis=1)
    return np.divide((left*right).sum(axis=1),den,out=np.full(len(a),np.nan),where=(n>=6)&(den>1e-12))

def heldout_metrics(train,heldout,meta,assignment,genes,out):
    scores,_,modules,basis=score_modules(train,assignment);rng=np.random.default_rng(24051);rows=[]
    for name,b in basis.groupby('module'):
        keys=b.event_key.tolist();g=genes.reindex(keys);known=g.map(single_gene).to_numpy();p=len(keys)
        if p*(p-1)//2<=20000:a,c=np.triu_indices(p,1)
        else:
            pairs=set()
            while len(pairs)<20000:
                aa=rng.integers(p,size=25000);bb=rng.integers(p,size=25000);pairs.update((int(min(i,j)),int(max(i,j))) for i,j in zip(aa,bb) if i!=j)
            pairs=sorted(pairs);chosen=rng.choice(len(pairs),20000,replace=False);a,c=np.array([pairs[i] for i in chosen]).T
        cross=known[a]&known[c]&(g.to_numpy()[a]!=g.to_numpy()[c]);a=a[cross];c=c[cross]
        query=heldout.reindex(keys);r=pair_correlations(query.to_numpy(),a,c);orient=np.sign(b.loading.to_numpy());oriented=r*orient[a]*orient[c]
        z=((query.T.fillna(b.set_index('event_key').center)-b.set_index('event_key').center)/b.set_index('event_key').scale)
        projected=z.to_numpy()@b.loading.to_numpy();coverage=query.notna().mean(axis=0).to_numpy();projected[coverage<.8]=np.nan
        md=meta.set_index('sample_id');cond_mean=scores[name].groupby(md.loc[scores.index,'condition']).mean()
        expected=md.loc[query.columns,'condition'].map(cond_mean).to_numpy();valid=np.isfinite(projected)&np.isfinite(expected)
        corr=float(np.corrcoef(expected[valid],projected[valid])[0,1]) if valid.sum()>=6 and np.std(expected[valid])>1e-12 and np.std(projected[valid])>1e-12 else np.nan
        rows.append({'module':name,'events':p,'sampled_cross_gene_pairs':len(a),'evaluable_pairs':int(np.isfinite(r).sum()),
            'mean_holdout_signed_correlation':float(np.nanmean(r)) if np.isfinite(r).any() else np.nan,
            'mean_holdout_PC1_oriented_correlation':float(np.nanmean(oriented)) if np.isfinite(oriented).any() else np.nan,
            'frozen_score_condition_replication_r':corr,'evaluable_holdout_samples':int(valid.sum())})
    d=pd.DataFrame(rows);d.to_csv(Path(out)/'heldout_modules.tsv',sep='\t',index=False)
    if not len(d):return {'modules':0,'evaluable_modules':0}
    return {'modules':len(d),'evaluable_modules':int(d.mean_holdout_PC1_oriented_correlation.notna().sum()),
        'median_heldout_coherence':float(d.mean_holdout_PC1_oriented_correlation.median()),'median_temporal_replication_r':float(d.frozen_score_condition_replication_r.median()),
        'modules_positive_replication':int(d.frozen_score_condition_replication_r.gt(0).sum()),'tested_temporal_modules':int(d.frozen_score_condition_replication_r.notna().sum()),
        'scope':'New culture replicates at the same eight conditions; full training-only QC and frozen PC1 basis; pairwise observed PSI, no heldout imputation for coherence'}
