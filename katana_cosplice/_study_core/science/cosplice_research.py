"""Independent CoSplice research CLI. No fixed event cap; native R engines."""
from __future__ import annotations
import argparse, hashlib, json, os, subprocess, time, warnings,shutil
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.decomposition import FactorAnalysis,PCA

ROOT=Path(__file__).resolve().parents[1]
def rscript_path():
    selected = os.environ.get('KATANA_RSCRIPT') or shutil.which('Rscript')
    if not selected or not Path(selected).is_file():
        raise RuntimeError('Rscript is unavailable; install the Katana native environment or select KATANA_RSCRIPT')
    return str(Path(selected).absolute())

def write_json(path,obj):
    path=Path(path);temp=path.with_name(path.name+'.tmp')
    temp.write_text(json.dumps(obj,indent=2,ensure_ascii=False,default=str),encoding='utf-8');temp.replace(path)

def prepare(matrix, metadata, *, value_kind='usage'):
    if matrix.empty or not matrix.index.is_unique or not matrix.columns.is_unique:
        raise ValueError('Nonempty matrix with unique event/sample IDs required')
    required={'sample_id','biological_sample_id','condition'}
    if not required.issubset(metadata): raise ValueError('Declare sample_id, biological_sample_id, condition')
    if metadata[list(required)].isna().any().any() or not metadata.sample_id.is_unique:
        raise ValueError('Missing or duplicate sample metadata')
    if set(matrix.columns)!=set(metadata.sample_id): raise ValueError('Matrix and sample metadata must match exactly')
    x=matrix.astype(float); vals=x.to_numpy()
    if value_kind not in ['usage','response']:raise ValueError('Unknown measurement kind')
    low,high=(0,1) if value_kind=='usage' else (-100,100)
    if np.isinf(vals).any() or ((vals[np.isfinite(vals)]<low)|(vals[np.isfinite(vals)]>high)).any():
        raise ValueError('Observed usage must lie in [0,1]' if value_kind=='usage' else 'Observed delta PSI percentage points must lie in [-100,100]')
    meta=metadata.set_index('sample_id').loc[x.columns].rename_axis('sample_id').reset_index()
    # Technical replicates contribute one observation per unit and condition.
    if meta.duplicated(['biological_sample_id','condition']).any():
        raise ValueError('Collapse technical replicates per biological unit/condition before this research entry')
    valid=x.notna().mean(axis=1);var=x.var(axis=1)
    audit=pd.DataFrame({'observed_fraction':valid,'variance':var})
    audit['eligible']=(valid>=.8)&(var>1e-12)
    audit['reason']=[';'.join(r for ok,r in [(a>=.8,'missingness'),(b>1e-12,'constant')] if not ok) or 'passed'
                     for a,b in zip(valid,var)]
    measured=x.loc[audit.eligible];means=measured.mean(axis=1)
    complete=measured.T.fillna(means).T
    return measured,complete,audit,means,meta

def fit(method, complete, meta, out,timeout_seconds=None):
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    complete.to_csv(out/'input.tsv',sep='\t',index_label='event_key')
    meta.to_csv(out/'samples.tsv',sep='\t',index=False)
    started=time.monotonic(); status={'method':method,'events':len(complete),'samples':complete.shape[1]}
    try:
        if method in ['joint','independent'] and meta.biological_sample_id.duplicated().any():
            raise ValueError('Joint-network benchmark requires independent units across conditions; paired estimation not implemented')
        if method in ['factor','pca']:
            mean=complete.mean(axis=1).to_numpy();scale=complete.std(axis=1).clip(lower=1e-8).to_numpy()
            z=(complete.to_numpy().T-mean)/scale
            k=max(1,min(8,z.shape[0]//5,z.shape[1]-1))
            model=FactorAnalysis(n_components=k,rotation='varimax',random_state=31,max_iter=1000) if method=='factor' else PCA(n_components=k,random_state=31)
            with warnings.catch_warnings(record=True) as caught:
                scores=model.fit_transform(z)
            status['warnings']=[str(w.message) for w in caught]
            status['iterations']=int(getattr(model,'n_iter_',0))
            loadings=model.components_.T if method=='factor' else model.components_.T*np.sqrt(model.explained_variance_)
            np.savez_compressed(out/'factor.npz',loadings=loadings,noise=np.broadcast_to(model.noise_variance_,len(mean)),mean=mean,scale=scale)
            pd.DataFrame(loadings,index=complete.index).to_csv(out/'loadings.tsv',sep='\t')
            pd.DataFrame(scores,index=complete.columns).to_csv(out/'scores.tsv',sep='\t')
            status['factors']=k
        elif method=='correlation':
            status['description']='Absolute Pearson, no module restriction; unthresholded descriptive baseline'
        else:
            threads=str(max(1,int(os.environ.get('KATANA_CPUS','2'))))
            env=os.environ.copy();env.update(OMP_NUM_THREADS=threads,OPENBLAS_NUM_THREADS=threads,MKL_NUM_THREADS=threads)
            engine=out/'engine_source.R'
            shutil.copyfile(ROOT/'science/engines.R',engine)
            with (out/'engine.log').open('w',encoding='utf-8') as log:
                proc=subprocess.Popen([rscript_path(),'--vanilla',str(engine),str(ROOT),method,
                    str(out/'input.tsv'),str(out),str(out/'samples.tsv')],stdout=log,stderr=subprocess.STDOUT,env=env,
                    creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
                try:
                    import psutil
                    peak=0
                    while proc.poll() is None:
                        if timeout_seconds is not None and time.monotonic()-started>timeout_seconds:raise subprocess.TimeoutExpired(proc.args,timeout_seconds)
                        try:
                            owned=psutil.Process(proc.pid);rss=sum(p.memory_info().rss for p in [owned]+owned.children(recursive=True) if p.is_running());peak=max(peak,rss)
                        except psutil.NoSuchProcess:pass
                        time.sleep(.5)
                    status['native_peak_rss_bytes']=peak
                except subprocess.TimeoutExpired:
                    # Windows Rscript is a wrapper; stopping only that process leaves
                    # its actual R child running beyond the declared budget.
                    import psutil
                    try:owned=psutil.Process(proc.pid);children=owned.children(recursive=True)
                    except psutil.NoSuchProcess:children=[]
                    for child in children:
                        try:child.terminate()
                        except psutil.NoSuchProcess:pass
                    proc.kill();proc.wait();psutil.wait_procs(children,timeout=10)
                    status['native_peak_rss_bytes']=peak
                    raise TimeoutError(f'Native {method} exceeded requested {timeout_seconds}s; owned R process tree stopped')
            if proc.returncode: raise RuntimeError(f'Native R exit {proc.returncode}; see engine.log')
        status['status']='completed'
    except Exception as exc:
        status.update(status='failed',error=str(exc))
    status['seconds']=time.monotonic()-started
    write_json(out/'status.json',status)
    return status

def edge_scores(method,out,keys,values):
    """Benchmark only: dense scores on the explicitly bounded common panel."""
    p=len(keys);out=Path(out)
    if p>2000: raise ValueError('Dense BENCHMARK scoring limited to <=2000; scalable discovery has no event cap')
    if method.startswith('wgcna'):
        labels=pd.read_csv(out/'assignments.tsv',sep='\t').set_index('event_key').loc[keys,'module'].to_numpy()
        corr=np.corrcoef(values)
        scores=np.maximum(corr,0)*((labels[:,None]==labels)&(labels[:,None]!=0))
    elif method in ['factor','pca']:
        a=np.load(out/'factor.npz');l=a['loadings'];cov=l@l.T
        sd=np.sqrt(np.diag(cov)+a['noise']);scores=np.abs(cov/sd[:,None]/sd)
    elif method=='correlation':
        scores=np.abs(np.corrcoef(values))
    elif method=='megena':
        tab=pd.read_csv(out/'memberships.tsv',sep='\t');scores=np.zeros((p,p));pos={str(x):i for i,x in enumerate(keys)}
        for _,g in tab.groupby('module'):
            ids=[pos[str(x)] for x in g.event_key if str(x) in pos]
            if len(ids)>1: scores[np.ix_(ids,ids)]=np.maximum(scores[np.ix_(ids,ids)],1/(len(ids)-1))
    else:
        scores=np.zeros((p,p));conds=pd.read_csv(out/'conditions.tsv',sep='\t')
        for k in conds['index']:
            theta=np.loadtxt(out/f'precision_{k}.tsv');sd=np.sqrt(np.diag(theta))
            scores+=np.abs(-theta/sd[:,None]/sd)/len(conds)
    np.fill_diagonal(scores,0)
    return np.nan_to_num(scores)

def within_module_support(baseline,refit):
    """Exact mean pairwise coassignment within frozen modules, O(p) storage."""
    base=np.asarray(baseline);lab=np.asarray(refit);out=np.full(len(base),np.nan)
    for module in np.unique(base):
        if module == 0 or str(module) in {'0','0.0','grey'}:continue
        ids=np.flatnonzero(base==module)
        if len(ids)<2:continue
        labels,counts=np.unique(lab[ids],return_counts=True);lookup=dict(zip(labels,counts))
        out[ids]=[(lookup[x]-1)/(len(ids)-1) if x != 0 and str(x) not in {'0','0.0','grey'} else 0. for x in lab[ids]]
    return out

def bootstrap(measured,meta,model_dir,out,iterations=50,seed=31,workers=4):
    """Resample biological blocks within condition profiles; no p-by-p matrix."""
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    baseline=pd.read_csv(Path(model_dir)/'assignments.tsv',sep='\t').set_index('event_key').loc[measured.index,'module'].to_numpy()
    meta=meta.set_index('sample_id').loc[measured.columns]
    unitcols={u:np.flatnonzero(meta.biological_sample_id.to_numpy()==u) for u in meta.biological_sample_id.unique()}
    strata={}
    for unit,cols in unitcols.items():strata.setdefault(tuple(sorted(meta.iloc[cols].condition)),[]).append(unit)
    from concurrent.futures import ThreadPoolExecutor
    rng=np.random.default_rng(seed);supports=[];records=[]
    draws=[[u for units in strata.values() for u in rng.choice(units,size=len(units),replace=True)] for _ in range(iterations)]
    write_json(out/'planned_draws.json',{'seed':seed,'workers':workers,'draws':draws})
    def refit(item):
        b,chosen=item
        cols=np.concatenate([unitcols[u] for u in chosen]);x=measured.iloc[:,cols].copy()
        x.columns=[f'b{b}_s{i}' for i in range(len(cols))]
        # Refit means within the draw. An unobservable/constant event contributes
        # zero support, without invalidating all other events in the same draw.
        mean=x.mean(axis=1)
        estimable=mean.notna()&x.var(axis=1).gt(1e-12)
        if estimable.sum()<40:
            return {'replicate':b,'status':'invalid_draw','drawn_units':chosen,'estimable_events':int(estimable.sum())},np.zeros(len(measured))
        x=x.loc[estimable].T.fillna(mean).T
        bm=meta.iloc[cols].reset_index(drop=True);bm['sample_id']=x.columns
        status=fit('wgcna',x,bm,out/f'refit_{b:03}')
        record={'replicate':b,'drawn_units':chosen,'unestimable_events':measured.index[~estimable].tolist(),**status}
        if status['status']=='completed':
            labels=pd.read_csv(out/f'refit_{b:03}/assignments.tsv',sep='\t').set_index('event_key').reindex(measured.index)['module'].fillna(0).to_numpy()
            support=within_module_support(baseline,labels)
        else:support=np.zeros(len(measured))
        return record,support
    with ThreadPoolExecutor(max_workers=workers) as executor:
        for record,support in executor.map(refit,enumerate(draws)):
            records.append(record);supports.append(support);write_json(out/'replicates.json',records)
            print(f'Biological-block refit {len(records)}/{iterations}: {record["status"]}',flush=True)
    write_json(out/'replicates.json',records)
    a=np.asarray(supports)
    count=np.isfinite(a).sum(axis=0)
    support_mean=np.divide(np.nansum(a,axis=0),count,out=np.full(len(baseline),np.nan),where=count>0)
    support_mean[np.isin(baseline.astype(str),['0','0.0','grey'])]=np.nan
    pd.DataFrame({'event_key':measured.index,'module':baseline,'within_module_coassignment':support_mean,
                  'requested_replicates':iterations,'successful_replicates':sum(r.get('status')=='completed' for r in records)}).to_csv(out/'stability.tsv',sep='\t',index=False)
    np.savez_compressed(out/'per_event_replicates.npz',support=a,event_key=measured.index.to_numpy(str))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--psi',required=True);parser.add_argument('--samples',required=True)
    parser.add_argument('--output',required=True);parser.add_argument('--method',default='wgcna',choices=['wgcna','megena','factor','joint','independent'])
    parser.add_argument('--bootstrap',type=int,default=50,help='WGCNA biological-block refits; 0 explicitly skips stability')
    args=parser.parse_args()
    if args.bootstrap < 0:parser.error('--bootstrap must be nonnegative')
    out=Path(args.output);out.mkdir(parents=True,exist_ok=False)
    x=pd.read_csv(args.psi,sep='\t',index_col=0);meta=pd.read_csv(args.samples,sep='\t')
    measured,complete,audit,means,meta=prepare(x,meta)
    audit.to_csv(out/'event_audit.tsv',sep='\t',index_label='event_key')
    means.to_csv(out/'training_means.tsv',sep='\t',index_label='event_key')
    write_json(out/'input_provenance.json',{'files':{str(Path(p).resolve()):hashlib.file_digest(open(p,'rb'),'sha256').hexdigest() for p in [args.psi,args.samples]},
        'version':'0.1.0-research.20261030','input_events':len(x),'eligible_events':len(complete),'feature_cap':None,
        'block_size':2000,'confidence':'within-module coassignment; not legacy composite confidence or hub validation',
        'limitations':'Research CLI; not production UI. Joint/MEGENA may require substantial resources.'})
    status=fit(args.method,complete,meta,out/'model');print(json.dumps(status))
    if status['status']=='completed' and args.method=='wgcna' and args.bootstrap>0:
        bootstrap(measured,meta,out/'model',out/'stability',iterations=args.bootstrap)
    raise SystemExit(0 if status['status']=='completed' else 1)

if __name__=='__main__':main()
