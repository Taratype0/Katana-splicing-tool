"""Direct Iancu co-splicing comparator, same 60 real samples, no feature cap.

The published gene-profile method and event-level production FA4 retain their
native representations. All resampling is of samples, never of distance pairs.
No statistical test treats the dependent distance pairs as replicates.
"""
from pathlib import Path
import argparse,hashlib,json,time,subprocess,traceback,concurrent.futures
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist,squareform

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'evidence/revision30';DATA=OUT/'data';FITS=OUT/'fits'
SOURCE=Path('E:/Project/RNASEQ/Benchmark/public_sources/GSE62669/GSE62669_RNA110714RH_normalizedExonCounts.txt.gz')
METHODS=['CoSplice_FA4','Iancu_CoSplicEx','Iancu_cosine']
R='C:/Program Files/R/R-4.5.0/bin/Rscript.exe'

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def write(p,obj):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(obj,indent=2,allow_nan=False),encoding='utf8')

def prepare():
    DATA.mkdir(parents=True,exist_ok=True)
    if (DATA/'manifest.json').exists():return
    x=pd.read_csv(SOURCE,sep='\t',index_col=0)
    assert x.shape==(114691,60) and x.index.is_unique and x.columns.is_unique
    assert np.isfinite(x.to_numpy()).all() and x.ge(0).all().all()
    gene=x.index.to_series().str.rsplit('_',n=1).str[0];size=gene.value_counts()
    x=x.loc[gene.map(size).ge(2)];gene=gene.loc[x.index]
    total=x.groupby(gene,sort=False).transform('sum')
    usage=x.div(total.where(total.ge(10)))
    np.save(DATA/'usage.npy',usage.to_numpy(float))
    pd.DataFrame({'event_key':usage.index,'gene_id':gene.to_numpy()}).to_csv(DATA/'events.tsv',sep='\t',index=False)
    meta=pd.DataFrame({'sample_id':x.columns,'biological_sample_id':x.columns,'condition':'HSCC_brain','cohort':'GSE62669'})
    meta.to_csv(DATA/'samples.tsv',sep='\t',index=False)
    vectors={m:[] for m in ['canberra','cosine']};genes=[];audit=[]
    for g,part in x.groupby(gene,sort=False):
        v=part.to_numpy().T;genes.append(g)
        for m in vectors:vectors[m].append(pdist(v,metric=m))
        audit.append(dict(gene_id=g,exons=len(part),samples_total_ge10=int((v.sum(axis=1)>=10).sum())))
    pd.DataFrame(audit).to_csv(DATA/'gene_audit.tsv',sep='\t',index=False)
    (DATA/'genes.txt').write_text('\n'.join(genes)+'\n')
    for m,values in vectors.items():np.save(DATA/f'{m}.npy',np.asarray(values))
    rng=np.random.default_rng(30031);draws=[sorted(rng.choice(60,48,replace=False).tolist()) for _ in range(20)]
    assert len({tuple(d) for d in draws})==20
    perm=rng.permutation(60);folds=[]
    for i in range(3):
        test=sorted(perm[20*i:20*i+20].tolist());folds.append(dict(train=[j for j in range(60) if j not in test],test=test))
    write(DATA/'splits.json',dict(seed=30031,deletions=draws,folds=folds,samples=x.columns.tolist()))
    # Numerical checks of distance construction and Mantel equivalence.
    a=np.asarray(vectors['canberra']);ok=np.isfinite(a).all(axis=1)&(np.nanstd(a,axis=1)>1e-12)
    z=a[ok][:4];corr=np.corrcoef(z)
    assert np.allclose(corr,np.array([[np.corrcoef(i,j)[0,1] for j in z] for i in z]))
    ix=draws[0];indices=pair_indices(ix)
    first=x.loc[gene.eq(genes[0])].to_numpy().T
    assert np.allclose(a[0,indices],pdist(first[ix],metric='canberra'))
    write(DATA/'manifest.json',dict(source=str(SOURCE),source_sha256=sha(SOURCE),original_exons=114691,
        original_genes=9066,multi_exon_genes=len(genes),exon_features=len(x),samples=60,
        representation='Normalized exon counts; CoSplice uses exon/gene-sum fractions, not rMATS PSI',
        abundance_mask='Gene total <10 normalized counts, only for fraction input; no integer read-support claim',
        cohort='GSE62669, Iancu original cohort; reused in prior Katana benchmarks, not new independent validation',
        input_scope='All available processed input genes; no additional connectivity-half filter or feature cap',
        phenotype='No verified high/low mapping; no condition tests',
        source_code={p.name:sha(p) for p in (ROOT/'vendor/Iancu_cosplicing_2015').glob('*.R')},
        source_repository='https://github.com/iancuo/cosplicingNetworks',numerical_checks=True))

def pair_indices(cols):
    n=60;a,b=np.triu_indices(len(cols),1);c=np.asarray(cols);i=c[a];j=c[b]
    assert (i<j).all()
    return n*i-i*(i+1)//2+j-i-1

def job(method,name,cols):
    out=FITS/method/name;out.mkdir(parents=True,exist_ok=True)
    if (out/'status.json').exists():
        status=json.loads((out/'status.json').read_text())
        if status['status']=='completed':return status
    start=time.monotonic()
    write(out/'requested.json',dict(method=method,split=name,sample_indices=cols))
    try:
        if method=='CoSplice_FA4':
            from unit_robust import unit_audit
            from algorithm28 import fit_latent
            e=pd.read_csv(DATA/'events.tsv',sep='\t');meta=pd.read_csv(DATA/'samples.tsv',sep='\t')
            raw=pd.DataFrame(np.load(DATA/'usage.npy',mmap_mode='r')[:,cols],index=e.event_key,columns=meta.sample_id.iloc[cols])
            meta=meta.set_index('sample_id').loc[raw.columns].reset_index()
            measured,audit=unit_audit(raw,meta,None)
            audit.to_csv(out/'event_audit.tsv.gz',sep='\t',index_label='event_key')
            measured=measured.loc[audit.eligible]
            # Equivalent to x.T.fillna(x.mean(axis=1)).T without 100k-column
            # pandas assignment overhead. No changes to QC, factors or gates.
            complete=impute(measured)
            status=fit_latent('factor_core',measured,complete,meta,out/'model')
            status.update(eligible_events=len(measured),input_events=len(raw),biological_units=len(meta))
            if status['status']!='completed':raise RuntimeError(str(status))
            a=pd.read_csv(out/'model/assignments.tsv',sep='\t').merge(e,on='event_key',validate='one_to_one')
            a.to_csv(out/'assignments.tsv',sep='\t',index=False)
            # Large diagnostic input tables are losslessly compressed, not retained twice.
            p=out/'model/input.tsv'
            if p.exists():
                import gzip,shutil
                with p.open('rb') as src,gzip.open(str(p)+'.gz','wb',compresslevel=1) as dst:shutil.copyfileobj(src,dst)
                p.unlink()
        else:
            metric='canberra' if method=='Iancu_CoSplicEx' else 'cosine'
            vals=np.load(DATA/f'{metric}.npy',mmap_mode='r')[:,pair_indices(cols)]
            valid=np.isfinite(vals).all(axis=1)&(np.nanstd(vals,axis=1)>1e-12)
            genes=np.asarray((DATA/'genes.txt').read_text().splitlines())
            pd.DataFrame({'gene_id':genes,'estimable':valid}).to_csv(out/'gene_audit.tsv',sep='\t',index=False)
            genes=genes[valid];vals=np.asarray(vals[valid],dtype='<f8',order='C')
            (out/'genes.txt').write_text('\n'.join(genes)+'\n');vals.tofile(out/'distances.bin')
            with (out/'R.log').open('w',encoding='utf8') as log:
                proc=subprocess.run([R,str(ROOT/'science/peer30.R'),str(out),str(len(genes)),str(vals.shape[1]),
                    str(ROOT/'vendor/Iancu_cosplicing_2015/analizeCoSplicEx.R')],stdout=log,stderr=subprocess.STDOUT,timeout=1800)
            if proc.returncode:raise RuntimeError(f'R exit {proc.returncode}; see R.log')
            a=pd.read_csv(out/'assignments.tsv',sep='\t');assert a.gene_id.tolist()==genes.tolist()
            status=dict(status='completed',method=method,eligible_genes=len(a),modules=int(a.loc[a.module.gt(0),'module'].nunique()),
                assigned_genes=int(a.module.gt(0).sum()),samples=len(cols),beta=6,metric=metric,
                execution='Official downstream R clustering expressions with documented input-stage adapter')
            # Deterministic temporary serialization; original distance arrays remain available.
            (out/'distances.bin').unlink()
        status.update(total_seconds=time.monotonic()-start,split=name)
        write(out/'status.json',status)
        print(json.dumps(dict(method=method,split=name,status=status['status'],seconds=round(status['total_seconds'],1),modules=status.get('modules'))),flush=True)
        return status
    except Exception as e:
        status=dict(status='failed',method=method,split=name,error=str(e),seconds=time.monotonic()-start,traceback=traceback.format_exc())
        write(out/'status.json',status);print(json.dumps(status),flush=True);return status

def run(phase,methods):
    prepare();splits=json.loads((DATA/'splits.json').read_text())
    jobs=[]
    if phase in ['full','all']:jobs += [(m,'full',list(range(60))) for m in methods]
    if phase in ['validation','all']:
        for i,c in enumerate(splits['deletions']):jobs += [(m,f'refit_{i:02}',c) for m in methods]
        for i,f in enumerate(splits['folds']):jobs += [(m,f'fold_{i:02}',f['train']) for m in methods]
    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as pool:
        results=list(pool.map(run_tuple,jobs))
    write(OUT/f'jobs_{phase}.json',results)

def run_tuple(args):return job(*args)

def impute(x):
    values=x.to_numpy(copy=True)
    values=np.where(np.isnan(values),x.mean(axis=1).to_numpy()[:,None],values)
    return pd.DataFrame(values,index=x.index,columns=x.columns)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--phase',choices=['prepare','full','validation','all'],default='all');ap.add_argument('--methods',nargs='+',default=METHODS);args=ap.parse_args()
    if args.phase=='prepare':prepare()
    else:run(args.phase,args.methods)
