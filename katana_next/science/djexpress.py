"""Source-pinned DJExpress integration. No Katana association-score substitution."""
from pathlib import Path
import hashlib,json,subprocess,gzip,shutil,re,os
import numpy as np,pandas as pd

def sources():
    root=Path(__file__).resolve().parents[1]/'vendor/DJExpress'
    expected=json.loads((root.parent/'DJExpress_SOURCE_SHA256.json').read_text(encoding='utf8'))
    for name,sha in expected['files'].items():
        if hashlib.sha256((root/name).read_bytes()).hexdigest()!=sha:
            raise ValueError('DJExpress author source integrity failed: '+name)
    return root,expected

def validate_inputs(counts,traits,gene_expression=None):
    if not counts.index.is_unique or not counts.columns.is_unique or counts.shape[1]<15 or counts.empty:
        raise ValueError('Unique feature/sample IDs and >=15 biological samples are required')
    y=counts.to_numpy(float)
    if not np.isfinite(y).all() or (y<0).any() or not np.equal(y,np.floor(y)).all():
        raise ValueError('DJExpress requires measured integer count evidence; PSI/TPM/DTU/DEXSeq statistics must not be converted into counts')
    if any(not re.fullmatch(r'[^:]+:\d+:\d+:[12]',str(f)) for f in counts.index):
        raise ValueError('Require junction-coordinate IDs chr:start:end:strand (1 or 2); label rMATS form evidence explicitly')
    if not traits.index.is_unique or not traits.columns.is_unique or set(traits.index)!=set(counts.columns):
        raise ValueError('Trait sample IDs must exactly match counts')
    t=traits.loc[counts.columns].to_numpy(float)
    if traits.shape[1]<2:
        raise ValueError('Pinned author JCNA1pass assumes >=2 traits; provide two real observed traits or >=2 observed condition indicators')
    if not np.isfinite(t).all() or np.any(np.var(t,axis=0)<=0):
        raise ValueError('Traits must be finite, observed and nonconstant; no fabricated dummy traits')
    if (y.sum(0)==0).any():raise ValueError('Zero-support samples cannot be treated as observed libraries')
    if gene_expression is not None:
        g=gene_expression
        if not g.index.is_unique or not g.columns.is_unique or set(g.columns)!=set(counts.columns):raise ValueError('Gene expression IDs must be unique and samples must match exactly')
        if not np.isfinite(g.to_numpy(float)).all():raise ValueError('Measured normalized gene expression must be finite')

def run_djexpress(counts,traits,out,*,rscript=None,gene_expression=None,gtf=None,input_kind='junction_counts'):
    if input_kind not in ['junction_counts','rmats_form_evidence']:raise ValueError('Unknown count evidence identity')
    validate_inputs(counts,traits,gene_expression)
    if (gene_expression is None)!=(gtf is None):raise ValueError('Second pass requires both measured normalized gene expression and the matching GTF')
    root,author=sources();out=Path(out);out.mkdir(parents=True,exist_ok=False);prepared=out/'prepared';prepared.mkdir()
    qc=pd.DataFrame(index=counts.index);qc['all_zero']=counts.sum(1).eq(0);qc['constant_count']=counts.var(1).eq(0)
    qc['passed']=~qc.all_zero;qc['reason']=np.where(qc.all_zero,'all_zero','passed')
    qc.rename_axis('event_key').to_csv(out/'EVERY_FEATURE_PREPARATION.tsv',sep='\t')
    retained=counts.loc[qc.passed];
    if len(retained)<4:raise ValueError('Fewer than four nonconstant count features remain; preparation audit retained')
    retained.to_csv(prepared/'counts.tsv',sep='\t');traits.loc[counts.columns].to_csv(prepared/'traits.tsv',sep='\t')
    if gene_expression is not None:
        gene_expression.loc[:,counts.columns].to_csv(prepared/'gene_expression.tsv',sep='\t')
        gtf=Path(gtf)
        if gtf.suffix=='.gz':shutil.copyfile(gtf,prepared/'reference.gtf.gz')
        else:
            with gtf.open('rb') as a,gzip.open(prepared/'reference.gtf.gz','wb') as b:shutil.copyfileobj(a,b)
    if not rscript:
        from katana_cosplice.benchmark import _rscript_with_wgcna
        rscript,diagnostic=_rscript_with_wgcna()
        if not rscript:raise RuntimeError(diagnostic)
    provenance=dict(engine='DJExpress',revision=author['revision'],author_files_verified=len(author['files']),source_unchanged=True,input_kind=input_kind,
        preparation='Katana prepared-object adapter: measured counts -> limma voom -> original JCNA1pass; optional original genePrepare/2pass',
        native_raw_entry_used=False,method_scope='author core integration, not an end-to-end JCNAprepare reproduction',primary_stage='JCNA1pass',
        second_pass='original per-trait outputs retained separately',novel_algorithm=False,
        prepared_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in prepared.iterdir()})
    (out/'ENGINE_PROVENANCE.json').write_text(json.dumps(provenance,indent=2),encoding='utf8')
    command=[str(rscript),'--vanilla',str(Path(__file__).with_name('djexpress_driver.R')),str(root),str(prepared.resolve()),str(out.resolve())]
    from katana_cosplice.benchmark import _r_subprocess_env
    environment=_r_subprocess_env(str(rscript))
    isolated_library=Path(__file__).resolve().parents[3]/'r_library'
    if isolated_library.is_dir():
        environment['R_LIBS']=str(isolated_library)+os.pathsep+environment.get('R_LIBS','')
    with (out/'R_STDOUT.log').open('w',encoding='utf8') as a,(out/'R_STDERR.log').open('w',encoding='utf8') as b:
        result=subprocess.run(command,stdout=a,stderr=b,env=environment)
    if result.returncode or not (out/'ENGINE_COMPLETE.json').is_file():
        (out/'RUN_FAILURE.json').write_text(json.dumps(dict(status='failed',R_exit_code=result.returncode,engine='DJExpress',fallback_run=False,logs=['R_STDOUT.log','R_STDERR.log']),indent=2),encoding='utf8')
        raise RuntimeError(f'DJExpress failed; original error/logs retained at {out}. No fallback method was run.')
    return out
