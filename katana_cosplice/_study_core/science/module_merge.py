"""Native WGCNA eigengene merging of a fitted sparse partition.

This is a borrowed standard hierarchical refinement, not a novel algorithm.
The external secondary comparison was added after poor initial deletion
stability was seen; it must not replace the prespecified primary benchmark.
"""
from pathlib import Path
import subprocess,time,json,shutil,os
from cosplice_research import ROOT,RSCRIPT,write_json

def run_owned(args,log,timeout=300):
    """Rscript may launch an R child on Windows; retain ownership on timeout."""
    proc=subprocess.Popen(args,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
    try:
        return proc.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        import psutil
        try:children=psutil.Process(proc.pid).children(recursive=True)
        except psutil.NoSuchProcess:children=[]
        for child in reversed(children):
            try:child.kill()
            except psutil.NoSuchProcess:pass
        if proc.poll() is None:proc.kill()
        proc.wait(timeout=10)
        _,alive=psutil.wait_procs(children,timeout=10)
        if alive:raise TimeoutError(f'Native merge exceeded {timeout}s; owned child cleanup incomplete: {[p.pid for p in alive]}')
        raise TimeoutError(f'Native merge exceeded {timeout}s; owned process tree stopped') from None


def merge(model,out):
    model=Path(model);out=Path(out)
    if (out/'status.json').exists():return json.loads((out/'status.json').read_text())
    out.mkdir(parents=True,exist_ok=True);started=time.monotonic()
    shutil.copyfile(ROOT/'science/merge_modules.R',out/'engine_source.R')
    with (out/'engine.log').open('w',encoding='utf-8') as log:
        try:
            code=run_owned([str(RSCRIPT),'--vanilla',str(out/'engine_source.R'),str(ROOT),str(model/'input.tsv'),str(model/'assignments.tsv'),str(out)],log,timeout=300)
            if code:raise RuntimeError(f'Native merge exit {code}')
            result={'status':'completed','method':'jackknife_wgcna_merge','cut_height':.25,'iterate':True,'native_function':'WGCNA::mergeCloseModules','base_model':str(model),'scope':'Secondary exploratory hierarchical refinement; original network and primary comparison retained'}
        except Exception as e:result={'status':'failed','method':'jackknife_wgcna_merge','error':str(e)}
    result['seconds']=time.monotonic()-started;write_json(out/'status.json',result);return result

if __name__=='__main__':
    import sys
    print(merge(sys.argv[1],sys.argv[2]))
