"""Data-selected shared-program rank followed by the native count-aware factor fit."""
import json,shutil
from pathlib import Path
from common32 import ROOT,EVID,sha
from parallel_rank33 import select
from run_native32 import fit as native_fit

def fit(method,input,out,timeout=21600):
 if method!='ebmf_counts':raise ValueError('This runner is the adaptive count-factor engine')
 input,out=map(Path,[input,out]);rankdir=out.parent/(out.name+'_rank')
 cached=rankdir/'parallel_rank.json';valid=False
 if cached.exists():
  r=json.loads(cached.read_text());valid=(r['input_sha256']==sha(input/'input.bin') and r['gene_mapping_sha256']==sha(input/'events.tsv')
    and r['samples_sha256']==sha(input/'samples.tsv') and r.get('selector_sha256')==sha(ROOT/'science/parallel_rank33.py') and r['repetitions']==999)
 if not valid:select(input,rankdir)
 prepared={n:rankdir/n for n in ['parallel_rank.tsv','parallel_rank.json','parallel_null_eigenvalues.npy']}
 return native_fit(method,input,out,timeout=timeout,fast_variance=True,engine_source='native33.R',helper_sources=['parallel_rank33.py'],prepared_files=prepared)
