"""Bounded owned process execution; no native model substitution."""
from pathlib import Path
import argparse,subprocess,os,time,json,shutil
import psutil
from common32 import ROOT,EVID,write,sha
def fit(method,input,out,timeout=21600,r_home=None,fast_variance=False,engine_source='native32.R',helper_sources=None,prepared_files=None):
 helper_sources=['rank32.R','export_tom32.R'] if helper_sources is None else helper_sources
 prepared_files={} if prepared_files is None else prepared_files
 input=Path(input);out=Path(out);out.mkdir(parents=True,exist_ok=True)
 if out.resolve()==ROOT.resolve() or not out.resolve().is_relative_to(ROOT.resolve()) or out.resolve()==input.resolve():raise ValueError('Native output must be an owned analysis subdirectory')
 if (out/'status.json').exists():
  s=json.loads((out/'status.json').read_text())
  if s['status']=='completed':
   if s.get('method')!=method or s.get('input_manifest_sha256')!=sha(input/'manifest.json') or s.get('engine_sha256')!=sha(ROOT/'science'/engine_source):
    raise RuntimeError('Completed fit cache does not match the current method, input or engine')
   if r_home and s.get('runtime',{}).get('rblas_sha256')!=sha(Path(r_home)/'bin/x64/Rblas.dll'):raise RuntimeError('Completed fit cache does not match the explicitly requested runtime')
   for name in helper_sources:
    if s.get('helper_sha256',{}).get(name)!=sha(ROOT/'science'/name):raise RuntimeError('Completed fit cache helper changed: '+name)
   for name,p in prepared_files.items():
    if s.get('prepared_sha256',{}).get(name)!=sha(p) or sha(out/name)!=sha(p):raise RuntimeError('Prepared selector evidence changed: '+name)
   if s.get('implementation_adapter'):
    for name,h in s['implementation_adapter']['files'].items():
     if sha(out/name)!=h:raise RuntimeError('Cached arithmetic adapter changed: '+name)
   if not (out/'COMPLETED').is_file():raise RuntimeError('Completed status is missing the native completion marker')
   required=['native_model.rds']+(['loadings.tsv','factor_scores.tsv','factors.tsv'] if method.startswith('ebmf') else ['assignments.tsv'])
   if any(not (out/n).is_file() or not (out/n).stat().st_size for n in required):raise RuntimeError('Completed cache is missing native model artifacts')
   return s
 if (out/'worker.json').exists():
  old=json.loads((out/'worker.json').read_text())
  try:
   p=psutil.Process(old['pid'])
   if abs(p.create_time()-old.get('create_time',0))<1:raise RuntimeError('Output directory already has a live owned worker')
  except psutil.NoSuchProcess:pass
 # Guard older Windows launcher orphans before reusing a partially written
 # output. Match this exact output argument, never terminate unrelated jobs.
 for proc in psutil.process_iter(['name','cmdline']):
  try:
   cmd=proc.info['cmdline'] or []
   if proc.info['name'].lower()=='rscript.exe' and any(Path(arg).name in ['native_engine.R','native_accelerated32.R'] for arg in cmd) and cmd and Path(cmd[-1]).resolve()==out.resolve():raise RuntimeError('Output still has a live native R process')
  except (psutil.NoSuchProcess,psutil.AccessDenied):continue
 old_files=[p for p in out.iterdir() if p.name!='attempts']
 if old_files:
  archive=out/'attempts'/time.strftime('%Y%m%d-%H%M%S');archive.mkdir(parents=True,exist_ok=False)
  for p in old_files:
   if p.is_symlink() or not p.resolve().is_relative_to(out.resolve()):raise RuntimeError('Refusing to move a linked native artifact')
   p.rename(archive/p.name)
  write(archive/'retry_preservation.json',dict(reason='Preserved incomplete previous native attempt before retry',previous_status='status.json' if (archive/'status.json').exists() else 'missing'))
 start=time.monotonic();result=dict(method=method,input_manifest_sha256=sha(input/'manifest.json'),engine_sha256=sha(ROOT/'science'/engine_source),engine_source=engine_source)
 engine=out/'native_engine.R';shutil.copyfile(ROOT/'science'/engine_source,engine)
 for name in helper_sources:shutil.copyfile(ROOT/'science'/name,out/name)
 result['helper_sha256']={name:sha(out/name) for name in helper_sources}
 for name,p in prepared_files.items():
  if Path(name).name!=name:raise ValueError('Prepared filenames must be simple names')
  shutil.copyfile(p,out/name)
 result['prepared_sha256']={name:sha(out/name) for name in prepared_files}
 runtime=Path(r_home or 'C:/Program Files/R/R-4.5.0')
 env=os.environ.copy()
 if r_home:env.update(R_HOME=str(runtime),OPENBLAS_NUM_THREADS='1')
 result['runtime']={'home':str(runtime),'r_dll_sha256':sha(runtime/'bin/x64/R.dll'),'rblas_sha256':sha(runtime/'bin/x64/Rblas.dll')}
 launch=engine
 if fast_variance:
  if method!='ebmf_counts' or r_home:raise ValueError('Verified arithmetic adapter supports only count-aware EBMF with the reference R runtime')
  gate=json.loads((EVID/'acceptance/accelerated/acceptance.json').read_text())
  if not gate['passed'] or gate['adapter_source_sha256']!=sha(ROOT/'science/variance32.R') or gate['adapter_dll_sha256']!=sha(ROOT/'runtime/variance32/variance32.dll'):raise RuntimeError('Accelerated arithmetic acceptance is missing or stale')
  for n in ['variance32.R','native_accelerated32.R','variance32.c']:shutil.copyfile(ROOT/'science'/n,out/n)
  shutil.copyfile(ROOT/'runtime/variance32/variance32.dll',out/'variance32.dll')
  result['implementation_adapter']=dict(kind='same native Brent optimizer and scalar likelihood; arithmetic acceleration only',files={n:sha(out/n) for n in ['variance32.R','native_accelerated32.R','variance32.c','variance32.dll']},acceptance_sha256=sha(EVID/'acceptance/accelerated/acceptance.json'),mkl_rt_sha256=sha('D:/anaconda3/Library/bin/mkl_rt.2.dll'))
  launch=out/'native_accelerated32.R';env.update(PATH='D:/anaconda3/Library/bin;'+env['PATH'],MKL_THREADING_LAYER='SEQUENTIAL',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
 with (out/'engine.log').open('w') as log:
  proc=subprocess.Popen([str(runtime/'bin/x64/Rscript.exe'),'--vanilla',str(launch),str(ROOT),method,str(input),str(out)],env=env,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
  write(out/'worker.json',dict(pid=proc.pid,create_time=psutil.Process(proc.pid).create_time()));peak=0;seen={}
  while proc.poll() is None:
   try:
    owned=psutil.Process(proc.pid);tree=[owned]+owned.children(recursive=True)
    for p in tree:seen[(p.pid,p.create_time())]=p
    peak=max(peak,sum(p.memory_info().rss for p in tree if p.is_running()))
   except psutil.NoSuchProcess:tree=[]
   if time.monotonic()-start>timeout or peak>12*2**30:
    for p in reversed(tree):
     try:p.terminate()
     except psutil.NoSuchProcess:pass
    proc.wait(timeout=20);result.update(status='failed',error=f'Owned native fit exceeded {timeout} seconds or 12 GiB');break
   time.sleep(1)
  else:result.update(status='completed' if proc.returncode==0 and (out/'COMPLETED').exists() else 'failed',exit_code=proc.returncode)
  # Windows Rscript can leave its native child alive when the launcher exits.
  # Only terminate processes observed as descendants of this exact invocation.
  for p in reversed(list(seen.values())):
   try:
    if p.is_running():p.terminate()
   except psutil.NoSuchProcess:pass
  psutil.wait_procs(list(seen.values()),timeout=10)
 result.update(seconds=time.monotonic()-start,peak_rss_bytes=peak);write(out/'status.json',result);print(json.dumps(result),flush=True);return result
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('method',choices=['ebmf','ebmf_adaptive','ebmf_counts','wgcna']);ap.add_argument('split',default='full',nargs='?');a=ap.parse_args()
 fit(a.method,EVID/'inputs'/a.split,EVID/'development'/a.method/a.split)
