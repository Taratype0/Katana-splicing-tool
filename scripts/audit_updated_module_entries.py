"""Read-only entry/packaging audit; never submits or starts an analysis."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import sys


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--release',type=Path,required=True)
    p.add_argument('--pipeline-repo',type=Path,required=True)
    p.add_argument('--outdir',type=Path,required=True)
    a=p.parse_args(); a.outdir.mkdir(parents=True,exist_ok=True)
    sha=lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    seal=json.loads((a.release/'SOURCE_SHA256.json').read_text('utf-8-sig'))
    mismatches=[name for name,h in seal.items() if not (a.release/name).is_file() or sha(a.release/name)!=h]
    if mismatches: raise ValueError('Selected sealed release changed: '+str(mismatches))
    probe="""import json,inspect
from katana_next.workflows.registry import WorkflowRegistry
from katana_motiftrace import scientific_release,__version__
import katana_cosplice.cli
r=WorkflowRegistry()
print(json.dumps({'motiftrace_version':__version__,'scientific_release':scientific_release,
 'cosplice_cli_source':inspect.getfile(katana_cosplice.cli),
 'workflows':{k:{'module':type(r.get(k)).__module__,'version':r.get(k).version,
 'defaults':r.get(k).normalize_parameters({})} for k in ['programs.cosplice','programs.psi_network','prioritize.motiftrace_consequences','reports.study']}},indent=2))
"""
    results={}
    for label,root in [('updated_local',a.release/'app'),('pipeline_branch',a.pipeline_repo)]:
        env={**os.environ,'PYTHONPATH':str(root),'PYTHONDONTWRITEBYTECODE':'1'}
        run=subprocess.run([sys.executable,'-B','-c',probe],cwd=a.outdir,env=env,text=True,capture_output=True)
        (a.outdir/(label+'.stdout.json')).write_text(run.stdout,encoding='utf8')
        (a.outdir/(label+'.stderr.log')).write_text(run.stderr,encoding='utf8')
        if run.returncode: raise RuntimeError(run.stderr)
        results[label]=json.loads(run.stdout)
    files={
      'planner':a.pipeline_repo/'katana_next/flow/planner.py',
      'network_adapter':a.pipeline_repo/'katana_next/flow/stages.py',
      'motif_workflow':a.pipeline_repo/'katana_next/workflows/motiftrace_consequences.py',
      'nf_downstream':a.pipeline_repo/'pipeline/workflows/katana_downstream.nf',
      'package':a.pipeline_repo/'pyproject.toml',
      'new_network_workflow':a.release/'app/katana_next/workflows/cosplice_study.py',
      'new_network_core':a.release/'cosplice_core/science/cosplice_research.py',
      'new_product':a.release/'cosplice_core/science/product.py',
      'new_report_backend':a.release/'FORMAL_REPORT_BACKEND.json',
      'new_standalone_network_cli':a.release/'app/katana_cosplice/cli.py'}
    tokens={
      'planner':['network_scope =','flow.network.fit','motiftrace_consequences'],
      'network_adapter':['from katana_cosplice','def _network('],
      'motif_workflow':['from katana_next.science.motiftrace_execution'],
      'nf_downstream':['NETWORK_FIT(', 'NETWORK_FOLLOWUP(', 'MOTIFTRACE(', 'KATANA_FLOW_REPORT('],
      'package':['katana-cosplice =','include =','katana_motiftrace'],
      'new_network_workflow':['cosplice_core','programs.cosplice','version=',"max_events",'production43.py'],
      'new_network_core':['RSCRIPT=','minimum_units=','8*2**30','timeout_seconds'],
      'new_product':['from cosplice_research','prepare(x'],
      'new_report_backend':['python_executable','release_relative_path'],
      'new_standalone_network_cli':['CoSplicePipeline','run_stage','config = load_config']}
    refs={key:{'path':str(path),'sha256':sha(path),'lines':[{'line':i,'text':line.strip()}
        for i,line in enumerate(path.read_text('utf-8-sig').splitlines(),1) if any(t in line for t in tokens[key])]} for key,path in files.items()}
    # Exercise exactly the updated preflight with six generated independent units.
    preflight="""import json,numpy as np,pandas as pd
from cosplice_research import prepare
x=pd.DataFrame(np.random.default_rng(433).uniform(.1,.9,(50,6)),columns=[f's{i}' for i in range(6)])
m=pd.DataFrame({'sample_id':x.columns,'biological_sample_id':x.columns,'condition':['A']*3+['B']*3})
try: prepare(x,m); print(json.dumps({'six_units_accepted':True}))
except ValueError as e: print(json.dumps({'six_units_accepted':False,'error':str(e)}))
"""
    env={**os.environ,'PYTHONPATH':str(a.release/'cosplice_core/science'),'PYTHONDONTWRITEBYTECODE':'1'}
    run=subprocess.run([sys.executable,'-B','-c',preflight],cwd=a.outdir,env=env,text=True,capture_output=True)
    if run.returncode: raise RuntimeError(run.stderr)
    results['six_unit_native_preflight']=json.loads(run.stdout)
    result={'schema':'katana.module_entry_audit/v1','status':'integration_gaps_confirmed',
        'release':str(a.release),'sealed_files_verified':len(seal),'source_mismatches':mismatches,
        'observed':results,'source_evidence':refs,'science_executed':False,
        'main_changed':False,'blocking_gaps':[
          'Pipeline branch MotifTrace workflow still imports the legacy execution adapter, not the selected EventGraph dispatcher.',
          'Nextflow network tasks call flow.network.fit / legacy CoSplicePipeline rather than programs.cosplice 4.3.3.',
          'katana-cosplice standalone CLI still invokes the legacy CoSplicePipeline; the new workflow is only registered through katana-next/browser.',
          'cosplice_core is a sibling of app in the local sealed release and is absent from pipeline wheel package discovery.',
          'New network native R command is a fixed Windows path; runtime memory/time limits and 15-unit check remain in that implementation.',
          'Local formal report workflow has a pinned external Python/renderer release configuration, not a portable installed backend.',
          'Updated workflow, packaged native dependencies, report collection and matching Linux container require joint acceptance before main promotion.']}
    (a.outdir/'AUDIT.json').write_text(json.dumps(result,indent=2),encoding='utf8')
    print(json.dumps({'status':result['status'],'sealed_files_verified':len(seal),'preflight':results['six_unit_native_preflight']}))


if __name__=='__main__': main()
