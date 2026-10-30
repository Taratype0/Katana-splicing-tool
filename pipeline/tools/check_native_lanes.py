#!/usr/bin/env python3
"""Test installed native lane grouping/merging with generated reads, including resume."""
import argparse
import csv
import gzip
import json
import os
from pathlib import Path
import subprocess

p=argparse.ArgumentParser()
group=p.add_mutually_exclusive_group(required=True)
group.add_argument('--nextflow')
group.add_argument('--nextflow-command-file',type=Path)
p.add_argument('--pipeline', type=Path, required=True)
p.add_argument('--outdir', type=Path, required=True)
a=p.parse_args()
root=a.outdir.resolve(); root.mkdir(parents=True, exist_ok=False)
rows=[]
for i in range(6):
    for lane in ['L2','L1']:
        row=dict(sample=f'S{i}',condition='WT' if i%2==0 else 'MUT',subject=f'mouse{i//2}',lane=lane,strandedness='auto',read_length='4')
        for mate in (1,2):
            f=root/f'S{i}_{lane}_{mate}.fq.gz'
            with gzip.open(f,'wb') as handle: handle.write(f'@S{i}_{lane}/{mate}\nACGT\n+\nIIII\n'.encode())
            row[f'fastq_{mate}']=str(f)
        rows.append(row)
with (root/'samples.csv').open('w',newline='') as handle:
    writer=csv.DictWriter(handle,list(rows[0])); writer.writeheader(); writer.writerows(rows)
(root/'contrasts.csv').write_text('contrast,treatment,control,paired\nMUT_vs_WT,MUT,WT,true\n')
(root/'ref.fa').write_text('>chr1\nACGT\n'); (root/'ref.gtf').write_text('# generated lane transport fixture\n')
(root/'test.config').write_text("process { executor='local'; cpus=1; memory='1 GB'; time='5m' }\n")
prefix=json.loads(a.nextflow_command_file.read_text()) if a.nextflow_command_file else [a.nextflow]
base=[*prefix,'run',str(a.pipeline.resolve()/'tests/lanes.nf'),'-c',str(root/'test.config'),'-ansi-log','false',
      '--input',str(root/'samples.csv'),'--contrasts',str(root/'contrasts.csv'),'--fasta',str(root/'ref.fa'),'--gtf',str(root/'ref.gtf'),'-work-dir',str(root/'work')]
results={}
for name,flags in [('first',[]),('resume',['-resume'])]:
    with (root/f'{name}.log').open('w') as stream:
        subprocess.run(base+flags+['-with-trace',str(root/f'{name}.tsv')],cwd=root,env={**os.environ,'NXF_OFFLINE':'true'},stdout=stream,stderr=subprocess.STDOUT,check=True,timeout=300)
    with (root/f'{name}.tsv').open() as handle: trace=list(csv.DictReader(handle,delimiter='\t'))
    expected='CACHED' if flags else 'COMPLETED'
    assert len(trace)==7 and all(r['status']==expected for r in trace),trace
    results[name]=dict(tasks=7,status=expected)
audits=list((root/'work').glob('*/*/fastq_staging_audit.json'))
assert len(audits)==6
for audit_file in audits:
    audit=json.loads(audit_file.read_text()); assert audit['lane_count']==2
    for mate in (1,2):
        with gzip.open(audit_file.parent/f"{audit['sample']}_R{mate}.fastq.gz",'rt') as handle:
            assert handle.read().splitlines()[::4]==[f"@{audit['sample']}_L1/{mate}",f"@{audit['sample']}_L2/{mate}"]
(root/'acceptance.json').write_text(json.dumps(dict(status='passed',samples=6,biological_pairs=3,native=results,paired_lane_order_verified=True),indent=2)+'\n')
print((root/'acceptance.json').read_text())
