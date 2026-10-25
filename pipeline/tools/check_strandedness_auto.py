#!/usr/bin/env python3
"""Run native auto-strand modules on generated alignments, never study data."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
from pathlib import Path


GENERATOR = r'''
import gzip,json,sys
from pathlib import Path
import pysam
root=Path(sys.argv[1])
cases=json.loads((root/'cases.json').read_text())
for case in cases:
    directory=root/case['id'];directory.mkdir()
    for i in range(4):
        sample='S'+str(i);single=case['single']
        orientation=case['orientation']
        if orientation=='mixed':orientation='forward' if i==0 else 'reverse'
        records=[]
        for n in range(600):
            gene_minus=n%2==1
            forward=orientation=='forward' or (orientation=='unstranded' and n%4<2) or (orientation=='ambiguous' and n%20<13)
            reverse_read=(not forward)!=gene_minus
            start=(3100 if gene_minus else 200)+(n%80)
            name=f'{sample}_{n}'
            for mate in range(1 if single else 2):
                read=pysam.AlignedSegment()
                read.query_name=name;read.query_sequence='A'*100
                reverse=reverse_read if mate==0 else not reverse_read
                read.flag=(16 if reverse else 0)
                if not single:
                    read.flag|=1|2|(64 if mate==0 else 128)|(32 if not reverse else 0)
                    read.next_reference_id=0
                read.reference_id=0
                read.reference_start=start+(200 if reverse else 0)
                read.next_reference_start=start+(0 if reverse else 200) if not single else -1
                read.template_length=(-300 if reverse else 300) if not single else 0
                read.mapping_quality=60;read.cigarstring='100M'
                read.query_qualities=pysam.qualitystring_to_array('I'*100)
                records.append(read)
        records.sort(key=lambda r:r.reference_start)
        with pysam.AlignmentFile(str(directory/(sample+'.bam')),'wb',header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':6000}]}) as f:
            for r in records:f.write(r)
        pysam.index(str(directory/(sample+'.bam')))
    (directory/'annotation.gtf').write_text('chr1\tfixture\texon\t101\t1900\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\nchr1\tfixture\texon\t3001\t4800\t.\t-\t.\tgene_id "g2"; transcript_id "t2";\n'.replace('\\t','\t').replace('\\n','\n'))
'''


def run(command, cwd, log, *, timeout=600):
    result = subprocess.run(command, cwd=cwd, text=True, capture_output=True, timeout=timeout)
    log.with_suffix('.stdout.txt').write_text(result.stdout, encoding='utf-8')
    log.with_suffix('.stderr.txt').write_text(result.stderr, encoding='utf-8')
    return result


def trace(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir',type=Path,required=True)
    parser.add_argument('--nextflow-command-json',default='["nextflow"]')
    parser.add_argument('--docker',default='docker')
    parser.add_argument('--rseqc-image',default='quay.io/biocontainers/rseqc:5.0.3--py39hf95cd2a_0')
    args=parser.parse_args()
    root=args.outdir.resolve();root.mkdir(parents=True,exist_ok=False)
    pipeline=Path(__file__).resolve().parents[1]
    runtime=json.loads(args.nextflow_command_json)
    assert isinstance(runtime,list) and runtime and all(isinstance(x,str) for x in runtime)
    cases=[{'id':f'{layout}-{direction}','single':layout=='se','orientation':direction,'expected':direction,'success':True}
           for layout in ('pe','se') for direction in ('forward','reverse','unstranded')]
    cases += [dict(id='explicit-only',single=False,orientation='forward',expected='forward',success=True),
              dict(id='ambiguous',single=False,orientation='ambiguous',expected='',success=False),
              dict(id='mixed',single=False,orientation='mixed',expected='',success=False)]
    (root/'cases.json').write_text(json.dumps(cases))
    (root/'generate.py').write_text(GENERATOR)
    generated=run([args.docker,'run','--rm','--network','none','--user',f'{os.getuid()}:{os.getgid()}',
                   '-v',f'{root}:/fixture',args.rseqc_image,'python','/fixture/generate.py','/fixture'],root,root/'generate')
    assert generated.returncode==0,(generated.stdout,generated.stderr)
    module=(pipeline/'subworkflows/local/resolve_strandedness').as_posix().replace("'","\\'")
    (root/'main.nf').write_text("""nextflow.enable.dsl=2
include { RESOLVE_STRANDEDNESS } from 'MODULE'
process CONSUMER {
    tag "${meta.id}"
    input: tuple val(meta), path(bam)
    output: path "${meta.id}.resolved.json"
    script:
    if (meta.strandedness != params.expected) error("Wrong propagated strand for ${meta.id}: ${meta.strandedness}")
    \"\"\"
    printf '%s\\n' '${groovy.json.JsonOutput.toJson(meta)}' > ${meta.id}.resolved.json
    \"\"\"
}
workflow {
    rows=Channel.fromPath(params.samples).splitCsv(header:true)
    bams=rows.map { row -> tuple([id:row.sample,source:row.source,condition:row.condition,
        strandedness:row.strandedness,single_end:row.single_end.toBoolean()],file("${params.inputs}/${row.sample}.bam")) }
    indices=rows.map { row -> tuple([id:row.sample],file("${params.inputs}/${row.sample}.bam.bai")) }
    RESOLVE_STRANDEDNESS(bams,indices,Channel.value(file(params.samples)),Channel.value(file(params.contrasts)),
        Channel.value(file("${params.inputs}/annotation.gtf")),Channel.value(file(params.resolver)),Channel.value(file(params.validator)))
    CONSUMER(RESOLVE_STRANDEDNESS.out.bam)
}
""".replace('MODULE',module))
    # Explicitly small resource requests are for this tiny generated fixture.
    # Production scheduling/resource parameters are unchanged by this check.
    fixture_config = "process.cpus=1\nprocess.memory='1 GB'\nprocess.time='5m'\nprocess.errorStrategy='terminate'\n" + \
        "docker.runOptions='--network none --user $(id -u):$(id -g)'\n" + \
        "process { withName: '(GFFREAD_STRAND_MODELS|RSEQC_INFEREXPERIMENT|RESOLVE_SAMPLE_STRANDEDNESS)' { cpus=1; memory='1 GB' } }\n"
    if args.rseqc_image != 'quay.io/biocontainers/rseqc:5.0.3--py39hf95cd2a_0':
        # Explicit developer diagnostic override, recorded in the receipt. CI
        # leaves this unset and qualifies the production RSeQC image itself.
        fixture_config += "process { withName: 'RSEQC_INFEREXPERIMENT' { container=" + json.dumps(args.rseqc_image) + " } }\n"
    (root/'fixture.config').write_text(fixture_config)
    result_cases=[]
    for case in cases:
        inputs=root/case['id'];case_root=root/('run-'+case['id']);case_root.mkdir()
        rows=[]
        for i in range(4):
            requested=case['expected'] if case['id']=='explicit-only' or (case['id']=='pe-forward' and i%2) else 'auto'
            rows.append(dict(sample='S'+str(i),condition='CONTROL' if i<2 else 'TREATED',
                source='fastq' if i%2==0 else 'genome_bam',assay='bulk_rna_seq',genome_bam=str(inputs/f'S{i}.bam'),
                fastq_1='',fastq_2='',transcriptome_bam='',strandedness=requested,
                single_end=str(case['single']).lower(),read_length='100',subject=f'animal{i}',batch=''))
        with (inputs/'samples.csv').open('w',newline='') as f:
            writer=csv.DictWriter(f,list(rows[0]));writer.writeheader();writer.writerows(rows)
        (inputs/'contrasts.csv').write_text('contrast,treatment,control\nTREATED_vs_CONTROL,TREATED,CONTROL\n')
        trace_file=case_root/'trace.tsv'
        command=runtime+['run',str(root/'main.nf'),'-c',str(pipeline/'nextflow.config'),'-c',str(root/'fixture.config'),
            '-profile','docker','--outdir',str(case_root/'outputs'),'--samples',str(inputs/'samples.csv'),
            '--contrasts',str(inputs/'contrasts.csv'),'--inputs',str(inputs),'--expected',case['expected'] or 'never',
            '--resolver',str(pipeline/'bin/katana_resolve_strandedness.py'),'--validator',str(pipeline/'bin/katana_validate_inputs.py'),
            '-ansi-log','false','-work-dir',str(case_root/'work'),'-with-trace',str(trace_file)]
        result=run(command,case_root,case_root/'first')
        assert (result.returncode==0)==case['success'],(case['id'],result.stdout,result.stderr)
        task_rows=trace(trace_file)
        calls=[r for r in task_rows if 'RSEQC_INFEREXPERIMENT' in r['name']]
        consumers=[r for r in task_rows if r['name'].startswith('CONSUMER')]
        if case['success']:
            assert len(consumers)==4 and all(r['status']=='COMPLETED' for r in task_rows),case['id']
            audits=list((case_root/'work').glob('*/*/strandedness/strandedness.audit.json'))
            assert len(audits)==1,case['id']
            audit=json.loads(audits[0].read_text())
            assert all(r['strandedness']==case['expected'] for r in audit['samples'])
            assert len(calls)==audit['auto_samples']
            assert all(r.get('usable_reads_sampled',1)>0 for r in audit['samples'] if r['requested']=='auto')
            if case['id']=='explicit-only':
                assert not any('GFFREAD_STRAND_MODELS' in r['name'] for r in task_rows)
        else:
            assert not consumers,'A failed strand gate released a downstream task'
            assert any('RESOLVE_SAMPLE_STRANDEDNESS' in r['name'] and r['status']=='FAILED' for r in task_rows)
        entry={'case':case['id'],'exit_code':result.returncode,'native_inference_tasks':len(calls),
               'downstream_consumers':len(consumers),'expected_success':case['success']}
        if case['id']=='pe-forward':
            resumed=run(command[:-2]+['-with-trace',str(case_root/'trace-resume.tsv'),'-resume'],case_root,case_root/'resume')
            assert resumed.returncode==0,(resumed.stdout,resumed.stderr)
            cached=trace(case_root/'trace-resume.tsv')
            assert len(cached)==len(task_rows) and all(r['status']=='CACHED' for r in cached)
            entry['actual_resume_cached_tasks']=len(cached)
        result_cases.append(entry)
        print(json.dumps(entry),flush=True)
    receipt={'schema':'katana.native_auto_strandedness_check/v1','passed':True,
             'scope':'Generated post-alignment FASTQ/BAM shared strand gate; no public/study dataset or complete upstream rerun',
             'actual_rseqc_image':args.rseqc_image,'cases':result_cases,
             'source_sha256':{str(p.relative_to(pipeline)):hashlib.sha256(p.read_bytes()).hexdigest()
                for p in [pipeline/'bin/katana_resolve_strandedness.py',pipeline/'subworkflows/local/resolve_strandedness.nf',
                          pipeline/'modules/local/infer_strandedness.nf',pipeline/'modules/local/resolve_sample_strandedness.nf']}}
    (root/'acceptance.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'passed':True,'cases':len(result_cases)}))


if __name__=='__main__':main()
