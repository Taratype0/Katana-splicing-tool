#!/usr/bin/env python3
"""Real Nextflow cache boundaries; generated transport inputs, no study fits."""
import argparse
import csv
import json
from pathlib import Path

from check_native_transport import command, make_case, run, tasks


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nextflow-command-file',type=Path,required=True)
    parser.add_argument('--outdir',type=Path,required=True)
    args=parser.parse_args()
    nextflow=json.loads(args.nextflow_command_file.read_text())
    root=args.outdir.resolve();root.mkdir(parents=True,exist_ok=False)
    case=make_case(root/'scope')
    config=case/'configuration/flow.json'
    document=json.loads(config.read_text());document['task_cache_scope']='ancestry-v1'
    config.write_text(json.dumps(document))
    run(command(nextflow),case,'first')
    def resume(label):
        previous=set((case/'output/pipeline_info').glob('trace-*.tsv'))
        run([*command(nextflow),'-resume'],case,label)
        paths=set((case/'output/pipeline_info').glob('trace-*.tsv'))-previous
        assert len(paths)==1
        with next(iter(paths)).open() as handle:
            rows=list(csv.DictReader(handle,delimiter='\t'))
        return {row['name'].rsplit(' (',1)[-1].rstrip(')'):row['status'] for row in rows if ' (' in row['name']}
    document['report_options']={'fixture_unused_selection':True}
    config.write_text(json.dumps(document))
    unrelated=resume('unrelated-option')
    assert all(unrelated[t['id']]=='CACHED' for t in tasks()),unrelated
    next(t for t in document['tasks'] if t['id']=='delta')['parameters']={'fixture_parameter':2}
    config.write_text(json.dumps(document))
    scientific=resume('changed-task-parameter')
    rerun={'delta','eta','theta'}
    assert all(scientific[t['id']]==('COMPLETED' if t['id'] in rerun else 'CACHED') for t in tasks()),scientific
    receipt={'schema':'katana.native_cache_scope/v1','status':'passed','scientific_execution':False,
             'unrelated_selection':unrelated,'changed_parameter':scientific,
             'expected_affected_tasks':sorted(rerun)}
    (root/'acceptance.receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt))


if __name__=='__main__': main()
