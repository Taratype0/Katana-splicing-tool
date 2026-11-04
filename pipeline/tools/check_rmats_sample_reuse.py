#!/usr/bin/env python3
"""Linux native DAG check: per-BAM reuse, pairing, single condition and resume.

The tool double checks channel/file identities, not rMATS numerical results.
"""
import argparse
import csv
import json
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nextflow-command-file', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    args = parser.parse_args()
    source = Path(__file__).resolve().parents[1]
    command = json.loads(args.nextflow_command_file.read_text())
    root = args.outdir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    proof = {}
    for case in ('contrasts', 'single', 'subset'):
        run = root / case
        (run / 'bin').mkdir(parents=True)
        inputs = run / 'inputs'
        inputs.mkdir()
        (inputs / 'ref.gtf').write_text('generated fixture\n')
        for name in ('A1', 'A2', 'A3', 'B1', 'B2'):
            (inputs / f'{name}.bam').write_text('generated fixture\n')
        (inputs / 'contrasts.csv').write_text('contrast,treatment,control,exclude_samples\n' +
            ('A_profile,A,A,\n' if case == 'single' else
             'A_vs_B,A,B,A3\nB_vs_A,B,A,\n' if case == 'subset' else 'A_vs_B,A,B,\nB_vs_A,B,A,\n'))
        fake = run / 'bin/rmats.py'
        fake.write_text('''#!/usr/bin/env python3
import json, pathlib, sys
argv=sys.argv[1:]
if '--version' in argv:
    print('v4.1.2'); sys.exit(0)
def value(key): return argv[argv.index(key)+1]
tmp=pathlib.Path(value('--tmp')); tmp.mkdir(parents=True, exist_ok=True)
out=pathlib.Path(value('--od')); out.mkdir(parents=True, exist_ok=True)
arms=[pathlib.Path(value(key)).read_text().strip().split(',') for key in ('--b1','--b2') if key in argv]
if value('--task')=='prep':
    assert len(arms)==1 and len(arms[0])==1, arms
    (tmp/'timestamp.rmats').write_text(json.dumps(arms[0]))
else:
    records=sorted(tmp.glob('*.rmats'))
    prepared=[bam for record in records for bam in json.loads(record.read_text())]
    expected=[bam for arm in arms for bam in arm]
    assert sorted(prepared)==sorted(expected), (prepared,expected)
    assert len(prepared)==len(set(prepared)), prepared
    (out/'arguments.json').write_text(json.dumps({'arms':arms, 'paired':'--paired-stats' in argv, 'statoff':'--statoff' in argv}))
''')
        fake.chmod(0o755)
        template = '''
include { RMATS } from '__SOURCE__/subworkflows/local/rmats'
workflow {
    p = { name -> file("${params.inputs}/${name}") }
    make = { condition ->
        def ids = __IDS__
        def metas = ids.collect { n -> [id:condition+n, subject:'mouse'+n, single_end:false, strandedness:'reverse'] }
        tuple(condition, metas, metas.collect { m -> p.call(m.id+'.bam') })
    }
    grouped = Channel.fromList(__GROUPS__.collect { make.call(it) })
    RMATS(Channel.value(p.call('contrasts.csv')), grouped, Channel.value(p.call('ref.gtf')),
        __SINGLE__,150,0.0001,false,50,500,[A_vs_B:true,B_vs_A:false,A_profile:false])
}
'''
        (run / 'main.nf').write_text(template.replace('__SOURCE__', source.as_posix())
            .replace('__IDS__', "condition == 'A' ? [1,2,3] : [1,2]" if case == 'subset' else '[1,2]')
            .replace('__GROUPS__', "['A']" if case == 'single' else "['A','B']")
            .replace('__SINGLE__', str(case == 'single').lower()))
        (run / 'nextflow.config').write_text('process.executor="local"\nprocess.cpus=1\n')
        trace_files = []
        for attempt in ('first', 'resume'):
            trace = run / f'{attempt}.tsv'
            invocation = command + ['run', str(run / 'main.nf'), '--inputs', str(inputs),
                '-ansi-log', 'false', '-with-trace', str(trace)]
            if attempt == 'resume': invocation.append('-resume')
            result = subprocess.run(invocation, cwd=run, text=True, capture_output=True, timeout=180)
            (run / f'{attempt}.log').write_text(result.stdout + result.stderr)
            assert result.returncode == 0, result.stdout + result.stderr
            rows = list(csv.DictReader(trace.open(), delimiter='\t'))
            assert rows and all(row['status'] == ('CACHED' if attempt == 'resume' else 'COMPLETED') for row in rows), rows
            preps = [row for row in rows if ':RMATS_PREP_BAM ' in row['name']]
            assert len(preps) == (2 if case == 'single' else 5 if case == 'subset' else 4), rows
            trace_files.append({'attempt':attempt,'tasks':len(rows),'prepared_bams':len(preps)})
        observations = [json.loads(path.read_text()) for path in (run / 'work').rglob('arguments.json')]
        if case == 'single':
            assert len(observations) == 1 and observations[0]['statoff'] and not observations[0]['paired']
        else:
            assert len(observations) == 2, observations
            paired = next(row for row in observations if row['paired'])
            assert paired['arms'] == [['A1.bam','A2.bam'],['B1.bam','B2.bam']], paired
            if case == 'subset':
                independent = next(row for row in observations if not row['paired'])
                assert independent['arms'] == [['B1.bam','B2.bam'],['A1.bam','A2.bam','A3.bam']], independent
        proof[case] = trace_files
    (root / 'receipt.json').write_text(json.dumps({'status':'passed','scope':'native DAG transport; no study data', 'cases':proof}, indent=2))


if __name__ == '__main__':
    main()
