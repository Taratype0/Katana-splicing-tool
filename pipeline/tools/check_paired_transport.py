#!/usr/bin/env python3
"""Exercise real Nextflow module argument transport using generated tool doubles.

This checks scheduling, paired flags and native module rendering, not statistical
models. pipeline/tests/paired_models.R separately runs the actual R models.
"""
import argparse
import csv
import json
import subprocess
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    runtime = ap.add_mutually_exclusive_group(required=True)
    runtime.add_argument('--nextflow-command-json')
    runtime.add_argument('--nextflow-command-file', type=Path)
    ap.add_argument('--outdir', required=True, type=Path)
    args = ap.parse_args()
    root = args.outdir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    source = Path(__file__).resolve().parents[1]
    bin_dir = root / 'bin'
    bin_dir.mkdir()
    scripts = {
        'rmats.py': '''#!/usr/bin/env python3
import json, pathlib, sys
if '--version' in sys.argv:
    print('v4.1.2'); sys.exit(0)
out=pathlib.Path(sys.argv[sys.argv.index('--od')+1])
out.mkdir(parents=True,exist_ok=True)
(out/'arguments.json').write_text(json.dumps(sys.argv[1:]))
''',
        'suppa.py': '''#!/usr/bin/env python3
import json, pathlib, sys
out=sys.argv[sys.argv.index('-o')+1]
pathlib.Path(out+'.dpsi').write_text(json.dumps(sys.argv[1:]))
pathlib.Path(out+'.psivec').write_text('generated transport check')
''',
        'python': '#!/bin/sh\nprintf "2.3\\n"\n',
    }
    for name, text in scripts.items():
        path = bin_dir / name
        path.write_text(text)
        path.chmod(0o755)
    inputs = root / 'inputs'
    inputs.mkdir()
    for name in ('ref.gtf', 'events.ioe', 'A.tpm', 'B.tpm', 'A.psi', 'B.psi', 'A.bam', 'B.bam', 'A.list', 'B.list', 'fixture.rmats'):
        (inputs / name).write_text('generated transport fixture\n')
    template = '''
include { RMATS_POST } from '__SOURCE__/modules/local/rmats_post'
include { DIFFSPLICE } from '__SOURCE__/modules/local/suppa_diffsplice'
workflow {
    p = { name -> file("${params.inputs}/${name}") }
    meta = [[subject:'mouse1',single_end:false,strandedness:'reverse']]
    flags = [paired:true, independent:false]
    rmats = Channel.fromList(['paired','independent']).map { id ->
        tuple(id, id+'A', id+'B', meta, meta, p.call('A.bam'),p.call('B.bam'),p.call('A.list'),p.call('B.list'),p.call('fixture.rmats'))
    }
    RMATS_POST(p.call('ref.gtf'),rmats,150,0.0001,false,50,500,flags)
    tpm = Channel.fromList(['paired','independent']).map { id -> tuple(id+'A',id+'B',p.call('A.tpm'),p.call('B.tpm')) }
    psi = Channel.fromList(['paired','independent']).map { id -> tuple(id+'A',id+'B',p.call('A.psi'),p.call('B.psi')) }
    suppaFlags = [('pairedA\\tpairedB'):true, ('independentA\\tindependentB'):false]
    DIFFSPLICE(p.call('events.ioe'),tpm,psi,'event','empirical',1000,0.05,0.05,1,0,true,suppaFlags,false)
}
'''
    (root / 'main.nf').write_text(template.replace('__SOURCE__', source.as_posix()))
    (root / 'nextflow.config').write_text('process.executor="local"\nprocess.cpus=1\n')
    runtime_json = args.nextflow_command_file.read_text() if args.nextflow_command_file else args.nextflow_command_json
    command = json.loads(runtime_json) + ['run', str(root / 'main.nf'),
        '--inputs', str(inputs), '-ansi-log', 'false', '-with-trace', str(root / 'trace.tsv')]
    run = subprocess.run(command, cwd=root, text=True, capture_output=True, timeout=180)
    (root / 'run.log').write_text(run.stdout + run.stderr)
    assert run.returncode == 0, run.stdout + run.stderr
    rows = list(csv.DictReader((root / 'trace.tsv').open(), delimiter='\t'))
    assert len(rows) == 4 and all(r['status'] == 'COMPLETED' for r in rows), rows
    records = []
    for pattern, flag in (('arguments.json', '--paired-stats'), ('*.dpsi', '-pa')):
        files = list((root / 'work').rglob(pattern))
        assert len(files) == 2, files
        for path in files:
            paired = 'independent' not in str(path.relative_to(root))
            argv = json.loads(path.read_text())
            assert (flag in argv) == paired, (path, argv)
            records.append({'file': str(path.relative_to(root)), 'paired': paired, 'flag': flag})
    (root / 'receipt.json').write_text(json.dumps({'status': 'passed', 'scope': 'native argument transport only',
        'statistical_model_validation': False, 'tasks': len(rows), 'records': records}, indent=2) + '\n')
    print('Paired/independent rMATS and SUPPA transport: 4 native tasks passed')


if __name__ == '__main__':
    main()
