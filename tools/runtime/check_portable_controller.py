"""Generated native handoff check of an installed, relocated controller."""
import csv
import json
import os
from pathlib import Path
import subprocess
import sys

from katana_next.flow.managed import _native_phase
from katana_next.flow.storage import create_managed_scope, register_stage, managed_run_lock

bundle, root = map(lambda x: Path(x).resolve(), sys.argv[1:3])
root.mkdir(parents=True, exist_ok=False)
scope = root / 'managed'
os.environ.update(NXF_HOME=str(root / 'nextflow-home'), NXF_OFFLINE='true',
                  PYTHONDONTWRITEBYTECODE='1')
source = root / 'generated.nf'
source.write_text('''nextflow.enable.dsl=2
process SEED {
 output:
 path 'seed.txt'
 script:
 """
 printf 'generated runtime check\\n' > seed.txt
 """
}
process FINISH {
 publishDir params.outdir, mode: 'copy'
 input:
 path seed
 output:
 path 'done.txt'
 script:
 """
 cat '${seed}' > done.txt
 """
}
workflow { SEED(); FINISH(SEED.out) }
''')
create_managed_scope(scope, protected_roots=(root / 'outputs',), run_id='independent-controller-check')
with managed_run_lock(scope):
    work = Path(register_stage(scope, 'upstream')['work_root'])
    command = [str(bundle / 'bin/nextflow'), 'run', str(source), '-ansi-log', 'false',
               '--outdir', str(root / 'outputs'), '-with-trace', str(root / 'first.trace.tsv')]
    first = _native_phase(command, work, scope, 'upstream')
    assert first['state'] == 'COMPLETED' and first['exit_code'] == 0
    command[-1] = str(root / 'resume.trace.tsv')
    second = _native_phase(command, work, scope, 'upstream', resume=True)
    assert second['state'] == 'COMPLETED' and second['exit_code'] == 0
for name, expected in [('first', 'COMPLETED'), ('resume', 'CACHED')]:
    rows = list(csv.DictReader((root / (name + '.trace.tsv')).open(), delimiter='\t'))
    assert len(rows) == 2 and all(r['status'] == expected for r in rows), rows
assert (root / 'outputs/done.txt').read_text() == 'generated runtime check\n'
assert Path(sys.executable).is_relative_to(bundle)
assert all(Path(p).resolve().is_relative_to(bundle) for p in sys.path)
subprocess.run([str(bundle / 'bin/katana-flow'), 'run', '--help'], check=True, stdout=subprocess.DEVNULL)
receipt = dict(status='passed', software_version='2.0.0b5.dev20261028',
               native_generated_tasks_completed=2, native_cached_tasks=2,
               installed_managed_implementation=True, relocated_python=True,
               old_environment_dependency=False, scientific_data_used=False,
               public_run_help_passed=True, slurm_job=os.environ.get('SLURM_JOB_ID'))
(root / 'acceptance.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps(receipt))
