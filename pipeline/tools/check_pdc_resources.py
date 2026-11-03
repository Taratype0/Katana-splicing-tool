#!/usr/bin/env python3
"""Render real Slurm wrappers using fake Slurm commands; never submit HPC jobs."""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
from check_upstream_resources import PIPELINE, fixture


def partitions_in(wrapper):
    # Slurm accepts several options on one #SBATCH line; clusterOptions uses
    # that form whereas the queue directive gets a dedicated line.
    directives = '\n'.join(line[len('#SBATCH'):].strip()
                           for line in wrapper.splitlines() if line.startswith('#SBATCH'))
    return re.findall(r'(?:^|\s)(?:-p\s+|--partition[=\s]+)(\S+)', directives)


def boundaries(root, command, env):
    proof = {}
    cases = [
        ('shared_edge', 54, '16h', [], 'shared'),
        ('main_edge', 56, '16h', [], 'main'),
        ('long_retry', 144, '32h', [], 'long'),
        ('unsupported_memory', 900, '16h', [], None),
        ('explicit_memory', 900, '16h', ['--slurm_partition', 'memory'], 'memory'),
        ('explicit_options', 900, '16h', ['--clusterOptions=--partition=custom-partition'], 'custom-partition'),
    ]
    for name, memory, duration, extra, partition in cases:
        run = root / name
        run.mkdir()
        (run/'main.nf').write_text('''process RESOURCE_BOUNDARY {
    output:
    path 'done.txt'
    script:
    """
    echo done > done.txt
    """
}
workflow { RESOURCE_BOUNDARY() }
''')
        (run/'probe.config').write_text(f"process {{ withName: 'RESOURCE_BOUNDARY' {{ cpus=12; memory={memory}.GB; time='{duration}'; container=null; publishDir=[enabled:false] }} }}\nexecutor.pollInterval='100ms'\nexecutor.queueStatInterval='100ms'\n")
        result = subprocess.run([*command, 'run', str(run/'main.nf'),
            '-c', str(PIPELINE/'nextflow.config'), '-c', str(run/'probe.config'),
            '-profile', 'pdc_kth', '--project', 'generated-account', *extra,
            '--output_layout', 'legacy', '--outdir', str(run/'results'),
            '-work-dir', str(run/'work'), '-ansi-log', 'false'],
            cwd=run, env=env, capture_output=True, text=True, timeout=120)
        (run/'stdout.txt').write_text(result.stdout)
        (run/'stderr.txt').write_text(result.stderr)
        if partition is None:
            assert result.returncode != 0
            assert 'Katana never selects memory automatically' in result.stdout+result.stderr
            assert not list((run/'work').glob('*/*/done.txt'))
        else:
            assert result.returncode == 0, result.stdout[-4000:]+result.stderr[-2000:]
            wrappers = list((run/'work').glob('*/*/.command.run'))
            assert len(wrappers) == 1
            text = wrappers[0].read_text()
            found = partitions_in(text)
            assert found == [partition], (name, found, text[:800])
        proof[name] = dict(memory_gib=memory, time=duration, partition=partition,
                           returncode=result.returncode, extra=extra)
    return proof


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nextflow-command-file', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    args = parser.parse_args()
    command = json.loads(args.nextflow_command_file.read_text())
    catalog = json.loads((PIPELINE / 'assets/resource-policy.json').read_text())
    root = args.outdir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    fake = root / 'bin'
    fake.mkdir()
    # Generated tasks contain only printf. Complete synchronously; the native
    # Nextflow executor still writes and submits each real Slurm wrapper.
    (fake / 'sbatch').write_text('''#!/bin/bash
bash "${@: -1}" >/dev/null 2>&1
echo "Submitted batch job $$"
exit 0
''')
    for name in ('squeue', 'sacct', 'scancel', 'module'):
        (fake / name).write_text('#!/bin/sh\nexit 0\n')
    for path in fake.iterdir():
        path.chmod(0o755)
    env = {**os.environ, 'PATH': str(fake) + os.pathsep + os.environ['PATH']}
    proof = {}
    for case in ('automatic', 'user_override'):
        run = root / case
        expected = fixture(run, catalog)
        config = run / 'probe.config'
        # No fake resource cap, no actual scheduler, no container/science.
        config.write_text("process { withName: '.*' { container = null; publishDir = [enabled: false] } }\nexecutor.pollInterval='100ms'\nexecutor.queueStatInterval='100ms'\n")
        extra = []
        if case == 'user_override':
            user = run / 'user.config'
            user.write_text("process { withName: 'STAR_ALIGN' { queue='user-star'; memory=80.GB; cpus=10 } }\n")
            extra = ['-c', str(user)]
        result = subprocess.run([*command, 'run', str(run/'main.nf'),
            '-c', str(PIPELINE/'nextflow.config'), '-c', str(config), *extra,
            '-profile', 'pdc_kth', '--project', 'generated-account',
            '--output_layout', 'legacy', '--outdir', str(run/'results'),
            '-work-dir', str(run/'work'), '-ansi-log', 'false'],
            cwd=run, env=env, capture_output=True, text=True, timeout=240)
        (run/'stdout.txt').write_text(result.stdout)
        (run/'stderr.txt').write_text(result.stderr)
        assert result.returncode == 0, result.stdout[-5000:]+result.stderr[-3000:]
        rows = {}
        for path in (run/'work').glob('*/*/allocation.tsv'):
            fields = path.read_text().strip().split('\t')
            name, cpu, memory, time, queue, attempt = fields
            wrapper = (path.parent/'.command.run').read_text()
            cpus, gib = int(cpu), float(memory.split()[0])
            # All inventory cases use GiB and <=24h initial/retry durations.
            inflated = max(cpus, __import__('math').ceil(gib*1024/888))
            desired = 'shared' if inflated <= 64 and gib <= 111 else ('main' if time != '1d 8h' else 'long')
            if name == 'KATANA_FLOW_CONFIGURE':
                assert '#SBATCH' not in wrapper, wrapper[:500]
                rows[name] = dict(executor='local')
                continue
            if case == 'user_override' and name == 'STAR_ALIGN':
                desired = 'user-star'
                assert cpus == 10 and gib == 80
            assert queue == desired, (case, name, fields, desired)
            partitions = partitions_in(wrapper)
            assert partitions == [desired], (name, partitions, wrapper[:700])
            assert '--account=generated-account' in wrapper
            rows.setdefault(name, []).append(dict(cpus=cpus, memory=memory, time=time,
                attempt=int(attempt), partition=queue, shared_equivalent_cpus=inflated))
        assert set(rows) == set(expected), set(rows)^set(expected)
        for name in ('STAR_ALIGN', 'NETWORK_FIT', 'RMATS_POST', 'STAR_GENOMEGENERATE'):
            assert rows[name][0]['partition'] == ('user-star' if case=='user_override' and name=='STAR_ALIGN' else 'main')
        for name in ('TRIMGALORE','FASTQC','DEXSEQ_COUNT','MOTIFTRACE','MAPP','QC','KATANA_FLOW_REPORT'):
            assert all(r['partition']=='shared' for r in rows[name]), (name, rows[name])
        proof[case] = rows
    receipt = dict(status='passed', scope='Native Nextflow Slurm wrappers; fake scheduler; no scientific jobs',
                   declared_processes=len(catalog['processes']), allocations=proof,
                   boundaries=boundaries(root, command, env))
    (root/'acceptance.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({'status':'passed','cases':list(proof),'tasks_per_case':len(expected)}))


if __name__ == '__main__':
    main()
