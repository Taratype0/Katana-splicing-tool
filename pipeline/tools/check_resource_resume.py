#!/usr/bin/env python3
"""Real Nextflow resource escalation/cache regression, generated shell tasks only."""
import argparse
import csv
import json
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nextflow-command-file', required=True, type=Path)
    parser.add_argument('--outdir', required=True, type=Path)
    args = parser.parse_args()
    root = args.outdir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    pipeline = Path(__file__).resolve().parents[1]
    prefix = json.loads(args.nextflow_command_file.read_text())
    base = [*prefix, 'run', str(pipeline / 'tests/resource_resume.nf'), '-c', str(pipeline / 'nextflow.config'),
            '-profile', 'local', '-ansi-log', 'false', '--fixture_marker', str(root / 'first-attempt.marker'),
            '-work-dir', str(root / 'work')]
    runs = []
    for name, options in [('initial', []), ('resume', ['-resume'])]:
        if name == 'resume':
            override = root / 'resource-only.config'
            override.write_text("process { withName: RESOURCE_CHECK { memory = '64 MB'; time = '2m' } }\n")
            options += ['-c', str(override)]
        trace = root / (name + '.trace.tsv')
        config = root / 'trace.config'
        config.write_text("trace.fields = 'task_id,hash,name,status,exit,attempt,cpus,memory,time'\n")
        run = subprocess.run([*base, '-c', str(config), '-with-trace', str(trace), *options],
                             cwd=root, capture_output=True, text=True, timeout=180)
        (root / (name + '.stdout.txt')).write_text(run.stdout)
        (root / (name + '.stderr.txt')).write_text(run.stderr)
        if run.returncode:
            raise AssertionError(run.stdout[-4000:] + run.stderr[-2000:])
        rows = list(csv.DictReader(trace.open(), delimiter='\t'))
        runs.append(rows)
    first, resumed = runs
    assert sum(r['status'] == 'COMPLETED' for r in first) == 2
    failed = [r for r in first if r['status'] == 'FAILED']
    assert len(failed) == 1 and failed[0]['exit'] == '137'
    retry = next(r for r in first if r['attempt'] == '2')
    assert retry['cpus'] == '1' and retry['memory'] == '32 MB'
    assert len(resumed) == 2 and all(r['status'] == 'CACHED' for r in resumed)
    proof = {'status': 'passed', 'scope': 'generated resource and resume tasks, not scientific validation',
             'completed': 2, 'failed_attempts': 1, 'actual_retry_memory': retry['memory'], 'cached_after_resource_change': 2}
    (root / 'acceptance.receipt.json').write_text(json.dumps(proof, indent=2))
    print(json.dumps(proof))


if __name__ == '__main__':
    main()
