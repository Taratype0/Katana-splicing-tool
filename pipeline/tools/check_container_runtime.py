#!/usr/bin/env python3
"""Run actual disconnected Nextflow/Docker installation checks, never science.

Requires an already built local image, Docker, and a standard Nextflow launcher.
Creates a new evidence directory outside the source checkout. No image is pulled,
no dependency is installed, and failed attempts are retained without resubmission.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


PIPELINE = Path(__file__).resolve().parents[1]
KINDS = {'network', 'mapp', 'motif'}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def inventory(root):
    return {path.relative_to(root).as_posix(): {'sha256': sha(path), 'size_bytes': path.stat().st_size}
            for path in sorted(root.rglob('*')) if path.is_file()}


def read_trace(path, expected_status, image):
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if len(rows) != 3:
        raise AssertionError(f'Expected exactly three independent image tasks; observed {len(rows)}')
    kinds = []
    for row in rows:
        matching = [kind for kind in KINDS if row.get('name', '').endswith('(' + kind + ')')]
        if len(matching) != 1:
            raise AssertionError('Unexpected image task name: ' + str(row))
        kinds.append(matching[0])
        if row.get('status') != expected_status or str(row.get('exit')) != '0':
            raise AssertionError('Image task did not reach the required successful state: ' + str(row))
        if row.get('container') != image or not Path(row.get('workdir', '')).is_dir():
            raise AssertionError('Image identity/work directory missing from actual Nextflow trace: ' + str(row))
    if set(kinds) != KINDS:
        raise AssertionError('Missing or repeated runtime probe kind: ' + str(kinds))
    return rows


def validate_products(root, uid, gid):
    reports = {}
    for kind in sorted(KINDS):
        report = json.loads((root / kind / 'acceptance.json').read_text(encoding='utf-8'))
        required = {'status': 'ready', 'probe_kind': kind, 'installed_prefix': '/opt/conda',
                    'task_uid': uid, 'task_gid': gid, 'host_runtime_rejected': True,
                    'existing_mounted_runtime_tested': True, 'host_runtime_executed': False,
                    'native_cli_runtime_check': True, 'scientific_computation': False}
        if any(report.get(key) != value for key, value in required.items()):
            raise AssertionError('Incomplete actual image/runtime evidence: ' + str(report))
        if 'external host path' not in ' '.join(report.get('host_runtime_error', [])):
            raise AssertionError('Mounted host interpreter was not rejected by the image boundary')
        if not (root / kind / 'katana-flow-help.stdout.txt').stat().st_size:
            raise AssertionError('Installed CLI help produced no output')
        reports[kind] = report
    return reports


def run(argv, root, label, timeout):
    with (root / (label + '.stdout.txt')).open('wb') as stdout, (root / (label + '.stderr.txt')).open('wb') as stderr:
        result = subprocess.run(argv, cwd=root, stdout=stdout, stderr=stderr, timeout=timeout, check=False)
    if result.returncode:
        raise RuntimeError(f'{label} exited {result.returncode}; retained stdout/stderr, task work and trace')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--image', required=True, help='Existing local Docker image name or sha256 ID; never pulled')
    parser.add_argument('--outdir', required=True, type=Path, help='New evidence directory outside the source checkout')
    parser.add_argument('--nextflow', default='nextflow', help='Standard Nextflow launcher executable/path (no shell string)')
    parser.add_argument('--uid', type=int, default=os.getuid() if hasattr(os, 'getuid') else None,
                        help='Non-root task UID; defaults to the Linux caller UID')
    parser.add_argument('--gid', type=int, default=os.getgid() if hasattr(os, 'getgid') else None,
                        help='Task GID; defaults to the Linux caller GID')
    parser.add_argument('--timeout', type=int, default=900, help='Per Nextflow invocation timeout in seconds')
    args = parser.parse_args()
    if os.name != 'posix' or args.uid is None or args.uid <= 0 or args.gid is None or args.gid < 0:
        parser.error('Use Linux/WSL with an explicit non-root UID and valid GID')
    if args.timeout <= 0:
        parser.error('--timeout must be positive')
    root = args.outdir.expanduser().resolve()
    if root == PIPELINE.parent or root.is_relative_to(PIPELINE.parent):
        parser.error('--outdir must be outside the source checkout')
    root.mkdir(parents=True, exist_ok=False)
    receipt = {'schema': 'katana.actual_nextflow_container_acceptance/v1', 'status': 'running',
               'started_utc': datetime.now(timezone.utc).isoformat(), 'requested_image': args.image,
               'task_uid': args.uid, 'task_gid': args.gid, 'container_network': 'none',
               'scientific_computation': False, 'study_inputs_used': False,
               'scope': 'Installed image runtime and actual Nextflow/Docker transport, not scientific FASTQ/BAM end-to-end acceptance'}
    def save():
        (root / 'acceptance.receipt.json').write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
    try:
        save()
        inspect = subprocess.run(['docker', 'image', 'inspect', '--format', '{{.Id}}', args.image],
                                 cwd=root, capture_output=True, text=True, timeout=30, check=False)
        (root / 'image-inspect.stdout.txt').write_text(inspect.stdout, encoding='utf-8')
        (root / 'image-inspect.stderr.txt').write_text(inspect.stderr, encoding='utf-8')
        image = inspect.stdout.strip()
        if inspect.returncode or not image.startswith('sha256:') or len(image) != 71:
            raise ValueError('A valid existing local image is required; no image will be pulled')
        receipt['actual_image_id'] = image
        fixture = root / 'fixture'
        shutil.copytree(PIPELINE / 'tests/container_runtime', fixture,
                        ignore=shutil.ignore_patterns('__pycache__', 'test_checker.py', 'README.md'))
        receipt['fixture_files'] = inventory(fixture)
        base = [args.nextflow, 'run', str(fixture / 'main.nf'), '--image', image,
                '--uid', str(args.uid), '--gid', str(args.gid), '--outdir', str(root / 'results'),
                '-ansi-log', 'false', '-work-dir', str(root / 'work')]
        receipt['first_command'] = [*base, '--trace_file', str(root / 'trace-first.tsv')]
        receipt['resume_command'] = [*base, '--trace_file', str(root / 'trace-resume.tsv'), '-resume']
        save()
        run(receipt['first_command'], root, 'first', args.timeout)
        first = read_trace(root / 'trace-first.tsv', 'COMPLETED', image)
        reports = validate_products(root / 'results', args.uid, args.gid)
        before = inventory(root / 'results')
        if not before:
            raise AssertionError('No actual image probe outputs were produced')
        work_before = {row['name']: inventory(Path(row['workdir']) / next(kind for kind in KINDS if row['name'].endswith('(' + kind + ')')))
                       for row in first}
        receipt.update(first_tasks=first, runtime_probes=reports, first_output_files=before)
        save()
        run(receipt['resume_command'], root, 'resume', args.timeout)
        resumed = read_trace(root / 'trace-resume.tsv', 'CACHED', image)
        if inventory(root / 'results') != before:
            raise AssertionError('Published output bytes changed on Nextflow resume')
        if {row['name']: inventory(Path(row['workdir']) / next(kind for kind in KINDS if row['name'].endswith('(' + kind + ')')))
            for row in resumed} != work_before:
            raise AssertionError('Task output bytes changed on Nextflow resume')
        validate_products(root / 'results', args.uid, args.gid)
        receipt.update(status='passed', resume_tasks=resumed, first_completed_count=3, resume_cached_count=3,
                       actual_container_execution=True, all_output_hashes_preserved=True,
                       actual_cli_and_native_runtime_probes=True, mounted_host_runtime_rejected=True)
    except Exception as exc:
        receipt.update(status='failed', error=f'{type(exc).__name__}: {exc}')
    finally:
        receipt['finished_utc'] = datetime.now(timezone.utc).isoformat()
        save()
    print(json.dumps({'status': receipt['status'], 'receipt': str(root / 'acceptance.receipt.json')}))
    return 0 if receipt['status'] == 'passed' else 1


if __name__ == '__main__':
    sys.exit(main())
