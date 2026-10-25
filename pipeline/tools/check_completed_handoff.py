#!/usr/bin/env python3
"""Verify production collection/report handoff of already-completed tasks.

This deliberately runs only KATANA_FLOW_COLLECT and KATANA_FLOW_REPORT. It does
not certify a full ``--mode complete`` run or execute an upstream/model task.
The input JSON has ``plan``, ``results`` (task-result directories), and optional
``models`` (external resource files). Paths may be relative to that JSON.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
from datetime import datetime, timezone


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')


def groovy(value):
    return "'" + str(value).replace('\\', '\\\\').replace("'", "\\'") + "'"


def bounded_inputs(path):
    path = Path(path).resolve()
    source = read(path)
    def located(value):
        value = Path(value)
        return (value if value.is_absolute() else path.parent / value).resolve()
    plan_path = located(source['plan'])
    plan = read(plan_path)
    results = sorted({located(p) for p in source['results']}, key=str)
    models = sorted({located(p) for p in source.get('models', [])}, key=str)
    expected = {task['id']: task for task in plan['tasks']}
    if len(expected) != len(plan['tasks']) or not expected:
        raise ValueError('The original sealed plan must have unique, nonempty task identities')
    observed, protected = {}, {path, plan_path, *models}
    for directory in results:
        receipt_path = directory / 'task.receipt.json'
        receipt = read(receipt_path)
        identity = receipt['task_id']
        if identity in observed or identity not in expected:
            raise ValueError('Duplicate or undeclared completed task')
        if receipt['plan_id'] != plan['plan_id'] or receipt['status'] not in {'completed', 'reused', 'not_applicable'}:
            raise ValueError('A source task has no matching successful completion')
        signed = {key: value for key, value in receipt.items() if key != 'receipt_sha256'}
        calculated = hashlib.sha256(json.dumps(signed, sort_keys=True, separators=(',', ':'), ensure_ascii=False, default=str).encode()).hexdigest()
        if calculated != receipt['receipt_sha256']:
            raise ValueError('Changed original task receipt: ' + str(receipt_path))
        observed[identity] = receipt_path
        protected.add(receipt_path)
        for record in receipt.get('inputs', []) + receipt.get('outputs', []) + receipt.get('producer_markers', []):
            candidate = Path(record['path'])
            # Image-owned /opt resources are checked by production collect;
            # never substitute a similarly named host resource here.
            if candidate.is_file():
                protected.add(candidate.resolve())
    if set(observed) != set(expected):
        raise ValueError('The original sealed plan is not fully covered')
    for model in models:
        if not model.is_file():
            raise FileNotFoundError(model)
    return {'plan': str(plan_path), 'results': list(map(str, results)), 'models': list(map(str, models))}, plan, protected


def inventory(paths):
    out = {}
    for p in sorted(paths, key=str):
        a = p.stat()
        sha = digest(p)
        b = p.stat()
        if (a.st_size, a.st_mtime_ns) != (b.st_size, b.st_mtime_ns):
            raise ValueError('Source changed during validation: ' + str(p))
        out[str(p)] = {'bytes': a.st_size, 'mtime_ns': a.st_mtime_ns, 'sha256': sha}
    return out


def trace(path):
    with Path(path).open(encoding='utf-8') as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    if len(rows) != 2 or {r['name'].split(' (', 1)[0] for r in rows} != {'KATANA_FLOW_COLLECT', 'KATANA_FLOW_REPORT'}:
        raise ValueError('The terminal handoff must contain exactly the two production processes')
    return rows


def work_products(work, rows):
    files = set()
    for row in rows:
        first, second = row['hash'].split('/')
        matches = [p for p in (work / first).glob(second + '*') if p.is_dir()]
        if len(matches) != 1:
            raise ValueError('Unable to identify one completed work directory')
        for product in ('collection', 'report'):
            root = matches[0] / product
            if root.is_dir():
                files.update(p for p in root.rglob('*') if p.is_file())
    return files


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', required=True, type=Path)
    parser.add_argument('--outdir', required=True, type=Path)
    parser.add_argument('--image', required=True, help='Already available full Katana Docker image, preferably exact ID')
    parser.add_argument('--pipeline-root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--nextflow-command-json', default='["nextflow"]', help='JSON argv; never shell-evaluated')
    parser.add_argument('--bind-readonly', action='append', default=[], type=Path,
                        help='Original absolute source roots needed by native task ledgers; mounted at unchanged paths')
    args = parser.parse_args(argv)
    out, pipeline = args.outdir.resolve(), args.pipeline_root.resolve()
    if out.exists():
        raise FileExistsError('Preserve an existing attempt and choose a new output directory')
    nextflow = json.loads(args.nextflow_command_json)
    if not isinstance(nextflow, list) or not nextflow or not all(isinstance(v, str) for v in nextflow):
        raise ValueError('Nextflow command must be a nonempty JSON string array')
    data, plan, protected = bounded_inputs(args.inputs)
    roots = [root.resolve() for root in args.bind_readonly]
    for root in roots:
        if not root.is_dir() or out == root or out.is_relative_to(root):
            raise ValueError('Original read-only roots must be independent of the new output')
    module = pipeline / 'modules/local/flow_control.nf'
    config = pipeline / 'nextflow.config'
    protected.update((module, config))
    before = inventory(protected)
    if shutil.disk_usage(out.parent).free < 1_000_000_000:
        raise ValueError('At least 1 GB headroom is required; size the retained products separately')
    out.mkdir(parents=True, exist_ok=False)
    write(out / 'source-inventory.before.json', before)
    write(out / 'inputs.json', data)
    source = """import groovy.json.JsonSlurper
include { KATANA_FLOW_COLLECT; KATANA_FLOW_REPORT } from %s
workflow {
    def data = new JsonSlurper().parse(file(params.completed_inputs).toFile())
    def original_plan = Channel.value(file(data.plan, checkIfExists: true))
    def products = Channel.value(data.results.collect { file(it, checkIfExists: true) })
    def resources = Channel.value(data.models.collect { file(it, checkIfExists: true) })
    KATANA_FLOW_COLLECT(original_plan, products, resources)
    KATANA_FLOW_REPORT(KATANA_FLOW_COLLECT.out.collection, resources)
}
""" % groovy(str(module.with_suffix('')))
    (out / 'main.nf').write_text(source, encoding='utf-8')
    mounts = []
    for root in roots:
        mounts.extend(('-v', shlex.quote(str(root) + ':' + str(root) + ':ro')))
    options = '--user $(id -u):$(id -g)' + (' ' + ' '.join(mounts) if mounts else '')
    (out / 'handoff.config').write_text(
        'includeConfig ' + groovy(config) + '\nprocess.executor = \'local\'\ndocker.writableInputMounts = false\ndocker.runOptions = ' + groovy(options) + '\n', encoding='utf-8')
    work = out / 'work'
    command = [*nextflow, 'run', str(out / 'main.nf'), '-c', str(out / 'handoff.config'),
               '-profile', 'docker', '--completed_inputs', str(out / 'inputs.json'),
               '--katana_container', args.image, '--outdir', str(out / 'results'),
               '-work-dir', str(work), '-ansi-log', 'false']
    lifecycle = {'schema': 'katana.completed_handoff_execution/v1', 'status': 'running',
                 'started_utc': datetime.now(timezone.utc).isoformat(), 'plan_id': plan['plan_id'],
                 'completed_source_tasks': len(plan['tasks']), 'scientific_tasks_submitted': 0,
                 'scope': 'Production COLLECT+REPORT only; not whole --mode complete acceptance', 'image': args.image}
    write(out / 'execution.json', lifecycle)
    completed_products = None
    try:
        for attempt in ('first', 'resume'):
            actual = [*command, '-with-trace', str(out / (attempt + '.trace.tsv'))]
            if attempt == 'resume':
                actual.append('-resume')
            write(out / (attempt + '.command.json'), actual)
            with (out / (attempt + '.stdout.log')).open('w') as stdout, (out / (attempt + '.stderr.log')).open('w') as stderr:
                child = subprocess.Popen(actual, cwd=out, stdout=stdout, stderr=stderr)
                lifecycle.update(attempt=attempt, pid=child.pid, command=actual)
                write(out / 'execution.json', lifecycle)
                code = child.wait()
            lifecycle[attempt + '_exit_code'] = code
            if code:
                raise RuntimeError('Nextflow ' + attempt + ' failed; preserve this attempt and inspect logs')
            rows = trace(out / (attempt + '.trace.tsv'))
            expected_status = 'COMPLETED' if attempt == 'first' else 'CACHED'
            if any(r['status'] != expected_status or r.get('exit', '0') not in ('0', '-') for r in rows):
                raise ValueError('Unexpected native terminal status: ' + str(rows))
            acceptance = read(out / 'results/report/flow.acceptance.json')
            if acceptance.get('status') != 'completed' or acceptance.get('plan_id') != plan['plan_id']:
                raise ValueError('Production report acceptance does not match the source plan')
            current = inventory(work_products(work, rows))
            if completed_products is None:
                completed_products = current
                write(out / 'completed-products.before-resume.json', current)
            elif current != completed_products:
                raise ValueError('Completed outputs changed during resume')
            if inventory(protected) != before:
                raise ValueError('Original source files or production modules changed')
        lifecycle.update(status='passed', finished_utc=datetime.now(timezone.utc).isoformat(),
                         first_completed_processes=2, resume_cached_processes=2,
                         original_inputs_and_completed_products_unchanged=True,
                         report=str(out / 'results/report/REPORT.html'),
                         report_acceptance=acceptance, protected_source_files=len(protected))
        write(out / 'acceptance.receipt.json', lifecycle)
    except Exception as exc:
        lifecycle.update(status='failed', finished_utc=datetime.now(timezone.utc).isoformat(), error=str(exc))
        write(out / 'failure.receipt.json', lifecycle)
        raise
    finally:
        write(out / 'execution.json', lifecycle)
    print(json.dumps({'status': lifecycle['status'], 'receipt': str(out / 'acceptance.receipt.json')}))


if __name__ == '__main__':
    main()
