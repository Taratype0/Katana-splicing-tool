#!/usr/bin/env python3
"""Real Nextflow acceptance using generated indexed BAMs, never study data."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys

REPO = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(REPO), str(REPO / 'tests/public_baseline')]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--nextflow-command-json', required=True)
    parser.add_argument('--outdir', required=True, type=Path)
    args = parser.parse_args()
    nf = json.loads(args.nextflow_command_json)
    root = args.outdir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    from test_report_workspace import report_fixture
    from katana_next.services.report_import import import_report
    from katana_next.targets.bundle import export_targets
    from katana_next.targets.nextflow import prepare_task
    from katana_next.targets.runner import inspect_results
    from katana_next.targets.restore import restore_results
    import pandas as pd

    report, fixture = report_fixture(root)
    source_hashes = {str(p.relative_to(report)): digest(p) for p in report.rglob('*') if p.is_file()}
    context = import_report(report, root / 'workspace')
    artifact = next(a for a in context.catalog.list_artifacts(context.manifest.project_id)
                    if a.name == 'All measured splice events')
    rows = pd.read_parquet(artifact.path)
    exported = export_targets(context, [
        {'source_artifact_id':artifact.artifact_id, 'row_ids':[rows.iloc[0].row_id], 'target_kind':'gene', 'filters':{}},
        {'source_artifact_id':artifact.artifact_id, 'row_ids':[rows.iloc[4].row_id], 'target_kind':'event', 'filters':{'include_nonsignificant':True}},
    ], root / 'request')
    prepared = prepare_task(exported['archive_path'], root / 'portable-task', fixture['bindings']['resources'])
    bindings = root / 'bindings.json'
    bindings.write_text(json.dumps(fixture['bindings'], indent=2))
    config = root / 'test.config'
    config.write_text('process.executor = "local"\nprocess.maxForks = 2\ntrace.overwrite = true\n')
    env = dict(os.environ, NXF_OFFLINE='true', MPLBACKEND='Agg', PYTHONDONTWRITEBYTECODE='1',
               PYTHONPATH=str(REPO), OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')

    def execute(label, main_nf, binding_file, out, work, resume=False):
        trace = root / f'{label}.trace.tsv'
        command = [*nf, 'run', str(main_nf), '-entry', 'TARGETS', '--target_request', exported['archive_path'],
                   '--target_bindings', str(binding_file), '--outdir', str(out), '-work-dir', str(work),
                   '-c', str(config), '-with-trace', str(trace), '-ansi-log', 'false']
        if resume:
            command.append('-resume')
        result = subprocess.run(command, cwd=root, env=env, capture_output=True, text=True, timeout=300)
        (root / f'{label}.stdout.txt').write_text(result.stdout)
        (root / f'{label}.stderr.txt').write_text(result.stderr)
        if result.returncode:
            raise AssertionError(f'{label} failed ({result.returncode})\n{result.stdout[-6000:]}\n{result.stderr[-3000:]}')
        with trace.open() as stream:
            tasks = list(csv.DictReader(stream, delimiter='\t'))
        assert len(tasks) == 3, tasks
        assert all('TARGET_REQUEST_INSPECT' in t['name'] or 'TARGET_DETAIL' in t['name'] for t in tasks), tasks
        return tasks

    outputs = root / 'native-results'
    first = execute('first', REPO / 'pipeline/main.nf', bindings, outputs, root / 'work')
    assert all(t['status'] == 'COMPLETED' for t in first), first
    targets = []
    for requested in exported['targets']:
        target_root = outputs / 'targets' / requested['target_id']
        result = inspect_results(target_root)
        assert result['verified'] is True, result
        target = next(t for t in result['targets'] if t['target_id'] == requested['target_id'])
        assert target['status'] == 'complete', target
        expected_events = 4 if requested['target_kind'] == 'gene' else 1
        assert len(target['sashimi']) == expected_events, target
        assert all(s['status'] == 'complete' for s in target['sashimi']), target
        assert list(target_root.rglob('*.png')), target_root
        targets.append({'target_id':requested['target_id'], 'kind':requested['target_kind'],
                        'status':target['status'], 'sashimi_events':expected_events,
                        'png_count':len(list(target_root.rglob('*.png'))), 'report':str(target_root / 'INDEX.html')})
    # Original result bytes must survive Nextflow resume and returned-result restore.
    before = {str(p.relative_to(outputs)):digest(p) for p in outputs.rglob('*') if p.is_file()}
    second = execute('resume', REPO / 'pipeline/main.nf', bindings, outputs, root / 'work', resume=True)
    assert all(t['status'] == 'CACHED' for t in second), second
    assert before == {str(p.relative_to(outputs)):digest(p) for p in outputs.rglob('*') if p.is_file()}
    returned = outputs / 'targets' / exported['targets'][1]['target_id']
    restored = restore_results(returned / 'REPORT.zip', root / 'restored', expected_project_id=context.manifest.project_id)
    assert inspect_results(root / 'restored')['verified'] is True, restored

    # The portable download uses the same module; missing indexes stay unavailable.
    missing = json.loads(bindings.read_text())
    index_id = next(r['resource_id'] for r in exported['resources'] if r['kind'] == 'bai')
    del missing['resources'][index_id]
    missing_file = root / 'missing-bindings.json'
    missing_file.write_text(json.dumps(missing))
    missing_out = root / 'missing-results'
    missing_tasks = execute('missing', Path(prepared['directory']) / 'main.nf', missing_file, missing_out, root / 'missing-work')
    for requested in exported['targets']:
        result = inspect_results(missing_out / 'targets' / requested['target_id'])
        target = next(t for t in result['targets'] if t['target_id'] == requested['target_id'])
        assert target['status'] == 'partial', target
        assert any(s['status'] == 'unavailable' for s in target['sashimi']), target
    # Rebinding the missing index creates fresh native work, then completed products.
    recovered = execute('recovered', Path(prepared['directory']) / 'main.nf', bindings, missing_out, root / 'missing-work', resume=True)
    assert sum(t['status'] == 'CACHED' for t in recovered) == 1, recovered
    for requested in exported['targets']:
        result = inspect_results(missing_out / 'targets' / requested['target_id'])
        assert next(t for t in result['targets'] if t['target_id'] == requested['target_id'])['status'] == 'complete'
    assert source_hashes == {str(p.relative_to(report)):digest(p) for p in report.rglob('*') if p.is_file()}
    summary = {'status':'passed','data':'generated synthetic indexed BAMs; no research data',
               'executor':'Nextflow local executor on Linux; not Slurm', 'targets':targets,
               'initial_completed':len(first),'resume_cached':len(second),
               'missing_index':'explicit partial/unavailable', 'recovered_completed':2,
               'source_report_unchanged':True,'returned_result_restore_verified':True}
    (root / 'acceptance.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
