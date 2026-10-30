#!/usr/bin/env python3
"""Native end-to-end delivery using generated analysis tables, never study data."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nextflow-command-file', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--pipeline', type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    from katana_next.full_fixture import create_full_fixture
    from katana_next.flow.storage_project import validate_retained_project
    root = args.outdir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    project = create_full_fixture(root / 'generated-input')
    # The general browser fixture has DTU but no DEXSeq input. Supply a tiny
    # explicit exon result for the full molecular report contract.
    data = json.loads(project.read_text())
    exon = project.parent / 'generated-exons.tsv'
    exon.write_text('groupID\tfeatureID\tpadj\tlog2fold\nG1\tE001\t0.01\t0.4\n')
    data['assets'].append(dict(asset_id='fixture-exons', kind='dexseq_result', path=str(exon),
                               metadata={'contrast_id': data['contrasts'][0]['contrast_id']}))
    project.write_text(json.dumps(data, indent=2))
    pipeline = args.pipeline.resolve()
    prefix = json.loads(args.nextflow_command_file.read_text())
    local = root / 'local.config'
    local.write_text("process { cpus=1; memory='1 GB'; time='10m'\n"
                     "withLabel:katana_control { cpus=1; memory='1 GB'; time='10m' }\n"
                     "withLabel:katana_report { cpus=1; memory='2 GB'; time='10m' }\n"
                     "withLabel:katana_downstream { cpus=1; memory='2 GB'; time='10m' }\n}\n")
    flags = ['--skip_' + name for name in ('network', 'mapp', 'motiftrace', 'motif_binding',
             'enrichment', 'rbp_identity', 'interaction', 'sequencing_qc', 'landscape', 'cross_comparison')]
    base = [*prefix, 'run', str(pipeline / 'main.nf'), '--mode', 'downstream-only',
            '--project', str(project), '--outdir', str(root / 'results'),
            '-c', str(local), '-ansi-log', 'false', *flags]
    records = {}
    for label, more in [('first', []), ('resume', ['-resume']), ('cleanup', ['--clean_work'])]:
        trace = root / (label + '.tsv')
        with (root / (label + '.log')).open('w') as stream:
            work = root / ('cleanup-work' if label == 'cleanup' else 'work')
            subprocess.run([*base, '-work-dir', str(work), '-with-trace', str(trace), *more], cwd=root,
                           stdout=stream, stderr=subprocess.STDOUT, check=True,
                           env={**os.environ, 'NXF_OFFLINE': 'true'}, timeout=600)
        rows = list(csv.DictReader(trace.open(), delimiter='\t'))
        expected = 'CACHED' if label == 'resume' else 'COMPLETED'
        assert rows and all(row['status'] == expected for row in rows), rows
        assert any(row['name'].endswith('KATANA_FLOW_RETAIN') for row in rows)
        assert any(row['name'].endswith('KATANA_FLOW_DELIVERY') for row in rows)
        records[label] = dict(tasks=len(rows), status=expected)
    receipt = json.loads((root / 'results/pipeline_info/delivery.receipt.json').read_text())
    assert receipt['status'] == 'completed'
    report = root / 'results/report/REPORT.html'
    digest = hashlib.sha256(report.read_bytes()).hexdigest()
    # Test-owned generated directories only: prove that the final client project
    # and report no longer depend on the source inputs or native work cache.
    assert not list((root / 'cleanup-work').rglob('.command.sh')), 'Native cleanup must remove unreferenced task work'
    shutil.rmtree(root / 'work')
    shutil.rmtree(root / 'generated-input')
    moved = root / 'relocated-results'
    (root / 'results').rename(moved)
    context = validate_retained_project(moved / 'katana/project')
    assert context.manifest.samples
    assert hashlib.sha256((moved / 'report/REPORT.html').read_bytes()).hexdigest() == digest
    assert (moved / 'report/REPORT.zip').stat().st_size > 0
    result = dict(status='passed', native_runs=records, relocated_project_readback=True,
                  generated_input_and_work_removed=True, nextflow_native_cleanup=True, report_bytes_unchanged=True,
                  scientific_study_rerun=False, optional_scientific_modules_tested=False)
    (root / 'acceptance.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
