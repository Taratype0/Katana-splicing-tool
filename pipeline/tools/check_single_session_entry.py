#!/usr/bin/env python3
"""Verify production entry wiring with an explicit non-biological upstream double."""
import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--nextflow-command-file', type=Path, required=True)
    p.add_argument('--outdir', type=Path, required=True)
    args = p.parse_args()
    source = Path(__file__).resolve().parents[1]
    root = args.outdir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    pipeline = root / 'pipeline'
    shutil.copytree(source, pipeline)
    # Only the upstream tool implementation is replaced in this test copy.
    # The public main, import and downstream scheduling code are unchanged.
    (pipeline / 'workflows/katana_core.nf').write_text('''process GENERATED_UPSTREAM {
    input:
    path samples
    path contrasts
    path fasta
    path gtf
    output:
    path 'katana_bundle'
    script:
    """
    test -s "$samples"
    test -s "$contrasts"
    test -s "$fasta"
    test -s "$gtf"
    mkdir katana_bundle
    printf 'generated transport only' > katana_bundle/fixture.txt
    """
}
workflow KATANA_PIPELINE {
    take:
    samples
    contrasts
    fasta
    gtf
    main:
    GENERATED_UPSTREAM(samples, contrasts, fasta, gtf)
    emit:
    bundle = GENERATED_UPSTREAM.out
}
''')
    prefix = json.loads(args.nextflow_command_file.read_text())
    flags = ['--skip_' + name for name in ('network', 'mapp', 'motiftrace', 'motif_binding', 'enrichment', 'rbp_identity')]
    cases = {}
    for entry in ('fastq', 'bam', 'bundle'):
        case = root / entry
        case.mkdir()
        command = [*prefix, 'run', str(pipeline / 'main.nf'), '-c', str(pipeline / 'tests/transport.config'),
                   '--flow_executable', str(pipeline / 'tests/protocol_fixture.py'),
                   '--outdir', str(case / 'results'), '-work-dir', str(case / 'work'), '-ansi-log', 'false', *flags]
        if entry == 'bundle':
            bundle = case / 'existing'
            bundle.mkdir()
            (bundle / 'fixture.txt').write_text('generated transport only')
            command += ['--mode', 'downstream-only', '--bundle', str(bundle)]
        else:
            for name in ('samples.csv', 'contrasts.csv', 'genome.fa', 'genes.gtf'):
                (case / name).write_text('generated ' + entry + ' transport fixture\n')
            command += ['--mode', 'complete', '--input', str(case / 'samples.csv'), '--contrasts', str(case / 'contrasts.csv'),
                        '--fasta', str(case / 'genome.fa'), '--gtf', str(case / 'genes.gtf')]
        for attempt in ('first', 'resume'):
            trace = case / (attempt + '.tsv')
            with (case / (attempt + '.log')).open('w') as log:
                subprocess.run([*command, '-with-trace', str(trace), *(['-resume'] if attempt == 'resume' else [])],
                               cwd=case, env={**os.environ, 'NXF_OFFLINE': 'true'},
                               stdout=log, stderr=subprocess.STDOUT, check=True, timeout=180)
            rows = list(csv.DictReader(trace.open(), delimiter='\t'))
            assert rows and all(r['status'] == ('CACHED' if attempt == 'resume' else 'COMPLETED') for r in rows)
            names = [r['name'] for r in rows]
            assert any(n.endswith('KATANA_FLOW_IMPORT') for n in names)
            assert any(n.endswith('KATANA_FLOW_REPORT') for n in names)
            assert any('GENERATED_UPSTREAM' in n for n in names) == (entry != 'bundle')
            cases[entry + '_' + attempt] = len(rows)
    result = dict(status='passed', cases=cases, single_nextflow_session_per_run=True,
                  scope='Production entry/import/downstream wiring only; upstream science replaced by an explicit generated fixture')
    (root / 'acceptance.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
