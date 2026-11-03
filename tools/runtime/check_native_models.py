#!/usr/bin/env python3
"""Exercise installed CoSplice, MotifTrace and delivery through native Nextflow.

Uses generated inputs only. Run with the relocated controller Python and a
matching local Docker image. The source checkout supplies test data generation,
never the analysis imports used inside task containers.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--image', required=True)
    parser.add_argument('--nextflow', required=True)
    args = parser.parse_args()
    root = args.outdir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    source = Path(__file__).resolve().parents[2]
    # Create all five event types, two contrasts and two known motif models.
    setup = '''import sys,json
from pathlib import Path
sys.path.insert(0,"/checks")
from check_unified_engines import fixture
p,m=fixture(Path.cwd()/"source")
cfg={"design_factors":["condition"],"contrast_ids":["perturbation-vs-baseline","rescue-vs-baseline"],"analyses":{"coverage":False,"mapp":False,"motif_binding":False,"interaction":False,"sequencing_qc":False},"parameters":{"prioritize.motiftrace_consequences":{"custom_motif_path":str(m),"cluster_permutations":99}}}
Path("config").mkdir()
Path("config/flow.json").write_text(json.dumps(cfg))
'''
    subprocess.run(['docker', 'run', '--rm', '--network=none', '--user',
        f'{os.getuid()}:{os.getgid()}', '-v', f'{root}:{root}', '-v',
        f'{source / "scripts"}:/checks:ro', '-w', str(root), args.image,
        'python', '-c', setup], check=True)
    from katana_next.flow.distribution import export_pipeline
    from katana_next.flow.storage_project import validate_retained_project
    export_pipeline(root / 'pipeline')
    (root / 'acceptance.config').write_text('executor.cpus = 4\nexecutor.memory = "6 GB"\n')
    base = [args.nextflow, 'run', str(root / 'pipeline/main.nf'), '-profile', 'docker',
        '--mode', 'downstream-only', '--project', str(root / 'source/katana.project.json'),
        '--flow_config', str(root / 'config/flow.json'), '--outdir', str(root / 'results'),
        '--storage_mode', 'lean', '--katana_container', args.image,
        '--max_cpus', '2', '--max_memory', '4 GB', '-c', str(root / 'acceptance.config'),
        '-work-dir', str(root / 'work'), '-ansi-log', 'false']
    records = {}
    for label, extra, expected in [('first', [], 'COMPLETED'), ('resume', ['-resume'], 'CACHED')]:
        log = root / f'{label}.log'
        trace = root / f'{label}.tsv'
        with log.open('w') as stream:
            run = subprocess.run([*base, '-with-trace', str(trace), *extra], cwd=root,
                stdout=stream, stderr=subprocess.STDOUT, timeout=3600)
        if run.returncode:
            print(log.read_text()[-24000:], flush=True)
            raise RuntimeError(f'{label} Nextflow failed: {run.returncode}')
        with trace.open() as stream:
            rows = list(csv.DictReader(stream, delimiter='\t'))
        assert rows and all(r['status'] == expected for r in rows), rows
        records[label] = {'tasks': len(rows), 'status': expected}
    results = root / 'results'
    receipt = json.loads((results / 'pipeline_info/delivery.receipt.json').read_text())
    assert receipt['status'] == 'completed', receipt
    validate_retained_project(results / 'katana/project')
    report = results / 'report'
    manifest = json.loads((report / 'report-manifest.json').read_text())
    assert manifest['presentation']['chapter_order'] == [1, 2, 3, 4, 5, 7, 6]
    assert manifest['network_browsers'][0]['whole_study'] is True
    assert (report / 'REPORT.html').stat().st_size > 0
    assert (report / 'REPORT.zip').stat().st_size > 0
    (root / 'acceptance.json').write_text(json.dumps(dict(status='passed',
        generated_data_only=True, native_runs=records, retained_project_readback=True,
        unchanged_report_format=True, network_and_motiftrace_enabled=True), indent=2) + '\n')


if __name__ == '__main__':
    main()
