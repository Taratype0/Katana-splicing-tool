#!/usr/bin/env python3
"""Exercise real resource selectors on tiny shell tasks; no biological data."""
import argparse
import importlib.util
import json
from pathlib import Path
import subprocess

# Resource rendering uses only the standard library. Loading it directly keeps
# this scheduler fixture independent of the scientific Python installation.
module_spec = importlib.util.spec_from_file_location('katana_resources',
    Path(__file__).resolve().parents[2] / 'katana_next/flow/resources.py')
resource_module = importlib.util.module_from_spec(module_spec)
module_spec.loader.exec_module(resource_module)
render_resources = resource_module.render_resources


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nextflow-command-file', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    args = parser.parse_args()
    root = args.outdir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    pipeline = Path(__file__).resolve().parents[1]
    command = json.loads(args.nextflow_command_file.read_text())
    proof = {}
    for name, options in [('portable', {}), ('site', {
        'defaults': {'queue': 'fallback'},
        'labels': {'katana_downstream': {'queue': 'shared-test'}},
        'processes': {'NETWORK_FIT': {'queue': 'parallel-test'}}
    })]:
        run_dir = root / name
        run_dir.mkdir()
        config = run_dir / 'resources.config'
        config.write_text(render_resources(options))
        run = subprocess.run([*command, 'run', str(pipeline / 'tests/resource_policy.nf'),
            '-c', str(pipeline / 'nextflow.config'), '-c', str(config),
            '-profile', 'local', '-ansi-log', 'false', '--outdir', str(run_dir / 'results'),
            '-work-dir', str(run_dir / 'work')], cwd=run_dir, capture_output=True, text=True, timeout=180)
        (run_dir / 'stdout.txt').write_text(run.stdout)
        (run_dir / 'stderr.txt').write_text(run.stderr)
        if run.returncode:
            raise AssertionError(run.stdout[-4000:] + run.stderr[-2000:])
        records = [p.read_text().strip().split('\t') for p in (run_dir / 'work').glob('*/*/allocation.tsv')]
        assert len(records) == 5, records
        rows = {r[0]: r[1:] for r in records}
        serial_queue = 'shared-test' if name == 'site' else 'unset'
        assert rows['gene-body'] == ['1', '4 GB', '12h', serial_queue], rows
        assert rows['sequencing-summary'] == ['1', '4 GB', '2h', serial_queue], rows
        assert rows['annotation'] == ['1', '4 GB', '2h', serial_queue], rows
        assert rows['parallel-fixture'] == ['2', '16 MB', '1m', 'parallel-test' if name == 'site' else 'unset'], rows
        assert rows['serial-upstream-fixture'] == ['1', '16 MB', '1m', serial_queue], rows
        proof[name] = rows
    receipt = {'status': 'passed', 'scope': 'native local allocation and queue-selector regression, not Slurm or scientific execution',
        'completed_tasks': 10, 'allocations': proof}
    (root / 'acceptance.receipt.json').write_text(json.dumps(receipt, indent=2))
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()
