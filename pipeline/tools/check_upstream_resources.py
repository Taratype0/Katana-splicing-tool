#!/usr/bin/env python3
"""Render the real allocation policy on empty native tasks, never study data.

The local executor's capacity is deliberately simulated: these tasks only print
their requests. No large allocation, tool, container or cluster job is started.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import re
import subprocess


PIPELINE = Path(__file__).resolve().parents[1]


def declared_processes():
    result = {}
    for path in (PIPELINE / 'modules').rglob('*.nf'):
        text = path.read_text()
        for match in re.finditer(r'process\s+(\w+)\s*\{', text):
            header = text[match.end():].split('input:', 1)[0]
            result[match[1]] = (path, re.findall(r"label\s+['\"]([^'\"]+)", header))
    return result


def declared_aliases():
    aliases = {}
    for folder in ('workflows', 'subworkflows'):
        for path in (PIPELINE / folder).rglob('*.nf'):
            for block in re.findall(r'include\s*\{(.*?)\}', path.read_text(), re.S):
                for original, alias in re.findall(r'(\w+)\s+as\s+(\w+)', block):
                    if original in declared_processes() and original != alias:
                        assert alias not in aliases or aliases[alias] == original, alias
                        aliases[alias] = original
    return aliases


def fixture(root, catalog):
    root.mkdir(parents=True, exist_ok=False)
    modules = root / 'probes'
    modules.mkdir()
    declared = declared_processes()
    assert set(declared) == {r['process'] for r in catalog['processes']}, 'Unreviewed or stale process inventory'
    include, calls, expected = [], [], {}
    for record in catalog['processes']:
        name = record['process']
        assert declared[name][1] == record['labels'], name
        labels = '\n'.join(f"    label '{label}'" for label in record['labels'])
        text = '''process NAME {
LABELS
    input:
    val spec
    output:
    path 'allocation.tsv'
    script:
    """
    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' '${spec.id}' '${task.cpus}' '${task.memory}' '${task.time}' '${task.queue ?: "unset"}' '${task.attempt}' > allocation.tsv
    if [ '${spec.retry ?: false}' = true ] && [ '${task.attempt}' = 1 ]; then exit 137; fi
    """
}
'''.replace('NAME', name).replace('LABELS', labels)
        if name == 'RMATS_POST':
            text = text.replace('    val spec\n', '    val spec\n    val contrast\n    val rmats_paired_stats\n')
        (modules / f'{name}.nf').write_text(text)
        include.append(f"include {{ {name} }} from './probes/{name}'")
        extras = ", Channel.value('probe'), Channel.value([probe:false])" if name == 'RMATS_POST' else ''
        calls.append(f"    {name}(Channel.of([id:'{name}', execution_group:'summary', workflow_id:'fixture', retry:{str(name in ('TRIMGALORE','DIFFSPLICE')).lower()}]){extras})")
        value = record['initial']
        expected[name] = ([2, '12 GB', '4h', 'unset', 1] if value['cpus'] == 'task-dependent' else
                          [int(value['cpus']), f"{value['memory_gib']} GB", f"{value['time_hours']}h", 'unset', 1])
    # Every currently wired alias, including all downstream execution groups.
    for alias, original in declared_aliases().items():
        if alias in expected:
            continue
        include.append(f"include {{ {original} as {alias} }} from './probes/{original}'")
        group = 'motiftrace' if alias == 'MOTIFTRACE_RESOURCE' else alias.lower()
        calls.append(f"    {alias}(Channel.of([id:'{alias}', execution_group:'{group}', workflow_id:'fixture']))")
        expected[alias] = expected[original].copy()
        if original in ('KATANA_FLOW_TASK', 'KATANA_FLOW_RESOURCE_TASK'):
            setting = catalog['downstream_groups'][group]
            expected[alias] = [setting['cpus'], f"{setting['memory_gib']} GB", f"{setting['time_hours']}h", 'unset', 1]
    include.append("include { RMATS_POST as RMATS_POST_PAIRED_PROBE } from './probes/RMATS_POST'")
    calls.append("    RMATS_POST_PAIRED_PROBE(Channel.of([id:'RMATS_POST_PAIRED_PROBE', retry:true]), Channel.value('matched'), Channel.value([matched:true]))")
    expected['RMATS_POST_PAIRED_PROBE'] = [12, '16 GB', '16h', 'unset', 1]
    (root / 'main.nf').write_text('nextflow.enable.dsl=2\n' + '\n'.join(include) +
        '\nworkflow KATANA_CORE {\n' + '\n'.join(calls) +
        '\n}\nworkflow KATANA_PIPELINE { KATANA_CORE() }\nworkflow { KATANA_PIPELINE() }\n')
    (root / 'probe.config').write_text('''// Capacity simulation for print-only tasks, not production settings.
executor.cpus = 1024
executor.memory = '4 TB'
process { withName: '.*' { container = null; publishDir = [enabled: false] } }
''')
    return expected


def memory_only_resume(root, command):
    """A resource reduction must not invalidate an unchanged successful task."""
    root.mkdir()
    (root / 'payload.txt').write_text('generated resource-cache fixture\n')
    (root / 'main.nf').write_text('''process DIFFSPLICE {
    label 'process_high'
    input:
    path payload
    output:
    path 'result.txt'
    script:
    """
    cat "$payload" > result.txt
    """
}
workflow { DIFFSPLICE(Channel.value(file("${projectDir}/payload.txt"))) }
''')
    (root / 'capacity.config').write_text("executor.cpus=1024\nexecutor.memory='4 TB'\n")
    (root / 'previous.config').write_text("process { withName: 'DIFFSPLICE' { memory=72.GB; time=16.h } }\n")
    records = {}
    for case in ('previous', 'measured'):
        trace = root / f'{case}.trace.tsv'
        extra = ['-c', str(root/'previous.config')] if case == 'previous' else ['-resume']
        result = subprocess.run([*command, 'run', str(root/'main.nf'),
            '-c', str(PIPELINE/'nextflow.config'), '-c', str(root/'capacity.config'), *extra,
            '-profile', 'local', '--output_layout', 'legacy', '--outdir', str(root/'results'),
            '-work-dir', str(root/'work'), '-with-trace', str(trace), '-ansi-log', 'false'],
            cwd=root, text=True, capture_output=True, timeout=120)
        (root/f'{case}.stdout.txt').write_text(result.stdout)
        (root/f'{case}.stderr.txt').write_text(result.stderr)
        assert result.returncode == 0, result.stdout[-3000:] + result.stderr[-2000:]
        with trace.open() as stream:
            rows = list(csv.DictReader(stream, delimiter='\t'))
        assert len(rows) == 1, rows
        assert rows[0]['status'] == ('COMPLETED' if case == 'previous' else 'CACHED'), rows
        records[case] = rows[0]
    assert records['previous']['hash'] == records['measured']['hash'], records
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nextflow-command-file', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    args = parser.parse_args()
    command = json.loads(args.nextflow_command_file.read_text())
    catalog = json.loads((PIPELINE / 'assets/resource-policy.json').read_text())
    root = args.outdir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    proof = {}
    for case in ('defaults', 'user_overrides'):
        run = root / case
        expected = fixture(run, catalog)
        expected['TRIMGALORE'] = [8, '16 GB', '16h', 'unset', 2]
        expected['DIFFSPLICE'] = [1, '8 GB', '4h', 'unset', 2]
        expected['RMATS_POST_PAIRED_PROBE'] = [12, '32 GB', '1d 8h', 'unset', 2]
        extra = []
        if case == 'user_overrides':
            user = run / 'user.config'
            # Synthetic configuration only; never reads the user's custom.config.
            user.write_text('''process {
    withName: 'TRIMGALORE' { cpus = 3; memory = '5 GB'; time = '3h'; queue = 'user-trim' }
    withName: 'KATANA_PIPELINE:KATANA_CORE:KATANA_DESEQ2_CONTRAST' {
        cpus = 2; memory = '9 GB'; time = '7h'; queue = 'user-statistics'
    }
    withName: 'NETWORK_FIT' { cpus = 4; memory = '24 GB'; time = '6h'; queue = 'user-network' }
}''')
            extra = ['-c', str(user)]
            expected['TRIMGALORE'] = [3, '5 GB', '3h', 'user-trim', 2]
            expected['KATANA_DESEQ2_CONTRAST'] = [2, '9 GB', '7h', 'user-statistics', 1]
            expected['NETWORK_FIT'] = [4, '24 GB', '6h', 'user-network', 1]
        result = subprocess.run([*command, 'run', str(run / 'main.nf'),
            '-c', str(PIPELINE / 'nextflow.config'), '-c', str(run / 'probe.config'), *extra,
            '-profile', 'local', '-ansi-log', 'false', '--output_layout', 'legacy',
            '--outdir', str(run / 'results'), '-work-dir', str(run / 'work')],
            cwd=run, text=True, capture_output=True, timeout=240)
        (run / 'stdout.txt').write_text(result.stdout)
        (run / 'stderr.txt').write_text(result.stderr)
        if result.returncode:
            raise AssertionError(result.stdout[-6000:] + result.stderr[-3000:])
        rows = {}
        for path in (run / 'work').glob('*/*/allocation.tsv'):
            fields = path.read_text().strip().split('\t')
            key = fields[0]
            value = [int(fields[1]), fields[2], fields[3], fields[4], int(fields[5])]
            rows.setdefault(key, []).append(value)
        assert set(rows) == set(expected), (set(rows) ^ set(expected))
        for key, value in expected.items():
            observed = sorted(rows[key], key=lambda row: row[-1])
            assert observed[-1] == value, (case, key, observed, value)
            if key == 'TRIMGALORE':
                assert len(observed) == 2, observed
                assert observed[0] == ([8, '8 GB', '8h', 'unset', 1] if case == 'defaults' else [3, '5 GB', '3h', 'user-trim', 1]), observed
            elif key == 'DIFFSPLICE':
                assert observed == [[1, '4 GB', '2h', 'unset', 1], [1, '8 GB', '4h', 'unset', 2]], observed
            elif key == 'RMATS_POST_PAIRED_PROBE':
                assert observed == [[12, '16 GB', '16h', 'unset', 1], [12, '32 GB', '1d 8h', 'unset', 2]], observed
            else:
                assert len(observed) == 1, (key, observed)
        proof[case] = rows
    cache = memory_only_resume(root/'memory-only-resume', command)
    receipt = {'status': 'passed', 'declared_processes': len(catalog['processes']),
        'scope': 'Native request rendering, aliases, retry and user overrides; no scientific or Slurm jobs',
        'local_capacity_simulated': True, 'allocations': proof, 'memory_only_resume': cache}
    (root / 'acceptance.receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'status': 'passed', 'declared_processes': len(catalog['processes']), 'cases': list(proof)}))


if __name__ == '__main__':
    main()
