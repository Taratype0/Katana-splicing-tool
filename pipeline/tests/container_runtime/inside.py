"""Exercise installed CLI/runtime probes inside one disconnected image task."""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys

from katana_next.flow.bundled_runtime import mapp_input_defaults
from katana_next.flow.model import seal_plan, write_json, file_sha256


def command(argv, directory, name, expected=0):
    result = subprocess.run(argv, cwd=directory, capture_output=True, text=True, check=False)
    (directory / (name + '.stdout.txt')).write_text(result.stdout, encoding='utf-8')
    (directory / (name + '.stderr.txt')).write_text(result.stderr, encoding='utf-8')
    if result.returncode != expected:
        raise AssertionError(f'{name} exited {result.returncode}, expected {expected}: {result.stderr[-2000:]}')
    return result


def main():
    kind = sys.argv[1]
    assert kind in {'network', 'mapp', 'motif'}
    assert sys.prefix == '/opt/conda', sys.prefix
    assert os.environ.get('KATANA_CONTAINER') == '1'
    assert os.getuid() != 0, 'Acceptance must exercise the non-root task user'
    out = Path(kind).resolve()
    out.mkdir(exist_ok=False)
    command(['katana-flow', '--help'], out, 'katana-flow-help')
    project = write_json(out / 'project.json', {'assets': [{'kind': 'fimo_runtime', 'path': 'fimo.json'}]})
    write_json(out / 'fimo.json', {'python': 'python'})
    adapter = {'network': 'network.fit', 'mapp': 'mapp.maei'}.get(kind)
    task = {'id': kind, 'workflow_id': 'flow.' + adapter if adapter else 'prioritize.rbp_mechanism',
            'execution_group': kind, 'parameters': {}, 'dependencies': [],
            'inputs': mapp_input_defaults({}) if kind == 'mapp' else {}}
    if adapter:
        task['adapter'] = adapter
    body = {'project_manifest': str(project), 'project_manifest_sha256': file_sha256(project), 'tasks': [task]}
    plan = write_json(out / 'plan.json', seal_plan(body))
    runtime = out / 'runtime.json'
    command(['katana-flow', 'runtime-check', '--plan', str(plan), '--project', str(project),
             '--task-id', kind, '--runtime-label', 'actual-nextflow-docker-image', '--out', str(runtime)], out, 'runtime-check')
    probe = json.loads(runtime.read_text(encoding='utf-8'))
    assert probe['status'] == 'ready' and probe['scientific_computation'] is False, probe
    if kind == 'mapp':
        assert probe.get('native_source', {}).get('native_source_provider') == 'katana.bundled_mapp/v1'

    # This existing executable is in mounted task storage, outside /opt/conda.
    # A rejection based only on a nonexistent path would not prove isolation.
    external = out / 'mounted-host-runtime' / 'python'
    external.parent.mkdir()
    external.write_text('#!/bin/sh\necho SHOULD_NOT_RUN > host-interpreter-executed\nexit 0\n', encoding='utf-8')
    external.chmod(0o755)
    bad_task = {**task, 'id': kind + '-reject-host', 'parameters': {'python': str(external)}}
    bad_plan = write_json(out / 'reject-host.plan.json', seal_plan({**body, 'tasks': [bad_task]}))
    bad_report = out / 'reject-host.runtime.json'
    command(['katana-flow', 'runtime-check', '--plan', str(bad_plan), '--project', str(project),
             '--task-id', bad_task['id'], '--out', str(bad_report)], out, 'reject-host', expected=2)
    rejected = json.loads(bad_report.read_text(encoding='utf-8'))
    assert rejected['status'] == 'failed' and 'external host path' in ' '.join(rejected['errors']), rejected
    assert not (out / 'host-interpreter-executed').exists()
    probe.update(probe_kind=kind, installed_prefix=sys.prefix, task_uid=os.getuid(), task_gid=os.getgid(),
                 host_runtime_rejected=True, existing_mounted_runtime_tested=True,
                 host_runtime_error=rejected['errors'], host_runtime_executed=False,
                 native_cli_runtime_check=True, container_network_disabled_by_docker=True,
                 acceptance_scope='Installed image runtime and Nextflow transport; no scientific model or study analysis')
    write_json(out / 'acceptance.json', probe)
    print(json.dumps({'probe_kind': kind, 'status': 'ready', 'host_runtime_rejected': True}))


if __name__ == '__main__':
    main()
