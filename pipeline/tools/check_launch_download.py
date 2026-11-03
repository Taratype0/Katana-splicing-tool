#!/usr/bin/env python3
"""Native launch-host download and offline reuse with a forbidden Slurm submitter."""
import argparse
import csv
import hashlib
import http.server
import json
import os
from pathlib import Path
import subprocess
import threading


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nextflow-command-file', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    args = parser.parse_args()
    runtime = json.loads(args.nextflow_command_file.read_text())
    pipeline = Path(__file__).resolve().parents[1]
    root = args.outdir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    quote = lambda value: "'" + str(value).replace('\\', '/').replace("'", "\\'") + "'"
    content = b'generated-pathway\tfixture\tGeneA\tGeneB\n'
    checksum = hashlib.sha256(content).hexdigest()
    requests = []

    class Handler(http.server.BaseHTTPRequestHandler):
        def do_GET(self):
            requests.append(self.path)
            self.send_response(200)
            self.end_headers()
            self.wfile.write(content)

        def log_message(self, *args):
            pass

    server = http.server.ThreadingHTTPServer(('127.0.0.1', 0), Handler)
    worker = threading.Thread(target=server.serve_forever, daemon=True)
    worker.start()
    cache = root / 'resource-cache'
    url = f'http://127.0.0.1:{server.server_port}/fixture.gmt'
    item = dict(asset_id='fixture-genes', parameter='gene_sets', filename='fixture.gmt',
                url=url, metadata={'sha256': checksum})
    request = dict(mode='complete', species='Mus_musculus', reference_build='GRCm39',
                   resource_names=[], resource_entries=[], resource_catalog=[item],
                   resource_cache=str(cache), skip_enrichment=False)
    for module in ('network', 'rbp_identity', 'motiftrace', 'motif_binding', 'mapp'):
        request['skip_' + module] = True
    request_path = root / 'request.json'
    script = "include { KATANA_FLOW_CONFIGURE } from " + quote(pipeline / 'modules/local/flow_control') + '\n'
    script += '''process SITE_EXECUTOR_PROBE {
    output:
    path 'probe.txt'
    script:
    "echo should-not-execute > probe.txt"
}
workflow {
    request = new groovy.json.JsonSlurper().parse(file(params.request).toFile())
    KATANA_FLOW_CONFIGURE(request, [], file(params.generator))
    KATANA_FLOW_CONFIGURE.out.configuration.view { "configuration=${it}" }
    if (params.probe) SITE_EXECUTOR_PROBE()
}
'''
    (root / 'main.nf').write_text(script)
    site = root / 'site.config'
    site.write_text('process.executor="slurm"\nprocess.cpus=1\nprocess.memory=128.MB\n'
                    'params.katana_container=null\nparams.probe=false\n')
    fake_bin = root / 'bin'
    fake_bin.mkdir()
    marker = root / 'sbatch-called.txt'
    sbatch = fake_bin / 'sbatch'
    sbatch.write_text('#!/bin/sh\nprintf "called\\n" >> "' + str(marker) + '"\nexit 44\n')
    sbatch.chmod(0o755)
    env = {**os.environ, 'PATH': str(fake_bin) + os.pathsep + os.environ['PATH'],
           'PYTHONPATH': str(pipeline.parent), 'NXF_OFFLINE': 'true'}
    evidence = {}

    def run(name, success=True, probe=False):
        request['study_name'] = name  # Force a new configure task, not a cache hit.
        request_path.write_text(json.dumps(request))
        trace = root / (name + '.tsv')
        result = subprocess.run([*runtime, 'run', str(root / 'main.nf'),
            '-c', str(pipeline / 'nextflow.config'), '-c', str(site),
            '--request', str(request_path), '--generator', str(pipeline / 'bin/katana_configure.py'),
            '--resource_cache', str(cache), '--digest_cache', str(root / 'digests'),
            '--outdir', str(root / 'output'), '--probe', str(probe).lower(),
            '-work-dir', str(root / 'work'), '-with-trace', str(trace), '-ansi-log', 'false'],
            cwd=root, env=env, text=True, capture_output=True, timeout=150)
        log = result.stdout + result.stderr
        (root / (name + '.log')).write_text(log)
        assert (result.returncode == 0) == success, log
        rows = list(csv.DictReader(trace.open(), delimiter='\t'))
        if success:
            assert len(rows) == 1 and rows[0]['status'] == 'COMPLETED', rows
            assert rows[0]['name'].endswith('KATANA_FLOW_CONFIGURE'), rows
            assert not marker.exists(), 'Download was submitted to Slurm'
        evidence[name] = dict(exit_code=result.returncode, tasks=rows)
        return log

    try:
        run('online')
        assert requests == ['/fixture.gmt'], requests
        saved = cache / checksum / 'fixture.gmt'
        assert saved.read_bytes() == content
        server.shutdown()
        server.server_close()
        worker.join(timeout=5)
        run('offline-cache')
        assert requests == ['/fixture.gmt'], requests
        request['resource_catalog'] = [{**item, 'filename': 'unavailable.gmt', 'url': url + '-missing'}]
        log = run('offline-missing', success=False)
        assert 'Unable to download resource unavailable.gmt' in log, log
        assert 'Nextflow launch host' in log, log
        assert saved.read_bytes() == content
        assert not marker.exists()
        request['resource_catalog'] = [item]
        run('site-executor', success=False, probe=True)
        calls = marker.read_text().splitlines()
        assert calls and all(line == 'called' for line in calls), 'Scientific task did not use site executor'
        receipt = dict(status='passed', launch_download_local=True, site_executor_preserved=True,
                       offline_cached_bytes_reused=True, missing_resource_error_identifies_asset=True,
                       cases=evidence, real_study_run=False)
        (root / 'acceptance.json').write_text(json.dumps(receipt, indent=2))
        print(json.dumps(receipt, indent=2))
    finally:
        server.shutdown()
        server.server_close()


if __name__ == '__main__':
    main()
