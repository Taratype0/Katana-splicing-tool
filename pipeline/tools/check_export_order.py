#!/usr/bin/env python3
"""Check actual export ordering under reversed arrival and native resume."""
import argparse
import csv
import hashlib
import json
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--nextflow', default='nextflow')
    parser.add_argument('--nextflow-command-json')
    args = parser.parse_args()
    runtime = json.loads(args.nextflow_command_json) if args.nextflow_command_json else [args.nextflow]
    if not isinstance(runtime, list) or not runtime or not all(isinstance(x, str) for x in runtime):
        raise ValueError('Nextflow command must be a nonempty string array')
    root = args.outdir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    source = Path(__file__).resolve().parents[1] / 'workflows/katana_core.nf'
    text = source.read_text(encoding='utf-8')
    helper = text[text.index('def exportArtifactInputs('):text.index('\ndef qcArtifacts(')]
    inputs = root / 'inputs'
    inputs.mkdir()
    for name in ('alpha', 'beta', 'gamma'):
        (inputs / (name + '.txt')).write_text(name + '\n')
    before = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs.iterdir()}
    script = r'''
import groovy.json.JsonOutput
process EXPORT_ORDER {
    input:
    tuple val(metadata), path(files, stageAs: 'artifacts/??????/*')
    output:
    path 'ordered.json'
    script:
    def rows = (0..<metadata.size()).collect { i -> [key: metadata[i].relative_path, file: files[i].toString()] }
    def encoded = JsonOutput.toJson(rows).bytes.encodeBase64().toString()
    """
    python3 - '${encoded}' <<'PY'
import base64,json,sys,pathlib
rows=json.loads(base64.b64decode(sys.argv[1]))
for row in rows:
    row['text']=pathlib.Path(row.pop('file')).read_text().strip()
pathlib.Path('ordered.json').write_text(json.dumps(rows))
PY
    """
}
workflow {
    rows = ['alpha','beta','gamma'].collect { key -> [[relative_path: "${key}.txt"], file("${params.inputs}/${key}.txt")] }
    if (params.reverse) rows = rows.reverse()
    EXPORT_ORDER(Channel.value(exportArtifactInputs(rows)))
}
'''
    (root / 'main.nf').write_text(helper + script)
    (root / 'nextflow.config').write_text('params.reverse=false\nprocess.executor="local"\n')
    evidence = []
    for name, extra in [('first', []), ('reverse', ['--reverse', 'true', '-resume'])]:
        trace = root / (name + '.tsv')
        command = runtime + ['run', str(root / 'main.nf'), '--inputs', str(inputs),
                   '-ansi-log', 'false', '-with-trace', str(trace), *extra]
        result = subprocess.run(command, cwd=root, text=True, capture_output=True, timeout=120)
        (root / (name + '.log')).write_text(result.stdout + result.stderr)
        assert result.returncode == 0, result.stdout + result.stderr
        with trace.open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        assert len(rows) == 1 and rows[0]['status'] == ('COMPLETED' if name == 'first' else 'CACHED'), rows
        evidence.append({'run': name, 'status': rows[0]['status'], 'hash': rows[0]['hash']})
    assert evidence[0]['hash'] == evidence[1]['hash']
    output = list((root / 'work').glob('*/*/ordered.json'))
    assert len(output) == 1
    assert json.loads(output[0].read_text()) == [{'key': name + '.txt', 'text': name} for name in ('alpha','beta','gamma')]
    assert before == {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs.iterdir()}
    receipt = {'status': 'passed', 'native_runs': evidence, 'metadata_file_pairs_preserved': True,
               'inputs_unchanged': True, 'production_helper_sha256': hashlib.sha256(helper.encode()).hexdigest()}
    (root / 'acceptance.json').write_text(json.dumps(receipt, indent=2))
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
