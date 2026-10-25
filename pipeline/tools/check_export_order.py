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
    version_source = source.parent.parent / 'subworkflows/local/dexseq_deu.nf'
    version_text = version_source.read_text(encoding='utf-8')
    version_helper = version_text[version_text.index('def stableCountVersion('):version_text.index('\nworkflow DEXSEQ_DEU')]
    inputs = root / 'inputs'
    inputs.mkdir()
    for name in ('alpha', 'beta', 'gamma'):
        (inputs / (name + '.txt')).write_text(name + '\n')
    # Directory order deliberately disagrees with sample identity. The old
    # first-arrival policy chooses B1 first, A2 after reversal; both are wrong.
    for sample, directory in [('B1', '40-b1'), ('A2', '00-a2'), ('A1', '90-a1')]:
        producer = inputs / directory
        producer.mkdir()
        (producer / (sample + '.clean.count.txt')).write_text(sample + '\n')
        (producer / 'versions.yml').write_text('fixture-count-tool: 1\n')
    before = {p.relative_to(inputs).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in inputs.rglob('*') if p.is_file()}
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
    producers = [[id:'B1', directory:'40-b1'], [id:'A2', directory:'00-a2'], [id:'A1', directory:'90-a1']]
    counts = producers.collect { row -> tuple([id:row.id], file("${params.inputs}/${row.directory}/${row.id}.clean.count.txt")) }
    versions = producers.collect { row -> file("${params.inputs}/${row.directory}/versions.yml") }
    versions = versions.tail() + [versions.head()]
    if (params.reverse) rows = rows.reverse()
    if (params.reverse) {
        counts = counts.reverse()
        versions = versions.reverse()
    }
    selected = stableCountVersion(Channel.fromList(counts), Channel.fromList(versions))
    selected.view { "selected_version=${it}" }
    with_version = Channel.value([artifacts: rows]).combine(selected).map { container, version ->
        exportArtifactInputs(container.artifacts + [[[relative_path: 'provenance/versions/' + version.parent.name + '.yml'], version]])
    }
    EXPORT_ORDER(with_version)
}
'''
    (root / 'main.nf').write_text(helper + version_helper + script)
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
    assert json.loads(output[0].read_text()) == [
        *[{'key': name + '.txt', 'text': name} for name in ('alpha','beta','gamma')],
        {'key': 'provenance/versions/90-a1.yml', 'text': 'fixture-count-tool: 1'}]
    assert before == {p.relative_to(inputs).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in inputs.rglob('*') if p.is_file()}
    receipt = {'status': 'passed', 'native_runs': evidence, 'metadata_file_pairs_preserved': True,
               'inputs_unchanged': True, 'version_representative': 'A1',
               'version_arrival_and_work_directory_order_ignored': True,
               'production_helper_sha256': hashlib.sha256(helper.encode()).hexdigest(),
               'production_version_helper_sha256': hashlib.sha256(version_helper.encode()).hexdigest()}
    (root / 'acceptance.json').write_text(json.dumps(receipt, indent=2))
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
