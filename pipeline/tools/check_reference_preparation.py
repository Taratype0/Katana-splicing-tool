#!/usr/bin/env python3
"""Exercise real reference preparation without modifying same-named inputs."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import subprocess
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--nextflow-command-json', default='["nextflow"]')
    args = parser.parse_args()
    root = args.outdir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    runtime = json.loads(args.nextflow_command_json)
    if not isinstance(runtime, list) or not runtime or not all(isinstance(x, str) for x in runtime):
        raise ValueError('Nextflow command must be a nonempty string array')
    module = Path(__file__).resolve().parents[1] / 'modules/local/prepare_references'
    escaped = module.as_posix().replace('\\', '\\\\').replace("'", "\\'")
    (root/'main.nf').write_text(
        "nextflow.enable.dsl = 2\ninclude { PREPARE_REFERENCES } from '" + escaped + "'\n"
        "workflow { PREPARE_REFERENCES(Channel.value(file(params.fasta, checkIfExists:true)), "
        "Channel.value(file(params.gtf, checkIfExists:true))) }\n")
    (root/'nextflow.config').write_text(
        "docker.enabled=true\ndocker.runOptions='--network none'\n"
        "process { cpus=1; memory='1 GB'; time='10m'; errorStrategy='terminate' }\n")
    fasta = b'>chr1\nACGTACGTACGTACGTACGT\n'
    gtf = b'chr1\tfixture\texon\t1\t20\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    cases = []
    for label, compressed, empty in [('same-basename', False, False), ('compressed', True, False), ('empty-output', False, True)]:
        case = root/label
        case.mkdir()
        inputs = case/'inputs'
        inputs.mkdir()
        suffix = '.gz' if compressed else ''
        fa, annotation = inputs/('genome.fasta'+suffix), inputs/('annotation.gtf'+suffix)
        encode = gzip.compress if compressed else lambda x: x
        fa.write_bytes(encode(fasta))
        annotation.write_bytes(encode(b'# no transcripts\n' if empty else gtf))
        before = {p.name: (digest(p), p.stat().st_mtime_ns) for p in (fa, annotation)}
        command = runtime + ['run', str(root/'main.nf'), '-c', str(root/'nextflow.config'),
                             '--fasta', str(fa), '--gtf', str(annotation), '-ansi-log', 'false',
                             '-work-dir', str(case/'work')]
        result = subprocess.run(command, cwd=case, text=True, capture_output=True, timeout=600)
        (case/'stdout.txt').write_text(result.stdout)
        (case/'stderr.txt').write_text(result.stderr)
        after = {p.name: (digest(p), p.stat().st_mtime_ns) for p in (fa, annotation)}
        assert before == after, f'{label}: a source file changed'
        assert (result.returncode == 0) != empty, (label, result.stdout, result.stderr)
        hashes = list((case/'work').rglob('reference.sha256'))
        if empty:
            assert not hashes, 'Failed reference preparation published a success hash list'
        else:
            assert len(hashes) == 1
            work = hashes[0].parent
            assert (work/'genome.fasta').read_bytes() == fasta
            assert (work/'annotation.gtf').read_bytes() == gtf
            assert 'ACGTACGTACGTACGTACGT' in (work/'transcripts.fasta').read_text()
            assert not (work/'annotation.gtf').is_symlink()
            assert not (work/'genome.fasta').is_symlink()
        cases.append({'case':label, 'exit_code':result.returncode, 'inputs_unchanged':True,
                      'source_sha256':{k:v[0] for k,v in before.items()}})
    receipt = {'schema':'katana.reference_preparation_check/v1','status':'passed',
               'actual_nextflow_and_gffread_container':True, 'cases':cases,
               'module_sha256':digest(module.with_suffix('.nf'))}
    (root/'acceptance.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()
