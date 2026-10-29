"""Execute both production import modules from the installed wheel (no export).

Only generated quantifications are used. The installed scripts are mounted
read-only into the same pinned upstream R image used by the real workflow.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess

from katana_next.flow.managed import _pipeline_source


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--nextflow-command-json', required=True)
    args = parser.parse_args()
    root = args.outdir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    pipeline = _pipeline_source().parent
    assert pipeline.name == '_pipeline', 'This check requires an installed wheel, not a checkout'
    scripts = ['katana_differential.R', 'tximport.R']
    before = {name: hashlib.sha256((pipeline / 'bin' / name).read_bytes()).hexdigest() for name in scripts}
    assert all(os.access(pipeline / 'bin' / name, os.X_OK) for name in scripts)
    # Nextflow discovers bin from its project directory. Use a read-only mounted
    # link to the installed resources, never export/chmod a replacement copy.
    (root / 'bin').symlink_to(pipeline / 'bin', target_is_directory=True)
    fixture = root / 'inputs'
    fixture.mkdir()
    for sample, counts in [('S1', [100, 80, 60, 40]), ('S2', [80, 60, 120, 20])]:
        quant = fixture / sample
        (quant / 'aux_info').mkdir(parents=True)
        with (quant / 'quant.sf').open('w') as handle:
            handle.write('Name\tLength\tEffectiveLength\tTPM\tNumReads\n')
            for i, count in enumerate(counts, 1):
                handle.write(f'T{i}\t1000\t800\t{count/sum(counts)*1e6}\t{count}\n')
        (quant / 'aux_info/meta_info.json').write_text(json.dumps({
            'salmon_version':'1.10.3', 'num_bootstraps':0, 'num_gibbs_samples':0,
            'num_valid_targets':4, 'num_processed':sum(counts), 'num_mapped':sum(counts)}))
    (fixture / 'samples.csv').write_text('sample,condition,subject\nS1,control,mouse1\nS2,treated,mouse1\n')
    (fixture / 'tx2gene.tsv').write_text(''.join(f'T{i}\tG{(i+1)//2}\tGene{(i+1)//2}\n' for i in range(1,5)))
    (fixture / 'annotation.gtf').write_text(''.join(
        f'chr1\tfixture\texon\t{i*2000}\t{i*2000+999}\t.\t+\t.\tgene_id "G{(i+1)//2}"; transcript_id "T{i}"; gene_name "Gene{(i+1)//2}"; gene_type "protein_coding";\n'
        for i in range(1,5)))
    module_dir = (pipeline / 'modules/local').as_posix()
    (root / 'main.nf').write_text(f'''nextflow.enable.dsl=2
include {{ KATANA_TXIMPORT }} from '{module_dir}/katana_differential'
include {{ TXIMPORT }} from '{module_dir}/tximport'
workflow {{
    quants = Channel.value([file('{fixture.as_posix()}/S1'),file('{fixture.as_posix()}/S2')])
    KATANA_TXIMPORT([id:'fixture'], quants, file('{fixture.as_posix()}/samples.csv'), file('{fixture.as_posix()}/annotation.gtf'), '')
    TXIMPORT(quants, file('{fixture.as_posix()}/tx2gene.tsv'))
}}
''')
    # Both processes use their production pinned container declarations.
    (root / 'nextflow.config').write_text(f'''process.executor='local'
process.cpus=1
process.memory='2 GB'
docker.enabled=true
docker.runOptions='--user {os.getuid()}:{os.getgid()} -v {pipeline.as_posix()}:{pipeline.as_posix()}:ro'
params.outdir='{(root / 'outputs').as_posix()}'
process {{ withName: TXIMPORT {{ publishDir = [path: '{(root / 'outputs/TXIMPORT').as_posix()}', mode: 'copy'] }} }}
''')
    prefix = json.loads(args.nextflow_command_json)
    for attempt, expected in [('first','COMPLETED'), ('resume','CACHED')]:
        trace = root / f'{attempt}.trace.tsv'
        command = [*prefix,'run',str(root / 'main.nf'),'-ansi-log','false','-with-trace',str(trace),'-work-dir',str(root / 'work')]
        if attempt == 'resume':
            command.append('-resume')
        with (root / f'{attempt}.log').open('w') as handle:
            subprocess.run(command, cwd=root, stdout=handle, stderr=subprocess.STDOUT, check=True)
        with trace.open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        assert len(rows) == 2 and all(r['status'] == expected for r in rows), rows
    with (root / 'outputs/DEG/gene.estimated_counts.tsv').open() as handle:
        counts = list(csv.DictReader(handle, delimiter='\t'))
    assert len(counts) == 2 and float(counts[0]['S1']) == 180 and float(counts[0]['S2']) == 140, counts
    assert (root / 'outputs/TXIMPORT/salmon.merged.gene_counts.tsv').stat().st_size > 0
    assert before == {name: hashlib.sha256((pipeline / 'bin' / name).read_bytes()).hexdigest() for name in scripts}
    receipt = dict(status='passed', installed_pipeline=str(pipeline), installed_scripts=before,
        native_import_tasks=2, resumed_cached_tasks=2, generated_samples=2, study_data_used=False,
        exported_pipeline_used=False, production_modules=['KATANA_TXIMPORT','TXIMPORT'])
    (root / 'acceptance.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()
