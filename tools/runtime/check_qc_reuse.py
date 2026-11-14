"""Generated QC logs -> real MultiQC -> installed Katana -> fixed report.

No raw reads, BAM scans, model fits, or production data are used.
"""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import zipfile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--multiqc', default='multiqc')
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--existing-multiqc', type=Path)
    parser.add_argument('--prepare-only', action='store_true')
    parser.add_argument('--nextflow', type=Path)
    args = parser.parse_args()
    out = args.outdir.resolve()
    out.mkdir(parents=True)
    raw = out/'native'
    raw.mkdir()
    for sample in ('S1', 'S2'):
        text = f'''##FastQC\t0.12.1
>>Basic Statistics\tpass
#Measure\tValue
Filename\t{sample}_R1.fastq.gz
File type\tConventional base calls
Encoding\tSanger / Illumina 1.9
Total Sequences\t1000
Sequences flagged as poor quality\t0
Sequence length\t150
%GC\t50
>>END_MODULE
>>Sequence Length Distribution\tpass
#Length\tCount
150\t1000
>>END_MODULE
'''
        with zipfile.ZipFile(raw/(sample+'_R1_fastqc.zip'), 'w') as archive:
            archive.writestr(sample+'_R1_fastqc/', '')
            archive.writestr(sample+'_R1_fastqc/fastqc_data.txt', text)
        directory = raw/sample/'raw_data_qualimapReport'
        directory.mkdir(parents=True)
        (directory/'coverage_profile_along_genes_(total).txt').write_text(
            '#Transcript position\tTranscript coverage profile\n' +
            ''.join(f'{i}\t{100+i}\n' for i in range(1,101)))
        (directory.parent/'rnaseq_qc_results.txt').write_text(
            f'bam file = {sample}.bam\ngff file = generated.gtf\nprotocol = non-strand-specific\n'
            "reads aligned = 1000\nexonic = 800\nintronic = 150\nintergenic = 50\n"
            "5' bias = 1.0\n3' bias = 1.2\n5'-3' bias = 0.83\n"
            "5'-3' bias region size = 100\n5'-3' bias number of top transcripts = 1000\n")
    if args.existing_multiqc:
        shutil.copytree(args.existing_multiqc, out/'multiqc')
    elif args.nextflow:
        import os
        module = Path(__file__).resolve().parents[2]/'pipeline/modules/local/multiqc.nf'
        nf = out/'check.nf'
        nf.write_text("include { KATANA_MULTIQC } from '"+module.as_posix()+"'\n"
                      "workflow { KATANA_MULTIQC(Channel.value(file(params.logs))) }\n")
        config = out/'nextflow.config'
        config.write_text("process { withName: KATANA_MULTIQC { publishDir = [path: '"+
                          (out/'multiqc').as_posix()+"', mode: 'copy'] } }\n")
        env = dict(os.environ, PATH=str(Path(args.multiqc).resolve().parent)+os.pathsep+os.environ['PATH'])
        command = [str(args.nextflow.resolve()), 'run', str(nf), '-c', str(config),
                   '--logs', str(raw), '-work-dir', str(out/'work')]
        subprocess.run(command+['-with-trace', str(out/'trace-first.tsv')], cwd=out, env=env, check=True)
        subprocess.run(command+['-resume', '-with-trace', str(out/'trace-resume.tsv')], cwd=out, env=env, check=True)
        assert '\tCACHED\t' in (out/'trace-resume.tsv').read_text()
    else:
        subprocess.run([args.multiqc, str(raw), '--outdir', str(out/'multiqc'),
                        '--filename', 'multiqc_report.html', '--data-format', 'json'], check=True)
    html = out/'multiqc/multiqc_report.html'
    data = out/'multiqc/multiqc_report_data/multiqc_data.json'
    if not data.exists() and (out/'multiqc/multiqc_data.json').is_file():
        data = out/'multiqc/multiqc_data.json'
    assert html.is_file() and data.is_file()
    saved = json.loads(data.read_text())['report_saved_raw_data']['qualimap_rnaseq_cov_hist']
    assert set(saved) == {'S1', 'S2'}, list(saved)
    if args.prepare_only:
        return
    from katana_next.flow.model import write_json
    from katana_next.flow.sequencing import execute_stage, readiness
    project = write_json(out/'project.json', {
        'schema_version':'2.0','project_id':'generated-qc-reuse','name':'Generated QC acceptance',
        'input_root':'.','output_root':'unused','annotation':{'genome_build':'generated'},
        'samples':[{'sample_id':s,'biological_sample_id':s,'condition':s} for s in ('S1','S2')],
        'contrasts':[], 'assets':[{'asset_id':str(i),'kind':'sequencing_qc','path':str(p),
            'metadata':{'qc_type':'multiqc'}} for i,p in enumerate((html,data))]})
    spec = {'adapter':'qc.sequencing', 'inputs':{'project':str(project)}}
    assert readiness(spec) == []
    result = execute_stage(spec, out/'qc')
    from katana_next.infrastructure.manifest_store import ManifestStore
    from katana_next.flow.native_coverage import profile_selection
    _, loaded = ManifestStore().load(project)
    selected, frames = profile_selection(loaded)
    assert set(selected) == {'S1', 'S2'}
    assert len(frames['1']) == 200
    assert frames['1'].coverage.min() == 101
    from katana_next.reports.engine import Builder
    from katana_next.reports.layout import render
    from katana_next.reports.engine import ASSET_ROOT
    report = out/'report'
    report.mkdir()
    shutil.copytree(ASSET_ROOT, report/'assets')
    builder = Builder({'project':json.loads(project.read_text()), 'artifacts':[],
                       'mode':'bulk','dummy':True}, report, 'stable_significant')
    index = json.loads((out/'qc/sequencing-qc-views.json').read_text())
    artifact = index['artifacts'][0]
    builder.add_views({'panels':[{'id':'native-multiqc', 'chapter':1, 'title':artifact['name'],
        'path':artifact['path'], 'metadata':artifact['metadata'], 'data_links':artifact['metadata']['data_links']}]})
    manifest = render(builder)
    assert (report/'REPORT.zip').is_file()
    text = (report/'REPORT.html').read_text()
    assert 'title="MultiQC sequencing quality"' in text
    assert manifest['panels'][0]['metadata']['native_multiqc']
    proof = {'status':'passed','scope':'generated logs only; no scientific data rerun',
             'native_multiqc':result['evidence'], 'report':str(report/'REPORT.html'),
             'original_html_bytes_preserved':(report/manifest['panels'][0]['image']).read_bytes()==html.read_bytes()}
    assert proof['original_html_bytes_preserved']
    write_json(out/'acceptance.json',proof)
    print(json.dumps(proof))


if __name__ == '__main__':
    main()
