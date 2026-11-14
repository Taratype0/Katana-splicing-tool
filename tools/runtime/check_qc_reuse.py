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
            archive.writestr(sample+'_R1_fastqc/fastqc_data.txt', text)
    subprocess.run([args.multiqc, str(raw), '--outdir', str(out/'multiqc'),
                    '--filename', 'multiqc_report.html', '--data-format', 'json'], check=True)
    html = out/'multiqc/multiqc_report.html'
    data = out/'multiqc/multiqc_data/multiqc_data.json'
    assert html.is_file() and data.is_file()
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
