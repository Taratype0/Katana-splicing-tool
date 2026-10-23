"""Real packaged DEXSeq CLI software fixture; not a study analysis."""
import hashlib
import json
import platform
import subprocess
import tempfile
from pathlib import Path

import HTSeq
import pkg_resources
import pysam

assert HTSeq.__version__ == '2.0.2'
assert pkg_resources.get_distribution('htseq').version == '2.0.2'
expected_source = {
    '/usr/local/bin/dexseq_prepare_annotation.py': '8faec3a78d15e8abd8465adb6ee4bd45d5577aa9e650bc2b09c8ea8107d952e6',
    '/usr/local/bin/dexseq_count.py': 'd0cb14bfddc2bc7cba95e0adcca56df363f073a82c19cff7937ecb854543e36e',
    '/usr/share/doc/kanata-dexseq-tools/LICENSE.GPL-3.0': '3972dc9744f6499f0f9b2dbf76696f2ae7ad8af9b23dde66d6af86c9dfb36986',
}
for name, wanted in expected_source.items():
    assert hashlib.sha256(Path(name).read_bytes()).hexdigest() == wanted, name

with tempfile.TemporaryDirectory(prefix='kanata-dexseq-software-fixture-') as temporary:
    root = Path(temporary)
    gtf = root/'annotation.gtf'
    gff = root/'flattened.gff'
    gtf.write_text(
        'X\tsoftware_fixture\texon\t101\t150\t.\t+\t.\tgene_id "fixture_gene"; transcript_id "fixture_transcript";\n'
        'X\tsoftware_fixture\texon\t201\t250\t.\t+\t.\tgene_id "fixture_gene"; transcript_id "fixture_transcript";\n'
    )
    header = {'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'X','LN':1000}]}
    observed = {}
    commands = []
    prepare = ['dexseq_prepare_annotation.py',str(gtf),str(gff)]
    subprocess.run(prepare,check=True)
    commands.append(prepare)
    for paired in (False,True):
        mode = 'paired_end' if paired else 'single_end'
        bam = root/f'{mode}.bam'
        with pysam.AlignmentFile(str(bam),'wb',header=header) as output:
            for index,start in enumerate((110,210)):
                read = pysam.AlignedSegment()
                read.query_name = 'fixture_pair' if paired else f'fixture_read{index}'
                read.query_sequence = 'A'*20
                read.query_qualities = pysam.qualitystring_to_array('I'*20)
                read.flag = (99 if index==0 else 147) if paired else 0
                read.reference_id = 0
                read.reference_start = start
                read.mapping_quality = 60
                read.cigar = [(0,20)]
                read.set_tag('NH',1)
                if paired:
                    read.next_reference_id = 0
                    read.next_reference_start = 210 if index==0 else 110
                    read.template_length = 120 if index==0 else -120
                output.write(read)
        pysam.index(str(bam))
        counts = root/f'{mode}.clean.count.txt'
        command = ['dexseq_count.py',str(gff)]
        if paired:
            command.extend(['-p','yes'])
        command.extend(['-f','bam',str(bam),'-r','pos',str(counts),'-a','10','-s','no'])
        subprocess.run(command,check=True)
        commands.append(command)
        values = dict((name,int(number)) for name,number in (line.split() for line in counts.read_text().splitlines()))
        assert {name:number for name,number in values.items() if not name.startswith('_')} == {'fixture_gene:001':1,'fixture_gene:002':1}
        assert all(number==0 for name,number in values.items() if name.startswith('_'))
        observed[mode] = values
    print(json.dumps({'status':'passed_real_dexseq_packaged_software_fixture',
                      'scientific_study_analysis':False,'full_pipeline_acceptance':'not_established',
                      'python':platform.python_version(),'architecture':platform.machine(),
                      'htseq':HTSeq.__version__,'pysam':pysam.__version__,
                      'source_and_license_sha256':expected_source,'counts':observed,
                      'command_count':len(commands)},indent=2))
