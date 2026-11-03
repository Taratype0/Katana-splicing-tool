#!/usr/bin/env python3
"""Compare grouped and per-BAM rMATS prep on generated junction reads (Linux)."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import pysam

IMAGE = 'quay.io/biocontainers/mulled-v2-8ea76ff0a6a4c7e5c818fd4281abf918f92eeeae:121e48ab4817ec619c157a346458efca1ccf3c0a-0'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, required=True)
    args = parser.parse_args()
    root = args.outdir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    gtf = []
    for transcript, exons in [('included',[(101,200),(301,400),(501,600)]), ('skipped',[(101,200),(501,600)])]:
        for start, end in exons:
            gtf.append(f'chr1\tfixture\texon\t{start}\t{end}\t.\t+\t.\tgene_id "gene1"; transcript_id "{transcript}"; gene_name "GENE1";')
    (root/'ref.gtf').write_text('\n'.join(gtf)+'\n')
    names=[]
    for group in ('A','B'):
        for replicate in range(3):
            name=f'{group}{replicate+1}.bam'; names.append(name)
            reads=[]
            for junction, start, gap, count in [('left',150,100,30 if group=='A' else 10),
                                               ('right',350,100,30 if group=='A' else 10),
                                               ('skip',150,300,10 if group=='A' else 30)]:
                for index in range(count+replicate):
                    read=pysam.AlignedSegment()
                    read.query_name=f'{group}{replicate}_{junction}_{index}'
                    read.query_sequence='A'*100
                    read.query_qualities=pysam.qualitystring_to_array('I'*100)
                    read.flag=0; read.reference_id=0; read.reference_start=start
                    read.mapping_quality=60; read.cigar=[(0,50),(3,gap),(0,50)]
                    read.set_tag('NH',1); reads.append(read)
            with pysam.AlignmentFile(str(root/name),'wb',header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':1000}]}) as handle:
                for read in sorted(reads,key=lambda r:r.reference_start): handle.write(read)
            pysam.index(str(root/name))
    (root/'A.txt').write_text(','.join(names[:3])+'\n')
    (root/'B.txt').write_text(','.join(names[3:])+'\n')
    base=['docker','run','--rm','-v',f'{root}:/data','-w','/data',IMAGE,'rmats.py']
    common=['--gtf','ref.gtf','--readLength','100','-t','single','--libType','fr-unstranded',
            '--variable-read-length','--allow-clipping','--nthread','1','--tstat','1','--cstat','0.0001']
    def execute(label, arguments):
        with (root/f'{label}.log').open('w') as stream:
            subprocess.run(base+common+arguments,stdout=stream,stderr=subprocess.STDOUT,check=True,timeout=180)
    execute('grouped-prep',['--task','prep','--b1','A.txt','--b2','B.txt','--tmp','grouped_tmp','--od','grouped_prep'])
    (root/'split_tmp').mkdir()
    for name in names:
        sample=name.removesuffix('.bam')
        (root/f'{sample}.txt').write_text(name+'\n')
        execute(sample+'-prep',['--task','prep','--b1',f'{sample}.txt','--tmp',sample+'_tmp','--od',sample+'_prep'])
        files=list((root/(sample+'_tmp')).glob('*.rmats'))
        assert len(files)==1, files
        (root/'split_tmp'/f'{sample}.rmats').write_bytes(files[0].read_bytes())
    for mode in ('grouped','split'):
        execute(mode+'-post',['--task','post','--b1','A.txt','--b2','B.txt',
            '--tmp',mode+'_tmp','--od',mode+'_post','--statoff'])
    compared={}
    def rows(path):
        with path.open() as stream:
            records=list(csv.DictReader(stream,delimiter='\t'))
        return sorted(tuple(sorted((key,value) for key,value in row.items() if key!='ID')) for row in records)
    for old in sorted((root/'grouped_post').glob('*.MATS.*.txt')):
        new=root/'split_post'/old.name
        assert new.is_file() and rows(old)==rows(new), old.name
        compared[old.name]=len(rows(old))
    assert len(compared)==10 and compared['SE.MATS.JC.txt']>0, compared
    (root/'acceptance.json').write_text(json.dumps({'status':'passed','scope':'generated reads; grouped versus per-BAM prep count equivalence',
        'real_study_run':False,'image':IMAGE,'files':compared},indent=2))


if __name__ == '__main__': main()
