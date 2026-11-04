"""Installed-package integration on generated data only; no study inputs.

Runs the public flow task executor, the standalone CoSplice command, and the
approved report builder. Native R and EventGraph calculations are not mocked.
"""
from pathlib import Path
import argparse
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd


def fixture(root, paired=False):
    from katana_next.full_fixture import create_full_fixture
    from katana_next.infrastructure.manifest_store import ManifestStore
    from katana_next.domain.models import InputAsset
    path = create_full_fixture(root)
    _, manifest = ManifestStore().load(path)
    if paired:
        names = {s.sample_id for c in manifest.contrasts[:1] for arm in (c.numerator,c.denominator) for s in manifest.samples if s.sample_id in arm.sample_ids[:3]}
        samples = tuple(s.model_copy(update={'biological_sample_id':f'unit{i%3}'}) for i,s in enumerate(s for s in manifest.samples if s.sample_id in names))
        c = manifest.contrasts[0]
        contrast = c.model_copy(update={k:getattr(c,k).model_copy(update={'sample_ids':getattr(c,k).sample_ids[:3]}) for k in ('numerator','denominator')})
        manifest = manifest.model_copy(update={'samples':samples,'contrasts':(contrast,), 'metadata':{**manifest.metadata,'comparison_designs':{contrast.contrast_id:{'paired':True}}}})
    rng=np.random.default_rng(20261102)
    ids=[s.sample_id for s in manifest.samples]
    n=len(ids)
    # Distinct, balanced latent profiles; generated values are not study evidence.
    latent=np.vstack([np.sin(np.arange(n)*2*np.pi/n),np.cos(np.arange(n)*2*np.pi/n)])
    if paired:
        latent=np.array([[-.9,-.7,-.5,.5,.7,.9],[-.5,-.9,-.7,.9,.5,.7]])
    values=np.clip(np.vstack([.5+.25*latent[i%2]+rng.normal(0,.018,n) for i in range(240)]),.01,.99)
    keys=[f'synthetic:{i}' for i in range(len(values))]
    inputs=path.parent/'inputs'
    # Reader-only fixtures need no splice geometry. The executable fixture
    # supplies valid five-type, both-strand forms before sequence inference.
    for asset in manifest.assets:
        if asset.kind != 'rmats_result':
            continue
        frame=pd.read_csv(asset.path,sep='\t')
        kind=asset.metadata['event_type']
        if kind == 'SE':
            # Supply 24 independent foreground/background gene pairs so the
            # native family actually estimates P/FDR, rather than testing only
            # the valid-but-unestimable small-family reporting branch.
            template=frame.iloc[0].to_dict(); rows=[]
            for i in range(24):
                offset=1000+2500*i
                for background in (False,True):
                    row=dict(template)
                    row.update(ID=2*i+int(background),GeneID=f'synthetic_{"bg" if background else "fg"}{i}',
                        geneSymbol=f'SYN_{"BG" if background else "FG"}{i}',chr='chr1',strand='+',
                        upstreamES=offset,upstreamEE=offset+40,downstreamES=offset+300,downstreamEE=offset+340,
                        exonStart_0base=offset+100+int(background),exonEnd=offset+135+int(background),
                        FDR=.8 if background else .001,IncLevelDifference=0. if background else .3,
                        IncLevel1=','.join(['.35' if background else '.65']*5),IncLevel2=','.join(['.35']*5),
                        IJC_SAMPLE_1=','.join(['35' if background else '65']*5),SJC_SAMPLE_1=','.join(['65' if background else '35']*5),
                        IJC_SAMPLE_2=','.join(['35']*5),SJC_SAMPLE_2=','.join(['65']*5))
                    rows.append(row)
            frame=pd.DataFrame(rows)
        if kind == 'RI':
            frame['riExonStart_0base']=frame['upstreamES']
            frame['riExonEnd']=frame['downstreamEE']
        elif kind in {'A3SS','A5SS'}:
            suffix=frame['strand'].eq('+') if kind=='A5SS' else frame['strand'].eq('-')
            frame.loc[suffix,'shortES']=frame.loc[suffix,'longExonStart_0base']
            frame.loc[suffix,'shortEE']=frame.loc[suffix,'longExonEnd']-45
            frame.loc[suffix,'flankingES']=frame.loc[suffix,'longExonEnd']+170
            frame.loc[suffix,'flankingEE']=frame.loc[suffix,'longExonEnd']+260
        frame.to_csv(asset.path,sep='\t',index=False)
    # The general reader fixture has only chr1; sequence inference requires a
    # reference and index covering every declared synthetic event chromosome.
    fasta=inputs/'fixture.genome.fa'
    sequence=fasta.read_text().splitlines()[1]
    chromosomes=sorted({str(chrom) for asset in manifest.assets if asset.kind=='rmats_result'
                        for chrom in pd.read_csv(asset.path,sep='\t')['chr'].unique()})
    blocks=[]; index=[]; offset=0
    for chrom in chromosomes:
        header='>'+chrom+'\n'; offset+=len(header)
        index.append(f'{chrom}\t{len(sequence)}\t{offset}\t{len(sequence)}\t{len(sequence)+1}\n')
        blocks.append(header+sequence+'\n'); offset+=len(sequence)+1
    fasta.write_text(''.join(blocks),encoding='ascii',newline='\n')
    Path(str(fasta)+'.fai').write_text(''.join(index),encoding='ascii',newline='\n')
    psi=inputs/'synthetic-psi.tsv'; events=inputs/'synthetic-events.tsv'
    pd.DataFrame(values,index=keys,columns=ids).to_csv(psi,sep='\t',index_label='event_key')
    pd.DataFrame({'event_key':keys,'gene_id':[f'ENSMUSG{i+1:011d}' for i in range(len(keys))], 'event_type':'SE'}).to_csv(events,sep='\t',index=False)
    # Generated DEXSeq-format results complete the four-layer input fixture.
    exon_assets=[]
    for c in manifest.contrasts:
        exon=inputs/(c.contrast_id+'.dexseq.tsv')
        pd.DataFrame({'groupID':[f'ENSMUSG{i:011d}.1' for i in range(1,41)],
                      'featureID':[f'E{i:03}' for i in range(1,41)],
                      'padj':[.005 if i<10 else .8 for i in range(40)],
                      'log2fold_'+c.denominator.label+'_'+c.numerator.label:[.6 if i<10 else .01 for i in range(40)]}).to_csv(exon,sep='\t',index=False)
        exon_assets.append(InputAsset(asset_id='synthetic-exons-'+c.contrast_id,kind='dexseq_result',path=exon,metadata={'contrast_id':c.contrast_id,'native_fold_change_direction':'control_over_treatment_native_header_authoritative'}))
    (inputs/'fixture.pathways.gmt').write_text(''.join(f'Synthetic_set{k}\tgenerated\t'+'\t'.join([f'GENE{i}' for i in range(k*20+1,k*20+21)]+[f'ENSMUSG{i:011d}' for i in range(k*20+1,k*20+21)])+'\n' for k in range(12)),encoding='utf8')
    meme=inputs/'synthetic-models.meme'
    meme.write_text('MEME version 4\nALPHABET= ACGU\n' + ''.join(f'MOTIF {name} Synthetic_{name}\nletter-probability matrix: alength= 4 w= 2\n{a}\n{b}\n' for name,a,b in [('AC','1 0 0 0','0 1 0 0'),('GU','0 0 1 0','0 0 0 1')]),encoding='utf8')
    Path(str(meme)+'.metadata.tsv').write_text('motif_id\tRBP\tspecies\tevidence\tsource\nAC\tR1\tMus_musculus\tuser_defined\tsynthetic\nGU\tR2\tMus_musculus\tuser_defined\tsynthetic\n',encoding='utf8')
    assets=tuple(a for a in manifest.assets if a.kind not in {'psi_matrix','event_annotation','cosplice_config'}) + tuple(exon_assets) + (
        InputAsset(asset_id='generated-psi',kind='psi_matrix',path=psi),InputAsset(asset_id='generated-events',kind='event_annotation',path=events))
    manifest=manifest.model_copy(update={'assets':assets,'metadata':{**manifest.metadata,'species':'Mus_musculus'}})
    ManifestStore().save(manifest,path)
    from qc_fixture import add_native_qc
    add_native_qc(path)
    return path,meme


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args(); root=args.out.resolve(); root.mkdir(parents=True,exist_ok=False)
    from katana_next.flow.planner import create_plan
    from katana_next.flow.model import write_json,read_json
    from katana_next.flow.runtime import execute_task,validate_receipt
    from katana_next.flow.collect import collect_results
    from katana_next.flow.cli import build_complete_report
    project,meme=fixture(root/'source')
    selected=['perturbation-vs-baseline','rescue-vs-baseline']
    settings={'design_factors':['condition'],'contrast_ids':selected,
        'analyses':{'coverage':False,'mapp':False,'motif_binding':False,'interaction':False,'sequencing_qc':True},
        'parameters':{'prioritize.motiftrace_consequences':{'custom_motif_path':str(meme),'cluster_permutations':99}}}
    plan=create_plan(project,settings); plan_path=write_json(root/'plan.json',plan)
    receipts={}
    for task in plan['tasks']:
        print('RUN',task['id'],flush=True)
        receipts[task['id']]=execute_task(plan_path,task['id'],root/'tasks'/task['id'],[receipts[d] for d in task['dependencies']])
    collection=collect_results(plan_path,list(receipts.values()),root/'collection')
    report=build_complete_report(collection.parent,root/'report')
    network=read_json(receipts['network-study'])
    browser=next(a for a in network['generated_artifacts'] if a.get('metadata',{}).get('report_role')=='whole_study_chapter')
    native=Path(browser['path']).parent
    assert read_json(native/'status.json')['status']=='completed'
    associations=pd.read_csv(native/'interpretations/default/associations.tsv',sep='\t')
    assert set(associations.loc[associations.kind.eq('contrast'),'contrast_id'])==set(selected)
    assert (native/'model/native_model.rds').stat().st_size>0
    assert read_json(native/'study/manifest.json')['samples']==15
    manifest=read_json(root/'report/report-manifest.json')
    assert manifest['presentation']['chapter_order']==[1,2,3,4,5,7,6]
    assert manifest['network_browsers'][0]['whole_study'] is True
    assert {p['contrast_id'] for p in manifest['panels'] if p['title'].startswith('DEXSeq | ')} == {'perturbation-vs-baseline','rescue-vs-baseline'}
    summary = next(t for t in manifest['complete_results'] if 'Resolved library strandedness' in t['name'])
    assert summary['rows'] == 15
    for analysis in ('DEG','DAS','DTU','DEXSeq'):
        assert {t['comparison_id'] for t in manifest['top10'] if t['analysis']==analysis} == {'perturbation-vs-baseline','rescue-vs-baseline'}

    assert (root/'report' / summary['source']).is_file()

    # Test actual public CLI, explicit three-pair design, and cache reuse.
    paired,_=fixture(root/'paired-source',paired=True)
    command=['katana-cosplice','--project',str(paired),'run']
    outcomes=[]
    for i in range(2):
        result=subprocess.run(command,check=True,capture_output=True,text=True)
        (root/f'paired-cli-{i}.json').write_text(result.stdout,encoding='utf8')
        outcomes.append(json.loads(result.stdout)['status'])
    assert outcomes==['complete','cache_hit'],outcomes
    from katana_cosplice.study_api import CoSpliceStudy
    with CoSpliceStudy(paired) as client:
        browser=next(a for a in client.artifacts() if a.metadata.get('report_role')=='whole_study_chapter')
        shared=read_json(browser.path.parent/'study/manifest.json')
        assert shared['samples']==6 and shared['biological_units']==3
    write_json(root/'acceptance.json',{'status':'passed','generated_data_only':True,'native_tasks':len(receipts),'paired_cli':outcomes,'paired_observations':6,'paired_units':3,'report':report})
    print('ACCEPTED',root,flush=True)


if __name__=='__main__':
    main()
