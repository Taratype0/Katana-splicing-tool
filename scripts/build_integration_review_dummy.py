"""Render a synthetic whole-report review with the selected installed renderer.

This is presentation acceptance, not a scientific pipeline run. No real data,
Slurm submission, database download, or inferred completion receipt is used.
Run with the candidate release's app directory on PYTHONPATH.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import shutil
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--outdir', type=Path, required=True)
    ap.add_argument('--network-demo', type=Path)
    ap.add_argument('--content-review', action='store_true', help='Four-comparison report-content review with a generated whole-study network and MAPP views.')
    args = ap.parse_args()
    root = args.outdir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    source = root / 'synthetic-inputs'; source.mkdir()
    rng = np.random.default_rng(20261102)
    blue, red = '#2863a5', '#c63f4a'
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
        'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
    views = {'panels':[], 'tables':[], 'top10':[]}
    warning = 'DUMMY DATA - illustrative values; not experimental results or computed inference.'
    conditions = ['Baseline','Treatment A','Treatment B','A + B','Recovery 1','Recovery 2']
    samples = [dict(sample_id=f'S{i+1:02}', biological_sample_id=f'BIO{i+1:02}',
                    condition=conditions[i//6]) for i in range(36)]
    contrast_specs = [('B_vs_A',1,0),('C_vs_A',2,0),('E_vs_A',4,0)]
    contrasts = [dict(contrast_id=cid, label=f'{conditions[n]} versus {conditions[d]}',
        numerator={'label':conditions[n],'sample_ids':[s['sample_id'] for s in samples if s['condition']==conditions[n]]},
        denominator={'label':conditions[d],'sample_ids':[s['sample_id'] for s in samples if s['condition']==conditions[d]]}) for cid,n,d in contrast_specs]
    if args.content_review:
        conditions = ['CONTROL_BASE','MUT_BASE','CONTROL_TREATED','MUT_TREATED']
        samples = [dict(sample_id=f'S{i+1:02}', biological_sample_id=f'{"BASE" if i<6 else "TREATED"}_PAIR{i%3+1}', condition=conditions[i//3]) for i in range(12)]
        contrast_specs = [('MUT_TREATED_vs_CONTROL_TREATED',3,2), ('MUT_BASE_vs_CONTROL_BASE',1,0),
                          ('CONTROL_TREATED_vs_CONTROL_BASE',2,0), ('MUT_TREATED_vs_MUT_BASE',3,1)]
        contrasts = [dict(contrast_id=cid, label=cid.replace('_vs_',' versus '),
            numerator={'label':conditions[n],'sample_ids':[s['sample_id'] for s in samples if s['condition']==conditions[n]]},
            denominator={'label':conditions[d],'sample_ids':[s['sample_id'] for s in samples if s['condition']==conditions[d]]}) for cid,n,d in contrast_specs]
    sample_count = len(samples)
    project = dict(project_id='katana-integration-dummy-20261102',
        name='Katana | SYNTHETIC DEMONSTRATION - NOT BIOLOGICAL RESULTS',
        samples=samples, contrasts=contrasts, assets=[],
        annotation={'genome_build':'GRCh38','annotation_release':'Dummy annotation; no FASTA or GTF analysed','species':'Homo sapiens'},
        metadata={'species':'Homo sapiens','synthetic':True,'gene_filter_policy':{'footnote_en':'Dummy pre-analysis exclusion: mitochondrial and pseudogene entries excluded; complete detection universe retained.'}})
    def table(name, frame, chapter, cid='all'):
        path = source/(name+'.tsv'); frame.to_csv(path,sep='\t',index=False)
        views['tables'].append(dict(name=name,path=str(path),chapter=chapter,contrast_id=cid,
            metadata={'synthetic':True}))
        return path
    def plot(name,title,chapter,frame,draw,cid='all',policy=None):
        data = table(name,frame,chapter,cid)
        fig, ax = plt.subplots(figsize=(9.8,4.4))
        draw(ax)
        ax.set_title(title,loc='left',pad=18)
        fig.text(.02,.015,'SYNTHETIC DEMONSTRATION - NOT BIOLOGICAL RESULTS',fontsize=8,color='#5b6470')
        fig.tight_layout(rect=(0,.05,1,1))
        path=source/(name+'.png'); fig.savefig(path,dpi=140); fig.savefig(path.with_suffix('.svg')); plt.close(fig)
        panel=dict(id=name,title=title,chapter=chapter,path=str(path),contrast_id=cid,
            policy=policy,caption=warning,metadata={'synthetic':True},
            data_links=[dict(role='plot_data',path=str(data),label='Exact displayed dummy data',unit='synthetic',contrast_id=cid,policy=policy)])
        views['panels'].append(panel)
    def heat(ax,values,labels=None):
        image=ax.imshow(values,aspect='auto',cmap='RdBu_r',vmin=-2,vmax=2)
        if labels is not None: ax.set_yticks(range(len(labels)),labels)
        ax.figure.colorbar(image,ax=ax,label='Illustrative standardized value')
    table('samples',pd.DataFrame(samples),1)
    table('contrasts',pd.DataFrame([{'contrast_id':x['contrast_id'],'numerator':x['numerator']['label'],
        'denominator':x['denominator']['label'],'numerator_samples':';'.join(x['numerator']['sample_ids']),
        'denominator_samples':';'.join(x['denominator']['sample_ids'])} for x in contrasts]),1)
    qc=pd.DataFrame({'cycle':range(1,151),'mean_phred':36-np.linspace(0,4,150)+rng.normal(0,.15,150)})
    plot('sequencing_quality','Sequencing quality by cycle',1,qc,lambda a:a.plot(qc.cycle,qc.mean_phred,color=blue))
    qc=pd.DataFrame({'sample_id':[s['sample_id'] for s in samples],'mapped_percent':rng.uniform(88,97,sample_count),'read_length':150})
    plot('mapping_quality','Read mapping and sample quality',1,qc,lambda a:(a.bar(range(sample_count),qc.mapped_percent,color=blue),a.set_ylabel('Mapped reads (%)')))
    cov=pd.DataFrame({'gene_body_percent':range(0,101),'relative_coverage':1+.1*np.sin(np.linspace(0,np.pi*2,101))})
    plot('gene_body','Gene-body coverage | 5-prime to 3-prime',1,cov,lambda a:(a.plot(cov.gene_body_percent,cov.relative_coverage,color=blue),a.set_xlabel('Gene body (%)')))
    pca=pd.DataFrame({'sample_id':[s['sample_id'] for s in samples],'condition':[s['condition'] for s in samples],
        'PC1':np.repeat(np.arange(len(conditions))-2,sample_count//len(conditions))+rng.normal(0,.3,sample_count),'PC2':rng.normal(0,1,sample_count)})
    plot('psi_pca',f'PSI PCA | all {sample_count} individual samples',2,pca,lambda a:[a.scatter(g.PC1,g.PC2,label=k) for k,g in pca.groupby('condition')]+[a.legend(fontsize=8)])
    psi=rng.uniform(.1,.9,(35,sample_count)); matrix=pd.DataFrame(psi,columns=[s['sample_id'] for s in samples]); matrix.insert(0,'event_id',[f'DEMO_EVENT_{i}' for i in range(35)])
    plot('psi_heatmap','PSI heatmap | individual samples',2,matrix,lambda a:heat(a,(psi-.5)*4))
    clusters=pd.DataFrame({'sample_index':range(sample_count),'pattern_1':np.sin(np.arange(sample_count)/6),'pattern_2':np.cos(np.arange(sample_count)/6)})
    plot('splicing_patterns','Splicing pattern clusters',2,clusters,lambda a:(a.plot(clusters.sample_index,clusters.pattern_1,color=blue,label='Pattern 1'),a.plot(clusters.sample_index,clusters.pattern_2,color=red,label='Pattern 2'),a.legend()))
    concord=pd.DataFrame({'sample_id':[s['sample_id'] for s in samples],'within_group_correlation':rng.uniform(.78,.98,sample_count)})
    plot('replicate_consistency','Within-group consistency and sample diagnostics',2,concord,lambda a:(a.bar(range(sample_count),concord.within_group_correlation,color=blue),a.set_ylim(0,1)))
    types=['SE','A5SS','A3SS','MXE','RI']
    genes=['PTBP1','SRSF1','HNRNPA1','RBFOX2','MBNL1','QKI','ELAVL1','SRSF3','U2AF2','YTHDC1']+[f'DEMO_GENE_{i:03}' for i in range(490)]
    table('detected_testable_gene_background',pd.DataFrame({'gene':genes,'detected':True,'mitochondrial':False,'pseudogene':False}),5)
    artifacts = []
    for ci,(cid,_,_) in enumerate(contrast_specs):
        comp=pd.DataFrame({'event_type':types,'raw':[900,230,280,160,210],'read_support':[650,170,190,125,160]})
        plot(cid+'_composition','Raw and read-support event types | '+cid,2,comp,lambda a:(a.bar(np.arange(5)-.18,comp.raw,.36,color=blue,label='Raw'),a.bar(np.arange(5)+.18,comp.read_support,.36,color=red,label='Read support'),a.set_xticks(range(5),types),a.legend()),cid)
        plot(cid+'_pie','Read-support event type distribution | '+cid,2,comp,lambda a:a.pie(comp.read_support,labels=types,autopct='%1.0f%%'),cid)
        for layer in ['DEG','DAS','DTU','DEXSeq']:
            effect=rng.normal(0,.2 if layer in ['DAS','DTU'] else 1.2,100)
            q=10**(-rng.uniform(.05,5,100)); order=np.argsort(q)
            frame=pd.DataFrame({'gene':genes[:100], 'object_id':[f'DEMO_{layer}_{i:03}' for i in range(100)],'effect':effect,'P':q/2,'FDR':q})
            if layer=='DAS': frame['event_type']=np.resize(types,100)
            for policy in (['stable_significant','significant'] if layer=='DAS' else [None]):
                subset=frame.iloc[:80] if policy=='stable_significant' else frame
                token=cid+'_'+layer+('_RS' if policy=='stable_significant' else '')
                label='RS-DAS' if policy=='stable_significant' else layer
                def volcano(a, f=subset):
                    color=np.where(f.FDR>.05,'#bac1c9',np.where(f.effect>0,red,blue))
                    a.scatter(f.effect,-np.log10(f.FDR),c=color,s=20)
                    a.axhline(-np.log10(.05),ls='--',color='#888',lw=.7)
                    a.set(xlabel='Effect (delta PSI / delta proportion or log2 fold change)',ylabel='-log10 FDR')
                    for r in f.nsmallest(5,'FDR').itertuples(): a.annotate(r.gene,(r.effect,-np.log10(r.FDR)),fontsize=8)
                plot(token+'_volcano',label+' volcano | '+cid,3,subset,volcano,cid,policy)
                ids=contrasts[ci]['numerator']['sample_ids']+contrasts[ci]['denominator']['sample_ids']
                hm=rng.normal(0,1,(10,len(ids))); hmf=pd.DataFrame(hm,columns=ids); hmf.insert(0,'gene',subset.nsmallest(10,'FDR').gene.to_list())
                plot(token+'_heatmap',label+' sample-level heatmap | '+cid,3,hmf,lambda a,v=hm,g=hmf.gene.to_list(),ids=ids:(heat(a,v,g),a.set_xticks(range(len(ids)),ids,rotation=45)),cid,policy)
                complete=source/(token+'_volcano.tsv')
                top=subset.nsmallest(10,'FDR').rename(columns={'gene':'gene_symbol','object_id':'event_id' if layer=='DAS' else 'feature_id','FDR':'fdr'})
                views['top10'].append(dict(id=token,layer=layer,kind=label,title=label+' Top 10 | '+cid,chapter=3,contrast_id=cid,
                    policy=policy,columns=list(top.columns),rows=top.to_dict('records'),
                    full_path=str(complete),total_rows=len(subset),ranking='Illustrative FDR ascending; effect retained'))
        for lib in ['GO_BP','GO_MF','GO_CC','KEGG']:
            terms=pd.DataFrame({'term':[lib+' demo term '+str(i+1) for i in range(6)],'FDR':[.0003,.001,.004,.008,.02,.04],
                'hit_genes':[';'.join(genes[i:i+8]) for i in range(6)],'input_genes':40,'background_genes':500,'comparison':cid})
            plot(cid+'_'+lib,lib.replace('_',' ')+' | RS-DAS genes | '+cid,5,terms,
                lambda a,t=terms:(a.barh(t.term,-np.log10(t.FDR),color=blue),a.set_xlabel('-log10 FDR')),cid,'stable_significant')
        foreground=table(cid+'_enrichment_input_genes',pd.DataFrame({'gene':genes[:40],'selection':'RS-DAS unique genes'}),5,cid)
        for panel in views['panels']:
            if panel.get('chapter')==5 and panel.get('contrast_id')==cid:
                panel['data_links'] += [dict(role='foreground',path=str(foreground),label='Input unique genes',contrast_id=cid),
                    dict(role='background',path=str(source/'detected_testable_gene_background.tsv'),label='Detected/testable gene background',contrast_id=cid)]
        if args.content_review:
            clip=pd.DataFrame({'position_nt':range(-100,101),'synthetic_overlap':np.exp(-np.linspace(-3,3,201)**2)})
            plot(cid+'_clip','Matched CLIP overlap | synthetic example | '+cid,7,clip,
                 lambda a,t=clip:(a.plot(t.position_nt,t.synthetic_overlap,color=blue),a.set(xlabel='Position relative to splice boundary (nt)',ylabel='Illustrative overlap')),cid,'stable_significant')
        rbp=pd.DataFrame({'set':['RBP','Splicing factors'],'P':[.001,.004],'FDR':[.002,.004],
            'hit_genes':['PTBP1;SRSF1;HNRNPA1;RBFOX2','SRSF1;SRSF3;U2AF2'],'background':'500 detected/testable genes'})
        if not args.content_review:
            plot(cid+'_rbp_identity','RBP gene-set enrichment | '+cid,7,rbp,lambda a:a.barh(rbp['set'],-np.log10(rbp.FDR),color=blue),cid,'stable_significant')
        # Deliberately supplied layout records. No native statistical receipt is invented.
        folder=source/cid; folder.mkdir()
        if args.content_review:
            from report_content_fixture import add_mapp_and_rbp
            add_mapp_and_rbp(project, cid, ci, folder, views, artifacts, rng)
        words=['UUUCUCU','UGCAUG','ACUAAC','GGAGGA','CUCUCU','UAAUAA','AGGAAG','CCUCCU','UGUGUG','AUAAUA','GCUAGC','CACACA']
        models={f'DEMO_PWM_{i:02}':{'kind':'PWM','matrix':[[.85 if b==c else .05 for b in 'ACGU'] for c in w],
            'RBPs':['PTBP1','RBFOX2','QKI','SRSF1'][i%4],
            'assignments':[{'RBP':['PTBP1','RBFOX2','QKI','SRSF1'][i%4],'species':'Homo sapiens','source':'Synthetic annotation fixture','native_id':f'DEMO_{i}'}]} for i,w in enumerate(words)}
        plan=dict(comparison_id=cid,models=models,synthetic_display_fixture=True,analysis_engine='eventgraph_20261030')
        (folder/'ANALYSIS_PLAN_BEFORE_SCORING.json').write_text(json.dumps(plan),encoding='utf8')
        rows=[]; windows=[]; events=[]; changes=[]; matches=[]
        for i,token in enumerate(models):
            q=.6 if ci>=2 else (.001*(i+1) if ci==0 else .4)
            p=.7 if ci>=2 else (.002*(i+1) if ci==0 else .01*(i+1))
            rows.append(dict(model_id=token,event_type=types[i%5],direction='inclusion_up' if i%2==0 else 'inclusion_down',
                layer='affected_sequence' if i<8 else 'regulatory_pre_mRNA',region='affected_inclusion' if i<8 else 'acceptor_flank',
                covered_region_start0=i*5,covered_region_end0=i*5+20,cluster_P_FWER_new_family=p,synthetic=True))
            windows.append(dict(model_id=token,event_type=types[i%5],direction=rows[-1]['direction'],region=rows[-1]['region'],
                window_start0=i*5,window_nt=20,raw_P_new_paired_swaps=q/2,q_new_window_family=q,status='estimable',method='synthetic_display_only',synthetic=True))
            for j in range(3):
                eid=f'{cid}_E{i:02}_{j}'
                events.append(dict(event_id=eid,gene_id=genes[i+j],gene_symbol=genes[i+j],event_type=types[i%5]))
                matches.append(dict(foreground_event=eid,status='matched',foreground_gene=genes[i+j],background_gene=f'DEMO_BG_{i}_{j}'))
                changes.append(dict(event_id=eid,model_id=token,instance_state='inclusion_only',consequence_direction='increased'))
        if args.content_review:
            if ci==3:
                rows=[]; windows=[]
                events=events[:4]; matches=matches[:4]; changes=changes[:4]
            (folder/'RESULTS.json').write_text(json.dumps({'comparison_id':cid,'eligible_new_window_tests':len(windows),
                'synthetic':True,'reason':'insufficient independent gene blocks' if ci==3 else 'synthetic completed tests'}),encoding='utf8')
        for name,data in [('NEW_CLUSTER_FWER.tsv',rows),('NEW_WINDOW_AND_FULL_SEQUENCE_TESTS.tsv.gz',windows),
            ('FROZEN_SELECTED_EVENTS.tsv',events),('MATCHING_ALL_SELECTED_WITH_LOSSES.tsv',matches),('PAIRED_FORM_INSTANCE_CHANGES.tsv.gz',changes)]:
            columns = (['model_id','cluster_P_FWER_new_family'] if name=='NEW_CLUSTER_FWER.tsv'
                       else ['model_id','q_new_window_family','status']) if not data else None
            pd.DataFrame(data,columns=columns).to_csv(folder/name,sep='\t',index=False)
        for path in folder.iterdir(): views['tables'].append(dict(name=path.name,path=str(path),chapter=7,contrast_id=cid,kind='file',metadata={'synthetic':True}))
    overlap_labels=['C1 only','C2 only','C3 only','C4 only','C1 & C2','C3 & C4','All four'] if args.content_review else ['B only','C only','E only','B & C','B & E','C & E','B & C & E']
    overlaps=pd.DataFrame({'intersection':overlap_labels,'events':[43,35,22,28,17,11,9]})
    if args.content_review:
        table('intersection_comparison_key',pd.DataFrame({'key':[f'C{i+1}' for i in range(4)],'contrast_id':[c[0] for c in contrast_specs]}),4)
    plot('exact_event_intersections','Exact-event cross-comparison intersections',4,overlaps,lambda a:a.bar(overlaps.intersection,overlaps.events,color=blue))
    dex=pd.DataFrame({'intersection':overlaps.intersection,'exon_bins':[32,24,21,18,12,8,5]})
    plot('exon_intersections','DEXSeq exon-bin cross-comparison intersections',4,dex,lambda a:a.bar(dex.intersection,dex.exon_bins,color=blue))
    interaction=pd.DataFrame({'gene':genes[:12],'interaction_effect':rng.normal(0,.2,12),'FDR':np.linspace(.001,.5,12)})
    plot('interaction','Declared-factor interaction | illustrative output',4,interaction,lambda a:a.barh(interaction.gene,interaction.interaction_effect,color=blue))
    # Import the user's current CoSplice dummy as a portable artifact. Its own
    # source files explicitly identify synthetic modules/statistics and no fit.
    network=source/'cosplice'
    if args.content_review:
        from report_content_fixture import build_network
        build_network(network, project)
    else:
        if not args.network_demo: ap.error('--network-demo is required outside --content-review')
        shutil.copytree(args.network_demo,network,ignore=shutil.ignore_patterns('*.zip','*.png','__pycache__'))
    pages=list(network.glob('*.html'))
    preferred=next((p for p in pages if p.name=='REPORT.html'),pages[0])
    # The copied page advertises its portable ZIP; keep that relative link valid.
    import zipfile
    with zipfile.ZipFile(network/'COSPLICE_DUMMY.zip','w',zipfile.ZIP_DEFLATED) as z:
        for p in network.rglob('*'):
            if p.is_file() and p.suffix!='.zip': z.write(p,p.relative_to(network))
    artifact=dict(artifact_id='dummy-whole-study-network',project_id=project['project_id'],analysis_id='programs.cosplice',
        name='CoSplice whole-study modules',kind='file',entity_type='cosplice_report',path=str(preferred),
        run_fingerprint='synthetic-display-only',fingerprint=hashlib.sha256(preferred.read_bytes()).hexdigest(),
        contrast_ids=[],metadata={'synthetic':True,'preserve_relative_report_links':True,'portable_root':str(network),
            'display_in_main':True,'report_role':'whole_study_chapter','contrast_id':'all'})
    views.update(project_id=project['project_id'],status='completed',
        source_project_manifest_sha256=hashlib.sha256(json.dumps(project,sort_keys=True,separators=(',',':'),default=str).encode()).hexdigest(),
        synthetic=True,audit={'scope':'generated presentation fixture; no scientific execution'})
    descriptor=dict(project=project,artifacts=[artifact,*artifacts],extra_views=[views],synthetic=True,
        report_contrast_groups=[{'group_id':'selected','label':'Declared demonstration comparisons','contrast_ids':[c[0] for c in contrast_specs]}])
    (root/'descriptor.json').write_text(json.dumps(descriptor,indent=2),encoding='utf8')
    from katana_next.reports.cli import main as render
    code=render(['--project',str(root/'descriptor.json'),'--outdir',str(root/'report')])
    if code: raise SystemExit(code)
    # Presentation-only treatment of a synthetic review. Keep the renderer's
    # untouched output and unavailable records; never invent native receipts.
    # Missing real scientific inputs are expected here, not failed demo panels.
    page_path=root/'report'/'REPORT.html'
    original_page=page_path.read_text(encoding='utf8')
    (root/'renderer-original.html').write_text(original_page,encoding='utf8')
    notices=re.findall(r'<div class="status-panel">.*?</div>',original_page,flags=re.S)
    page=re.sub(r'<div class="status-panel">.*?</div>','',original_page,flags=re.S)
    page=page.replace('COMPLETED RESULTS · BULK','SYNTHETIC DISPLAY REVIEW · BULK')
    banner=('<p class="policy-explanation"><strong>Synthetic demonstration.</strong> '
        'All values, statistical tests and network memberships shown here are generated display fixtures. '
        'No biological analysis or end-to-end pipeline acceptance is claimed. '
        + ('MAPP activities and CLIP tracks are synthetic display examples; no fitting or binding evidence is claimed. ' if args.content_review else 'Optional MAPP and matched CLIP calculations are not represented. ')+
        'See the appendix for the unmodified renderer availability records.</p>')
    page=page.replace('<div class="scopebar">',banner+'<div class="scopebar">',1)
    details=('<details><summary>Real-analysis availability records (expected for this dummy)</summary>'
        '<p>The fixture supplies illustrations, not completed native scientific inputs. '
        'The following original renderer records are retained unchanged.</p>'+''.join(notices)+'</details>')
    page=page.replace('</section></main>',details+'</section></main>',1)
    page_path.write_text(page,encoding='utf8')
    (root/'report'/'DUMMY_REVIEW.json').write_text(json.dumps(dict(
        synthetic=True,science_run=False,renderer_original_sha256=hashlib.sha256(original_page.encode()).hexdigest(),
        display_only_changes=['synthetic header and banner','availability notices moved to appendix'],
        original_availability_notices=notices,production_renderer_modified=False),indent=2),encoding='utf8')
    with zipfile.ZipFile(root/'report'/'REPORT.zip','w',zipfile.ZIP_DEFLATED) as z:
        for p in sorted((root/'report').rglob('*')):
            if p.is_file() and p != root/'report'/'REPORT.zip': z.write(p,p.relative_to(root/'report'))
    (root/'DUMMY_SCOPE.json').write_text(json.dumps(dict(synthetic=True,science_run=False,
        purpose='Seven-chapter presentation review; not evidence of native pipeline integration',
        renderer_file=__import__('katana_next.reports.engine',fromlist=['x']).__file__,
        samples=sample_count,contrasts=len(contrasts),network='Synthetic whole-study module fixture',
        motif='12 supplied dummy PWM models per contrast: positive, regional-only and null examples',
        main_changed=False,real_data_read=False),indent=2),encoding='utf8')


if __name__=='__main__': main()
