"""Small synthetic display fixtures. No scientific fitting or real input reads."""
from pathlib import Path
import json
import runpy
import zipfile
import numpy as np
import pandas as pd


def add_mapp_and_rbp(project, cid, ci, folder, views, artifacts, rng):
    from katana_next.reports.mapp_activity import render_activity
    contrast = project['contrasts'][ci]
    ids = contrast['numerator']['sample_ids'] + contrast['denominator']['sample_ids']
    frame = pd.DataFrame({'site': ['3ss', '5ss'] * 10,
                          'native_model_id': [f'DEMO_PWM_{i:02}' for i in range(20)]})
    for sid in ids:
        frame['activities_' + sid] = rng.normal(0, .65, 20)
    data = folder/'mapp-activities.tsv'; frame.to_csv(data, sep='\t', index=False)
    path = folder/'mapp-activities.png'
    axis = render_activity(frame, project, cid, path)
    views['panels'].append(dict(id=cid+'-mapp', title='MAEI motif activity | '+cid,
        chapter=7, contrast_id=cid, path=str(path),
        caption='SYNTHETIC DISPLAY FIXTURE. Six individual samples; condition labels and group separator. No MAPP fit performed.',
        metadata={'synthetic': True, 'activity_axis': axis},
        data_links=[dict(role='plot_data', path=str(data), label='All displayed synthetic activities', contrast_id=cid, unit='sample')]))
    rows=[]
    for i, genes in enumerate(['PTBP1;SRSF1;HNRNPA1;RBFOX2','SRSF1;SRSF3;U2AF2','QKI;MBNL1']):
        rows.append(dict(term=['RBP identity','Splicing factors','Custom declared RBP set'][i],
            genes=genes, overlap_genes=len(genes.split(';')), pvalue=[.001,.008,.4][i] if ci<2 else .4+i*.1,
            fdr=[.003,.012,.4][i] if ci<2 else .6+i*.1, fold_enrichment=3-i*.7,
            selected_genes=40, background_genes=500, comparison_id=cid, fdr_family='Synthetic three-set family'))
    path=folder/'rbp_enrichment.parquet'; pd.DataFrame(rows).to_parquet(path,index=False)
    artifacts.append(dict(artifact_id=cid+'-rbp', project_id=project['project_id'],
        analysis_id='prioritize.rbp_support', name='Synthetic RBP identity enrichment', kind='table',
        path=str(path), contrast_ids=[cid], run_fingerprint='synthetic-display-only',
        metadata={'synthetic':True,'foreground_source':'DEG unique genes (synthetic)',
                  'background_definition':'500 detected/testable genes in this comparison (synthetic)',
                  'resource_semantics':'SYNTHETIC DATA. RBP identity is not RBP-target regulation.'}))


def build_network(directory, project):
    directory.mkdir(); tables=directory/'report_tables'; tables.mkdir()
    rng=np.random.default_rng(20261104)
    samples=project['samples']; modules=[]; members=[]; genes=[]; associations=[]
    graphs={}; usage={}; files={}; pathways=[]; libraries=[]
    labels=['PTBP1','SRSF1','HNRNPA1','RBFOX2','MBNL1','QKI','ELAVL1','SRSF3','U2AF2','YTHDC1']
    for i in range(1,5):
        m=f'M{i}'; n=20+5*i
        modules.append(dict(module=m,events=n,pc1_variance_fraction=.6-i*.04,
                            latent_factor='Synthetic shared programme', loading_direction='positive'))
        for j in range(n):
            key=f'DEMO_{m}_SE_{j+1:03}'; symbol=labels[j%len(labels)]
            members.append(dict(module=m,event_key=key,gene_id=symbol,event_type='SE',
                kME=round(.95-j*.01,3),within_module_weight=n-j,absolute_loading=.8,
                loading_margin=.4,condition_mean_psi_range=.2,observed_samples=12))
            usage[key]=rng.uniform(.15,.85,12).round(3).tolist()
        mm=[r for r in members if r['module']==m]
        for symbol in labels:
            gg=[r for r in mm if r['gene_id']==symbol]
            genes.append(dict(module=m,gene_id=symbol,events=len(gg),max_abs_kME=max(x['kME'] for x in gg),
                              within_module_weight=sum(x['within_module_weight'] for x in gg)))
        edges=[[j,(j+1)%n,.8] for j in range(n)]+[[j,(j+5)%n,.45] for j in range(0,n,3)]
        graphs[m]=dict(nodes=[r['event_key'] for r in mm],edges=edges,total_nodes=n,total_edges=len(edges))
        files[m]=dict(genes=m+'-genes.tsv', events=m+'-events.tsv')
        pd.DataFrame([r for r in genes if r['module']==m]).to_csv(tables/files[m]['genes'],sep='\t',index=False)
        pd.DataFrame(mm).to_csv(tables/files[m]['events'],sep='\t',index=False)
        pd.DataFrame(mm).to_csv(tables/(m+'-nodes.tsv'),sep='\t',index=False)
        pd.DataFrame([dict(source=mm[a]['event_key'],target=mm[b]['event_key'],weight=w) for a,b,w in edges]).to_csv(tables/(m+'-edges.tsv'),sep='\t',index=False)
        scopes=[dict(kind='omnibus',factor='condition')]+[dict(kind='contrast',factor='condition',a=c['numerator']['label'],b=c['denominator']['label']) for c in project['contrasts']]
        for k, scope in enumerate(scopes):
            missing=i==3
            q=.01 if k==1 and i==1 else .12+.15*i
            associations.append(dict(**scope,module=m,status='NA' if missing else 'estimated',
                p_value=None if missing else q/4,q_value=None if missing else q,
                effect=None if missing else .4,ci_low=.1,ci_high=.7,statistic=2,
                qualified_candidate=False,qualification_reason='Stability not assessed',
                reason='Synthetic insufficient score coverage example' if missing else '',
                median_member_support=None,median_match_jaccard=None))
        for lib in ['GO','KEGG','Reactome']:
            libraries.append(dict(module=i,library=lib,module_genes=10,all_module_genes=10,
                background_genes=50,network_genes=50,status='tested'))
        for cat in ['GO BP','GO MF','GO CC','KEGG','Reactome']:
            for j in range(3):
                pathways.append(dict(module=m,category=cat,pathway=f'{cat} synthetic term {j+1}',
                    p_value=.001*(j+1) if i<3 else .4,q_value=.01*(j+1) if i<3 else .6,
                    overlap=3+j,module_genes=10,genes=';'.join(labels[:3+j])))
    scores=[dict(sample_id=s['sample_id'],**{f'M{i}':round(float(rng.normal()),3) if not(i==3 and j==11) else None for i in range(1,5)}) for j,s in enumerate(samples)]
    result=dict(input_events=150,eligible_events=140,assigned_events=len(members),network_events=140,
        model='SYNTHETIC DISPLAY FIXTURE; no network fit',bootstrap_requested=0,
        pathway_enrichment=dict(background_genes=50,network_events=140,tests=len(pathways),tested_modules=4,libraries=libraries))
    data=dict(modules=modules,members=members,genes=genes,associations=associations,graphs=graphs,
        event_usage=usage,module_files=files,samples=samples,scores=scores,pathways=pathways,
        manifest=dict(study_id='SYNTHETIC DISPLAY REVIEW',samples=12,biological_units=6,network_id='synthetic',assignment_sha256='synthetic'),
        result=result,analysis=dict(inference_revision='Synthetic P/FDR supplied for presentation only; no inference performed',
            limits='No scientific acceptance is claimed'),assembly={},qualification={'status':'not_measured'})
    for name, rows in [('pathways',pathways),('module_qualification',associations),('qualified_modules',[]),
                       ('sample_metadata',samples),('module_scores',scores),('module_gene_summary',genes),('gene_event_membership',members)]:
        pd.DataFrame(rows,columns=None if rows else ['module','q_value']).to_csv(tables/(name+'.tsv'),sep='\t',index=False)
    pd.DataFrame.from_dict(usage,orient='index',columns=[s['sample_id'] for s in samples]).to_csv(tables/'module_event_psi.tsv.gz',sep='\t',index_label='event_id')
    for name in ['all_module_lists.zip','qualified_module_lists.zip']:
        with zipfile.ZipFile(tables/name,'w',zipfile.ZIP_DEFLATED) as z:
            for path in sorted(tables.glob('M*.tsv')): z.write(path,path.name)
    data['downloads']=[p.name for p in sorted(tables.iterdir())]
    base=Path(__file__).resolve().parents[1]/'katana_cosplice/_study_core/science'
    page=(base/'study_report39.html').read_text('utf8')
    page=page.replace('__KATANA_CSS__',(base/'katana_report.css').read_text('utf8')).replace('__STUDY_TITLE__','Chapter 06 · SYNTHETIC DISPLAY REVIEW')
    page=page.replace('__REPORT_DATA__',json.dumps(data))
    for i in range(1,5):
        page=page.replace(f'<b>0{i}</b>',f'<b>6.{i}</b>').replace(f'<span class="chapter-number">0{i}</span>',f'<span class="chapter-number">6.{i}</span>')
    page=runpy.run_path(str(base/'no_refit_report.py'))['adapt'](page,result)
    page=page.replace('<section class="cover">','<section class="cover"><p class="signal"><strong>SYNTHETIC DATA.</strong> Illustrative modules, edges, scores and P/FDR; no scientific analysis was run.</p>')
    (directory/'REPORT.html').write_text(page,encoding='utf8')
