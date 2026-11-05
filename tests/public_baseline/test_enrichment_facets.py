import pandas as pd
import pytest
from katana_next.reports.enrichment_facets import render_enrichment_facets


def test_facets_use_saved_counts_ratios_and_adjusted_p(tmp_path):
    rows=pd.DataFrame([dict(category=c,term=f'{c} term {i}',Count=3+i,GeneRatio=(3+i)/40,
        **{'p.adjust':.001*(i+1) if i<6 else .6}) for c in ['BP','CC','MF','KEGG'] for i in range(7)])
    original=rows.copy(deep=True)
    selected=render_enrichment_facets(rows,tmp_path/'facets.png')
    assert len(selected)==20
    assert selected.category.drop_duplicates().tolist()==['BP','CC','MF','KEGG']
    assert (selected['p.adjust']<=.05).all()
    pd.testing.assert_frame_equal(rows,original)
    for row in selected.to_dict('records'):
        source=original.loc[(original.category==row['category'])&(original.term==row['term'])].iloc[0]
        assert all(row[key]==source[key] for key in ['Count','GeneRatio','p.adjust'])
    svg=(tmp_path/'facets.svg').read_text('utf8')
    assert all(x in svg for x in ['Gene Ratio','p.adjust','Counts','BP','CC','MF','KEGG'])


def test_unknown_ontology_is_not_guessed(tmp_path):
    frame=pd.DataFrame([{'category':'Unclassified','term':'RNA processing','Count':4,'GeneRatio':.1,'p.adjust':.01}])
    with pytest.raises(ValueError,match='explicitly recorded'):
        render_enrichment_facets(frame,tmp_path/'bad.png')


def test_empty_facet_is_reported_without_inventing_terms(tmp_path):
    frame=pd.DataFrame([{'category':'BP','term':'nonsignificant','Count':4,'GeneRatio':.1,'p.adjust':.8}])
    result=render_enrichment_facets(frame,tmp_path/'empty.png')
    assert result.empty
    assert 'No significant terms' in (tmp_path/'empty.svg').read_text('utf8')


def test_public_report_cli_uses_native_full_tables_and_preserves_statistics(tmp_path):
    import json
    import subprocess
    import sys
    from pathlib import Path
    artifacts=[]
    originals={}
    for foreground in ('deg','das','dtu','dexseq'):
        for cid in ('B_vs_A','D_vs_C'):
            folder=tmp_path/foreground/cid; folder.mkdir(parents=True)
            rows=pd.DataFrame([dict(comparison_id=cid,term=f'{c} term {i}',ontology_category=c,
                overlap_genes=3+i,selected_genes=40,background_genes=500,
                pvalue=.0001*(i+1),fdr=.001*(i+1) if i<6 else .9,
                fdr_family='original-family-'+c,genes='G1;G2;G3')
                for c in ('GO BP','GO CC','GO MF','KEGG') for i in range(7)])
            for entity,name,frame in [('gene_set','Gene-set enrichment',rows),
                ('gene_input_audit','Enrichment gene inputs by comparison',pd.DataFrame({'gene_id':['G1','G2'],'selected_input':[True,False]})),
                ('background_audit','Enrichment background eligibility',pd.DataFrame({'gene_id':['G1','G2','G3']}))]:
                frame['comparison_id'] = cid
                path=folder/(name+'.parquet'); frame.to_parquet(path,index=False)
                originals[str(path)]=path.read_bytes()
                artifacts.append(dict(artifact_id=foreground+cid+name,analysis_id='prioritize.enrichment',
                    kind='table',entity_type=entity,name=name,path=str(path),
                    run_fingerprint=foreground+cid,metadata=dict(completion_verified=True,
                    foreground_source=foreground,background_source=foreground,native_background_provenance={'native_inputs':[]},
                    das_selection={'event_set':'stable_significant','read_support_policy':'bulk-all-samples-counts-v1'})))
    project=tmp_path/'source.json'
    project.write_text(json.dumps({'project':{'project_id':'native-report-fixture','name':'Native saved ORA fixture',
        'samples':[],'contrasts':[{'contrast_id':cid,'label':cid,'numerator':{'label':cid.split('_vs_')[0],'sample_ids':[]},'denominator':{'label':cid.split('_vs_')[1],'sample_ids':[]}} for cid in ('B_vs_A','D_vs_C')]},
        'artifacts':artifacts,'report_require_plot_data':True}),encoding='utf8')
    out=tmp_path/'report'
    result=subprocess.run([sys.executable,'-m','katana_next.reports.cli','--project',str(project),'--outdir',str(out)],capture_output=True,text=True)
    assert result.returncode==0, result.stderr+result.stdout
    manifest=json.loads((out/'report-manifest.json').read_text())
    panels=[p for p in manifest['panels'] if p.get('metadata',{}).get('renderer')=='go-kegg-facets-v1']
    assert len(panels)==8
    for panel in panels:
        assert panel['chapter']==5 and not panel.get('output_only')
        assert {l['role'] for l in panel['data_links']}=={'plot_data','full_results','foreground','background'}
        link=next(l for l in panel['data_links'] if l['role']=='plot_data')
        selected=pd.read_csv(out/link['path'],sep='\t')
        assert len(selected)==20
        assert set(selected.comparison_id)=={panel['contrast_id']}
        assert (selected['p.adjust']==selected.fdr).all()
        assert (selected.GeneRatio==selected.overlap_genes/selected.selected_genes).all()
        assert set(selected.fdr_family)=={'original-family-'+c for c in ('GO BP','GO CC','GO MF','KEGG')}
    assert all(Path(path).read_bytes()==content for path,content in originals.items())
    assert (out/'REPORT.zip').stat().st_size>0
    # Older portable reports omitted scope metadata on complete table entries.
    # Recover it from exact artifact provenance and the saved comparison column.
    for table in manifest['complete_results']:
        table.pop('completed_enrichment', None)
        table.pop('contrast_id', None)
    for panel in manifest['panels']: panel['output_only'] = True
    (out/'report-manifest.json').write_text(json.dumps(manifest))
    report_before = (out/'REPORT.html').read_bytes()
    refreshed = tmp_path/'refreshed'
    result = subprocess.run([sys.executable,'-m','katana_next.reports.cli','--from-report',str(out),
                             '--outdir',str(refreshed),'--no-archive'],capture_output=True,text=True)
    assert result.returncode == 0, result.stderr+result.stdout
    current = json.loads((refreshed/'report-manifest.json').read_text())
    visible = [p for p in current['panels'] if p.get('metadata',{}).get('renderer')=='go-kegg-facets-v1' and not p.get('output_only')]
    assert len(visible) == 8
    assert {p['contrast_id'] for p in visible} == {'B_vs_A','D_vs_C'}
    for panel in visible:
        assert {r['role'] for r in panel['data_links']} == {'plot_data','full_results','foreground','background'}
    assert (out/'REPORT.html').read_bytes() == report_before
    assert all(Path(path).read_bytes()==content for path,content in originals.items())


def test_native_unclassified_and_missing_tests_are_not_negative_findings(tmp_path):
    from katana_next.reports.native_enrichment import build_native_enrichment_facets
    frame=pd.DataFrame([dict(term='Unclassified term',ontology_category='Unclassified',overlap_genes=4,selected_genes=40,fdr=.001)])
    path=tmp_path/'native.parquet';frame.to_parquet(path)
    a=dict(path=str(path),analysis_id='prioritize.enrichment',kind='table',entity_type='gene_set',
           contrast_ids=['B_vs_A'],metadata=dict(foreground_source='deg',native_background_provenance={'source':'recorded'},report_completion_evidence='completed_workflow_manifest'))
    result=build_native_enrichment_facets([a],tmp_path/'out')
    from pathlib import Path
    svg=Path(result['panels'][0]['path']).with_suffix('.svg').read_text('utf8')
    assert 'No recorded eligible tests' in svg
    assert 'No significant terms' not in svg
    assert result['audit'][0]['unclassified_terms']==1


def test_completed_go_ids_get_offline_namespaces_without_changing_native_statistics(tmp_path):
    from katana_next.reports.native_enrichment import build_native_enrichment_facets
    from katana_next.reports.go_categories import annotate_completed_categories
    from pathlib import Path
    frame=pd.DataFrame([dict(term='GO / '+gid+' saved term',ontology_category='Unclassified',
        overlap_genes=3,selected_genes=20,background_genes=500,fdr=.01,pvalue=.001,
        comparison_id='B_vs_A',library='GO',genes='G1;G2;G3')
        for gid in ('GO:0008150','GO:0005575','GO:0003674')])
    path=tmp_path/'native.parquet';frame.to_parquet(path,index=False);original=path.read_bytes()
    artifact=dict(path=str(path),analysis_id='prioritize.enrichment',kind='table',entity_type='gene_set',
        contrast_ids=['B_vs_A'],metadata=dict(foreground_source='deg',native_background_provenance={'source':'recorded'},
        report_completion_evidence='completed_workflow_manifest'))
    result=build_native_enrichment_facets([artifact],tmp_path/'out')
    assert result['audit'][0]['unclassified_terms']==0
    panel=result['panels'][0];data=pd.read_csv(next(x['path'] for x in panel['data_links'] if x['role']=='plot_data'),sep='\t')
    assert data.category.tolist()==['BP','CC','MF']
    assert data.ontology_category.tolist()==['Unclassified']*3
    assert data.fdr.tolist()==[.01]*3 and data.pvalue.tolist()==[.001]*3
    assert panel['metadata']['go_category_annotation']['rows_resolved']==3
    assert panel['metadata']['go_category_annotation']['ontology_version']=='releases/2026-07-26'
    assert path.read_bytes()==original
    assert 'GO category metadata unavailable' not in Path(panel['path']).with_suffix('.svg').read_text('utf8')
    frame.loc[0,'ontology_category']='GO MF'
    annotated,_=annotate_completed_categories(frame)
    assert annotated.loc[0,'category']=='MF'  # Do not overwrite a recorded category.
