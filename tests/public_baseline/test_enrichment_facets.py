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
            for entity,name,frame in [('gene_set','tests',rows),
                ('gene_input_audit','inputs',pd.DataFrame({'gene_id':['G1','G2'],'selected_input':[True,False]})),
                ('background_audit','background',pd.DataFrame({'gene_id':['G1','G2','G3']}))]:
                path=folder/(name+'.parquet'); frame.to_parquet(path,index=False)
                originals[str(path)]=path.read_bytes()
                artifacts.append(dict(artifact_id=foreground+cid+name,analysis_id='prioritize.enrichment',
                    kind='table',entity_type=entity,name=name,path=str(path),contrast_ids=[cid],
                    run_fingerprint=foreground+cid,metadata=dict(completion_verified=True,
                    foreground_source=foreground,background_source=foreground,native_background_provenance={'native_inputs':[]},
                    das_selection={'event_set':'stable_significant','read_support_policy':'bulk-all-samples-counts-v1'})))
    project=tmp_path/'source.json'
    project.write_text(json.dumps({'project':{'project_id':'native-report-fixture','name':'Native saved ORA fixture',
        'samples':[],'contrasts':[{'contrast_id':cid,'label':cid} for cid in ('B_vs_A','D_vs_C')]},
        'artifacts':artifacts,'report_require_plot_data':True}),encoding='utf8')
    out=tmp_path/'report'
    subprocess.run([sys.executable,'-m','katana_next.reports.cli','--project',str(project),'--outdir',str(out)],check=True,capture_output=True,text=True)
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
