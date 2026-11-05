#!/usr/bin/env python3
"""Installed report CLI acceptance from generated native tables; no scientific fit."""
import pandas as pd
from pathlib import Path
import argparse
def verify(tmp_path):
    from katana_next.reports.go_categories import annotate_completed_categories
    categorized, provenance = annotate_completed_categories(pd.DataFrame({
        'term':['GO:0008150 biological process','GO:0005575 cellular component','GO:0003674 molecular function'],
        'ontology_category':['Unclassified']*3, 'fdr':[.01,.02,.03]}))
    assert categorized.category.tolist()==['BP','CC','MF']
    assert categorized.fdr.tolist()==[.01,.02,.03]
    assert provenance['rows_resolved']==3
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


if __name__ == "__main__":
    parser=argparse.ArgumentParser(); parser.add_argument("--outdir",type=Path,required=True)
    args=parser.parse_args(); args.outdir.mkdir(parents=True,exist_ok=False)
    verify(args.outdir)
    print("Approved report CLI: 8 native panels, original P/FDR and data links verified")
