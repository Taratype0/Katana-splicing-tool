#!/usr/bin/env python3
"""Installed report CLI acceptance from generated native tables; no scientific fit."""
import pandas as pd
from pathlib import Path
import argparse
def verify(tmp_path):
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



if __name__ == "__main__":
    parser=argparse.ArgumentParser(); parser.add_argument("--outdir",type=Path,required=True)
    args=parser.parse_args(); args.outdir.mkdir(parents=True,exist_ok=False)
    verify(args.outdir)
    print("Approved report CLI: 8 native panels, original P/FDR and data links verified")
