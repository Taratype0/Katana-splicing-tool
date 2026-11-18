import json
import pandas as pd
from katana_next.api.report_workspace import BrowseReport
from katana_next.services.report_import import import_report, sha
from katana_next.services.report_query import browse, plot, select_rows
from test_report_workspace import report_fixture


def test_report_names_preserve_source_and_reject_ambiguous_aliases(tmp_path):
    root, _ = report_fixture(tmp_path)
    manifest=json.loads((root/'report-manifest.json').read_text())
    original=pd.DataFrame({'gene_id':['ENSMUSG001','ENSMUSG002','ENSMUSG003'],
        'gene_symbol':['ENSMUSG001.2','ENSMUSG002.4','Native'], 'log2fc':[-2.,5.,0.], 'deg_fdr':[.01,.2,None]})
    original.to_parquet(root/'full.parquet',index=False)
    pd.DataFrame({'alias':['Fos','A','B'],'canonical_gene_ids':['ENSMUSG001','ENSMUSG002','ENSMUSG002'],
                  'mapping_status':['mapped']*3}).to_parquet(root/'gene_identifier_mapping.parquet',index=False)
    manifest['complete_results']=[{'source':n,'name':n,'analysis_id':'differential.deg'} for n in ['full.parquet','gene_identifier_mapping.parquet']]
    manifest['source_files']=[]
    (root/'report-manifest.json').write_text(json.dumps(manifest))
    context=import_report(root,tmp_path/'workspace')
    artifact=next(a for a in context.catalog.list_artifacts(context.manifest.project_id) if a.name=='full.parquet')
    filters=BrowseReport(fdr_column='deg_fdr',effect_column='log2fc')
    rows=browse(artifact,filters)['rows']
    assert [r['gene_name'] for r in rows]==['Fos',None,'Native']
    assert [r['gene_symbol'] for r in rows]==original.gene_symbol.tolist()
    filters.identifiers=['Fos']
    assert select_rows(artifact,filters)['rows'][0]['gene_id']=='ENSMUSG001'
    filters.identifiers=[]
    chart=plot(artifact,filters)
    assert chart['x_min']==-chart['x_max'] and chart['plotted']==2 and chart['unplottable']==1
    assert chart['points'][0][3]=='Fos'
    assert sha(artifact.path)==artifact.metadata['imported_table_sha256']
