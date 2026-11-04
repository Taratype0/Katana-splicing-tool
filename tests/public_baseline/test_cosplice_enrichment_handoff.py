"""Real-shaped ID/symbol handoff, complete library families and empty outcomes."""
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pandas as pd
import pytest
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests


@pytest.fixture
def engine(monkeypatch):
    from katana_cosplice import _study_core
    root=Path(_study_core.__file__).parent/'science'
    monkeypatch.syspath_prepend(str(root))
    spec=importlib.util.spec_from_file_location('enrichment_handoff_test',root/'enrichment42.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def inputs():
    keys=[f'e{i}' for i in range(40)]
    annotation=pd.DataFrame({'gene_id':[f'ENSMUSG{i:011d}' for i in range(40)],
                             'gene_symbol':[f'Gene{i}' for i in range(40)]},index=keys)
    assignment=pd.DataFrame({'module':[1]*15+[2]*15+[0]*10},index=keys)
    return assignment,annotation


def test_two_symbol_libraries_map_to_native_ids_and_keep_joint_family(engine,tmp_path):
    assignment,annotation=inputs()
    libraries=[]
    for name,starts in [('GO',[0,20]),('KEGG',[10,25])]:
        path=tmp_path/(name+'.gmt')
        path.write_text(''.join(f'{name}_{start}\turl\t'+'\t'.join(f'Gene{i}' for i in range(start,start+15))+'\n' for start in starts))
        libraries.append({'path':str(path),'library':name,'category':name})
    out=tmp_path/'out'
    summary=engine.enrich_libraries(assignment,annotation,libraries,out)
    actual=pd.read_csv(out/'pathway_enrichment.tsv',sep='\t')
    assert summary['status']=='completed' and summary['tests']==8
    assert set(actual.library)=={'GO','KEGG'}
    assert summary['background_genes']==40 # includes unassigned events
    assert set(actual.background_genes)=={30}
    expected=hypergeom.sf(actual.overlap-1,actual.background_genes,actual.pathway_background_genes,actual.module_genes)
    np.testing.assert_allclose(actual.p_value,expected)
    np.testing.assert_allclose(actual.q_value,multipletests(expected,method='fdr_bh')[1])
    assert all(g.startswith('ENSMUSG') for text in actual.genes.dropna() for g in text.split(';'))
    # Legacy exact-ID inputs produce identical tests, without adding a second
    # copy of every gene to the GMT to disguise missing identifier mapping.
    terms=[]
    for item in libraries:
        for term in engine.read_gmt(item['path'],item['library'],item['category']):
            terms.append({**term,'genes':{f'ENSMUSG{int(g[4:]):011d}' for g in term['genes']}})
    engine.calculate(assignment,annotation.gene_id,terms,tmp_path/'legacy')
    legacy=pd.read_csv(tmp_path/'legacy/pathway_enrichment.tsv',sep='\t')
    pd.testing.assert_frame_equal(actual,legacy)


@pytest.mark.parametrize('mode,reason',[
    ('mismatch','no_mapped_background'),('small_terms','no_eligible_terms'),
    ('small_modules','insufficient_mapped_genes'),('unassigned','no_assigned_modules')])
def test_not_estimable_is_a_recorded_outcome_not_nonsignificance(engine,tmp_path,mode,reason):
    assignment,annotation=inputs()
    genes=set(annotation.gene_symbol)
    if mode=='mismatch':genes={'Other'+str(i) for i in range(40)}
    if mode=='small_terms':genes=set(annotation.gene_symbol[:3])
    if mode=='small_modules':assignment['module']=[i//3+1 for i in range(40)]
    if mode=='unassigned':assignment['module']=0
    summary=engine.calculate(assignment,annotation,[dict(library='GO',category='GO',pathway='P',url='',genes=genes)],tmp_path)
    assert summary['status']=='not_estimable' and reason in summary['reasons']
    assert summary['tests']==summary['tested_modules']==0
    result=pd.read_csv(tmp_path/'pathway_enrichment.tsv',sep='\t')
    assert result.empty and {'p_value','q_value','genes','module'}<=set(result)
    assert json.loads((tmp_path/'enrichment_manifest.json').read_text())==summary
    assert (tmp_path/'enrichment_backgrounds.tsv').stat().st_size>0


def test_ambiguous_symbols_not_arbitrarily_assigned(engine,tmp_path):
    assignment,annotation=inputs()
    annotation.loc[['e0','e1'],'gene_symbol']='Shared'
    terms=[dict(library='GO',category='GO',pathway='P',url='',genes={'Shared',*annotation.gene_symbol[2:]})]
    engine.calculate(assignment,annotation,terms,tmp_path)
    mapping=pd.read_csv(tmp_path/'enrichment_identifier_mapping.tsv',sep='\t')
    assert mapping.loc[mapping.token.eq('Shared'),'status'].tolist()==['ambiguous']
    assert mapping.loc[mapping.token.eq('Shared'),'gene_id'].isna().all()


@pytest.mark.parametrize('content',['','broken\n','P\turl\n','P\turl\tA\nP\turl\tB\n'])
def test_invalid_gmt_remains_an_error(engine,tmp_path,content):
    path=tmp_path/'bad.gmt';path.write_text(content)
    with pytest.raises(ValueError):engine.read_gmt(path,'GO','GO')


def test_report_distinguishes_unestimable_and_preserves_legacy_go_categories():
    from katana_cosplice import _study_core
    node=shutil.which('node')
    if not node:pytest.skip('Node is required for JavaScript report behavior')
    source=(Path(_study_core.__file__).parent/'science/enrichment42.js').read_text('utf8')
    harness=r'''
const vm=require('vm'),assert=require('assert');
function run(info){
 const elements={pathwayContext:{},pathwayTable:{replaceChildren(){},append(x){this.child=x}}};
 const select={value:'GO BP'},body={};
 const context={D:{result:{pathway_enrichment:info},pathways:[]},significant:[{module:'M1'}],
  $:id=>elements[id],E:x=>x,document:{createElement:()=>({dataset:{},querySelector:s=>s==='select'?select:body})}};
 vm.runInNewContext(SOURCE+'\nrenderEnrichment39()',context);
 return {elements,body};
}
const base={background_genes:40,network_events:40,tests:0,tested_modules:0,libraries:[]};
const empty=run({...base,status:'not_estimable',reasons:['no_eligible_terms']});
assert(empty.elements.pathwayContext.textContent.includes('not estimable'));
assert(empty.elements.pathwayContext.textContent.includes('No P values or FDR were computed'));
assert(!empty.elements.pathwayContext.textContent.includes('No significant'));
const legacy=run({...base,tests:1,tested_modules:1,libraries:[{module:1,library:'GO',status:'tested',module_genes:10,all_module_genes:15,background_genes:30,network_genes:40}]});
assert(legacy.elements.pathwayTable.child.innerHTML.includes('<option>GO BP</option>'));
assert(legacy.elements.pathwayTable.child.innerHTML.includes('<option>GO MF</option>'));
assert(legacy.body.innerHTML.includes('No significant terms'));
'''
    subprocess.run([node,'-e','const SOURCE='+json.dumps(source)+';\n'+harness],check=True,capture_output=True,text=True)
