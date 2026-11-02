"""Gene ORA conditional on the actual network input, including unassigned events."""
from pathlib import Path
import json,hashlib
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
from validation24 import single_gene

ROOT=Path(__file__).resolve().parents[1]

def read_gmt(path,library,category):
    rows=[]
    for line in Path(path).read_text(encoding='utf8').splitlines():
        name,url,*genes=line.split('\t')
        rows.append(dict(library=library,category=category,pathway=name,url=url,genes=set(genes)))
    return rows

def resources():
    p=ROOT/'data/revision42/annotation';terms=[]
    for c in ['bp','mf','cc']:
        terms+=read_gmt(p/f'c5.go.{c}.v2025.1.Hs.symbols.gmt','GO','GO '+c.upper())
    terms+=read_gmt(ROOT/'data/revision24/annotation/c2.cp.reactome.v2025.1.Hs.symbols.gmt','Reactome','Reactome')
    # Canonical KEGG symbol only: no ambiguous or fuzzy alias expansion.
    symbols={}
    for line in (p/'kegg_genes.txt').read_text(encoding='utf8').splitlines():
        cols=line.split('\t'); description=cols[-1]
        if ';' in description:
            symbol=description.split(';',1)[0].split(',',1)[0].strip()
            if single_gene(symbol):symbols[cols[0]]=symbol
    names=dict(line.split('\t',1) for line in (p/'kegg_pathways.txt').read_text(encoding='utf8').splitlines())
    sets={key.removeprefix('path:'):set() for key in names}
    for line in (p/'kegg_links.txt').read_text(encoding='utf8').splitlines():
        gene,term=line.split('\t');term=term.removeprefix('path:')
        if gene in symbols and term in sets:sets[term].add(symbols[gene])
    for key,genes in sets.items():
        name=names.get(key,names.get('path:'+key,key)).removesuffix(' - Homo sapiens (human)')
        terms.append(dict(library='KEGG',category='KEGG',pathway=key+' | '+name,url='https://www.kegg.jp/pathway/'+key,genes=genes))
    return terms

def calculate(assignment,annotation,terms,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    if not assignment.index.is_unique or not annotation.index.is_unique:raise ValueError('Duplicate event IDs')
    network=annotation.reindex(assignment.index)
    background=set(network[network.map(single_gene)])
    pd.DataFrame({'gene_id':sorted(background)}).to_csv(out/'enrichment_network_genes.tsv',sep='\t',index=False)
    modules={int(m):set(network[assignment.module.eq(m)&network.map(single_gene)]) for m in sorted(assignment.module.dropna().unique()) if m>0}
    libraries=sorted({t['library'] for t in terms}); audits=[];rows=[]
    mapping=pd.DataFrame({'gene_id':sorted(background)})
    for library in libraries:
        collection=[t for t in terms if t['library']==library]
        annotated=set().union(*(t['genes'] for t in collection));universe=background&annotated
        mapping[library+'_mapped']=mapping.gene_id.isin(universe)
        eligible=[(t,t['genes']&universe) for t in collection if 10<=len(t['genes']&universe)<=500]
        for m,all_genes in modules.items():
            foreground=all_genes&universe
            audits.append(dict(module=m,library=library,network_genes=len(background),background_genes=len(universe),unmapped_network_genes=len(background-universe),all_module_genes=len(all_genes),module_genes=len(foreground),eligible_terms=len(eligible),status='tested' if len(foreground)>=5 else 'insufficient_mapped_genes'))
            if len(foreground)<5:continue
            for term,genes in eligible:
                hits=foreground&genes
                rows.append(dict(module=m,library=library,category=term['category'],pathway=term['pathway'],url=term['url'],module_genes=len(foreground),all_module_genes=len(all_genes),background_genes=len(universe),network_genes=len(background),pathway_background_genes=len(genes),overlap=len(hits),genes=';'.join(sorted(hits)),p_value=float(hypergeom.sf(len(hits)-1,len(universe),len(genes),len(foreground)))))
    d=pd.DataFrame(rows)
    if d.empty:raise ValueError('No estimable enrichment tests')
    # Include all discovered modules, nonsignificant terms and zero-overlap tests.
    d['q_value']=multipletests(d.p_value,method='fdr_bh')[1]
    d['fdr_family']='all discovered modules x all eligible terms across supplied libraries'
    d.to_csv(out/'pathway_enrichment.tsv',sep='\t',index=False)
    mapping.to_csv(out/'enrichment_gene_mapping.tsv',sep='\t',index=False)
    pd.DataFrame(audits).to_csv(out/'enrichment_backgrounds.tsv',sep='\t',index=False)
    summary=dict(tests=len(d),tested_modules=len(modules),background_genes=len(background),network_events=len(assignment),assigned_events=int(assignment.module.gt(0).sum()),gene_set_collection='GO BP/MF/CC + KEGG human pathways + Reactome',background_rule='unique unambiguous genes from every selected network event, including unassigned events, intersected with each library annotation',foreground_rule='all unique unambiguous module genes intersected with the same library universe; not Top 10',fdr_family=d.fdr_family.iloc[0],term_size_rule='10 to 500 genes after intersection with library-mapped network background',mapping_rule='exact symbols; KEGG canonical symbols only; no ambiguous aliases',libraries=audits,significant_counts=d[d.q_value.le(.05)].groupby(['module','category']).size().rename('terms').reset_index().to_dict('records'),interpretation='Functional annotation relative to selected network genes. Not causal regulation, splice-specific mechanism, or independent biological validation; selection and gene event multiplicity can still bias enrichment.')
    summary['significant_modules']=int(d[d.q_value.le(.05)].module.nunique())
    summary['significant_pathways']=int(d[d.q_value.le(.05)].pathway.nunique())
    (out/'enrichment_manifest.json').write_text(json.dumps(summary,indent=2),encoding='utf8')
    return summary

def enrich_single(assignment,annotation,path,out):
    return calculate(assignment,annotation,read_gmt(path,Path(path).stem,Path(path).stem),out)

def main():
    import shutil
    run=ROOT/'runs/20261102-revision38-wgcna5000-donor-inference';evidence=ROOT/'evidence/revision42';evidence.mkdir(exist_ok=True)
    cfg=json.loads((run/'config.json').read_text())
    assignment=pd.read_csv(run/'model/assignments.tsv',sep='\t',index_col=0)
    audit=pd.read_csv(run/'event_audit.tsv',sep='\t',index_col=0)
    assert len(assignment)==5000 and set(assignment.index)==set(audit.index[audit.network_selected])
    annotation=pd.read_csv(cfg['annotation'],sep='\t',dtype=str).set_index('event_key').gene_id
    for name in ['pathway_enrichment.tsv','summary.json','report.html']:
        backup=evidence/(name+'.before')
        if not backup.exists():shutil.copyfile(run/name,backup)
    summary=calculate(assignment,annotation,resources(),run)
    provenance=json.loads((ROOT/'data/revision42/annotation/sources.json').read_text())
    reactome=Path(cfg['gene_sets']);provenance['sources'].append(dict(file=reactome.name,url='https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/'+reactome.name,sha256=hashlib.sha256(reactome.read_bytes()).hexdigest()))
    (run/'enrichment_sources.json').write_text(json.dumps(provenance,indent=2))
    original=json.loads((run/'summary.json').read_text());original['pathway_enrichment']=summary
    (run/'summary.json').write_text(json.dumps(original,indent=2),encoding='utf8')
    print(json.dumps({k:v for k,v in summary.items() if k!='libraries'},indent=2))
if __name__=='__main__':main()
