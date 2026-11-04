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
    seen=set()
    for number,line in enumerate(Path(path).read_text(encoding='utf8').splitlines(),1):
        if not line.strip():continue
        fields=line.split('\t')
        if len(fields)<3 or not fields[0].strip() or not any(g.strip() for g in fields[2:]):
            raise ValueError(f'Malformed GMT row {number}: {path}')
        name,url,*genes=fields
        if name in seen:raise ValueError(f'Duplicate GMT term {name}: {path}')
        seen.add(name)
        genes=[g.strip() for g in genes if g.strip()]
        rows.append(dict(library=library,category=category,pathway=name,url=url,genes=set(genes)))
    if not rows:raise ValueError(f'Empty GMT library: {path}')
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
    # Keep native gene identities. Symbols are exact aliases supplied by the
    # event annotation, never fuzzy matches or an expanded background.
    if isinstance(annotation,pd.Series):annotation=annotation.to_frame('gene_id')
    if 'gene_id' not in annotation:raise ValueError('Event annotation requires gene_id')
    aliases={}
    for row in annotation.to_dict('records'):
        gene=row['gene_id']
        if not single_gene(gene):continue
        for token in (gene,row.get('gene_symbol')):
            if single_gene(token):aliases.setdefault(token,set()).add(gene)
    unique={token:next(iter(ids)) for token,ids in aliases.items() if len(ids)==1}
    mapped_terms=[]; token_audit=[]; seen_tokens=set()
    for term in terms:
        mapped=set()
        for token in sorted(term['genes']):
            gene=unique.get(token)
            if gene is not None:mapped.add(gene)
            identity=(term['library'],token)
            if identity not in seen_tokens:
                seen_tokens.add(identity)
                token_audit.append(dict(library=term['library'],token=token,gene_id=gene or '',status='mapped' if gene is not None else 'ambiguous' if token in aliases else 'unmapped'))
        mapped_terms.append({**term,'genes':mapped})
    terms=mapped_terms
    pd.DataFrame(token_audit,columns=['library','token','gene_id','status']).drop_duplicates().to_csv(out/'enrichment_identifier_mapping.tsv',sep='\t',index=False)
    network=annotation.gene_id.reindex(assignment.index)
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
            state='no_mapped_background' if not universe else 'no_eligible_terms' if not eligible else 'insufficient_mapped_genes' if len(foreground)<5 else 'tested'
            audits.append(dict(module=m,library=library,network_genes=len(background),background_genes=len(universe),unmapped_network_genes=len(background-universe),all_module_genes=len(all_genes),module_genes=len(foreground),eligible_terms=len(eligible),status=state))
            if len(foreground)<5:continue
            for term,genes in eligible:
                hits=foreground&genes
                rows.append(dict(module=m,library=library,category=term['category'],pathway=term['pathway'],url=term['url'],module_genes=len(foreground),all_module_genes=len(all_genes),background_genes=len(universe),network_genes=len(background),pathway_background_genes=len(genes),overlap=len(hits),genes=';'.join(sorted(hits)),p_value=float(hypergeom.sf(len(hits)-1,len(universe),len(genes),len(foreground)))))
    d=pd.DataFrame(rows,columns=['module','library','category','pathway','url','module_genes','all_module_genes','background_genes','network_genes','pathway_background_genes','overlap','genes','p_value'])
    # Include all discovered modules, nonsignificant terms and zero-overlap tests.
    d['q_value']=multipletests(d.p_value,method='fdr_bh')[1] if len(d) else pd.Series(dtype=float)
    d['fdr_family']='all discovered modules x all eligible terms across supplied libraries'
    d.to_csv(out/'pathway_enrichment.tsv',sep='\t',index=False)
    mapping.to_csv(out/'enrichment_gene_mapping.tsv',sep='\t',index=False)
    pd.DataFrame(audits,columns=['module','library','network_genes','background_genes','unmapped_network_genes','all_module_genes','module_genes','eligible_terms','status']).to_csv(out/'enrichment_backgrounds.tsv',sep='\t',index=False)
    summary=dict(status='completed' if len(d) else 'not_estimable',reasons=[] if len(d) else sorted({a['status'] for a in audits}) if audits else ['no_assigned_modules' if not modules else 'no_supplied_libraries'],tests=len(d),tested_modules=int(d.module.nunique()),discovered_modules=len(modules),background_genes=len(background),network_events=len(assignment),assigned_events=int(assignment.module.gt(0).sum()),gene_set_collection=', '.join(libraries),supplied_libraries=[dict(library=lib,categories=sorted({t['category'] for t in terms if t['library']==lib})) for lib in libraries],background_rule='unique unambiguous genes from every selected network event, including unassigned events, intersected with each library annotation',foreground_rule='all unique unambiguous module genes intersected with the same library universe; not Top 10',fdr_family='all discovered modules x all eligible terms across supplied libraries',term_size_rule='10 to 500 genes after intersection with library-mapped network background',mapping_rule='exact unambiguous gene IDs and gene symbols from event annotation; no fuzzy or ambiguous aliases',libraries=audits,significant_counts=d[d.q_value.le(.05)].groupby(['module','category']).size().rename('terms').reset_index().to_dict('records'),interpretation='Functional annotation relative to selected network genes. Not causal regulation, splice-specific mechanism, or independent biological validation; selection and gene event multiplicity can still bias enrichment.')
    summary['significant_modules']=int(d[d.q_value.le(.05)].module.nunique())
    summary['significant_pathways']=int(d[d.q_value.le(.05)].pathway.nunique())
    (out/'enrichment_manifest.json').write_text(json.dumps(summary,indent=2),encoding='utf8')
    return summary

def enrich_single(assignment,annotation,path,out):
    return calculate(assignment,annotation,read_gmt(path,Path(path).stem,Path(path).stem),out)

def enrich_libraries(assignment,annotation,libraries,out):
    terms=[];seen=set()
    for item in libraries:
        path=Path(item['path']);library=item.get('library') or path.stem
        if library in seen:raise ValueError(f'Duplicate enrichment library identity: {library}')
        seen.add(library)
        terms.extend(read_gmt(path,library,item.get('category') or library))
    return calculate(assignment,annotation,terms,out)

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
