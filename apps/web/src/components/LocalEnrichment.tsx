import {useState} from 'react';
import {api} from '../api';
import {VegaFigure} from './ReportFigures';
import {PathPicker} from './PathPicker';

export type EnrichmentResult={audit:Record<string,unknown>;columns:string[];rows:Record<string,unknown>[];files:Record<string,string>};
export function LocalEnrichment({projectId,artifactId,columns,filters,selectedIds,onGenes}:{projectId:string;artifactId:string;columns:string[];filters:Record<string,unknown>;selectedIds:string[];onGenes:(genes:string[])=>void}) {
  const [library,setLibrary]=useState('');const [species,setSpecies]=useState('');
  const [gene,setGene]=useState('gene_id');const [test,setTest]=useState('');const [scope,setScope]=useState('filtered');
  const [busy,setBusy]=useState(false);const [error,setError]=useState('');const [result,setResult]=useState<EnrichmentResult>();
  const [maximum,setMaximum]=useState(2000);const [cutoff,setCutoff]=useState('0.05');const [page,setPage]=useState(0);
  const validTest=test||['pvalue','PValue','p_value','pval','deg_fdr','padj','fdr','FDR'].find(c=>columns.includes(c))||'';
  async function run(){setBusy(true);setError('');setResult(undefined);try{setResult(await api.enrichReport(projectId,artifactId,{filters:{...filters,...(scope==='selected'?{row_ids:selectedIds}:{})},library_path:library,library_species:species,gene_column:gene,test_column:validTest,maximum_term_size:maximum}));setPage(0);}catch(e){setError(String(e));}finally{setBusy(false);}}
  const rows=(result?.rows||[]).filter(r=>cutoff===''||Number(r.fdr)<=Number(cutoff));
  return <section className="local-enrichment"><h3>Local pathway enrichment</h3>
    <p>Foreground: matching genes or selected rows. Background: all testable genes in the same comparison and cell/event scope, before browsing thresholds. Event rows count once per gene.</p>
    <PathPicker label="Pathway library GMT / TSV / CSV" value={library} onChange={setLibrary} kind="file" extensions={['.gmt','.tsv','.csv']} placeholder="Choose GO, KEGG, Reactome or a custom gene-set library"/>
    <div className="report-filters">
      <label>Library species<select aria-label="Library species" value={species} onChange={e=>setSpecies(e.target.value)}><option value="">Choose species</option><option value="human">Human</option><option value="mouse">Mouse</option><option value="other">Other / project species</option></select></label>
      <label>Identifiers<select aria-label="Enrichment gene column" value={gene} onChange={e=>setGene(e.target.value)}>{['gene_id','gene_symbol'].filter(c=>columns.includes(c)).map(c=><option key={c}>{c}</option>)}</select></label>
      <label>Original P / FDR<select aria-label="Background test column" value={validTest} onChange={e=>setTest(e.target.value)}><option value="">Choose column</option>{columns.map(c=><option key={c}>{c}</option>)}</select></label>
      <label>Foreground<select aria-label="Enrichment foreground" value={scope} onChange={e=>setScope(e.target.value)}><option value="filtered">All matching rows</option><option value="selected">Selected rows ({selectedIds.length})</option></select></label>
      <label>Max pathway size<input aria-label="Maximum pathway size" type="number" min={1} value={maximum} onChange={e=>setMaximum(Number(e.target.value))}/></label>
    </div>
    <button disabled={busy||!library||!species||!validTest||!columns.includes(gene)||(scope==='selected'&&!selectedIds.length)||maximum<1} onClick={()=>void run()}>{busy?'Calculating…':'Run local enrichment'}</button>
    {error&&<p role="alert">{error}</p>}
    {result&&<><p>{String(result.audit.mapped_selected_genes)} mapped foreground genes / {String(result.audit.annotated_background_genes)} background genes · {String(result.audit.tested_terms)} tested pathways. BH includes zero-overlap pathways.</p>
      <div className="target-actions">{Object.entries(result.files).map(([name,url])=><a key={name} href={url}>{name}</a>)}</div>
      <details><summary>Mapping, background and method audit</summary><pre>{JSON.stringify(result.audit,null,2)}</pre></details>
      <label>Display FDR ≤ <input aria-label="Enrichment display FDR" type="number" min={0} max={1} step={.01} placeholder="All terms" value={cutoff} onChange={e=>{setCutoff(e.target.value);setPage(0);}}/></label>
      <VegaFigure title="Pathway enrichment" spec={{$schema:'https://vega.github.io/schema/vega-lite/v5.json',width:'container',height:Math.max(120,Math.min(rows.length,20)*24),data:{values:rows.slice(0,20)},mark:{type:'circle',opacity:.85},encoding:{y:{field:'term',type:'nominal',sort:'-x',axis:{title:null,labelLimit:360}},x:{field:'fold_enrichment',type:'quantitative',axis:{title:'Fold enrichment'}},size:{field:'overlap_genes',type:'quantitative',legend:{title:'Hit genes'}},color:{field:'fdr',type:'quantitative',scale:{scheme:'blues',reverse:true},legend:{title:'FDR'}},tooltip:[{field:'term'},{field:'overlap_genes'},{field:'fold_enrichment'},{field:'fdr'}]}}} onPoint={r=>r.genes&&onGenes(String(r.genes).split(',').map(g=>g.trim()))}/>
      <div className="report-table-wrap"><table><thead><tr>{['Pathway','Overlap','Fold enrichment','P','FDR','Hit genes'].map(c=><th key={c}>{c}</th>)}</tr></thead><tbody>{rows.slice(page*50,page*50+50).map(r=><tr key={String(r.term)}><td>{String(r.term)}</td><td>{String(r.overlap_genes)} / {String(r.term_genes_in_background)}</td><td>{Number(r.fold_enrichment).toPrecision(3)}</td><td>{Number(r.pvalue).toExponential(2)}</td><td>{Number(r.fdr).toExponential(2)}</td><td><button disabled={!r.genes} onClick={()=>onGenes(String(r.genes).split(',').map(g=>g.trim()))}>Inspect {String(r.overlap_genes)} genes</button></td></tr>)}</tbody></table></div>
      <div className="target-actions"><button disabled={!page} onClick={()=>setPage(page-1)}>Previous pathways</button><span>{rows.length} matching pathways</span><button disabled={(page+1)*50>=rows.length} onClick={()=>setPage(page+1)}>Next pathways</button></div></>}
  </section>;
}
