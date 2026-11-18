import {useState,useRef} from 'react';
import {useQuery} from '@tanstack/react-query';
import {api} from '../api';
import {makeTargetEntry} from '../targetBasket';
import {PopulationPlot, type PlotRange} from './PopulationPlot';
import {SelectedEvidence} from './SelectedEvidence';
import {LocalEnrichment} from './LocalEnrichment';
import type {Artifact, TargetBasketEntry, TablePage} from '../types';

type ResultPage=TablePage & {source_rows?:number;scopes?:Record<string,string[]>};
export function ReportExplorer({projectId, artifacts, onAdd, imported=true, eventType="", initialSearch="", initialView="plot"}: {projectId:string; artifacts:Artifact[]; onAdd:(entry:TargetBasketEntry)=>void;imported?:boolean;eventType?:string;initialSearch?:string;initialView?:"plot"|"table"|"enrichment"}) {
  const priority=(a:Artifact)=>{
    const type=String(a.metadata.source_analysis_id||a.analysis_id);
    const subset=/significant|top.?10|volcano|plot.data/i.test(a.name)?5:0;
    return (type==='differential.deg'?0:type==='differential.as'?1:type==='single_cell.differential'?2:10)+subset;
  };
  const tables=artifacts.filter(a=>a.analysis_id==='imported.report'||(a.kind==='table'&&a.path.endsWith('.parquet'))).sort((a,b)=>priority(a)-priority(b)||Number(b.row_count||0)-Number(a.row_count||0));
  const [id,setId]=useState('');const artifact=tables.find(a=>a.artifact_id===id)||tables[0];
  // A fresh table gets independent query state; selections already added to the basket persist.
  return artifact?<ResultExplorer key={`${projectId}:${artifact.artifact_id}`} projectId={projectId} artifact={artifact} tables={tables} onTable={setId} onAdd={onAdd} imported={imported} eventType={eventType} initialSearch={initialSearch} initialView={initialView}/>:<p>No complete result tables are available. Import a report or complete a project analysis.</p>;
}

function ResultExplorer({projectId,artifact,tables,onTable,onAdd,imported,eventType,initialSearch,initialView}:{projectId:string;artifact:Artifact;tables:Artifact[];onTable:(id:string)=>void;onAdd:(entry:TargetBasketEntry)=>void;imported:boolean;eventType:string;initialSearch:string;initialView:"plot"|"table"|"enrichment"}) {
  const [search,setSearch]=useState(initialSearch);const [identifiers,setIdentifiers]=useState('');const [comparison,setComparison]=useState('');
  const [dasMode,setDasMode]=useState<"ps"|"ordinary"|"legacy10">("ps");
  const readSupport=dasMode!=="ordinary";
  const [fdr,setFdr]=useState<string|null>(null);const [effect,setEffect]=useState<string|null>(null);
  const [fdrMax,setFdrMax]=useState('0.05');const [effectMin,setEffectMin]=useState(/differential.as|splicing/i.test(String(artifact.metadata.source_analysis_id))?'0.1':'1');const [direction,setDirection]=useState('both');
  const [scopes,setScopes]=useState<Record<string,string>>(eventType?{event_type:eventType}:{});const [range,setRange]=useState<PlotRange>();
  const [offset,setOffset]=useState(0);const [sort,setSort]=useState('');const [descending,setDescending]=useState(false);
  const [selected,setSelected]=useState<Record<string,Record<string,unknown>>>({});const [error,setError]=useState('');const [busy,setBusy]=useState('');
  const [view,setView]=useState<'plot'|'table'|'enrichment'>(initialView);const [showColumns,setShowColumns]=useState(false);const [compact,setCompact]=useState<{archive_url:string;directory:string}>();
  const schema=useQuery({queryKey:['report-schema',projectId,artifact.artifact_id],queryFn:()=>api.browseReport(projectId,artifact.artifact_id,{limit:1})});
  const columns=schema.data?.columns||[];
  const fdrColumn=fdr??['deg_fdr','padj','fdr','FDR','qvalue','q_value','adjusted_pvalue'].find(c=>columns.includes(c))??'';
  const effectColumn=effect??['log2FoldChange','log2fc','log2FC','dpsi','delta_psi','IncLevelDifference','effect'].find(c=>columns.includes(c))??'';
  const settings={search,identifiers:identifiers.split(/[\s,;]+/).filter(Boolean),comparison,scopes:{...scopes,...(dasMode==="ps"&&columns.includes("passes_read_support")?{passes_read_support:"true"}:dasMode==="legacy10"&&columns.includes("passes_legacy_all10")?{passes_legacy_all10:"true"}:{})},
    fdr_column:fdrColumn,fdr_max:!fdrColumn||fdrMax===''?null:Number(fdrMax),effect_column:effectColumn,effect_min:!effectColumn||effectMin===''?null:Number(effectMin),direction,...range};
  const data=useQuery({queryKey:['report-table',projectId,artifact.artifact_id,settings,offset,sort,descending],
    queryFn:()=>api.browseReport(projectId,artifact.artifact_id,{...settings,offset,limit:100,sort_column:sort,sort_desc:descending}),enabled:!!schema.data});
  const plotSettings={comparison,scopes,search,identifiers:settings.identifiers,fdr_column:fdrColumn,effect_column:effectColumn};
  const population=useQuery({queryKey:['report-population',projectId,artifact.artifact_id,plotSettings],queryFn:()=>api.plotReport(projectId,artifact.artifact_id,plotSettings),enabled:!!fdrColumn&&!!effectColumn});
  const hasPlot=!!fdrColumn&&!!effectColumn;const activeView=!hasPlot&&view==='plot'?'table':view;
  const effectCutoff=effectMin===''?(effectColumn.toLowerCase().includes('psi')||effectColumn==='IncLevelDifference'?.1:1):Number(effectMin);
  const fdrCutoff=fdrMax===''?.05:Number(fdrMax);
  const result=data.data as ResultPage|undefined;const rows=data.isFetching?[]:(result?.rows||[]);
  const scopeOptions=(schema.data as ResultPage|undefined)?.scopes||{};
  const priority=['gene_name','deg_fdr','log2fc','gene_symbol','gene_id','event_id','event_type','comparison_id','cell_type','cluster_id','donor_id','fdr','FDR','padj','qvalue','dpsi','IncLevelDifference','log2FoldChange'];
  const primary=['gene_name','gene_id','event_id','event_type','comparison_id',effectColumn,fdrColumn,'baseMean','cell_type'];
  const visibleColumns=showColumns||!primary.some(c=>columns.includes(c))?[...priority.filter(c=>columns.includes(c)),...columns.filter(c=>!priority.includes(c)&&!['row_id','source_row_id'].includes(c))]:[...new Set(primary.filter(c=>columns.includes(c)))];
  const selectionEpoch=useRef(0);
  const reset=()=>{selectionEpoch.current++;setOffset(0);setSelected({});setError('');};
  const add=(kind:'gene'|'event')=>{try{onAdd(makeTargetEntry(artifact,Object.values(selected),kind,settings));setError('');}catch(e){setError(String(e));}};
  async function all(){const epoch=selectionEpoch.current;setBusy('selection');setError('');try{const allRows=await api.selectReport(projectId,artifact.artifact_id,settings);if(epoch===selectionEpoch.current)setSelected(Object.fromEntries(allRows.rows.map(r=>[String(r.row_id),r])));}catch(e){setError(String(e));}finally{setBusy('');}}
  async function point(rowId:string){const epoch=selectionEpoch.current;try{const response=await api.browseReport(projectId,artifact.artifact_id,{row_ids:[rowId],limit:1});if(epoch===selectionEpoch.current)setSelected(s=>({...s,...Object.fromEntries(response.rows.map(r=>[String(r.row_id),r]))}));}catch(e){setError(String(e));}}
  async function download(){setBusy('export');try{await api.exportReportRows(projectId,artifact.artifact_id,settings);}catch(e){setError(String(e));}finally{setBusy('');}}
  async function compactReport(){setBusy('compact');setError('');try{setCompact(await api.compactReport(projectId));}catch(e){setError(String(e));}finally{setBusy('');}}
  const clearFilters=()=>{setDasMode("ordinary");setSearch('');setIdentifiers('');setFdrMax('');setEffectMin('');setDirection('both');setRange(undefined);setScopes({});setComparison('');reset();};
  return <section className="result-card report-explorer" aria-label="Report explorer">
    <div className="artifact-nav report-toolbar"><label>Result<select aria-label="Result table" value={artifact.artifact_id} onChange={e=>onTable(e.target.value)}>{tables.map(a=><option key={a.artifact_id} value={a.artifact_id}>{a.name}{a.contrast_ids?.length?` · ${a.contrast_ids.join(', ')}`:''} · {(a.row_count||0).toLocaleString()} rows</option>)}</select></label>
      {imported&&<><a href={`/api/projects/${encodeURIComponent(projectId)}/report-source/REPORT.html`} target="_blank" rel="noreferrer">Original report ↗</a><button disabled={!!busy} onClick={()=>void compactReport()}>{busy==='compact'?'Compacting…':'Export compact report'}</button></>}
    </div>
    {compact&&<p>New report: {compact.directory} · <a href={compact.archive_url}>Download deduplicated REPORT.zip</a></p>}
    <details className="report-filter-panel" open><summary>Filters and target selection</summary><div className="report-filters">
      <label className="search-field"><input aria-label="Search" placeholder="Search results" value={search} onChange={e=>{setSearch(e.target.value);reset();}}/></label>
      {(scopeOptions.comparison_id||[]).length>1&&<label>Comparison<select aria-label="Comparison ID" value={comparison} onChange={e=>{setComparison(e.target.value);reset();}}><option value="">All comparisons</option>{(scopeOptions.comparison_id||[]).map(c=><option key={c}>{c}</option>)}</select></label>}
      <details className="identifier-filter"><summary>Gene / event list{settings.identifiers.length?` · ${settings.identifiers.length}`:''}</summary><textarea aria-label="Gene IDs, symbols or event IDs" value={identifiers} onChange={e=>{setIdentifiers(e.target.value);reset();}} placeholder="Paste exact IDs or gene symbols"/></details>
      <label>FDR ≤<input aria-label="Maximum FDR" placeholder="All results" type="number" min="0" max="1" step="0.01" disabled={!fdrColumn} value={fdrMax} onChange={e=>{setFdrMax(e.target.value);setRange(undefined);reset();}}/></label>
      <label>{effectColumn.toLowerCase().includes('psi')?'|ΔPSI| ≥':'|log2 fold change| ≥'}<input aria-label="Minimum absolute effect" placeholder="All effects" type="number" min="0" step="0.01" disabled={!effectColumn} value={effectMin} onChange={e=>{setEffectMin(e.target.value);setRange(undefined);reset();}}/></label>
      {columns.includes("passes_read_support")&&<label>DAS event set<select aria-label="DAS event set" value={dasMode} onChange={e=>{setDasMode(e.target.value as "ps"|"ordinary"|"legacy10");reset();}}><option value="ps">PS-DAS · precision & quantile support (default)</option><option value="ordinary">Ordinary DAS · unfiltered native statistics</option>{columns.includes("passes_legacy_all10")&&<option value="legacy10">Legacy All-Sample ≥10 (strict benchmark)</option>}</select></label>}
      <label>Direction<select aria-label="Effect direction" disabled={!effectColumn} value={direction} onChange={e=>{setDirection(e.target.value);reset();}}><option value="both">Both</option><option value="positive">Positive</option><option value="negative">Negative</option></select></label>
      {Object.entries(scopeOptions).filter(([c])=>!['comparison_id','passes_read_support','passes_precision_support','passes_legacy_all10','nmd_asymmetry_rescued'].includes(c)).map(([c,options])=><label key={c}>{c}<select aria-label={`Scope ${c}`} value={scopes[c]||''} onChange={e=>{setScopes(s=>({...s,[c]:e.target.value}));reset();}}><option value="">All</option>{options.map(v=><option key={v}>{v}</option>)}</select></label>)}
      <details className="column-bindings"><summary>Column bindings</summary><label>FDR column<select aria-label="FDR column" value={fdrColumn} onChange={e=>{setFdr(e.target.value);setFdrMax('0.05');setRange(undefined);reset();}}><option value="">Choose column</option>{columns.map(c=><option key={c}>{c}</option>)}</select></label><label>Effect column<select aria-label="Effect column" value={effectColumn} onChange={e=>{setEffect(e.target.value);setEffectMin(/psi|IncLevelDifference/i.test(e.target.value)?'0.1':'1');setRange(undefined);reset();}}><option value="">Choose column</option>{columns.map(c=><option key={c}>{c}</option>)}</select></label></details><button onClick={()=>{clearFilters();setDasMode("ps");setFdrMax("0.05");setEffectMin(effectColumn.toLowerCase().includes("psi")||effectColumn==="IncLevelDifference"?"0.1":"1");setScopes(eventType?{event_type:eventType}:{});}}>Restore defaults</button><button onClick={clearFilters}>Include nonsignificant results</button><button onClick={()=>{setFdrMax(String(fdrCutoff));setEffectMin(String(effectCutoff));reset();}}>Apply significance cutoffs</button>
    </div></details>
    {range&&<p>Plot range: {range.effect_low.toPrecision(3)} ≤ effect ≤ {range.effect_high.toPrecision(3)}; {range.fdr_min.toExponential(2)} ≤ FDR ≤ {range.fdr_max.toExponential(2)} <button onClick={()=>{setRange(undefined);reset();}}>Clear plot range</button></p>}
    {(error||data.error||schema.error||population.error)&&<p role="alert">{error||data.error?.message||schema.error?.message||population.error?.message}</p>}
    <div className="artifact-tabs report-view-tabs" role="tablist" aria-label="Result views">{(hasPlot?['plot','table','enrichment']:['table','enrichment']).map(v=><button role="tab" aria-selected={activeView===v} className={activeView===v?'active':''} key={v} onClick={()=>setView(v as typeof view)}>{v==='plot'?'Volcano plot':v==='table'?'Complete result table':'Pathway enrichment'}</button>)}</div>
    {activeView==='plot'&&(population.isFetching?<p>Updating full-population plot…</p>:population.data&&<PopulationPlot data={population.data} effect={effectColumn} fdr={fdrColumn} fdrCutoff={fdrCutoff} effectCutoff={effectCutoff} readSupport={readSupport&&columns.includes("passes_read_support")} onRange={r=>{setRange(r);reset();}} onSelect={id=>void point(id)}/>)}
    {activeView==='plot'&&<SelectedEvidence rows={Object.values(selected)} effect={effectColumn} fdr={fdrColumn}/>}
    <div className="target-actions"><button disabled={data.isFetching||!rows.length||!!busy} onClick={()=>setSelected(s=>({...s,...Object.fromEntries(rows.map(r=>[String(r.row_id),r]))}))}>Select this page</button>
      <button disabled={data.isFetching||!result?.total||!!busy} onClick={()=>void all()}>{busy==='selection'?'Selecting…':`Select all ${result?.total||0} matching rows`}</button>
      <button disabled={!Object.keys(selected).length} onClick={()=>{selectionEpoch.current++;setSelected({});}}>Clear</button>
      <button disabled={!Object.keys(selected).length||!!busy} onClick={()=>add('gene')} aria-label="Add genes to targets">Add genes</button><button disabled={!Object.keys(selected).length||!!busy} onClick={()=>add('event')} aria-label="Add events to targets">Add events</button>
      <button disabled={!!busy||!result?.total} onClick={()=>void download()}>Download all matching rows</button>
      <button onClick={()=>setView('enrichment')}>Enrich matching / selected genes</button></div>
    <div className="report-table-status"><span>{data.isFetching?'Loading…':`${result?.total??0} matching rows`} / {artifact.row_count} source rows · {Object.keys(selected).length} selected</span><details><summary>Source</summary><p>{String(artifact.metadata.source_analysis_id||artifact.analysis_id)} · {artifact.metadata.report_hash_verified?'Checksum verified':'See recorded source verification'}</p><p>All delivered rows remain available, including nonsignificant results. Browsing never changes original statistics or test families.</p></details></div>
    {activeView==='table'&&<><label className="report-column-toggle"><input type="checkbox" checked={showColumns} onChange={e=>setShowColumns(e.target.checked)}/> Show all original columns and mapping provenance</label><div className="report-table-wrap"><table><thead><tr><th>Select</th>{visibleColumns.map(c=><th key={c}><button onClick={()=>{setSort(c);setDescending(sort===c?!descending:false);setOffset(0);}}>{c==='gene_name'?'Gene name':c}{sort===c?(descending?' ↓':' ↑'):''}</button></th>)}</tr></thead><tbody>{rows.map(r=><tr key={String(r.row_id)}><td><input aria-label={`Select ${r.row_id}`} type="checkbox" checked={!!selected[String(r.row_id)]} onChange={e=>setSelected(s=>{const next={...s};if(e.target.checked)next[String(r.row_id)]=r;else delete next[String(r.row_id)];return next;})}/></td>{visibleColumns.map(c=><td key={c}>{r[c]===null||r[c]===undefined?'NA':String(r[c])}</td>)}</tr>)}</tbody></table></div>
    <div className="target-actions"><button disabled={!offset} onClick={()=>setOffset(Math.max(0,offset-100))}>Previous</button><span>{rows.length?offset+1:0}–{offset+rows.length} / {result?.total??0}</span><button disabled={offset+rows.length>=(result?.total||0)} onClick={()=>setOffset(offset+100)}>Next</button></div></>}
    {activeView==='enrichment'&&<LocalEnrichment projectId={projectId} artifactId={artifact.artifact_id} columns={columns} filters={settings} selectedIds={Object.keys(selected)} onGenes={genes=>{setIdentifiers(genes.join(' '));setFdrMax('');setEffectMin('');setDirection('both');setSearch('');setRange(undefined);reset();}}/>}
  </section>;
}
