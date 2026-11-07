import {useState} from 'react';
import {useQuery} from '@tanstack/react-query';
import {api} from '../api';
import {makeTargetEntry} from '../targetBasket';
import type {Artifact, TargetBasketEntry} from '../types';

export function ReportExplorer({projectId, artifacts, onAdd}: {projectId:string; artifacts:Artifact[]; onAdd:(entry:TargetBasketEntry)=>void}) {
  const tables = artifacts.filter(a=>a.analysis_id==='imported.report');
  const [id,setId]=useState('');
  const artifact=tables.find(a=>a.artifact_id===id) || tables[0];
  const [search,setSearch]=useState('');
  const [identifiers,setIdentifiers]=useState('');
  const [comparison,setComparison]=useState('');
  const [fdr,setFdr]=useState('');
  const [effect,setEffect]=useState('');
  const [fdrMax,setFdrMax]=useState('');
  const [effectMin,setEffectMin]=useState('');
  const [offset,setOffset]=useState(0);
  const [selected,setSelected]=useState<Record<string,Record<string,unknown>>>({});
  const [error,setError]=useState('');
  const settings={search,identifiers:identifiers.split(/[\s,;]+/).filter(Boolean),comparison,
    fdr_column:fdr,fdr_max:fdrMax===''?null:Number(fdrMax),effect_column:effect,effect_min:effectMin===''?null:Number(effectMin)};
  const data=useQuery({queryKey:['report-table',projectId,artifact?.artifact_id,settings,offset],
    queryFn:()=>api.browseReport(projectId,artifact!.artifact_id,{...settings,offset,limit:100}),enabled:!!artifact});
  const columns=data.data?.columns || [];
  const rows=data.data?.rows || [];
  const reset=()=>{setOffset(0);setSelected({});};
  const add=(kind:'gene'|'event')=>{try {onAdd(makeTargetEntry(artifact!,Object.values(selected),kind,settings));setError('');}catch(e){setError(String(e));}};
  return <section className="target-workbench report-explorer" aria-label="Report explorer">
    <h2>Explore completed report</h2><a href={`/api/projects/${encodeURIComponent(projectId)}/report-source/REPORT.html`} target="_blank" rel="noreferrer">Open original seven-section report and figures</a><p>Browse complete result tables. Filters preserve original statistics and test families. Add exact rows to the Gene / Event workspace for a new target pipeline task.</p>
    <label>Result table <select value={artifact?.artifact_id || ''} onChange={e=>{setId(e.target.value);reset();setFdr('');setEffect('');setFdrMax('');setEffectMin('');setComparison('');}}>{tables.map(a=><option key={a.artifact_id} value={a.artifact_id}>{a.name} · {a.row_count} rows</option>)}</select></label>
    {artifact && <p>Source: {String(artifact.metadata.source_analysis_id || '')} · checksum {artifact.metadata.report_hash_verified?'verified against report':'recorded at import; original hash unavailable'}</p>}
    <div className="report-filters">
      <label>Search <input value={search} onChange={e=>{setSearch(e.target.value);reset();}}/></label>
      <label>Comparison ID <input value={comparison} onChange={e=>{setComparison(e.target.value);reset();}}/></label>
      <label>Gene IDs, symbols or event IDs <textarea value={identifiers} onChange={e=>{setIdentifiers(e.target.value);reset();}} placeholder="Paste one or many exact identifiers"/></label>
      <label>FDR column <select value={fdr} onChange={e=>{setFdr(e.target.value);setFdrMax('');reset();}}><option value="">No FDR filter</option>{columns.map(c=><option key={c}>{c}</option>)}</select><input aria-label="Maximum FDR" type="number" min="0" max="1" step="0.01" disabled={!fdr} value={fdrMax} onChange={e=>{setFdrMax(e.target.value);reset();}}/></label>
      <label>Effect column <select value={effect} onChange={e=>{setEffect(e.target.value);setEffectMin('');reset();}}><option value="">No effect filter</option>{columns.map(c=><option key={c}>{c}</option>)}</select><input aria-label="Minimum absolute effect" type="number" min="0" step="0.01" disabled={!effect} value={effectMin} onChange={e=>{setEffectMin(e.target.value);reset();}}/></label>
    </div>
    {(error || data.error) && <p role="alert">{error || data.error?.message}</p>}
    {fdr && effect && <><p>Effect versus −log10({fdr}), current page only. Click a point to select its exact row.</p><svg viewBox="0 0 640 250" role="img" aria-label="Interactive effect and FDR plot" style={{maxWidth:800,background:'#f5f7fa'}}>
      <line x1="320" y1="10" x2="320" y2="225" stroke="#aaa"/>
      {rows.filter(r=>r[fdr]!==null && r[fdr]!=='' && r[effect]!==null && r[effect]!=='' && Number.isFinite(Number(r[fdr])) && Number(r[fdr])>=0 && Number(r[fdr])<=1 && Number.isFinite(Number(r[effect]))).map(r=>{
        const max=Math.max(1,...rows.map(v=>Math.abs(Number(v[effect]))||0));const x=320+Number(r[effect])/max*290;
        const y=225-Math.min(20,-Math.log10(Math.max(1e-20,Number(r[fdr]))))*10;
        return <circle key={String(r.row_id)} cx={x} cy={y} r={selected[String(r.row_id)]?6:4} fill={Number(r[effect])>=0?'#C63F4A':'#2863A5'} onClick={()=>setSelected(s=>({...s,[String(r.row_id)]:r}))}><title>{String(r.gene_symbol||r.gene_id||r.event_id||r.row_id)}</title></circle>;
      })}<text x="280" y="245">{effect}</text></svg></>}
    <div className="target-actions"><button disabled={data.isFetching || !rows.length} onClick={()=>setSelected(s=>({...s,...Object.fromEntries(rows.map(r=>[String(r.row_id),r]))}))}>Select this page</button><button onClick={()=>setSelected({})}>Clear selection</button><button disabled={!Object.keys(selected).length} onClick={()=>add('gene')}>Add genes to targets</button><button disabled={!Object.keys(selected).length} onClick={()=>add('event')}>Add events to targets</button></div>
    <p>{data.isFetching?'Loading…':`${data.data?.total ?? 0} matching rows`} · {Object.keys(selected).length} selected across pages</p>
    <div style={{overflow:'auto',maxHeight:520}}><table><thead><tr><th>Select</th>{columns.map(c=><th key={c}>{c}</th>)}</tr></thead><tbody>{rows.map(r=><tr key={String(r.row_id)}><td><input aria-label={`Select ${r.row_id}`} type="checkbox" checked={!!selected[String(r.row_id)]} onChange={e=>setSelected(s=>{const next={...s};if(e.target.checked)next[String(r.row_id)]=r;else delete next[String(r.row_id)];return next;})}/></td>{columns.map(c=><td key={c}>{r[c]===null?'NA':String(r[c])}</td>)}</tr>)}</tbody></table></div>
    <div className="target-actions"><button disabled={!offset} onClick={()=>setOffset(Math.max(0,offset-100))}>Previous</button><span>{offset+1}–{offset+rows.length}</span><button disabled={offset+rows.length>=(data.data?.total||0)} onClick={()=>setOffset(offset+100)}>Next</button></div>
  </section>;
}
