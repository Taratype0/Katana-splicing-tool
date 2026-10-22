import { useRef, useState, type Ref } from "react";
import { keepPreviousData, useQuery } from "@tanstack/react-query";
import { api, type HubEvent, type ModuleDashboard as Dashboard } from "../api";
import { moduleFileStem, saveBlob, saveNetwork, tableTSV } from "./moduleDownloads";

const colors = ["#287a98", "#b85362", "#3a8c78", "#a88337", "#646c9f"];
function zColor(value: number|null) {
  if (value == null) return "#d8dcde";
  const t = Math.min(Math.abs(value) / 2, 1);
  const end = value < 0 ? [122,171,205] : [221,151,158];
  return `rgb(${end.map(c => Math.round(250+(c-250)*t)).join(",")})`;
}

function Network({data, select, svgRef}: {data: Dashboard; select: (event: string)=>void; svgRef: Ref<SVGSVGElement>}) {
  const points = new Map(data.nodes.map(n => [n.event_key, n]));
  const x = (v: number)=>190+v*145;
  const y = (v: number)=>170+v*130;
  const max = Math.max(1, ...data.nodes.map(n=>n.connectivity));
  const selected = data.nodes.find(n=>n.event_key===data.selected_event);
  const selectedLabel = selected ? `${selected.gene || "Unannotated"} | ${selected.event_key}` : "";
  return <svg ref={svgRef} viewBox="0 0 380 340" role="img" aria-label="Splice-event association network">
    <text x="10" y="15" fontSize="9" fill="#243f4b"><title>{selectedLabel}</title>Selected: {selectedLabel.length>60?selectedLabel.slice(0,57)+"...":selectedLabel}</text>
    {data.edges.map(e=>{const a=points.get(e.event_a)!, b=points.get(e.event_b)!; return <line key={`${e.event_a}:${e.event_b}`} x1={x(a.x)} y1={y(a.y)} x2={x(b.x)} y2={y(b.y)} stroke={e.same_gene?"#bf9095":"#9bb6bd"} strokeOpacity=".5" strokeWidth={.5+e.weight*2}><title>{e.same_gene?"Same gene; possible structural coupling. ":""}Association weight: {e.weight.toFixed(3)}</title></line>})}
    {data.nodes.map(n=><g key={n.event_key} onClick={()=>select(n.event_key)} tabIndex={0} role="button" aria-label={`Select event ${n.event_key}`} onKeyDown={e=>{if(e.key==="Enter")select(n.event_key)}}>
      <circle cx={x(n.x)} cy={y(n.y)} r={3+6*Math.sqrt(n.connectivity/max)} fill={n.event_key===data.selected_event?"#b85362":"#287a98"} stroke="white"><title>{n.gene||"Unannotated gene"}: {n.event_key}; weighted connectivity {n.connectivity.toFixed(2)}</title></circle>
    </g>)}
  </svg>;
}

function Samples({data, scores=false}: {data: Dashboard; scores?: boolean}) {
  const rows = scores ? data.scores.map(r=>({...r, psi:r.value})) : data.event_profile;
  const levels = data.conditions.filter(c=>rows.some(r=>r.condition===c));
  const values = rows.flatMap(r=>r.psi==null?[]:[r.psi]);
  const low = scores ? Math.min(0,...values) : 0, high = scores ? Math.max(1e-3,...values) : 1;
  const y = (v:number)=>205-(v-low)/(high-low)*160;
  return <svg viewBox="0 0 390 260" role="img" aria-label={scores?"Module score by condition":"Event PSI by condition"}>
    <line x1="42" y1="45" x2="42" y2="205" stroke="#93a7ad"/><line x1="42" y1="205" x2="380" y2="205" stroke="#93a7ad"/>
    {[low,(low+high)/2,high].map(v=><text key={v} x="37" y={y(v)+4} textAnchor="end" fontSize="10">{v.toFixed(2)}</text>)}
    {levels.map((c,j)=>{const group=rows.filter(r=>r.condition===c); const x=65+(j+.5)*300/levels.length; return <g key={c}>
      {group.map((r,i)=>r.psi==null?null:<circle key={r.sample_id} cx={x+((i*17)%31-15)} cy={y(r.psi)} r="3" fill={colors[data.conditions.indexOf(c)%colors.length]} fillOpacity=".55"><title>{r.sample_id}: {r.psi.toFixed(3)}</title></circle>)}
      <text x={x} y="229" textAnchor="middle" fontSize="10">{c}</text><text x={x} y="245" textAnchor="middle" fontSize="9">{group.length} measurements</text>
    </g>})}
  </svg>;
}

export function ModuleDashboard({projectId,artifactId,module}: {projectId:string;artifactId:string;module:string}) {
  const networkRef = useRef<SVGSVGElement>(null);
  const [exportStatus,setExportStatus] = useState("");
  const [metricSort,setMetricSort] = useState<keyof HubEvent>("rank_all_edges");
  const [descending,setDescending] = useState(false);
  const [condition,setCondition] = useState("all"), [event,setEvent] = useState<string>(), [hide,setHide] = useState(false);
  const query = useQuery({queryKey:["module-dashboard",projectId,artifactId,module,condition,event,hide],queryFn:()=>api.moduleDashboard(projectId,artifactId,module,condition,event,hide),placeholderData:keepPreviousData,retry:false});
  const d=query.data;
  async function exportFigure(format: "svg" | "png") {
    if (!d || !networkRef.current) return;
    try { await saveNetwork(networkRef.current, d, format, hide); setExportStatus(`${format.toUpperCase()} network exported: ${d.module} / ${d.condition}`); }
    catch (error) { setExportStatus(error instanceof Error ? error.message : "Export failed"); }
  }
  function exportTable(suffix: string, rows: Array<object>, columns: string[]) {
    if (!d) return;
    const scoped = rows.map(row=>({module:d.module,condition:d.condition,discovery_run:d.discovery_run,hide_same_gene_in_drawing:hide,...row}));
    saveBlob(new Blob([tableTSV(scoped,["module","condition","discovery_run","hide_same_gene_in_drawing",...columns])],{type:"text/tab-separated-values;charset=utf-8"}), `${moduleFileStem(d)}-${hide&&suffix.startsWith("displayed-")?"cross-gene-view-":""}${suffix}.tsv`);
    setExportStatus(`Exported ${rows.length} rows: ${d.module} / ${d.condition}`);
  }
  const hubColumns = ["event_key","gene","kWithin","kWithin_no_same_gene","kME","rank_all_edges","rank_no_same_gene","rank_change_no_same_gene","preset_selection_frequency_all_edges","preset_selection_frequency_no_same_gene","bootstrap_frequency_all_edges","bootstrap_frequency_no_same_gene","hub_candidate","label_status","hub_label_status","estimation_status","estimation_reason","selection_reason"];
  const fullMetrics = [...(d?.hub_comparison || [])].sort((a,b)=>{
    const av=a[metricSort], bv=b[metricSort];
    if (av==null) return bv==null ? a.event_key.localeCompare(b.event_key) : 1;
    if (bv==null) return -1;
    const order=typeof av==="number" && typeof bv==="number" ? av-bv : String(av).localeCompare(String(bv));
    return (descending?-order:order) || a.event_key.localeCompare(b.event_key);
  });
  function metricHeading(key: keyof HubEvent, label: string) {
    return <th aria-sort={metricSort===key ? descending?"descending":"ascending" : "none"}><button type="button" onClick={()=>{
      setDescending(metricSort===key?!descending:!key.startsWith("rank")); setMetricSort(key);
    }}>{label}</button></th>;
  }
  return <section className="module-dashboard" aria-label="Linked event module views">
    <h4>Event Network &amp; Module Activity</h4>
    <p>One node is one splice event. Click a node to inspect its PSI measurements.</p>
    {query.isLoading && <p>Loading the completed discovery run...</p>}
    {query.error && <p role="alert">{query.error.message}</p>}
    {query.isFetching && !query.isLoading && <p role="status">Updating views. Currently showing: {d?.condition}.</p>}
    {d && <>
      <div className="module-controls"><label>Network condition<select value={condition} onChange={e=>setCondition(e.target.value)}><option value="all">All samples (discovery network)</option>{d.conditions.map(c=><option key={c}>{c}</option>)}</select></label>
        <label><input type="checkbox" checked={hide} onChange={e=>setHide(e.target.checked)}/> Hide same-gene edges in drawing</label></div>
      <div className="module-counts"><strong>{d.total_events} events</strong><span>{d.biological_samples} biological samples</span><span>{d.total_edges} association edges</span><span>{d.same_gene_edges} same-gene edges</span></div>
      {d.status!=="available" && <p className="module-notice">{d.method}. Measurements remain available below.</p>}
      <div className="module-downloads" role="group" aria-label="Export current module results">
        <button disabled={query.isFetching} onClick={()=>exportFigure("svg")}>Network SVG</button>
        <button disabled={query.isFetching} onClick={()=>exportFigure("png")}>Network PNG</button>
        <button disabled={query.isFetching} onClick={()=>exportTable("displayed-nodes",d.nodes,["event_key","gene","connectivity","x","y"])}>Displayed nodes TSV</button>
        <button disabled={query.isFetching} onClick={()=>exportTable("displayed-edges",d.edges,["event_a","event_b","weight","correlation","pair_n","same_gene"])}>Displayed edges TSV</button>
        <button disabled={query.isFetching || !d.hub_comparison?.length} onClick={()=>exportTable("hub-all-events",d.hub_comparison || [],hubColumns)}>All hub metrics TSV</button>
        <button disabled={query.isFetching || !d.hub_events.length} onClick={()=>exportTable("hub-candidates",d.hub_events,hubColumns)}>Optional labelled hub candidates TSV</button>
      </div>
      {exportStatus && <p role="status">{exportStatus}</p>}
      <p>Exports follow the displayed module and sample scope. Drawing tables contain only shown nodes/edges; hub metrics include all assessed module events.</p>
      <div className="module-view-grid"><article><h5>Event Association Network</h5><Network data={d} select={setEvent} svgRef={networkRef}/><p>{d.method}</p></article>
        <article><h5>Sample Heatmap: Scaled PSI</h5><div className="module-heat-scroll"><svg viewBox={`0 0 ${Math.max(360,d.heatmap.samples.length*7)} ${d.heatmap.events.length*7+22}`} role="img" aria-label="Row-scaled PSI heatmap">
          {d.heatmap.z.flatMap((row,i)=>row.map((v,j)=><rect key={`${i}:${j}`} x={j*7} y={i*7+15} width="7" height="7" fill={zColor(v)}><title>{d.heatmap.events[i]} / {d.heatmap.samples[j]}: {v==null?"Missing":`z=${v.toFixed(2)}`}</title></rect>))}
          {d.heatmap.conditions.map((c,j)=><rect key={j} x={j*7} y="0" width="7" height="9" fill={colors[d.conditions.indexOf(c)%colors.length]}><title>{c}</title></rect>)}
        </svg></div><p>{d.heatmap.scale}. Shown: {d.nodes.length}/{d.total_events} events and {d.heatmap.samples.length} samples.</p></article>
        <article><h5>Module Score by Condition</h5>{d.scores.length?<Samples data={d} scores/>:<p>No cached module scores in this run.</p>}<p>{d.score_definition}</p></article>
        <article><h5>Selected Event: Inclusion by Condition</h5><select aria-label="Selected splice event" value={d.selected_event} onChange={e=>setEvent(e.target.value)}>{d.event_options.map(e=><option key={e.event_key} value={e.event_key}>{e.gene?`${e.gene} | `:""}{e.event_key}</option>)}</select><Samples data={d}/><p>PSI is the fraction of transcripts including this event (0-1). Each dot is a measurement, not an independent replication claim.</p></article></div>
      <h5>Hub metrics and optional labels</h5><p>{d.hub_method || "Run the explicit hub-metric analysis on this completed discovery result. Drawing connectivity does not establish a hub."}</p>
      {d.hub_label_status==="not_applied" || d.hub_status==="labels_not_applied" ? <p role="status">Optional threshold labels were not applied. All continuous metrics and complete rankings remain available; unlabelled events have not failed a rule.</p>
        : !d.hub_events.length && <p role="status">No optional labelled hub candidates for this scope. Status: {d.hub_status || "not_computed"}. Event measurements remain available.</p>}
      <p>Label status: {d.hub_label_status || "not recorded"}. Preset-selection frequency is descriptive sensitivity to an optional rule, not a probability, P value or FDR.</p>
      {!!d.hub_events.length && <table className="module-hubs"><thead><tr><th>Gene label</th><th>Event identifier</th><th>Full kWithin</th><th>kME</th><th>Bootstrap selection</th><th>Rank without same-gene edges</th></tr></thead><tbody>{d.hub_events.map(h=><tr key={h.event_key}><td>{h.gene||"Not annotated"}</td><td><button onClick={()=>setEvent(h.event_key)}>{h.event_key}</button></td><td>{(h.kWithin??h.connectivity).toFixed(2)}</td><td>{h.kME?.toFixed(3)??"Not computed"}</td><td>{h.bootstrap_frequency_all_edges?.toFixed(2)??"Not computed"} / {h.bootstrap_frequency_no_same_gene?.toFixed(2)??"Not computed"} (all / cross-gene)</td><td>{h.rank_no_same_gene??"Not computed"}</td></tr>)}</tbody></table>}
      {!!d.hub_comparison?.length && <section aria-label="Complete hub metrics"><h5>All events: same-gene sensitivity and complete rankings ({d.hub_comparison.length})</h5><p>Click a column heading to sort every assessed event. Full connectivity uses all estimable pairs, without the drawing threshold or edge cap. No top-rank display cutoff is applied. Unavailable measurements and optional label status are kept explicitly.</p>
        <table><thead><tr>{metricHeading("event_key","Event")}{metricHeading("rank_all_edges","Full rank")}{metricHeading("rank_no_same_gene","Cross-gene rank")}{metricHeading("rank_change_no_same_gene","Rank change")}{metricHeading("kWithin","Full kWithin")}{metricHeading("kWithin_no_same_gene","Cross-gene kWithin")}{metricHeading("kME","kME")}{metricHeading("preset_selection_frequency_all_edges","Preset-selection frequency: all / cross-gene")}<th>Label / estimation status</th></tr></thead><tbody>{fullMetrics.map(h=><tr key={h.event_key}><td>{h.gene} | {h.event_key}</td><td>{h.rank_all_edges??"Unavailable"}</td><td>{h.rank_no_same_gene??"Unavailable"}</td><td>{h.rank_change_no_same_gene??"Unavailable"}</td><td>{h.kWithin?.toFixed(3)??"Unavailable"}</td><td>{h.kWithin_no_same_gene?.toFixed(3)??"Unavailable"}</td><td>{h.kME?.toFixed(3)??"Unavailable"}</td><td>{(h.preset_selection_frequency_all_edges??h.bootstrap_frequency_all_edges)?.toFixed(2)??"Unavailable"} / {(h.preset_selection_frequency_no_same_gene??h.bootstrap_frequency_no_same_gene)?.toFixed(2)??"Unavailable"}</td><td>{h.label_status??h.hub_label_status??"not recorded"}; {h.estimation_status??h.selection_reason??"not recorded"}</td></tr>)}</tbody></table></section>}
      <h5>Module Functions &amp; Contributing Genes</h5><p>{d.enrichment_scope}</p>
      {d.enrichment.length ? <><svg viewBox={`0 0 720 ${d.enrichment.length*28+28}`} role="img" aria-label="Module functional enrichment bubble plot">
        {d.enrichment.map((e,i)=><g key={`${e.term}:${i}`}><text x="4" y={24+i*28} fontSize="11">{e.term.length>65?e.term.slice(0,62)+"...":e.term}</text><circle cx={540+Math.min(8,-Math.log10(Math.max(e.FDR??1,1e-300)))*18} cy={20+i*28} r={3+Math.sqrt(e.overlap_count)} fill={(e.FDR??1)<=.05?"#b85362":"#83abc3"}><title>{e.term}; {e.overlap_count} genes; FDR {e.FDR}; {e.overlap_genes}</title></circle><text x="694" y={24+i*28} fontSize="10">{e.overlap_count}</text></g>)}
      </svg><p>Bubble size: contributing genes. Horizontal position: -log10(FDR). Red: FDR≤0.05; blue: other terms. Shown: {d.enrichment.length}/{d.enrichment_total} terms.</p><details><summary>Contributing gene lists and statistics</summary>{d.enrichment.map((e,i)=><p key={i}><strong>{e.term}</strong>: {e.overlap_genes||"No contributing genes"}; P={e.pvalue?.toPrecision(3)??"missing"}, FDR={e.FDR?.toPrecision(3)??"missing"}</p>)}</details></> : <p>No cached functional-enrichment result for this module. Configure a species-matched pathway resource and run module interpretation.</p>}
      <details><summary>How this view was built</summary><p>{d.display_rule}</p><p>{d.interpretation}</p><p>Discovery run: {d.discovery_run}</p>{d.sources.map(s=><a key={s.artifact_id} href={api.artifactExportUrl(projectId,s.artifact_id,"tsv")}>{s.name}</a>)}</details>
    </>}
  </section>;
}
