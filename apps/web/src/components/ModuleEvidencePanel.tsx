import { useQuery } from "@tanstack/react-query";
import { useState } from "react";
import { createPortal } from "react-dom";
import { api } from "../api";
import { ModuleDashboard } from "./ModuleDashboard";
import { SPLICE_MODULE_NAME } from "../branding";

export function ModuleEvidencePanel({projectId, artifactId, module}: {projectId: string; artifactId: string; module: string}) {
  const [expanded,setExpanded] = useState(false);
  const result = useQuery({queryKey:["module-evidence",projectId,artifactId,module],queryFn:()=>api.moduleEvidence(projectId,artifactId,module)});
  return <section aria-label="Module evidence">
    <h3>Module {module}</h3>
    <button className="module-open" onClick={()=>setExpanded(true)}>Explore event network, heatmap &amp; hubs</button>
    {expanded && createPortal(<div className="module-modal-backdrop"><section role="dialog" aria-modal="true" aria-label={`Module ${module} explorer`} className="module-modal" onKeyDown={e=>{if(e.key==="Escape")setExpanded(false)}}><header><h3>{SPLICE_MODULE_NAME} / Module {module}</h3><button autoFocus onClick={()=>setExpanded(false)}>Close</button></header><ModuleDashboard key={`${artifactId}:${module}`} projectId={projectId} artifactId={artifactId} module={module}/></section></div>,document.body)}
    {result.isLoading && <p>Reading cached module evidence...</p>}
    {result.error && <p role="alert">{result.error.message}</p>}
    {result.data && <>
      <h4>Functional Evidence &amp; Contributing Genes</h4>
      <p>Read the full module-linked tables below: member events, hub confidence, condition associations, and enrichment contributions. Missing resources are not negative biological results.</p>
      {!result.data.tables.length && <p>No module-linked tables in this completed run.</p>}
      {result.data.tables.map((table)=><details key={table.artifact_id}>
        <summary>{table.name} ({table.row_count.toLocaleString()} rows)</summary>
        <p>Source run: {table.run_fingerprint.slice(0,12)}</p>
        <p>Preview: first {Math.min(table.preview_limit,table.row_count)} matching records. Full TSV retains all records.</p>
        <a href={api.artifactExportUrl(projectId,table.artifact_id,"tsv")} download>Full source TSV (all modules)</a>
        {table.rows.map((row,index)=><dl className="detail-list" key={String(row.row_id ?? index)}>{Object.entries(row).map(([key,value])=><div key={key}><dt>{key.replaceAll("_"," ")}</dt><dd>{value == null ? "Not available" : typeof value === "object" ? JSON.stringify(value) : String(value)}</dd></div>)}</dl>)}
      </details>)}
    </>}
  </section>;
}
