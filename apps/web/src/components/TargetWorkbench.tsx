import { useEffect, useState } from "react";
import { useMutation, useQuery } from "@tanstack/react-query";
import { api } from "../api";
import type { ManifestContrast, TargetBasketEntry, TargetBatch, TargetExport, TargetRun, TargetPlanMessage } from "../types";
import "./TargetWorkbench.css";

type Props = {projectId: string; projectFingerprint: string; entries: TargetBasketEntry[];
  contrasts?: ManifestContrast[]; storageError?: string; onRemove: (id: string) => void; onClear: () => void};

function planMessage(value: TargetPlanMessage) {
  if (typeof value === "string") return value;
  const identity = value.target_id || value.resource_id;
  return `${identity ? `${identity}: ` : ""}${value.message || JSON.stringify(value)}`;
}

export function TargetWorkbench({projectId, projectFingerprint, entries, contrasts = [], storageError, onRemove, onClear}: Props) {
  const [modules, setModules] = useState(["gene_card", "event_card", "sashimi"]);
  const [exported, setExported] = useState<TargetExport>();
  const [bindings, setBindings] = useState<Record<string, string>>({});
  const [workers, setWorkers] = useState("");
  const [manifestPath, setManifestPath] = useState("");
  const [comparisonScope, setComparisonScope] = useState<"source" | "all" | "custom">("source");
  const [customComparisons, setCustomComparisons] = useState<string[]>([]);
  const runKey = `katana:target-run:v1:${projectId}:${projectFingerprint}`;
  const [runId, setRunId] = useState(() => localStorage.getItem(runKey) || "");
  const [returned, setReturned] = useState<TargetRun>();
  const signature = JSON.stringify({entries, modules, comparisonScope, customComparisons});
  useEffect(() => {setExported(undefined); setBindings({});}, [signature]);

  function acceptRun(run: TargetRun) {
    setRunId(run.run_id); setReturned(run);
    try {localStorage.setItem(runKey, run.run_id);} catch { /* The returned links remain available in this session. */ }
  }
  const exportRequest = useMutation({mutationFn: () => {
    if (comparisonScope === "custom" && !customComparisons.length) throw new Error("Choose at least one target comparison; an empty selection never means all.");
    return api.exportTargets(projectId, {
    selections: entries.map(({source_artifact_id, row_ids, target_kind, filters}) =>
      ({source_artifact_id, row_ids, target_kind, filters})), modules,
    ...(comparisonScope === "all" ? {contrast_ids: contrasts.map(contrast => contrast.contrast_id)} :
      comparisonScope === "custom" ? {contrast_ids: [...customComparisons]} : {}),
  });}, onSuccess: setExported});
  function executionRequest(): TargetBatch {
    if (!exported) throw new Error("Export and verify the target request first.");
    if (workers && (!Number.isInteger(Number(workers)) || Number(workers) < 1)) throw new Error("Workers must be a positive integer or blank for the available allocation.");
    return {request_id: exported.request_id,
      resources: Object.fromEntries(Object.entries(bindings).filter(([, value]) => value.trim()).map(([key, value]) => [key, value.trim()])),
      ...(workers ? {workers: Number(workers)} : {})};
  }
  const plan = useMutation({mutationFn: () => api.planTargets(projectId, executionRequest())});
  const run = useMutation({mutationFn: () => api.runTargets(projectId, executionRequest()), onSuccess: acceptRun});
  const imported = useMutation({mutationFn: () => api.importTargetResults(projectId, manifestPath.trim()), onSuccess: acceptRun});
  const status = useQuery({queryKey: ["target-run", projectId, runId], queryFn: () => api.targetRun(projectId, runId),
    enabled: Boolean(runId), refetchInterval: query => ["queued", "running"].includes(query.state.data?.status || "queued") ? 1000 : false});
  const current = status.data || returned;
  const active = current && ["queued", "running"].includes(current.status);
  const errors = [storageError, exportRequest.error?.message, plan.error?.message, run.error?.message,
    imported.error?.message, status.error?.message, current?.error].filter(Boolean);

  return <section className="target-workbench" aria-label="Gene and event target workspace">
    <header><span className="target-eyebrow">Selected-target analysis</span><h2>Gene / Event workspace</h2>
      <p>Choose exact genes or events in result tables, then add them to this basket. Cards reuse completed evidence. Optional Sashimi reads only the requested events from BAMs on the execution machine.</p>
      <p className="target-note">This workspace produces separate target reports and interactive data. The approved seven-section study report is unchanged.</p></header>
    {errors.map((message, i) => <p className="target-error" role="alert" key={i}>{message}</p>)}
    <section aria-label="Persistent target basket">
      <div className="target-heading"><h3>Target basket · {entries.reduce((sum, item) => sum + item.row_ids.length, 0)} source rows</h3>
        <button type="button" disabled={!entries.length} onClick={onClear}>Clear basket</button></div>
      {!entries.length ? <p>No targets selected. Open a result table, select rows and choose “Add genes to targets” or “Add events to targets”.</p> :
        entries.map(entry => <article className="target-selection" key={entry.basket_id}>
          <div><strong>{entry.target_kind === "gene" ? "Gene targets" : "Exact event targets"}</strong> · {entry.row_ids.length} source rows · {entry.source_name}
            <button type="button" aria-label={`Remove ${entry.source_name} ${entry.target_kind} targets`} onClick={() => onRemove(entry.basket_id)}>Remove</button></div>
          <details><summary>Selected identities and browsing context</summary><ul>{entry.labels.map((label, i) => <li key={`${label}:${i}`}>{label} <code>{entry.row_ids[i] || ""}</code></li>)}</ul>
            <pre>{JSON.stringify(entry.filters, null, 2)}</pre></details>
        </article>)}
      <fieldset><legend>Requested products</legend>{[["gene_card", "Gene cards"], ["event_card", "Event cards"], ["sashimi", "Sashimi read evidence"]].map(([id, label]) =>
        <label key={id}><input type="checkbox" checked={modules.includes(id)} onChange={event => setModules(current => event.target.checked ? [...current, id] : current.filter(value => value !== id))}/>{label}</label>)}</fieldset>
      <label className="target-resource">Target comparison scope<select value={comparisonScope} onChange={event => setComparisonScope(event.target.value as "source" | "all" | "custom")}>
        <option value="source">Selected source comparisons</option><option value="all" disabled={!contrasts.length}>All project comparisons</option>
        <option value="custom" disabled={!contrasts.length}>Custom comparisons</option></select></label>
      {comparisonScope === "source" && <p>By default, the batch uses the union of comparison IDs recorded in the exact selected source rows; if those rows have no comparison column, it uses the declared project comparisons. Browsing thresholds and comparison controls are saved separately and do not silently expand this scope.</p>}
      {comparisonScope === "custom" && <fieldset><legend>Include target comparisons</legend>{contrasts.map(contrast =>
        <label key={contrast.contrast_id}><input type="checkbox" checked={customComparisons.includes(contrast.contrast_id)} onChange={event => setCustomComparisons(current => event.target.checked ? [...current, contrast.contrast_id] : current.filter(id => id !== contrast.contrast_id))}/>{contrast.label}</label>)}</fieldset>}
      {comparisonScope === "custom" && !customComparisons.length && <p role="status">Choose at least one target comparison. An empty selection never means all.</p>}
      {comparisonScope !== "source" && <p>The requested comparisons must include those represented by the selected source rows. Exact source identities are checked during export; mismatches are rejected.</p>}
      <button className="target-primary" type="button" disabled={!entries.length || !modules.length || exportRequest.isPending || (comparisonScope === "custom" && !customComparisons.length)} onClick={() => exportRequest.mutate()}>
        {exportRequest.isPending ? "Verifying and exporting…" : "Export portable target request"}</button>
      {exported && <div className="target-export" role="status"><strong>{exported.target_count} verified targets</strong> · request <code>{exported.request_id}</code>
        <a href={exported.download_url}>Download request ZIP</a>
        <ul>{exported.targets.map((target, index) => <li key={String(target.target_id || index)}>
          {String(target.target_kind === "gene" ? target.gene_symbol || target.gene_id || target.target_id : target.event_id || target.target_id)}:
          {" "}{Array.isArray(target.comparison_ids) ? target.comparison_ids.map(String).join(", ") : "Comparison scope not recorded"}</li>)}</ul>
        <p>Run this request locally or on a server with the same Katana package. BAMs and large references are excluded.</p>
        <pre>katana-targets plan --request /path/to/request.zip --bindings /path/to/bindings.json{"\n"}katana-targets run --request /path/to/request.zip --outdir /path/to/target-results --bindings /path/to/bindings.json</pre>
      </div>}
    </section>
    <section aria-label="Target execution"><h3>Run on this service host</h3>
      <p>Paths below belong to the machine running the Katana Python service. This starts ordinary Python work on that host; it does not submit a Slurm job. Use the CLI within your scheduler allocation for cluster execution.</p>
      {exported?.resources.map(resource => <label className="target-resource" key={resource.resource_id}>
        <span>{resource.kind} · {resource.sample_id || resource.resource_id}</span><code>{resource.resource_id}</code>
        <input aria-label={`Resource path ${resource.resource_id}`} value={bindings[resource.resource_id] || ""} placeholder="Path on this service host" onChange={event => setBindings(current => ({...current, [resource.resource_id]: event.target.value}))}/></label>)}
      <label className="target-resource">Workers <input type="number" min={1} step={1} value={workers} placeholder="Available allocation" onChange={event => setWorkers(event.target.value)}/></label>
      <div className="target-actions"><button type="button" disabled={!exported || plan.isPending} onClick={() => plan.mutate()}>Check target readiness</button>
        <button type="button" className="target-primary" disabled={!exported || run.isPending || Boolean(active)} onClick={() => run.mutate()}>Run requested targets</button></div>
      {plan.data && <div role="status"><strong>{plan.data.ready ? "Ready" : "Resource or input attention required"}</strong>
        {[...plan.data.errors, ...plan.data.warnings].map((message, i) => <p key={i}>{planMessage(message)}</p>)}
        <details><summary>Tasks and reuse plan</summary><pre>{JSON.stringify(plan.data, null, 2)}</pre></details></div>}
      {current && <div className="target-run" role="status"><h4>Target run: {current.status}</h4><p>{current.execution}</p>
        {current.progress && <pre>{JSON.stringify(current.progress, null, 2)}</pre>}
        {!current.result && current.target_statuses && <ul>{Object.entries(current.target_statuses).map(([id, item]) =>
          <li key={id}><code>{id}</code>: {item.status || item.type || "observed"}{item.type === "target_reused" ? " · verified cache reused" : ""}</li>)}</ul>}
        <div className="target-actions">{current.report_url && <a href={current.report_url} target="_blank" rel="noreferrer">Open target index</a>}
          {current.archive_url && <a href={current.archive_url}>Download target result ZIP</a>}
          {current.manifest_url && <a href={current.manifest_url}>Verified result manifest</a>}</div>
        {current.result?.targets.map(target => <article className="target-result" key={target.target_id}>
          <strong>{target.target_kind === "gene" ? target.gene_symbol || target.gene_id || target.target_id : target.event_id || target.target_id}</strong> · {target.target_kind || "target"} · {target.status}
          {(target.error || target.reason) && <p>{target.error || target.reason}</p>}
          {target.sashimi?.filter(item => !["complete", "not_requested"].includes(item.status)).map((item, i) => <p key={i}>Sashimi {item.status}: {item.reason || "See target diagnostics"}</p>)}
          {target.report && <a href={api.targetFileUrl(projectId, current.run_id, target.report)} target="_blank" rel="noreferrer">Open report and data</a>}</article>)}
      </div>}
    </section>
    <section aria-label="Import completed target results"><h3>Reopen returned results</h3>
      <p>Extract a completed target result ZIP on this service host, then enter its result manifest or directory. Import verifies project identity and member hashes before exposing reports or data.</p>
      <label className="target-resource">Result manifest or directory <input value={manifestPath} onChange={event => setManifestPath(event.target.value)} placeholder="Path on this service host"/></label>
      <button type="button" disabled={!manifestPath.trim() || imported.isPending} onClick={() => imported.mutate()}>{imported.isPending ? "Verifying results…" : "Verify and import target results"}</button>
    </section>
  </section>;
}
