import { useState } from "react";
import { useMutation, useQueryClient } from "@tanstack/react-query";
import { api } from "../api";
import type { ProjectPayload } from "../types";

export type DesignSample = {sample_id: string; biological_sample_id: string; condition: string; batch: string | null; attributes: Record<string, string>};
export function DesignEditor({projectId, onSaved}: {projectId: string; onSaved?: (project: ProjectPayload) => void}) {
  const client = useQueryClient();
  const [rows, setRows] = useState<DesignSample[] | null>(null);
  const [revision, setRevision] = useState("");
  const [reason, setReason] = useState("");
  const [factor, setFactor] = useState("");
  const [validated, setValidated] = useState(false);
  const load = useMutation({mutationFn: () => api.editableDesign(projectId), onSuccess: data => {
    setRows(data.samples); setRevision(data.revision); setValidated(false);
  }});
  const save = useMutation({mutationFn: (apply: boolean) => api.editDesign(projectId, {revision, samples: rows, reason}, apply), onSuccess: (data, apply) => {
    setValidated(!apply);
    if (apply) { setRevision(data.revision!); if (data.project) onSaved?.(data.project); void client.invalidateQueries(); }
  }});
  const factors = [...new Set((rows ?? []).flatMap(row => Object.keys(row.attributes)))];
  function change(index: number, key: string, value: string, attribute = false) {
    setValidated(false);
    setRows(current => current!.map((r, i) => i !== index ? r : attribute ? {...r, attributes: {...r.attributes, [key]: value}} : {...r, [key]: value}));
  }
  return <details className="design-editor"><summary>edit experimental design</summary>
    <p>Sample IDs, input files and imported contrast directions stay fixed. Saving preserves the previous manifest and invalidates old analysis caches. No analysis runs automatically.</p>
    <button className="secondary-action" onClick={() => load.mutate()} disabled={load.isPending}>Load editable design</button>
    {rows && <>
      <div className="design-table-scroll"><table><thead><tr><th>Sample (fixed)</th><th>Biological unit</th><th>Condition</th><th>Batch (optional)</th>{factors.map(f => <th key={f}>{f}</th>)}</tr></thead>
        <tbody>{rows.map((row, index) => <tr key={row.sample_id}><td>{row.sample_id}</td>
          {(["biological_sample_id", "condition", "batch"] as const).map(key => <td key={key}><input aria-label={`${row.sample_id} ${key}`} value={row[key] ?? ""} onChange={e => change(index, key, e.target.value)}/></td>)}
          {factors.map(key => <td key={key}><input aria-label={`${row.sample_id} ${key}`} value={row.attributes[key] ?? ""} onChange={e => change(index, key, e.target.value, true)}/></td>)}
        </tr>)}</tbody></table></div>
      <label>New factor<input value={factor} onChange={e => setFactor(e.target.value)}/></label>
      <button className="secondary-action" disabled={!factor.trim() || factors.includes(factor.trim())} onClick={() => {setRows(rows.map(r => ({...r, attributes: {...r.attributes, [factor.trim()]: ""}}))); setFactor(""); setValidated(false);}}>Add factor column</button>
      <label>Reason for change<input value={reason} onChange={e => {setReason(e.target.value); setValidated(false);}}/></label>
      <button className="secondary-action" disabled={!reason.trim() || save.isPending} onClick={() => save.mutate(false)}>Validate changes</button>
      <button className="secondary-action" disabled={!validated || save.isPending} onClick={() => save.mutate(true)}>Save design revision</button>
      {save.data && <p role="status">{save.data.applied ? "Saved. Previous results are retained but stale; rebuild only the analyses you need." : `Validated. ${save.data.review.batch_design.reason}`}</p>}
    </>}
    {(load.error || save.error) && <p role="alert">{load.error?.message ?? save.error?.message}</p>}
  </details>;
}
