import { useMutation } from "@tanstack/react-query";
import { api } from "../api";
import { DesignEditor } from "./DesignEditor";
import type { ProjectPayload } from "../types";

export function RunPlanPanel({projectId, workflowId, parameters, onDesignSaved}: {projectId: string; workflowId: string; parameters: Record<string, unknown>; onDesignSaved?: (project: ProjectPayload) => void}) {
  const plan = useMutation({mutationFn: () => api.workflowPlan(projectId, workflowId, parameters)});
  const design = useMutation({mutationFn: () => api.projectDesign(projectId)});
  return <section className="run-plan-panel" aria-label="Prerequisites and experimental design">
    <DesignEditor projectId={projectId} key={projectId} onSaved={project => {plan.reset(); design.reset(); onDesignSaved?.(project);}}/>
    <div className="run-plan-actions">
      <button className="secondary-action" disabled={plan.isPending} onClick={() => plan.mutate()}>{plan.isPending ? "Checking inputs..." : "Check prerequisites"}</button>
      <button className="secondary-action" disabled={design.isPending} onClick={() => design.mutate()}>Review experimental design</button>
      <small>Read-only checks. No analysis or download.</small>
    </div>
    {(plan.error || design.error) && <p role="alert">{plan.error?.message ?? design.error?.message}</p>}
    {plan.data && <div aria-live="polite">
      <strong>{plan.data.ready ? "Declared prerequisites available" : "Prerequisites need attention"}</strong>
      {[...plan.data.errors, ...plan.data.warnings].map((text) => <p key={text}>{text}</p>)}
      <details><summary>Source dependencies ({plan.data.parent_artifacts.length})</summary>
        {plan.data.parent_artifacts.map((p) => <p key={p.artifact_id}>{p.name} · run {p.run_fingerprint.slice(0,10)}</p>)}
        <p>Input files: {(plan.data.input_bytes/1024**3).toFixed(3)} GiB. {plan.data.memory_estimate}</p><p>{plan.data.scope}</p>
      </details>
    </div>}
    {design.data && <details open><summary>Declared comparisons and shared controls</summary>
      <p>{design.data.factor_policy}</p><p>{design.data.batch_design.reason}</p>
      <div className="design-table-scroll"><table><thead><tr><th>Comparison</th><th>Effect definition</th><th>Biological units (numerator / denominator)</th><th>Pairing</th></tr></thead>
        <tbody>{design.data.contrasts.map((c) => <tr key={c.contrast_id}><td>{c.label}<small>{c.contrast_id}</small></td><td>{c.effect_definition}</td><td>{c.numerator_biological_units} / {c.denominator_biological_units}</td><td>{c.pairing}</td></tr>)}</tbody></table></div>
      {design.data.comparison_pairs.map((p) => <p key={`${p.contrast_A}:${p.contrast_B}`}><strong>{p.label_A} compared with {p.label_B}</strong><br/>
        Same control samples: {p.matched_denominator ? "Yes" : "No"}. Shared biological units: {p.shared_biological_units.length}. {p.interpretation}</p>)}
    </details>}
  </section>;
}
