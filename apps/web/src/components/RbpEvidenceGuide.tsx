import { RBP_EVIDENCE, rbpEvidenceLayer } from "../rbpEvidence";
import type { Artifact, WorkflowDefinition } from "../types";

type Props = {
  workflows: WorkflowDefinition[];
  artifacts: Artifact[];
  activeWorkflowId: string;
  parameters: Record<string, unknown>;
  onSelect: (workflowId: string, parameters?: Record<string, unknown>) => void;
};

export function RbpEvidenceGuide({workflows, artifacts, activeWorkflowId, parameters, onSelect}: Props) {
  const active = rbpEvidenceLayer(activeWorkflowId, parameters);
  return <section className="rbp-evidence-guide" aria-label="RBP evidence layers">
    <h3>{RBP_EVIDENCE.title_en}</h3>
    <p>The original RBP features are retained. Annotation, candidate regulation and RNA consequences answer different questions.</p>
    <div className="rbp-evidence-grid">
      {RBP_EVIDENCE.layers.map(layer => {
        const count = artifacts.filter(a => rbpEvidenceLayer(a.analysis_id, a.parameters)?.id === layer.id).length;
        const available = workflows.some(w => w.workflow_id === layer.workflow_id && w.enabled);
        return <article key={layer.id} className={active?.id === layer.id ? "active" : ""}>
          <h4>{layer.number}. {layer.title_en}</h4>
          <p lang="zh">{layer.title_zh}</p>
          <p>{layer.question_en}</p><p>{layer.value_en}</p>
          <p className="rbp-evidence-boundary"><strong>Boundary: </strong>{layer.boundary_en}</p>
          <small>{count ? `${count} saved outputs; not gene/event counts.` : "No saved output loaded for this layer; not a negative finding."}</small>
          <button type="button" className="text-action" disabled={!available}
            onClick={() => onSelect(layer.workflow_id, layer.id === "regulation" ? {sequence_scope:"boundary_association"} : undefined)}>
            View {layer.number}. {layer.title_en}
          </button>
        </article>;
      })}
    </div>
    <p>{RBP_EVIDENCE.shared_boundary_en}</p>
  </section>;
}

export function RbpLayerNote({workflowId, parameters}: {workflowId:string; parameters?:Record<string,unknown>}) {
  const layer = rbpEvidenceLayer(workflowId, parameters);
  return layer ? <aside className="rbp-current-layer" aria-label="Current RBP evidence layer">
    <strong>{layer.number}. {layer.title_en}</strong>
    <p>{layer.unit_en}</p><p>{layer.boundary_en}</p>
  </aside> : null;
}
