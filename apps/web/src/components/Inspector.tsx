import { CircleDotDashed, Fingerprint, FlaskConical, Route, Rows3, TriangleAlert } from "lucide-react";
import type { Artifact, WorkflowDefinition } from "../types";
import { displayModuleName, displayArtifactName } from "../branding";
import { ModuleEvidencePanel } from "./ModuleEvidencePanel";
import { LineagePanel } from "./LineagePanel";

type Props = { workflow?: WorkflowDefinition; artifact?: Artifact; row?: Record<string, unknown> };

function display(value: unknown) {
  if (value === null || value === undefined || value === "") return "Not available";
  if (typeof value === "boolean") return value ? "Yes" : "No";
  if (typeof value === "number") return Number.isInteger(value) ? String(value) : value.toPrecision(4);
  return typeof value === "object" ? JSON.stringify(value, null, 2) : String(value);
}

export function Inspector({ workflow, artifact, row }: Props) {
  const module = row?.analysis_module ?? row?.module ?? row?.consensus_module ?? row?.module_id;
  return (
    <aside className="inspector" id="scientific-detail" aria-label="Scientific detail">
      <div className="panel-heading">
        <div><span className="eyebrow">Detail</span><h3>{row ? "Selected record" : "Scientific context"}</h3></div>
        <CircleDotDashed size={19} />
      </div>
      {row && workflow && <details className="method-summary"><summary>Method and column definitions</summary><p>{workflow.method}</p><p>Unit: {workflow.analysis_unit}</p></details>}
      {module != null && artifact && ["programs.cosplice", "programs.psi_network", "programs.module_functions", "programs.diagnostics", "programs.external_validation", "programs.hub_events"].includes(artifact.analysis_id) && <ModuleEvidencePanel projectId={artifact.project_id} artifactId={artifact.artifact_id} module={String(module)} />}
      {row ? (
        <dl className="detail-list">
          {Object.entries(row).map(([key, value]) => (
            <div key={key}><dt>{key.replaceAll("_", " ")}</dt><dd>{display(value)}</dd>{Boolean((artifact?.metadata.column_definitions as Record<string, string> | undefined)?.[key]) && <small>{(artifact?.metadata.column_definitions as Record<string, string>)[key]}</small>}</div>
          ))}
        </dl>
      ) : (
        <>
          {workflow && (
            <div className="method-summary">
              <span>Analysis purpose</span><p>{workflow.question}</p>
              <span>Method</span><p>{workflow.method}</p>
              <span>Analysis unit</span><p>{workflow.analysis_unit}</p>
              <span>Outputs</span><p>{workflow.output_summary}</p>
            </div>
          )}
          {artifact && (
            <div className="artifact-summary">
              <strong>{displayArtifactName(artifact.name)}</strong>
              <dl>
                <div><dt><Rows3 size={15} />Rows</dt><dd>{artifact.row_count ?? "n/a"}</dd></div>
                <div><dt><Route size={15} />Workflow</dt><dd>{artifact.analysis_id}</dd></div>
                <div><dt><Fingerprint size={15} />Run</dt><dd>{artifact.run_fingerprint.slice(0, 12)}</dd></div>
                <div><dt><FlaskConical size={15} />Entity</dt><dd>{artifact.entity_type}</dd></div>
              </dl>
            </div>
          )}
        </>
      )}
      {workflow?.limitations && (
        <div className="claim-boundary">
          <span><TriangleAlert size={13} />Interpretation boundary</span>
          <p>{workflow.limitations}</p>
        </div>
      )}
      {artifact && <details className="method-summary"><summary>Source, scope and effective parameters</summary><pre style={{whiteSpace:"pre-wrap", overflowWrap:"anywhere"}}>{JSON.stringify({comparisons:artifact.contrast_ids, sources:artifact.parent_artifact_ids, parameters:artifact.metadata.effective_evidence_parameters ?? artifact.parameters}, null, 2)}</pre></details>}
      {artifact && <LineagePanel key={artifact.artifact_id} projectId={artifact.project_id} artifactId={artifact.artifact_id}/>}
    </aside>
  );
}
