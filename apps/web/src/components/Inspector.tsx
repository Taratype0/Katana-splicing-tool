import { CircleDotDashed, Fingerprint, FlaskConical, Route, Rows3, TriangleAlert } from "lucide-react";
import type { Artifact, WorkflowDefinition } from "../types";

type Props = { workflow?: WorkflowDefinition; artifact?: Artifact; row?: Record<string, unknown> };

function display(value: unknown) {
  if (value === null || value === undefined || value === "") return "Not available";
  if (typeof value === "boolean") return value ? "Yes" : "No";
  if (typeof value === "number") return Number.isInteger(value) ? String(value) : value.toPrecision(4);
  return String(value);
}

export function Inspector({ workflow, artifact, row }: Props) {
  return (
    <aside className="inspector">
      <div className="panel-heading">
        <div><span className="eyebrow">Detail</span><h3>{row ? "Selected record" : "Scientific context"}</h3></div>
        <CircleDotDashed size={19} />
      </div>
      {row ? (
        <dl className="detail-list">
          {Object.entries(row).map(([key, value]) => (
            <div key={key}><dt>{key.replaceAll("_", " ")}</dt><dd>{display(value)}</dd></div>
          ))}
        </dl>
      ) : (
        <>
          {workflow && (
            <div className="method-summary">
              <span>Question</span><p>{workflow.question}</p>
              <span>Method</span><p>{workflow.method}</p>
              <span>Analysis unit</span><p>{workflow.analysis_unit}</p>
              <span>Outputs</span><p>{workflow.output_summary}</p>
            </div>
          )}
          {artifact && (
            <div className="artifact-summary">
              <strong>{artifact.name}</strong>
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
    </aside>
  );
}
