import { CircleDotDashed, Fingerprint, Route, Rows3 } from "lucide-react";
import type { Artifact } from "../types";

type Props = {
  artifact?: Artifact;
  row?: Record<string, unknown>;
};

function display(value: unknown) {
  if (value === null || value === undefined || value === "") return "Not available";
  if (typeof value === "number") return Number.isInteger(value) ? String(value) : value.toPrecision(4);
  return String(value);
}

export function Inspector({ artifact, row }: Props) {
  return (
    <aside className="inspector">
      <div className="panel-heading">
        <div><span className="eyebrow">Inspector</span><h3>{row ? "Selected event" : "Artifact context"}</h3></div>
        <CircleDotDashed size={19} />
      </div>
      {row ? (
        <dl className="detail-list">
          {Object.entries(row).map(([key, value]) => (
            <div key={key}><dt>{key.replaceAll("_", " ")}</dt><dd>{display(value)}</dd></div>
          ))}
        </dl>
      ) : artifact ? (
        <div className="artifact-summary">
          <strong>{artifact.name}</strong>
          <p>This panel explains provenance and valid downstream actions without changing the result.</p>
          <dl>
            <div><dt><Rows3 size={15} />Rows</dt><dd>{artifact.row_count ?? "n/a"}</dd></div>
            <div><dt><Route size={15} />Workflow</dt><dd>{artifact.analysis_id}</dd></div>
            <div><dt><Fingerprint size={15} />Run</dt><dd>{artifact.run_fingerprint.slice(0, 12)}</dd></div>
          </dl>
        </div>
      ) : (
        <div className="inspector-empty">Select a cached artifact or table row to inspect its evidence trail.</div>
      )}
      <div className="claim-boundary">
        <span>Claim boundary</span>
        <p>Phase 1 fixture output validates the platform only. It is not a biological result.</p>
      </div>
    </aside>
  );
}
