import { Database, Fingerprint, GitCompareArrows, HardDrive } from "lucide-react";
import type { ProjectPayload } from "../types";

type Props = {
  project: ProjectPayload;
  contrastId: string;
  onContrastChange: (value: string) => void;
};

function shortHash(value: string) {
  return `${value.slice(0, 8)}...${value.slice(-6)}`;
}

export function ContextBar({ project, contrastId, onContrastChange }: Props) {
  const contrasts = project.manifest.contrasts;
  return (
    <section className="context-bar" aria-label="Active analysis context">
      <div className="context-identity">
        <span className="status-dot" />
        <div>
          <small>Active project</small>
          <strong>{project.snapshot.name}</strong>
        </div>
      </div>
      <label className="context-field">
        <GitCompareArrows size={16} />
        <span>
          <small>Contrast</small>
          <select value={contrastId} onChange={(event) => onContrastChange(event.target.value)}>
            <option value="">Project-wide</option>
            {contrasts.map((contrast) => (
              <option key={contrast.contrast_id} value={contrast.contrast_id}>{contrast.label}</option>
            ))}
          </select>
        </span>
      </label>
      <div className="context-field">
        <Database size={16} />
        <span><small>Reference</small><strong>{project.snapshot.genome_build}</strong></span>
      </div>
      <div className="context-field wide">
        <HardDrive size={16} />
        <span><small>Output root</small><strong title={project.snapshot.output_root}>{project.snapshot.output_root}</strong></span>
      </div>
      <div className="context-field fingerprint">
        <Fingerprint size={16} />
        <span><small>Input fingerprint</small><strong>{shortHash(project.snapshot.manifest_fingerprint)}</strong></span>
      </div>
    </section>
  );
}
