import { Database, Fingerprint, GitCompareArrows, HardDrive } from "lucide-react";
import type { ProjectPayload } from "../types";

type Props = { project: ProjectPayload; selectedContrastCount: number; workflowLabel: string };

function shortHash(value: string) {
  return `${value.slice(0, 8)}...${value.slice(-6)}`;
}

export function ContextBar({ project, selectedContrastCount, workflowLabel }: Props) {
  return (
    <section className="context-bar" aria-label="Active analysis context">
      <div className="context-identity">
        <span className="status-dot" />
        <div><small>Active project</small><strong>{project.snapshot.name}</strong></div>
      </div>
      <div className="context-field">
        <GitCompareArrows size={16} />
        <span><small>Analysis context</small><strong>{workflowLabel}</strong></span>
      </div>
      <div className="context-field">
        <Database size={16} />
        <span><small>Comparisons</small><strong>{selectedContrastCount || project.snapshot.contrast_count} selected</strong></span>
      </div>
      <div className="context-field wide">
        <HardDrive size={16} />
        <span><small>External output root</small><strong title={project.snapshot.output_root}>{project.snapshot.output_root}</strong></span>
      </div>
      <div className="context-field fingerprint">
        <Fingerprint size={16} />
        <span><small>{project.snapshot.genome_build}</small><strong>{shortHash(project.snapshot.manifest_fingerprint)}</strong></span>
      </div>
    </section>
  );
}
