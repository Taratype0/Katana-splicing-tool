import {
  Activity,
  FileChartColumn,
  FolderKanban,
  GitBranch,
  Microscope,
  Network,
  ScanLine,
  ShieldCheck,
  SlidersHorizontal,
  Target,
} from "lucide-react";
import type { Workspace } from "../types";
import { Fragment } from "react";

const workflowStages = [
  { title: "1 / Inputs and QC", ids: ["project", "quality", "landscape"] },
  { title: "2 / Differential analysis", ids: ["differential"] },
  { title: "3 / Networks and functions", ids: ["programs", "prioritize"] },
  { title: "4 / Validation and reports", ids: ["validate", "reports"] },
  { title: "Selected-target follow-up", ids: ["targets"] },
  { title: "Single-cell analysis", ids: ["single-cell"] },
];

const icons = {
  folder: FolderKanban,
  shield: ShieldCheck,
  activity: Activity,
  split: GitBranch,
  target: Target,
  scan: ScanLine,
  network: Network,
  cells: Microscope,
  report: FileChartColumn,
} as const;

type Props = {
  workspaces: Workspace[];
  activeId: string;
  onChange: (workspaceId: string) => void;
};

export function WorkspaceNav({ workspaces, activeId, onChange }: Props) {
  return (
    <aside className="workspace-nav" aria-label="Analysis workspaces">
      <div className="brand-lockup">
        <img className="brand-mark" src="/katana_icon.png" alt="Katana logo" width={48} height={48} />
        <div>
          <strong>Katana</strong>
          <small>RNA splicing analysis</small>
        </div>
      </div>

      <nav>
        {workflowStages.map((stage) => <Fragment key={stage.title}>
          <h2 style={{ fontSize: 11, letterSpacing: ".04em", margin: "20px 10px 8px", opacity: .7 }}>{stage.title}</h2>
          {stage.ids.flatMap((id) => workspaces.filter((w) => w.workspace_id === id)).map((workspace) => {
          const Icon = icons[workspace.icon as keyof typeof icons] ?? SlidersHorizontal;
          return (
            <button
              type="button"
              key={workspace.workspace_id}
              className={activeId === workspace.workspace_id ? "nav-item active" : "nav-item"}
              disabled={!workspace.enabled}
              aria-label={workspace.label}
              onClick={() => onChange(workspace.workspace_id)}
              title={workspace.enabled ? workspace.purpose : `${workspace.purpose} Planned phase ${workspace.phase}`}
            >
              <Icon size={18} strokeWidth={1.8} />
              <span>{workspace.label}</span>
              {!workspace.enabled && <small>Phase {workspace.phase}</small>}
            </button>
          );
        })}</Fragment>)}
      </nav>

      <div className="local-note">
        <ShieldCheck size={17} />
        <div>
          <strong>Local by design</strong>
          <span>Inputs remain on this workstation.</span>
        </div>
      </div>
    </aside>
  );
}
