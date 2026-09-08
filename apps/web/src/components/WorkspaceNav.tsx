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
        <div className="brand-mark" aria-hidden="true">
          <span />
          <span />
          <span />
        </div>
        <div>
          <strong>Katana</strong>
          <small>Evidence workspace</small>
        </div>
      </div>

      <nav>
        {workspaces.map((workspace) => {
          const Icon = icons[workspace.icon as keyof typeof icons] ?? SlidersHorizontal;
          return (
            <button
              type="button"
              key={workspace.workspace_id}
              className={activeId === workspace.workspace_id ? "nav-item active" : "nav-item"}
              disabled={!workspace.enabled}
              onClick={() => onChange(workspace.workspace_id)}
              title={workspace.enabled ? workspace.purpose : `${workspace.purpose} Planned phase ${workspace.phase}`}
            >
              <Icon size={18} strokeWidth={1.8} />
              <span>{workspace.label}</span>
              {!workspace.enabled && <small>Phase {workspace.phase}</small>}
            </button>
          );
        })}
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
