import { Download, ExternalLink, FileArchive, FileText, Image as ImageIcon } from "lucide-react";
import { api } from "../api";
import type { Artifact } from "../types";
import { ArtifactTable } from "./ArtifactTable";
import { PlotViewer } from "./PlotViewer";

type Props = {
  projectId: string;
  artifact?: Artifact;
  selectedIds: Set<string>;
  onToggle: (rowId: string) => void;
  onFocus: (row: Record<string, unknown>) => void;
};

function extension(path: string) {
  return path.split(".").pop()?.toLowerCase() ?? "";
}

export function ArtifactViewer({ projectId, artifact, selectedIds, onToggle, onFocus }: Props) {
  if (!artifact) {
    return (
      <div className="viewer-empty">
        <FileText size={30} />
        <strong>No cached result for this workflow</strong>
        <span>Opening a page never launches analysis. Use Run when you are ready.</span>
      </div>
    );
  }
  if (artifact.kind === "table" || artifact.kind === "selection") {
    return (
      <ArtifactTable
        projectId={projectId}
        artifact={artifact}
        selectedIds={selectedIds}
        onToggle={onToggle}
        onFocus={onFocus}
      />
    );
  }
  if (artifact.kind === "plot_spec") {
    return <PlotViewer projectId={projectId} artifact={artifact} />;
  }

  const ext = extension(artifact.path);
  const source = api.artifactFileUrl(projectId, artifact.artifact_id);
  const download = api.artifactFileUrl(projectId, artifact.artifact_id, true);
  if (["png", "jpg", "jpeg", "webp", "svg"].includes(ext)) {
    return (
      <div className="file-preview image-preview">
        <div className="preview-toolbar"><ImageIcon size={15} /><span>{artifact.name}</span><a href={download}><Download size={15} />Export</a></div>
        <div className="image-stage"><img src={source} alt={artifact.name} /></div>
      </div>
    );
  }
  if (ext === "pdf") {
    return (
      <div className="file-preview">
        <div className="preview-toolbar"><FileText size={15} /><span>{artifact.name}</span><a href={download}><Download size={15} />Export</a></div>
        <object className="document-frame" data={source} type="application/pdf"><a href={source}>Open PDF</a></object>
      </div>
    );
  }
  if (["html", "htm"].includes(ext)) {
    return (
      <div className="file-preview">
        <div className="preview-toolbar"><FileText size={15} /><span>{artifact.name}</span><a href={source} target="_blank" rel="noreferrer"><ExternalLink size={15} />Open</a><a href={download}><Download size={15} />Export</a></div>
        <iframe className="document-frame" src={source} title={artifact.name} />
      </div>
    );
  }
  return (
    <div className="download-card">
      <FileArchive size={34} />
      <strong>{artifact.name}</strong>
      <span>{artifact.path}</span>
      <a className="primary-action" href={download}><Download size={16} />Download artifact</a>
    </div>
  );
}
