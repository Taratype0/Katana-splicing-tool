import { useState } from "react";
import { Image as ImageIcon } from "lucide-react";
import { api } from "../api";
import { displayModuleName } from "../branding";
import type { Artifact } from "../types";

export function ImagePreview({ projectId, artifact }: { projectId: string; artifact: Artifact }) {
  const [mode, setMode] = useState("fit");
  const exports = artifact.metadata?.exports as Record<string, string> | undefined;
  return <div className="file-preview image-preview">
    <div className="preview-toolbar image-toolbar">
      <ImageIcon size={15} /><span>{displayModuleName(artifact.name)}</span>
      {[["fit", "Fit page"], ["width", "Fit width"], ["actual", "100%"]].map(([key, label]) =>
        <button type="button" key={key} aria-pressed={mode === key} onClick={() => setMode(key)}>{label}</button>)}
      {exports?.png && <a href={api.artifactExportUrl(projectId, artifact.artifact_id, "png")}>PNG</a>}
      {exports?.svg && <a href={api.artifactExportUrl(projectId, artifact.artifact_id, "svg")}>SVG</a>}
      {!exports?.png && <a href={api.artifactFileUrl(projectId, artifact.artifact_id, true)}>Export</a>}
    </div>
    <div className={`image-stage image-stage-${mode}`} aria-label="Image viewport">
      <img src={api.artifactFileUrl(projectId, artifact.artifact_id)} alt={displayModuleName(artifact.name)} />
    </div>
  </div>;
}
