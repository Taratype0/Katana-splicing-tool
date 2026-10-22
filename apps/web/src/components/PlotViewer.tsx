import { useEffect, useRef, useState } from "react";
import { useQuery } from "@tanstack/react-query";
import type { Result } from "vega-embed";
import { ChartNoAxesColumnIncreasing, LoaderCircle } from "lucide-react";
import { api } from "../api";
import { displayModuleName } from "../branding";
import type { Artifact } from "../types";

type Props = { projectId: string; artifact?: Artifact };

export function PlotViewer({ projectId, artifact }: Props) {
  const host = useRef<HTMLDivElement>(null);
  const embedded = useRef<Result | null>(null);
  const [renderError, setRenderError] = useState("");
  const query = useQuery({
    queryKey: ["plot", projectId, artifact?.artifact_id],
    queryFn: () => api.artifactContent(projectId, artifact!.artifact_id),
    enabled: Boolean(artifact),
  });

  useEffect(() => {
    const element = host.current;
    if (!element || !query.data) return;
    let disposed = false;
    setRenderError("");
    element.replaceChildren();
    void import("vega-embed").then(({ default: embed }) =>
      embed(element, query.data!, { actions: { export: true, source: false, compiled: false, editor: false }, renderer: "svg" })
        .then((result) => {
          if (disposed) result.finalize();
          else { embedded.current = result; window.requestAnimationFrame(() => void result.view.resize().runAsync()); }
        }),
    ).catch((error: unknown) => { if (!disposed) setRenderError(String(error)); });
    const observer = new ResizeObserver(() => {
      window.requestAnimationFrame(() => void embedded.current?.view.resize().runAsync());
    });
    observer.observe(element);
    return () => {
      disposed = true;
      observer.disconnect();
      embedded.current?.finalize();
      embedded.current = null;
      element.replaceChildren();
    };
  }, [query.data]);

  if (!artifact) {
    return <div className="viewer-empty"><ChartNoAxesColumnIncreasing size={30} /><strong>No cached plot</strong></div>;
  }
  if (query.isLoading) {
    return <div className="viewer-empty"><LoaderCircle className="spin" size={26} />Loading plot</div>;
  }
  if (query.error) return <div className="viewer-empty error">{query.error.message}</div>;
  return <div className="plot-stage">
    <div className="preview-toolbar"><span>{displayModuleName(artifact.name)}</span>{["png", "svg"].map((format) => <a key={format} href={api.artifactExportUrl(projectId, artifact.artifact_id, format)}>Download {format.toUpperCase()}</a>)}</div>
    {renderError && <div className="warning-strip">Interactive preview failed: {renderError}. The cached PNG/SVG remains available.</div>}
    <div ref={host} className="plot-host" aria-label={displayModuleName(artifact.name)} />
  </div>;
}
