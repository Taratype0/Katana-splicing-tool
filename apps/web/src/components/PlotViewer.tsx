import { useEffect, useRef } from "react";
import { useQuery } from "@tanstack/react-query";
import type { Result } from "vega-embed";
import { ChartNoAxesColumnIncreasing, LoaderCircle } from "lucide-react";
import { api } from "../api";
import type { Artifact } from "../types";

type Props = { projectId: string; artifact?: Artifact };

export function PlotViewer({ projectId, artifact }: Props) {
  const host = useRef<HTMLDivElement>(null);
  const embedded = useRef<Result | null>(null);
  const query = useQuery({
    queryKey: ["plot", projectId, artifact?.artifact_id],
    queryFn: () => api.artifactContent(projectId, artifact!.artifact_id),
    enabled: Boolean(artifact),
  });

  useEffect(() => {
    const element = host.current;
    if (!element || !query.data) return;
    let disposed = false;
    element.replaceChildren();
    void import("vega-embed").then(({ default: embed }) =>
      embed(element, query.data!, { actions: { export: true, source: false, compiled: false, editor: false }, renderer: "svg" })
        .then((result) => {
          if (disposed) result.finalize();
          else embedded.current = result;
        }),
    );
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
  return <div className="plot-stage"><div ref={host} className="plot-host" aria-label={artifact.name} /></div>;
}
