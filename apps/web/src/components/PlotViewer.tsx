import { useEffect, useRef } from "react";
import { useQuery } from "@tanstack/react-query";
import type { Result } from "vega-embed";
import { ChartNoAxesColumnIncreasing, LoaderCircle } from "lucide-react";
import { api } from "../api";
import type { Artifact } from "../types";

type Props = {
  projectId: string;
  artifact?: Artifact;
};

export function PlotViewer({ projectId, artifact }: Props) {
  const host = useRef<HTMLDivElement>(null);
  const view = useRef<Result | null>(null);
  const query = useQuery({
    queryKey: ["plot", projectId, artifact?.artifact_id],
    queryFn: () => api.artifactContent(projectId, artifact!.artifact_id),
    enabled: Boolean(artifact),
  });

  useEffect(() => {
    if (!host.current || !query.data) return;
    let disposed = false;
    void import("vega-embed").then(({ default: embed }) =>
      embed(host.current!, query.data!, { actions: false, renderer: "svg" }).then((result) => {
        if (disposed) result.finalize();
        else view.current = result;
      }),
    );
    const observer = new ResizeObserver(() => void view.current?.view.resize().runAsync());
    observer.observe(host.current);
    return () => {
      disposed = true;
      observer.disconnect();
      view.current?.finalize();
      view.current = null;
    };
  }, [query.data]);

  if (!artifact) {
    return (
      <div className="viewer-empty">
        <ChartNoAxesColumnIncreasing size={30} />
        <strong>No cached plot</strong>
        <span>Run the current workspace explicitly to create one.</span>
      </div>
    );
  }
  if (query.isLoading) {
    return <div className="viewer-empty"><LoaderCircle className="spin" size={26} />Loading plot specification</div>;
  }
  if (query.error) {
    return <div className="viewer-empty error">{query.error.message}</div>;
  }
  return <div ref={host} className="plot-host" aria-label={artifact.name} />;
}
