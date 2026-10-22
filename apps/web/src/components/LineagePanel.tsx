import { useMutation } from "@tanstack/react-query";
import { api } from "../api";

export function LineagePanel({projectId, artifactId}: {projectId:string;artifactId:string}) {
  const lineage = useMutation({mutationFn: () => api.lineage(projectId,artifactId)});
  return <details className="method-summary"><summary>Analysis path and filtering provenance</summary>
    <button className="secondary-action" disabled={lineage.isPending} onClick={() => lineage.mutate()}>Read recorded analysis path</button>
    {lineage.error && <p role="alert">{lineage.error.message}</p>}
    {lineage.data && <><p>{lineage.data.scope}</p>{lineage.data.nodes.map(node => <div key={node.artifact_id}>
      <strong>{node.name}</strong><p>{node.workflow}: {node.rows ?? "unavailable"} rows ({node.entity_type}); {node.state}</p>
      {Boolean(node.parents?.length) && <small>Derived from: {node.parents!.map(id => lineage.data!.nodes.find(n => n.artifact_id === id)?.name ?? id).join("; ")}</small>}
    </div>)}</>}
  </details>;
}
