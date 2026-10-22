"""Read-only ancestry, with counts kept in their original analysis units."""
from katana_next.workflows.sources import source_complete


def artifact_lineage(context, artifact_id):
    root = context.catalog.get_artifact(artifact_id)
    if root is None or root.project_id != context.manifest.project_id:
        raise ValueError("Artifact is not in the open project")
    queue, seen, nodes, edges = [artifact_id], set(), [], []
    while queue:
        identity = queue.pop(0)
        if identity in seen:
            continue
        if len(seen) >= 1000:
            raise ValueError("More than 1000 ancestor artifacts; narrow the inspected result")
        seen.add(identity)
        artifact = context.catalog.get_artifact(identity)
        if artifact is None:
            nodes.append({"artifact_id":identity,"name":"Unavailable ancestor","state":"missing","rows":None,"entity_type":"unknown","workflow":"unknown"})
            continue
        if artifact.project_id != context.manifest.project_id:
            raise ValueError("Cross-project ancestry is not implicitly transferable")
        state = "missing_file" if not artifact.path.is_file() else "old_project_revision" if artifact.metadata.get("project_snapshot_fingerprint") != context.snapshot.manifest_fingerprint else "complete" if source_complete(context,artifact) else "incomplete_or_old_method"
        nodes.append({"artifact_id":identity,"name":artifact.name,"state":state,"rows":artifact.row_count,
            "entity_type":artifact.entity_type,"workflow":artifact.analysis_id,"run_fingerprint":artifact.run_fingerprint,
            "comparisons":list(artifact.contrast_ids),"parameters":artifact.parameters,"parents":list(artifact.parent_artifact_ids),
            "scope":artifact.metadata.get("interpretation",artifact.metadata.get("limitations",""))})
        for parent in artifact.parent_artifact_ids:
            edges.append({"source":parent,"target":identity})
            queue.append(parent)
    return {"root":artifact_id,"nodes":nodes,"edges":edges,
        "scope":"Exact recorded ancestry only. Row counts are not subtracted across gene, event, sample, module or comparison units. A missing count is unavailable, not zero. Open the Quality attrition table for cumulative filter losses."}
