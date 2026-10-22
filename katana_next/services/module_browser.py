"""Read only module evidence from one frozen discovery lineage."""
import pyarrow.dataset as ds
import pyarrow.parquet as pq

from katana_next.workflows.sources import validate_source, source_complete


def module_evidence(context, artifact_id, module, *, limit=30):
    if not module or len(module) > 200 or not 1 <= limit <= 100:
        raise ValueError("A module label and a page limit of 1-100 are required")
    source = context.catalog.get_artifact(artifact_id)
    validate_source(context, source)
    if source.analysis_id not in {"programs.djexpress", "programs.cosplice", "programs.psi_network", "programs.module_functions", "programs.diagnostics", "programs.external_validation", "programs.hub_events"}:
        raise ValueError("Module browser requires a CoSplice or fixed-module diagnostic artifact")
    run = source.metadata.get("source_discovery_run", source.run_fingerprint)
    all_sources = context.catalog.list_artifacts(context.manifest.project_id)
    related = [a for a in all_sources if a.kind == "table" and
        a.metadata.get("project_snapshot_fingerprint") == context.snapshot.manifest_fingerprint and source_complete(context, a) and (
        a.analysis_id in {"programs.djexpress", "programs.cosplice", "programs.psi_network"} and a.run_fingerprint == run or
        a.analysis_id in {"programs.diagnostics", "programs.external_validation", "programs.module_functions", "programs.hub_events"} and a.metadata.get("source_discovery_run") == run)]
    tables = []
    for artifact in related:
        validate_source(context, artifact)
        columns = pq.read_schema(artifact.path).names
        labels = [name for name in ("module", "module_id", "analysis_module", "consensus_module", "observed_module") if name in columns]
        if not labels:
            continue
        # Prefer analysis membership over original/backbone labels when both exist.
        name = next((c for c in ("analysis_module", "module", "consensus_module", "module_id", "observed_module") if c in labels))
        dataset = ds.dataset(artifact.path, format="parquet")
        expression = ds.field(name).cast("string") == module
        count = dataset.count_rows(filter=expression)
        if not count:
            continue
        scanner = dataset.scanner(filter=expression, batch_size=limit)
        head = scanner.head(limit).to_pandas()
        import json
        rows = json.loads(head.to_json(orient="records"))
        tables.append({"artifact_id": artifact.artifact_id, "name": artifact.name, "entity_type": artifact.entity_type,
            "matching_column": name, "row_count": count, "preview_limit": limit, "rows": rows,
            "workflow": artifact.analysis_id, "run_fingerprint": artifact.run_fingerprint})
    return {"module": module, "discovery_run": run, "tables": tables,
        "scope": "Cached evidence from one discovery lineage. Empty sections mean unavailable evidence, not biological absence. Associations, motifs and enrichment do not prove causality."}
