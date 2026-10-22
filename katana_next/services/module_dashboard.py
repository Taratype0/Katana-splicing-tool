"""Build linked module views only from immutable, completed discovery artifacts."""
from pathlib import Path
import pandas as pd
import pyarrow.parquet as pq

from katana_next.workflows.sources import validate_source, source_complete
from katana_next.science.module_dashboard import build_module_dashboard


def module_dashboard(context, artifact_id, module, *, include_measurements=False, **options):
    source = context.catalog.get_artifact(artifact_id)
    validate_source(context, source)
    if source.analysis_id not in {"programs.cosplice", "programs.psi_network", "programs.module_functions", "programs.diagnostics", "programs.external_validation", "programs.hub_events"}:
        raise ValueError("Select a completed cosplice module artifact")
    run = source.metadata.get("source_discovery_run", source.run_fingerprint)
    artifacts = [a for a in context.catalog.list_artifacts(context.manifest.project_id)
                 if a.kind == "table" and (a.analysis_id in {"programs.cosplice", "programs.psi_network"} and a.run_fingerprint == run
                    or a.analysis_id in {"programs.module_functions", "programs.hub_events"} and a.metadata.get("source_discovery_run") == run)
                 and (a.analysis_id != "programs.hub_events" or source.analysis_id != "programs.hub_events" or a.run_fingerprint == source.run_fingerprint)
                 and a.metadata.get("project_snapshot_fingerprint") == context.snapshot.manifest_fingerprint and source_complete(context, a)]
    hub_anchor = next((a for a in artifacts if a.analysis_id == 'programs.hub_events' and a.path.stem == 'hub_events__metrics'), None)
    if hub_anchor:
        artifacts = [a for a in artifacts if a.analysis_id != 'programs.hub_events' or a.run_fingerprint == hub_anchor.run_fingerprint]
    used = []

    def read(*stems, required=False):
        for stem in stems:
            match = next((a for a in artifacts if Path(a.path).stem == stem), None)
            if match:
                validate_source(context, match)
                used.append({"artifact_id": match.artifact_id, "name": match.name, "run_fingerprint": match.run_fingerprint})
                return pq.read_table(match.path).to_pandas().drop(columns=["row_id"], errors="ignore")
        if required:
            raise ValueError(f"Completed discovery run has no cached {stems[0]} table")
        return None

    assignments = read("consensus_modules__assignments", "network_observed__assignments", required=True)
    matrix = read("matrix__psi", "stable_qc__matrix", required=True)
    key = next((c for c in ("event_key", "Unnamed: 0") if c in matrix), None)
    if key is None:
        raise ValueError("Cached matrix has no event identifier column")
    matrix = matrix.set_index(key)
    samples = pd.DataFrame([{"sample_id": s.sample_id, "biological_sample_id": s.biological_sample_id, "condition": s.condition} for s in context.manifest.samples])
    cached_samples = read('matrix__samples')
    if cached_samples is not None:
        if not {'sample_id','condition','biological_sample_id'} <= set(cached_samples):
            raise ValueError('Cached discovery samples must declare biological identity and condition')
        samples = cached_samples
    if not set(matrix.columns) <= set(samples.sample_id):
        raise ValueError('Discovery sample metadata does not cover the PSI matrix')
    samples = samples.loc[samples.sample_id.isin(matrix.columns)].copy()
    metadata = read("matrix__events")
    edges = read("network_observed__edges")
    scores = read("consensus_modules__scores", "network_observed__scores")
    hub_edges = read('hub_events__display_edges')
    hubs = read('hub_events__metrics')
    if hubs is not None:
        selected_hubs = hubs.loc[hubs.module.astype(str).eq(str(module)) & hubs.condition.astype(str).eq(str(options.get('condition','all')))]
        if len(selected_hubs):
            options['node_metrics'] = selected_hubs
    if hub_edges is not None:
        condition = str(options.get('condition','all'))
        edges = hub_edges.loc[hub_edges.module.astype(str).eq(str(module)) & hub_edges.condition.astype(str).eq(condition)].copy()
        options['edges_are_scope_specific'] = True
    result = build_module_dashboard(matrix, assignments, samples, module, metadata=metadata, edges=edges, scores=scores, **options)
    result["discovery_run"] = run
    result["sources"] = used
    if hubs is not None:
        selected = hubs.loc[hubs.module.astype(str).eq(str(module)) & hubs.condition.astype(str).eq(result['condition'])].copy()
        result['hub_method'] = 'Full observed-PSI module adjacency and PC1/kME; biological-block bootstrap; candidates must pass both all-edge and no-same-gene criteria. All ranks and rejection reasons remain exportable.'
        result['hub_status'] = 'scope_not_computed' if selected.empty else 'no_hubs_passed'
        if not selected.empty:
            selected = selected.assign(connectivity=selected.kWithin).sort_values(['rank_all_edges','event_key'],kind='stable')
            result['total_edges'] = int(selected.supported_degree.sum()//2)
            result['connectivity_definition'] = 'Complete module adjacency sum (kWithin) from the explicit hub analysis; no edge threshold/cap in connectivity. Node sizes use this value. The drawing alone is thinned.'
            # Do not turn an unsupported/underpowered network into a hub list.
            eligible = selected.loc[selected.hub_candidate.eq(True) & selected.supported_degree.gt(0) & selected.cross_gene_supported_degree.gt(0)]
            if result['status'] != 'available' or selected.supported_degree.sum() == 0:
                eligible = eligible.iloc[:0]
                result['hub_status'] = 'no_supported_edges'
            elif len(eligible):
                result['hub_status'] = 'available'
            eligible = eligible.assign(connectivity=eligible.kWithin)
            import json
            result['hub_events'] = json.loads(eligible.to_json(orient='records'))
            result['hub_comparison'] = json.loads(selected.to_json(orient='records'))
            result['hub_analysis_run'] = next((a.run_fingerprint for a in artifacts if a.path.stem == 'hub_events__metrics'), '')
    functions = read("interpret_modules__gene_sets", "module_functions__enrichment")
    result["enrichment"] = []
    result["enrichment_total"] = 0
    if functions is not None and "module" in functions:
        import json
        functions = functions.loc[functions.module.astype(str).eq(module)]
        result["enrichment_total"] = len(functions)
        if "overlap_count" in functions:
            functions = functions.loc[functions.overlap_count.gt(0)]
        if "FDR" in functions:
            order = ["FDR"] + (["pvalue"] if "pvalue" in functions else [])
            functions = functions.sort_values(order, kind="stable")
        result["enrichment"] = json.loads(functions.head(20).to_json(orient="records"))
    result["enrichment_scope"] = "Genes from this frozen discovery module; not a condition-specific DEG list. All terms, raw P, FDR and contributing genes remain in the source table."
    if include_measurements:
        import json
        keys = [e["event_key"] for e in result["event_options"]]
        selected = samples if result["condition"] == "all" else samples.loc[samples.condition.eq(result["condition"])]
        result["event_measurements"] = {"events": keys, "samples": json.loads(selected.to_json(orient="records")),
            "psi": json.loads(matrix.loc[keys, selected.sample_id].to_json(orient="values"))}
    return result
