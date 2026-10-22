"""Event-level module views; display thinning never changes module membership."""
from __future__ import annotations

import json
import numpy as np
import pandas as pd
import networkx as nx


def _records(frame):
    return json.loads(frame.to_json(orient="records"))


def pairwise_pearson(values: np.ndarray, minimum: int = 15):
    """Pairwise-complete correlations, not zero-imputed PSI correlations."""
    valid = np.isfinite(values).astype(float)
    x = np.where(np.isfinite(values), values, 0.0)
    n = valid @ valid.T
    sums = x @ valid.T
    squares = (x * x) @ valid.T
    safe_n = np.maximum(n, 1)
    cov = x @ x.T - sums * sums.T / safe_n
    var = np.maximum(squares - sums * sums / safe_n, 0)
    denom = np.sqrt(var * var.T)
    corr = np.divide(cov, denom, out=np.full_like(cov, np.nan), where=denom > 1e-12)
    corr[n < minimum] = np.nan
    np.fill_diagonal(corr, np.nan)
    return np.clip(corr, -1, 1), n.astype(int)


def build_module_dashboard(matrix, assignments, samples, module, *, metadata=None,
                           edges=None, scores=None, condition="all", selected_event=None,
                           hide_same_gene=False, node_limit=80, minimum_samples=15, edges_are_scope_specific=False, node_metrics=None):
    if not module or not 10 <= node_limit <= 200:
        raise ValueError("A module and display node limit between 10 and 200 are required")
    if not matrix.index.is_unique or not matrix.columns.is_unique:
        raise ValueError("Event and sample IDs must be unique; gene names cannot replace event IDs")
    if not assignments.event_key.is_unique or not samples.sample_id.is_unique:
        raise ValueError("Assignments and sample metadata must have unique IDs")
    label = next((c for c in ("analysis_module", "consensus_module", "module", "observed_module") if c in assignments), None)
    if label is None:
        raise ValueError("Module membership is absent")
    keys = assignments.loc[assignments[label].astype(str).eq(module), "event_key"].astype(str).tolist()
    if not keys:
        raise ValueError("This module is absent from the selected discovery run")
    missing = set(keys) - set(matrix.index.astype(str))
    if missing:
        raise ValueError(f"The cached PSI matrix is missing {len(missing)} module events")
    missing_samples = set(matrix.columns) - set(samples.sample_id)
    if missing_samples:
        raise ValueError("Sample metadata does not cover the cached PSI matrix")
    design = samples.set_index("sample_id").loc[matrix.columns].copy()
    if "condition" not in design or design.condition.isna().any():
        raise ValueError("Every matrix sample requires an explicit condition")
    identity_assumed = "biological_sample_id" not in design
    if identity_assumed:
        design["biological_sample_id"] = design.index
    if design.biological_sample_id.isna().any() or design.biological_sample_id.astype(str).eq("").any():
        raise ValueError("Biological sample identity must not be missing")
    conditions = sorted(design.condition.astype(str).unique())
    if condition != "all" and condition not in conditions:
        raise ValueError("Unknown condition; choose one declared in the project")
    data = matrix.loc[keys].astype(float)
    if np.any((data.to_numpy()[np.isfinite(data.to_numpy())] < 0) | (data.to_numpy()[np.isfinite(data.to_numpy())] > 1)):
        raise ValueError("Module views require measured PSI in [0,1], not a residual matrix")
    if condition != "all":
        design = design.loc[design.condition.astype(str).eq(condition)]
        data = data.loc[:, design.index]
    genes = {}
    host_genes = {}
    if metadata is not None and "event_key" in metadata:
        if not metadata.event_key.is_unique:
            raise ValueError("Event annotation must have unique event IDs")
        gene_column = next((c for c in ("gene_symbol", "gene_name", "gene_id", "GeneID") if c in metadata), None)
        if gene_column:
            genes = metadata.set_index("event_key")[gene_column].fillna("").astype(str).to_dict()
        host_column = next((c for c in ('gene_id','GeneID','gene_symbol','gene_name') if c in metadata), None)
        if host_column:
            host_genes = metadata.set_index('event_key')[host_column].fillna('').astype(str).to_dict()
    node_frame = pd.DataFrame({"event_key": keys, "gene": [genes.get(k, "") for k in keys]})
    status, method = "available", "Cached discovery edges; no new module fit"
    independent = int(design.biological_sample_id.nunique())
    cached = (condition == "all" or edges_are_scope_specific) and edges is not None and {"event_a", "event_b"}.issubset(edges) and independent >= minimum_samples
    if cached:
        edge_frame = edges.loc[edges.event_a.isin(keys) & edges.event_b.isin(keys)].copy()
        edge_frame["weight"] = pd.to_numeric(edge_frame.get("adjacency", edge_frame.get("correlation", 0)), errors="coerce")
        edge_frame["pair_n"] = None
        if edges_are_scope_specific:
            edge_frame['pair_n'] = edges.loc[edge_frame.index,'pair_n'] if 'pair_n' in edges else None
            method = 'Cached scope-specific official WGCNA hub-analysis edges; top eight per event for drawing only. Full connectivity uses every estimable module pair.'
    elif independent < minimum_samples:
        status = "insufficient_biological_samples"
        method = f"Within-condition association requires at least {minimum_samples} biological samples"
        edge_frame = pd.DataFrame(columns=["event_a", "event_b", "weight", "correlation", "pair_n"])
    elif len(keys) > 2000:
        status = "module_exceeds_interactive_budget"
        method = "More than 2,000 events: retain all measurements; use exported matrix for an explicit large-network job"
        edge_frame = pd.DataFrame(columns=["event_a", "event_b", "weight", "correlation", "pair_n"])
    elif condition == "all" and design.biological_sample_id.duplicated().any():
        status = "select_condition_for_repeated_units"
        method = "Choose a condition; repeated measures are not independent samples"
        edge_frame = pd.DataFrame(columns=["event_a", "event_b", "weight", "correlation", "pair_n"])
    else:
        # Technical repeats are averaged within a biological unit, never counted twice.
        collapsed = data.T.groupby(design.biological_sample_id, sort=True).mean().T
        corr, pair_n = pairwise_pearson(collapsed.to_numpy(), minimum_samples)
        adjacency = ((1 + corr) / 2) ** 6
        row, col = np.triu_indices(len(keys), 1)
        use = np.isfinite(adjacency[row, col]) & (adjacency[row, col] >= .10)
        edge_frame = pd.DataFrame({"event_a": np.asarray(keys)[row[use]], "event_b": np.asarray(keys)[col[use]],
            "weight": adjacency[row[use], col[use]], "correlation": corr[row[use], col[use]], "pair_n": pair_n[row[use], col[use]]})
        method = "Descriptive pairwise-complete Pearson association; signed adjacency=((1+r)/2)^6, >=0.10; fixed discovery membership"
    edge_frame["same_gene"] = [bool(host_genes.get(a)) and host_genes.get(a) == host_genes.get(b) for a, b in zip(edge_frame.event_a, edge_frame.event_b)]
    coupling_count = int(edge_frame.same_gene.sum())
    degree = pd.concat([edge_frame[["event_a", "weight"]].rename(columns={"event_a": "event_key"}),
                       edge_frame[["event_b", "weight"]].rename(columns={"event_b": "event_key"})]).groupby("event_key").weight.sum()
    node_frame["connectivity"] = node_frame.event_key.map(degree).fillna(0)
    if node_metrics is not None:
        if not node_metrics.event_key.is_unique:
            raise ValueError('Full hub metrics require unique event identities within a scope')
        full = node_metrics.set_index('event_key').kWithin
        node_frame['connectivity'] = node_frame.event_key.map(full).fillna(0)
    node_frame = node_frame.sort_values(["connectivity", "event_key"], ascending=[False, True], kind="stable")
    shown = node_frame.head(node_limit).copy()
    visible = edge_frame.loc[edge_frame.event_a.isin(shown.event_key) & edge_frame.event_b.isin(shown.event_key)].copy()
    if hide_same_gene:
        visible = visible.loc[~visible.same_gene]
    graph = nx.Graph()
    graph.add_nodes_from(shown.event_key)
    for e in visible.itertuples():
        graph.add_edge(e.event_a, e.event_b, weight=float(e.weight))
    # Thin edges for legibility, not for a new biological or significance filter.
    display_edges = set()
    if graph.number_of_edges():
        display_edges.update(tuple(sorted(e)) for e in nx.maximum_spanning_tree(graph).edges)
        for key in graph:
            display_edges.update(tuple(sorted((key, other))) for other in sorted(graph[key], key=lambda other: (-graph[key][other]["weight"], other))[:2])
    layout = nx.spring_layout(graph.edge_subgraph(display_edges).copy() if display_edges and len(set(sum((list(e) for e in display_edges), []))) == len(graph) else graph, seed=17, iterations=40)
    shown["x"] = [float(layout[k][0]) for k in shown.event_key]
    shown["y"] = [float(layout[k][1]) for k in shown.event_key]
    visible = visible.loc[[tuple(sorted((a, b))) in display_edges for a, b in zip(visible.event_a, visible.event_b)]]
    order = design.sort_values(["condition", "biological_sample_id"], kind="stable").index.tolist()
    display_samples = order if len(order) <= 120 else [order[i] for i in np.unique(np.linspace(0, len(order)-1, 120).astype(int))]
    heat = data.loc[shown.event_key, display_samples]
    center = data.loc[shown.event_key].mean(axis=1)
    scale = data.loc[shown.event_key].std(axis=1, ddof=0).replace(0, np.nan)
    z = heat.sub(center, axis=0).div(scale, axis=0)
    z.loc[scale.isna()] = heat.loc[scale.isna()].where(heat.loc[scale.isna()].isna(), 0)
    if selected_event is None:
        selected_event = str(shown.event_key.iloc[0])
    if selected_event not in keys:
        raise ValueError("Selected event is outside this module")
    profile = design.reset_index().rename(columns={"index": "sample_id"})
    profile["psi"] = data.loc[selected_event].reindex(design.index).to_numpy()
    score_records = []
    if scores is not None and "sample_id" in scores and module in scores:
        score_data = scores.set_index("sample_id")[module]
        score_records = [{"sample_id": k, "condition": str(design.loc[k, "condition"]), "value": None if pd.isna(score_data.get(k)) else float(score_data[k])} for k in design.index]
    return {"module": module, "condition": condition, "conditions": conditions, "node_unit": "splice_event", "status": status,
        "connectivity_definition": "Sum of exported edge adjacency weights; discovery exporter may cap edges. Not full WGCNA kWithin or regulatory influence." if cached else "Sum of signed Pearson adjacency weights >=0.10 within the fixed module, using pairwise-complete biological units.",
        "method": method, "biological_samples": independent, "minimum_biological_samples": minimum_samples,
        "identity_assumed": identity_assumed, "total_events": len(keys), "total_edges": len(edge_frame), "same_gene_edges": coupling_count,
        "nodes": _records(shown), "edges": _records(visible[["event_a", "event_b", "weight", "correlation", "pair_n", "same_gene"]]),
        "hub_events": [], "hub_status": "not_computed", "hub_comparison": [],
        "hub_method": "Run cosplice: validate hub events on this completed discovery run. Drawing connectivity is not a validated hub ranking.",
        "event_options": _records(node_frame[["event_key", "gene"]]),
        "heatmap": {"events": shown.event_key.tolist(), "samples": display_samples, "conditions": [str(design.loc[k, "condition"]) for k in display_samples],
                    "z": json.loads(z.to_json(orient="values")), "scale": "Row z-score from PSI; blue=lower, red=higher; missing=grey"},
        "selected_event": selected_event, "event_profile": _records(profile[["sample_id", "biological_sample_id", "condition", "psi"]]),
        "scores": score_records, "score_definition": "Cached module eigengene/program score; its sign is a convention, not inclusion direction",
        "display_rule": "Top connectivity events; maximum spanning forest plus two strongest edges per visible node. All module events remain in source tables. Layout is not biological distance.",
        "interpretation": "Edges are splice-usage associations, not PPI or proven regulation. Same-gene edges can reflect shared reads or event geometry. Condition views do not refit or relabel modules and have no edge P/FDR."}
