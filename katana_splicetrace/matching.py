"""Layer-consistent two-stage global bipartite caliper matching for MotifTrace v2.

Extends EventGraph v1 matching with an optional 5th caliper on proximal intronic
flank GC (`flank_GC_caliper`, default 0.08) so genomic isochore GC differences
do not manufacture false-positive G/C-rich motifs in regulatory intronic flanks.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching, min_weight_full_bipartite_matching

BASE_FEATURES = ("GC", "length", "baseline_PSI", "depth")
BASE_CALIPERS = np.array([0.05, 0.20, 0.15, 1.00], dtype=float)


def match_events_v2(
    foreground: pd.DataFrame,
    background: pd.DataFrame,
    *,
    flank_gc_caliper: float | None = 0.08,
) -> tuple[pd.DataFrame, dict[str, float | int], pd.DataFrame]:
    """Globally assign foreground events to unaffected control genes under exact calipers."""
    use_flank = flank_gc_caliper is not None and "flank_GC" in foreground.columns and "flank_GC" in background.columns
    features = (*BASE_FEATURES, "flank_GC") if use_flank else BASE_FEATURES
    calipers = np.append(BASE_CALIPERS, float(flank_gc_caliper)) if use_flank else BASE_CALIPERS.copy()

    required = {"event_id", "gene_id", "event_type", "ambiguous_fraction", *BASE_FEATURES}
    for frame in (foreground, background):
        if not required <= set(frame.columns):
            raise ValueError(f"Missing required columns: {sorted(required - set(frame.columns))}")
        if frame["event_id"].duplicated().any():
            raise ValueError("Duplicate event_id detected")
        if frame["gene_id"].isna().any() or frame["gene_id"].astype(str).str.strip().eq("").any():
            raise ValueError("Non-empty gene_id required to define gene-block independence")
    if set(foreground["gene_id"].astype(str)) & set(background["gene_id"].astype(str)):
        raise ValueError("Foreground and background gene sets overlap")

    fg = foreground.sort_values("event_id").reset_index(drop=True).copy()
    bg = background.sort_values("event_id").reset_index(drop=True).copy()
    genes = np.array(sorted(bg["gene_id"].astype(str).unique()))
    gene_index = {g: i for i, g in enumerate(genes)}
    bg_codes = bg["gene_id"].astype(str).map(gene_index).to_numpy()

    bg_vals = bg[list(features)].to_numpy(dtype=float)
    bg_kinds = bg["event_type"].to_numpy()
    bg_valid = (
        np.isfinite(bg_vals[:, :4]).all(axis=1)
        & (bg_vals[:, 1] > 0)
        & (bg_vals[:, 3] > 0)
        & (bg["ambiguous_fraction"].to_numpy(dtype=float) <= 0.05)
    )
    if use_flank:
        bg_valid &= np.isfinite(bg_vals[:, 4])

    edges: list[tuple] = []
    cand_counts: list[int] = []
    for i, row in enumerate(fg.itertuples(index=False)):
        x = np.array([getattr(row, c) for c in features], dtype=float)
        ok_fg = np.isfinite(x[:4]).all() and x[1] > 0 and x[3] > 0 and float(row.ambiguous_fraction) <= 0.05
        if use_flank:
            ok_fg = ok_fg and np.isfinite(x[4])
        if not ok_fg:
            cand_counts.append(0)
            continue
        with np.errstate(divide="ignore", invalid="ignore"):
            cols = [
                np.abs(bg_vals[:, 0] - x[0]),
                np.abs(bg_vals[:, 1] - x[1]) / x[1],
                np.abs(bg_vals[:, 2] - x[2]),
                np.abs(np.log2(bg_vals[:, 3] / x[3])),
            ]
            if use_flank:
                cols.append(np.abs(bg_vals[:, 4] - x[4]))
            distances = np.column_stack(cols)
        idx = np.flatnonzero(
            bg_valid
            & (bg_kinds == row.event_type)
            & np.isfinite(distances).all(axis=1)
            & (distances <= calipers).all(axis=1)
        )
        cand_counts.append(len(idx))
        if not len(idx):
            continue
        cost = 1.0 + np.square(distances[idx] / calipers).sum(axis=1)
        order = np.lexsort((idx, cost, bg_codes[idx]))
        prev_gene = -1
        for k in order:
            j = int(idx[k])
            g_code = int(bg_codes[j])
            if g_code == prev_gene:
                continue
            prev_gene = g_code
            edges.append((i, g_code, j, float(cost[k]), *distances[j].tolist()))

    dist_cols = ["GC_distance", "length_distance", "PSI_distance", "depth_distance"] + (
        ["flank_GC_distance"] if use_flank else []
    )
    edge_df = pd.DataFrame(edges, columns=["foreground_index", "background_gene_index", "background_index", "cost", *dist_cols])
    if len(edge_df):
        edge_df[["foreground_index", "background_gene_index", "background_index"]] = edge_df[
            ["foreground_index", "background_gene_index", "background_index"]
        ].astype(int)
    fg["candidate_background_events"] = cand_counts

    n, g_count = len(fg), len(genes)
    if n == 0:
        raise ValueError("Empty foreground event table")
    ri = edge_df["foreground_index"].to_numpy(dtype=int)
    ci = edge_df["background_gene_index"].to_numpy(dtype=int)
    graph = csr_matrix((np.ones(len(edge_df)), (ri, ci)), shape=(n, g_count))
    optimum = int((maximum_bipartite_matching(graph, perm_type="column") >= 0).sum())

    max_edge_cost = float(len(calipers) + 1)
    penalty = max_edge_cost * n + 1.0
    dummy_rows = np.arange(n, dtype=int)
    costs = csr_matrix(
        (np.r_[edge_df["cost"].to_numpy(dtype=float), np.full(n, penalty)], (np.r_[ri, dummy_rows], np.r_[ci, g_count + dummy_rows])),
        shape=(n, g_count + n),
    )
    rows, cols = min_weight_full_bipartite_matching(costs)
    lookup = {(int(r.foreground_index), int(r.background_gene_index)): r for r in edge_df.itertuples(index=False)}
    assignments = {int(r_i): lookup[(int(r_i), int(c_j))] for r_i, c_j in zip(rows, cols) if c_j < g_count}
    assert len(assignments) == optimum

    out_rows: list[dict] = []
    for i, row in enumerate(fg.itertuples(index=False)):
        rec = {
            "foreground_event": row.event_id,
            "foreground_gene": row.gene_id,
            "event_type": row.event_type,
            "direction": getattr(row, "direction", ""),
            "dpsi": float(getattr(row, "dpsi", np.nan)),
            "status": "unmatched",
            "background_event": "",
            "background_gene": "",
            "candidate_background_events": int(row.candidate_background_events),
            "match_cost": np.nan,
        }
        if i in assignments:
            edge = assignments[i]
            ctrl = bg.iloc[int(edge.background_index)]
            rec.update(
                status="matched",
                background_event=ctrl.event_id,
                background_gene=ctrl.gene_id,
                match_cost=float(edge.cost),
                reason="globally_assigned_v2",
            )
            for dc in dist_cols:
                rec[dc] = float(getattr(edge, dc))
        else:
            rec["reason"] = "no_compatible_background" if row.candidate_background_events == 0 else "background_gene_capacity_exhausted"
        out_rows.append(rec)

    res = pd.DataFrame(out_rows)
    matched = res[res["status"].eq("matched")]
    assert not matched["background_gene"].duplicated().any()
    cert = {
        "maximum_cardinality_certificate": optimum,
        "candidate_gene_edges": int(len(edge_df)),
        "used_flank_gc_caliper": bool(use_flank),
        "flank_gc_caliper": float(flank_gc_caliper) if use_flank else np.nan,
        "total_matched_cost": float(matched["match_cost"].sum()) if len(matched) else 0.0,
    }
    return res, cert, edge_df
