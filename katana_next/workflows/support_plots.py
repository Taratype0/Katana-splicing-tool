from __future__ import annotations

from itertools import combinations
import hashlib

import numpy as np
import pandas as pd

from katana_next.science.enrichment import overrepresentation
from katana_next.science.statistics import bh_adjust
from katana_next.workflows.bulk import _plot


def rbp_support_plots(writer, candidates, matched, selected_genes, gene_sets, labels, parent_ids, enrichment=None):
    if "comparison_id" not in candidates:
        return
    rows = []
    for comparison, frame in ([] if enrichment is not None else candidates.groupby("comparison_id", sort=False)):
        eligible = frame["selected_input"] if "selected_input" in frame else frame["canonical_gene_id"].isin(selected_genes)
        if "selected_input" not in frame and "evidence_dimensions" in frame:
            eligible &= frame["evidence_dimensions"].fillna(0).gt(0)
        tested = set(frame["canonical_gene_id"]) - {""}
        selected = set(frame.loc[eligible, "canonical_gene_id"]) - {""}
        result = overrepresentation(selected, tested, gene_sets, minimum_overlap=1)
        result["comparison_id"] = comparison
        result["comparison_label"] = labels.get(comparison, comparison)
        rows.append(result)
    enrichment = enrichment.copy() if enrichment is not None else (pd.concat(rows, ignore_index=True) if rows else pd.DataFrame())
    if not enrichment.empty:
        enrichment["fdr"] = bh_adjust(enrichment["pvalue"])
        enrichment["fdr_family"] = "all eligible comparison-by-RBP-set tests in this run"
        enrichment["minus_log10_fdr"] = -np.log10(enrichment.fdr.clip(lower=np.finfo(float).tiny))
        enrichment["row_id"] = enrichment.comparison_id + "::" + enrichment.term
        enrichment["comparison_label"] = enrichment.comparison_id.map(labels).fillna(enrichment.comparison_id)
    table = writer.table("Cross RBP enrichment support", "cross_rbp_enrichment_support.parquet", enrichment, entity_type="rbp_gene_set", parents=parent_ids)
    _heatmap(writer, "Cross RBP enrichment support", "cross_rbp_enrichment_support", enrichment, "term", "minus_log10_fdr", "-log10 FDR", (table.artifact_id,), size="overlap_genes")
    active = matched.copy()
    if "selected_input" not in active and "evidence_dimensions" in active:
        active = active[active.evidence_dimensions.fillna(0).gt(0)]
    active["comparison_label"] = active.comparison_id.map(labels).fillna(active.comparison_id)
    active["gene_label"] = active.get("gene_symbol", active.canonical_gene_id).replace("", np.nan).fillna(active.canonical_gene_id)
    active["rank"] = pd.to_numeric(active.get("rank_within_comparison", pd.Series(np.nan, index=active.index)), errors="coerce")
    matrix = writer.table("Cross RBP matched gene matrix", "cross_rbp_matched_gene_matrix.parquet", active, entity_type="gene", parents=parent_ids)
    _heatmap(writer, "RBP candidate membership and rank", "cross_rbp_matched_gene_matrix", active.assign(present=1), "gene_label", "present", "Candidate membership", (matrix.artifact_id,))
    pairwise_support_plots(writer, active, "canonical_gene_id", "rank", "RBP candidate rank", "pairwise_rbp_rank_support", labels, (matrix.artifact_id,), comparisons=list(candidates.comparison_id.unique()))


def isoform_support_plots(writer, contrasts, switches, labels, parent_ids):
    if contrasts.empty:
        return
    aggregations = {"delta_tpm": ("delta_tpm", "sum"), "gene_symbol": ("gene_symbol", "first")}
    if "best_dtu_fdr" in contrasts:
        aggregations["best_dtu_fdr"] = ("best_dtu_fdr", "min")
    genes = contrasts.groupby(["comparison_id", "gene_id"], as_index=False).agg(**aggregations)
    genes["row_id"] = genes.comparison_id + "::" + genes.gene_id
    genes["comparison_label"] = genes.comparison_id.map(labels).fillna(genes.comparison_id)
    genes["gene_label"] = genes.gene_symbol.replace("", np.nan).fillna(genes.gene_id)
    genes["dtu_supported"] = pd.to_numeric(genes.get("best_dtu_fdr", pd.Series(np.nan, index=genes.index)), errors="coerce").le(.05)
    table = writer.table("Cross isoform gene support matrix", "cross_isoform_gene_support_matrix.parquet", genes, entity_type="gene", parents=parent_ids)
    _heatmap(writer, "Gene TPM change and DTU support", "cross_isoform_gene_support_heatmap", genes, "gene_label", "delta_tpm", "Total gene TPM difference", (table.artifact_id,))
    switches = switches.copy()
    switches["comparison_label"] = switches.comparison_id.map(labels).fillna(switches.comparison_id)
    switches["gene_label"] = switches.gene_symbol.replace("", np.nan).fillna(switches.gene_id)
    _heatmap(writer, "Descriptive dominant-transcript switches", "cross_isoform_dominant_isoform_switches", switches.assign(switched=1), "gene_label", "switched", "Mean-TPM dominant transcript changed", parent_ids)
    pairwise_support_plots(writer, contrasts, "transcript_id", "delta_tpm", "Transcript delta TPM", "pairwise_isoform_delta_tpm", labels, parent_ids)


def _heatmap(writer, name, stem, frame, row_key, value_key, value_title, parents, size=None):
    values = frame.copy()
    # Only the view is capped; source tables retain all genes and terms.
    if not values.empty:
        ordered = values.assign(magnitude=pd.to_numeric(values[value_key], errors="coerce").abs()).sort_values("magnitude", ascending=False)
        values = values[values[row_key].isin(ordered[row_key].drop_duplicates().head(50))]
    encoding = {"x": {"field": "comparison_label", "type": "nominal", "axis": {"title": None, "labelAngle": -90, "labelLimit": 500}},
                "y": {"field": row_key, "type": "nominal", "axis": {"title": None, "labelLimit": 350}},
                "color": {"field": value_key, "type": "quantitative", "legend": {"title": value_title}},
                "tooltip": [{"field": c} for c in values.columns if c not in {"row_id", "source_path"}]}
    if size:
        encoding["size"] = {"field": size, "type": "quantitative", "legend": {"title": "Overlapping genes"}}
    writer.plot(name, f"{stem}.vl.json", _plot(values.to_dict("records"), "circle" if size else "rect", encoding,
                title=name, description="Support-prioritization only, not a causal model. At most 50 rows displayed; all rows are exported.",
                height=max(250, min(1100, 21 * values[row_key].nunique()))), entity_type="support_plot", parents=parents)


def pairwise_support_plots(writer, frame, identity, value, title, stem, labels, parents, comparisons=None):
    comparisons = comparisons or list(frame.comparison_id.drop_duplicates())
    for left, right in combinations(comparisons, 2):
        a = frame[frame.comparison_id.eq(left)].drop_duplicates(identity).set_index(identity)[[value]].rename(columns={value: "value_a"})
        b = frame[frame.comparison_id.eq(right)].drop_duplicates(identity).set_index(identity)[[value]].rename(columns={value: "value_b"})
        joined = a.join(b, how="outer").reset_index()
        joined["row_id"] = joined[identity].astype(str)
        joined["present_a"] = joined[identity].isin(a.index)
        joined["present_b"] = joined[identity].isin(b.index)
        joined["membership"] = np.select([joined.present_a & joined.present_b, joined.present_a],
                                         ["shared", f"only {labels.get(left, left)}"], default=f"only {labels.get(right, right)}")
        joined["comparison_a"], joined["comparison_b"] = left, right
        suffix = hashlib.sha256(f"{left}\0{right}".encode()).hexdigest()[:16]
        table = writer.table(f"{title}: {labels.get(left, left)} compared with {labels.get(right, right)}", f"{stem}_{suffix}.parquet", joined, entity_type="support_pair", parents=parents)
        measurable = joined.dropna(subset=["value_a", "value_b"])
        plotted, mark, encoding, mode = pairwise_display(measurable, identity)
        encoding["x"]["axis"] = {"title": f"{labels.get(left, left)}: {title}"}
        encoding["y"]["axis"] = {"title": f"{labels.get(right, right)}: {title}"}
        writer.plot(table.name, f"{stem}_{suffix}.vl.json", _plot(plotted, mark, encoding,
                     title=f"{title}: comparison A versus comparison B",
                     description=f"Support-prioritization only. {len(measurable)} measured on both sides; {mode}. {len(joined)-len(measurable)} rows lacking paired measurements retained in the table, not plotted as zero.", height=460),
                    metadata={"display_mode": mode, "measured_pair_count": len(measurable), "display_cells": len(plotted)},
                    entity_type="support_plot", parents=(table.artifact_id,))


def pairwise_display(measurable, identity, threshold=5000):
    if len(measurable) <= threshold:
        return (measurable[[identity, "value_a", "value_b", "membership"]].to_dict("records"),
                {"type": "point", "filled": True, "opacity": .65},
                {"x": {"field": "value_a", "type": "quantitative"},
                 "y": {"field": "value_b", "type": "quantitative"},
                 "tooltip": [{"field": c} for c in (identity, "value_a", "value_b", "membership")]},
                "individual scatter points")
    counts, xs, ys = np.histogram2d(measurable.value_a, measurable.value_b, bins=60)
    rows = [{"x0": xs[i], "x1": xs[i+1], "y0": ys[j], "y1": ys[j+1], "count": int(counts[i, j])}
            for i, j in zip(*np.nonzero(counts))]
    return (rows, "rect",
            {"x": {"field": "x0", "type": "quantitative"}, "x2": {"field": "x1"},
             "y": {"field": "y0", "type": "quantitative"}, "y2": {"field": "y1"},
             "color": {"field": "count", "type": "quantitative", "scale": {"type": "log", "scheme": "tealblues"},
                       "legend": {"title": "Transcript count"}},
             "tooltip": [{"field": c} for c in ("x0", "x1", "y0", "y1", "count")]},
            "60-by-60 density bins containing every measured pair; full individual values in TSV")
