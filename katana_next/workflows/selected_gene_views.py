"""Selection-preserving candidate heatmaps and gene cards, using saved upstream results."""
from pathlib import Path
import hashlib
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist

from katana_next.workflows.base import WorkflowBase, WorkflowDefinition, ParameterDefinition, ArtifactWriter
from katana_next.workflows.sources import selected_source, selection_rows, resolve_table
from katana_next.workflows.figure_exports import _RENDER_LOCK


def cluster_order(matrix):
    if len(matrix) < 2:
        return list(range(len(matrix)))
    # Imputation is for ordering only; the displayed/exported PSI stays missing.
    values = matrix.to_numpy(float)
    means = np.nanmean(values, axis=1, keepdims=True)
    filled = np.where(np.isfinite(values), values, np.nan_to_num(means, nan=0.0))
    distance = pdist(filled, metric="euclidean")
    return leaves_list(linkage(distance, method="average", optimal_ordering=True)).tolist()


def save_figure(writer, figure, stem, name, parents):
    exports = {}
    directory = writer.directory / "plots"
    directory.mkdir(exist_ok=True)
    for suffix in ("png", "svg"):
        path = directory / f"{stem}.{suffix}"
        figure.savefig(path, dpi=160, bbox_inches="tight", facecolor="white")
        exports[suffix] = str(path)
    figure.clear()
    writer.existing_file(name, Path(exports["png"]), entity_type="gene_plot",
                         metadata={"mime_type": "image/png", "exports": exports}, parents=parents)


class SelectedGeneViewsWorkflow(WorkflowBase):
    version = "2.1.4"
    definition = WorkflowDefinition(
        workflow_id="prioritize.selected_gene_views", workspace_id="prioritize",
        label="Selected-gene AS heatmap and gene cards",
        question="How do all selected genes splice across samples, and how does their expression compare?",
        method="Saved PSI matrix with average-linkage Euclidean row/column clustering; gene cards use cached sample TPM and rMATS dPSI.",
        output_summary="Selection coverage audit, clustered PSI/order tables, paginated heatmap and one expression/AS card per selected gene.",
        analysis_unit="AS event across biological samples; gene for expression cards",
        supports_contrast_selection=True,
        parameters=(ParameterDefinition("selection_id", "Saved gene selection", "selection", "", "Select rows from the candidate table; every selected gene is audited, including genes with no measured AS."),
                    ParameterDefinition("events_per_gene", "Events per selected gene", "number", 3, "Highest PSI variance per gene; 0 includes every measured event. No global Top-N gene truncation.", minimum=0, maximum=100, step=1),
                    ParameterDefinition("rows_per_page", "Heatmap rows per page", "number", 40, "More events create additional pages, never dropped rows.", minimum=5, maximum=100, step=5)),
        limitations="Clustering orders observed PSI patterns, not new biological classes. TPM is descriptive; no-data cells stay blank. Requires explicitly completed upstream workflows.",
    )

    def parent_artifacts(self, context, parameters):
        source = selected_source(context, parameters)
        if source is None:
            raise ValueError("Create a saved gene selection first; planning never runs candidate selection.")
        values = [source, resolve_table(context, "differential.as", "Differential AS events", parameters),
                  resolve_table(context, "landscape.bulk", "Sample-level PSI matrix", parameters),
                  resolve_table(context, "prioritize.isoform_support", "Sample gene abundance", parameters)]
        if values[1] is None or values[2] is None:
            raise ValueError("Run DAS and Landscape explicitly first; gene views do not launch them.")
        return [item for item in values if item]

    def input_paths(self, context, parameters):
        return [a.path for a in self.parent_artifacts(context, parameters)]

    def run(self, context, job, run_fingerprint, parameters):
        source = selected_source(context, parameters)
        if not source:
            raise ValueError("Create a saved candidate gene selection first.")
        selected = selection_rows(context, parameters, source, pd.read_parquet(source.path))
        if "gene_id" not in selected:
            raise ValueError("Selected rows require canonical gene_id annotation.")
        selected["gene_symbol"] = selected.get("gene_symbol", selected.gene_id).fillna(selected.gene_id)
        selected = selected.drop_duplicates("gene_id")
        das = resolve_table(context, "differential.as", "Differential AS events", parameters)
        psi_source = resolve_table(context, "landscape.bulk", "Sample-level PSI matrix", parameters)
        if not das or not psi_source:
            raise ValueError("Run DAS and Landscape explicitly first; gene views do not launch them.")
        events = pd.read_parquet(das.path)
        ids = parameters.get("contrast_ids") or [c.contrast_id for c in context.manifest.contrasts]
        events = events[events.comparison_id.isin(ids) & events.gene_id.isin(selected.gene_id)]
        psi = pd.read_parquet(psi_source.path).set_index("event_id").drop(columns="row_id", errors="ignore")
        samples = [s.sample_id for s in context.manifest.samples if any(s.sample_id in (*c.numerator.sample_ids, *c.denominator.sample_ids) for c in context.manifest.contrasts if c.contrast_id in ids)]
        missing_samples = set(samples) - set(psi.columns)
        if missing_samples:
            raise ValueError(f"Cached PSI matrix misses selected samples: {sorted(missing_samples)}")
        psi = psi[samples]
        candidate_events = events.drop_duplicates("event_id").set_index("event_id")
        candidate_events = candidate_events[candidate_events.index.isin(psi.index)]
        candidate_events["variance"] = psi.var(axis=1)
        ordered = candidate_events.sort_values(["variance"], ascending=False)
        limit = int(parameters["events_per_gene"])
        if limit:
            ordered = ordered.groupby("gene_id", sort=False).head(limit)
        matrix = psi.reindex(ordered.index).dropna(how="all")
        audit = selected[["gene_id", "gene_symbol"]].copy()
        audit["measured_AS_events"] = audit.gene_id.map(candidate_events.gene_id.value_counts()).fillna(0).astype(int)
        audit["displayed_AS_events"] = audit.gene_id.map(ordered.loc[matrix.index].gene_id.value_counts()).fillna(0).astype(int)
        audit["status"] = np.where(audit.displayed_AS_events.gt(0), "displayed", "no measured AS in cached PSI; gene retained in audit/card")
        audit["row_id"] = audit.gene_id
        parents = tuple(a.artifact_id for a in self.parent_artifacts(context, parameters))
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version,
                                common_metadata={"method": self.definition.method, "limitations": self.definition.limitations,
                                                 "selected_unique_genes": len(selected), "missing_PSI": "gray; row-mean fill for clustering only",
                                                 "heatmap_display":"row z-score; ddof=1; pastel blue-white-red; limits +/-2.5; raw PSI retained in matrix"})
        writer.table("Selected gene coverage audit", "selected_gene_audit.parquet", audit, entity_type="gene", parents=parents)
        writer.table("Selected-gene AS event detail", "selected_gene_events.parquet", events, entity_type="splice_event", parents=parents)
        if not matrix.empty:
            matrix = matrix.iloc[cluster_order(matrix), cluster_order(matrix.T)]
            display_ids = {event: f"E{i + 1:03d}" for i, event in enumerate(matrix.index)}
            export = matrix.reset_index()
            export["row_id"] = export.event_id
            export["display_event_id"] = export.event_id.map(display_ids)
            writer.table("Clustered selected-gene PSI matrix", "selected_gene_psi.parquet", export, entity_type="splice_event", parents=parents)
            order = pd.DataFrame([{"row_id": f"{axis}::{identifier}", "axis": axis, "identifier": identifier, "cluster_order": i}
                                  for axis, values in (("event", matrix.index), ("sample", matrix.columns)) for i, identifier in enumerate(values)])
            writer.table("Heatmap clustering order", "clustering_order.parquet", order, entity_type="cluster_order", parents=parents)
            step = int(parameters["rows_per_page"])
            with _RENDER_LOCK:
                for start in range(0, len(matrix), step):
                    chunk = matrix.iloc[start:start + step]
                    figure = Figure(figsize=(max(9, .38 * len(samples) + 5), max(5, .24 * len(chunk) + 3)), layout="constrained")
                    axis = figure.subplots()
                    from matplotlib import colormaps
                    from katana_next.science.heatmap_style import row_scaled, pastel_diverging
                    art = axis.imshow(np.ma.masked_invalid(row_scaled(chunk.to_numpy(float))), cmap=pastel_diverging(), vmin=-2.5, vmax=2.5, aspect="auto")
                    axis.set_xticks(range(len(chunk.columns)), chunk.columns, rotation=90, fontsize=8)
                    names = [f"{ordered.loc[e, 'gene_symbol']} | {ordered.loc[e, 'event_type']} | {display_ids[e]}" for e in chunk.index]
                    axis.set_yticks(range(len(chunk)), names, fontsize=7)
                    axis.set_title(f"Selected-gene AS | row-scaled PSI | page {start // step + 1}")
                    figure.colorbar(art, ax=axis, label="Row z score (+/-2.5 display limits); gray = missing")
                    save_figure(writer, figure, f"selected_gene_heatmap_{start // step + 1:03d}", f"Clustered PSI heatmap page {start // step + 1}", parents)
                if len(matrix) > 1:
                    from scipy.cluster.hierarchy import dendrogram
                    filled = matrix.T.fillna(matrix.mean(axis=1)).T.fillna(0)
                    tree = linkage(pdist(filled), method="average", optimal_ordering=True)
                    figure = Figure(figsize=(12, 5), layout="constrained")
                    axis = figure.subplots()
                    dendrogram(tree, labels=[display_ids[e] for e in matrix.index], ax=axis,
                               truncate_mode="lastp", p=min(40, len(matrix)), leaf_rotation=90, leaf_font_size=8)
                    axis.set_title("Event clustering dendrogram | average linkage, Euclidean PSI distance")
                    axis.set_ylabel("PSI distance")
                    axis.set_xlabel("E codes map to exact event IDs in the PSI table; parenthesized labels count collapsed leaves")
                    save_figure(writer, figure, "selected_gene_dendrogram", "AS event clustering dendrogram", parents)
        expression = resolve_table(context, "prioritize.isoform_support", "Sample gene abundance", parameters)
        abundance = pd.read_parquet(expression.path) if expression else pd.DataFrame()
        if not abundance.empty:
            abundance = abundance[abundance.sample_id.isin(samples) & abundance.gene_id.isin(selected.gene_id)]
            writer.table("Selected gene expression", "selected_gene_expression.parquet", abundance, entity_type="gene", parents=parents)
        with _RENDER_LOCK:
            for i, gene in enumerate(selected.itertuples()):
                job.report(.2 + .7 * i / max(1, len(selected)), f"Gene card {i + 1}/{len(selected)}: {gene.gene_symbol}")
                figure = Figure(figsize=(11, 4.8), layout="constrained")
                left, right = figure.subplots(1, 2)
                part = abundance[abundance.gene_id.eq(gene.gene_id)] if not abundance.empty else pd.DataFrame()
                if part.empty:
                    left.text(.5, .5, "No cached sample TPM for this gene\nRun Isoform support explicitly", ha="center", transform=left.transAxes)
                else:
                    means = part.groupby("condition", sort=False).total_tpm.mean()
                    left.bar(range(len(means)), means, color="#0072B2")
                    for j, condition in enumerate(means.index):
                        values = part.loc[part.condition.eq(condition), "total_tpm"]
                        left.scatter(j + np.linspace(-.10, .10, len(values)), values, c="black", s=12, zorder=3)
                    left.set_xticks(range(len(means)), means.index, rotation=45, ha="right")
                left.set_title("Expression by experimental condition")
                left.set_ylabel("Total gene TPM; dots = samples")
                gene_events = events[events.gene_id.eq(gene.gene_id) & events.stable_significant.fillna(False)].copy()
                gene_events["abs_dpsi"] = gene_events.dpsi.abs()
                dominant = gene_events.sort_values("abs_dpsi", ascending=False).drop_duplicates("comparison_id").set_index("comparison_id")
                values = [dominant.loc[c, "dpsi"] if c in dominant.index else np.nan for c in ids]
                right.bar(range(len(ids)), values, color=["#009E73" if v >= 0 else "#D55E00" for v in values])
                for j, value in enumerate(values):
                    right.text(j, value if np.isfinite(value) else 0, f"{value:+.3f}" if np.isfinite(value) else "no strict AS", ha="center", va="bottom", fontsize=8)
                labels = {c.contrast_id: c.label for c in context.manifest.contrasts}
                right.set_xticks(range(len(ids)), [labels[c] for c in ids], rotation=45, ha="right")
                right.set_ylabel("Standardized dPSI")
                right.set_title("Largest |dPSI| among strict AS events")
                right.axhline(0, color="#68727a", linewidth=.7)
                for axis in (left, right):
                    axis.spines[["top", "right"]].set_visible(False)
                figure.suptitle(f"{gene.gene_symbol} | {gene.gene_id}")
                save_figure(writer, figure, "gene_card_" + hashlib.sha256(str(gene.gene_id).encode()).hexdigest()[:12], f"Gene card: {gene.gene_symbol}", parents)
        return writer.finish()
