import pandas as pd

from katana_next.science.cell_embedding import build_cell_embedding, cell_feature_values
from katana_next.science.blind_clustering import cluster_latent, reveal_partition
from katana_next.science.single_cell import load_single_cell_counts, load_single_cell_metadata
from katana_next.workflows.base import ArtifactWriter, ParameterDefinition, WorkflowBase, WorkflowDefinition
from katana_next.workflows.figure_style import readable_plot as _plot


class CellEmbeddingWorkflow(WorkflowBase):
    version = "1.3.0"
    definition = WorkflowDefinition(
        "single_cell.cell_embedding", "single-cell", "Single-cell add-on: expression / event UMAP",
        "Do cells share expression or splice-event patterns, and do they separate by the supplied cell types?",
        "Label-blind variable features, measured-event median imputation, scaling, PCA and seeded UMAP. Cell-type prediction is evaluated separately with training-only features and held-out donors.",
        "Cell-type/donor/event-support/missingness UMAPs, event distributions, coverage audit and held-out-donor prediction results.",
        "cell for descriptive visualization; donor held out for prediction",
        required_asset_kinds=("single_cell_metadata",),
        optional_asset_kinds=("single_cell_junction_counts", "single_cell_expression_counts"),
        parameters=(
            ParameterDefinition("mode", "Feature source", "select", "splicing", "Splicing uses inclusion/exclusion event counts; expression uses long-form cell_id/gene_id/count.", ("splicing", "expression")),
            ParameterDefinition("min_reads", "Reads per cell-event", "number", 10, "Below this depth PSI is missing, not zero; this is a visualization gate, not a DAS test.", minimum=1, step=1),
            ParameterDefinition("min_feature_fraction", "Cells measuring a feature", "number", .2, "Label-blind feature coverage; all exclusions are retained.", minimum=.01, maximum=1, step=.05),
            ParameterDefinition("min_observed_features", "Measured features per cell", "number", 5, "Cells without enough features are reported separately.", minimum=1, step=1),
            ParameterDefinition("max_features", "Maximum variable features", "number", 500, "Variance-selected features; no DAS significance or cell-type labels are used.", minimum=2, maximum=3000, step=50),
            ParameterDefinition("n_neighbors", "UMAP / prediction neighbors", "number", 15, "Fixed before looking at cell-type separation.", minimum=2, maximum=100, step=1),
            ParameterDefinition("min_dist", "UMAP point spacing", "number", .3, "Lower values pack points more tightly; this changes the display, not the Leiden cluster resolution.", minimum=0, maximum=1, step=.05),
            ParameterDefinition("pca_components", "Pre-UMAP PCA dimensions", "number", 30, "Used for UMAP and optional cell clustering; capped by available features/cells.", minimum=2, maximum=100, step=1),
            ParameterDefinition("cluster_cells", "Label-blind Leiden clustering", "boolean", False, "Cluster pre-UMAP PCA features without cell-type or donor labels; annotations are revealed only for separate evaluation."),
            ParameterDefinition("clustering_neighbors", "Leiden graph neighbors", "number", 15, "The cell graph is built in pre-UMAP PCA space, not from the two-dimensional map.", minimum=2, maximum=100, step=1),
            ParameterDefinition("clustering_resolution", "Leiden cluster resolution", "number", 1., "Higher resolution usually produces more groups. Declare it before revealing annotations; do not tune to the desired number of cell types.", minimum=.01, maximum=10, step=.1),
            ParameterDefinition("random_seed", "Random seed", "number", 17, "Reproducible descriptive embedding.", minimum=0, step=1),
            ParameterDefinition("gene_filter", "Optional splicing gene ID / name", "text", "", "Restrict the splicing feature matrix to one gene; at least two variable events are needed. For expression, use the overlay option instead."),
            ParameterDefinition("expression_overlay_gene", "Color by one gene's expression", "text", "", "Optional expression-count overlay on the same coordinates; requires a complete single_cell_expression_counts matrix. This does not recluster cells or select features."),
            ParameterDefinition("evaluate_cell_types", "Evaluate held-out donors", "boolean", True, "Predictions use the feature matrix, not UMAP coordinates. Labels never select embedding features."),
        ),
        limitations="UMAP is not a significance test. Missingness or sequencing depth can separate cells. Donor-held-out prediction evaluates supplied labels, not independent truth; expression-derived labels favor the expression baseline.",
    )

    def input_paths(self, context, parameters):
        expected = "single_cell_junction_counts" if parameters.get("mode", "splicing") == "splicing" else "single_cell_expression_counts"
        if not any(a.kind == expected and a.path.is_file() for a in context.manifest.assets):
            raise ValueError(f"Register a {expected} input; no expression-to-splicing substitution")
        if parameters.get("mode", "splicing") == "expression" and parameters.get("gene_filter"):
            raise ValueError("Use expression_overlay_gene to show one gene on expression UMAP; do not generate a 2-D map from one expression feature")
        if parameters.get("expression_overlay_gene"):
            expression_assets = [a for a in context.manifest.assets if a.kind == "single_cell_expression_counts" and a.path.is_file()]
            if len(expression_assets) != 1:
                raise ValueError("A gene expression overlay requires exactly one single_cell_expression_counts matrix; intron/event support is not gene expression")
        return super().input_paths(context, parameters)

    def run(self, context, job, run_fingerprint, parameters):
        metadata = load_single_cell_metadata(context)
        mode = parameters.get("mode", "splicing")
        if mode == "splicing":
            counts = load_single_cell_counts(context)
        else:
            assets = [a for a in context.manifest.assets if a.kind == "single_cell_expression_counts"]
            if len(assets) != 1:
                raise ValueError("Declare exactly one expression-count matrix")
            path = assets[0].path
            counts = pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path, sep="\t")
        job.report(.1, "Selecting label-blind cell features; retaining missingness and donor identity")
        result = build_cell_embedding(counts, metadata, **{k: parameters[k] for k in (
            "mode", "min_reads", "min_feature_fraction", "min_observed_features", "max_features", "n_neighbors",
            "min_dist", "pca_components", "random_seed", "gene_filter", "evaluate_cell_types") if k in parameters})
        cluster_tables = []
        if parameters.get("cluster_cells", False):
            assignments, metrics = cluster_latent(result.latent_coordinates,
                n_neighbors=parameters.get("clustering_neighbors", 15),
                resolution=parameters.get("clustering_resolution", 1.), random_seed=parameters.get("random_seed", 17))
            result.coordinates["cluster_id"] = result.coordinates.cell_id.map(assignments.set_index("cell_id").cluster_id)
            cluster_tables = [("Label-blind Leiden assignments", "leiden_assignments", assignments),
                ("Label-free Leiden diagnostics", "leiden_diagnostics", pd.DataFrame([metrics]))]
            if parameters.get("evaluate_cell_types", True):
                counts, agreement = reveal_partition(assignments, result.coordinates)
                cluster_tables.extend([("Post-clustering annotation composition", "cluster_cell_types", counts),
                    ("Post-clustering annotation agreement", "cluster_agreement", pd.DataFrame([agreement]))])
        overlay_gene = parameters.get("expression_overlay_gene", "").strip()
        if overlay_gene:
            asset = next(a for a in context.manifest.assets if a.kind == "single_cell_expression_counts")
            expression = pd.read_parquet(asset.path) if asset.path.suffix == ".parquet" else pd.read_csv(asset.path, sep="\t")
            overlay = expression_gene_overlay(expression, result.coordinates.cell_id, overlay_gene)
            result.coordinates["gene_expression_logCP10K"] = overlay.to_numpy()
            result.coordinates["expression_overlay_gene"] = overlay_gene
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version,
            common_metadata={"method": self.definition.method, "limitations": self.definition.limitations,
                             "significance_test": "none; UMAP is descriptive", "feature_selection": "no cell-type labels or DAS filtering",
                             "support_metric": "sum of per-event evidence; overlapping events can share reads, not unique library depth"})
        for label, name, frame in (
            ("Cell embedding coverage audit", "cells", result.cells),
            ("Cell embedding feature audit", "features", result.features),
            ("Cell-level UMAP coordinates", "coordinates", result.coordinates),
            ("Pre-UMAP PCA coordinates", "latent_coordinates", result.latent_coordinates),
            ("Held-out-donor cell-type prediction", "donor_validation", result.donor_validation),
            ("Embedding depth and missingness diagnostics", "diagnostics", result.diagnostics),
            ("Event / expression distributions", "feature_distributions", result.distributions),
        ):
            writer.table(label, name + ".parquet", frame, entity_type="cell_embedding")
        for label, name, frame in cluster_tables:
            writer.table(label, name + ".parquet", frame, entity_type="cell_clustering")
        if parameters.get("cluster_cells", False) and parameters.get("evaluate_cell_types", True):
            composition = next(frame for _, name, frame in cluster_tables if name == "cluster_cell_types")
            writer.plot("Blind clusters: cell types revealed afterward", "cluster_cell_types.vl.json", _plot(
                composition.to_dict("records"), "rect", {
                    "x": {"field": "cell_type", "type": "nominal", "title": "Provided cell type (evaluation only)"},
                    "y": {"field": "cluster_id", "type": "nominal", "title": "Label-blind Leiden cluster"},
                    "color": {"field": "fraction_of_cluster", "type": "quantitative", "title": "Fraction of this cluster", "scale": {"scheme": "blues", "domain": [0, 1]}},
                    "tooltip": [{"field": c} for c in ("cluster_id", "cell_type", "cells", "fraction_of_cluster")],
                }, title="What cell types make up each blind cluster?", description="Labels do not enter feature selection, graph construction or clustering. Agreement with supplied annotations is not independent biological truth.", height=380), entity_type="cell_clustering_plot")
        colors = ["cell_type", "donor_id", "condition", "total_evidence_counts", "missing_feature_fraction"]
        if parameters.get("cluster_cells", False):
            colors.insert(0, "cluster_id")
        if overlay_gene:
            colors.append("gene_expression_logCP10K")
        for color in colors:
            continuous = color in {"total_evidence_counts", "missing_feature_fraction", "gene_expression_logCP10K"}
            color_label = "Summed event support (not unique reads)" if color == "total_evidence_counts" else color.replace("_", " ")
            if color == "gene_expression_logCP10K":
                color_label = f"{overlay_gene}: log1p(CP10K)"
            writer.plot(f"{mode.title()} UMAP: {color}", f"umap_{color}.vl.json", _plot(result.coordinates.to_dict("records"),
                {"type": "circle", "size": 18, "opacity": .75}, {
                    "x": {"field": "umap_1", "type": "quantitative", "title": "UMAP 1", "scale": {"zero": False}},
                    "y": {"field": "umap_2", "type": "quantitative", "title": "UMAP 2", "scale": {"zero": False}},
                    "color": {"field": color, "type": "quantitative" if continuous else "nominal", "title": color_label,
                              "scale": {"scheme": "viridis" if continuous else "tableau20"}},
                    "tooltip": [{"field": c} for c in ("cell_id", "donor_id", "cell_type", "condition", "observed_features", "missing_feature_fraction")],
                }, title=f"Cell-level {mode} UMAP: {color_label}", description=self.definition.limitations, height=430), entity_type="cell_embedding_plot")
        writer.plot("Measured feature distributions by cell type", "feature_distributions.vl.json", _plot(
            result.distributions.loc[result.distributions.value.notna()].to_dict("records"), "boxplot", {
                "x": {"field": "cell_type", "type": "nominal", "title": "Provided cell type"},
                "y": {"field": "value", "type": "quantitative", "title": "Posterior PSI" if mode == "splicing" else "log1p(CP10K)"},
                "facet": {"field": "feature_id", "type": "nominal", "columns": 2},
                "color": {"field": "cell_type", "type": "nominal", "scale": {"scheme": "tableau20"}},
            }, title="Measured feature distributions (missing PSI excluded)", description="First eight variance-ranked features, not selected by significance or cell-type separation.", height=280), entity_type="cell_distribution")
        return writer.finish()


def expression_gene_overlay(expression, cells, gene):
    values = cell_feature_values(expression, mode="expression")
    chosen = values.feature_id.eq(gene)
    if "gene_symbol" in values:
        chosen |= values.gene_symbol.eq(gene)
    identifiers = values.loc[chosen, "feature_id"].unique()
    if len(identifiers) != 1:
        raise ValueError("Gene expression overlay needs one unambiguous gene present in the count matrix")
    totals = values.groupby("cell_id").depth.sum()
    gene_values = values.loc[chosen].set_index("cell_id").value
    result = cells.map(gene_values)
    measured = cells.map(totals).gt(0)
    # Sparse zeros are defined only for cells with an observed expression library.
    return result.mask(measured & result.isna(), 0.)
