"""One conservative source contract for web, Python and command-line routes."""
import pyarrow.parquet as pq


SELECTION_TARGETS = {
    "prioritize.enrichment", "prioritize.rbp_support", "prioritize.selected_gene_views",
    "prioritize.event_motifs", "prioritize.binding_evidence", "validate.native_sashimi",
    "single_cell.sashimi",
    "validate.perturbation",
    "prioritize.motiftrace_consequences",
    "prioritize.binding_enrichment", "prioritize.rbp_mechanism", "prioritize.rbp_sequence_results", "prioritize.se_rbp_impact", "differential.exploratory_event_interaction",
}


def compatible_workflows(artifact):
    if artifact.kind != "table" or not artifact.path.is_file() or artifact.path.suffix != ".parquet":
        return []
    columns = set(pq.read_schema(artifact.path).names)
    routes = []
    if "row_id" not in columns:
        return routes
    if columns & {"gene_id", "gene_symbol"}:
        routes.extend(["prioritize.enrichment", "prioritize.rbp_support"])
    if "gene_id" in columns and not artifact.analysis_id.startswith("single_cell."):
        routes.append("prioritize.selected_gene_views")
        routes.append("validate.perturbation")
    if {"event_id", "event_type", "chr", "strand"} <= columns:
        routes.extend(["prioritize.event_motifs", "prioritize.binding_evidence"])
        if {"gene_id", "comparison_id"} <= columns:
            routes.append("prioritize.binding_enrichment")
            routes.append("prioritize.rbp_mechanism")
            routes.append("prioritize.rbp_sequence_results")
            routes.append("prioritize.se_rbp_impact")
            routes.append("prioritize.motiftrace_consequences")
        if {"comparison_id", "IncLevel1", "IncLevel2"} <= columns:
            routes.append("validate.native_sashimi")
    if artifact.analysis_id == "differential.as" and {"event_id", "comparison_id", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncFormLen", "SkipFormLen"} <= columns:
        routes.append("differential.exploratory_event_interaction")
    if artifact.analysis_id.startswith("single_cell.") and {"event_id", "cluster_id"} <= columns:
        routes.append("single_cell.sashimi")
    return routes


def validate_route(context, workflow_id, parameters):
    if not parameters.get("selection_id"):
        return
    from katana_next.workflows.sources import selected_source
    source = selected_source(context, parameters)
    if workflow_id not in compatible_workflows(source):
        raise ValueError(f"Selection from {source.name} cannot be sent to {workflow_id}. Keep the original gene/event/assay identity; select annotated member rows rather than a module summary.")
    # Selection-driven MotifTrace kernels use the source event comparison. An
    # explicit interface comparison must not relabel that computation. The
    # completed-map 'all frozen' option remains a separate display-only scope.
    motif_selection = workflow_id in {
        "prioritize.rbp_sequence_results", "prioritize.se_rbp_impact",
        "prioritize.motiftrace_consequences",
    } or (workflow_id == "prioritize.rbp_mechanism" and
          parameters.get("sequence_scope", "boundary_association") in {
              "boundary_association", "se_sequence_impact", "completed_sequence_maps", "event_consequences"})
    if motif_selection:
        from katana_next.workflows.sources import selection_rows
        selected = selection_rows(context, parameters, source,
            pq.read_table(source.path, columns=["row_id", "comparison_id"]).to_pandas())
        cases = selected["comparison_id"]
        if selected.empty or cases.isna().any() or cases.nunique() != 1:
            raise ValueError("MotifTrace requires events selected from exactly one declared comparison")
        case = str(cases.iloc[0])
        if case not in {c.contrast_id for c in context.manifest.contrasts}:
            raise ValueError("The selected event comparison is not declared in this project")
        explicit = parameters.get("contrast_ids")
        if explicit is not None:
            if isinstance(explicit, str):
                requested = {explicit}
            elif isinstance(explicit, (list, tuple, set)):
                requested = {str(value) for value in explicit}
            else:
                raise ValueError("contrast_ids must be a list of declared contrast IDs")
            if requested != {case}:
                raise ValueError("The selected events and explicit contrast_ids refer to different comparisons")
        elif workflow_id == "prioritize.rbp_mechanism":
            # This shared entry exposes comparison selection; record the exact
            # source scope rather than leaving an apparent all-project request.
            parameters["contrast_ids"] = [case]
