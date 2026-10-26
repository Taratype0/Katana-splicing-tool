"""Public optional analysis switches, shared by planning and reporting."""
DEFAULTS = dict.fromkeys(("landscape", "cross_comparison", "enrichment", "rbp_identity",
                          "motiftrace", "motif_binding", "mapp", "coverage", "network",
                          "interaction", "sequencing_qc"), True)
# Extra transcript-by-transcript BAM traversal is opt-in. It is not part of
# nf-core/rnaseq's default RSeQC module list. Existing explicit selections stay
# authoritative, including previously sealed plans that selected this analysis.
DEFAULTS["coverage"] = False


def analysis_selection(configuration):
    values = configuration.get("analyses", {})
    if not isinstance(values, dict) or set(values) - set(DEFAULTS):
        raise ValueError("analyses accepts: " + ", ".join(DEFAULTS))
    if any(not isinstance(value, bool) for value in values.values()):
        raise ValueError("Each analyses switch must be true or false")
    return {**DEFAULTS, **values}


def resource_kinds(configuration):
    enabled = analysis_selection(configuration)
    result = set()
    if enabled["enrichment"] or enabled["network"]:
        result.add("gene_sets_gmt")
    if enabled["rbp_identity"]:
        result.add("rbp_gene_sets")
    if enabled["motiftrace"] or enabled["motif_binding"]:
        result.update(("rbp_motif_database", "fimo_runtime"))
    if enabled["motif_binding"]:
        result.add("rbp_binding_peaks")
    return result


def component_enabled(component, enabled):
    group = {"gene_body_coverage": "coverage", "network_bootstrap": "network",
             "network_dendrogram": "network", "motif_enrichment": "motiftrace",
             "motiftrace_sequence": "motiftrace", "clip": "motif_binding",
             "mapp_native": "mapp", "mapp_parallel": "mapp"}.get(component, component)
    return enabled.get(group, True)


def active_chapters(enabled):
    chapters = {1, 3}
    if enabled["landscape"]: chapters.add(2)
    if enabled["cross_comparison"] or enabled["interaction"]: chapters.add(4)
    if enabled["enrichment"]: chapters.add(5)
    if enabled["network"]: chapters.add(6)
    if any(enabled[key] for key in ("rbp_identity", "motiftrace", "motif_binding", "mapp")): chapters.add(7)
    return chapters


def artifact_enabled(record, enabled):
    workflow = record.get("analysis_id", record.get("metadata", {}).get("analysis_id", ""))
    group = {"landscape.bulk": "landscape", "differential.cross_comparison": "cross_comparison",
             "prioritize.enrichment": "enrichment", "prioritize.rbp_support": "rbp_identity",
             "prioritize.motiftrace_consequences": "motiftrace", "prioritize.rbp_mechanism": "motif_binding",
             "prioritize.motiftrace_parallel": "mapp", "quality.read_coverage_bias": "coverage",
             "differential.replicate_interaction": "interaction"}.get(workflow)
    if workflow.startswith(("flow.network.", "programs.")): group = "network"
    if workflow.startswith("flow.mapp."): group = "mapp"
    if workflow == "flow.qc.sequencing": group = "sequencing_qc"
    if workflow.startswith(("flow.qc.genebody", "flow.reference.")): group = "coverage"
    component = record.get("component", record.get("metadata", {}).get("component", ""))
    return enabled.get(group, True) and component_enabled(component, enabled) and (
        not record.get("chapter") or int(record["chapter"]) in active_chapters(enabled))
