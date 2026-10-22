from __future__ import annotations

from katana_cosplice.architecture import ArtifactKind, ArtifactSpec, StageSpec


def _artifact(key: str, path: str, kind: ArtifactKind, required: bool = True) -> ArtifactSpec:
    return ArtifactSpec(key=key, relative_path=path, kind=kind, required=required)


STAGE_SPECS: tuple[StageSpec, ...] = (
    StageSpec(
        key="audit",
        title="1. Audit inputs",
        order=10,
        description="Validate raw input paths, sample identity, and reference availability.",
        cli_command="audit",
        accepts_force=False,
        artifacts=(
            _artifact("checks", "01_audit/input_audit.tsv", ArtifactKind.TABLE),
            _artifact("samples", "01_audit/sample_metadata_resolved.tsv", ArtifactKind.TABLE),
            _artifact("parameters", "01_audit/parameters.json", ArtifactKind.PARAMETERS),
        ),
    ),
    StageSpec(
        key="matrix",
        title="2. Build PSI matrix",
        order=20,
        description=(
            "Build sample-level splice-usage and read-count matrices from either rMATS events "
            "or annotation-independent BAM junction clusters."
        ),
        cli_command="build-matrix",
        artifacts=(
            _artifact("psi", "02_matrix/sample_level_psi_matrix.tsv", ArtifactKind.MATRIX),
            _artifact("ijc", "02_matrix/sample_level_ijc_matrix.tsv", ArtifactKind.MATRIX),
            _artifact("sjc", "02_matrix/sample_level_sjc_matrix.tsv", ArtifactKind.MATRIX),
            _artifact("events", "02_matrix/event_metadata.tsv", ArtifactKind.TABLE),
            _artifact("samples", "02_matrix/sample_metadata_resolved.tsv", ArtifactKind.TABLE),
            _artifact("sources", "02_matrix/matrix_source_selection.tsv", ArtifactKind.TABLE),
            _artifact(
                "legacy_sources",
                "02_matrix/rmats_source_selection.tsv",
                ArtifactKind.TABLE,
                False,
            ),
            _artifact("integrity", "02_matrix/matrix_integrity_summary.tsv", ArtifactKind.TABLE),
            _artifact("parameters", "02_matrix/parameters.json", ArtifactKind.PARAMETERS),
            _artifact("junction_counts", "02_matrix/junction_count_matrix.tsv", ArtifactKind.MATRIX, False),
            _artifact("junction_clusters", "02_matrix/junction_cluster_metadata.tsv", ArtifactKind.TABLE, False),
            _artifact("rmats_mapping", "02_matrix/rmats_junction_annotation.tsv", ArtifactKind.TABLE, False),
        ),
    ),
    StageSpec(
        key="stable_qc",
        title="3. Network input QC",
        order=30,
        description="Require valid nonconstant PSI and at least ten informative counts in every network sample; condition-response statistics are diagnostic only.",
        dependencies=("matrix",),
        cli_command="qc",
        artifacts=(
            _artifact("metrics", "03_stable_as_qc/stable_as_event_qc.tsv", ArtifactKind.TABLE),
            _artifact("condition_metrics", "03_stable_as_qc/stable_as_condition_qc.tsv", ArtifactKind.TABLE),
            _artifact("summary", "03_stable_as_qc/stable_as_qc_summary.tsv", ArtifactKind.TABLE),
            _artifact("matrix", "03_stable_as_qc/stable_psi_matrix.tsv", ArtifactKind.MATRIX),
            _artifact(
                "posterior_matrix",
                "03_stable_as_qc/stable_posterior_psi_matrix.tsv",
                ArtifactKind.MATRIX,
                False,
            ),
            _artifact(
                "posterior_alpha",
                "03_stable_as_qc/stable_posterior_alpha_matrix.tsv",
                ArtifactKind.MATRIX,
                False,
            ),
            _artifact(
                "posterior_beta",
                "03_stable_as_qc/stable_posterior_beta_matrix.tsv",
                ArtifactKind.MATRIX,
                False,
            ),
            _artifact(
                "uncertainty_events",
                "03_stable_as_qc/beta_binomial_event_uncertainty.tsv",
                ArtifactKind.TABLE,
                False,
            ),
            _artifact(
                "uncertainty_samples",
                "03_stable_as_qc/beta_binomial_sample_uncertainty.tsv",
                ArtifactKind.TABLE,
                False,
            ),
            _artifact("parameters", "03_stable_as_qc/parameters.json", ArtifactKind.PARAMETERS),
        ),
    ),
    StageSpec(
        key="network_observed",
        title="4A. Observed network",
        order=40,
        description="Discover co-splicing modules from the observed stable sample-level PSI matrix.",
        dependencies=("stable_qc",),
        cli_command="fit-network",
        artifacts=(
            _artifact("assignments", "04_network/observed/event_module_assignments.tsv", ArtifactKind.TABLE),
            _artifact("edges", "04_network/observed/network_edges.tsv", ArtifactKind.TABLE),
            _artifact("scores", "04_network/observed/module_program_scores.tsv", ArtifactKind.MATRIX),
            _artifact("diagnostics", "04_network/observed/network_diagnostics.tsv", ArtifactKind.TABLE),
            _artifact("sensitivity", "04_network/observed/module_resolution_sensitivity.tsv", ArtifactKind.TABLE),
            _artifact("selected", "04_network/observed/selected_event_matrix.tsv", ArtifactKind.MATRIX),
            _artifact("parameters", "04_network/observed/parameters.json", ArtifactKind.PARAMETERS),
            _artifact("posterior_correlation", "04_network/observed/posterior_ensemble_correlation.npz", ArtifactKind.OTHER, False),
            _artifact("gene_prior_matrix", "04_network/observed/gene_mantel_prior_matrix.npz", ArtifactKind.OTHER, False),
            _artifact("gene_prior_genes", "04_network/observed/gene_mantel_prior_gene_audit.tsv", ArtifactKind.TABLE, False),
            _artifact("gene_prior_events", "04_network/observed/gene_mantel_prior_event_audit.tsv", ArtifactKind.TABLE, False),
            _artifact("module_sizes_plot", "04_network/observed/plots/observed_module_sizes.png", ArtifactKind.PLOT, False),
            _artifact("score_heatmap", "04_network/observed/plots/observed_program_score_heatmap.png", ArtifactKind.PLOT, False),
            _artifact("dendrogram", "04_network/observed/plots/observed_event_dendrogram.png", ArtifactKind.PLOT, False),
        ),
    ),
    StageSpec(
        key="network_residual",
        title="4B. Residual network",
        order=41,
        description="Discover modules after removing explicitly configured sample covariates.",
        dependencies=("stable_qc",),
        cli_command="fit-network",
        artifacts=(
            _artifact("assignments", "04_network/residual/event_module_assignments.tsv", ArtifactKind.TABLE),
            _artifact("edges", "04_network/residual/network_edges.tsv", ArtifactKind.TABLE),
            _artifact("scores", "04_network/residual/module_program_scores.tsv", ArtifactKind.MATRIX),
            _artifact("diagnostics", "04_network/residual/network_diagnostics.tsv", ArtifactKind.TABLE),
            _artifact("sensitivity", "04_network/residual/module_resolution_sensitivity.tsv", ArtifactKind.TABLE),
            _artifact("selected", "04_network/residual/selected_event_matrix.tsv", ArtifactKind.MATRIX),
            _artifact("parameters", "04_network/residual/parameters.json", ArtifactKind.PARAMETERS),
            _artifact("posterior_correlation", "04_network/residual/posterior_ensemble_correlation.npz", ArtifactKind.OTHER, False),
            _artifact("gene_prior_matrix", "04_network/residual/gene_mantel_prior_matrix.npz", ArtifactKind.OTHER, False),
            _artifact("gene_prior_genes", "04_network/residual/gene_mantel_prior_gene_audit.tsv", ArtifactKind.TABLE, False),
            _artifact("gene_prior_events", "04_network/residual/gene_mantel_prior_event_audit.tsv", ArtifactKind.TABLE, False),
            _artifact("module_sizes_plot", "04_network/residual/plots/residual_module_sizes.png", ArtifactKind.PLOT, False),
            _artifact("score_heatmap", "04_network/residual/plots/residual_program_score_heatmap.png", ArtifactKind.PLOT, False),
            _artifact("dendrogram", "04_network/residual/plots/residual_event_dendrogram.png", ArtifactKind.PLOT, False),
        ),
    ),
    StageSpec(
        key="consensus_modules",
        title="4C. Bootstrap module confidence",
        order=45,
        description=(
            "Annotate the fixed WGCNA modules with sample-stratified bootstrap "
            "co-assignment, fixed-module kME, split-half confirmation and optional "
            "read-reliability confidence; this does not change membership and is "
            "not external preservation."
        ),
        dependencies=("network_observed",),
        cli_command="consensus-modules",
        artifacts=(
            _artifact("assignments", "04c_consensus_modules/consensus_event_module_assignments.tsv", ArtifactKind.TABLE),
            _artifact("modules", "04c_consensus_modules/consensus_module_summary.tsv", ArtifactKind.TABLE),
            _artifact("scores", "04c_consensus_modules/consensus_module_program_scores.tsv", ArtifactKind.MATRIX),
            _artifact("bootstrap", "04c_consensus_modules/coassignment_bootstrap_runs.tsv", ArtifactKind.TABLE),
            _artifact("edges", "04c_consensus_modules/bootstrap_coassignment_edges.tsv", ArtifactKind.TABLE),
            _artifact("matrix", "04c_consensus_modules/bootstrap_coassignment_matrix.npz", ArtifactKind.OTHER),
            _artifact("diagnostics", "04c_consensus_modules/consensus_diagnostics.tsv", ArtifactKind.TABLE),
            _artifact("parameter_search", "04c_consensus_modules/consensus_parameter_search.tsv", ArtifactKind.TABLE),
            _artifact("selection_audit", "04c_consensus_modules/multiobjective_selection_audit.tsv", ArtifactKind.TABLE),
            _artifact("parameters", "04c_consensus_modules/parameters.json", ArtifactKind.PARAMETERS),
            _artifact("heatmap", "04c_consensus_modules/plots/bootstrap_coassignment_heatmap.png", ArtifactKind.PLOT, False),
            _artifact("confidence", "04c_consensus_modules/plots/event_consensus_confidence.png", ArtifactKind.PLOT, False),
            _artifact("sizes", "04c_consensus_modules/plots/consensus_module_sizes.png", ArtifactKind.PLOT, False),
            _artifact("search", "04c_consensus_modules/plots/consensus_parameter_search.png", ArtifactKind.PLOT, False),
            _artifact("multiobjective_tradeoff", "04c_consensus_modules/plots/multiobjective_candidate_tradeoff.png", ArtifactKind.PLOT, False),
        ),
    ),
    StageSpec(
        key="hybrid_programs",
        title="Legacy: tight and diffuse programs",
        order=46,
        description=(
            "Deprecated benchmark-only compatibility stage. Production CoSplice uses the single "
            "WGCNA partition with confidence annotations from 4C."
        ),
        dependencies=("stable_qc", "network_observed"),
        optional_dependencies=("network_residual", "consensus_modules"),
        cli_command="hybrid-programs",
        visible=False,
        artifacts=(
            _artifact("events", "04d_hybrid_programs/hybrid_event_program_assignments.tsv", ArtifactKind.TABLE),
            _artifact("programs", "04d_hybrid_programs/hybrid_program_summary.tsv", ArtifactKind.TABLE),
            _artifact("scores", "04d_hybrid_programs/hybrid_program_scores.tsv", ArtifactKind.MATRIX),
            _artifact("bootstrap", "04d_hybrid_programs/diffuse_bootstrap_coassignment_runs.tsv", ArtifactKind.TABLE),
            _artifact("overlap", "04d_hybrid_programs/tight_diffuse_program_overlap.tsv", ArtifactKind.TABLE),
            _artifact("engine", "04d_hybrid_programs/diffuse_engine_audit.tsv", ArtifactKind.TABLE),
            _artifact("preservation", "04d_hybrid_programs/external_program_preservation.tsv", ArtifactKind.TABLE),
            _artifact("mapping", "04d_hybrid_programs/external_event_mapping.tsv", ArtifactKind.TABLE),
            _artifact("null", "04d_hybrid_programs/external_preservation_null.tsv", ArtifactKind.TABLE),
            _artifact("preservation_audit", "04d_hybrid_programs/external_preservation_audit.tsv", ArtifactKind.TABLE),
            _artifact("parameters", "04d_hybrid_programs/parameters.json", ArtifactKind.PARAMETERS),
            _artifact("confidence_plot", "04d_hybrid_programs/plots/tight_diffuse_program_confidence.png", ArtifactKind.PLOT, False),
            _artifact("overlap_plot", "04d_hybrid_programs/plots/tight_diffuse_program_overlap.png", ArtifactKind.PLOT, False),
            _artifact("preservation_plot", "04d_hybrid_programs/plots/external_module_preservation_zsummary.png", ArtifactKind.PLOT, False),
        ),
    ),
    StageSpec(
        key="associate_observed",
        title="5A. Observed associations",
        order=50,
        description="Associate observed-network program scores with configured sample metadata.",
        dependencies=("network_observed",),
        cli_command="associate",
        artifacts=(
            _artifact("table", "05_association/observed/module_trait_association.tsv", ArtifactKind.TABLE),
            _artifact("factorial", "05_association/observed/factorial_module_association.tsv", ArtifactKind.TABLE),
            _artifact("coefficients", "05_association/observed/factorial_module_coefficients.tsv", ArtifactKind.TABLE),
            _artifact("group_summary", "05_association/observed/module_program_score_group_summary.tsv", ArtifactKind.TABLE),
            _artifact("design", "05_association/observed/association_design.tsv", ArtifactKind.TABLE),
            _artifact("parameters", "05_association/observed/parameters.json", ArtifactKind.PARAMETERS),
            _artifact("plot", "05_association/observed/module_trait_association.png", ArtifactKind.PLOT, False),
            _artifact("factorial_plot", "05_association/observed/factorial_module_association.png", ArtifactKind.PLOT, False),
        ),
    ),
    StageSpec(
        key="associate_residual",
        title="5B. Residual associations",
        order=51,
        description="Associate residual-network program scores with configured sample metadata.",
        dependencies=("network_residual",),
        cli_command="associate",
        artifacts=(
            _artifact("table", "05_association/residual/module_trait_association.tsv", ArtifactKind.TABLE),
            _artifact("factorial", "05_association/residual/factorial_module_association.tsv", ArtifactKind.TABLE),
            _artifact("coefficients", "05_association/residual/factorial_module_coefficients.tsv", ArtifactKind.TABLE),
            _artifact("group_summary", "05_association/residual/module_program_score_group_summary.tsv", ArtifactKind.TABLE),
            _artifact("design", "05_association/residual/association_design.tsv", ArtifactKind.TABLE),
            _artifact("parameters", "05_association/residual/parameters.json", ArtifactKind.PARAMETERS),
            _artifact("plot", "05_association/residual/module_trait_association.png", ArtifactKind.PLOT, False),
            _artifact("factorial_plot", "05_association/residual/factorial_module_association.png", ArtifactKind.PLOT, False),
        ),
    ),
    StageSpec(
        key="validate_modules",
        title="6. Validate modules",
        order=60,
        description="Measure bootstrap recovery, parameter sensitivity, and observed/residual preservation.",
        dependencies=("network_observed",),
        optional_dependencies=("network_residual",),
        cli_command="validate-network",
        artifacts=(
            _artifact("events", "06_module_validation/event_membership_stability.tsv", ArtifactKind.TABLE),
            _artifact("modules", "06_module_validation/module_bootstrap_stability.tsv", ArtifactKind.TABLE),
            _artifact("bootstrap", "06_module_validation/bootstrap_partition_stability.tsv", ArtifactKind.TABLE),
            _artifact("sensitivity", "06_module_validation/network_parameter_sensitivity.tsv", ArtifactKind.TABLE),
            _artifact("preservation_pairs", "06_module_validation/observed_residual_preservation_pairs.tsv", ArtifactKind.TABLE),
            _artifact("preservation_summary", "06_module_validation/observed_residual_preservation_summary.tsv", ArtifactKind.TABLE),
            _artifact("program_score_confidence", "06_module_validation/module_program_score_confidence.tsv", ArtifactKind.TABLE),
            _artifact("report", "06_module_validation/module_validation_report.html", ArtifactKind.REPORT),
            _artifact("parameters", "06_module_validation/parameters.json", ArtifactKind.PARAMETERS),
            _artifact("module_stability_plot", "06_module_validation/plots/module_bootstrap_stability.png", ArtifactKind.PLOT, False),
            _artifact("ari_plot", "06_module_validation/plots/bootstrap_adjusted_rand_index.png", ArtifactKind.PLOT, False),
            _artifact("sensitivity_plot", "06_module_validation/plots/parameter_sensitivity_heatmap.png", ArtifactKind.PLOT, False),
            _artifact("preservation_plot", "06_module_validation/plots/observed_residual_module_preservation.png", ArtifactKind.PLOT, False),
        ),
    ),
    StageSpec(
        key="annotate_modules",
        title="6B. Annotate module evidence",
        order=65,
        description=(
            "Add cached DEG, rMATS, DEXSeq, DTU and RBP gene-set support to observed modules; "
            "annotations never redefine network membership."
        ),
        dependencies=("network_observed",),
        cli_command="annotate-modules",
        artifacts=(
            _artifact("events", "06b_module_annotation/event_module_annotation.tsv", ArtifactKind.TABLE),
            _artifact("genes", "06b_module_annotation/module_gene_annotation.tsv", ArtifactKind.TABLE),
            _artifact("summary", "06b_module_annotation/module_evidence_summary.tsv", ArtifactKind.TABLE),
            _artifact("rbp", "06b_module_annotation/module_rbp_enrichment.tsv", ArtifactKind.TABLE),
            _artifact("sources", "06b_module_annotation/support_source_audit.tsv", ArtifactKind.TABLE),
            _artifact("parameters", "06b_module_annotation/parameters.json", ArtifactKind.PARAMETERS),
            _artifact("plot", "06b_module_annotation/plots/module_evidence_support_heatmap.png", ArtifactKind.PLOT, False),
        ),
    ),
    StageSpec(
        key="benchmark_methods",
        title="6C. Benchmark methods",
        order=67,
        description=(
            "Compare CoSplice with neutral hierarchical, graphical-lasso, and optional exact WGCNA "
            "baselines on one predeclared stable-AS event universe."
        ),
        dependencies=("stable_qc",),
        cli_command="benchmark-methods",
        artifacts=(
            _artifact("universe", "06c_method_benchmark/common_event_universe.tsv", ArtifactKind.TABLE),
            _artifact("assignments", "06c_method_benchmark/event_method_assignments.tsv", ArtifactKind.TABLE),
            _artifact("summary", "06c_method_benchmark/method_benchmark_summary.tsv", ArtifactKind.TABLE),
            _artifact("bootstrap", "06c_method_benchmark/bootstrap_partition_stability.tsv", ArtifactKind.TABLE),
            _artifact("matching", "06c_method_benchmark/module_matching.tsv", ArtifactKind.TABLE),
            _artifact("audit", "06c_method_benchmark/method_availability.tsv", ArtifactKind.TABLE),
            _artifact("metric_ranking", "06c_method_benchmark/method_metric_ranking.tsv", ArtifactKind.TABLE),
            _artifact("interpretation", "06c_method_benchmark/method_strengths_limitations.tsv", ArtifactKind.TABLE),
            _artifact("report", "06c_method_benchmark/method_benchmark_report.html", ArtifactKind.REPORT),
            _artifact("parameters", "06c_method_benchmark/parameters.json", ArtifactKind.PARAMETERS),
            _artifact("quality_plot", "06c_method_benchmark/plots/method_quality_heatmap.png", ArtifactKind.PLOT, False),
            _artifact("structure_plot", "06c_method_benchmark/plots/method_module_structure.png", ArtifactKind.PLOT, False),
            _artifact("agreement_plot", "06c_method_benchmark/plots/cross_method_module_agreement.png", ArtifactKind.PLOT, False),
            _artifact("scorecard_plot", "06c_method_benchmark/plots/head_to_head_scorecard.png", ArtifactKind.PLOT, False),
        ),
    ),
    StageSpec(
        key="preserve_modules",
        title="6D. External module preservation",
        order=68,
        description=(
            "Project bootstrap-consensus modules into configured independent cohorts and quantify "
            "density, connectivity, membership and program preservation with permutation calibration."
        ),
        dependencies=("consensus_modules",),
        cli_command="preserve-modules",
        artifacts=(
            _artifact("modules", "06d_external_preservation/external_module_preservation.tsv", ArtifactKind.TABLE),
            _artifact("mapping", "06d_external_preservation/external_event_mapping.tsv", ArtifactKind.TABLE),
            _artifact("null", "06d_external_preservation/external_preservation_null_summary.tsv", ArtifactKind.TABLE),
            _artifact("audit", "06d_external_preservation/external_preservation_target_audit.tsv", ArtifactKind.TABLE),
            _artifact("parameters", "06d_external_preservation/parameters.json", ArtifactKind.PARAMETERS),
            _artifact("zsummary_plot", "06d_external_preservation/plots/external_module_preservation_zsummary.png", ArtifactKind.PLOT, False),
            _artifact("components_plot", "06d_external_preservation/plots/external_module_preservation_components.png", ArtifactKind.PLOT, False),
        ),
    ),
    StageSpec(
        key="interpret_modules",
        title="6E. Module functional interpretation",
        order=69,
        description=(
            "Prioritize hub events, hub genes and candidate RBP regulators using cached module membership, "
            "DEG/DTU support, FASTA/PWM motif scans, RBP expression, GO/KEGG and isoform-switch evidence."
        ),
        dependencies=("network_observed",),
        optional_dependencies=("consensus_modules",),
        cli_command="interpret-modules",
        artifacts=(
            _artifact("hub_events", "06e_module_interpretation/module_hub_events.tsv", ArtifactKind.TABLE),
            _artifact("hub_genes", "06e_module_interpretation/module_hub_genes.tsv", ArtifactKind.TABLE),
            _artifact("joint", "06e_module_interpretation/module_joint_evidence.tsv", ArtifactKind.TABLE),
            _artifact("motif", "06e_module_interpretation/module_rbp_motif_enrichment.tsv", ArtifactKind.TABLE),
            _artifact("motif_hits", "06e_module_interpretation/event_rbp_motif_hits.tsv", ArtifactKind.TABLE),
            _artifact("rbp_expression", "06e_module_interpretation/rbp_module_expression_correlation.tsv", ArtifactKind.TABLE),
            _artifact("rbp_expression_matrix", "06e_module_interpretation/rbp_gene_expression_matrix.tsv", ArtifactKind.TABLE),
            _artifact("gene_sets", "06e_module_interpretation/module_go_kegg_enrichment.tsv", ArtifactKind.TABLE),
            _artifact("isoform_switches", "06e_module_interpretation/dominant_isoform_switches.tsv", ArtifactKind.TABLE),
            _artifact("isoform_enrichment", "06e_module_interpretation/module_isoform_switch_enrichment.tsv", ArtifactKind.TABLE),
            _artifact("regulators", "06e_module_interpretation/module_regulator_prioritization.tsv", ArtifactKind.TABLE),
            _artifact("audit", "06e_module_interpretation/interpretation_source_audit.tsv", ArtifactKind.TABLE),
            _artifact("motif_parameters", "06e_module_interpretation/motif_scan_parameters.json", ArtifactKind.PARAMETERS),
            _artifact("isoform_parameters", "06e_module_interpretation/isoform_switch_parameters.json", ArtifactKind.PARAMETERS),
            _artifact("expression_parameters", "06e_module_interpretation/rbp_expression_parameters.json", ArtifactKind.PARAMETERS),
            _artifact("parameters", "06e_module_interpretation/parameters.json", ArtifactKind.PARAMETERS),
        ),
    ),
    StageSpec(
        key="fit_module_models",
        title="6F. Replicate-aware module statistics",
        order=70,
        description=(
            "Fit beta-binomial condition-effect models to sample IJC/SJC and robust program-score models; "
            "batch is included only when explicitly configured."
        ),
        dependencies=("matrix", "network_observed"),
        optional_dependencies=("consensus_modules",),
        cli_command="fit-module-models",
        artifacts=(
            _artifact("events", "06f_module_statistics/event_beta_binomial_condition_effect.tsv", ArtifactKind.TABLE),
            _artifact("modules", "06f_module_statistics/module_program_condition_effect.tsv", ArtifactKind.TABLE),
            _artifact("variance", "06f_module_statistics/module_variance_partition.tsv", ArtifactKind.TABLE),
            _artifact("benchmark", "06f_module_statistics/beta_binomial_vs_simple_delta_dpsi.tsv", ArtifactKind.TABLE),
            _artifact("diagnostics", "06f_module_statistics/beta_binomial_model_diagnostics.tsv", ArtifactKind.TABLE),
            _artifact("parameters", "06f_module_statistics/parameters.json", ArtifactKind.PARAMETERS),
        ),
    ),
    StageSpec(
        key="meta_analyze_modules",
        title="6G. Cross-cohort module meta-analysis",
        order=71,
        description=(
            "Harmonize exact canonical AS events across independent cohorts, match independently discovered modules "
            "when supplied, and meta-analyze module program-score effects."
        ),
        dependencies=("consensus_modules",),
        cli_command="meta-analyze-modules",
        artifacts=(
            _artifact("alignment", "06g_cross_cohort_meta/exact_event_alignment.tsv", ArtifactKind.TABLE),
            _artifact("matching", "06g_cross_cohort_meta/cross_cohort_module_matching.tsv", ArtifactKind.TABLE),
            _artifact("states", "06g_cross_cohort_meta/cross_cohort_module_state.tsv", ArtifactKind.TABLE),
            _artifact("effects", "06g_cross_cohort_meta/cohort_module_effects.tsv", ArtifactKind.TABLE),
            _artifact("meta", "06g_cross_cohort_meta/random_effects_module_meta_analysis.tsv", ArtifactKind.TABLE),
            _artifact("audit", "06g_cross_cohort_meta/meta_analysis_cohort_audit.tsv", ArtifactKind.TABLE),
            _artifact("consensus", "06g_cross_cohort_meta/multi_cohort_module_consensus.tsv", ArtifactKind.TABLE),
            _artifact("parameters", "06g_cross_cohort_meta/parameters.json", ArtifactKind.PARAMETERS),
        ),
    ),
    StageSpec(
        key="factor_contrasts",
        title="6H. Configured factor contrasts and phenocopy",
        order=72,
        description=(
            "Test predeclared target-reference, perturbation and difference-in-differences contrasts on "
            "stable sample-level PSI and cached CoSplice program scores. Phenocopy results prioritize "
            "mechanistic follow-up and do not by themselves prove loss or gain of function."
        ),
        dependencies=("stable_qc",),
        optional_dependencies=("network_observed", "network_residual"),
        cli_command="contrast-factors",
        artifacts=(
            _artifact("events", "06h_factor_contrasts/event_factor_contrasts.tsv", ArtifactKind.TABLE),
            _artifact("modules", "06h_factor_contrasts/module_factor_contrasts.tsv", ArtifactKind.TABLE),
            _artifact("phenocopy", "06h_factor_contrasts/phenocopy_occlusion_summary.tsv", ArtifactKind.TABLE),
            _artifact("evidence", "06h_factor_contrasts/factor_evidence_summary.tsv", ArtifactKind.TABLE),
            _artifact("notes", "06h_factor_contrasts/presentation_notes.md", ArtifactKind.REPORT),
            _artifact("parameters", "06h_factor_contrasts/parameters.json", ArtifactKind.PARAMETERS),
            _artifact("dashboard", "06h_factor_contrasts/plots/factor_phenocopy_dashboard.png", ArtifactKind.PLOT, False),
            _artifact("scatter", "06h_factor_contrasts/plots/factor_phenocopy_signature_scatter.png", ArtifactKind.PLOT, False),
            _artifact("modules_plot", "06h_factor_contrasts/plots/module_perturbation_occlusion_heatmap.png", ArtifactKind.PLOT, False),
        ),
    ),
    StageSpec(
        key="differential_network",
        title="6I. Differential co-splicing network",
        order=73,
        description=(
            "Compare module density and edge topology across condition-specific "
            "networks while retaining one fixed consensus event universe."
        ),
        dependencies=("network_observed", "consensus_modules"),
        cli_command="differential-network",
        artifacts=(
            _artifact(
                "modules",
                "06i_differential_network/module_differential_network.tsv",
                ArtifactKind.TABLE,
            ),
            _artifact(
                "events",
                "06i_differential_network/event_network_rewiring.tsv",
                ArtifactKind.TABLE,
            ),
            _artifact(
                "edges",
                "06i_differential_network/condition_network_edges.tsv",
                ArtifactKind.TABLE,
            ),
            _artifact(
                "parameters",
                "06i_differential_network/parameters.json",
                ArtifactKind.PARAMETERS,
            ),
            _artifact("density_plot", "06i_differential_network/plots/differential_module_density_forest.png", ArtifactKind.PLOT, False),
            _artifact("rewiring_plot", "06i_differential_network/plots/differential_module_rewiring_heatmap.png", ArtifactKind.PLOT, False),
        ),
    ),
    StageSpec(
        key="multiview_fusion",
        title="6J. AS-primary multi-view graph fusion",
        order=74,
        description=(
            "Reweight existing AS edges with configured RBP, DEG, DTU or isoform "
            "features without allowing auxiliary evidence to create new AS edges."
        ),
        dependencies=("network_observed",),
        cli_command="multiview-fusion",
        artifacts=(
            _artifact(
                "edges",
                "06j_multiview_fusion/fused_network_edges.tsv",
                ArtifactKind.TABLE,
            ),
            _artifact(
                "assignments",
                "06j_multiview_fusion/multiview_module_assignments.tsv",
                ArtifactKind.TABLE,
            ),
            _artifact(
                "audit",
                "06j_multiview_fusion/view_audit.tsv",
                ArtifactKind.TABLE,
            ),
            _artifact(
                "parameters",
                "06j_multiview_fusion/parameters.json",
                ArtifactKind.PARAMETERS,
            ),
        ),
    ),
    StageSpec(
        key="bayesian_preservation",
        title="6K. Bayesian cross-cohort preservation",
        order=75,
        description=(
            "Combine cached density, connectivity, membership and program-score "
            "preservation across independent cohorts into posterior probabilities."
        ),
        dependencies=("preserve_modules",),
        cli_command="bayesian-preservation",
        artifacts=(
            _artifact(
                "modules",
                "06k_bayesian_preservation/module_preservation_posterior.tsv",
                ArtifactKind.TABLE,
            ),
            _artifact(
                "metrics",
                "06k_bayesian_preservation/metric_preservation_posterior.tsv",
                ArtifactKind.TABLE,
            ),
            _artifact(
                "parameters",
                "06k_bayesian_preservation/parameters.json",
                ArtifactKind.PARAMETERS,
            ),
        ),
    ),
    StageSpec(
        key="report",
        title="7. Build report",
        order=80,
        description="Summarize existing cached outputs and provenance without recomputing analyses.",
        optional_dependencies=(
            "audit",
            "matrix",
            "stable_qc",
            "network_observed",
            "network_residual",
            "consensus_modules",
            "associate_observed",
            "associate_residual",
            "validate_modules",
            "annotate_modules",
            "benchmark_methods",
            "preserve_modules",
            "interpret_modules",
            "fit_module_models",
            "meta_analyze_modules",
            "factor_contrasts",
            "differential_network",
            "multiview_fusion",
            "bayesian_preservation",
        ),
        cli_command="report",
        accepts_force=False,
        artifacts=(
            _artifact("report", "07_report/katana_cosplice_report.html", ArtifactKind.REPORT),
            _artifact("pdf", "07_report/katana_cosplice_validation_report.pdf", ArtifactKind.REPORT),
            _artifact("browser", "07_report/module_browser.html", ArtifactKind.REPORT),
        ),
    ),
)

STAGE_BY_KEY = {spec.key: spec for spec in STAGE_SPECS}


def stage_specs(*, visible_only: bool = False) -> tuple[StageSpec, ...]:
    specs = STAGE_SPECS
    if visible_only:
        specs = tuple(spec for spec in specs if spec.visible and spec.implemented)
    return tuple(sorted(specs, key=lambda spec: spec.order))


def stage_spec(key: str) -> StageSpec:
    try:
        return STAGE_BY_KEY[key]
    except KeyError as exc:
        raise KeyError(f"Unknown CoSplice stage: {key}") from exc


def validate_stage_registry() -> None:
    keys = [spec.key for spec in STAGE_SPECS]
    if len(keys) != len(set(keys)):
        raise ValueError("CoSplice stage keys must be unique.")
    if len({spec.order for spec in STAGE_SPECS}) != len(STAGE_SPECS):
        raise ValueError("CoSplice stage order values must be unique.")
    known = set(keys)
    for spec in STAGE_SPECS:
        missing = set(spec.dependencies + spec.optional_dependencies) - known
        if missing:
            raise ValueError(f"Stage {spec.key} references unknown dependencies: {sorted(missing)}")


validate_stage_registry()
