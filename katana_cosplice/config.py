from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import yaml


STABILITY_PRESETS: dict[str, dict[str, float]] = {
    "lenient": {
        "min_valid_fraction": 0.70,
        "min_read_supported_fraction": 0.50,
        "min_total_junction_reads": 8.0,
        "min_psi_variance": 0.0010,
        "min_repeatability": 0.15,
        "min_condition_supported_fraction": 0.70,
        "min_supported_replicates_per_condition": 2.0,
        "min_condition_mean_psi_range": 0.20,
        "min_alternative_fraction": 0.05,
        "alternative_psi_low": 0.20,
        "alternative_psi_high": 0.80,
        "min_sample_psi_range": 0.40,
    },
    "default": {
        "min_valid_fraction": 0.80,
        "min_read_supported_fraction": 0.70,
        "min_total_junction_reads": 10.0,
        "min_psi_variance": 0.0025,
        "min_repeatability": 0.30,
        "min_condition_supported_fraction": 0.80,
        "min_supported_replicates_per_condition": 2.0,
        "min_condition_mean_psi_range": 0.25,
        "min_alternative_fraction": 0.10,
        "alternative_psi_low": 0.20,
        "alternative_psi_high": 0.80,
        "min_sample_psi_range": 0.50,
    },
    "strict": {
        "min_valid_fraction": 0.90,
        "min_read_supported_fraction": 0.80,
        "min_total_junction_reads": 20.0,
        "min_psi_variance": 0.0050,
        "min_repeatability": 0.50,
        "min_condition_supported_fraction": 0.90,
        "min_supported_replicates_per_condition": 3.0,
        "min_condition_mean_psi_range": 0.30,
        "min_alternative_fraction": 0.15,
        "alternative_psi_low": 0.20,
        "alternative_psi_high": 0.80,
        "min_sample_psi_range": 0.60,
    },
}


@dataclass
class StableASConfig:
    preset: str = "default"
    min_valid_fraction: float | None = None
    min_read_supported_fraction: float | None = None
    min_total_junction_reads: float | None = None
    min_psi_variance: float | None = None
    min_repeatability: float | None = None
    min_condition_supported_fraction: float | None = None
    min_supported_replicates_per_condition: float | None = None
    min_condition_mean_psi_range: float | None = None
    min_alternative_fraction: float | None = None
    alternative_psi_low: float | None = None
    alternative_psi_high: float | None = None
    min_sample_psi_range: float | None = None

    def resolved(self) -> dict[str, float]:
        if self.preset not in STABILITY_PRESETS:
            raise ValueError(f"Unknown stable-AS preset: {self.preset}")
        values = dict(STABILITY_PRESETS[self.preset])
        for key in values:
            override = getattr(self, key)
            if override is not None:
                values[key] = float(override)
        return values


@dataclass
class NetworkConfig:
    minimum_network_samples: int = 15
    engine: str = "tom_hierarchical"
    module_detection: str = "wgcna_blockwise"
    method: str = "bicor"
    signed: bool = True
    soft_power: int = 6
    edge_threshold: float = 0.10
    module_cut_height: float = 0.80
    min_module_size: int = 20
    max_events: int | None = None
    tom_dense_max_events: int = 3_000
    tom_cut_heights: list[float] = field(
        default_factory=lambda: [0.70, 0.75, 0.80, 0.85, 0.90]
    )
    tom_cut_selection_bootstraps: int = 10
    dynamic_cut_height: float = 0.995
    dynamic_deep_split: int = 2
    dynamic_deep_splits: list[int] = field(default_factory=lambda: [1, 2, 3, 4])
    dynamic_pam_respects_dendro: bool = False
    wgcna_rscript: str | None = None
    wgcna_merge_cut_height: float = 0.25
    knn_neighbors: int = 50
    edge_smoothing: str = "shared_neighbor"
    shared_neighbor_weight: float = 0.50
    shared_neighbor_top_k: int = 50
    multi_k_consensus: bool = True
    consensus_knn_neighbors: list[int] = field(default_factory=lambda: [25, 50, 75])
    consensus_edge_min_fraction: float = 0.66
    consensus_bootstrap_iterations: int = 50
    consensus_bootstrap_weight: float = 0.75
    consensus_bootstrap_min_support: float = 0.70
    community_consensus_runs: int = 20
    community_consensus_weight: float = 0.75
    community_resolution: float = 1.0
    module_refinement: bool = True
    min_module_membership: float = 0.45
    min_module_membership_margin: float = 0.0
    adaptive_membership: bool = True
    adaptive_membership_floor: float = 0.08
    adaptive_membership_mad_scale: float = 2.0
    module_merge_correlation: float = 0.85
    refinement_iterations: int = 2
    max_output_edges: int = 100_000
    bootstrap_iterations: int = 50
    bootstrap_edge_threshold: float = 0.70
    random_seed: int = 17


@dataclass
class UncertaintyConfig:
    enabled: bool = True
    network_mode: str = "posterior"
    prior_strength: float = 1.0
    reliability_read_scale: float = 20.0
    reliability_floor: float = 0.10
    reliability_power: float = 0.50
    filter_low_reliability: bool = True
    min_network_reliability: float = 0.15
    posterior_draws: int = 25
    posterior_random_seed: int = 17
    min_edge_sign_probability: float = 0.80
    posterior_edge_diagnostics: bool = True


@dataclass
class GenePriorConfig:
    enabled: bool = False
    exon_count_matrix: Path | None = None
    feature_gene_map: Path | None = None
    feature_id_column: str = "feature_id"
    gene_id_column: str = "gene_id"
    encoded_gene_rule: str = "prefix_before_last_underscore"
    min_exons: int = 2
    max_genes: int | None = None
    soft_power: float = 6.0
    weight: float = 0.10
    prior_mode: str = "scaffold_reinforcement"
    scaffold_knn_neighbors: int = 25
    scaffold_min_module_size: int = 10
    scaffold_resolution: float = 1.0
    random_seed: int = 73


@dataclass
class DifferentialNetworkConfig:
    enabled: bool = True
    condition_column: str = "condition"
    min_samples_per_condition: int = 3
    shrinkage: float = 5.0
    permutations: int = 200
    bootstrap_iterations: int = 200
    max_output_edges_per_module_pair: int = 5_000
    random_seed: int = 71


@dataclass
class MultiViewConfig:
    enabled: bool = False
    event_feature_table: Path | None = None
    auxiliary_columns: list[str] = field(default_factory=list)
    as_weight: float = 0.80
    auxiliary_weights: dict[str, float] = field(default_factory=dict)
    edge_threshold: float = 0.10
    min_module_size: int = 20
    resolution: float = 1.0


@dataclass
class BayesianPreservationConfig:
    enabled: bool = True
    prior_sd: float = 5.0
    tau_scale: float = 2.0
    grid_size: int = 301


@dataclass
class ValidationConfig:
    bootstrap_iterations: int = 20
    robust_jaccard: float = 0.75
    moderate_jaccard: float = 0.50
    knn_neighbors: list[int] = field(default_factory=lambda: [25, 50, 75])
    community_resolutions: list[float] = field(default_factory=lambda: [0.75, 1.0, 1.25])
    random_seed: int = 29


@dataclass
class ConsensusConfig:
    bootstrap_iterations: int = 20
    min_coassignment: float = 0.70
    min_module_size: int | None = None
    clustering_method: str = "backbone_confidence_only"
    dynamic_cut_height: float = 0.995
    dynamic_deep_splits: list[int] = field(default_factory=lambda: [1, 2, 3, 4])
    dynamic_pam_respects_dendro: bool = False
    hierarchical_thresholds: list[float] = field(
        default_factory=lambda: [0.50, 0.60, 0.70, 0.80]
    )
    graph_edge_thresholds: list[float] = field(
        default_factory=lambda: [0.50, 0.60, 0.70, 0.80]
    )
    graph_resolutions: list[float] = field(
        default_factory=lambda: [0.50, 0.80, 1.10, 1.40]
    )
    graph_top_k: int = 50
    min_candidate_modules: int = 2
    max_candidate_modules: int = 30
    core_membership_min: float = 0.70
    peripheral_membership_min: float = 0.45
    peripheral_membership_margin: float = 0.05
    core_kme_mean_min: float = 0.70
    core_kme_lower_quantile_min: float = 0.50
    core_kme_sign_consistency_min: float = 0.80
    peripheral_kme_mean_min: float = 0.50
    peripheral_kme_lower_quantile_min: float = 0.30
    peripheral_kme_sign_consistency_min: float = 0.70
    max_grey_fraction: float = 0.25
    max_largest_module_fraction: float = 0.60
    dense_max_events: int = 8_000
    output_edge_threshold: float = 0.50
    max_output_edges: int = 200_000
    random_seed: int = 37
    apply_threshold_labels: bool = False


@dataclass
class MultiObjectiveSelectionConfig:
    """Leakage-safe validation objectives; legacy modes may use them for selection."""

    enabled: bool = False
    tuning_cohorts: list[str] = field(default_factory=list)
    holdout_cohorts: list[str] = field(default_factory=list)
    functional_gene_sets_gmt: Path | None = None
    stability_weight: float = 0.40
    preservation_weight: float = 0.35
    functional_weight: float = 0.25
    functional_weight_mode: str = "auto_significance"
    min_significant_functional_modules: int = 1
    min_shared_events_per_module: int = 10
    min_module_genes: int = 3
    functional_fdr_cutoff: float = 0.10
    score_floor: float = 0.01


@dataclass
class ExternalCohortConfig:
    name: str
    psi_matrix: Path
    sample_metadata: Path
    module_assignments: Path | None = None
    effect_table: Path | None = None
    event_metadata: Path | None = None
    event_mapping_table: Path | None = None
    preservation_table: Path | None = None
    genome_build: str = ""
    annotation_version: str = ""


@dataclass
class ModuleInterpretationConfig:
    rbp_expression_matrix: Path | None = None
    rbp_motif_hits: Path | None = None
    rbp_motif_database: Path | None = None
    rbp_motif_organism: str = ""
    gene_sets_gmt: Path | None = None
    motif_window_bp: int = 250
    motif_score_threshold: float = 0.80
    gprofiler_enabled: bool = False
    gprofiler_organism: str = "mmusculus"
    isoform_min_gene_tpm: float = 1.0
    isoform_min_abs_delta_if: float = 0.10
    min_module_genes: int = 3
    fdr_cutoff: float = 0.05


@dataclass
class ModuleStatisticsConfig:
    condition_column: str = "condition"
    min_total_reads: int = 10
    min_samples_per_condition: int = 2
    fdr_cutoff: float = 0.05
    max_events: int | None = None
    optimizer_max_iterations: int = 200
    random_effect_column: str | None = None
    quadrature_points: int = 9


@dataclass
class MetaAnalysisConfig:
    min_shared_events: int = 10
    min_cohorts: int = 2
    fdr_cutoff: float = 0.05
    event_matching: str = "annotation_aware"
    min_supported_cohort_fraction: float = 0.50


@dataclass
class FactorContrastConfig:
    """Predeclared factor contrasts for phenotype and perturbation follow-up."""

    enabled: bool = False
    primary_factor: str = ""
    target_level: str = ""
    reference_level: str = ""
    context_factors: list[str] = field(default_factory=list)
    perturbation_factor: str = ""
    baseline_level: str = ""
    perturbation_levels: list[str] = field(default_factory=list)
    network_source: str = "residual"
    fdr_cutoff: float = 0.05
    min_abs_delta_psi: float = 0.10
    min_samples_per_group: int = 2


@dataclass
class PreservationConfig:
    permutations: int = 100
    min_shared_events: int = 20
    dense_max_events: int = 8_000
    preserved_zsummary: float = 10.0
    partial_zsummary: float = 2.0
    random_seed: int = 43


@dataclass
class HybridProgramConfig:
    """Settings for retaining both tight sparse and diffuse TOM programs."""

    matrix_source: str = "auto"
    diffuse_engine: str = "auto"
    bootstrap_iterations: int = 10
    min_module_size: int | None = None
    python_tom_cut_height: float = 0.95
    python_max_events: int = 8_000
    robust_confidence: float = 0.75
    moderate_confidence: float = 0.55
    random_seed: int = 53


@dataclass
class BenchmarkConfig:
    methods: list[str] = field(
        default_factory=lambda: [
            "katana_sparse",
            "pearson_hierarchical",
            "graphical_lasso",
            "wgcna_exact",
        ]
    )
    max_events: int = 150
    bootstrap_iterations: int = 5
    graphical_lasso_alpha: float = 0.10
    graphical_lasso_max_iter: int = 100
    neutral_knn_neighbors: int = 15
    wgcna_merge_cut_height: float = 0.25
    random_seed: int = 41


@dataclass
class MatrixInputConfig:
    """Controls the sample-level splice-usage matrix source."""

    input_mode: str = "rmats_events"
    backend: str = "auto"
    min_mapping_quality: int = 1
    min_anchor_bp: int = 8
    max_intron_length: int = 500_000
    min_total_junction_reads: int = 20
    min_samples_with_reads: int = 2
    min_cluster_fraction: float = 0.001
    exclude_duplicates: bool = True
    max_parallel_samples: int = 1


@dataclass
class CoSpliceConfig:
    project_name: str
    rmats_root: Path | None
    sample_metadata: Path
    output_dir: Path
    quant_root: Path | None = None
    gtf: Path | None = None
    genome_fasta: Path | None = None
    transcript_fasta: Path | None = None
    bam_root: Path | None = None
    support_root: Path | None = None
    rbp_gene_set: Path | None = None
    genome_build: str = ""
    annotation_version: str = ""
    gene_exclusion: dict[str, Any] = field(default_factory=dict)
    rmats_mode: str = "JCEC"
    selected_conditions: list[str] = field(default_factory=list)
    residual_covariates: list[str] = field(default_factory=list)
    association_factors: list[str] = field(default_factory=lambda: ["condition"])
    factorial_factors: list[str] = field(default_factory=list)
    batch_column: str | None = None
    matrix_input: MatrixInputConfig = field(default_factory=MatrixInputConfig)
    stable_as: StableASConfig = field(default_factory=StableASConfig)
    network: NetworkConfig = field(default_factory=NetworkConfig)
    uncertainty: UncertaintyConfig = field(default_factory=UncertaintyConfig)
    gene_prior: GenePriorConfig = field(default_factory=GenePriorConfig)
    differential_network: DifferentialNetworkConfig = field(
        default_factory=DifferentialNetworkConfig
    )
    multiview: MultiViewConfig = field(default_factory=MultiViewConfig)
    bayesian_preservation: BayesianPreservationConfig = field(
        default_factory=BayesianPreservationConfig
    )
    validation: ValidationConfig = field(default_factory=ValidationConfig)
    consensus: ConsensusConfig = field(default_factory=ConsensusConfig)
    multiobjective_selection: MultiObjectiveSelectionConfig = field(
        default_factory=MultiObjectiveSelectionConfig
    )
    external_cohorts: list[ExternalCohortConfig] = field(default_factory=list)
    preservation: PreservationConfig = field(default_factory=PreservationConfig)
    hybrid_programs: HybridProgramConfig = field(default_factory=HybridProgramConfig)
    benchmark: BenchmarkConfig = field(default_factory=BenchmarkConfig)
    interpretation: ModuleInterpretationConfig = field(default_factory=ModuleInterpretationConfig)
    statistics: ModuleStatisticsConfig = field(default_factory=ModuleStatisticsConfig)
    meta_analysis: MetaAnalysisConfig = field(default_factory=MetaAnalysisConfig)
    factor_contrasts: FactorContrastConfig = field(default_factory=FactorContrastConfig)

    def validate(self) -> None:
        minimum_samples = self.network.minimum_network_samples
        if isinstance(minimum_samples, bool) or not isinstance(minimum_samples, int) or minimum_samples < 3:
            raise ValueError("network.minimum_network_samples must be an integer of at least 3")
        if self.rmats_mode not in {"JC", "JCEC"}:
            raise ValueError("rmats_mode must be JC or JCEC")
        if not self.project_name.strip():
            raise ValueError("project_name is required")
        if self.matrix_input.input_mode not in {"rmats_events", "bam_junctions"}:
            raise ValueError("matrix_input.input_mode must be rmats_events or bam_junctions")
        if self.matrix_input.backend not in {
            "auto",
            "pysam",
            "wsl_pysam",
            "bamnostic",
            "precomputed_junctions",
        }:
            raise ValueError(
                "matrix_input.backend must be auto, pysam, wsl_pysam, bamnostic "
                "or precomputed_junctions"
            )
        if self.matrix_input.min_mapping_quality < 0:
            raise ValueError("matrix_input.min_mapping_quality cannot be negative")
        if self.matrix_input.min_anchor_bp < 1:
            raise ValueError("matrix_input.min_anchor_bp must be positive")
        if self.matrix_input.max_intron_length < 1:
            raise ValueError("matrix_input.max_intron_length must be positive")
        if self.matrix_input.min_total_junction_reads < 1:
            raise ValueError("matrix_input.min_total_junction_reads must be positive")
        if self.matrix_input.min_samples_with_reads < 1:
            raise ValueError("matrix_input.min_samples_with_reads must be positive")
        if not 0.0 <= self.matrix_input.min_cluster_fraction < 1.0:
            raise ValueError("matrix_input.min_cluster_fraction must be within [0, 1)")
        if self.matrix_input.max_parallel_samples < 1:
            raise ValueError("matrix_input.max_parallel_samples must be positive")
        for name, values in (
            ("selected_conditions", self.selected_conditions),
            ("residual_covariates", self.residual_covariates),
            ("association_factors", self.association_factors),
            ("factorial_factors", self.factorial_factors),
        ):
            if any(not str(value).strip() for value in values):
                raise ValueError(f"{name} cannot contain blank values")
            if len(values) != len(set(values)):
                raise ValueError(f"{name} cannot contain duplicate values")
        stable = self.stable_as.resolved()
        for key in (
            "min_valid_fraction",
            "min_read_supported_fraction",
            "min_repeatability",
            "min_condition_supported_fraction",
            "min_alternative_fraction",
            "alternative_psi_low",
            "alternative_psi_high",
        ):
            if not 0.0 <= stable[key] <= 1.0:
                raise ValueError(f"stable_as.{key} must be within [0, 1]")
        if stable["alternative_psi_low"] >= stable["alternative_psi_high"]:
            raise ValueError("stable_as alternative_psi_low must be lower than alternative_psi_high")
        if stable["min_total_junction_reads"] < 0:
            raise ValueError("stable_as.min_total_junction_reads cannot be negative")
        if stable["min_supported_replicates_per_condition"] < 1:
            raise ValueError("stable_as.min_supported_replicates_per_condition must be at least 1")
        if stable["min_psi_variance"] < 0:
            raise ValueError("stable_as.min_psi_variance cannot be negative")
        if not 0.0 <= stable["min_condition_mean_psi_range"] <= 1.0:
            raise ValueError("stable_as.min_condition_mean_psi_range must be within [0, 1]")
        if not 0.0 <= stable["min_sample_psi_range"] <= 1.0:
            raise ValueError("stable_as.min_sample_psi_range must be within [0, 1]")
        if self.network.method not in {"bicor", "spearman", "pearson"}:
            raise ValueError("network.method must be bicor, spearman or pearson")
        if self.network.engine not in {
            "tom_hierarchical",
            "sparse_knn",
            "exact_hierarchical",
        }:
            raise ValueError(
                "network.engine must be tom_hierarchical, sparse_knn or "
                "exact_hierarchical"
            )
        if self.network.soft_power < 1:
            raise ValueError("network.soft_power must be positive")
        if self.network.min_module_size < 2:
            raise ValueError("network.min_module_size must be at least 2")
        if self.network.engine == "exact_hierarchical" and self.network.max_events is None:
            raise ValueError("exact_hierarchical requires an explicit network.max_events memory cap")
        if self.network.max_events is not None and self.network.max_events < self.network.min_module_size * 2:
            raise ValueError("network.max_events must be at least twice min_module_size")
        if self.network.tom_dense_max_events < self.network.min_module_size * 2:
            raise ValueError(
                "network.tom_dense_max_events must be at least twice min_module_size"
            )
        if not self.network.tom_cut_heights or any(
            not 0.0 < value < 1.0 for value in self.network.tom_cut_heights
        ):
            raise ValueError(
                "network.tom_cut_heights must contain values within (0, 1)"
            )
        if self.network.tom_cut_selection_bootstraps < 2:
            raise ValueError(
                "network.tom_cut_selection_bootstraps must be at least 2"
            )
        if self.network.module_detection not in {
            "wgcna_blockwise",
            "dynamic_tree_cut",
            "static_cut",
        }:
            raise ValueError(
                "network.module_detection must be wgcna_blockwise, "
                "dynamic_tree_cut or static_cut"
            )
        if not 0.0 < self.network.wgcna_merge_cut_height < 1.0:
            raise ValueError(
                "network.wgcna_merge_cut_height must be within (0, 1)"
            )
        if not 0.0 < self.network.dynamic_cut_height <= 1.0:
            raise ValueError(
                "network.dynamic_cut_height must be within (0, 1]"
            )
        if self.network.dynamic_deep_split not in {0, 1, 2, 3, 4}:
            raise ValueError("network.dynamic_deep_split must be within 0..4")
        if not self.network.dynamic_deep_splits or any(
            value not in {0, 1, 2, 3, 4}
            for value in self.network.dynamic_deep_splits
        ):
            raise ValueError(
                "network.dynamic_deep_splits must contain values within 0..4"
            )
        if self.network.knn_neighbors < 2:
            raise ValueError("network.knn_neighbors must be at least 2")
        if self.network.edge_smoothing not in {"none", "shared_neighbor"}:
            raise ValueError("network.edge_smoothing must be none or shared_neighbor")
        if not 0.0 <= self.network.shared_neighbor_weight <= 1.0:
            raise ValueError("network.shared_neighbor_weight must be within [0, 1]")
        if self.network.shared_neighbor_top_k < 2:
            raise ValueError("network.shared_neighbor_top_k must be at least 2")
        if not self.network.consensus_knn_neighbors or any(
            value < 2 for value in self.network.consensus_knn_neighbors
        ):
            raise ValueError("network.consensus_knn_neighbors must contain values >= 2")
        if not 0.0 < self.network.consensus_edge_min_fraction <= 1.0:
            raise ValueError("network.consensus_edge_min_fraction must be within (0, 1]")
        if self.network.consensus_bootstrap_iterations < 0:
            raise ValueError("network.consensus_bootstrap_iterations must be non-negative")
        if not 0.0 <= self.network.consensus_bootstrap_weight <= 1.0:
            raise ValueError("network.consensus_bootstrap_weight must be within [0, 1]")
        if not 0.0 <= self.network.consensus_bootstrap_min_support <= 1.0:
            raise ValueError("network.consensus_bootstrap_min_support must be within [0, 1]")
        if self.network.community_consensus_runs < 1:
            raise ValueError("network.community_consensus_runs must be positive")
        if not 0.0 <= self.network.community_consensus_weight <= 1.0:
            raise ValueError("network.community_consensus_weight must be within [0, 1]")
        if self.network.community_resolution <= 0:
            raise ValueError("network.community_resolution must be positive")
        if not 0.0 <= self.network.min_module_membership <= 1.0:
            raise ValueError("network.min_module_membership must be within [0, 1]")
        if not 0.0 <= self.network.min_module_membership_margin <= 1.0:
            raise ValueError("network.min_module_membership_margin must be within [0, 1]")
        if not 0.0 <= self.network.adaptive_membership_floor <= self.network.min_module_membership:
            raise ValueError(
                "network.adaptive_membership_floor must be within [0, min_module_membership]"
            )
        if self.network.adaptive_membership_mad_scale < 0.0:
            raise ValueError("network.adaptive_membership_mad_scale must be non-negative")
        if not 0.0 < self.network.module_merge_correlation <= 1.0:
            raise ValueError("network.module_merge_correlation must be within (0, 1]")
        if self.network.refinement_iterations < 1:
            raise ValueError("network.refinement_iterations must be positive")
        if self.network.max_output_edges < 1:
            raise ValueError("network.max_output_edges must be positive")
        if not 0.0 < self.network.module_cut_height < 2.0:
            raise ValueError("network.module_cut_height must be between 0 and 2")
        if not 0.0 <= self.network.edge_threshold <= 1.0:
            raise ValueError("network.edge_threshold must be between 0 and 1")
        if not 0.0 <= self.network.bootstrap_edge_threshold <= 1.0:
            raise ValueError("network.bootstrap_edge_threshold must be between 0 and 1")
        if self.network.bootstrap_iterations < 0:
            raise ValueError("network.bootstrap_iterations cannot be negative")
        if self.uncertainty.prior_strength <= 0:
            raise ValueError("uncertainty.prior_strength must be positive")
        if self.uncertainty.network_mode not in {
            "diagnostic_only",
            "raw_reliability_weighted",
            "posterior",
            "posterior_reliability_weighted",
            "posterior_ensemble",
        }:
            raise ValueError(
                "uncertainty.network_mode must be diagnostic_only, "
                "raw_reliability_weighted, posterior, "
                "posterior_reliability_weighted or posterior_ensemble"
            )
        if self.uncertainty.reliability_read_scale <= 0:
            raise ValueError("uncertainty.reliability_read_scale must be positive")
        if not 0.0 <= self.uncertainty.reliability_floor <= 1.0:
            raise ValueError("uncertainty.reliability_floor must be within [0, 1]")
        if self.uncertainty.reliability_power <= 0:
            raise ValueError("uncertainty.reliability_power must be positive")
        if not 0.0 <= self.uncertainty.min_network_reliability <= 1.0:
            raise ValueError(
                "uncertainty.min_network_reliability must be within [0, 1]"
            )
        if self.uncertainty.posterior_draws < 2:
            raise ValueError("uncertainty.posterior_draws must be at least 2")
        if not 0.5 <= self.uncertainty.min_edge_sign_probability <= 1.0:
            raise ValueError(
                "uncertainty.min_edge_sign_probability must be within [0.5, 1]"
            )
        if self.gene_prior.enabled and self.gene_prior.exon_count_matrix is None:
            raise ValueError("Enabled gene_prior requires exon_count_matrix")
        if self.gene_prior.encoded_gene_rule not in {
            "prefix_before_last_underscore",
            "exact_feature_id",
        }:
            raise ValueError(
                "gene_prior.encoded_gene_rule must be "
                "prefix_before_last_underscore or exact_feature_id"
            )
        if self.gene_prior.min_exons < 2:
            raise ValueError("gene_prior.min_exons must be at least 2")
        if self.gene_prior.max_genes is not None and self.gene_prior.max_genes < 2:
            raise ValueError("gene_prior.max_genes must be at least 2 when configured")
        if self.gene_prior.soft_power <= 0:
            raise ValueError("gene_prior.soft_power must be positive")
        if not 0.0 <= self.gene_prior.weight <= 0.5:
            raise ValueError(
                "gene_prior.weight must be within [0, 0.5] so AS remains primary"
            )
        if self.gene_prior.prior_mode not in {
            "continuous_adjacency",
            "scaffold_reinforcement",
        }:
            raise ValueError(
                "gene_prior.prior_mode must be continuous_adjacency or "
                "scaffold_reinforcement"
            )
        if self.gene_prior.scaffold_knn_neighbors < 2:
            raise ValueError("gene_prior.scaffold_knn_neighbors must be at least 2")
        if self.gene_prior.scaffold_min_module_size < 2:
            raise ValueError("gene_prior.scaffold_min_module_size must be at least 2")
        if self.gene_prior.scaffold_resolution <= 0:
            raise ValueError("gene_prior.scaffold_resolution must be positive")
        if not self.gene_prior.feature_id_column.strip():
            raise ValueError("gene_prior.feature_id_column is required")
        if not self.gene_prior.gene_id_column.strip():
            raise ValueError("gene_prior.gene_id_column is required")
        if self.differential_network.min_samples_per_condition < 3:
            raise ValueError(
                "differential_network.min_samples_per_condition must be at least 3"
            )
        if self.differential_network.shrinkage < 0:
            raise ValueError("differential_network.shrinkage cannot be negative")
        if self.differential_network.permutations < 19:
            raise ValueError("differential_network.permutations must be at least 19")
        if self.differential_network.bootstrap_iterations < 20:
            raise ValueError(
                "differential_network.bootstrap_iterations must be at least 20"
            )
        if self.differential_network.max_output_edges_per_module_pair < 1:
            raise ValueError(
                "differential_network.max_output_edges_per_module_pair must be positive"
            )
        if not self.differential_network.condition_column.strip():
            raise ValueError("differential_network.condition_column is required")
        if not 0.5 <= self.multiview.as_weight <= 1.0:
            raise ValueError("multiview.as_weight must be within [0.5, 1.0]")
        if len(self.multiview.auxiliary_columns) != len(
            set(self.multiview.auxiliary_columns)
        ):
            raise ValueError("multiview.auxiliary_columns cannot contain duplicates")
        if self.multiview.enabled and (
            self.multiview.event_feature_table is None
            or not self.multiview.auxiliary_columns
        ):
            raise ValueError(
                "Enabled multiview fusion requires event_feature_table and auxiliary_columns"
            )
        if any(value < 0 for value in self.multiview.auxiliary_weights.values()):
            raise ValueError("multiview auxiliary weights cannot be negative")
        if not 0.0 <= self.multiview.edge_threshold <= 1.0:
            raise ValueError("multiview.edge_threshold must be within [0, 1]")
        if self.multiview.min_module_size < 2:
            raise ValueError("multiview.min_module_size must be at least 2")
        if self.multiview.resolution <= 0:
            raise ValueError("multiview.resolution must be positive")
        if self.bayesian_preservation.prior_sd <= 0:
            raise ValueError("bayesian_preservation.prior_sd must be positive")
        if self.bayesian_preservation.tau_scale <= 0:
            raise ValueError("bayesian_preservation.tau_scale must be positive")
        if self.bayesian_preservation.grid_size < 51:
            raise ValueError("bayesian_preservation.grid_size must be at least 51")
        if self.validation.bootstrap_iterations < 1:
            raise ValueError("validation.bootstrap_iterations must be positive")
        if not 0.0 <= self.validation.moderate_jaccard <= self.validation.robust_jaccard <= 1.0:
            raise ValueError("validation Jaccard thresholds must satisfy 0 <= moderate <= robust <= 1")
        if not self.validation.knn_neighbors or any(value < 2 for value in self.validation.knn_neighbors):
            raise ValueError("validation.knn_neighbors must contain values >= 2")
        if not self.validation.community_resolutions or any(value <= 0 for value in self.validation.community_resolutions):
            raise ValueError("validation.community_resolutions must contain positive values")
        if not isinstance(self.consensus.apply_threshold_labels, bool):
            raise ValueError("consensus.apply_threshold_labels must be true or false")
        if self.consensus.bootstrap_iterations < 2:
            raise ValueError("consensus.bootstrap_iterations must be at least 2")
        if self.consensus.bootstrap_iterations > 65_535:
            raise ValueError("consensus.bootstrap_iterations exceeds uint16 co-assignment storage")
        if not 0.0 < self.consensus.min_coassignment <= 1.0:
            raise ValueError("consensus.min_coassignment must be within (0, 1]")
        if self.consensus.min_module_size is not None and self.consensus.min_module_size < 2:
            raise ValueError("consensus.min_module_size must be at least 2 when provided")
        if self.consensus.clustering_method not in {
            "backbone_confidence_only",
            "backbone_refined_consensus",
            "dynamic_tree_consensus",
            "hierarchical_consensus",
            "leiden_graph",
            "hierarchical",
        }:
            raise ValueError(
                "consensus.clustering_method must be backbone_confidence_only, "
                "backbone_refined_consensus, "
                "dynamic_tree_consensus, hierarchical_consensus, leiden_graph "
                "or hierarchical"
            )
        if not 0.0 < self.consensus.dynamic_cut_height <= 1.0:
            raise ValueError(
                "consensus.dynamic_cut_height must be within (0, 1]"
            )
        if not self.consensus.dynamic_deep_splits or any(
            value not in {0, 1, 2, 3, 4}
            for value in self.consensus.dynamic_deep_splits
        ):
            raise ValueError(
                "consensus.dynamic_deep_splits must contain values within 0..4"
            )
        if not self.consensus.hierarchical_thresholds or any(
            not 0.0 < value <= 1.0
            for value in self.consensus.hierarchical_thresholds
        ):
            raise ValueError(
                "consensus.hierarchical_thresholds must contain values within (0, 1]"
            )
        if not self.consensus.graph_edge_thresholds or any(
            not 0.0 <= value <= 1.0 for value in self.consensus.graph_edge_thresholds
        ):
            raise ValueError("consensus.graph_edge_thresholds must contain values within [0, 1]")
        if not self.consensus.graph_resolutions or any(
            value <= 0.0 for value in self.consensus.graph_resolutions
        ):
            raise ValueError("consensus.graph_resolutions must contain positive values")
        if self.consensus.graph_top_k < 2:
            raise ValueError("consensus.graph_top_k must be at least 2")
        if not 1 <= self.consensus.min_candidate_modules <= self.consensus.max_candidate_modules:
            raise ValueError(
                "consensus module limits must satisfy 1 <= min_candidate_modules <= max_candidate_modules"
            )
        if not 0.0 <= self.consensus.peripheral_membership_min <= self.consensus.core_membership_min <= 1.0:
            raise ValueError(
                "consensus membership thresholds must satisfy 0 <= peripheral <= core <= 1"
            )
        if not 0.0 <= self.consensus.peripheral_membership_margin <= 1.0:
            raise ValueError("consensus.peripheral_membership_margin must be within [0, 1]")
        for name in (
            "core_kme_mean_min",
            "core_kme_lower_quantile_min",
            "core_kme_sign_consistency_min",
            "peripheral_kme_mean_min",
            "peripheral_kme_lower_quantile_min",
            "peripheral_kme_sign_consistency_min",
        ):
            if not 0.0 <= getattr(self.consensus, name) <= 1.0:
                raise ValueError(f"consensus.{name} must be within [0, 1]")
        if self.consensus.peripheral_kme_mean_min > self.consensus.core_kme_mean_min:
            raise ValueError(
                "consensus.peripheral_kme_mean_min cannot exceed core_kme_mean_min"
            )
        if (
            self.consensus.peripheral_kme_lower_quantile_min
            > self.consensus.core_kme_lower_quantile_min
        ):
            raise ValueError(
                "consensus.peripheral_kme_lower_quantile_min cannot exceed "
                "core_kme_lower_quantile_min"
            )
        if (
            self.consensus.peripheral_kme_sign_consistency_min
            > self.consensus.core_kme_sign_consistency_min
        ):
            raise ValueError(
                "consensus.peripheral_kme_sign_consistency_min cannot exceed "
                "core_kme_sign_consistency_min"
            )
        if not 0.0 <= self.consensus.max_grey_fraction < 1.0:
            raise ValueError("consensus.max_grey_fraction must be within [0, 1)")
        if not 0.0 < self.consensus.max_largest_module_fraction <= 1.0:
            raise ValueError("consensus.max_largest_module_fraction must be within (0, 1]")
        if self.consensus.dense_max_events < 2:
            raise ValueError("consensus.dense_max_events must be at least 2")
        if not 0.0 <= self.consensus.output_edge_threshold <= 1.0:
            raise ValueError("consensus.output_edge_threshold must be within [0, 1]")
        if self.consensus.max_output_edges < 1:
            raise ValueError("consensus.max_output_edges must be positive")
        cohort_names = [cohort.name.strip() for cohort in self.external_cohorts]
        if any(not name for name in cohort_names):
            raise ValueError("external_cohorts names cannot be blank")
        if len(cohort_names) != len(set(cohort_names)):
            raise ValueError("external_cohorts names must be unique")
        multiobjective = self.multiobjective_selection
        for label, values in (
            ("tuning_cohorts", multiobjective.tuning_cohorts),
            ("holdout_cohorts", multiobjective.holdout_cohorts),
        ):
            if any(not str(value).strip() for value in values):
                raise ValueError(f"multiobjective_selection.{label} cannot contain blanks")
            if len(values) != len(set(values)):
                raise ValueError(f"multiobjective_selection.{label} cannot contain duplicates")
        unknown_selection_cohorts = (
            set(multiobjective.tuning_cohorts)
            | set(multiobjective.holdout_cohorts)
        ).difference(cohort_names)
        if unknown_selection_cohorts:
            raise ValueError(
                "multiobjective_selection references unknown external cohorts: "
                f"{sorted(unknown_selection_cohorts)}"
            )
        overlap = set(multiobjective.tuning_cohorts) & set(
            multiobjective.holdout_cohorts
        )
        if overlap:
            raise ValueError(
                "multiobjective_selection tuning and holdout cohorts must be disjoint: "
                f"{sorted(overlap)}"
            )
        objective_weights = (
            multiobjective.stability_weight,
            multiobjective.preservation_weight,
            multiobjective.functional_weight,
        )
        if any(value <= 0 for value in objective_weights):
            raise ValueError(
                "multiobjective_selection objective weights must all be positive"
            )
        if multiobjective.functional_weight_mode not in {
            "auto_significance",
            "always",
            "disabled",
        }:
            raise ValueError(
                "multiobjective_selection.functional_weight_mode must be "
                "auto_significance, always, or disabled"
            )
        if multiobjective.min_significant_functional_modules < 1:
            raise ValueError(
                "multiobjective_selection.min_significant_functional_modules must be at least 1"
            )
        if multiobjective.min_shared_events_per_module < 3:
            raise ValueError(
                "multiobjective_selection.min_shared_events_per_module must be at least 3"
            )
        if multiobjective.min_module_genes < 2:
            raise ValueError(
                "multiobjective_selection.min_module_genes must be at least 2"
            )
        if not 0.0 < multiobjective.functional_fdr_cutoff <= 1.0:
            raise ValueError(
                "multiobjective_selection.functional_fdr_cutoff must be within (0, 1]"
            )
        if not 0.0 < multiobjective.score_floor < 1.0:
            raise ValueError(
                "multiobjective_selection.score_floor must be within (0, 1)"
            )
        if multiobjective.enabled:
            if not multiobjective.tuning_cohorts:
                raise ValueError(
                    "Enabled multiobjective_selection requires at least one tuning cohort"
                )
            if multiobjective.functional_gene_sets_gmt is None:
                raise ValueError(
                    "Enabled multiobjective_selection requires functional_gene_sets_gmt"
                )
        if self.preservation.permutations < 10:
            raise ValueError("preservation.permutations must be at least 10")
        if self.preservation.min_shared_events < 3:
            raise ValueError("preservation.min_shared_events must be at least 3")
        if self.preservation.dense_max_events < self.preservation.min_shared_events:
            raise ValueError("preservation.dense_max_events must be >= min_shared_events")
        if self.preservation.partial_zsummary >= self.preservation.preserved_zsummary:
            raise ValueError("preservation partial_zsummary must be lower than preserved_zsummary")
        if self.hybrid_programs.matrix_source not in {"auto", "observed", "residual"}:
            raise ValueError("hybrid_programs.matrix_source must be auto, observed or residual")
        if self.hybrid_programs.diffuse_engine not in {"auto", "wgcna_exact", "python_soft_tom"}:
            raise ValueError(
                "hybrid_programs.diffuse_engine must be auto, wgcna_exact or python_soft_tom"
            )
        if self.hybrid_programs.bootstrap_iterations < 2:
            raise ValueError("hybrid_programs.bootstrap_iterations must be at least 2")
        if (
            self.hybrid_programs.min_module_size is not None
            and self.hybrid_programs.min_module_size < 2
        ):
            raise ValueError("hybrid_programs.min_module_size must be at least 2 when provided")
        if not 0.0 < self.hybrid_programs.python_tom_cut_height < 1.0:
            raise ValueError("hybrid_programs.python_tom_cut_height must be within (0, 1)")
        if self.hybrid_programs.python_max_events < self.network.min_module_size * 2:
            raise ValueError(
                "hybrid_programs.python_max_events must be at least twice network.min_module_size"
            )
        if not (
            0.0
            <= self.hybrid_programs.moderate_confidence
            <= self.hybrid_programs.robust_confidence
            <= 1.0
        ):
            raise ValueError(
                "hybrid confidence thresholds must satisfy 0 <= moderate <= robust <= 1"
            )
        allowed_benchmark_methods = {
            "katana_sparse",
            "pearson_hierarchical",
            "graphical_lasso",
            "wgcna_exact",
        }
        unknown_benchmark_methods = set(self.benchmark.methods).difference(allowed_benchmark_methods)
        if unknown_benchmark_methods:
            raise ValueError(f"Unknown benchmark methods: {sorted(unknown_benchmark_methods)}")
        if not self.benchmark.methods or len(self.benchmark.methods) != len(set(self.benchmark.methods)):
            raise ValueError("benchmark.methods must contain unique method names")
        if self.benchmark.max_events < self.network.min_module_size * 2:
            raise ValueError("benchmark.max_events must be at least twice network.min_module_size")
        if self.benchmark.bootstrap_iterations < 1:
            raise ValueError("benchmark.bootstrap_iterations must be positive")
        if self.benchmark.graphical_lasso_alpha <= 0:
            raise ValueError("benchmark.graphical_lasso_alpha must be positive")
        if self.benchmark.graphical_lasso_max_iter < 10:
            raise ValueError("benchmark.graphical_lasso_max_iter must be at least 10")
        if self.benchmark.neutral_knn_neighbors < 2:
            raise ValueError("benchmark.neutral_knn_neighbors must be at least 2")
        if not 0.0 < self.benchmark.wgcna_merge_cut_height < 1.0:
            raise ValueError("benchmark.wgcna_merge_cut_height must be within (0, 1)")
        if self.interpretation.min_module_genes < 2:
            raise ValueError("interpretation.min_module_genes must be at least 2")
        if not 0.0 < self.interpretation.fdr_cutoff <= 1.0:
            raise ValueError("interpretation.fdr_cutoff must be within (0, 1]")
        if self.interpretation.motif_window_bp < 20:
            raise ValueError("interpretation.motif_window_bp must be at least 20")
        if not 0.0 < self.interpretation.motif_score_threshold <= 1.0:
            raise ValueError("interpretation.motif_score_threshold must be within (0, 1]")
        if self.interpretation.isoform_min_gene_tpm < 0:
            raise ValueError("interpretation.isoform_min_gene_tpm cannot be negative")
        if not 0.0 <= self.interpretation.isoform_min_abs_delta_if <= 1.0:
            raise ValueError("interpretation.isoform_min_abs_delta_if must be within [0, 1]")
        if self.statistics.min_total_reads < 0:
            raise ValueError("statistics.min_total_reads cannot be negative")
        if self.statistics.min_samples_per_condition < 2:
            raise ValueError("statistics.min_samples_per_condition must be at least 2")
        if self.statistics.max_events is not None and self.statistics.max_events < 1:
            raise ValueError("statistics.max_events must be positive when provided")
        if self.statistics.optimizer_max_iterations < 20:
            raise ValueError("statistics.optimizer_max_iterations must be at least 20")
        if self.statistics.quadrature_points < 5:
            raise ValueError("statistics.quadrature_points must be at least 5")
        if not 0.0 < self.statistics.fdr_cutoff <= 1.0:
            raise ValueError("statistics.fdr_cutoff must be within (0, 1]")
        if self.meta_analysis.min_shared_events < 3:
            raise ValueError("meta_analysis.min_shared_events must be at least 3")
        if self.meta_analysis.min_cohorts < 2:
            raise ValueError("meta_analysis.min_cohorts must be at least 2")
        if not 0.0 < self.meta_analysis.fdr_cutoff <= 1.0:
            raise ValueError("meta_analysis.fdr_cutoff must be within (0, 1]")
        if self.meta_analysis.event_matching not in {"exact", "annotation_aware"}:
            raise ValueError("meta_analysis.event_matching must be exact or annotation_aware")
        if not 0.0 < self.meta_analysis.min_supported_cohort_fraction <= 1.0:
            raise ValueError("meta_analysis.min_supported_cohort_fraction must be within (0, 1]")
        contrast = self.factor_contrasts
        if contrast.enabled:
            required = {
                "primary_factor": contrast.primary_factor,
                "target_level": contrast.target_level,
                "reference_level": contrast.reference_level,
                "perturbation_factor": contrast.perturbation_factor,
                "baseline_level": contrast.baseline_level,
            }
            missing = [name for name, value in required.items() if not str(value).strip()]
            if missing:
                raise ValueError(f"factor_contrasts missing required settings: {missing}")
            if contrast.target_level == contrast.reference_level:
                raise ValueError("factor_contrasts target_level and reference_level must differ")
            if not contrast.perturbation_levels:
                raise ValueError("factor_contrasts perturbation_levels cannot be empty when enabled")
            if len(contrast.perturbation_levels) != len(set(contrast.perturbation_levels)):
                raise ValueError("factor_contrasts perturbation_levels cannot contain duplicates")
            if len(contrast.context_factors) != len(set(contrast.context_factors)):
                raise ValueError("factor_contrasts context_factors cannot contain duplicates")
            if contrast.network_source not in {"observed", "residual"}:
                raise ValueError("factor_contrasts network_source must be observed or residual")
            if not 0.0 < contrast.fdr_cutoff <= 1.0:
                raise ValueError("factor_contrasts fdr_cutoff must be within (0, 1]")
            if not 0.0 <= contrast.min_abs_delta_psi <= 1.0:
                raise ValueError("factor_contrasts min_abs_delta_psi must be within [0, 1]")
            if contrast.min_samples_per_group < 2:
                raise ValueError("factor_contrasts min_samples_per_group must be at least 2")

    def to_mapping(self) -> dict[str, Any]:
        mapping = asdict(self)
        for key in [
            "rmats_root",
            "sample_metadata",
            "output_dir",
            "quant_root",
            "gtf",
            "genome_fasta",
            "transcript_fasta",
            "bam_root",
            "support_root",
            "rbp_gene_set",
        ]:
            value = mapping.get(key)
            mapping[key] = str(value) if value is not None else None
        mapping["multiview"]["event_feature_table"] = (
            str(self.multiview.event_feature_table)
            if self.multiview.event_feature_table is not None
            else None
        )
        mapping["gene_prior"]["exon_count_matrix"] = (
            str(self.gene_prior.exon_count_matrix)
            if self.gene_prior.exon_count_matrix is not None
            else None
        )
        mapping["gene_prior"]["feature_gene_map"] = (
            str(self.gene_prior.feature_gene_map)
            if self.gene_prior.feature_gene_map is not None
            else None
        )
        mapping["external_cohorts"] = [
            {
                "name": cohort.name,
                "psi_matrix": str(cohort.psi_matrix),
                "sample_metadata": str(cohort.sample_metadata),
                "module_assignments": str(cohort.module_assignments) if cohort.module_assignments else None,
                "effect_table": str(cohort.effect_table) if cohort.effect_table else None,
                "event_metadata": str(cohort.event_metadata) if cohort.event_metadata else None,
                "event_mapping_table": str(cohort.event_mapping_table) if cohort.event_mapping_table else None,
                "preservation_table": str(cohort.preservation_table) if cohort.preservation_table else None,
                "genome_build": cohort.genome_build,
                "annotation_version": cohort.annotation_version,
            }
            for cohort in self.external_cohorts
        ]
        mapping["multiobjective_selection"]["functional_gene_sets_gmt"] = (
            str(self.multiobjective_selection.functional_gene_sets_gmt)
            if self.multiobjective_selection.functional_gene_sets_gmt is not None
            else None
        )
        for section in ("interpretation",):
            for key, value in mapping[section].items():
                if key.endswith("_matrix") or key.endswith("_hits") or key.endswith("_gmt") or key.endswith("_database"):
                    mapping[section][key] = str(value) if value is not None else None
        return mapping

    def write(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(yaml.safe_dump(self.to_mapping(), sort_keys=False), encoding="utf-8")
        return target


def _path(value: Any, base: Path) -> Path | None:
    if value in (None, ""):
        return None
    path = Path(str(value)).expanduser()
    return path if path.is_absolute() else (base / path).resolve()


def load_config(path: str | Path) -> CoSpliceConfig:
    config_path = Path(path).resolve()
    raw = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    base = config_path.parent
    stable = StableASConfig(**(raw.get("stable_as") or {}))
    matrix_input = MatrixInputConfig(**(raw.get("matrix_input") or {}))
    network = NetworkConfig(**(raw.get("network") or {}))
    uncertainty = UncertaintyConfig(**(raw.get("uncertainty") or {}))
    gene_prior_raw = dict(raw.get("gene_prior") or {})
    gene_prior_raw["exon_count_matrix"] = _path(
        gene_prior_raw.get("exon_count_matrix"), base
    )
    gene_prior_raw["feature_gene_map"] = _path(
        gene_prior_raw.get("feature_gene_map"), base
    )
    gene_prior = GenePriorConfig(**gene_prior_raw)
    differential_network = DifferentialNetworkConfig(
        **(raw.get("differential_network") or {})
    )
    multiview_raw = dict(raw.get("multiview") or {})
    multiview_raw["event_feature_table"] = _path(
        multiview_raw.get("event_feature_table"), base
    )
    multiview = MultiViewConfig(**multiview_raw)
    bayesian_preservation = BayesianPreservationConfig(
        **(raw.get("bayesian_preservation") or {})
    )
    validation = ValidationConfig(**(raw.get("validation") or {}))
    consensus = ConsensusConfig(**(raw.get("consensus") or {}))
    multiobjective_raw = dict(raw.get("multiobjective_selection") or {})
    multiobjective_raw["functional_gene_sets_gmt"] = _path(
        multiobjective_raw.get("functional_gene_sets_gmt"), base
    )
    multiobjective_selection = MultiObjectiveSelectionConfig(
        **multiobjective_raw
    )
    external_cohorts = [
        ExternalCohortConfig(
            name=str(item.get("name", "")).strip(),
            psi_matrix=_path(item.get("psi_matrix"), base) or Path(),
            sample_metadata=_path(item.get("sample_metadata"), base) or Path(),
            module_assignments=_path(item.get("module_assignments"), base),
            effect_table=_path(item.get("effect_table"), base),
            event_metadata=_path(item.get("event_metadata"), base),
            event_mapping_table=_path(item.get("event_mapping_table"), base),
            preservation_table=_path(item.get("preservation_table"), base),
            genome_build=str(item.get("genome_build", "") or ""),
            annotation_version=str(item.get("annotation_version", "") or ""),
        )
        for item in (raw.get("external_cohorts") or [])
    ]
    preservation = PreservationConfig(**(raw.get("preservation") or {}))
    hybrid_programs = HybridProgramConfig(**(raw.get("hybrid_programs") or {}))
    benchmark = BenchmarkConfig(**(raw.get("benchmark") or {}))
    interpretation_raw = dict(raw.get("interpretation") or {})
    for key in ("rbp_expression_matrix", "rbp_motif_hits", "rbp_motif_database", "gene_sets_gmt"):
        interpretation_raw[key] = _path(interpretation_raw.get(key), base)
    interpretation = ModuleInterpretationConfig(**interpretation_raw)
    statistics = ModuleStatisticsConfig(**(raw.get("statistics") or {}))
    meta_analysis = MetaAnalysisConfig(**(raw.get("meta_analysis") or {}))
    factor_contrasts = FactorContrastConfig(**(raw.get("factor_contrasts") or {}))
    config = CoSpliceConfig(
        project_name=str(raw.get("project_name", "")).strip(),
        rmats_root=_path(raw.get("rmats_root"), base),
        sample_metadata=_path(raw.get("sample_metadata"), base) or Path(),
        output_dir=_path(raw.get("output_dir"), base) or (base / "katana_cosplice_output"),
        quant_root=_path(raw.get("quant_root"), base),
        gtf=_path(raw.get("gtf"), base),
        genome_fasta=_path(raw.get("genome_fasta"), base),
        transcript_fasta=_path(raw.get("transcript_fasta"), base),
        bam_root=_path(raw.get("bam_root"), base),
        support_root=_path(raw.get("support_root"), base),
        rbp_gene_set=_path(raw.get("rbp_gene_set"), base),
        genome_build=str(raw.get("genome_build", "") or ""),
        annotation_version=str(raw.get("annotation_version", "") or ""),
        gene_exclusion=dict(raw.get("gene_exclusion") or {}),
        rmats_mode=str(raw.get("rmats_mode", "JCEC")).upper(),
        selected_conditions=[str(value) for value in raw.get("selected_conditions", [])],
        residual_covariates=[str(value) for value in raw.get("residual_covariates", [])],
        association_factors=[str(value) for value in raw.get("association_factors", ["condition"])],
        factorial_factors=[str(value) for value in raw.get("factorial_factors", [])],
        batch_column=str(raw["batch_column"]) if raw.get("batch_column") else None,
        matrix_input=matrix_input,
        stable_as=stable,
        network=network,
        uncertainty=uncertainty,
        gene_prior=gene_prior,
        differential_network=differential_network,
        multiview=multiview,
        bayesian_preservation=bayesian_preservation,
        validation=validation,
        consensus=consensus,
        multiobjective_selection=multiobjective_selection,
        external_cohorts=external_cohorts,
        preservation=preservation,
        hybrid_programs=hybrid_programs,
        benchmark=benchmark,
        interpretation=interpretation,
        statistics=statistics,
        meta_analysis=meta_analysis,
        factor_contrasts=factor_contrasts,
    )
    config.validate()
    return config
