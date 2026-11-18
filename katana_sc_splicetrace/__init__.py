"""Katana-SC-SpliceTrace 2.4.1 — Smart-seq3 RBP discovery, mature-mRNA motif, Module-RBP bridge & GTF/UniProt ORF/NMD engine."""
from katana_sc_splicetrace.contrast_discovery import (
    SCMotifContrastResultV2,
    compute_sc_splicing_contrast_v2,
    compute_sc_stratified_meta_contrast_v2,
    discover_sc_rbp_programs_v2,
)
from katana_sc_splicetrace.extensible_db import (
    export_katana_motif_pack,
    iupac_to_pwm,
    load_external_motif_entries,
    merge_motif_entries_into_library,
)
from katana_sc_splicetrace.fast_accel import HAS_NUMBA, configure_server_threads, fast_score_mean_log_odds
from katana_sc_splicetrace.library import CuratedMotifLibrary, load_bundled_human_rbp_library
from katana_sc_splicetrace.module_bridge import bridge_splicenet_modules_to_rbps
from katana_sc_splicetrace.protein_consequences import (
    load_bundled_uniprot_domains,
    load_gtf_coding_transcripts,
    trace_protein_and_nmd_consequences,
)
from katana_sc_splicetrace.rna_forms import IndexedFasta, extract_event_rna_v2
from katana_sc_splicetrace.sequence_and_consequence import (
    build_sc_directional_motif_tensor_v2,
    extract_sc_event_sequences_v2,
    load_default_curated_library_v2,
    trace_sc_mature_motif_consequences_v2,
)
from katana_sc_splicetrace.target_and_psi_models import (
    estimate_sc_psi_eb_baseline_v2,
    gene_equal_psi_mse_v2,
    predict_sc_psi_precision_v2,
    rank_rbp_target_events_v2,
)

__all__ = [
    "__version__",
    "HAS_NUMBA",
    "configure_server_threads",
    "fast_score_mean_log_odds",
    "IndexedFasta",
    "CuratedMotifLibrary",
    "SCMotifContrastResultV2",
    "load_bundled_human_rbp_library",
    "load_default_curated_library_v2",
    "load_external_motif_entries",
    "merge_motif_entries_into_library",
    "export_katana_motif_pack",
    "iupac_to_pwm",
    "extract_event_rna_v2",
    "extract_sc_event_sequences_v2",
    "build_sc_directional_motif_tensor_v2",
    "trace_sc_mature_motif_consequences_v2",
    "compute_sc_splicing_contrast_v2",
    "compute_sc_stratified_meta_contrast_v2",
    "discover_sc_rbp_programs_v2",
    "rank_rbp_target_events_v2",
    "estimate_sc_psi_eb_baseline_v2",
    "predict_sc_psi_precision_v2",
    "gene_equal_psi_mse_v2",
    "load_bundled_uniprot_domains",
    "load_gtf_coding_transcripts",
    "trace_protein_and_nmd_consequences",
    "bridge_splicenet_modules_to_rbps",
]

try:
    from katana_sc_splicetrace.protein_viz import (
        classify_proteome_fate,
        plot_causal_alluvial_flow_v24,
        plot_consequence_landscape_v24,
        plot_dual_track_nmd_architecture_v24,
        plot_lollipop_protein_splicing_map_v24,
        plot_side_by_side_inframe_domains_v24,
    )
    __all__.extend([
        "classify_proteome_fate",
        "plot_consequence_landscape_v24",
        "plot_dual_track_nmd_architecture_v24",
        "plot_side_by_side_inframe_domains_v24",
        "plot_lollipop_protein_splicing_map_v24",
        "plot_causal_alluvial_flow_v24",
    ])
except ImportError:
    pass

__version__ = "2.4.1"
