"""Katana-SpliceTrace 2.4.3 — High-speed bulk AS RBP discovery, mature-mRNA motif, Module-RBP bridge & GTF/UniProt ORF/NMD engine."""
from katana_splicetrace.continuous_rescue import MotifRescueResultV22, run_bulk_continuous_rescue_v22
from katana_splicetrace.extensible_db import (
    export_katana_motif_pack,
    iupac_to_pwm,
    load_external_motif_entries,
    merge_motif_entries_into_library,
)
from katana_splicetrace.fast_accel import HAS_NUMBA, configure_server_threads, fast_score_mean_log_odds
from katana_splicetrace.library import (
    CuratedMotifLibrary,
    load_bundled_human_rbp_library,
    load_curated_human_rbp_library,
)
from katana_splicetrace.matching import match_events_v2
from katana_splicetrace.module_bridge import bridge_splicenet_modules_to_rbps
from katana_splicetrace.protein_consequences import (
    load_bundled_uniprot_domains,
    load_gtf_coding_transcripts,
    trace_protein_and_nmd_consequences,
)
from katana_splicetrace.rmats_io import load_rmats_measurable_events
from katana_splicetrace.rna_forms import IndexedFasta, extract_event_rna_v2, trace_motif_instances_v2
from katana_splicetrace.scanner import run_motiftrace_v2_scan

__all__ = [
    "__version__",
    "HAS_NUMBA",
    "configure_server_threads",
    "fast_score_mean_log_odds",
    "CuratedMotifLibrary",
    "MotifRescueResultV22",
    "IndexedFasta",
    "load_bundled_human_rbp_library",
    "load_curated_human_rbp_library",
    "load_external_motif_entries",
    "merge_motif_entries_into_library",
    "export_katana_motif_pack",
    "iupac_to_pwm",
    "load_rmats_measurable_events",
    "match_events_v2",
    "extract_event_rna_v2",
    "trace_motif_instances_v2",
    "load_bundled_uniprot_domains",
    "load_gtf_coding_transcripts",
    "trace_protein_and_nmd_consequences",
    "bridge_splicenet_modules_to_rbps",
    "run_motiftrace_v2_scan",
    "run_bulk_continuous_rescue_v22",
]

try:
    from katana_splicetrace.protein_viz import (
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

__version__ = "2.4.3"
