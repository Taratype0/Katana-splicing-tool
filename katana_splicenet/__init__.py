"""Katana-SpliceNet 2.3.2 — Pure-Python Excess-BioVar Signed TOM Co-Splicing Network."""
from katana_splicenet.evaluation import enrich_module_genes_v2, orthogonalize_host_depth, score_modules_v2
from katana_splicenet.network import (
    compute_signed_tom_from_bicor,
    compute_tukey_bicor_matrix,
    consensus_core_audit,
    fit_wgcna_v2,
    pick_soft_threshold_v2,
)
from katana_splicenet.qc_and_selection import select_events_v2, unit_robust_audit_v2

__version__ = "2.3.2"

__all__ = [
    "__version__",
    "unit_robust_audit_v2",
    "select_events_v2",
    "compute_tukey_bicor_matrix",
    "compute_signed_tom_from_bicor",
    "fit_wgcna_v2",
    "pick_soft_threshold_v2",
    "consensus_core_audit",
    "score_modules_v2",
    "orthogonalize_host_depth",
    "enrich_module_genes_v2",
]
