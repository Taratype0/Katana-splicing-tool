# Katana v1.0.2 Release Notes

Release date: 2026-08-31

Katana v1.0.2 hardens the CoSplice methods workflow and closes the locked public-validation path.

## CoSplice validation

- Adds leakage-safe discovery/tuning/holdout ablation with pre-holdout fingerprints and access audit.
- Adds exact same-event-universe R WGCNA comparison and a standalone aggregate panel report.
- Adds adaptive functional weighting: GO/KEGG can influence selection only after a formal FDR gate.
- Penalizes degenerate one-module candidates during selection rather than rewarding their high ARI.
- Supports TCGA SpliceSeq ES, RI, AA, AD, AP, AT and ME inputs plus incremental manifest merging.
- Records negative enrichment and preservation results without converting them into a superiority claim.

## Current evidence boundary

In 11 frozen TCGA settings, CoSplice recovered valid multi-module partitions more often than exact
WGCNA and had modestly higher mean program-score preservation, while WGCNA retained higher internal
bootstrap ARI. Traditional Zsummary and functional replication do not establish universal CoSplice
superiority. Cross-study, cross-platform holdouts remain publication-critical.

## Reproducibility

- Adds `katana-cosplice-locked-ablation`.
- Adds `katana-cosplice-locked-panel-report`.
- Every locked child run records inputs, parameters, assignment hashes, random seed and generated files.
