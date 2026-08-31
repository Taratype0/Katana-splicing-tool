# Katana v1.1.0 Release Notes

## CoSplice method architecture

- Replaced the production sparse kNN/Leiden-first network with official R WGCNA signed adjacency, exact TOM, Dynamic Tree Cut and eigengene merging.
- Enabled beta-binomial replicate/read reliability QC by default (`network_reliability >= 0.15`) and use posterior PSI without multiplying event reliability directly into TOM adjacency.
- Changed production bootstrap consensus to confidence-only annotation: the full-data WGCNA partition is preserved exactly, and weak or ambiguous events remain assigned but are flagged.
- Kept module-program merging, kME-style membership refinement, and explicit core/peripheral/uncertain/grey assignments.
- Added a configurable dense-TOM memory cap that fails explicitly rather than silently truncating Stable-AS events.
- Retained sparse/Leiden and tight/diffuse implementations only for historical cache and benchmark compatibility.

## Leakage and objective safeguards

- Holdout runs receive one pre-frozen TOM height and do not perform holdout-side cut selection.
- Multi-objective selection now disables and renormalizes an unevaluable tuning-preservation objective instead of collapsing every candidate score.
- Functional coherence remains inactive unless the predeclared tuning-GMT FDR gate is passed.
- Selection audits record TOM pilot candidates, effective objective weights, activation state, fingerprints, and interpretation limits.

## Frozen public validation

The final 11-run TCGA SpliceSeq panel compared WGCNA-backbone CoSplice with exact R WGCNA on identical event universes:

- mean bootstrap ARI among non-empty partitions: 0.867 versus 0.871;
- mean holdout program-score preservation: 0.497 versus 0.472;
- non-empty module recovery: 11/11 for both methods;
- strict valid multi-module recovery: 2/11 for both methods;
- functional-objective activation: 0/11.

This establishes near-WGCNA internal stability while retaining CoSplice uncertainty, confidence, differential-network and interpretation outputs. It does not establish universal superiority or replicated biological function.
Zero-module/all-grey results are excluded from ARI summaries because a constant partition can otherwise receive a misleading perfect score.

## Real-count uncertainty ablation

A matched three-arm Setd2 ablation used 3,000 Stable-AS events, 36 samples,
replicate IJC/SJC counts, seed 37 and the same 20 condition-stratified
bootstrap samples in every arm:

- raw-PSI WGCNA: 5 modules and mean bootstrap ARI 0.548;
- posterior-PSI WGCNA after read/repeat reliability filtering: 6 modules and
  mean bootstrap ARI 0.564;
- posterior-PSI WGCNA plus conservative co-assignment boundary rejection:
  6 modules and mean bootstrap ARI 0.540.

The read-uncertainty layer therefore produced a small reproducibility gain on
this real-count dataset. Because boundary rejection reduced bootstrap ARI, it
was removed from the production path and retained only for historical
reproduction. Direct reliability multiplication into
TOM is not a production mode because the negative control ablation collapsed
the partition to all grey.

## Independent validation boundary

- In the independent GSE62669 co-splicing challenge, CoSplice and exact WGCNA
  both recovered two valid modules. Mean bootstrap ARI was 0.500 versus 0.497,
  minimum ARI was 0.347 for both, and within-module co-assignment was 0.784
  versus 0.611. CoSplice was substantially slower.
- Two additional boundary datasets (GSE141633 bulk PRAD PSI and GSE85908
  single-cell PSI) returned all-grey partitions for both methods under the
  frozen 500-event settings. These are retained as negative stress tests.
- A post-lock human TCGA KIRC SETD2 analysis tested module preservation,
  GO/KEGG coherence and phenotype association. Eight unadjusted module-score
  comparisons passed Mann-Whitney FDR, but zero survived adjustment for tumor
  purity and stage; zero GO/KEGG terms passed FDR 0.05. The release therefore
  makes no independent SETD2 biological-replication claim.

## Compatibility

- Existing sparse/Leiden configuration values remain accepted for reproducibility of historical runs.
- New projects use `network.module_detection: wgcna_blockwise` and `consensus.clustering_method: backbone_confidence_only`.
- The production dense TOM default is intended for up to 3,000 Stable-AS events unless the user explicitly changes the memory-safety cap.
