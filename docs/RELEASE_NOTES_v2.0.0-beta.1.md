# Katana Next 2.0.0 beta 1

Katana 2.0 is an isolated rewrite of the legacy `v1.2.3-legacy` desktop tool.
It replaces numbered pages and implicit dependencies with nine question-led
workspaces, typed artifacts, explicit jobs, and a shared responsive viewer.

## Product changes

- Added a React/TypeScript evidence workspace hosted by FastAPI and the native
  PySide6 QtWebEngine desktop launcher.
- Added project manifests, explicit numerator-minus-denominator contrasts,
  DuckDB/Parquet artifacts, named selections, full input/parameter
  fingerprints, and compatible-cache resolution.
- Added isolated progress/cancel jobs. Opening a page never runs analysis and a
  Run cannot launch unrelated workflows.
- Added responsive Vega, PNG/SVG, HTML, PDF and paginated table viewers with
  readable wrapping comparison controls.
- Added legacy-project migration with automatic local Rscript discovery and a
  strict external output-root policy.

## Analysis changes

- Added shared stable-AS QC using replicate IJC/SJC support, PSI uncertainty,
  pairwise consistency, and leave-one-replicate-out direction checks.
- Added bulk PSI Landscape, DAS, DEG, stageR/DTU/DEXSeq import, exact
  cross-comparison overlap, transparent candidate evidence, enrichment, RBP,
  and isoform support.
- Added Native Sashimi with real BAM coverage, count-scaled non-overlapping
  arcs, per-group/replicate PSI, genomic alignment tracks, and event schematics.
- Locked CoSplice production discovery to signed WGCNA/TOM with bootstrap
  co-assignment rejection of unstable boundary events to grey.
- Added single-cell assay audit, CB/UB-aware BAM junction choices, GTF junction
  annotation, donor-cell-type pseudobulk, beta-binomial testing, and FDR.
- Added HTML/PDF/ZIP evidence packages with unique artifact names, parameters,
  checksums, methods, interpretations, and explicit claim limits.

## Performance and correctness

- Normalized rMATS event/replicate data are cached by input and comparison
  fingerprint, so Landscape and DAS no longer reparse the same large tables.
- On the real Ythdc1 project, the first cache build took 156.0 seconds and a
  forced downstream Landscape rebuild reused it in 7.65 seconds.
- Vectorized dominant-transcript selection reduced the real 12-sample isoform
  support rebuild from 145.0 to 17.6 seconds with identical switch identities.
- File fingerprints are cached in memory by resolved path, size, and mtime for
  responsive compatible-cache checks against large BAM collections.
- UI results, parent artifacts, saved selections, and reports are all bound to
  the current project snapshot, preventing stale cross-project reuse.

## Validation boundary

The beta software gate passes on a complete deterministic fixture and a real
Ythdc1 project. CoSplice and single-cell publication-performance claims remain
pending locked independent-cohort benchmarks. See
[`v2/V2_BETA_ACCEPTANCE.md`](v2/V2_BETA_ACCEPTANCE.md).
