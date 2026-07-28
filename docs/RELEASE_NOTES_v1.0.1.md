# Katana v1.0.1 Release Notes

Release date: 2026-07-28

Katana v1.0.1 consolidates the comparison-first RNA-seq and alternative
splicing workflow into independently runnable, cache-aware modules.

## Highlights

- Stable comparison identity, canonical direction audit and parameter-fingerprinted caches.
- Shared stable-AS filtering from rMATS significance, replicate read support and replicate direction consistency.
- Per-comparison transcript/splicing integration, unique-gene candidate ranking and independent candidate selection.
- Event-level cross-comparison overlap, direction-aware Venn/similarity, dPSI pattern analysis and replicate-aware interaction follow-up.
- Selected-event Native Sashimi with BAM/GTF validation, cacheable event catalogs and PNG preview.
- Native 5.6 sample-PSI heatmaps, PCA, replicate QC, pattern summaries and confidence-aware comparison similarity.
- RBP and isoform candidate support with per-comparison and cross-comparison TSV/PNG outputs.
- Standalone and embedded Katana-CoSplice stages for stable co-splicing program discovery, preservation, interpretation and reporting.
- Explicit external output roots plus guarded one-click migration of recognized legacy Katana outputs.
- Project-level artifact registry and cache-only report/provenance inventory.

## Jutils validation

A strict same-input Setd2 run matched the filtered event count for all six
comparisons. Native 5.6 completed in 56.736 seconds versus 346.371 seconds for
bundled Jutils, about 6.1x faster. See
`docs/JUTILS_HEAD_TO_HEAD_VALIDATION_v1.0.1.md` for scope and limitations.

## Execution contract

Opening a project or page reads lightweight state and existing cache only.
Analysis is started only through an explicit Run, Generate, Apply or
Recalculate action. A downstream module may consume an existing upstream
artifact but does not silently execute the upstream module.

## Upgrade note

Legacy projects whose output path points inside the Katana source checkout are
prompted to choose an external output root. Recognized Katana output trees can
be moved through `Settings > Migrate Existing Outputs`; arbitrary project/input
directories are never treated as migration sources.

## Compatibility and limits

- Native 5.6 is not a pixel-identical replacement for every Jutils CLI option.
- DEXSeq, DTU, RBP, isoform, enrichment and CoSplice interpretation remain optional and depend on project inputs.
- GO/KEGG online enrichment requires network access.
- Manual GUI validation across the user's target monitor/DPI combinations remains a release checklist item.
