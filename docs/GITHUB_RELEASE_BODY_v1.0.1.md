# Katana v1.0.1

Katana v1.0.1 is a comparison-first, direction-aware desktop workflow for
RNA-seq, alternative splicing candidate screening and reproducible follow-up.

## What is included

- independent, parameter-fingerprinted analysis modules;
- stable-AS filtering with replicate read-support checks;
- DEG + rMATS integration with optional DEXSeq/DTU support;
- unique-gene ranking, candidate selection and gene/event follow-up;
- event-level cross-comparison overlap, direction and interaction analysis;
- selected-event Native Sashimi;
- Native 5.6 PSI heatmaps, PCA, replicate QC and comparison similarity;
- RBP and isoform support tables and plots;
- standalone and embedded Katana-CoSplice;
- explicit external output directories and guarded legacy-output migration.

## Validation

- Automated Python test suite and offscreen desktop startup smoke test.
- Strict same-input Setd2 Jutils comparison: exact filtered-event agreement
  across six comparisons and about 6.1x lower wall time for Native 5.6.
- Cache reopening reads prior TSV/PNG/JSON artifacts without automatically
  rerunning analysis.

Full Jutils validation:
`docs/JUTILS_HEAD_TO_HEAD_VALIDATION_v1.0.1.md`

## Upgrade

If an older project writes outputs inside the Katana checkout, select an
external directory and use `Settings > Migrate Existing Outputs`. Katana only
migrates recognized output trees and refuses arbitrary input/project folders.

## Known limits

Optional evidence modules require their corresponding project resources.
Native 5.6 is functionally validated for the documented rMATS scope but does
not reproduce every legacy Jutils CLI/MDS option.
