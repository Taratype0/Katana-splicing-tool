# Katana Next Rewrite

This directory is the executable design contract for the Katana rewrite. The
legacy application remains preserved at tag `v1.2.3-legacy`; development of the
new platform takes place on branch `codex/v2-platform-rewrite` in
`E:\Project\katana-next`.

The rewrite has three non-negotiable goals:

1. Build one continuous analysis workspace instead of a collection of numbered
   pages.
2. Make statistical validity, stability, and provenance visible in every
   result.
3. Support bulk and single-cell splicing without treating cells as biological
   replicates or forcing one assay chemistry into an unsuitable model.

The implementation is governed by:

- [PRODUCT_ARCHITECTURE.md](PRODUCT_ARCHITECTURE.md): product structure, UI,
  runtime, artifacts, caching, and migration.
- [SCIENTIFIC_METHODS.md](SCIENTIFIC_METHODS.md): bulk AS, CoSplice 2.0,
  single-cell splicing, validation, and claim boundaries.
- [DELIVERY_GATES.md](DELIVERY_GATES.md): phased delivery and objective release
  criteria.

No legacy analysis is copied into production merely because it already exists.
Each module must have a defined biological question, valid experimental unit,
explicit output contract, test fixture, and downstream destination.

## Current status

The beta implementation now registers all nine workspaces and 17 isolated
workflows. It includes typed manifests and contrasts, stable-AS QC, bulk
landscape, DAS/DEG/DTU, exact cross-comparison analysis, transparent candidate
evidence, enrichment/RBP/isoform support, real-BAM Native Sashimi, WGCNA-backed
CoSplice programs, donor-aware single-cell junction analysis, and reproducible
HTML/PDF/ZIP reporting.

Bootstrap and launch the development build:

```powershell
cd E:\Project\katana-next
.\scripts\bootstrap_katana_next.ps1
.\scripts\launch_katana_next.ps1
```

The complete fixture validates software and scientific contracts only and is
intentionally not a biological result. Real Ythdc1 rMATS, stageR DTU, quant.sf,
GTF, BAM and BAI inputs are used separately for adapter and Native Sashimi smoke
acceptance. Publication-level CoSplice and single-cell performance claims still
require locked independent cohorts; green software tests do not replace those
benchmarks.
