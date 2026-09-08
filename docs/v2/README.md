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

Phase 1 is implemented as a runnable vertical slice. It includes the typed
manifest and contrast domain, file/workflow fingerprints, a DuckDB artifact and
selection catalog, Parquet result paging, explicit background jobs, cache-hit
semantics, FastAPI endpoints, and the modern React workspace shell.

Bootstrap and launch the development build:

```powershell
cd E:\Project\katana-next
.\scripts\bootstrap_katana_next.ps1
.\scripts\launch_katana_next.ps1
```

The Phase 1 fixture validates platform behavior only and is intentionally not a
biological analysis. Phase 2 replaces that fixture with the first verified bulk
splicing vertical slice.
