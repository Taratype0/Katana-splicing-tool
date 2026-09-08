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
