# Katana Next Delivery Gates

## Phase 0: Preserve and specify

- Preserve legacy release and create an isolated rewrite branch/worktree.
- Freeze domain vocabulary, artifact schema, claim vocabulary, and methods.
- Create representative bulk and single-cell fixtures with immutable checksums.

Exit gate: the legacy checkout is unchanged and every proposed production
module has a question, input, output, analysis unit, and validation plan.

## Phase 1: Platform foundation

- Implement typed project manifest, contrast model, artifact registry, cache
  fingerprint, job queue, and progress/event API.
- Implement DuckDB/Parquet result storage and user-selected output roots.
- Build the React workspace shell, context bar, selection basket, table viewer,
  plot viewer, inspector, and job drawer.

Exit gate: a fixture project can be loaded, fingerprinted, cached, filtered,
selected, and reopened without running an analysis implicitly.

## Phase 2: Bulk analysis vertical slice

- Implement input/QC, event catalog, Landscape, stable-event QC, differential
  AS, exact cross-comparison, Native Sashimi, and report export.
- Add DEG and DTU adapters and transparent candidate evidence views.

Exit gate: one real project runs end to end from manifest to report, every plot
links to its source rows, and changing one parameter invalidates only dependent
artifacts.

## Phase 3: CoSplice 2.0

- Implement exact blockwise WGCNA baseline without duplicate dense TOM.
- Implement posterior/bootstrap consensus candidate and grey-event handling.
- Add module-trait tests, preservation, interpretation, and feasibility gates.
- Run locked same-input ablation on simulations, Setd2/Ythdc1, and independent
  public cohorts.

Exit gate: the production engine is selected from evidence. It must not have a
worse calibrated false-positive rate than WGCNA, and any claimed advantage must
replicate on untouched data.

## Phase 4: Single-cell splicing

- Implement assay audit, barcode/UMI junction extraction, junction confidence,
  droplet and full-length quantification tracks, donor-aware pseudobulk tests,
  cell-state analysis, interpretation, and pseudobulk Sashimi.
- Add adapters/comparators for SICILIAN, SpliZ, BRIE2, scQuint and MARVEL where
  license and runtime permit.

Exit gate: simulated truth, matched bulk/pseudobulk, and at least two public
datasets demonstrate junction precision, PSI/usage accuracy, calibrated DAS,
replicate concordance, and reproducible biological interpretation.

## Phase 5: Product acceptance

- Complete keyboard, high-DPI, resizing, multi-monitor, cancellation, crash
  recovery, and large-table tests.
- Generate a full reproducibility bundle and migration tool for legacy projects.
- Remove production imports of legacy pages and experimental engines.

Exit gate: clean-machine installation and a full bulk plus single-cell run pass
without source-tree outputs, hidden jobs, stale-cache reuse, or manual window
resizing.

## Release policy

Each phase is independently shippable but does not inherit a scientific claim
from an earlier prototype. Feature completion requires unit, schema-contract,
integration, scientific-regression, and UI acceptance evidence. A green unit
test suite alone is not evidence that the biological method works.
