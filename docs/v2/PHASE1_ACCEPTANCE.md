# Phase 1 Platform Acceptance

## Scope

Phase 1 implements the platform foundation specified in
`DELIVERY_GATES.md`. Its fixture is deliberately non-biological. Passing this
gate means the application runtime is ready to receive validated scientific
workflows; it does not validate a DAS or CoSplice method.

## Implemented contracts

- Immutable Pydantic project, sample, contrast, artifact, selection, and job
  models.
- Manifest path resolution and source-tree output rejection.
- Full SHA-256 for ordinary inputs and an explicitly labeled sampled hash for
  files larger than 64 MiB.
- Complete workflow fingerprints over workflow version, project fingerprint,
  input fingerprints, and normalized parameters.
- Project-local DuckDB catalog and Parquet result tables under the user-selected
  output root.
- Explicit asynchronous jobs with progress, cancellation state, cache hit,
  rebuild, failure, and completion semantics.
- Artifact deduplication across forced rebuilds.
- FastAPI project, artifact, table-page, plot, selection, job, and SSE routes.
- React workspace navigation, context bar, responsive plot/table viewer,
  server-side search/sort/paging, virtualized rows, inspector, selection basket,
  and job drawer.
- Automatic recent-project reopen that reads cached artifacts but creates no
  analysis job.

## Automated evidence

Executed on 2026-09-08:

```text
Katana Next Python tests: 5 passed
Katana Next React tests: 3 passed
TypeScript strict build: passed
Vite production build: passed
Legacy Python tests in lightweight venv: 154 passed; 3 correctly required R/WGCNA
Those 3 tests in katana-splicing-tool Conda environment: 3 passed
Python compileall for katana_next: passed
```

The Starlette test client emits one upstream deprecation warning. It does not
affect runtime behavior; migration to the replacement HTTP test client is a
dependency-maintenance task rather than a suppressed test failure.

## Real-process acceptance

The acceptance fixture was created outside the repository at:

```text
E:\Project\katana-next-runs\phase1-acceptance-20260908
```

Observed behavior:

1. Opening the manifest created zero jobs and initially read zero artifacts.
2. Explicit Run generated one Parquet table and one Vega-Lite plot artifact.
3. Searching `Afg2a` returned one of eight event rows.
4. The selected event persisted as a named selection in DuckDB.
5. Reloading the browser automatically reopened the manifest, recovered two
   artifacts and one selection, and showed no job in the new runtime session.
6. Force rebuild retained exactly two artifacts rather than accumulating stale
   catalog rows.
7. The next ordinary Run returned `cache_hit` without rebuilding.
8. At the 1280 x 720 acceptance viewport, the plot, inspector, and activity bar
   remained visible without window resizing.

Generated acceptance artifacts remain outside source control. They include the
DuckDB catalog, exported artifact manifest, compressed Parquet table, and
Vega-Lite JSON specification.

## Phase 1 release decision

Phase 1 passes its foundation gate. Phase 2 may begin with the canonical bulk
event catalog and stable-event QC vertical slice. Planned navigation remains
disabled until its corresponding scientific workflow passes its own gate.
