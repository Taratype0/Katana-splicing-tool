# Katana Next 2.0.0 beta acceptance

## Decision

The Katana Next beta software gate passes. The rewritten platform is usable as
one continuous, artifact-driven analysis product. This decision covers runtime
isolation, cache identity, output provenance, the implemented statistical
contracts, and verified adapters. It is not a claim that CoSplice or the
single-cell method already outperforms published methods on independent data.

## Deterministic full-product run

The final acceptance run is stored outside the repository:

```text
E:\Project\katana-next-runs\v2-final-acceptance-20260908-r9
```

The fixture contains 15 bulk samples, three declared contrasts, five rMATS
event classes, DEG, DTU, quant.sf, GTF, genome/transcript FASTA, pathway/RBP
sets, and donor-aware single-cell junction counts. It is synthetic truth for
software and method-contract verification, not a biological dataset.

Observed acceptance results:

- nine workspaces, 17 registered workflows, and 101 typed artifacts;
- all 16 fixture-compatible workflows completed independently;
- no Run created artifacts for another workflow;
- every output and `parameters.json` file exists;
- closing and reopening the project created no job;
- all 16 unchanged reruns returned `cache_hit`;
- CoSplice resolves to `backbone_refined_consensus`;
- the report ZIP contains unique paths and no raw inputs;
- every generated artifact is bound to the current project snapshot fingerprint.

Native Sashimi is excluded from the synthetic fixture because fabricated BAM
coverage would not validate a BAM adapter.

## Real Ythdc1 acceptance

The migrated real project is:

```text
E:\Project\katana-next-runs\ythdc1-v2\katana.project.json
```

Verified inputs include four rMATS comparisons, stageR DTU output, sample-level
`quant.sf`, 12 indexed BAM files, and GTF/FASTA references. Project audit,
stable-AS QC, bulk Landscape, DAS, DEG, DTU, exact cross-comparison, transparent
candidate evidence, isoform support, Native Sashimi, and the evidence report run
without hidden upstream execution. The real Sashimi smoke uses one selected
Ptprz1 A5SS event and produces BAM-derived coverage, junction arcs, replicate
PSI, aligned transcript tracks, and an rMATS evidence preview.

On the real project, the first normalized rMATS event-cache build took 156.0
seconds. A forced Landscape rebuild then reused that same input cache in 7.65
seconds (`source_cache_status=hit`), without rerunning another analysis module.
Vectorizing dominant-transcript selection reduced the complete 12-sample
isoform-support rebuild from 145.0 to 17.6 seconds while preserving all 13,512
switch identities.

CoSplice discovery is intentionally rejected for this project because 12
independent samples do not pass the declared 15-sample network-discovery gate.
This is a scientific feasibility result, not a runtime defect.

## Automated checks

```text
Python compileall: passed
Full Python regression suite: 165 passed
Ruff (katana_next, tests/v2, acceptance runner): passed
React/Vitest: 3 passed
TypeScript strict build: passed
Vite production build: passed
Desktop embedded API check: passed
```

The remaining Starlette TestClient warning is emitted by the upstream
FastAPI/httpx compatibility layer and does not alter runtime results.

## UI acceptance

Browser-mode inspection at 1081 x 897 confirms that navigation, scientific
detail, run configuration, horizontally wrapped comparison chips, artifacts,
responsive Vega plots, and real PNG previews stay within the workspace without
manual window refresh. Table pages are limited to 50 rows and use server-side
search, sorting, and pagination.

The PySide6 QtWebEngine desktop host starts the loopback API, verifies health,
loads the same production web bundle, and shuts the server down with the
window. A clean-machine packaged executable and physical multi-monitor sweep
remain release-candidate checks rather than beta blockers.

## Scientific claim boundary

The beta validates implementations and guards. Publication claims still need
predeclared same-input WGCNA ablations, independent biological cohorts, and
assay-matched public single-cell benchmarks. In particular, a synthetic
recovery result cannot establish external module preservation, biological
replication, or superiority over WGCNA, SpliZ, BRIE2, scQuint, or MARVEL.
