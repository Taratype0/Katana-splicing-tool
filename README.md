# Katana Next

Katana Next is a local-first RNA-seq evidence workspace for bulk expression,
alternative splicing, transcript usage, read validation, co-splicing programs,
and donor-aware single-cell splicing. It replaces the numbered-page Katana
interface with one artifact-driven analysis journey:

`Project -> Quality -> Landscape -> Differential -> Prioritize -> Validate -> Programs -> Single Cell -> Reports`

The previous application is preserved at tag `v1.2.3-legacy`. Katana Next is
developed independently on branch `codex/v2-platform-rewrite`; project results
must be written to a user-selected directory outside this repository.

## Production inventory

| Workspace | Workflows | Scientific purpose |
| --- | --- | --- |
| Project | Input and design audit | Validate files, sample identity, contrasts, references, and direction contracts. |
| Quality | Stable-AS evidence | Separate event measurability and replicate consistency from statistical significance. |
| Landscape | Bulk splicing landscape | Review event composition, PSI distributions, sample similarity, and PCA. |
| Differential | AS, DEG, DTU, cross-comparison | Test or import effects at the correct event, gene, and transcript units. |
| Prioritize | Evidence, enrichment, RBP, isoform | Keep evidence dimensions explicit while prioritizing genes and mechanisms. |
| Validate | Native Sashimi | Inspect real BAM coverage, junction reads, replicate PSI, and aligned exon models. |
| Programs | CoSplice | Discover signed WGCNA/TOM co-splicing programs and annotate bootstrap confidence. |
| Single Cell | Audit, quantify, differential | Quantify junction choices and test donor-pseudobulk effects without cell-level pseudoreplication. |
| Reports | Evidence package | Export HTML, PDF, Parquet, plots, parameters, checksums, provenance, and ZIP. |

All 17 workflows are registered through one runtime. Opening or switching a
workspace only reads the catalog. `Run analysis` launches only the selected
workflow; it never runs another module implicitly. Unchanged inputs and
parameters return an explicit cache hit.

## Install and launch

```powershell
cd E:\Project\katana-next
.\scripts\bootstrap_katana_next.ps1
.\scripts\launch_katana_next_desktop.ps1
```

Browser mode is also available:

```powershell
.\scripts\launch_katana_next.ps1
```

Then open `http://127.0.0.1:8765/`. The desktop and browser modes render the
same React interface and use the same typed FastAPI workflow API.

## Project migration

Legacy `.katana_project.yaml` projects are migrated explicitly. The command
requires an external output directory and never reuses the source checkout as a
result cache:

```powershell
.\.venv-next\Scripts\python.exe -m katana_next.migration `
  E:\path\to\.katana_project.yaml `
  E:\path\to\katana-next-results `
  --destination E:\path\to\katana.project.json
```

The generated `katana.project.json` records sample/contrast direction,
available BAM/BAI and quantification files, references, source assets, and
detected external runtimes.

## Artifact contract

Every run records its normalized parameters, parent artifacts, input and
manifest fingerprints, workflow version, row count, content checksum, warnings,
and claim boundary. Tables are stored as compressed Parquet, indexed by a
project-local DuckDB catalog. Every run directory includes `parameters.json`
with the generated artifact list.

The UI provides responsive Vega or image preview, searchable and sortable
50-row table pages, named gene/event selections, export links, progress and
cancel state, and cache-hit versus rebuild messages. Long comparison lists wrap
horizontally instead of shrinking controls below readable width.

## Scientific boundaries

- rMATS FDR and dPSI define caller-specific differential evidence; IJC/SJC,
  replicate PSI, leave-one-replicate-out consistency, and posterior uncertainty
  independently describe measurement support.
- Exact cross-comparison overlap uses the same normalized event structure.
  Direction concordance and effect similarity are separate endpoints.
- CoSplice uses a signed WGCNA soft-adjacency/TOM backbone. Bootstrap
  co-assignment can reject unstable boundary events to grey but cannot
  repartition stable modules; phenotype, GO, RBP, DEG, and DTU evidence do not
  create module membership.
- RBP and isoform layers prioritize plausible support and never establish
  causal regulation or protein-domain disruption on their own.
- Single-cell confirmatory tests use biological donor/sample pseudobulk. Cells
  are measurement units, not independent biological replicates.

See [Scientific methods](docs/v2/SCIENTIFIC_METHODS.md) and
[Product architecture](docs/v2/PRODUCT_ARCHITECTURE.md) for the complete
contracts.

## Verification

Run the publication acceptance fixture in a new external directory:

```powershell
.\.venv-next\Scripts\python.exe -m scripts.run_v2_acceptance `
  E:\Project\katana-next-runs\my-clean-acceptance
```

The runner executes every workflow that can be validated with deterministic
synthetic inputs, checks that no Run creates artifacts for another workflow,
reopens the project in a fresh runtime, verifies every unchanged workflow is a
cache hit, and writes `acceptance_summary.json`. Native Sashimi is validated
separately with a real indexed BAM.

Developer checks:

```powershell
.\.venv-next\Scripts\python.exe -m compileall -q katana_next
.\.venv-next\Scripts\python.exe -m pytest tests\v2 -q
pnpm --dir apps/web test
pnpm --dir apps/web build
.\.venv-next\Scripts\python.exe -m katana_next.desktop --check
```

Passing software acceptance demonstrates reproducible execution and valid
analysis contracts. It does not by itself establish a novel biological or
method-performance claim; those require locked external benchmark cohorts.
