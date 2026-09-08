# Katana Next Product Architecture

## 1. Audit conclusion

Katana 1.2.3 contains valuable analysis implementations, but its product model
is page-first rather than question-first. More than forty technical navigation
entries expose implementation details to the user. Large coordinating classes
(`ProjectService` and `ResultsPage`) combine discovery, execution, caching,
plotting, UI state, and interpretation. This makes apparently independent
modules feel disconnected while also allowing hidden coupling between them.

Changing colors or widget spacing would not solve this. Katana Next therefore
uses a new application structure and imports only verified legacy components.

## 2. User journey

The primary workflow is:

`Create project -> Audit data -> Explore -> Test -> Prioritize -> Validate -> Interpret -> Report`

The permanent top-level workspaces are:

| Workspace | Question answered | Main outputs |
| --- | --- | --- |
| Project | What data, samples, references, and contrasts are available? | validated manifest, design matrix, fingerprints |
| Quality | Are counts, junctions, replicates, and annotations usable? | QC gates, exclusions, coverage and replicate diagnostics |
| Landscape | What broad expression and splicing patterns exist? | distributions, sample similarity, PCA, event composition |
| Differential | Which genes, transcripts, and AS events change? | DEG, DTU, DAS tables with effect sizes, CI and FDR |
| Prioritize | Which genes/events have convergent evidence? | transparent evidence matrix and selection basket |
| Validate | Does read coverage support the selected event? | Native Sashimi and junction evidence |
| Programs | Which reproducible co-splicing programs exist? | modules, associations, stability, preservation, interpretation |
| Single Cell | Which cell types/states show reproducible splicing changes? | junction QC, cell-type DAS, state scores and programs |
| Reports | What can be communicated and reproduced? | HTML/PDF bundle, tables, figures, methods and provenance |

Numbered pages such as `5.4.3` are removed. A result is identified by analysis
type, dataset, contrast, parameter fingerprint, and run ID.

## 3. Interaction model

Every analysis workspace uses the same composition:

- A compact context bar shows dataset, contrast, analysis unit, thresholds,
  cache state, and provenance.
- A responsive central canvas links plots and a virtualized result table.
- A right inspector explains the selected gene, event, module, or cell type and
  shows valid next actions.
- A persistent selection basket can send artifacts to Sashimi, candidate
  evidence, enrichment, RBP, isoform, program, or report workflows.
- A bottom job drawer shows queued/running/completed/failed jobs, progress,
  cache hit versus rebuild, logs, cancel, and retry.

Opening a workspace never starts an analysis. It reads compatible cached
artifacts only. Running an analysis creates one explicit job and never triggers
unrelated modules.

## 4. Modern UI direction

Katana Next is a local-first web desktop application:

- React and TypeScript for the interface.
- TanStack Table for searchable, sortable, filterable, virtualized tables.
- Plotly or Vega-Lite for responsive scientific plots and linked brushing.
- Cytoscape for module and regulator networks.
- A Python API and worker process for analysis.
- FastAPI for typed local IPC and server-sent job progress.
- DuckDB plus Parquet/Arrow for large tables and cached matrices.
- Tauri as the Windows desktop shell after browser-mode acceptance.

This is preferred over extending the current PySide widget tree because it
supports responsive previews, high-DPI rendering, reusable viewers, and large
interactive tables without embedding Matplotlib independently on every page.
Streamlit-style frameworks are rejected because they do not provide the state,
job isolation, or desktop packaging control required here.

Visual direction:

- Light scientific workspace by default, with a neutral warm background,
  high-contrast ink, and one restrained accent family.
- Typography, spacing, and color tokens are shared across every module.
- Tables and plots use semantic colors consistently: direction, significance,
  support, uncertainty, and missingness never reuse the same encoding.
- Plot titles state the biological comparison; subtitles state the unit and
  filter; captions state the denominator and test.
- PNG is an export format, not the primary viewer. Interactive SVG/canvas plots
  resize immediately without requiring a window refresh.

## 5. Runtime architecture

```text
katana-next/
  apps/
    desktop/                 # Tauri shell
    web/                     # React/TypeScript UI
    api/                     # FastAPI composition root
  packages/
    domain/                  # immutable project, sample, contrast, artifact models
    workflows/               # typed analysis DAGs and feasibility gates
    bulk_splicing/           # event catalog, QC, DAS and comparison methods
    cosplice/                # network discovery and validation engines
    single_cell_splicing/    # junction, cell-type and state workflows
    expression/              # DEG and enrichment adapters
    isoform/                 # DTU and isoform evidence
    interpretation/          # RBP, motif, pathway and report evidence
    visualization/           # plot specifications, not GUI widgets
  tests/
    unit/
    contract/
    integration/
    scientific/
    ui/
```

Domain packages do not import UI code. Workflows declare inputs and outputs;
workers execute them; viewers render their artifacts. External tools are
adapters with version checks and normalized output contracts.

## 6. Artifact and cache contract

Every artifact records:

- project and run ID;
- parent artifact IDs;
- input file size, modification time, and content fingerprint;
- genome build and annotation version;
- sample and contrast IDs;
- analysis unit and statistical model;
- normalized parameters and random seed;
- software and external-tool versions;
- schema version, row count, checksums, and creation status.

Cache compatibility is determined from the complete artifact fingerprint, not
from a filename or page. Changing a threshold invalidates only dependent
artifacts. Results can be moved to a user-selected output root without writing
large caches to the source tree or system drive.

## 7. Cross-module continuity

The universal objects are `Gene`, `Transcript`, `SpliceEvent`, `Junction`,
`Module`, `Sample`, `CellType`, `Contrast`, `GeneSet`, and `SelectionSet`.
Actions operate on these objects rather than page-specific rows. Examples:

- A DAS event can be sent directly to Native Sashimi.
- A module's hub events can be sent to event validation.
- A candidate gene set can be sent to GO/KEGG or RBP enrichment.
- A cell-type-specific event can be compared with its bulk evidence.
- Any filtered table can become a named, versioned selection set.

The inspector always distinguishes observation, statistical result,
supporting evidence, interpretation, and unsupported causal claim.

## 8. Existing functions to keep, redesign, or retire

| Area | Decision | Reason |
| --- | --- | --- |
| Project scanning and comparison identity | Redesign and keep | strong foundation, but move to manifest-first typed models |
| Stable-AS filtering | Rebuild as shared QC service | one definition must feed all downstream analyses |
| DEG, DTU, rMATS import | Keep through adapters | valuable evidence, but never merge p-values into an opaque score |
| Candidate ranking | Redesign | expose evidence dimensions and Pareto/tier logic |
| Native Sashimi | Keep and modularize | useful event-level read validation |
| Cross-comparison overlap | Rebuild | explicitly separate exact event, direction concordance, locus, and gene |
| RBP/motif and enrichment | Keep as interpretation | matched background, FDR, and non-causal language required |
| Isoform support | Keep and strengthen | distinguish DTU, dominant switch, sequence, ORF and domain evidence |
| Jutils/native Jutils pages | Retire | replace with unified Landscape and similarity workflows |
| Legacy duplicate/benchmark engines | Exclude from production | retain only in benchmark package or archive |
| CoSplice | Rebuild behind a stable API | current implementation has duplicated dense computation and mixed engines |

## 9. Reporting

Reports are assembled from immutable artifacts, not screenshots of pages. Each
figure carries its source table and plot specification. A report bundle contains
HTML, PDF, SVG/PNG, TSV/Parquet, methods, parameter JSON, software versions,
input fingerprints, warnings, and a machine-readable manifest.
