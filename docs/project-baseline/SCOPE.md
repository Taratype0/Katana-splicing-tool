# Adopted scope and baseline: T01–T08

Status: adopted development scope, 22 October 2026. This is a design and source
inventory checkpoint, not clean-install or integrated scientific acceptance.

| Task | Decision |
|---|---|
| T01 — inventory | Catalogue current analysis implementations, readers, interfaces and dependencies; distinguish new computation from reading existing results. |
| T02 — consolidate | Converge duplicate entry points onto shared cores; exclude superseded historical and unadopted experimental branches from the new delivery. Keep dependencies needed by retained modules. |
| T03 — method boundaries | Keep all adopted functions formal. Track method/source versions and statistical units; default/optional selection is independent of migration, packaging and acceptance. |
| T04 — assays | Bulk upstream accepts conventional short-read FASTQ or eligible genome BAM. Katana downstream has separately qualified bulk, scRNA and snRNA routes. Compatible result import, new calculation and raw-data processing are distinct capabilities. |
| T05 — activation | Separate core defaults, defaults eligible when inputs/design/resources support them, and opt-in calculation. Reusable existing optional results can remain visible without rerunning them. |
| T06 — ownership and dependencies | One Katana analysis core with MotifTrace and Katana-Splice-Network modules; Nextflow upstream, reports, CLI/API and UI have explicit boundaries. |
| T07 — delivery and community | Emphasise downstream interpretation and reuse of upstream components. Separate existing source implementations from the package, containers and Nextflow work still needed. |
| T08 — baseline and recovery | Develop in isolation; preserve the original application, frozen methods and previous results. Use task commits, separately versioned code/method/results/presentation and recoverable baselines. |

## Functional scope

The inventory has 54 formal primary workflow entries. The wider inventory has
134 records: 119 adopted formal records, including supporting functions, and
15 unadopted historical/experimental records. These are inventory units, not
claims of 119 independent algorithms or 54 accepted package workers.

- Core evidence: QC and landscape, differential event/gene/transcript evidence,
  same-comparison integration, candidates, enrichment and isoform interpretation.
- **MotifTrace**: local sequence and motif consequences, usage-aware evidence,
  candidate RBP and CLIP layers, and separately identified native MAPP evidence.
- **Katana-Splice-Network**: discovery, confidence, fixed-module diagnostics,
  condition effects, hub/functional interpretation and external validation.
- Qualified single-cell/nucleus splicing routes: retain assay, donor, cell type,
  aggregation and native inference identity; no implicit single-cell DEG step.
- English report generation and interactive review consume existing results.

WGCNA/TOM, DJExpress, MAPP and other third-party methods keep their original
method identities and attribution. Sequence, binding, statistical association
and mechanistic causation are different evidence levels.

## Defaults and incomplete inputs

The standard profile classifies the 54 primary entries as 5 core defaults,
28 eligible defaults and 21 opt-in entries. Eligibility still depends on the
actual route, available inputs, design and resources; these are not all run on
every project. A complete run resolves the selected profile, not every option.

Network new fits, native MAPP new fits, deep CLIP, complex/calibration models,
Sashimi and advanced single-cell calculations are explicit options. Importing
their compatible results does not request a new fit. Missing eligible defaults
produce explicit omissions; missing required or explicitly selected inputs block
the affected work. Failures preserve successful upstream results. Missing,
failed, cancelled, not estimable and zero significant findings stay distinct.

## Execution boundaries

Prepare/count a sample once and reuse compatible evidence. Schedule independent
contrasts and legal blocks separately; executor/Slurm capacity determines
concurrency, with no arbitrary four-task ceiling. Preserve joint-model scope and
complete correction families. Discover a network once per declared cohort/method/
event universe, then test effects on the fixed network. Native MAPP fitting keeps
its joint sample design. Do not split a scientific model merely to increase jobs.

Opening a project or rendering an English report never starts missing analysis.
The original one-click report chapter structure is retained as a product
requirement. A UI/rendering change reuses valid results instead of realigning reads.

## What remains engineering work

Existing local scientific implementations do not yet constitute the generic new
distribution. Shared headless workers, actual format adapters, resource/dependency
packaging, clean installs, the selected downstream Nextflow DAG and route-specific
acceptance are still required. snRNA support requires explicit qualification;
it is not automatically inferred from a label. Raw scRNA/snRNA FASTQ and raw
long-read reconstruction are outside the initial upstream scope.

Kanata is the provisional Nextflow proposal name; Katana is the package/product.
The current proposal and upstream naming remain subject to community discussion.
No formal nf-core approval or release is claimed by this checkpoint.
