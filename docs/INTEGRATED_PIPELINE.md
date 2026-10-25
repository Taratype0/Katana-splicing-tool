# Integrated bulk analysis and report

Start with [Choose a Katana entry point](ENTRY_POINTS.md) for the concise input,
command and current acceptance map.

The workflow connects existing upstream modules, Katana computations and the
approved seven-section report. It does not change the report layout, frozen
MotifTrace algorithms or existing study results. Generated bulk FASTQ/BAM
complete execution and unchanged cache repeats passed the
[recorded checks](validation/entrypoints-runtime-20261025.md). Site deployment
and biological suitability for a new experiment remain separate. It is not an
approved nf-core pipeline.

## Execution

Use the accepted workflow on `main`. Further development continues on
`codex/package-development`. For container execution, the host
needs Nextflow 25.10.4 or newer, compatible Java and a container engine. Native
acceptance uses Linux x86-64 Docker; Apptainer/Singularity profiles need site
qualification. Build or load the matching Katana downstream
image using the [container instructions](../containers/downstream/README.md).
That image supplies Python, R, RSeQC, FIMO and the other native dependencies;
upstream nf-core tool modules keep their own tool images. The host does not need
a standalone Katana Conda environment for containerized tasks. References and
analysis databases remain explicit inputs.

For direct execution without containers, install the [standalone Linux
environment](RUNTIME_ENVIRONMENTS.md#standalone-linux-environment). A source-only
`pip install -e '.[flow,next]'` does not install R or the native toolchain.
The runtime and actual-image transport have [passed installation checks](validation/downstream-runtime-20261023.md);
full scientific workflow acceptance has its own scope and evidence.

Use a prepared project with completed upstream outputs to avoid repeating
sequencing processing. Set `KATANA_IMAGE` to the newly built, matching image tag
or immutable image reference. The previous runtime-check image does not include
later source changes:

```sh
nextflow run pipeline/main.nf -profile docker \
  --mode downstream-only --project /data/study/project.json \
  --flow_config /data/study/flow.json \
  --katana_container "$KATANA_IMAGE" \
  --outdir /data/study/new-results -resume
```

Installed wheels contain the same workflow sources. Export once using
`katana-flow export-pipeline --outdir /data/tools/katana-pipeline`, then use
`nextflow run /data/tools/katana-pipeline/main.nf` with the same arguments.
An export never replaces an existing pipeline directory.

For a new bulk study:

First [build the separate DEXSeq tools image](../containers/dexseq-tools/README.md).
The upstream annotation/count steps use unchanged GPL tools that are not included
in the Katana downstream image. No public DEXSeq image is assumed to exist.
Set `DEXSEQ_IMAGE` to the built local image or verified site-provided image/SIF.

```sh
nextflow run pipeline/main.nf -profile docker \
  --mode complete --input /data/study/samples.csv \
  --dexseq_tools_container "$DEXSEQ_IMAGE" \
  --contrasts /data/study/contrasts.csv \
  --fasta /references/genome.fa --gtf /references/annotation.gtf \
  --project_template /data/study/project-template.json \
  --flow_config /data/study/flow.json \
  --katana_container "$KATANA_IMAGE" \
  --outdir /data/study/results -resume
```

`--mode upstream-only` also requires `--dexseq_tools_container` with a container
profile. It produces the verified native bundle and does not claim
to have completed Katana analysis or its report. The examples use a locally
built or loaded image; they do not refer to a public registry. For a cluster,
use `-profile slurm,apptainer` or `slurm,singularity` and set
`--katana_container` to the matching image/SIF accessible to compute nodes.
Slurm account/partition settings belong in site configuration. The executor
controls admission; there is no hidden four-task concurrency cap. Standalone
local execution without containers uses `-profile local` with the native
environment already available.
Do not use a local executor on a shared login node for heavy analyses.
The Docker profile runs tasks with the invoking Linux user's UID/GID so tools
such as STAR can create private directories that Nextflow can still collect.
If site configuration replaces `docker.runOptions`, retain
`--user $(id -u):$(id -g)` alongside any additional options. Apptainer and
Singularity retain their native caller-identity behavior.

### Resume an unchanged downstream plan

For a report-only presentation update, retain the original plan and cache and
add `--flow_report_revision <new-version-label>` to the resumed Nextflow command.
This changes only the report task identity, preserving earlier reports and
reusing completed calculations and collection. The label accepts letters,
digits, dots, underscores and hyphens. It does not change scientific parameters
or authorize changed inputs. The report task requires `flow.acceptance.json`.

RBP and boundary-motif enrichment are shown per comparison as compact
significant-result tables (native FDR ≤ 0.05). Every row exposes all recorded
hit genes, P/FDR and its native effect; repeated support plots and all tested
rows remain in outputs. Motif contributors are joined from completed matched
foreground pairs and recorded FIMO sites, never inferred from an RBP name.
MotifTrace sequence Top 10 and complete model/test files remain separate.

Keep the original Nextflow launch directory, `.nextflow` cache, work directory,
input files and original absolute `plan.json` path. Validate that plan and its
project before resuming:

```sh
katana-flow validate-plan \
  --plan /data/study/original-plan/plan.json \
  --project /data/study/project.json \
  --out /data/study/reuse-validation-01.json

nextflow run pipeline/main.nf -profile docker \
  --mode downstream-only --project /data/study/project.json \
  --flow_plan /data/study/original-plan/plan.json \
  --katana_container "$KATANA_IMAGE" \
  --outdir /data/study/new-results -resume
```

Use the exact plan path passed to the original tasks, not a new copy. The
validator checks plan identity, unchanged project bytes and bound resources;
`ready: true` means that reuse is ready, not that the analysis is complete.
`--flow_plan` bypasses planning and does not need `--flow_config`. If both are
provided, the workflow logs that the original plan is used and the config is
not reapplied. Changing comparisons, scientific parameters or resources requires
a new explicit plan. Existing successful work is reused only when Nextflow's
task fingerprint still matches; a new image or executable can invalidate that
fingerprint. Do not edit cache records to force a hit.

### Recover a completed MotifTrace calculation whose display failed

There is a narrow recovery route for the recorded error
`Portable event report exceeds 64MiB; complete tables retained, no truncation`.
It verifies the failed task, comparison, parameters, dependencies and full
retained statistical family, then generates its presentation in a new task
directory. It does not recompute statistics, change thresholds, truncate the
complete tables or rewrite the original failure as success.

Set `TASK_ID` to the exact failed task ID in the existing plan. Supply all of
that task's original dependency receipts, as recorded in its task command:

```sh
katana-flow task \
  --plan /data/study/original-plan/plan.json \
  --project /data/study/project.json --task-id "$TASK_ID" \
  --dependency-receipts /data/study/dependency-one/task.receipt.json \
                        /data/study/dependency-two/task.receipt.json \
  --presentation-recovery /data/study/failed-task/failure.receipt.json \
  --outdir /data/study/recovered-presentation-01
```

The dependency paths above illustrate two dependencies; use the complete list
for the actual task. The new result has its own receipt and is subject to normal
collection and report checks. This command does not inject a result into the
Nextflow cache. Other errors, incomplete numeric outputs, changed inputs or
already completed task families are rejected by this special recovery route.

## Calculation and result contracts

- Resource checks precede new upstream execution. A complete imported project
  is checked again before downstream tasks are scheduled.
- Only declared, selected contrast IDs are scheduled. The plan does not invent
  all pairwise comparisons or infer experimental purposes from sample names.
- Each task has its own manifest, output directory and catalog. Parallel jobs
  do not write concurrently to one DuckDB file.
- AS read support, landscape, DEG/DAS/DTU/DEXSeq integration, exact overlaps,
  four evidence-specific enrichments and RBP identity use the existing
  scientific workflows. Each comparison has independent MotifTrace and
  boundary/CLIP tasks. Backgrounds match detected/testable evidence.
- Gene-body QC is calculated from BAM and matched BED12; it is not inferred
  from FastQC. BED12 is generated from the matching GTF when no explicit BED
  is supplied. Coverage summaries wait for all declared sample measurements.
- For BAM input, an independent FastQC task measures each supplied BAM's
  retained records. The report labels this as **Supplied BAM records**;
  it cannot recover discarded reads or the original trimming/alignment QC.
  FASTQ inputs retain separate raw, trimmed and alignment measurements.
- By default each selected comparison gets its own network discovery,
  associations and annotation tasks, with explicit dependencies. Hub events
  are ranked by `|kME|`, retaining signed membership and complete event/gene
  tables. Bootstrap is disabled by default; explicitly setting
  `network.bootstrap.enabled: true` adds the real resampling stage. Disabled
  resampling is recorded as a user-selected omission, not a failed calculation
  or a completed confidence analysis.
  Its input includes only the comparison's biological samples. The network
  tables, figures and enrichment keep this identity throughout collection and
  reporting. Numerical linkage and leaf IDs are saved before
  plotting, using the unchanged fit. Edge-support bootstrap and module
  bootstrap remain distinct.
- Network configuration can be generated from the imported project's native
  rMATS sources and sample-vector identities. An explicit existing YAML may
  override configuration. JC/JCEC selection uses `network.parameters.rmats_mode`
  first, then `project.metadata.rmats_mode`, then the sole mode declared by the
  selected network's rMATS assets. Supplying both modes requires an explicit
  choice; a missing mode or absent requested source fails readiness. The resolved
  mode is recorded in new task plans. An execution-only compatibility wrapper restores
  integral ordinal-rank indices for dynamicTreeCut; ranks, ordering, fitting
  parameters and frozen source files remain unchanged.
- Interaction tasks require explicit designs and real counts. Simple contrast
  differences are not relabeled as fitted interactions.
- MAEI prepares exact annotated SE inclusion/skipping transcript TPM and the
  official 5,440 kmer counts in the existing 14-window family, with one
  independent task per comparison/window. It then executes the official
  external MAPP model, Z-score and native mixture
  significance chain. It requires hash-bound source scripts and matching
  reference annotation, genomic FASTA and existing transcript quantification.
  An explicit already-prepared quantity/count input remains supported. It
  must describe one selected comparison; use project-derived inputs to
  prepare independent families for multiple comparisons. The workflow
  does not run MAPP's other eight
  subworkflows, create TPM from rMATS PSI, or claim calibrated MotifTrace/MAPP
  parallel evidence from an uncalibrated native fit. New-project preparation,
  native fitting and report panels have component tests. Generated FASTQ/BAM entries have passed all 14-window families and the
  complete report handoff; real-study and target-site validation remain distinct.
- Completed task receipts bind parameters, method versions, input and output
  SHA256 and dependency identities. Unknown/failed outputs cannot be adopted
  as success. `-resume` reuses the exact Nextflow task cache; prior independent
  studies require an explicit verified adapter and are not automatically
  relabeled as this run.
- A failed downstream task makes the run fail. Its already submitted independent
  tasks finish normally so their successful outputs can be reused with
  `-resume`; the workflow does not automatically retry or ignore failures.
  A report depending on a failed required task cannot run as complete. If the
  final report gate fails after rendering, its ZIP is named
  `REPORT.diagnostic.zip`, and no successful `flow.acceptance.json` is issued.
- Collection waits for every required task. A partial/test collection cannot
  become a complete report. The report acceptance gate rejects missing fixed
  computations and missing required chapter figures. An explicitly declared
  single-contrast design without interactions retains chapter 4 with a
  not-applicable explanation instead of a fabricated overlap plot.

## Configuration

`project-template.json` supplies `metadata.species`, resource `assets`, optional
`sample_attributes`, and explicit JC/JCEC selection. Sample identities,
contrast direction and rMATS vector order come from the verified upstream
bundle. Resources retain their own species, build, biosample and source
provenance. Reference CLIP is not a measurement in the current experiment.

`flow.json` uses schema `katana.flow.config/v1`:

```json
{
  "schema": "katana.flow.config/v1",
  "contrast_ids": ["B_vs_A"],
  "event_set": "stable_significant",
  "parameters": {"differential.deg": {"fdr": 0.05, "min_abs_log2fc": 1.0}},
  "qc": {"reference_bed": "reference.bed", "rseqc_executable": "geneBody_coverage.py"},
  "network": {"scope": "contrast"},
  "design_factors": ["condition"],
  "interactions": [],
  "mapp": {"adapter": "mapp.maei", "inputs": {}, "parameters": {}}
}
```

With empty MAPP `inputs`, a new plan binds the package's unchanged bundled
author source and recorded hashes, then resolves those files inside each
worker. A host-specific `native_root` is not required. Preparation still needs
the project's matching transcript quantification, event identities, GTF and
genomic FASTA to produce included/total TPM and motif counts. Explicit external
author-source inputs remain supported with their exact source hashes; they
are not silently replaced. This example is a configuration shape, not a
resource-complete study. Missing reference or quantification inputs fail
readiness rather than producing an empty successful result.
Declared multifactor studies need explicit `interactions` entries with `id`,
`design_path`, `reference_a` and `reference_b`. Previously approved thresholds
remain in [THRESHOLDS.md](THRESHOLDS.md). Scientific applicability is recorded
separately from a failed or unavailable computation.

Inspect without executing:

```sh
katana-flow preflight --project project.json --config flow.json
katana-flow plan --project project.json --config flow.json --out plan.json
```

The fixed report stays read-only. `katana-flow report` is the strict completion
entrypoint; `katana-report` may still render honestly labeled historical partial
collections. Neither command changes source statistics or invents missing data.

Every planned native MAEI task is required by the completion gate. The approved
report also supports importing compatible, already calibrated MAPP/MotifTrace
parallel results; this is conditional on supplying that separate evidence.
Register a `motiftrace_parallel_results` asset, selecting its asset ID through
`parameters.prioritize.motiftrace_parallel.parallel_bundle_asset_id` if needed.
The planner then schedules an independent verified import for each selected
comparison. Missing, changed, ambiguous or mismatched supplied bundles fail
readiness or acceptance. Without such a bundle, the report states that the
optional calibrated comparison was not supplied; this does not suppress native
MAPP fitting or relabel native statistics as calibrated evidence. The historical
method-comparison calibration is not automatically rerun for a new study.
