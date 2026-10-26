# Storage management

`katana-flow run` runs the native Nextflow workflow with verified storage
handoffs. It removes the completed upstream work before downstream analysis,
then removes downstream work after the final report and retained analysis
project have passed verification. The scientific modules, selected comparisons
and seven-section report format are unchanged.

Use a matching updated Katana installation and workflow. Existing wheels,
exported workflows, container images and running jobs are not updated by a
change to GitHub `main`.

## Run with automatic work cleanup

Run on a Linux execution host with native Nextflow, compatible Java and the
chosen tool environments or containers. Windows users can run these commands
inside an already configured WSL/Linux environment; native Windows Nextflow
execution is not claimed.

Run from a separate study directory, outside the software checkout. Managed
storage deliberately refuses to create reclaimable work inside the installed
source tree.

For FASTQ or BAM, use the same samplesheet, contrasts, references, resources and
scientific settings as the [Nextflow usage guide](NEXTFLOW_USAGE.md):

```sh
katana-flow run --mode complete \
  --input inputs/samples.csv --contrasts inputs/contrasts.csv \
  --fasta references/genome.fa --gtf references/annotation.gtf \
  --project_template inputs/project-template.json \
  --flow_config inputs/flow.json \
  --outdir results \
  --profile slurm,apptainer --site-config site.config \
  --katana_container "$KATANA_IMAGE" \
  --dexseq_tools_container "$DEXSEQ_IMAGE"
```

Use `--profile local` for an installed native local environment, or a matching
container profile such as `docker`. Account, partition, CPU, memory and wall
time remain site/executor settings; cleanup does not change the scheduler's
admission policy or scientific parameters. Repeat `--site-config` to supply
more than one native Nextflow configuration file.

An existing prepared project can enter downstream analysis directly:

```sh
katana-flow run --mode downstream-only \
  --project inputs/project.json --flow_config inputs/flow.json \
  --outdir downstream-results \
  --profile slurm,apptainer --site-config site.config \
  --katana_container "$KATANA_IMAGE"
```

Use `--flow_plan /absolute/original/plan.json` instead of `--flow_config` to
reuse an exact compatible sealed downstream plan. This is not an automatic
adoption of arbitrary older Nextflow work directories. Completed results from
another run still need their supported import/adoption or project mapping.

`--outdir` is required. The default private run directory is its sibling
`results.run`; choose another with `--run-dir`. Keep input project/configuration
directories separate from both output and run directories. New runs do not
overwrite nonempty completed deliveries.

Use `--params-file inputs/scientific-params.yaml` for supported native
Nextflow scientific/resource parameters. Keep managed routing fields, output
locations and cleanup policy in the `katana-flow run` arguments. A scientific
params file for this command is not the complete mode/input/output params file
shown for a direct `nextflow run` invocation.

## What is deleted, and when

| Point in the run | Required evidence | Cleanup |
| --- | --- | --- |
| Upstream handoff | Native upstream execution succeeded; the complete upstream bundle and every available member pass size/path/SHA checks. | Delete only this managed run's upstream work. The verified bundle remains available to downstream analysis. |
| Final downstream handoff | Native downstream execution succeeded; all planned calculations, the fixed report, final delivery and retained project pass their verification gates. | Delete only this managed run's downstream work, then verify the final delivery again. |
| Failure or interrupted execution | Success or engine termination is absent or unverified. | Preserve the affected work and diagnostic records. No success cleanup is inferred from an empty queue or a partially written report. |

These are whole-stage handoffs. A process finishing does not authorize deletion
of its inputs: other comparisons, QC tasks or downstream consumers may still
need them. Katana does not delete a task's files while its declared consumers
are using them.

Cleaned upstream intermediates include trimmed-read copies, unsorted and
name-grouped BAMs, task-local temporary files and regenerable indexes. The
retained project keeps the formal input evidence needed for subsequent target
analysis, including available BAM/index and reference resources. Original raw
data, source projects and shared reference/database locations are not deleted.

## What remains

```text
results/
  project/
    project.json              # browser/API/target-analysis entry
    inputs/                   # BAM/indexes, references, resources and input tables
    collection/               # complete scientific results and original proofs
    retention.receipt.json
  report/
    REPORT.html
    REPORT.zip
    report-manifest.json
    flow.acceptance.json
  upstream-support/           # complete-mode supporting native outputs
  pipeline_info/
  runner-logs/
  delivery.receipt.json
  managed-run.receipt.json

results.run/
  managed-run.json
  logs/
  cleanup-*.receipt.json

.katana-runtime-cache/          # shared tool-cache default when no site setting exists
```

The retained project has its own verified resource bindings. It must not depend
on soft links into deleted work. Available original producer products and
proofs are preserved; relocating resource/catalogue paths does not refit the
statistics. The report is copied without redesigning or rewriting its format.

`upstream-support` preserves native results not included in the canonical
bundle: fitted RDS objects, DEXSeq diagnostic PDFs, SUPPA definitions/statistics,
Salmon `aux_info`, library-format statistics, command provenance, QC and logs.
They are scientific outputs, not disposable cache. A Salmon directory can also
retain its small `quant.sf` table because its declared output is the complete
native directory. Already bundled large BAM/reference copies are not separately
published again in lean mode. Rebuildable STAR indexes need not be retained as
final results; user-supplied shared indexes remain untouched.

This controls temporary storage and unnecessary duplication. Completed
scientific tables, native models and result matrices remain available even
when they are large.

On a filesystem supporting hard links, identical retained products can share
one owned data allocation. The manager can also hand off immutable files from
its locked, registered downstream work without making another full BAM copy,
but only after proving their existing links belong to that work. Borrowed user
inputs are copied independently; unsupported or cross-filesystem linking falls
back to verified copies. File names and complete scientific evidence remain
available in either case.

## Keep work or clean a completed run later

Add `--keep-work` to the original command to retain native stage work after
success. Lean publication still applies, but automatic deletion is disabled:

```sh
katana-flow run --mode downstream-only \
  --project inputs/project.json --flow_config inputs/flow.json \
  --outdir downstream-results --keep-work --profile local
```

Preview later cleanup of that completed managed run:

```sh
katana-flow clean-run --run-dir downstream-results.run
```

Apply the verified, bounded plan:

```sh
katana-flow clean-run --run-dir downstream-results.run --apply
```

The default is a dry run. `clean-run` accepts a completed managed scope with
its ownership records, native completion evidence and verified delivery. It
does not accept an arbitrary `work` directory, running/failed scope or guessed
filesystem target. Cleanup receipts record eligible paths, logical bytes and
estimated reclaimable allocation. They do not claim a measured return of
filesystem/project quota; retained hard links are not counted as freed space.

## Resume and failures

Repeat the same managed command with `--resume` and the same `--run-dir` after
resolving an observed failure. Inputs, software, parameters, configuration
identities and the original cleanup choice must still match. Changed settings
require a new run directory rather than editing an old completion record.

A verified completed handoff is reused. If the affected native work remains,
Nextflow can reuse eligible matching completed tasks with its normal `-resume`
rules. Work deleted after a successful handoff no longer supports that stage's
native cache resume; the independent checkpoint is retained instead. Neither
mode promises that unrelated historical results are automatically cache hits.

An interrupted engine whose termination is unknown is not automatically
submitted again. Check its native processes/jobs and recorded attempt before
resuming; preserve the original attempt and its logs. The manager does not
silently cancel or replace unknown-running scientific work.

## Tool caches and site configuration

Container and Conda caches are separate from scientific task work. When no
external setting is provided, managed execution puts the applicable tool
caches under `.katana-runtime-cache` beside the private run directory, allowing
neighbouring runs to reuse them. Existing `NXF_*_CACHEDIR` environment
settings and explicit site cache configuration remain authoritative.

These tool caches are not automatically removed: another run or user may need
them. `clean-run` also does not remove original inputs, shared references,
databases, installed environments or caches outside its registered stage work.
Configure tool caches outside reclaimable work before enabling cleanup.

For a direct native invocation, `--storage_mode lean` reduces duplicate
publication only; it does not perform managed stage handoff or delete work.
Direct Nextflow remains available with its ordinary work/resume behaviour.
`--storage_mode retain` preserves the legacy publication layout.

## Retain an already completed collection

An existing verified collection can be made independent of its source work
without rerunning completed analysis:

```sh
katana-flow retain-project \
  --collection completed/collection --project inputs/project.json \
  --outdir retained-project

katana-flow verify-retained --directory retained-project
```

Keep the source project and resources available until retention succeeds.
`verify-retained` verifies resource/product bytes and reloads the project at
its present location. Retention does not by itself authorize deletion of an
unmanaged old run directory.

Open the retained `project.json` through the Katana browser or use the same
Python/CLI project entry. Saved results and gene/event selections can feed the
[target-analysis workflow](TARGET_COMMANDS.md). New Sashimi still needs the
retained usable BAM/index resources; a static HTML report cannot execute a
target job on its own.

## Validation boundary

Storage tests use generated software fixtures to exercise ownership, copying,
publication, resume guards and cleanup after verified handoff. They do not
rerun an existing study or establish fresh whole-study scientific validation
for every platform or scheduler. Matching installation/container checks and
runtime qualifications remain part of the release's
[acceptance record](validation/entrypoints-runtime-20261025.md).
