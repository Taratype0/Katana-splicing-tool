# Choose analyses and manage execution

Katana runs the selected analyses and then builds their report and browsing
files. By default it selects the full seven-chapter report. To omit an analysis,
set its switch to `false` in your existing `flow.json`:

```json
{
  "design_factors": ["condition"],
  "analyses": {
    "network": false,
    "mapp": false,
    "coverage": false
  }
}
```

Merge these fields into your study's analysis configuration. A disabled analysis
has no calculation task, does not require its dedicated resources, and is
recorded as `not_selected` in the report manifest. Enabled analyses still need
their inputs and completed outputs. A failed selected calculation is not treated
as a user opt-out.

| Switch (default `true` except `coverage`) | What it includes |
|---|---|
| `landscape` | AS landscape, PSI summaries and associated report views |
| `cross_comparison` | Event overlap and direction comparisons, where at least two contrasts are selected |
| `enrichment` | GO/pathway analysis for the selected DEG, DAS, DTU and exon-use gene sets |
| `rbp_identity` | RBP identity-set enrichment and hit genes |
| `motiftrace` | Mature-RNA MotifTrace; significant motifs, sequence views and complete output tables |
| `motif_binding` | Boundary motif and CLIP analysis |
| `mapp` | Native MAPP MAEI fits and any explicitly supplied calibrated parallel results |
| `sequencing_qc` | Imported sequencing QC summaries (FASTQ/BAM caller outputs) |
| `coverage` (default `false`) | Optional extra RSeQC gene-body BAM scan and summary; set `true` explicitly to calculate |
| `network` | Independent per-comparison networks, associations, annotation and browsers |
| `interaction` | Explicitly configured factorial interaction models |

The study design, input audit, read-support classification and core molecular
results remain the report foundation. Interaction models still require an
explicit design; the switch does not invent one. The default approved chapter
order, colours, comparison selector and figure formats are unchanged. A wholly
disabled optional chapter is omitted without renumbering the remaining ones.
This selection applies to downstream calculations; it does not disable the
upstream expression/splicing calls needed to produce the core input bundle.

The extra RSeQC `geneBody_coverage.py` scan is disabled for new plans unless
`analyses.coverage` is explicitly `true`. The default RSeQC module lists in
[nf-core/rnaseq 3.26.0](https://github.com/nf-core/rnaseq/blob/3.26.0/nextflow.config)
and [3.27.0](https://github.com/nf-core/rnaseq/blob/3.27.0/nextflow.config) do not
include this tool. Routine BAM processing, strand inference and sequencing QC
are separate and remain enabled. Skipping this optional measurement does not
change the seven-chapter layout or claim that a coverage curve was computed.
Existing explicit selections and sealed plans are unchanged. To stop this extra
scan in a new plan based on an older config, set `analyses.coverage` to `false`.

## Resources and partitions

The downstream process defaults follow the task-size pattern used by
[nf-core/rnaseq](https://github.com/nf-core/rnaseq/blob/master/conf/base.config)
and [nf-core/rnasplice](https://github.com/nf-core/rnasplice/blob/master/conf/base.config):

| Task | CPUs | Memory | Scheduler time |
|---|---:|---:|---:|
| Network fit | 12 | 72 GB | 16 h |
| MotifTrace, MAPP preparation and MAPP fit | 6 | 36 GB | 8 h |
| RSeQC gene-body coverage (one BAM per task) | 1 | 4 GB | 12 h |
| Sequencing-QC collection, reference BED and project audit | 1 | 4 GB | 2 h |
| Other downstream analysis tasks | 2 | 12 GB | 4 h |

These are starting requests, not dataset-size estimates or guarantees of fit.
Selected resource/signal terminations receive one retry, with the default
memory and time requests doubled. CPU count stays fixed: more threads do not
help a serial scan. Other errors finish with failure so they can be diagnosed.
Optional native `max_cpus`, `max_memory` and `max_time` settings bound requests.
No resource value causes a switch to a `memory` partition. The user/site selects
the queue; Slurm controls admission and actual concurrency.

Use an ordinary Nextflow site configuration for executor/container/account
setup. The managed command also accepts a separate execution-only JSON:

```sh
katana-flow run --mode downstream-only \
  --project imported-project/project.json --flow_config inputs/flow.json \
  --profile slurm,apptainer --site-config site.config \
  --resource-config resources.json --outdir results
```

Copy [resources.json](../examples/nextflow/resources.json) and select appropriate
task requests. The generic example does not choose a queue or account.
`defaults` optionally sets queue/account; `labels` overrides a Nextflow process
label; `processes` overrides exact uppercase Nextflow aliases, with higher
priority than label/default settings. Explicit resource values
are fixed overrides; default automatic scaling does not change those explicit
values. An optional `maxRetries` value changes the retry count for that alias.
The `account` field generates Slurm's `--account`; configure other schedulers
with their native site configuration instead.

## Configure your own cluster once

The pipeline declares CPU, memory and wall-time requests. The site configuration
chooses an executor and containers; the resource file chooses queues if needed.
There are no built-in `main`, `memory` or `shared` queue names, automatic queue
switches, exclusive-node requests or manual dependency-submission loops.
Nextflow submits ready tasks and monitors completion; the scheduler allocates
resources. Keep `--slurm_partition` unset when using per-task queue mappings.
A global partition selection intentionally applies to everything without a more
specific override.

For Dardel, copy [resources.dardel.json](../examples/nextflow/resources.dardel.json),
replace `YOUR_ACCOUNT`, and use it with `--profile slurm,apptainer`. It requests
`shared` by default, including RSeQC and light processing. It does not force even
threaded tasks onto whole nodes. If a particular task requires a different
partition under your allocation, add that exact alias to `processes`. Confirm
partition limits with the site; no universal CPU/memory threshold identifies a
shared partition across clusters.

For another Slurm cluster, replace the queue and account names. For local use,
use `--profile local,docker` and the queue-free generic example. PBS/SGE/cloud
users can set `process.executor` and scheduler-specific options in a standard
Nextflow `--site-config`; the JSON `account` convenience is Slurm-specific.
Do not copy Slurm account directives into other executors.

For example, a site with queues named `short` and `compute` can use:

```json
{
  "defaults": {"queue": "short"},
  "labels": {"process_high": {"queue": "compute"}},
  "processes": {
    "QC": {"queue": "short", "time": "24h"},
    "NETWORK_FIT": {"queue": "compute", "memory": "96 GB"}
  }
}
```

These names are examples, not built-in partitions. `QC` includes per-sample
RSeQC coverage and sequencing-summary collection; their native time defaults
are different. A `QC` override deliberately applies to both. For finer control,
use a dynamic selector in your site config keyed by `spec.workflow_id`.

Upstream DEXSeq counting/annotation, FASTQ staging, reference preparation and
FASTA indexing also request one CPU for their current serial commands. Their
existing memory/time defaults remain intact; parallel aligners, quantifiers and
model fits retain their tool-specific allocations.

RSeQC `geneBody_coverage.py` scans alignments for each transcript and has no
thread-count argument in this adapter. Its initial request is now one CPU and
4 GB, with 12 hours for the potentially long scan. This allocation reduces wasted
reservation; it does not make a serial algorithm faster. Different samples are
independent tasks. Large BAMs or annotations may need more memory/time, configured
without changing the measurement or forcing a different queue.

Review `pipeline_info/trace-*.tsv`, execution reports and your scheduler's usage
statistics before adjusting resources. CPU budgets are passed to the Katana
worker/BLAS environment; a budget is a ceiling, not a promise of full utilization
in every phase. Defaults for network/MotifTrace/MAPP remain initial requests,
not a claim of measured efficiency on every dataset.

Reuse completed QC through existing task receipts/cache when inputs, annotation,
method/version and parameters match. A historical Qualimap coverage result does
not automatically replace an RSeQC gene-body profile: their measurements differ.
This resource update neither reruns nor relabels existing study results.

References: [Nextflow resource allocation and retry](https://training.nextflow.io/latest/config_exec/02_resources_and_retries/),
[Nextflow site configuration](https://www.nextflow.io/docs/latest/config.html),
and [RSeQC gene-body coverage](https://rseqc.sourceforge.net/#genebody-coverage-py).

## Resume without repeating completed calculations

If a run fails, keep its run directory and use the same command with `--resume`.
You may change `resources.json` first. Each resource revision is recorded
separately from the scientific settings; changing only that file no longer
invalidates the managed study identity.

Nextflow performs dependency scheduling and cache reuse. There is no manual
submit/check/submit loop between comparisons. A change to task inputs, code,
container or the generated command can require recalculation. For example,
changing CPU count changes a task hash when that count appears in the command.
Memory/time-only adjustments normally preserve otherwise identical task hashes.
See [Nextflow cache and resume](https://www.nextflow.io/docs/latest/cache-and-resume.html).

Keep scientific settings, inputs and arbitrary Groovy site configuration
separate from execution-only JSON. Site configuration can contain scientific
overrides, so changes to that file still require a new managed run. Editing
scientific parameters is not a resource-only resume.

Completed upstream handoffs can release upstream work before downstream runs.
After successful final delivery and reload checks, managed runs release their
registered work. Failed/unknown work stays available for recovery. `--keep-work`
retains native work after success. Input files and shared tool/container caches
are not part of this cleanup.

## Network and MAPP execution

Networks still discover modules using the selected correlation, TOM and module
algorithm. kME ranks event membership/Hubs after discovery; it does not replace
module discovery. Bootstrap remains off by default. The legacy
`tom_dense_max_events`, `consensus.dense_max_events`,
`preservation.dense_max_events` and `wgcna_timeout_seconds` protection fields
are accepted for compatibility but do not reject/cancel a calculation. There is
no internal event-count stop or R-WGCNA elapsed-time deadline. User-requested
event selection and scientific filters remain explicit. Scheduler wall time and
available physical memory still apply.

MAPP now reads the comparison's annotation/TPM and builds transcript mappings
once. Its 14 site/window preparations share that context; the 14 native model
fits remain separate tasks. Windows, model families, TPM eligibility and
statistics are unchanged. Preparation is still independent across comparisons.

## File verification and retention

New task receipts store file size and filesystem change identity alongside the
original SHA-256. Unchanged task inputs/products use a metadata check; changed
or transported files are hashed. Older receipts without that identity fall back
to content verification. SHA-256 calculations are also reused within one
process for unchanged files, including hard-linked aliases. Set
`KATANA_VERIFY_CONTENT=1` to request content verification at receipt boundaries.

This fast path assumes the filesystem reports changes correctly, as native
Nextflow caching does. For archives or filesystems with unreliable timestamps,
use content verification. Final delivery still checks content and reloads the
retained project before deleting work.

Completed, owned task products are retained with hard links on the same
filesystem and copied across filesystems. Hard links survive removal of work;
they are not symlinks into deleted directories. These products are immutable:
make a separate copy before editing their content. Borrowed study inputs are
not converted into writable shared outputs. Supporting outputs now go directly
to their final staging tree instead of being copied through another full tree.
