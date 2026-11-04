# Task resource allocation

Katana uses a version-pinned **per-tool** mapping from nf-core/rnaseq 3.27.0,
nf-core/rnasplice 2.0.0 and nf-core/differentialabundance 2.0.0. Copying the generic
`process_high` label alone is not equivalent to inheriting a pipeline's effective
settings: rnaseq has a separate Trim Galore override.

This integrates upstream resource requests and per-BAM rMATS preparation. The existing tool versions, model formulas,
paired design, analysis arguments and approved report format are unchanged.
These defaults are requests, not measurements of every study's actual CPU/RAM
usage. A complete scientific run was not restarted to collect resource metrics.

## Initial upstream requests

| Task | CPU | Memory | Time | Pinned source or exception |
|---|---:|---:|---:|---|
| Trim Galore | 8 | 8 GB | 8 h | rnaseq workers; measured 4.1 GiB peak requires more than the upstream 4 GiB starting request |
| FastQC, including BAM alias | 2 | 4 GB | 4 h | Measured 1.1 GiB peak in 12 mouse libraries, with headroom |
| DEXSeq annotation / per-sample counting | 1 | 6 GB | 4 h | rnasplice 2.0.0 single label |
| DEXSeq exon / DTU models | 6 | 16 / 8 GB | 8 h | Measured overrides of the rnasplice medium label |
| DESeq2 contrast / normalization | 1 | 6 GB | 4 h | differentialabundance single label; normalization is serial |
| FASTA index / lane staging | 1 | 6 GB | 4 h | rnaseq single label; serial Katana lane adapter |
| tximport / DRIMSeq filtering / stageR | 1 | 4 GB | 2 h | Completed native peaks up to 1.02 GiB; explicit per-tool overrides |
| SUPPA differential tests / event generation / PSI | 1 | 4 GB | 2 h | Completed native peaks below 0.85 GiB; no inherited 72 GiB request |
| SUPPA table splitting / BAM header audit | 1 | 2 GB | 2 h | Small completed table/header operations |
| Reference/metadata streaming | 1 | 4 GB | inherited | Observed reference tasks below 0.8 GiB; no extra full input copies |
| BAM evidence streaming | 2 | 6 GB | 8 h | Samtools decompression plus a serial AWK accumulator |
| Salmon projection | 6 | 12 GB | 4 h | Completed native peak 3.70 GiB |
| BAM collation | 6 | 8 GB | 4 h | Completed native peak 1.55 GiB |
| samtools sort | 6 | 36 GB | 8 h | Buffer budget derives from task.memory; larger user-buffer peaks are not transferred to this default |
| STAR alignment | 12 | 72 GB | 16 h | rnaseq high label retained |
| STAR genome index | 16 | 96 GB | 16 h | Existing GRCm39 measured-memory exception; 128 GB on retry |
| rMATS prep, per BAM | 1 | 6 GB | 4 h | rnasplice 2.0.0 single label and per-BAM granularity |
| rMATS post, unpaired contrast | 12 | 8 GB | 4 h | Completed post-only peak 1.40 GiB |
| rMATS post, paired contrast | 12 | 16 GB | 16 h | Engineering starting request, not a measured paired peak; existing OOM retry requests 32 GB / 32 h, with unchanged workers and model |

DEXSeq exon models now start at 16 GiB (observed 6.59 GiB), and DTU models at
8 GiB (observed 2.28 GiB), both retaining six workers and the eight-hour limit.
These named overrides replace the corresponding 36 GiB medium allocation.

The [measurement record](../pipeline/assets/measured-resources.json) combines
completed Slurm job-step MaxRSS with available Nextflow peak RSS and excludes
unfinished, failed and tiny public validation tasks from calibration. Batch
sampling can miss short application peaks; a nearly empty shell sample is not
evidence that an R model needs only a few megabytes. Counts are distinct completed
jobs, not independent biological experiments. No completed analysis was rerun.
Most older task directories have been cleaned; retained accounting observations
remain available, but do not reconstruct their deleted inputs or claim new
same-version scientific acceptance.

The [complete inventory](resource-allocation-audit.csv) covers all 61 declared
processes, including optional and retained historical declarations. Each row records its upstream source or
states that the task is Katana-specific. The machine-readable
[policy](../pipeline/assets/resource-policy.json) pins repository commits and
source hashes. Runtime aliases and downstream execution groups are additionally
covered by the native allocation probe.

### Preserve the unit of work when borrowing a request

Katana follows rnasplice 2.0.0: **one BAM per prep task**, reused by every contrast
that includes that sample. Post remains one task per contrast, with explicit
subject order preserved for paired models. The existing rMATS 4.1.2 / PAIRADISE
runtime and scientific arguments are retained. Historical group-prep modules
remain in the source snapshot but are not called by the current workflow.

DEXSeq annotation and counting use the published Katana DEXSeq tools image
(HTSeq 2.0.2 and the unchanged upstream GPL scripts). Nextflow downloads its
pinned Docker image or published SIF automatically. No local build is needed.
The original scripts are staged explicitly from `pipeline/vendor/dexseq`;
statistical fits still use the existing Bioconductor DEXSeq container.

`KATANA_FLOW_CONFIGURE` runs on the Nextflow launch host, including database
download and cache lookup. Only this small configuration task uses the local
executor; alignment, counting, models and other scientific tasks retain the
site executor. Generic `process.executor = 'slurm'` settings do not move the
download to network-isolated compute nodes. Launch where containers can run
and outbound network is available, or supply a pre-populated resource cache.

The existing STAR-index request is also explicit: recorded GRCm39 index peak
memory was about 68.3 GB, so the already approved 96 GB first attempt is retained.
It is not presented as a universal requirement for every reference genome.

## Katana downstream

| Work | Initial CPU / memory / time | Basis |
|---|---|---|
| CoSplice network fit | 12 / 72 GB / 16 h | Existing Katana policy; threads passed to BLAS/OpenMP |
| MotifTrace | 6 / 36 GB / 8 h | Existing Katana policy; not an nf-core upstream module |
| MAPP preparation | 6 / 12 GB / 8 h | Existing native preparation peak 4.45 GiB; preparation remains separate |
| MAPP fitting | 6 / 8 GB / 8 h | Existing native fitting peak 1.89 GiB |
| Network follow-up | 2 / 12 GB / 4 h | Reuses fitted network |
| Optional gene-body coverage | 1 / 4 GB / 12 h | Serial traversal |
| Sequencing-QC summary / reference BED / audit | 1 / 4 GB / 2 h | Existing serial policy |
| Other downstream summaries | 2 / 12 GB / 4 h | Existing Katana starting request |
| Report | 2 / 8 GB / 8 h | Existing native report peak 1.44 GiB; unchanged rendering and packaging |

These modules have no interchangeable resource configuration in the three
upstream pipelines. Their allocations are identified separately rather than
claimed as inherited nf-core defaults. Memory needs vary with event count,
model count, annotation and cohort size. Downstream automatic retries remain
disabled by default. If enabled explicitly, retry requests increase memory and
time without automatically increasing CPU count. The selected site profile can
choose a different partition for the larger request; scientific task arguments stay unchanged.

The current whole-study CoSplice and newer MotifTrace engine have different
execution scopes from the historical jobs. The old network reached 174.5 GiB;
its numbers cannot justify reducing the new engine's 72 GiB default. MotifTrace
likewise retains 36 GiB until comparable current-engine measurements exist.
Optional SUPPA clustering (`CLUSTEREVENTS`, disabled by default) has no matching
completed memory record and retains 72 GiB. These are explicit unresolved
calibration cases, not claims that every module has been measured or minimized.
Unused legacy group-prep/single-post declarations do not reserve any resources
in the current workflow.

## Site and user configuration

The generic profiles leave partition selection to the site. For Dardel select
the packaged `pdc_kth` profile:

```bash
nextflow -bg run Taratype0/Katana-splicing-tool -r v2.0.0-beta.5.dev2026110402 \
  -profile singularity,pdc_kth -c custom.config --project YOUR_ALLOCATION \
  ...
```

This replaces the separate `-c pdc_kth.config`; do not load both site policies.
Your existing `custom.config` stays unchanged. Named CPU, RAM, time and queue
overrides remain supported. `--slurm_partition` explicitly selects a partition;
an explicit `--partition` in `--clusterOptions` is preserved without a second
partition directive. In imported-project mode use `--slurm_account` for the
allocation, since `--project` identifies the project file.

The profile accounts for the observed Dardel shared memory policy:
`shared-equivalent CPUs = max(tool CPUs, ceil(memory MiB / 888))`.
It prefers shared up to 64 such CPUs and 111 GiB. This 64-CPU preference is
configurable, not a scheduler limit or a limit on the number of concurrent tasks.
It leaves moderate jobs on shared while moving larger reservations to main.
For example, STAR's 12-worker/72-GiB request needs 84 shared CPUs and is routed
to main. CoSplice has the same initial reservation. FastQC,
trimming, per-BAM rMATS prep, paired and unpaired rMATS post, SUPPA differential tests,
DEXSeq models, 36-GiB MotifTrace and the reduced MAPP requests
remain eligible for shared. A user-specified 256-GiB STAR index also goes to main.

Main is a whole-node allocation, limited here to 24 h/894 GiB. Longer jobs use
long when within 7 days/222 GiB. Larger or incompatible requests need an explicit
site override; Katana never selects memory automatically. Shared routing can be
adjusted with `--pdc_shared_max_cpus` (64) and
`--pdc_shared_memory_per_cpu_mb` (888 MiB). Other servers do not load this policy.
Partition choice cannot guarantee a shorter queue; admission remains Slurm's decision.

CPU worker counts remain fixed during memory retries, including named FastQC,
trimming and DEXSeq overrides. Existing retry eligibility/counts and scientific
commands are unchanged. Trace output records requested CPU/RAM/time, partition,
attempt, timestamps and peak RSS so future adjustments can use measurements.
No live jobs are cancelled or resubmitted by a profile update.

The native print-only probe checks all 61 declared processes and every wired
alias, including downstream execution groups and synthetic user overrides.
A separate test uses Nextflow's real Slurm executor with fake scheduler commands
to inspect the generated `#SBATCH` headers, routing boundaries and explicit
partition overrides. These checks allocate no HPC resources and run no scientific
data. Retained-cache and installed-package checks run separately in CI.

Memory requests for noncomparable or unmeasured heavy tasks are explicitly
identified above. The recorded measurements come from existing mouse studies, not a claim
that every species, depth or future dataset needs the same memory. Compare the
trace with Slurm allocation and peak RSS for each new workload. `Priority` means
waiting for scheduling priority; it is not proof of quota exhaustion.

## Sources

- [rnaseq 3.27.0 Trim Galore override](https://github.com/nf-core/rnaseq/blob/a1fcdddd3b826fe46eb46f0479f2ff8a7815af05/conf/modules/trimgalore.config)
- [rnasplice 2.0.0 modules](https://github.com/nf-core/rnasplice/tree/aceb47ab77cf0edfb290a44035328ed5b7eb0be8/modules)
- [differentialabundance 2.0.0 DESeq2](https://github.com/nf-core/differentialabundance/blob/30ed7741fc392127156c2fb10cfa3d69d216b54b/modules/nf-core/deseq2/differential/main.nf)
- [Nextflow selector precedence](https://www.nextflow.io/docs/latest/config.html#process-selectors)

- [PDC job scheduling and node capacities](https://pdc-web.eecs.kth.se/doc/run_jobs/job_scheduling/)
