# Task resource allocation

Katana uses a version-pinned **per-tool** mapping from nf-core/rnaseq 3.27.0,
nf-core/rnasplice 2.0.0 and nf-core/differentialabundance 2.0.0. Copying the generic
`process_high` label alone is not equivalent to inheriting a pipeline's effective
settings: rnaseq has a separate Trim Galore override.

This is a resource-only integration. The existing tool versions, model formulas,
paired design, analysis arguments and approved report format are unchanged.
These defaults are requests, not measurements of every study's actual CPU/RAM
usage. A complete scientific run was not restarted to collect resource metrics.

## Initial upstream requests

| Task | CPU | Memory | Time | Pinned source or exception |
|---|---:|---:|---:|---|
| Trim Galore | 8 | 4 GB | 8 h | rnaseq 3.27.0 explicit CPU/memory override and medium time |
| FastQC, including BAM alias | 2 | 12 GB | 4 h | rnaseq 3.27.0 / rnasplice 2.0.0 low label |
| DEXSeq annotation / per-sample counting | 1 | 6 GB | 4 h | rnasplice 2.0.0 single label |
| DEXSeq exon / DTU models | 6 | 36 GB | 8 h | rnasplice 2.0.0 medium label |
| DESeq2 contrast / normalization | 1 | 6 GB | 4 h | differentialabundance single label; normalization is serial |
| FASTA index / lane staging | 1 | 6 GB | 4 h | rnaseq single label; serial Katana lane adapter |
| tximport | 6 | 36 GB | 8 h | rnasplice medium request retained |
| Salmon projection / samtools sort | 6 | 36 GB | 8 h | Corresponding rnaseq medium requests retained |
| STAR alignment | 12 | 72 GB | 16 h | rnaseq high label retained |
| STAR genome index | 16 | 96 GB | 16 h | Existing GRCm39 measured-memory exception; 128 GB on retry |
| rMATS prep / post | 12 | 72 GB | 16 h | Current contrast-level implementation; see below |

The [complete inventory](resource-allocation-audit.csv) covers all 60 declared
processes, including optional branches. Each row records its upstream source or
states that the task is Katana-specific. The machine-readable
[policy](../pipeline/assets/resource-policy.json) pins repository commits and
source hashes. Runtime aliases and downstream execution groups are additionally
covered by the native allocation probe.

### Preserve the unit of work when borrowing a request

Katana's current rMATS prep still processes both groups for one contrast, as in
rnasplice 1.0.4. rnasplice 2.0.0 uses `process_single` for **one BAM per prep task**.
Applying the newer single-BAM reservation to the older multi-BAM task would not
be an equivalent resource mapping. The contrast-level allocation is retained
explicitly; this patch does not silently change preparation granularity or
rerun completed results.

The existing STAR-index request is also explicit: recorded GRCm39 index peak
memory was about 68.3 GB, so the already approved 96 GB first attempt is retained.
It is not presented as a universal requirement for every reference genome.

## Katana downstream

| Work | Initial CPU / memory / time | Basis |
|---|---|---|
| CoSplice network fit | 12 / 72 GB / 16 h | Existing Katana policy; threads passed to BLAS/OpenMP |
| MotifTrace | 6 / 36 GB / 8 h | Existing Katana policy; not an nf-core upstream module |
| MAPP preparation / fitting | 6 / 36 GB / 8 h | Existing Katana policy; separate tasks |
| Network follow-up | 2 / 12 GB / 4 h | Reuses fitted network |
| Optional gene-body coverage | 1 / 4 GB / 12 h | Serial traversal |
| Sequencing-QC summary / reference BED / audit | 1 / 4 GB / 2 h | Existing serial policy |
| Other downstream summaries | 2 / 12 GB / 4 h | Existing Katana starting request |
| Report | 2 / 12 GB / 8 h | Rendering and packaging existing products |

These modules have no interchangeable resource configuration in the three
upstream pipelines. Their allocations are identified separately rather than
claimed as inherited nf-core defaults. Memory needs vary with event count,
model count, annotation and cohort size. Downstream retries increase memory and
time without automatically increasing CPU count or changing partition.

## Site and user configuration

Pass site/user settings with normal Nextflow `-c` files. Katana does not choose
a cluster partition or edit a user's configuration. Exact `withName` selectors
can override CPU, memory, time and queue; fully qualified names can distinguish
a specific invocation. Generic label overrides have lower priority than
process-name overrides. No Dardel-specific memory-per-CPU rule is built in.

The native print-only probe loads the production configuration and checks all
60 process declarations, the BAM FastQC alias, downstream execution groups,
retry requests and independent synthetic user overrides. It simulates local
executor capacity because it prints requests without allocating those resources.
It is not a Slurm benchmark or a real-data performance measurement.

For real runs, retain Nextflow trace/report and compare Slurm requested vs
allocated CPU/memory and peak RSS. A scheduler may allocate extra CPUs to satisfy
its memory-per-CPU policy; pipeline worker threads and charged allocation can
differ. `Priority` is a queue reason, not proof of quota exhaustion.

## Sources

- [rnaseq 3.27.0 Trim Galore override](https://github.com/nf-core/rnaseq/blob/a1fcdddd3b826fe46eb46f0479f2ff8a7815af05/conf/modules/trimgalore.config)
- [rnasplice 2.0.0 modules](https://github.com/nf-core/rnasplice/tree/aceb47ab77cf0edfb290a44035328ed5b7eb0be8/modules)
- [differentialabundance 2.0.0 DESeq2](https://github.com/nf-core/differentialabundance/blob/30ed7741fc392127156c2fb10cfa3d69d216b54b/modules/nf-core/deseq2/differential/main.nf)
- [Nextflow selector precedence](https://www.nextflow.io/docs/latest/config.html#process-selectors)
