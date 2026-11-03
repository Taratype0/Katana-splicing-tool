# Task resource allocation

The upstream module scripts and nf-core base resource labels are retained.
Katana does not choose a cluster partition, increase an account allocation, or
rewrite a user's configuration. Pass site and user configuration with `-c`.
A user-supplied exact `withName` selector overrides the corresponding pipeline
request; a broad label selector has lower priority than an exact process selector.
The resource-policy integration test verifies queue, CPU, memory and time overrides.

The [complete task inventory](resource-allocation-audit.csv) lists every declared
module and its initial allocation **before site/user overrides**. It is a code
inventory, not a benchmark or proof of peak memory on a particular study.
Aliases retain their original module labels. Downstream allocations are selected
by workflow/execution group in `pipeline/conf/task_resources.config`.

| Inherited label | Initial CPUs | Memory | Time |
|---|---:|---:|---:|
| process_single | 1 | 6 GB | 4 h |
| process_low | 2 | 12 GB | 4 h |
| process_medium | 6 | 36 GB | 8 h |
| process_high | 12 | 72 GB | 16 h |

Trim Galore retains `process_high`: **12 CPU, 72 GB, 16 h** before overrides.
These are generic nf-core reservations, not a measurement that trimming needs
72 GB. FastQC, sorting and statistical fitting likewise inherit their labels.
This revision does not silently replace those allocations with new guesses.
STAR index generation and the existing serial-tool overrides are identified in
the inventory. Scientific command arguments and imported module source stay intact.

## Katana downstream

| Work | Initial CPU / memory / time | Review note |
|---|---|---|
| CoSplice network fit | 12 / 72 GB / 16 h | CPU budget passed to BLAS/OpenMP; memory scales with the event matrix |
| MotifTrace | 6 / 36 GB / 8 h | Complete model family retained; requesting cores does not make serial traversal parallel |
| MAPP preparation and fit | 6 / 36 GB / 8 h | Tasks retain model inputs; no automatic switch to a high-memory partition |
| Network follow-up | 2 / 12 GB / 4 h | Separate task group, no network rediscovery |
| Gene-body coverage (optional) | 1 / 4 GB / 12 h | Serial traversal; disabled unless selected |
| Sequencing-QC summary, reference BED, audit | 1 / 4 GB / 2 h | Reuses existing products where available |
| Other downstream summaries | 2 / 12 GB / 4 h | Generic starting request, not study-specific peak data |
| Report | 2 / 12 GB / 8 h | Rendering and packaging; scientific results reused |

Downstream retries increase memory/time while keeping CPU count fixed. Upstream
retry policy is inherited from the upstream configuration. Task-level overrides
remain available; increasing a request does not guarantee scheduler capacity.

## Check allocation rather than only requested threads

Slurm may allocate additional CPUs to satisfy a partition's memory-per-CPU
constraint. For example, 72 GiB / 888 MiB is approximately 84 CPUs even when the
program requests 12 worker threads. That is a scheduler allocation, not 84 threads
used by Trim Galore. `Priority` denotes queue priority; it is not a quota error.

Compare Nextflow's trace/report with Slurm `ReqCPUS`, `AllocCPUS`, `ReqMem`, elapsed
runtime and step-level `MaxRSS` after a representative run. Do not lower memory
solely because CPU use is low. Limits and partitions belong to the site's config,
so no Dardel-specific memory ratio is built into the portable workflow.

Sources: [nf-core rnaseq base configuration](https://github.com/nf-core/rnaseq/blob/master/conf/base.config),
[Nextflow process selector priority](https://www.nextflow.io/docs/latest/config.html#process-selectors),
[Slurm memory-per-CPU behavior](https://slurm.schedmd.com/sbatch.html#OPT_mem-per-cpu).
