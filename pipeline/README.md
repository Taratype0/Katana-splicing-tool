# Katana bulk workflow

This entry point connects the existing bulk FASTQ/BAM analysis to Katana's
downstream calculations and the approved seven-section English report. It is
a Katana development workflow built with pinned nf-core modules; it is not an
official nf-core release. It does not run single-cell FASTQ processing.

## Run

Choose the [tested downstream image or standalone environment](../docs/RUNTIME_ENVIRONMENTS.md).
Container runs require Nextflow and the selected container engine on the host;
the downstream Python and R tools are supplied by the image. Nextflow requires Java 17 or newer and a Linux
execution environment (a Linux server, Linux workstation or an explicitly
configured Linux environment on Windows). The Python package's native Windows
interfaces are separate from this Nextflow execution requirement.

For a new study:

```bash
nextflow run ./pipeline \
  --mode complete \
  --input samples.csv --contrasts contrasts.csv \
  --fasta genome.fa --gtf annotation.gtf \
  --project_template configuration/project-template.json --flow_config configuration/flow.json \
  --katana_container /path/to/katana-downstream.sif \
  --dexseq_tools_container /path/to/kanata-dexseq-tools.sif \
  --outdir results -profile slurm,singularity
```

For completed upstream results, import the verified bundle once with
`katana-flow import`, then run only the missing downstream work:

```bash
nextflow run ./pipeline \
  --mode downstream-only --project project/project.json \
  --flow_config configuration/flow.json --outdir results -profile local
```

`--mode upstream-only` takes the same samples, contrasts and reference inputs
and generates the native upstream bundle without starting Katana downstream.
Samplesheets explicitly declare FASTQ or genome-BAM inputs. Reuse the existing
upstream schemas; do not generate every possible contrast automatically.

Add `-resume` to resume this run's Nextflow work/cache. Katana task receipts
also verify matching completed results before reuse. Preserve the original
inputs and prior run directories. Changing scientific inputs or a model
setting creates a different result identity.

The `slurm` profile has no application queue-size limit. Account, partition,
QoS and cluster policies determine admission. Set `--slurm_account` and
`--slurm_partition` for the target institution, or provide a site config with
`-c`. No Dardel paths, account names or private data are included. `local`,
`docker`, `singularity`, `apptainer` and `conda` profiles are independent choices;
a container profile requires an available Katana image via
`--katana_container`. No unpublished image is assumed to exist.
Container-enabled downstream/complete runs fail before submission if that full
runtime image is unspecified; they never silently use the host Python.
Use the [container build helper](../containers/downstream/README.md) with the
accepted wheel and dependency locks, then pass the actual image (a Docker tag,
registry digest or converted SIF as appropriate for the selected engine). The development
tag is a local build name, not a claim that an image has been published.

All Katana stages use the full image by default. An advanced site configuration
can override an execution group with
`params.katana_stage_containers = [network_fit: 'your-tested-image', mapp: 'your-tested-maei-image']`.
Control/preflight still require the full runtime; each scientific task also runs
`katana-flow runtime-check` inside its selected image before calculation. This
checks Python imports, the relevant R namespaces and declared worker runtime,
and writes `runtime-<task-id>.json`. A failed probe prevents that calculation.
The probe does not run a model or establish scientific completion.
Container tasks clear inherited Python/R library overrides and restrict selected
Python, Rscript, RSeQC executables and explicit R library directories to the
fixed image environment `/opt/conda`. The container launcher must be `katana-flow`
or `/opt/conda/bin/katana-flow`; an interpreter cannot declare its own host prefix
to be the image runtime. Absolute paths to mounted host Conda tools
are rejected before calculation, including project metadata and FIMO worker
configuration. Singularity/Apptainer Katana tasks use `--cleanenv`; the same
checks apply when calling `katana-flow task` directly in the image. Input data
and reference paths remain ordinary mounted resources. A standalone installation
can still explicitly select a separate installed native runtime.
Upstream module containers remain unchanged. Standalone downstream users can
install the documented Conda environment without using Nextflow or Docker.

Keep the configuration and its relative resource files in a dedicated
configuration directory. The workflow stages that directory and the project
directory so relative resources remain accessible inside task containers.
Keep Nextflow work and generated outputs outside both input directories.

## Execution and completion

Full mode validates required configured resources before upstream execution.
Upstream channel outputs are passed directly to bundle import. Downstream
tasks have separate output directories/catalogues. Each contrast or module
starts when its own named dependencies finish; network or MAPP fitting does
not hold unrelated comparisons behind a shared barrier. There is no fixed
four-task limit.

The report renderer reads validated results. It does not start a hidden refit.
Collection checks every planned task and its receipt before rendering. A failed
required analysis, missing binding resource or unavailable required fitting
backend blocks complete acceptance. An empty placeholder is not a completed
analysis. Historical partial reports remain readable through the existing
renderer, without being promoted to complete workflow runs.

New upstream analyses remove annotated pseudogenes and mitochondrial-encoded
genes before Salmon quantification/counting and native differential models.
STAR uses the unchanged full alignment annotation. The unfiltered GTF,
per-gene decisions and filtering receipt are retained in bundle provenance;
the filtered analysis GTF and matching transcript FASTA are exported. Existing
completed models are never silently refitted or relabelled.

Gene/event cards and Sashimi remain separate, user-selected follow-up analyses.
This entry point does not change the fixed seven-section report structure.

## Source identity and validation scope

`UPSTREAM_SOURCE.json` records the copied executable source closure, its base
Git commit and per-file SHA-256 values, including the already developed
per-contrast differential split. `INTEGRATION_CHANGES.json` records the explicit
upstream integration changes. Original source checkouts are not modified.
`UPSTREAM_LICENSE`, `UPSTREAM_CITATIONS.md` and native module metadata retain
upstream attribution. External tool/container licences still apply.

`tests/protocol_fixture.py` is an explicit transport-only fixture. Actual
Nextflow execution verifies task fan-out, dependency directories, skipped
groups, failure propagation and resuming. It calculates no biology and is
never selected by a production profile. The reference-filter tests separately
exercise real GTF filtering. Passing these checks is not a substitute for
real FASTQ/BAM end-to-end scientific acceptance or for the independently
tracked native dependency/container distribution work.

Run the native transport checks with an installed standard Nextflow launcher:

```bash
python -B pipeline/tools/check_native_transport.py --outdir /tmp/katana-transport-check
```

The output directory must be new. The Linux CI pins Nextflow 25.10.4 and its
distribution SHA-256. Its direct Java invocation includes the Java module
opens required by Nextflow cache serialization; calling a distribution JAR
without the launcher's JVM options is not an equivalent runtime check.
