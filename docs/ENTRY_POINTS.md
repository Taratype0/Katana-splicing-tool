# Choose a Katana entry point

**Use completed results when you already have them. Use the bulk Nextflow
workflow when reads still need processing, and the target workflow after
choosing genes or events.** These entries share the package but perform
different work; successful report assembly does not certify a new raw-data run.

For full input/contrast formats and execution choices, use the
[Nextflow usage guide](NEXTFLOW_USAGE.md), [parameter reference](NEXTFLOW_PARAMETERS.md)
or [Chinese quick start](NEXTFLOW_QUICKSTART_ZH.md).

| Starting material | Command or interface | Work performed |
| --- | --- | --- |
| Bulk FASTQ | Nextflow `--mode complete` | Upstream processing → Katana inputs → selected downstream analyses → report. |
| Genome-aligned bulk BAM | Nextflow `--mode complete`, with a BAM samplesheet | BAM/reference audits, quantification and differential results → downstream analyses → report. |
| FASTQ/BAM, with only upstream output needed | Nextflow `--mode upstream-only` | Produces the native input bundle and stops before Katana downstream calculations. |
| Existing native bundle/project, with analyses still to calculate | Import if needed, then Nextflow `--mode downstream-only` | Computes the selected downstream analyses. |
| Completed analysis files or a portable report | `katana-report` | Assembles or refreshes the seven-section report without statistical fits. |
| Selected genes or exact events | Browser target basket or `katana-targets` | Cards, requested Sashimi, independent reports and returned interactive products. |

**Verified now:** generated bulk FASTQ and BAM complete execution and
all-cached repeats; completed-result reporting; installed Windows/Linux target
restore and reuse.
See the [current versions and acceptance boundaries](validation/entrypoints-runtime-20261025.md).

The raw-read workflow is for **bulk RNA-seq**. Existing single-cell/single-nucleus
results can use their supported downstream adapters, but that does not make the
bulk FASTQ/BAM command a raw scRNA/snRNA pipeline. Assay- and method-specific
eligibility still applies.

## Install one matching software set

The [downloadable release](RELEASE_INSTALL.md) supplies the matching wheel,
Docker archives, exported workflow and dependency locks.

Use `main` and the [installation guide](INSTALL.md). Further development
continues on `codex/package-development`. The
Nextflow host needs Java/Nextflow and its container engine. Katana Python, R
and external tools run in the downstream image; upstream tools use their
declared images. Standalone downstream users can instead use the
[locked Linux Conda environment](RUNTIME_ENVIRONMENTS.md).

Set `KATANA_IMAGE` to the matching built image or SIF, and `DEXSEQ_IMAGE` to the
[separately built upstream DEXSeq tools image](../containers/dexseq-tools/README.md).
Match **source revision, wheel SHA256, image identity and dependency locks**;
the package version string alone is insufficient. The current accepted
installed artifacts are bound to source revision `098ceae`; their full hashes
are in the [runtime record](validation/entrypoints-runtime-20261025.md#software-and-runtime).
Earlier wheels/images keep their original identities. No public
container-registry or PyPI release is implied.

An installed wheel can export the workflow without a source checkout:

```sh
katana-flow export-pipeline --outdir tools/katana-pipeline
```

The examples below use that exported `tools/katana-pipeline/main.nf`. Keep
configuration/input directories separate from work and output directories.

For the development addition `strandedness=auto`, use a
`codex/package-development` checkout's `pipeline/main.nf` in place of that export;
see [the development checkout command](NEXTFLOW_USAGE.md#2-installation-and-directory-layout).
The published `v2.0.0-beta.5.dev20261025` wheel and original workflow export remain
frozen and still require manual strandedness. No update to those assets is
implied by the new source support.

## 1. New bulk FASTQ or BAM

Both use the same command; the samplesheet declares the actual files:

```sh
nextflow run tools/katana-pipeline/main.nf -profile docker \
  --mode complete --input inputs/samples.csv \
  --contrasts inputs/contrasts.csv \
  --fasta references/genome.fa --gtf references/annotation.gtf \
  --project_template inputs/project-template.json \
  --flow_config inputs/flow.json \
  --katana_container "$KATANA_IMAGE" \
  --dexseq_tools_container "$DEXSEQ_IMAGE" \
  --outdir results -work-dir work
```

The raw-entry sheet requires unique sample IDs, `condition` and `strandedness`,
plus `fastq_1`/optional `fastq_2` **or** `genome_bam`. Supply a
positive `read_length` for each sample or `--read_length` for the study. BAM
input also requires explicit `single_end=true/false`. `source`, when supplied,
is `fastq` or `genome_bam`; a transcriptome BAM does not substitute for it.
The matched-reference and BAM-content checks still run. File paths must already
be local or mounted. Do not use a downstream metadata-adapter example as a
replacement for this raw-entry schema.

The development workflow accepts `auto` alongside explicit
`forward`/`reverse`/`unstranded` rows. Only `auto` rows add RSeQC inference after
the shared BAM audit/sort/index, using the matched unfiltered GTF converted to
BED12. FASTQ rows use their STAR alignments; supplied BAMs are not realigned.
All selected samples must resolve and each selected comparison must retain
consistent actual strandedness before Salmon/DEXSeq/rMATS. Undetermined calls
stop execution without an unstranded fallback. Native evidence and requested/
resolved settings are retained in pipeline QC and bundle provenance; see
[the inference rules and defaults](NEXTFLOW_USAGE.md#automatic-strandedness-in-the-development-workflow).

The contrast sheet uses `contrast,treatment,control`: include only the
comparisons you want, with that direction. The project template declares
species and matched reference/database resources; `flow.json` configures
downstream modules. See [integrated execution and resource requirements](INTEGRATED_PIPELINE.md).

For upstream output only, use `--mode upstream-only` and omit
`--project_template`/`--flow_config`. Keep samples, contrasts, FASTA/GTF and
the container arguments. It stops after the native bundle, not after a report.

## 2. Existing upstream results, with downstream computation

A completed Katana native bundle can be imported into a new workspace:

```sh
katana-flow import --bundle upstream/katana_bundle \
  --project-template inputs/project-template.json --outdir imported-project

nextflow run tools/katana-pipeline/main.nf -profile docker \
  --mode downstream-only --project imported-project/project.json \
  --flow_config inputs/flow.json --katana_container "$KATANA_IMAGE" \
  --outdir downstream-results -work-dir downstream-work
```

An existing valid `project.json` can be used directly. Unadapted nf-core output
directories are not interchangeable with that manifest: first bind their
actual files, samples, comparisons, references and methods through the
[existing-result interfaces](EXISTING_RESULTS.md). Importing completed evidence
does not certify an analysis that was never performed.

For unchanged work, preserve the original plan path, inputs, launch directory,
Nextflow cache and work directory. Validate the plan, then use `--flow_plan`
and `-resume`; see [the exact commands](INTEGRATED_PIPELINE.md#resume-an-unchanged-downstream-plan).
A new executable/image can change Nextflow fingerprints. Do not rewrite a
receipt or force a cache hit to hide that change.

## 3. Completed results: report only

```sh
katana-report --project results.json --outdir exports/study-report \
  --compact --no-archive

katana-report --from-report portable-report --outdir exports/refreshed-report \
  --compact --no-archive
```

`results.json` can bind several completed sources explicitly. The first route
assembles their files; the second refreshes an existing portable report. Add
`--path-map paths.json` when the completed inputs moved to another machine.
Use `--archive downloads/study-report.zip` instead of `--no-archive` for a
download package. Both routes preserve native statistics and the approved
seven-section format. They do not calculate missing networks, MotifTrace or
MAPP families. See [completed sources and missing-analysis planning](EXISTING_RESULTS.md).

## 4. Selected genes/events and returned products

Start the browser with `katana-open`, browse completed results, add genes or
exact events to the target basket, and export its request ZIP. The same request
can run locally or in a server allocation:

```sh
katana-targets inspect --request targets.zip
katana-targets plan --request targets.zip --bindings bindings.json
katana-targets run --request targets.zip --bindings bindings.json --outdir target-results
katana-targets results --outdir target-results
```

Bindings map resource IDs to that machine's BAM/index/reference files. Cards
can use saved evidence; requested Sashimi requires its verified read resources.
The request does not embed large BAMs or assume a particular cluster. Return
the result ZIP to the browser for verified import, or restore it for subsequent
CLI use without the original project database:

```sh
katana-targets restore --results completed-targets.zip --outdir restored-targets
katana-targets plan --request restored-targets/TARGET_REQUEST.zip \
  --bindings bindings.json --outdir restored-targets
```

Restore verifies the archive and completed products. Reuse still requires the
same verified methods/evidence/resources; restoration alone cannot promise a
cache hit under a different dependency environment. See [target commands](TARGET_COMMANDS.md)
for manual target selection, per-target scheduler commands and result reuse.

## Validation and platform scope

The [current entry/runtime record](validation/entrypoints-runtime-20261025.md)
separates complete-mode execution, cache reuse, target analysis and retained
real-study evidence. The [completed-result report record](validation/completed-results-runtime-20261025.md)
describes the full portable study report, while the
[preserved raw-entry record](validation/raw-entry-preserved-evidence-20261025.md)
retains earlier execution and recovery scopes.

Windows supports the Python/browser/report route. Native tools were checked on
Linux x86-64 with standalone Conda and Docker. Site deployment is separate:
these checks do not certify every Slurm/Apptainer installation, macOS/ARM or
external BAM/reference combination. The workflow is not yet an approved nf-core
pipeline. Installing a runtime does not supply references or validate a new
experiment automatically.
