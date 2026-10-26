# Quick start

Katana produces an RNA-splicing study report and a project you can explore in
the browser. Choose your starting files below. You select the comparisons;
Katana does not generate every possible group combination.

## Install

For the downloadable beta, follow [release installation](RELEASE_INSTALL.md).
For current source features, including automatic stage cleanup and the
configurable integrated-network sample minimum, install current `main` and
build a matching downstream image. The older beta downloads keep their
original behavior.

For the Python commands and browser, with Python 3.11 or newer:

```sh
git clone --branch main https://github.com/Taratype0/Katana-splicing-tool.git
cd Katana-splicing-tool
python -m venv .venv
. .venv/bin/activate
python -m pip install -e ".[next,flow,contracts]"
katana-open
```

On Windows, activate with `.venv\Scripts\Activate.ps1` instead. Run Nextflow on
Linux or a configured WSL environment with Java 17+ and Nextflow 25.10.4+.
The Python installation alone supports browsing/reporting; new analyses use
the [downstream container](../containers/downstream/README.md) or the
[standalone Conda environment](RUNTIME_ENVIRONMENTS.md).

For the full container workflow, build the downstream image from this checkout
and the [DEXSeq tools image](../containers/dexseq-tools/README.md). Set
`KATANA_IMAGE` and `DEXSEQ_IMAGE` to the resulting image names or SIF paths.
The upstream modules supply their other tool images automatically.
References and organism-specific analysis databases are configured separately.

Create your study directory outside the software checkout. Keep input tables
under `inputs/`, reads under `reads/`, and references under `references/`.
Commands below run from that study directory with Katana on `PATH`.

## 1. Start from FASTQ or BAM

Copy the [FASTQ sample sheet](../examples/nextflow/samples_fastq.csv) or
[BAM sample sheet](../examples/nextflow/samples_bam.csv) and replace its paths
and groups. Write one row per biological sample. Currently, multiple lanes
need to be merged before this entry.

FASTQ sheet columns:

```csv
sample,condition,fastq_1,fastq_2,strandedness,read_length,subject
Control_1,Control,../reads/Control_1_R1.fastq.gz,../reads/Control_1_R2.fastq.gz,auto,150,animal_01
Treated_1,Treated,../reads/Treated_1_R1.fastq.gz,../reads/Treated_1_R2.fastq.gz,auto,150,animal_04
```

The two rows illustrate the format; include all biological replicates.
For BAM, use `genome_bam` and `single_end` instead of the FASTQ columns, as in
the BAM template. Paths resolve relative to the sample sheet. `auto` measures
strandedness; known values can be `forward`, `reverse` or `unstranded`.

List only your intended comparisons in `inputs/contrasts.csv`:

```csv
contrast,treatment,control
Treated_vs_Control,Treated,Control
```

Effects follow treatment relative to control. The raw-data route currently
supports independent groups with `~ condition`; paired/batch-adjusted models
need a separate compatible analysis route. See [experimental design](NEXTFLOW_USAGE.md#4-contrastsheet-input-and-direction).

Prepare the [project resource template and analysis settings](NEXTFLOW_USAGE.md#6-downstream-analysis-configuration).
These identify the organism, reference, motif/pathway/CLIP resources and
analysis parameters. The full default report uses these resources; downloading
FASTQ alone does not provide them. You can [turn off selected analyses](ANALYSIS_AND_EXECUTION_OPTIONS.md)
and omit their calculation and dedicated resource requirements.

```sh
katana-flow run --mode complete \
  --input inputs/samples.csv --contrasts inputs/contrasts.csv \
  --fasta references/genome.fa --gtf references/annotation.gtf \
  --project_template inputs/project-template.json --flow_config inputs/flow.json \
  --profile docker --katana_container "$KATANA_IMAGE" \
  --dexseq_tools_container "$DEXSEQ_IMAGE" --outdir results
```

The same command handles a BAM sheet without realigning its BAMs. It runs
the two native Nextflow stages automatically, verifies each handoff, and
cleans their temporary work after success. No manual approval between tasks
is needed. Add `--keep-work` if you want all Nextflow work retained.

### Small studies

The upstream differential route requires at least two independent samples per
group. The integrated network defaults to **six total samples per comparison**,
not six per group. To use a two-plus-two study, set the explicit network
minimum in `flow.json`:

```json
{"network": {"parameters": {"minimum_network_samples": 4}}}
```

Merge this field into your full configuration. The network engine accepts a
minimum of at least three; the upstream two-per-group requirement still applies.
This changes the eligibility setting, not the number of observations or the
method's statistical power. The selected minimum and actual sample count are
recorded. This setting alone does not guarantee that other analyses are
estimable in a small study.

### On a cluster

Use `--profile slurm,apptainer --site-config site.config` instead of
`--profile docker`. Configure the partition, account, threads, memory and time
for your cluster. Katana does not choose a `memory` partition for you.
See the [site configuration example](../examples/nextflow/site.config).
Use an allocated compute environment for local execution, not a shared login
node. Slurm manages actual admission and concurrency.

## 2. Start from existing results

For a Katana upstream results bundle:

```sh
katana-flow import --bundle upstream-results/katana_bundle \
  --project-template inputs/project-template.json --outdir imported-project
katana-flow run --mode downstream-only \
  --project imported-project/project.json --flow_config inputs/flow.json \
  --profile docker --katana_container "$KATANA_IMAGE" --outdir results
```

For a prepared Katana project, start at the second command. For external
nf-core output, first map its result files, samples, comparisons and reference
using the [existing-results guide](EXISTING_RESULTS.md). Automatic discovery
from an arbitrary rnaseq/rnasplice/differentialabundance directory is not yet
available. A single pipeline's outputs may not cover every report analysis.

## 3. Open results or regenerate a report

After a successful managed run:

| Path | Use |
| --- | --- |
| `results/report/REPORT.html` | Open the study report in a browser. |
| `results/report/REPORT.zip` | Share the offline report and attachments. |
| `results/project/project.json` | Open the reusable analysis project in Katana. |
| `results/project/collection/` | Complete downstream results. |
| `results/pipeline_info/` | Execution traces, resource usage and logs. |

Start the interactive application with `katana-open` and open the project.
For a server, use `katana-open --no-open` and the documented
[server connection](INSTALL.md#server-access).

To regenerate a report from a completed collection whose source resources
remain available:

```sh
katana-flow report --collection completed-collection --outdir refreshed-report
```

This command assembles existing results without fitting new models. The
seven-section format and comparison selector stay the same. A managed run
already delivers its independent offline report; use that report after source
work cleanup. Retained-project browsing and target analysis are separate from
regenerating every report view from an archived collection's original paths.

## 4. Inspect selected genes or events

Select targets in the browser and export `targets.zip`. On your local analysis
machine or server, bind any required BAM/reference locations and run:

```sh
katana-targets inspect --request targets.zip
katana-targets run --request targets.zip --bindings bindings.json --outdir target-results
katana-targets results --outdir target-results
```

The result includes Gene/Event cards, available network and motif evidence,
requested Sashimi, individual reports and browser-readable data. Existing
differential analyses are reused. See [bindings and target commands](TARGET_COMMANDS.md).

## If a run stops

Keep its run directory. Fix the reported input or environment problem and
check [recovery rules](STORAGE_MANAGEMENT.md) before `--resume`: the managed
runner preserves the scientific configuration and source identity. Execution-only
CPU/memory/time/queue changes can be supplied through `--resource-config`; see
[analysis and recovery options](ANALYSIS_AND_EXECUTION_OPTIONS.md).
Native Nextflow `-resume` is available for direct runs with their work retained.
Removing work saves storage but removes that work's native resume cache.

For full formats and advanced commands, use [Nextflow usage](NEXTFLOW_USAGE.md),
[parameters](NEXTFLOW_PARAMETERS.md) and [storage management](STORAGE_MANAGEMENT.md).
