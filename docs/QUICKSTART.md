# Quick start

Katana produces an RNA-splicing study report and a project you can explore in
the browser. Choose your starting files below. You select the comparisons;
Katana does not generate every possible group combination.

## Install

For Nextflow, use the [direct launch command](NEXTFLOW_SIMPLE.md). Nextflow,
Java and a container engine are sufficient; matching containers download
automatically. Mouse GRCm39 analysis databases also download automatically.
The Python installation below is for standalone analysis and the browser.

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

Create your study directory outside the software checkout. Keep sample and
contrast tables separate from the work and results directories.

## 1. Start from FASTQ or BAM

Copy the [FASTQ sample sheet](../examples/nextflow/samples_fastq.csv) or
[BAM sample sheet](../examples/nextflow/samples_bam.csv) and replace its paths
and groups. For multiple lanes, repeat the biological sample with a distinct
`lane` ID and each pair of original FASTQs. Nextflow performs the merge. CSV
and TSV are accepted.

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

Effects follow treatment relative to control. Independent comparisons use
`~ condition`. Paired comparisons declare `paired=true,blocking=subject` in
the contrast sheet and matching `subject` identities in the sample sheet; they
use `~ subject + condition`. See [experimental design](NEXTFLOW_USAGE.md#4-contrastsheet-input-and-direction).

Supply sample/contrast sheets, reference files and output paths. Internal JSON
is generated automatically. The selected modules use the versioned
[mouse databases](RESOURCE_DATABASES.md), or your explicit resource overrides.
Optional analyses can be disabled with the documented `--skip_*` switches.

```sh
nextflow -bg run Taratype0/Katana-splicing-tool -r v2.0.0-beta.5.dev20261103 \
  --input inputs/samples.csv --contrasts inputs/contrasts.csv \
  --fasta references/genome.fa --gtf references/annotation.gtf \
  --species Mus_musculus --reference_build GRCm39 \
  -profile docker --outdir results -work-dir work
```

The same command handles a BAM sheet without realigning its BAMs. It runs
one Nextflow DAG through report and client-project delivery. No manual stage
handoff or external coordinator is needed. Completed task outputs are published
as they finish; `-resume` reuses matching work after interruption. See the
[output layout and recovery guide](NATIVE_OUTPUTS.md).

### Small studies

The upstream differential route requires at least two independent samples per
group. CoSplice builds one shared network across all comparable study samples,
with no minimum biological-unit gate. Twelve independent units is a project
recommendation. Pairing is optional: declare matching subjects and paired
contrasts when appropriate. Six observations from three pairs remain three
independent biological units. Model terms that cannot be estimated are recorded
as unavailable. See [CoSplice](COSPLICE.md).

### On a cluster

Use `-profile slurm,apptainer -c site.config` instead of
`-profile docker`. Configure the partition, account, threads, memory and time
for your cluster. Katana does not choose a `memory` partition for you.
See the [site configuration example](../examples/nextflow/site.config).
Use an allocated compute environment for local execution, not a shared login
node. Slurm manages actual admission and concurrency.

## 2. Start from existing results

For a Katana upstream results bundle:

```sh
nextflow -bg run Taratype0/Katana-splicing-tool -r v2.0.0-beta.5.dev20261103 \
  --mode downstream-only --bundle upstream-results/katana_bundle \
  --species Mus_musculus --reference_build GRCm39 \
  -profile docker --outdir results -work-dir work
```

Resource options and analysis switches are the same as in the full-run example. For a prepared Katana project, replace `--bundle` with `--project
/data/project.json`. For external
nf-core output, first map its result files, samples, comparisons and reference
using the [existing-results guide](EXISTING_RESULTS.md). Automatic discovery
from an arbitrary rnaseq/rnasplice/differentialabundance directory is not yet
available. A single pipeline's outputs may not cover every report analysis.

## 3. Open results or regenerate a report

After a successful native run:

| Path | Use |
| --- | --- |
| `results/report/REPORT.html` | Open the study report in a browser. |
| `results/report/REPORT.zip` | Share the offline report and attachments. |
| `results/katana/project/project.json` | Open the reusable analysis project in Katana. |
| `results/katana/project/collection/` | Complete downstream results. |
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
seven-section format and comparison selector stay the same. A complete run
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
