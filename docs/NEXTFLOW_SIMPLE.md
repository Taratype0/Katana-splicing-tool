# Run Katana with Nextflow

Use the [matching package, controller and containers](INDEPENDENT_RUNTIME.md).
The 28 October release includes these direct options. The older 25 October
downloads retain their previous interface.

Supply samples, the comparisons you want, references, and an output directory.
Katana generates its internal project template and downstream settings. You do
not need to create `study.yml`, `project-template.json` or `flow.json`.

Use the native Nextflow entry for new studies:

```bash
nextflow run pipeline/main.nf -params-file examples/nextflow/params.complete.yaml \
  -profile slurm,apptainer -c examples/nextflow/kth.config \
  --slurm_account YOUR_ALLOCATION --outdir /data/study/results -work-dir /data/study/work
```

Edit the template paths first. A single Nextflow DAG connects all stages and
publishes the report and independent client project. Use `-resume` with the same
launch directory and work directory after interruption. No Python coordinator
is needed. Historical `katana-flow run` sessions retain their old recovery
command and installation until they finish; do not switch an active run's code.

## 1. FASTQ or BAM to analysis and report

From the repository checkout:

```bash
nextflow run pipeline/main.nf \
  -profile slurm,apptainer \
  -c examples/nextflow/kth.config \
  --slurm_account YOUR_ALLOCATION \
  --input /data/study/samples.csv \
  --contrasts /data/study/contrasts.csv \
  --fasta /data/reference/genome.fa \
  --gtf /data/reference/annotation.gtf \
  --reference_build GRCm39 \
  --annotation_release YOUR_ANNOTATION_RELEASE \
  --species 'Mus musculus' \
  --study_name 'My study' \
  --gene_sets /data/resources/pathways.gmt \
  --rbp_gene_sets /data/resources/rbp_gene_sets.gmt \
  --motif_database /data/resources/motifs.json \
  --clip_peaks /data/resources/matched_clip_peaks.tsv \
  --katana_container /data/images/katana-matching-version.sif \
  --dexseq_tools_container /data/images/dexseq-tools.sif \
  --outdir /data/study/results \
  -work-dir /data/study/work \
  -resume
```

Replace paths, allocation, species/build and release with your own. Resource
formats are the existing Katana resource formats: a filename alone does not
convert an arbitrary motif or CLIP database. See the resource sections of [the parameter reference](NEXTFLOW_PARAMETERS.md).

For a workstation, use `-profile local,docker`, Docker image references instead
of SIF paths, and omit `-c kth.config` / `--slurm_account`. For another cluster,
use your site's Slurm config. Queue choices are site configuration; analyses do
not automatically select a memory partition. The supplied `kth.config` uses
`shared`. It is a Katana example, not an automatic import of nf-core's `pdc_kth`.

The sample sheet selects the input type. FASTQ:

```csv
sample,condition,source,fastq_1,fastq_2,strandedness,read_length
control_1,CONTROL,fastq,/data/reads/control_1_R1.fastq.gz,/data/reads/control_1_R2.fastq.gz,auto,150
treated_1,TREATED,fastq,/data/reads/treated_1_R1.fastq.gz,/data/reads/treated_1_R2.fastq.gz,auto,150
```

Genome-aligned BAM:

```csv
sample,condition,source,genome_bam,strandedness,single_end,read_length
control_1,CONTROL,genome_bam,/data/bam/control_1.bam,auto,false,150
treated_1,TREATED,genome_bam,/data/bam/treated_1.bam,auto,false,150
```

These two-row fragments illustrate columns, not an adequate replicated study.
Use one biological sample identity across its lanes, and distinct identities
for independent samples. Complete six-sample examples are in
[samples_fastq.csv](../examples/nextflow/samples_fastq.csv) and
[samples_bam.csv](../examples/nextflow/samples_bam.csv). Detailed lane, batch and
donor rules remain in [the input guide](NEXTFLOW_USAGE.md).

Choose comparisons explicitly:

```csv
contrast,treatment,control
TREATED_vs_CONTROL,TREATED,CONTROL
```

Effects are oriented treatment minus control (expression: treatment/control).
The workflow does not enumerate every group combination. `--contrast_ids`
optionally selects a subset for downstream work; the upstream comparison CSV
still controls upstream tests. To avoid unwanted upstream tests, put only the
desired comparisons in that CSV.

## 2. Use one optional Nextflow parameter file

Download/copy [params.complete.yaml](../examples/nextflow/params.complete.yaml),
edit paths and optional settings, then run:

```bash
nextflow run pipeline/main.nf \
  -params-file /data/study/params.complete.yaml \
  -profile slurm,apptainer \
  -c examples/nextflow/kth.config \
  --slurm_account YOUR_ALLOCATION \
  -work-dir /data/study/work -resume
```

YAML is optional. Every scalar option can instead be passed as `--option value`;
CLI values override the parameter file. Use `-c` for executor, queues and task
resource requests. Use `-params-file` or `--...` for analysis inputs/settings.
The [complete native option table](NEXTFLOW_NATIVE_OPTIONS.md) lists every
declared parameter and its default, including STAR, Salmon, rMATS and SUPPA.

## 3. Select analyses

| Option | Default | Effect |
| --- | --- | --- |
| `--skip_network` | false | Omit independent per-comparison network calculations and their report results. |
| `--skip_motiftrace` | false | Omit MotifTrace sequence consequences and its motif summaries. |
| `--skip_motif_binding` | false | Omit motif/CLIP binding analysis; useful when matched CLIP is unavailable. |
| `--skip_mapp` | false | Omit MAPP fitting and its report results. |
| `--skip_enrichment` | false | Omit general gene-set enrichment; networks still need gene sets. |
| `--skip_rbp_identity` | false | Omit RBP gene-set annotation/enrichment. |
| `--skip_landscape` | false | Omit the optional AS landscape analysis. |
| `--skip_cross_comparison` | false | Omit cross-comparison intersections; per-comparison molecular analysis remains. |
| `--skip_interaction` | false | Disable explicitly configured interaction fits. No design or biological purpose is guessed. |
| `--skip_sequencing_qc` | false | Omit downstream sequencing-QC summaries. Upstream FastQC has its own `--skip_fastqc`. |
| `--with_genebody_coverage` | false | Opt in to extra RSeQC gene-body BAM traversal. |
| `--with_network_bootstrap` | false | Opt in to bootstrap; the normal network uses continuous kME without bootstrap. |

Disabled components are recorded as not selected. The approved seven-chapter
layout stays unchanged; disabled analyses do not launch placeholder computation.
DEG/DAS/DTU/DEXSeq import and the core molecular summary remain part of the
downstream workflow. Missing resources for an enabled analysis cause an early
message identifying the required option; they do not silently disable it.

For example, without matched CLIP, append `--skip_motif_binding`. Without any
motif analysis resources, explicitly select `--skip_motiftrace
--skip_motif_binding --skip_mapp`. To omit all gene-set-dependent analyses,
select `--skip_enrichment --skip_network --skip_rbp_identity`.

## 4. Thresholds and resources

| Option | Default | Meaning |
| --- | --- | --- |
| `--deg_fdr` | 0.05 | Imported DEG adjusted-P cutoff. Does not refit the upstream model. |
| `--deg_log2fc` | 1.0 | Minimum absolute expression log2 fold change. |
| `--das_fdr` | 0.05 | Caller FDR cutoff in event QC/classification. |
| `--das_dpsi` | 0.10 | Minimum absolute delta PSI. |
| `--dtu_fdr` | 0.05 | Imported DTU/DEXSeq adjusted-P cutoff. |
| `--event_set` | stable_significant | RS-DAS; `significant` selects ordinary DAS where supported. RS-required modules retain that boundary. |
| `--rmats_mode` | JCEC | Native result universe imported from the upstream bundle: JC or JCEC. Existing projects keep their declared mode. |
| `--network_min_module_size` | 20 | Minimum events per discovered module. |
| `--network_merge_cut_height` | 0.25 | Module eigengene merge cut height, strictly between 0 and 1. |
| `--contrast_ids` | all declared | Comma-separated downstream comparison IDs; never creates new comparisons. |
| `--gene_sets` | unset | GO/KEGG/custom GMT used by selected enrichment and networks. |
| `--rbp_gene_sets` | unset | RBP identity gene sets, using the project's gene identifier convention. |
| `--motif_database` | unset | Katana-compatible motif models, including multiple models. |
| `--clip_peaks` | unset | Matched CLIP resource for selected binding analysis. |
| `--fimo_runtime` | generated | Optional existing worker descriptor. Normally generated as `python3` in the selected runtime. |

Read support remains at least 10 reads in **every** sample of the selected
comparison. Network defaults remain 20 events per module, merge height 0.25,
per-comparison fitting, kME ranking and no bootstrap. Module size and merge height have direct custom options. This entrypoint change
does not alter the algorithms. Explicit factorial design settings can still
use the existing advanced configuration route documented below.

Resource bytes are hashed once during configuration generation. Without
provider metadata the recorded version is `sha256:<digest>` (content identity),
not an invented database release. Supply real provider metadata in the same
optional Nextflow YAML, for example:

```yaml
resource_metadata:
  motif_database:
    source_url: https://YOUR_DATABASE/EXACT_RELEASE
    version: YOUR_RELEASE
    species: Mus musculus
```

Conflicting declared hashes/species/builds are rejected. A database or FASTA/GTF
is a data resource; installing Python packages does not supply it.

## 5. Existing results and advanced settings

An already imported Katana project needs no new study/configuration JSON:

```bash
nextflow run pipeline/main.nf \
  --mode downstream-only --project /data/existing-inputs/project.json \
  -profile slurm,apptainer -c examples/nextflow/kth.config \
  --slurm_account YOUR_ALLOCATION \
  --katana_container /data/images/katana-matching-version.sif \
  --outdir /data/new-run/report -work-dir /data/new-run/work -resume
```

Its samples, resources, reference and JC/JCEC mode come from the imported
project. Flat thresholds and analysis switches generate its execution settings.
New resource paths are not replacements for assets in an existing project.
An arbitrary nf-core output folder still needs the appropriate import adapter;
this interface does not guess results, sample identities or effect directions.

Existing explicit `--project_template` plus `--flow_config` remain supported
for advanced raw runs (including explicitly specified interaction designs).
`--flow_config` alone is supported for an existing project. In these routes,
the explicit configuration controls downstream settings; flat downstream
switches do not override it. A reused `--flow_plan` retains its original settings.
Do not mix these routes expecting two configurations to be merged.

Generated settings are published in `OUTDIR/pipeline_info/configuration/`.
The final report, result tables and client-readable outputs follow the existing
[standard output layout](NATIVE_OUTPUTS.md). `-resume` reuses retained matching Nextflow work; it cannot reuse
deleted cache. `--storage_mode lean` reduces publication duplication but does
not itself delete work. Keep work until the report and portable result package
are verified, then use the existing documented cleanup route.

Target gene/event deep dive remains a separate user-selected follow-up; see
[the target commands](TARGET_COMMANDS.md). It is not automatically run for every
gene during the standard report.
