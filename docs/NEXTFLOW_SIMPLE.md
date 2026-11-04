# Run Katana with Nextflow

Supply sample and contrast tables, a genome FASTA, an annotation GTF, and output
paths. One Nextflow run connects read processing, statistical analysis, Katana,
and the report. No separate Python installation, coordinator, JSON or YAML is
needed for this entry.

The server needs Nextflow >=25.10.4, Java >=17, and Apptainer/Singularity or
Docker. Nextflow downloads the versioned tool images automatically. For mouse
GRCm39, Katana also downloads the [versioned analysis databases](RESOURCE_DATABASES.md)
required by your selected analyses. FASTA and GTF remain explicit inputs.

## 1. FASTQ or BAM to analysis and report

```bash
nextflow -bg run Taratype0/Katana-splicing-tool -r v2.0.0-beta.5.dev2026110403 \
  -profile singularity -c /path/to/site.config \
  --input /data/study/samples.csv \
  --contrasts /data/study/contrasts.csv \
  --fasta /data/reference/genome.fa \
  --gtf /data/reference/annotation.gtf \
  --reference_build GRCm39 --annotation_release 'GENCODE vM34' \
  --species Mus_musculus --study_name 'My study' \
  --resource_cache /data/cache/katana/resources \
  --digest_cache /data/cache/katana/checksums \
  --outdir /data/study/results -work-dir /data/study/work
```

Your site configuration selects the executor, allocation and queues. On Dardel,
use `-profile singularity,pdc_kth -c custom.config --project YOUR_ALLOCATION`
in place of the generic profile/site options above. The packaged profile accounts
for memory-driven shared CPU allocation; do not also load an older `pdc_kth.config`.
Your custom configuration remains yours. See [task resources](RESOURCE_ALLOCATION.md).
Katana does not automatically switch tasks to a memory partition.
On a workstation use `-profile docker` and omit the site config/account options.
The same command works with the released pipeline directory in place of the
GitHub repository name.

`-bg` backgrounds the Nextflow driver. Follow `.nextflow.log` and the reports in
`results/pipeline_info`. After interruption, rerun the same command with
`-resume`, from the same launch directory and using the same release and work
path. Add `--clean_work true` when you want Nextflow to remove task caches after
successful final delivery; failed runs keep their cache for recovery. Formal
results and BAMs are retained. Shared image/database caches remain reusable.

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
for independent samples. For multiple lanes, repeat the sample row with a
distinct `lane` column and the original FASTQ paths; Nextflow merges the lanes
automatically before QC. Complete six-sample examples are in
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

## 2. Existing results and optional configuration

For a Katana upstream bundle, replace the raw input arguments with
`--mode downstream-only --bundle /data/katana_bundle`; keep the species,
reference-build and resource settings. For an imported project, use
`--mode downstream-only --project /data/project/katana.project.json`; its existing
resources and identities are reused. These entries do not redo alignment.

All routine settings use command-line options. A Nextflow `-params-file` is
optional. The advanced `--flow_config` and `--project_template` inputs remain
available for specialized designs, but are not prerequisites for a normal run.

## 3. Select analyses

| Option | Default | Effect |
| --- | --- | --- |
| `--skip_network` | false | Omit whole-study CoSplice and its condition/contrast associations. |
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
whole-study fitting with contrast associations, kME ranking and no bootstrap. Module size and merge height have direct custom options. This entrypoint change
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
