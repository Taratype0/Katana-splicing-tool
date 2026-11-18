# Katana

Start with the [Nextflow quick start](docs/NEXTFLOW_SIMPLE.md) or [standalone installation](docs/INDEPENDENT_RUNTIME.md).

**RNA splicing analysis, co-splicing networks and MotifTrace, with an interactive study report.**

Katana connects RNA-seq processing with expression, splicing, transcript usage and exon usage analysis. Start from bulk FASTQ, genome-aligned BAM, or existing analysis results. Explore the results in a browser and select genes or events for a closer look.

## Start here — Two Independent Yet Converging Workflows

Katana v2.7 provides two independent entrypoints that converge on the same statistical engine, publication figures, and portable `REPORT.html` / `REPORT.zip` archive:

| Entrypoint | Starting Inputs | Command |
| --- | --- | --- |
| **Route 1: Standalone Python Package** *(no Nextflow, Docker, or R needed)* | Completed `rMATS` (+ optional `DESeq2`, `MultiQC`, `FASTA`, `GTF`, `BAM`) or Single-Cell Smart-seq3 junction/PSI tables | `katana run-bulk --rmats ./rmats --samples samples.tsv -o ./report`<br>`katana run-sc --events sc_events.tsv --cells sc_cells.tsv -o ./sc_report` |
| **Route 2: End-to-End Nextflow Pipeline** *(raw reads or aligned BAMs)* | Bulk `FASTQ` (`fastq_1,fastq_2`) **or** genome-aligned `BAM` (`genome_bam`) + `contrasts.csv` + `FASTA` + `GTF` | `nextflow run pipeline/main.nf --mode complete --input samples.csv --contrasts contrasts.csv --fasta genome.fa --gtf genes.gtf --outdir results` *(or `katana pipeline --suite-outdir ./report -- ...`)* |
| **Convergence Bridge** | Any completed Nextflow upstream bundle (`katana_bundle`) | `katana run-bulk --bundle ./results/katana_bundle -o ./report` |
| **Interactive Local Compute & Sashimi** | Any `REPORT.zip` or report directory | `katana view ./report/REPORT.zip` |

### Python Quick Start (Standalone Package)

```bash
pip install -e ".[standalone]"

# 1. Bulk splicing suite directly from rMATS (+ optional DESeq2, MultiQC 5'-3', FASTA & GTF)
katana run-bulk \
  --rmats ./rmats_contrasts \
  --samples ./samples.tsv \
  --multiqc ./multiqc_data \
  --fasta ./genome.fa --gtf ./annotation.gtf \
  --outdir ./katana_bulk_report

# 2. Single-cell splicing suite (Precision-Kernel Event-PSI UMAP, Δ_kNN, SC-SpliceNet & ORA)
katana run-sc \
  --events ./sc_differential_splicing.tsv \
  --cells ./cell_metadata.tsv \
  --psi-matrix ./cell_psi.parquet \
  --expression ./gene_logtpm.parquet \
  --outdir ./katana_sc_report

# 3. Open interactive local workspace (Live ORA, PSI Landscape Lab, 170-PWM Scanner, 4-Track GeneCard & Sashimi)
katana view ./katana_bulk_report/REPORT.zip
```

Read the **[quick start](docs/QUICKSTART.md)** for full input table schemas. The [full Nextflow usage guide](docs/NEXTFLOW_USAGE.md) and [parameter reference](docs/NEXTFLOW_PARAMETERS.md) cover cluster profiles (`docker`, `singularity`, `apptainer`, `slurm`).

## What Katana analyzes

| Analysis | Results |
| --- | --- |
| Design and quality | Sample groups, selected comparisons, sequencing QC, read support and gene-body coverage. |
| Splicing landscape | Event-type distributions, PSI heatmaps, PCA, sample patterns and replicate consistency. PSI is the fraction of transcripts including a splice feature. |
| Molecular differences | Differential gene expression (DEG), alternative splicing (DAS), transcript usage (DTU) and exon usage (DEXSeq), with effects, significance, gene labels and Top 10 tables. |
| Comparisons and pathways | Event overlaps and directions across selected comparisons; GO/KEGG and other supplied gene sets, with detected-gene backgrounds and contributing genes. Explicit interaction designs support additional fitted models. |
| Co-splicing networks | One network across comparable study samples: modules, complete gene/event membership, activity, selected condition associations, network graphs, enrichment and Hub rankings. |
| MotifTrace and RBP analysis | RBP identity enrichment, splice-boundary motif/CLIP evidence, and motif sites gained, lost or retained in reconstructed RNA forms. |
| Selected targets | Gene cards, Event cards, transcript and predicted protein consequences, read-backed Sashimi plots and individual HTML reports. |
| Processed single-cell/nucleus data | Supported donor-aware splicing, transcript, network and visualization workflows, plus compatible MARVEL and long-read result imports. The raw-read pipeline targets bulk RNA-seq. |

### MotifTrace

For each selected comparison, MotifTrace reconstructs local RNA forms for five splice-event types and scans supplied motif models. It links sites to sequences, genomic positions, events and genes, and calculates the supported matched-background and PSI-weighted analyses. Boundary/CLIP evidence and mature-RNA sequence consequences are shown separately.

The report displays up to ten significant motif models per comparison, including their sequences, original P values, regional FWER or window FDR, and hit-gene links. Full model and site tables remain downloadable. Matching MAPP MAEI inputs can be prepared and fitted in the integrated workflow; compatible completed MAPP results can also be imported.

See [MotifTrace commands and Python API](docs/MOTIFTRACE.md).

### Co-splicing networks

Katana groups events with related PSI patterns using its WGCNA/TOM workflow. Discovery uses all comparable study samples; the selected comparisons describe how the same modules vary with condition. The browser provides network graphs, searchable member lists, module activity and significant enrichment with contributing genes.

Hub ranking uses absolute eigengene membership, `|kME|`, while retaining the signed values and complete ranking. Bootstrap is optional and off by default. Pairing is optional; repeated observations retain their biological-unit identities. Twelve independent units are recommended, without a hard sample-count gate. Module size, merge height, event count and execution resources are configurable; see [CoSplice commands and settings](docs/COSPLICE.md).

### Reports and browser

The default runs the full report. [Choose which optional analyses to run](docs/ANALYSIS_AND_EXECUTION_OPTIONS.md); disabled modules require no calculation or dedicated resources.

The study report has seven sections:

1. Experimental design, references and QC.
2. Alternative-splicing landscape.
3. DEG, DAS, DTU and DEXSeq results.
4. Cross-comparison results.
5. GO/KEGG and pathway enrichment.
6. RBP and MotifTrace.
7. Co-splicing networks across the study.

Sequencing QC uses the original interactive MultiQC report. [Import completed QC](docs/QC_REUSE.md)
to reuse Qualimap coverage and other measured metrics without another BAM scan.

A comparison selector keeps the report focused. Tables link to complete results, and the splicing display can switch between Precision-Supported DAS (`PS-DAS`, the default: effective-length normalized replicate-quantile + AS-NMD asymmetry rescue + Jeffreys Binomial Fisher precision) and unfiltered Ordinary DAS.

A whole-study CoSplice result stays the same when the comparison selector changes. Its module composition, activities, associations and complete gene/event lists close the report. Imported older networks keep their recorded scope. [CoSplice](docs/COSPLICE.md) and [MotifTrace](docs/MOTIFTRACE.md) also have standalone commands and Python APIs.

The interactive browser adds adjustable thresholds, result filtering, gene/event selection and new analyses through the same Python service used by the command line. Selected-target reports are separate from the seven-section study report. Offline HTML can be shared without installing Katana.

## Run a study

Supply the sample sheet, comparison sheet, references and analysis resources. One Nextflow run schedules read processing, differential analysis, Katana, the report and the reusable browser project. Internal configuration files are generated automatically:

```sh
nextflow -bg run Taratype0/Katana-splicing-tool -r v2.0.0-beta.5.dev2026111401 \
  -profile docker --input inputs/samples.csv --contrasts inputs/contrasts.csv \
  --fasta references/genome.fa --gtf references/annotation.gtf \
  --species Mus_musculus --reference_build GRCm39 \
  --outdir results -work-dir work
```

Nextflow fetches the matching tool images. Mouse GRCm39 defaults also fetch versioned GO, KEGG, RBP, motif and reference CLIP resources; see [databases](docs/RESOURCE_DATABASES.md) for sources and overrides. The BAM entry uses the same command with a BAM sample sheet. Nextflow schedules dependencies automatically; cluster accounts, partitions and resources belong in your site configuration.

Nextflow passes results between tasks through channels; no external Python coordinator or study-specific monitoring script is required. Add `-resume` to reuse this run's completed tasks. Work is separate from deliverables and remains available for recovery. See [storage management](docs/STORAGE_MANAGEMENT.md).

The same [parameter template](examples/nextflow/params.complete.yaml) also works
with `nextflow run pipeline/main.nf -params-file params.yaml -profile slurm,apptainer -c site.config`.
Explicit command-line parameters override the file. Use `-resume` to recover
the same run. The older `katana-flow run` command remains for resuming historical
managed sessions; it is not the entry for new studies.

Results are published as stages finish:

```text
results/
  rnaseq/           # FastQC, trimming/alignment diagnostics, BAMs, Salmon, matrices
  rnasplicing/      # Native rMATS and SUPPA outputs
  differential/     # DEG, DTU and exon-usage results
  katana/           # Per-task results and retained browser/target-analysis project
  report/           # REPORT.html, REPORT.zip, figures and complete tables
  pipeline_info/    # Nextflow execution records and final delivery receipt
```

For an existing upstream bundle, select `--mode downstream-only --bundle /data/katana_bundle`.
For an imported Katana project, use `--mode downstream-only --project /data/project.json`.
These entries do not launch alignment or upstream differential fitting.

## What you receive

- An offline HTML report and downloadable figures, complete tables and report archive.
- A reusable analysis project for the Katana browser, CLI and Python API.
- Comparison-specific network and motif results, including member and hit-gene tables.
- Method settings, software/reference versions and execution logs.
- BAM/index resources retained separately for later read-backed target analysis.

## Installation and current status

Use the [release installation guide](docs/RELEASE_INSTALL.md) for downloadable packages, or [build from current source](docs/QUICKSTART.md#install). The browser and reporting tools work from the Python package; new native analyses also need the documented R/tools environment. The integrated workflow uses containers, while standalone downstream users can use [Conda](docs/RUNTIME_ENVIRONMENTS.md).

The source includes FASTQ/BAM workflows, existing-result adapters, downstream
calculations, reports and selected-target commands. Install the matching release
workflow and runtime described in the installation guide. [Native delivery
validation](docs/validation/native-entry-20261030.md) distinguishes generated
integration tests from existing scientific validation; no completed study was
rerun for this release.

Katana is an independently maintained Nextflow workflow that reuses nf-core modules. Official nf-core inclusion is pending. This does not require manual task scheduling: Nextflow manages execution, dependencies and caching.

## Documentation

[Quick start](docs/QUICKSTART.md) 路 [Input and comparison tables](docs/NEXTFLOW_USAGE.md) 路 [Parameters](docs/NEXTFLOW_PARAMETERS.md) 路 [Existing results](docs/EXISTING_RESULTS.md) 路 [Target reports](docs/TARGET_COMMANDS.md) 路 [Methods](docs/METHODS.md) 路 [Installation](docs/INSTALL.md)

## License

Katana code is available under the [MIT License](LICENSE). Bundled third-party components retain their original licenses; see [third-party notices](THIRD_PARTY_NOTICES.md).
