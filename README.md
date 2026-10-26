# Katana

**RNA splicing analysis, co-splicing networks and MotifTrace, with an interactive study report.**

Katana connects RNA-seq processing with expression, splicing, transcript usage and exon usage analysis. Start from bulk FASTQ, genome-aligned BAM, or existing analysis results. Explore the results in a browser and select genes or events for a closer look.

## Start here

| What do you have? | What to do |
| --- | --- |
| FASTQ or genome-aligned BAM | Run the [complete workflow](docs/QUICKSTART.md#1-start-from-fastq-or-bam). |
| Existing RNA-seq/splicing results | [Connect your results](docs/QUICKSTART.md#2-start-from-existing-results) and run the remaining analyses. |
| Completed Katana results | [Generate or open the report](docs/QUICKSTART.md#3-open-results-or-regenerate-a-report). |
| Genes or events you want to inspect | [Create a target report](docs/QUICKSTART.md#4-inspect-selected-genes-or-events), including Sashimi when BAMs are available. |

Read the **[quick start](docs/QUICKSTART.md)** for installation, example input tables and commands. The [full usage guide](docs/NEXTFLOW_USAGE.md) and [parameter reference](docs/NEXTFLOW_PARAMETERS.md) cover advanced settings.

## What Katana analyzes

| Analysis | Results |
| --- | --- |
| Design and quality | Sample groups, selected comparisons, sequencing QC, read support and gene-body coverage. |
| Splicing landscape | Event-type distributions, PSI heatmaps, PCA, sample patterns and replicate consistency. PSI is the fraction of transcripts including a splice feature. |
| Molecular differences | Differential gene expression (DEG), alternative splicing (DAS), transcript usage (DTU) and exon usage (DEXSeq), with effects, significance, gene labels and Top 10 tables. |
| Comparisons and pathways | Event overlaps and directions across selected comparisons; GO/KEGG and other supplied gene sets, with detected-gene backgrounds and contributing genes. Explicit interaction designs support additional fitted models. |
| Co-splicing networks | An independent network for each comparison: modules, complete gene/event membership, activity, condition associations, network graphs, enrichment and Hub rankings. |
| MotifTrace and RBP analysis | RBP identity enrichment, splice-boundary motif/CLIP evidence, and motif sites gained, lost or retained in reconstructed RNA forms. |
| Selected targets | Gene cards, Event cards, transcript and predicted protein consequences, read-backed Sashimi plots and individual HTML reports. |
| Processed single-cell/nucleus data | Supported donor-aware splicing, transcript, network and visualization workflows, plus compatible MARVEL and long-read result imports. The raw-read pipeline targets bulk RNA-seq. |

### MotifTrace

For each selected comparison, MotifTrace reconstructs local RNA forms for five splice-event types and scans supplied motif models. It links sites to sequences, genomic positions, events and genes, and calculates the supported matched-background and PSI-weighted analyses. Boundary/CLIP evidence and mature-RNA sequence consequences are shown separately.

The report displays up to ten significant motif models per comparison, including their sequences, original P/FDR and hit-gene links. Full model and site tables remain downloadable. Matching MAPP MAEI inputs can be prepared and fitted in the integrated workflow; compatible completed MAPP results can also be imported.

See [MotifTrace commands and Python API](docs/MOTIFTRACE.md).

### Co-splicing networks

Katana groups events with related PSI patterns using its WGCNA/TOM workflow. Each comparison has its own modules and matched downstream results. The browser provides network graphs, searchable member lists, module activity and significant enrichment with contributing genes.

Hub ranking uses absolute eigengene membership, `|kME|`, while retaining the signed values and complete ranking. Bootstrap is optional and off by default. Module size, merge height, sample minimum and execution resources are configurable. The current dense network method can require substantial memory for large event sets; see [network settings](docs/NETWORK_EXECUTION.md).

### Reports and browser

The default runs the full report. [Choose which optional analyses to run](docs/ANALYSIS_AND_EXECUTION_OPTIONS.md); disabled modules require no calculation or dedicated resources.

The study report has seven sections:

1. Experimental design, references and QC.
2. Alternative-splicing landscape.
3. DEG, DAS, DTU and DEXSeq results.
4. Cross-comparison results.
5. GO/KEGG and pathway enrichment.
6. Co-splicing networks.
7. RBP and MotifTrace.

A comparison selector keeps the report focused. Tables link to complete results, and the splicing display can switch between ordinary significant events and read-supported significant events (RS-DAS, the default).

The interactive browser adds adjustable thresholds, result filtering, gene/event selection and new analyses through the same Python service used by the command line. Selected-target reports are separate from the seven-section study report. Offline HTML can be shared without installing Katana.

## Run a study

After preparing the sample sheet, comparison sheet, reference and analysis-resource configuration, the storage-managed command runs Nextflow and delivers the report plus a reusable analysis project:

```sh
katana-flow run --mode complete \
  --input inputs/samples.csv --contrasts inputs/contrasts.csv \
  --fasta references/genome.fa --gtf references/annotation.gtf \
  --project_template inputs/project-template.json --flow_config inputs/flow.json \
  --profile docker --katana_container "$KATANA_IMAGE" \
  --dexseq_tools_container "$DEXSEQ_IMAGE" --outdir results
```

Use the matching current package, workflow and images described in the [quick start](docs/QUICKSTART.md). The BAM entry uses the same command with a BAM sample sheet. Nextflow schedules dependencies automatically; cluster accounts, partitions and resources belong in your site configuration.

Successful stages hand off their results before temporary work is removed. Failed work stays available for recovery. Add `--keep-work` to retain work for later Nextflow resume. See [storage management](docs/STORAGE_MANAGEMENT.md).

## What you receive

- An offline HTML report and downloadable figures, complete tables and report archive.
- A reusable analysis project for the Katana browser, CLI and Python API.
- Comparison-specific network and motif results, including member and hit-gene tables.
- Method settings, software/reference versions and execution logs.
- BAM/index resources retained separately for later read-backed target analysis.

## Installation and current status

Use the [release installation guide](docs/RELEASE_INSTALL.md) for downloadable packages, or [build from current source](docs/QUICKSTART.md#install). The browser and reporting tools work from the Python package; new native analyses also need the documented R/tools environment. The integrated workflow uses containers, while standalone downstream users can use [Conda](docs/RUNTIME_ENVIRONMENTS.md).

The source includes FASTQ/BAM workflows, existing-result adapters, downstream calculations, reports and selected-target commands. The current source has newer storage and configuration changes than the 25 October beta downloads; build matching images to use those changes throughout a run. [Validation records](docs/validation/entrypoints-runtime-20261025.md) and the [current integration review](docs/development/PIPELINE_REVIEW_20261026.md) describe tested routes and remaining work.

Katana is an independently maintained Nextflow workflow that reuses nf-core modules. Official nf-core inclusion is pending. This does not require manual task scheduling: Nextflow manages execution, dependencies and caching.

## Documentation

[Quick start](docs/QUICKSTART.md) · [Input and comparison tables](docs/NEXTFLOW_USAGE.md) · [Parameters](docs/NEXTFLOW_PARAMETERS.md) · [Existing results](docs/EXISTING_RESULTS.md) · [Target reports](docs/TARGET_COMMANDS.md) · [Methods](docs/METHODS.md) · [Installation](docs/INSTALL.md)

## License

Katana code is available under the [MIT License](LICENSE). Bundled third-party components retain their original licenses; see [third-party notices](THIRD_PARTY_NOTICES.md).
