# Katana

**RNA splicing analysis, MotifTrace and co-splicing networks.**

Katana integrates RNA-seq splicing, expression and transcript-usage results to help researchers identify candidate genes and inspect their sequence changes, regulatory context and read support. Analyses can be explored in a browser, run from the command line, or called from Python. Completed results can be assembled into study reports.

## Inputs

A project records sample groups, biological replicates, donors, reference versions and the comparisons to analyze. Provide the files required by the selected analysis:

| Input | Used for |
| --- | --- |
| rMATS JC or JCEC and compatible DEG, DTU or DEXSeq results | Differential splicing, expression, transcript/exon usage and cross-comparison analysis. |
| Salmon `quant.sf`, compatible count, PSI or junction matrices | Transcript composition, sample comparisons and eligible statistical or network analyses. |
| GTF, genome/transcript FASTA and relevant annotation resources | Event and transcript annotation, RNA reconstruction, enrichment, motif and protein-consequence analysis. |
| BAM files and indexes, when needed | Sashimi, coverage and junction inspection, target-gene checks, and supported network or single-cell counting routes. |
| Supported processed single-cell/nucleus data, MARVEL or long-read results | The corresponding downstream analyses or completed-result imports. |

Analyses based on existing results do not require BAM files. Upstream raw-data processing continues to use established nf-core workflows; Katana connects their outputs with the experimental design. Input schemas and compatible references are specified for each analysis.

## Analysis features

| Area | Analyses |
| --- | --- |
| Samples and experimental design | Check files, sample identities, biological replicates, donors, groups and effect direction; select the comparisons of interest. |
| Quality and overview | Inspect read support, replicate agreement, coverage bias and uncertainty; explore PSI, event-type distributions, PCA, correlations and heatmaps. |
| Differential and cross-condition analysis | Integrate DAS, DEG, DTU and DEXSeq while retaining source statistics; compare exact-event overlaps, concordant effects and direction reversals. Additional expression and exon-usage interaction calculations require suitable counts, designs and R tools. |
| Candidate genes and events | Combine supporting results, browse complete candidate sets, save selections, and inspect gene cards and event details. |
| Functional enrichment | Analyze declared gene sets against appropriate backgrounds using GO/pathway resources; perform ranked pathway analysis where the design and resources support it. |
| Transcript and protein consequences | Examine transcript composition and usage changes; inspect annotated ORFs, UTRs, predicted NMD and optional Pfam/domain or protein-feature evidence. |
| Read validation | Generate replicate-aware Sashimi plots and inspect target-gene coverage, exons, junctions and RNA-level perturbation evidence. |
| Single-cell/nucleus downstream analysis | For eligible data, aggregate by donor and cell type and analyze splicing measurements, differences, embeddings, networks and transcript usage. Import supported completed MARVEL and long-read results. |

### MotifTrace

MotifTrace links splicing changes to RNA sequence and RBP evidence at three levels:

1. **RBP identity and target sets:** annotate candidate RBPs, splicing factors and declared target sets.
2. **Regulatory context around splice events:** examine neighborhood motifs, matched-background enrichment and compatible reference CLIP evidence.
3. **RNA sequence consequences:** reconstruct local reference RNA forms for SE, A3SS, A5SS, MXE and RI; inspect gained, lost or retained motif sites and new splice junctions, and assess motif availability using PSI.

Browse results by gene, event or motif and export RNA FASTA, site tables, genomic coordinates, statistics and figures. Custom motif libraries are supported. Eligible existing transcript/exon-usage and expression results can also inform motif-availability analysis.

MotifTrace supports new sequence analyses and imports of compatible completed sequence or parallel MAPP results. Imported method identities and P/FDR values are retained. General one-command MAPP fitting for arbitrary new projects is not currently provided. Local RNA reconstruction and motif predictions do not establish measured RBP binding.

The development package provides `from katana_motiftrace import MotifTrace`,
`katana-motiftrace` commands, and a shared MotifTrace mode selector in the
browser interface. See [MotifTrace inputs, API and commands](docs/MOTIFTRACE.md)
for the six analysis/import modes, requirements and examples.

### Katana-Splice-Network

Network analysis accepts eligible rMATS events, competing junctions extracted from BAMs, or external PSI matrices, retaining the units and read evidence of each input route.

The existing WGCNA/TOM route identifies co-splicing modules and produces memberships, sample activity scores, condition associations, hubs and network views. Supporting analyses include bootstrap stability and membership confidence, module function, structural and sensitivity diagnostics, projection onto existing modules, and preservation in independent cohorts. Optional routes, including the original DJExpress author core, require their own inputs.

Discovery uses an eligible sample cohort. Condition associations and external validation retain their respective designs; membership, association, stability and functional evidence are presented separately.

### One-click reports

The approved report follows seven sections: design and QC; AS landscape; DEG/DAS/DTU/DEXSeq; cross-comparison results; GO/KEGG and pathways; co-splicing networks; and RBP/MotifTrace. It includes per-comparison Top 10 details, full-result links and a report-wide ordinary DAS / RS-DAS choice, with RS-DAS as the default. See the [approved report specification](docs/REPORT_SPECIFICATION.md).

The development package provides `katana-report` and `from katana_next.reports import build_report` for completed projects, and the same renderer is available as `reports.study` in the application. See [report inputs, commands and outputs](docs/REPORTS.md). Two additional export routes remain available:

- **Analysis and reproducibility package:** an HTML/PDF summary, figures, parameters, file inventory and ZIP attachments.
- **Completed study-report import:** package an existing supported report and its attachments while retaining its original scope and results.

Reports retain significant and non-significant results and document actual filtering and exclusions. Unrun analyses and missing matching downstream results are identified. Report generation reads completed analyses; it does not infer experimental purpose or rescue, generate a validation programme, or rerun upstream processing. Selected gene/event dossiers remain separate interactive analyses.

### Integrated analysis workflow (development)

The [Nextflow workflow](pipeline/main.nf) connects bulk FASTQ/BAM processing,
Katana calculations and the fixed seven-section report. Expensive comparisons
and network stages are separate executor tasks. It supports upstream-only,
downstream-only and complete execution; completed work is retained with
`-resume`. The `katana-flow` commands validate resources, prepare independent
tasks, verify their products and enforce report completion. The existing
`katana-report` command remains a renderer for already completed results.

This integration is under execution acceptance, and is not an official nf-core
release. A missing required calculation or resource blocks a complete report.
Native MAPP fitting needs explicitly prepared, matching transcript quantities
and motif-count inputs; these are not substituted with event PSI. See
[integrated execution and current limits](docs/INTEGRATED_PIPELINE.md).

### Browser interaction

The browser interface lets users open projects, inspect samples and designs, select comparisons and parameters, follow task status, filter result tables, save candidate genes or events, inspect gene cards and figures, and export results. MotifTrace, networks and reports are accessible through their analysis pages.

The browser connects to a Katana Python service running locally or on the analysis server. The service accesses project files and runs calculations; the browser displays results and controls jobs. Offline report HTML remains available for reading completed results and does not start new analyses. The CLI and Python interface support scripted and server-side use. [MotifTrace](docs/MOTIFTRACE.md) has dedicated API and command-line examples; other modules remain available through the shared workflow commands as their packaged interfaces are finalized.

The browser is the current interface and release target. Existing desktop and mode-chooser entry points remain available for compatibility; a new standalone desktop distribution is not planned.

The **Gene / Event workspace** keeps selected targets across result pages and exports a portable request. `katana-targets` reads the same request locally or on an analysis server to create Gene cards, exact Event cards, native Sashimi, independent HTML reports and complete interactive tables. Completed expression, splicing, transcript, network and MotifTrace evidence is reused with its original statistics. BAMs and references are bound separately on the execution machine. See [commands, portable inputs and outputs](docs/TARGET_COMMANDS.md) and [selection and evidence rules](docs/TARGET_REPORTS.md).

## Outputs

Depending on the analysis, Katana produces browsable result tables, interactive views and PNG/SVG/PDF figures. Sequence analyses also export FASTA and site files; network analyses export memberships, activities, associations, hubs and stability results. Reports and ZIP packages retain parameters, methods and input sources.

## Install and use

The current version is installed from a source checkout. Use Python 3.11 or newer in a separate environment:

```sh
git clone https://github.com/Taratype0/Katana-splicing-tool.git
cd Katana-splicing-tool
python -m venv .venv
```

Activate the environment with `.venv\Scripts\Activate.ps1` in Windows PowerShell or `. .venv/bin/activate` on Linux. Then install and open the browser interface:

```sh
python -m pip install -e ".[next]"
katana-open --runtime-root ../katana-runtime
```

`katana-open` opens the browser by default. On a server, add `--no-open` to start or reuse the Python service and print its URL without opening a window; see [server access](docs/INSTALL.md#server-access). Open an existing Katana project or migrate a compatible legacy project, check the samples and inputs, then select the analysis. Explicit `--mode desktop` and `--mode choose` retain the existing compatibility modes. The command line and Python interface call the same analyses; for example, export a reproducibility package from an existing project:

```sh
python -m katana_next.cli --project path/to/katana.project.json run reports.bundle
```

Inputs are currently registered in a project configuration. Sample/contrast-sheet onboarding and upstream-result adapters are being unified for new users. See [installation](docs/INSTALL.md) and [usage](docs/USAGE.md) for migration, formats and resources. R, motif-scanning environments, references and databases are configured for the selected modules.

## Development status

The analyses described above are implemented in the existing Katana application and include both new calculations and completed-result imports. The source installation provides interface, CLI and Python access.

Current development packages these capabilities for general distribution, including independent installation, dependency separation, input adapters, local/server execution and reports for new projects. Release packaging and complete Nextflow downstream integration still require acceptance work; a published PyPI installation command is not yet provided here.

Further technical details are in [methods](docs/METHODS.md) and [development documentation](docs/development/ARCHITECTURE.md).

## License

Katana's own code uses the [MIT License](LICENSE). Third-party components retain their original licenses; see [third-party notices](THIRD_PARTY_NOTICES.md).
