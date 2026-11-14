# Sequencing QC without duplicate BAM scans

The first chapter keeps the approved Katana report structure and experimental
design tables. Sequencing-quality plots are displayed in the original interactive
MultiQC report, with its complete numerical data available beside it. MultiQC
aggregates completed tool outputs; it does not rerun alignment or scan BAM files.

New FASTQ/BAM runs aggregate the available FastQC, trimming and alignment outputs
in one Nextflow MultiQC task. BAM input cannot recover raw-read QC or reads removed
before the supplied BAM was created. Only measured modules are displayed.

## Reuse completed QC

Add `--qc_results /path/to/completed/qc` to the usual Nextflow command. This is also
available with `--mode downstream-only --bundle /path/to/katana_bundle`.
No hand-written JSON is needed. Supply one study's QC directory:

```text
qc/
  multiqc_report.html
  multiqc_data/
    multiqc_data.json
  sample_1/
    rnaseq_qc_results.txt
    raw_data_qualimapReport/
      coverage_profile_along_genes_(total).txt
  sample_2/
    ...
```

MultiQC HTML and JSON must belong to the same completed report. Its measured
sample names must match the declared samples (normal FastQC mate suffixes are
accepted). A software-version-only MultiQC report is not sequencing QC. Choose
the appropriate original report when an upstream tool produced multiple reports.
Per-sample Qualimap directories must use the exact sample IDs in the sample sheet.
Do not substitute a different alignment, reference or earlier experimental batch
just because its sample names are similar. These files are small; BAMs and FASTQs
are not copied by the QC importer.

Existing project users can declare `gene_body_coverage` assets with an explicit
`sample_id` and optional `qualimap_rnaseq_metrics` assets in the project manifest.
Declare the MultiQC HTML and `multiqc_data.json` as `sequencing_qc` assets with
`metadata.qc_type = "multiqc"`. The same installed `katana-flow plan` and
`katana-flow task` commands use this policy in standalone and Nextflow execution.

## Coverage policy

For each sample, Katana prefers a valid completed Qualimap coverage profile,
then Qualimap's saved raw curve in `multiqc_data.json`, then an existing
RSeQC/declared coverage profile. HTML alone can display the plots but cannot
provide a machine-readable profile; supply its matching JSON as shown above.
It never schedules another RSeQC
scan for a sample with a valid completed profile. Existing curves are included
even when `--with_genebody_coverage false` (the default) disables new scans.

`--with_genebody_coverage true` authorizes RSeQC for **missing samples only**.
With scanning disabled, missing samples are recorded as unavailable, not passed.
Advanced standalone configuration can use `qc.coverage_source: existing_native`
to forbid all new coverage scans. Duplicate sources of the same priority,
invalid profiles and unknown sample IDs require correction rather than being
silently combined. No sample is excluded based on the shape of the curve.

Qualimap native terminal-window metrics and RSeQC profile-derived end ratios are
different measurements. Source names, sample mapping and metric definitions are
preserved; they are not relabelled as interchangeable values.

The seven chapters, colours, molecular statistics and completed analysis results
are unchanged. Existing releases and reports remain intact; installing a new
release does not rewrite an earlier report or submit a scientific run.
