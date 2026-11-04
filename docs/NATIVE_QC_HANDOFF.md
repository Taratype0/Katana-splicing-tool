# Native QC handoff

The pipeline exports FastQC (raw, trimmed or supplied-BAM records), Trim Galore,
STAR alignment logs and the all-sample strandedness summary. The report reader
accepts each of these producer formats. The strandedness table keeps every
sample, its requested and resolved orientation, inference method and recorded
measurements. Explicitly supplied orientations retain empty measurement fields.
No inference is repeated when creating the report.

The table is available with the chapter-one QC downloads. The approved report
chapters, scientific thresholds, analysis methods and sample contrasts are unchanged.
Unknown QC formats and incomplete or conflicting sample identities still produce
an actionable input error.

Release acceptance now includes this table in both the installed standalone and
native Nextflow generated-data runs, through report delivery and cache reuse.
Unit regressions also cover malformed summaries, mixed FASTQ/BAM projects and
relocated report data. Previously the complete generated-data runtime fixture
explicitly disabled sequencing QC, leaving this producer/consumer gap untested.

A real project's successful plan preflight establishes input readiness only;
it does not establish that its downstream analyses or final report are complete.

DEXSeq report handoff also accepts the native-header-authoritative direction
marker emitted by the pipeline importer. Original effect columns, direction and
P/FDR are retained in the existing volcano, Top 10 and downloadable tables.

Independently scheduled comparisons are selected by their own result scope when
building molecular panels and Top 10 tables. Combined legacy tables remain
supported through their comparison column. A result from another comparison
is never substituted for missing data. The runtime check verifies all four
molecular layers for both generated comparisons.
