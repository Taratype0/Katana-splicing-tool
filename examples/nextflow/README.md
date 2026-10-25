# Raw-entry sheet templates

These CSV files are format templates, not sequencing data or completed results.
Copy one sample sheet and `contrasts.csv` into your study's `inputs/` directory,
then replace every sample ID, condition, file path, strandedness and read length.
Paths in sample sheets resolve from the sheet's directory: `../reads` and
`../alignments` refer to sibling study directories when the sheet is in `inputs/`.

The two sample sheets demonstrate paired-end bulk data from independent
biological samples, with three samples per arm. BAM must be genome-aligned and
match the declared reference. `reverse` and 150 are examples, not inferred
library settings. FASTQ lane rows must be merged into one file/pair per sample
before this raw-entry schema; do not duplicate a sample ID for each lane.

`contrasts_factorial.csv` illustrates twelve selected simple comparisons for
eight independent background/condition groups. It requires a corresponding
sample sheet containing those exact condition names; it is not compatible with
the two-group sheets unchanged. A condition label containing two factors does
not fit an interaction model. Only list the comparisons needed for the study.

Use the [Nextflow usage guide](../../docs/NEXTFLOW_USAGE.md) for required files,
validation commands, supported designs and effect direction. The
[parameter reference](../../docs/NEXTFLOW_PARAMETERS.md) separates Nextflow
options from project resources and downstream settings.
