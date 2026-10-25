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

### Development workflow: unknown strandedness

In a `codex/package-development` source checkout, any chosen sample's
`strandedness` cell may be `auto`. The existing manual examples remain valid;
declared `forward`, `reverse` and `unstranded` rows keep their values and can
coexist with `auto` rows. Use the checkout's `pipeline/main.nf`. The frozen
`v2.0.0-beta.5.dev20261025` wheel/export does not support `auto`.

Only `auto` samples run RSeQC 5.0.3 inference on the shared audited, sorted and
indexed genome BAM, using BED12 derived from the matched unfiltered GTF.
FASTQ samples use their STAR BAM; supplied BAMs are not realigned. Defaults are
`--strandedness_sample_size 200000`, `--strandedness_mapq 30`,
`--stranded_threshold 0.8` and `--unstranded_threshold 0.1`. Fractions are
normalized over informative sense/antisense evidence. Undetermined calls fail,
and every selected comparison must have consistent resolved strandedness before
Salmon/DEXSeq/rMATS. Pipeline QC and bundle provenance retain the requested
values, decisions and native RSeQC evidence. See
[the detailed inference rules](../../docs/NEXTFLOW_USAGE.md#automatic-strandedness-in-the-development-workflow).

Use the [Nextflow usage guide](../../docs/NEXTFLOW_USAGE.md) for required files,
validation commands, supported designs and effect direction. The
[parameter reference](../../docs/NEXTFLOW_PARAMETERS.md) separates Nextflow
options from project resources and downstream settings.
