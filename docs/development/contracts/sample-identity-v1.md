# Sample identity and upstream metadata adapter v1

Katana retains the upstream samplesheet and adds a small experiment metadata
table. Established nf-core workflows remain responsible for read preparation,
lane concatenation, alignment, quantification and their own input/QC checks.
The sample-identity validator checks identity and mappings; it does not replace those workflows or
introduce another read-processing implementation.

Users select the experimental comparisons they need. Recording conditions does
not generate every pairwise comparison. Selected-contrast preflight validates the explicitly selected contrasts and their designs; sample-identity validation generates no contrasts or statistical models.

## Two tables, different responsibilities

The upstream table retains one row per input run/lane. The Katana metadata table
has exactly one row per canonical `sample` ID, regardless of lane count.

| Identity | Meaning | Counting rule |
| --- | --- | --- |
| `sample` | Stable ID linking upstream inputs and downstream sample outputs | Several lane rows may share it; a label is not its identity |
| Input row | One read pair, single-end read file, BAM record or processed-sample record | Never an independent biological replicate by itself |
| `biological_unit_id` | Experimental unit from which the sample originates | Shared identity is preserved; unknown values remain unknown |
| `donor_id` | Person, animal or other explicitly recorded donor | Repeated conditions/specimens retain the same donor ID |
| `technical_replicate_of` | Canonical parent sample for a separately named technical replicate | Does not authorize an extra independent replicate or automatic merging |
| `display_name` | Human-readable experimental name | May change without changing the stable sample ID |

Use stable, explicit IDs; avoid spaces in upstream sample IDs. Do not infer
conditions, donors or replicate numbers from underscore-separated file names.
For independent animals, assign distinct biological-unit IDs. For aliquots or
repeated measurements from the same organism/experimental unit, retain their
shared unit and donor identities. A condition, time point or specimen label
alone does not create a new biological unit. Specimen identity is not separately
modeled by this thin adapter. Genuinely independently established cultures can
have different experimental-unit IDs while sharing their original donor; this
still does not establish independent donors. Selected-contrast preflight determines whether a paired or
repeated-measures design is valid.

### Upstream adapter columns

| `--format` | Required columns | Optional columns and scope |
| --- | --- | --- |
| `fastq` | `sample,fastq_1,fastq_2,strandedness` | `fastq_2` may be blank for single-end input; existing upstream columns are preserved |
| `bam` | `sample` and either `genome_bam` or compatibility alias `bam` | `bai` may be supplied; extra native columns are retained |
| `sample_ids` | `sample` | Processed-input identity mapping only; no raw-read processing |

The FASTQ column convention follows the [official nf-core/rnaseq usage
documentation](https://nf-co.re/rnaseq/latest/docs/usage/). That workflow groups
repeated sample IDs for read concatenation. The exact accepted schema and merging
behavior must still be qualified against the pinned upstream pipeline/version.
Sample-identity validation does not certify every nf-core release or pipeline using these columns.

The BAM profile accepts native `genome_bam` and the Katana compatibility alias
`bam`. When both are nonempty they must identify the same input; conflicting
values are rejected. `transcriptome_bam` never substitutes for genome-aligned
input. The current rnaseq reprocessing sheet also retains its other native
columns and alignment-skipping configuration. Its supported provenance and
quantifier restrictions still apply. This profile does not certify a complete
runnable sheet for every source pipeline/version. Successful sample-identity validation does
not establish arbitrary BAM suitability, genome/reference compatibility or
complete Katana input generation.

### Katana experiment metadata columns

Required headers are `sample,condition,biological_unit_id,donor_id`.

- `sample` must uniquely identify the same canonical sample as the upstream table.
- `condition` is an explicit group/condition label, not a filename-derived guess.
- `biological_unit_id` may be blank when unknown; validation reports that
  uncertainty and does not invent a biological replicate count.
- `donor_id` may be blank for bulk data. Retain it whenever known. scRNA/snRNA
  donor-level inference needs donor mapping and later method/design qualification;
  a metadata check does not confer donor-level statistical evidence.
- Optional fields are `display_name`, `batch`, `cell_type`,
  `technical_replicate_of` and named `covariate_*` columns. These record declared
  information. They do not fit a model or automatically select covariates.

Put no lane-specific information in a duplicated metadata row. If a sample has
two lanes, reference its single metadata record from both upstream rows.
Metadata rows with duplicate canonical IDs, conflicting mappings or missing
counterparts are diagnostics to resolve explicitly.
If native input columns also declare a condition/group, donor or biological
unit, their sample-level values must agree with the Katana metadata mapping.
Validation preserves the original samplesheet bytes and does not rewrite it
into a new format.

## Multiple lanes and technical replicates

For re-sequencing of the same upstream sample, repeat its `sample` ID in the
upstream sheet and retain one metadata row. Each read path must identify the
intended file without duplicate reuse or mate collisions. The adapter returns
read-only `upstream_groups`, retaining input row and mate order. nf-core performs
the actual concatenation according to its qualified workflow/version.

A separately identified library or technical replicate retains a distinct
`sample` ID and may name its parent in `technical_replicate_of`. Do not silently
convert several libraries into independent biological replicates. Conversely,
do not merge different IDs simply because they share a donor or biological unit.
Sample-identity validation performs no cross-ID read merge and no count-table summation.

The BAM adapter records sample/input identity. It does not execute `samtools
merge`, deduplicate alignments, establish the relationship between multiple BAMs
or make a quantification route scientifically eligible. Multiple BAM rows for
one sample are rejected as an unresolved mapping; supply the appropriate
qualified upstream representation rather than applying FASTQ merge rules.

## What is checked, and what is not

Sample-identity validation checks table structure, sample identity, metadata coverage, consistent lane
grouping, duplicate/mate input references, explicit technical-replicate links and
the requested assay scope. An optional exchange-contract project mapping check checks that the
manifest and table identities agree. It does not create or rewrite the manifest.

It does not open or validate FASTQ/BAM contents, inspect read groups, infer donor
identity, determine experimental independence from names, fetch data, run QC,
merge files, aggregate counts or submit Nextflow/Slurm jobs. Keep upstream and
method-specific eligibility checks in their existing responsible workflow.
Duplicate input checks compare lexically normalised declared locators, not file
contents, symlink targets or all relative/absolute aliases. They do not establish
that differently named files contain different biological material.

A successful result means the declared metadata mapping is internally
consistent. It is not a scientific acceptance receipt or execution approval.
Reports retain `execution_authorized: false`, `merges_executed: 0`,
`scientific_tools_executed: 0` and `contrasts_generated: 0`.

## Run the synthetic examples

The examples use invented names and placeholder paths only. There are no FASTQ,
BAM or quantitative data files, no invented completion evidence and no file
hash claims. Run the commands below from the repository root. The contracts share
the existing `katana_next` package rather than a second nested source tree.

Linux/macOS:

```bash
export PYTHONPATH="$PWD"
python -B -m katana_next.contracts samplesheets \
  --input examples/contracts/sample-identity-v1/fastq-multilane/upstream.csv \
  --metadata examples/contracts/sample-identity-v1/fastq-multilane/metadata.csv \
  --format fastq --assay bulk_rna
```

Windows PowerShell:

```powershell
$env:PYTHONPATH = (Get-Location).Path
python -B -m katana_next.contracts samplesheets `
  --input examples/contracts/sample-identity-v1/fastq-multilane/upstream.csv `
  --metadata examples/contracts/sample-identity-v1/fastq-multilane/metadata.csv `
  --format fastq --assay bulk_rna
```

Use `--project path/to/project.json` to check against an existing exchange-contract v1 manifest.
That optional manifest must describe the same canonical sample IDs; this option
does not transform an arbitrary older Katana project schema.

The API is also read-only:

```python
from katana_next.contracts import validate_samplesheets

report = validate_samplesheets(
    "examples/contracts/sample-identity-v1/fastq-multilane/upstream.csv",
    "examples/contracts/sample-identity-v1/fastq-multilane/metadata.csv",
    input_format="fastq",
    assay_mode="bulk_rna",
)
```

Assay labels are `bulk_rna`, `single_cell_rna` and `single_nucleus_rna`.
The latter two use `sample_ids` for processed-input identity mapping at this input boundary;
they do not enable a new raw single-cell/nucleus FASTQ pipeline.
