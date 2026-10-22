# Synthetic sample-identity examples

These CSVs contain invented sample IDs and placeholder paths. They contain no
real experiment information or read/alignment contents. Validation is metadata
only; no example establishes upstream pipeline acceptance or biological power.

| Directory | Adapter format / assay | What it illustrates |
| --- | --- | --- |
| `fastq-multilane` | `fastq` / `bulk_rna` | Two samples, two lane rows each, one metadata row per sample |
| `bam` | `bam` / `bulk_rna` | Two canonical native `genome_bam` records; metadata only, not a complete runnable rnaseq BAM samplesheet |
| `repeated-donor` | `sample_ids` / `bulk_rna` | Two donors measured in two conditions, retaining the same organism/unit ID across conditions |
| `technical-replicate` | `sample_ids` / `bulk_rna` | A distinct technical-library sample linked to a parent, sharing its unit/donor |
| `single-cell-processed` | `sample_ids` / `single_cell_rna` | Processed donor/cell-type identities retained without claiming raw processing or new donor inference |

No directory provides contrasts. Users choose the contrasts needed; Selected-contrast preflight will
validate their designs. This validator neither expands all pairs nor invokes a statistical
model. The two-sample examples are intentionally small identity illustrations,
not a recommendation for an adequately replicated differential analysis.

See [the sample identity contract](../../../docs/development/contracts/sample-identity-v1.md) for
columns, scope and CLI/API commands. Run those commands from the package source
directory (the repository root), replacing the example
directory and `--format`/`--assay` values as needed.
