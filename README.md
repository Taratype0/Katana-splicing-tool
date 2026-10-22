# Katana

Katana connects differential RNA-splicing evidence with candidate review,
isoform interpretation, regulatory hypotheses and coordinated splicing networks.

**Development checkpoint — 22 October 2026.** T01–T08 define the adopted scope
and engineering baseline. T09 adds tested exchange contracts. The complete
reusable analysis package and integrated Nextflow distribution are being built.
This checkpoint is not a new desktop release or an approved nf-core pipeline.

## Agreed architecture

| Layer | Responsibility |
|---|---|
| Bulk upstream | FASTQ or eligible genome BAM → qualified analysis inputs, using Nextflow and reusable nf-core modules |
| Katana Python analysis core | Shared QC, differential evidence, candidates, enrichment, isoforms, MotifTrace and network analysis |
| Reports and interfaces | English reports, CLI/API and interactive views consuming the same analysis results |

Katana MotifTrace connects sequence/motif changes, transcript usage and RBP/CLIP
evidence. Katana-Splice-Network discovers coordinated splicing modules and
supports fixed-module condition effects, hubs and validation. These remain
distinct evidence layers; neither establishes causal regulation on its own.

The upstream scope is short-read **bulk RNA-seq**. Qualified scRNA/snRNA routes
belong to Katana downstream analysis of compatible existing results, aligned BAM
or counts; raw single-cell FASTQ processing is outside the first upstream scope.

The adopted inventory has **54 formal primary workflow entries**. Formal scope,
default/optional selection and migration/packaging/acceptance status are separate.
An optional feature remains formal. Inclusion in scope is not a claim that its
new package interface has passed acceptance.

## Available in this public checkpoint

- [T01–T08 decisions and delivery boundaries](docs/project-baseline/SCOPE.md).
- [Development milestones](docs/project-baseline/MILESTONES.md).
- [T09 read-only contracts, synthetic examples and tests](development/README.md).
- Preserved legacy application code and its existing Git history. Root-level
  application/build files describe that legacy application; use `development/`
  for the new contract checkpoint.

The new T09 source files are byte-preserving copies of the isolated development
baseline. Historical experimental scripts, local configurations, datasets,
databases and internal acceptance logs are excluded from the new publication.
Scientific workers will be published with their dependency and licence review
as each package migration is completed.

The pre-checkpoint main commit is
[`f19453c`](https://github.com/Taratype0/Katana-splicing-tool/tree/f19453c3cad403eb9dc326dc65b2c21a7e91cfcc).
Existing local Katana installations and completed results are not modified.

MIT applies to repository-owned code; third-party methods and future bundled
resources retain their own attribution and licences.
