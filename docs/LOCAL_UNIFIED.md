# Unified local Katana — 12 November 2026

Start `launch_katana.bat` or `Start-Local-Katana.ps1`. The installed Windows
entry uses port 18775 and verifies the application identity before reusing a
service. The earlier sealed port-18773 application remains available through
the separately retained legacy launcher; its source and data were not replaced.

## One entry for projects, report folders and ZIPs

Choose **Open report** and select the extracted folder (recommended),
`report-manifest.json`, or `REPORT.zip`. ZIPs are unpacked into a checked local
cache; folders are read directly. Import creates a separate browsing workspace
with complete Parquet table copies. Keep the source folder available for its
original report and attachments. A lone HTML without complete tables/design
metadata cannot recreate the analysis. The supported report schema remains
`katana.approved_report/v1`; target-result ZIPs use the target-results importer.

**Existing project** opens the normal analysis workspaces in the same app.
**Results / DEG / DAS** also queries complete native project Parquet tables.
No raw FASTQ workflow or remote job starts when a report is opened.

## Full-population results

Filters cover comparison, exact gene IDs/symbols/event IDs, free text, original
FDR, absolute effect, direction, event type and recorded cell/donor/cluster
columns. Thresholds are unset initially, so nonsignificant rows remain visible.
Tables are paginated, but filtering, sorting, export and **Select all matching
rows** query the entire table. Original statistics and testing families are
never refitted by browsing.

The effect/FDR plot displays all valid matching rows, independently of the
table page. Above 30,000 points it aggregates the complete population into a
density grid; it does not sample or silently discard points. Drag a rectangle
to filter, then inspect individual points as the population becomes smaller.
Invalid/missing coordinates are counted separately and remain in the source
table. Zero FDR has a labelled plotting floor of 1e-300, not a changed source FDR.

## Local enrichment and target reports

**Pathway enrichment** accepts a local GMT or term/gene TSV/CSV library for
GO, KEGG, Reactome or custom pathways. Declare its species and matching ID
column; no database is silently downloaded. Foreground is all matching rows or
selected rows. The background uses all rows with valid original P/FDR in the
same comparison and biological scope, before significance/search/target-list
filters, intersected with the library. Events are collapsed to unique genes.
Multiple comparisons or cell types must be separated before ORA. Declared
hit-only/subset/Top-10 tables are rejected as backgrounds; undisclosed upstream
truncation cannot be independently certified.

The method is one-sided hypergeometric ORA, with BH across every eligible term,
including zero-hit terms. The UI can filter displayed FDR afterward. Exports
retain every tested pathway, original source/library hashes, foreground,
background and mapping losses. **Inspect genes** returns contributors to the
result explorer for target selection. Choosing a new test column is explicit;
use original P or adjusted-P values, never an effect size as testability evidence.

Select genes/events, choose **Add genes / Add events**, then open **Gene / Event
workspace**. **Generate cards locally** exports the exact selection internally
and creates saved-evidence cards without BAMs. Gene targets expand to available
measured events, including nonsignificant ones. Comparison scopes and browsing
filters are recorded. Cards retain gene/event evidence across their requested
comparisons, including recorded cell-type columns; this does not refit or pool
single-cell tests. Exact event cards require canonical identities and geometry.
Sashimi is optional and needs accessible BAM/index/reference bindings. The
existing independent Nextflow target task remains an advanced execution option;
preparing it is not submitting it.

## Bulk and single-cell use the same entry

Assay identity comes from the report `mode`, project `metadata.assay_mode` /
`metadata.assay`, or explicitly typed single-cell assets. Missing identity is
shown as **Assay not recorded**, never inferred from filenames. The same report
browser displays both. Single-cell metadata/columns retain donor, cell type,
cluster and condition information where supplied; an import cannot invent
missing metadata. Upstream single-cell models require biological donor/sample
replication; cells are not counted as independent biological replicates.

## Lossless storage by default

Fresh `katana-report`, Python report builds and presentation refreshes now
default to compaction. `--no-compact` is the explicit storage opt-out;
`--no-archive` still creates the full report directory without a second ZIP.
**Export compact report** makes a new storage-only copy of an imported report,
preserving its presentation and source. It deduplicates only SHA-identical
supported payloads, redirects portable links and losslessly encodes full HTML
table data. Native table bytes, nonsignificant rows, NA, P/FDR and source
provenance remain intact. `storage-audit.json` records the changes. Old delivery
receipts remain source provenance; storage transformation is not a new scientific
execution certification. Existing reports/archives are never compacted in place.

## Verification boundaries

Regression coverage includes full-population counts beyond pagination, density
accounting, original NA/nonsignificant values, filters/export/selection agreement,
ORA against an independent SciPy calculation, full-background and zero-hit-term
correction, ZIP traversal/tampering rejection, source preservation, target export,
native generated-BAM fixtures, and browser interactions. Synthetic fixtures are
software validation only. Local packaging does not establish biological validity
or a new Slurm/Nextflow cluster acceptance run.
