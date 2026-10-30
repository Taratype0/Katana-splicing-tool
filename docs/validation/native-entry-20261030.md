# Native entry and delivery validation

This change makes the existing single-session `pipeline/main.nf` the documented
new-study entry, adds native client-project retention and final delivery tasks,
and gives native upstream products a consistent output layout. The scientific
modules and approved report renderer are unchanged.

## Checks

- Production main-script/import/downstream wiring with explicitly generated
  upstream test doubles: FASTQ-shaped and BAM-shaped entries each schedule ten
  transport tasks; existing-bundle entry schedules eight. Each repeat is fully
  cached. These checks do not run an aligner or claim new FASTQ/BAM scientific
  acceptance.
- Real installed Katana commands on generated result tables: twenty native
  tasks produce the retained client project and accepted seven-section report.
  A second run tests cache reuse. A separate generated-data session tests Nextflow-owned cleanup.
  The project is moved and read back after generated source inputs are removed.
- Standard publication tests check native result paths, unchanged file content,
  cached republication, absence of extra trimmed-read/index copies and absence
  of output symlinks into work. Historical publication tests remain separate.
- The final gate rejects mismatched report/project plans, changed report ZIPs,
  missing HTML and failed report acceptance. Retained-project tests exercise
  actual catalogue readback, relocation and target-data access.
- CI installs both Linux-built and Windows-built wheels on Linux, tests packaged
  R importer execution and resume, and runs the new single-session delivery
  check from the exported installed pipeline.

The release manifest records the exact candidate/main CI and artifact hashes.
Optional network, MotifTrace, MAPP and upstream scientific implementations retain
their earlier validation; generated delivery checks do not replace those checks
or claim that any in-progress study has completed. Existing studies were not
rerun for this change. Active managed sessions retain their original software
and cache identity rather than being overwritten during execution.
