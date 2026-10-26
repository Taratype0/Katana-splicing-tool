# Execution options and recovery validation

This change removes internal network event-count/time guards, adds optional
analysis selection, shares MAPP preparation within each comparison, reduces
repeated file reads/copies, and separates scheduler resources from scientific
settings during managed resume. Default report content and scientific settings
are retained. See [usage](../ANALYSIS_AND_EXECUTION_OPTIONS.md).

## Checks completed

- Public regression suite, targeted import/managed-run checks and network
  extensions: 1,097 distinct tests passed and 12 skipped across the recorded
  runs, taking the latest result for each test after repairs. No unresolved
  test failures remained. This is not a claim of a single uninterrupted run.
- Real Nextflow 25.10.4 scheduling on generated shell fixtures: dependency
  ordering, parallel dispatch, required-failure handling and cached resume.
- Real resource retry: first attempt deliberately exited 137; retry used two
  CPUs and 32 MB instead of one CPU and 16 MB. Both successful tasks were
  cached after a memory/time-only configuration change.
- Generated study tables ran through six selected core analysis tasks,
  collection, the report renderer and report acceptance with optional modules
  disabled. Missing disabled motif resources were omitted at import; missing
  core resources still failed validation.
- Shared MAPP preparation covered all 14 windows. The bundled author's kmer
  counter ran on generated sequences; matched shared and single-window tables
  and quantity provenance agreed. Native model fitting was not repeated.
- File checks detected same-size content changes even with restored mtime on
  Windows. Retained hard links survived removal of the work-side link; the
  cross-filesystem copy fallback preserved bytes.
- The full Nextflow entry compiled. All 89 copied upstream source files and
  37 declared integration changes passed the source check.
- The approved real report's content hash remained unchanged.

## Scope

These checks used generated fixtures and existing software tests. No Primary
culture or public biological dataset was rerun, and no HPC jobs were submitted.
They do not measure full-study speedup, peak memory, queue waiting time, or
validate a newly published container. The full TOM remains quadratic in event
count; removing an internal guard does not change that algorithm.
