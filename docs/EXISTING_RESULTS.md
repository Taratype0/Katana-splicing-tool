# Reusing completed results

The existing-results entry verifies and copies previously completed Katana
results. It does not repeat their statistical models or label an old job as a
new calculation. Native tables, parameters, effect directions, P values, FDR,
sample vectors, scientific family manifests and completion receipts are kept.

First export a read-only artifact catalog and a report descriptor for the old
project. A descriptor identifies the original `project_manifest` and an
`artifacts` list or `artifact_catalog` JSON file. The import never writes the old
project database. All referenced result files must lie inside explicitly named
allowed roots.

```bash
katana-flow seal-existing \
  --descriptor /data/old-results/report-project.json \
  --allowed-root /data/old-results \
  --completion-evidence /data/old-results/native-completion.receipt.json \
  --family-root /data/old-results/motiftrace/comparison-A \
  --include-file /data/old-results/other-required-table.tsv \
  --out /work/new-run/historical-source.json

katana-flow adopt-existing \
  --source /work/new-run/historical-source.json \
  --outdir /work/new-run/adopted-results
```

Repeat `--allowed-root`, `--completion-evidence`, `--family-root` and
`--include-file` for additional explicit sources. Every member of each supplied
`COMPLETE.json` family is checked. Local browser attachments are copied with
their existing relative layout. Missing, altered, cross-project or unbounded
files fail the import. Unknown attempts and previous sealed snapshots are never
overwritten.

The adoption receipt has status `adopted`,
`scientific_computation_performed: false` and `native_jobs_submitted: 0`. Its
native evidence remains available for method-specific numerical and scientific
acceptance. An inventory or a historical report build receipt alone does not
prove that a missing analysis was performed.

## Compute only missing analyses

A new compute project may map input files and resources to a different machine.
An explicit study binding checks the original project ID, ordered biological
samples, conditions, batches, sample attributes, selected contrast arms and
genome build. It records the distinct original and compute-manifest hashes.
Annotation, sequence and database file contents must also be bound in the
subsequent computational tasks; equality of a genome-build name alone does not
prove resource identity.

```bash
katana-flow bind-existing \
  --adoption /work/new-run/adopted-results \
  --project /work/new-run/compute-project.json \
  --contrast-id comparison-A \
  --contrast-id comparison-B \
  --out /work/new-run/study-identity.json
```

A flow plan can use a transport task with `workflow_id:
"flow.historical.adopt"`, `adapter: "historical.adopt"` and inputs `adoption`
and `identity_projection`. Its dependencies are empty. Optional
`parameters.artifact_ids` selects exact old artifact IDs. Downstream tasks depend
on this transport task and retain the imported artifacts' original analysis IDs,
fingerprints and parameters. Completing the transport task does not claim that
the imported scientific workflows ran again.

Large human-readable TSV duplicates need not travel to a compute worker that
only consumes the original Parquet analysis input. This requires an explicit
analysis-inputs snapshot:

```bash
katana-flow seal-existing \
  --descriptor /data/old-results/analysis-input-descriptor.json \
  --allowed-root /data/old-results \
  --purpose analysis-inputs --exclude-exports \
  --out /work/new-run/analysis-input-source.json
```

This option retains the original export pointers as provenance and explicitly
does not claim to supply those exports or a portable report. A complete report
uses the full result snapshot. In particular, an old read-support policy is not
silently upgraded to the current policy: any new eligibility projection must be
recorded separately and must preserve the old flags and native statistics.

The portable adoption owns its result files and report views. Its original
project JSON remains byte-identical, including historical raw-read and reference
paths. Remapping those external computational resources requires the explicit
compute project and study binding above; adoption alone does not make the old
raw-data project runnable on another machine.
