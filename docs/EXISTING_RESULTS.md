# Reusing completed results

## Generate one report from several completed sources

Use `katana-report` to assemble already completed results through the installed
package. It preserves the approved seven-section format and does not require a
previously generated report, the source checkout, or a shared artifact database.
This is report assembly: it submits no upstream, network, MotifTrace or MAPP
calculation. Missing analyses remain missing; see [compute only missing
analyses](#compute-only-missing-analyses) for a separate execution plan.

Install a matching wheel, standalone environment or downstream image using the
[installation guide](INSTALL.md#installed-package-report-commands). Keep its
wheel SHA, source revision and dependency/image identities together. No PyPI or
public container-registry installation is assumed.

A result-set descriptor names the target project and each completed source.
For example, arrange a workspace as follows:

```text
workspace/
  results.json
  study/katana.project.json
  sources/base/report-project.json
  sources/downstream-collection/collection.receipt.json
  sources/downstream-collection/...
```

The example below has two sources: existing differential/sequence results and
a completed downstream collection. The latter must contain the selected
comparisons' separate network results; it is not a pooled replacement.

```json
{
  "schema": "katana.completed_result_set/v1",
  "project_manifest": "study/katana.project.json",
  "mode": "bulk",
  "report_contrast_ids": ["B_vs_A", "C_vs_A"],
  "report_contrast_groups": [
    {
      "group_id": "selected_comparisons",
      "label": "Selected comparisons",
      "contrast_ids": ["B_vs_A", "C_vs_A"]
    }
  ],
  "report_native_figures_in_main": false,
  "report_require_plot_data": true,
  "sources": [
    {
      "source_id": "base_results",
      "kind": "descriptor",
      "path": "sources/base/report-project.json",
      "target_project_manifest_sha256": "<target-manifest-sha256>",
      "exclude_analysis_ids": ["programs.cosplice", "prioritize.mapp_native"]
    },
    {
      "source_id": "completed_downstream",
      "kind": "collection",
      "path": "sources/downstream-collection",
      "target_project_manifest_sha256": "<target-manifest-sha256>"
    }
  ]
}
```

Replace both SHA placeholders with the SHA256 of the exact target manifest
bytes. From `workspace`, this command prints that value without editing it:

```sh
python -c "from pathlib import Path; import hashlib; print(hashlib.sha256(Path('study/katana.project.json').read_bytes()).hexdigest())"
```

`source_id` values must be unique. Supported source kinds are:

| Kind | Input | Meaning |
| --- | --- | --- |
| `descriptor` | A report-source JSON naming its project, catalog and completed assets | Reuses the named artifacts and their original statistical scope. |
| `collection` | A completed Katana flow collection directory | Verifies its collection receipt and completed task/output identities. |
| `views` | A completed extra-views JSON | Adds explicitly bound figures/tables, such as previously completed native MAPP families; it does not establish that an absent calculation ran. |

The source list can contain a base descriptor, one network collection per
comparison and separate MAPP views. There is no fixed source count. Add rows for
the actual completed sources rather than combining experiments or rerunning
completed fits. `exclude_analysis_ids` removes superseded copies from that
source's report contribution; it does not delete the original files. Declare
overlapping source contributions explicitly; conflicting artifact identities
are rejected rather than selected by arrival order.

Assembly checks the ordered samples, biological units, conditions, batches,
sample attributes, complete contrast definitions/directions and genome build.
A different manifest revision requires the explicit exact
`target_project_manifest_sha256` binding shown above. A changed project ID also
requires `allow_project_rebind: true` on that source; this is presentation-only
rebinding and does not waive the biological-design checks. Native project
identities, file hashes and statistics remain in provenance.
This design comparison checks the genome-build declaration, not equality of
FASTA/GTF bytes. Resource equivalence needs the original plan/input fingerprints
or separately bound reference-file hashes. A target-manifest SHA authorizes the
stated presentation binding; it does not replace those resource checks.

For `views` containing native MAPP catalogs, declare `source_collection_root`
for the completed family directory. Optional `completion_evidence` entries have
`path` and `sha256` fields pointing to completed/verified evidence. Resource and
family identities still need their actual numerical/output records. An old
unavailable status is superseded only through an explicit
`superseded_statuses` rule naming its exact source/title and verified replacement
sources covering the required comparisons. The original status is retained in
assembly provenance; a new figure alone cannot erase a failed calculation.

Generate a compact, complete report directory without making an additional ZIP:

```sh
katana-report --project results.json --outdir exports/study-report \
  --compact --no-archive --result-store payload-store
```

The output directory must be new. The command prints JSON with the report,
manifest, provenance and archive location, if one was requested. Open
`exports/study-report/REPORT.html`; the comparison selector shows one declared
comparison at a time. The approved report layout and source statistics are
unchanged.

To deliver a ZIP at a separate location instead:

```sh
katana-report --project results.json --outdir exports/study-report-for-delivery \
  --compact --archive downloads/study-report.zip --result-store payload-store
```

`--no-archive` and `--archive` are mutually exclusive. Without either flag, the
archive is `REPORT.zip` in the new report directory. A custom archive location
is reported by the CLI; the complete report directory is retained in all cases.

### Moving completed inputs to another machine

Use relative paths where possible. If existing descriptors contain old absolute
paths, provide explicit consumer-side prefix bindings in `paths.json`:

```json
{
  "path_mappings": [
    {
      "source_prefix": "/previous-machine/study-results",
      "target_prefix": "copied-study-results"
    }
  ]
}
```

The `target_prefix` above resolves relative to `paths.json`. Prefix matching
respects path-component boundaries; it does not find resources by basename.
Use Windows prefixes such as `Z:/study-results` when that is the original
producer location. Do not change source manifests merely to replace paths.

```sh
katana-report --project results.json --path-map paths.json \
  --outdir exports/moved-study-report --compact --no-archive
```

The mappings affect this fresh report build only. They do not copy missing
inputs, change biological identities, rewrite original JSON files, or make a
historical raw-data project executable on another machine. This renderer needs
the completed evidence files it reads; it does not download BAMs by default.

### Storage and complete results

`--compact` consolidates byte-identical payload copies and losslessly encodes
complete table-browser data. It keeps every row, native P/FDR and original
source-table bytes, and records the transformations in `storage-audit.json`.
It does not discard large but unique matrices or scientific result families.

The optional `--result-store` is a separate, immutable payload store shared by
fresh report builds. Same-filesystem hardlinks can save disk space; unsupported
or cross-filesystem links fall back to copies. Treat stored payloads and their
linked report files as immutable. The report contains real files, rather than
links requiring the store on a recipient's machine. The store must be outside
source results and report output. It cannot be combined with `--from-report`.

For presentation-only refresh of an already portable report, use
[`--from-report`](REPORTS.md#selected-comparison-view-and-read-only-refresh)
instead. Both routes preserve the original report and completed analyses.

## Adopt completed scientific results for later tasks

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
