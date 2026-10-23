# Target workflow implementation contract

Internal interface agreement for the portable target workflow. Public user
behaviour is specified in `../TARGET_REPORTS.md`. No scheduler, hostname,
reference species or workstation path is built into these interfaces.

## Package and ownership

New code lives in `katana_next.targets`. It must not modify the seven-section
renderer or frozen numerical/Sashimi algorithms. Importing this package does
not open a project, service, GUI or job.

## Export API: `bundle.py`

```python
export_targets(context, selections, output_dir, *, modules=None,
               contrast_ids=None) -> dict
open_bundle(path)  # context manager yielding Bundle(request: dict, root: Path)
validate_request(request, root, *, verify_files=True) -> dict
```

`context` is the existing open ProjectContext. Each selection is a dictionary
with `source_artifact_id`, `row_ids`, `target_kind` (`gene` or `event`), and
`filters` (the browsing snapshot). Validate source completion, current project
identity and exact selected rows; do not resolve the latest unrelated table.
Multiple selections from different source tables are accepted in one request.
`modules=None` means `gene_card`, `event_card`, `sashimi`.

Export response has `request_id`, `request_path`, `bundle_directory`,
`archive_path`, `target_count`, `targets` and `resources`. Return absolute local
output paths in this response, but never embed machine paths in the portable
request. The directory and ZIP include `request.json` and SHA-covered evidence.

Portable request schema: `katana.target_request/v1`.

```text
schema, request_id, created_at
project:
  project_id, name, assay, source_manifest_sha256, source_snapshot_fingerprint
  references: genome_build, coordinate_system, gtf_resource_id (if available)
  samples[]: sample_id, biological_sample_id, condition, attributes,
             bam_resource_id / bai_resource_id (if available)
  contrasts[]: contrast_id, label, numerator{label,sample_ids}, denominator{...}
targets[]:
  target_id, target_kind, gene_id, gene_symbol
  event_id/event_type/event_geometry (event targets)
  comparison_ids[], source_selections[]
evidence[]:
  evidence_id, analysis_id, name, entity_type, path (bundle-relative)
  sha256, bytes, source_sha256, source_fingerprint
  parameters, metadata (portable), target_ids[]
resources[]:
  resource_id, kind (bam/bai/gtf/...), sample_id when applicable
  sha256 and bytes when available, otherwise explicit verification metadata
modules[], members{relative_path:{sha256,bytes}}, execution
```

Gene targets expand from complete measured event evidence, not just significant
rows or three high-variance events. Event targets keep exact identities and
geometry. Source table columns, dtypes, effects and statistics are retained.
Export only matching genes/events and necessary context (samples/module scores,
etc.); do not export unrelated cohorts, BAMs or large references. Members and
ZIP paths are strictly relative, bounded, unique, non-symlink regular files.
Opening a ZIP verifies its member set/content in a fresh temporary directory;
opening a directory/JSON verifies the same ledger. Reject tampering, foreign
sample/contrast identities and unsafe archives before calculation.

Bindings are independent JSON: `katana.target_bindings/v1`,
`resources: {logical_resource_id: local_path}`. Relative bindings resolve from
the bindings file. Source-export paths may be returned separately as a local
bindings template; they are not part of `request.json` or its portable ZIP.

## Execution API: `runner.py`

```python
plan_targets(request_path, resources=None, *, outdir=None) -> dict
run_targets(request_path, output_dir, resources=None, *, workers=None,
            target_ids=None, progress=None) -> dict
inspect_results(output_dir) -> dict
```

`resources` is a bindings JSON path or dictionary. `target_ids` is an explicit
subset, not a hidden truncation. `progress` may be a callable accepting one
JSON-safe event dictionary. Plan reports `ready`, `errors`, `warnings`,
`tasks`, resource requirements and reusable products. No job or source write
occurs in planning.

Run produces `INDEX.html`, per-target `REPORT.html`, cards, complete selected
tables and interoperable JSON/Parquet data, a result manifest and `REPORT.zip`.
Return `schema`, `status` (`complete`, `partial`, `failed` or `cache_hit`),
`request_id`, `report`, `archive`, `manifest`, `targets` and cache summary.
Original P/FDR and source test-family identity remain intact. No DESeq2,
rMATS, network fit or MotifTrace rescan is implicitly started.

Use the existing NativeSashimiAdapter, preserve its exact event geometry and
biological replicate ordering, and record the actual reader. Missing BAM/index
or unsupported resources are explicit per-target failures/unavailability, never
zero-read results. Verify bindings before invoking the adapter. Resource checks
must distinguish a full SHA from a sampled fingerprint.

Every target has an independent input/method/product identity and verified
completion manifest. Repeated matching completed targets reuse verified files
without changing their mtimes. Partial targets retain diagnostics and can be
retried independently without replacing completed ones. Protect concurrent
writers with explicit owned claims; never delete unknown live claims. A batch
retains all target statuses even if one fails.

No arbitrary four-job cap. Local worker choice is explicit or derived from the
available allocation. The same CLI must run inside any scheduler allocation;
an optional Nextflow example can split target IDs using normal executor config.
Do not advertise a local thread pool as generic scheduler submission.

## Browser/API

A cross-page target basket stores exact selected row IDs plus browsing context.
Add gene/event choices, remove/clear, module choices, export ZIP and run controls.
Use existing page selection and keep the seven-section report files untouched.
Frontend saving must retain actual source parameters and current browse settings
separately. Imported result manifests can be reopened for table/figure browsing.

Use a new router `katana_next.api.targets` registered in `create_app`:

```text
POST /api/projects/{project_id}/targets/export
  {selections, modules, contrast_ids} -> export response + download_url
POST /api/projects/{project_id}/targets/run
  {selections, modules, contrast_ids, resources?} -> explicit local job/status
GET  /api/projects/{project_id}/targets/exports/{request_id}/download
GET  /api/projects/{project_id}/targets/runs/{run_id}
GET  /api/projects/{project_id}/targets/runs/{run_id}/files/{path}
POST /api/projects/{project_id}/targets/results/import
  {manifest_path} -> verified result summary and browse URL
```

Only serve files bounded by the owned export/result root. Server execution is
ordinary Python on the selected machine, not tied to Dardel. CLI portability
does not require that the destination has the original artifact database.

## CLI and acceptance

Root integrates `katana-targets` with export/inspect/plan/run/tasks commands,
distribution metadata, user documentation and installed-package acceptance.
Keep examples generic. Test changed machine roots and moved bundles with source
database removed; gene and exact-event identity; thresholds/reference/contrast
scope; corrupted/missing data; resource binding; real small-BAM Sashimi;
per-target failure/cache/concurrency; browser selections/download/import; and
unchanged original inputs/approved report files. CI must exercise Linux and
Windows independently; do not claim untested systems are verified.
