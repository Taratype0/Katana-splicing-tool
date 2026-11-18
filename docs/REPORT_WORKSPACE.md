# Current unified local entry

Use [the unified local guide](LOCAL_UNIFIED.md) for report ZIP/folder import, full-population exploration, local enrichment and target cards (port 18775). The notes below describe the preserved 8 November base acceptance.

# Report exploration and target pipeline follow-up

This version adds a separate local report workspace. Start with
`Start-Report-Workspace.ps1` (port 18774), then **Open extracted report**.
Supply the extracted directory containing `report-manifest.json` and a new,
empty workspace directory. The source report and sealed release are preserved.

The importer supports `katana.approved_report/v1` with recorded display-project
sample and comparison metadata. It imports complete TSV/CSV/Parquet tables,
including compressed tables and nonsignificant rows. It verifies recorded hashes
where present, records explicit missing-hash status, preserves source metadata,
and writes its own Parquet copies and import audit. A single HTML file without
its data and design metadata is insufficient.

In **Report explorer**, choose a complete table, comparison, identifier list,
FDR and effect columns/thresholds. Pagination covers the full table. The
interactive effect/FDR plot explicitly represents the current page; it does not
claim to show the whole population. Missing values are not zeroes. Select rows
across pages or paste exact gene IDs/symbols/event IDs to filter, then add genes
or events to the persistent target basket. Original report figures remain
available via the seven-section report link. Browse thresholds do not rerun or
relabel enrichment, network discovery or original statistical families.

## New follow-up task

In **Gene / Event workspace**, export the selected request, then choose
**Create Nextflow target task**. The package is a new independent task with the
same `KATANA_TARGETS` module used by the main Katana pipeline:

```sh
nextflow run pipeline/main.nf -entry TARGETS \
  --target_request request.zip --target_bindings bindings.json \
  --outdir target-results -work-dir target-work -c site.config
```

The downloaded task has its own `main.nf`, module, config, request ZIP, logical
resource requirements and binding template. Fill the template on the execution
machine and save as `bindings.json`; BAMs remain on that machine. Bindings are
native Nextflow `path` inputs and are staged for executor/container visibility.
Configure a container containing this Katana revision or provide this Python
package in the task environment. Site executor/resources belong in `site.config`.
Resume with `-resume`. This entry never launches FASTQ/alignment, differential
models or whole-study network/MotifTrace fits.

The browser marks the task **prepared, not submitted**. It does not claim a
remote scheduler submission or completion. Each target publishes its own result
directory under `results/targets/<target_id>`, with cards, Sashimi when available,
tables, report ZIP and verified result manifest. Partial read evidence remains
partial in the manifest even when the diagnostic-producing process finishes.
Returned target directories can be verified and reopened through the existing
target-results import. Local Python execution remains an explicit alternate
execution option, not a scheduler adapter.

## Validation boundary

Windows tests exercise report import, original-file preservation, corruption
rejection, complete table querying, exact target export, resource binding, real
generated indexed BAM drawing, target result verification and task packaging.
Browser acceptance covers report import, batch selection, target export and
Nextflow task download. Actual PTBP1 report acceptance imported 44 complete
tables containing 1,690,979 rows. Those rows include different evidence types;
this count is not a count of unique genes/events.

Native Nextflow 25.10.4 execution on the Linux local executor passed on
8 November 2026 ([acceptance run](https://github.com/Taratype0/Katana-splicing-tool/actions/runs/42395577378)).
A selected gene expands into four measured events; a second exact event retains
its nonsignificant FDR. Three native processes complete and all three are cached
on resume, with identical target product bytes. Removing a BAM index produces
explicit partial/unavailable records. Rebinding it reruns the two detail tasks
and recovers complete products. Returned result ZIP restoration verifies the
original request and product ledger. The original report remains unchanged.

This verifies software behavior on generated indexed BAMs, not biological
validity of a study or Slurm scheduling. The generated BAM counts and saved
statistical fixture values are independent test inputs. No formal dataset or
cluster job was launched. The sealed port-18773 default remains preserved.

## Interface

Imported-report navigation shows only Results and the Gene / Event workspace.
Result selection, compact filters and table actions occupy separate rows; gene,
event and statistical columns precede technical coordinates. Original statistics
and checksum details remain under Source. Target choices and outputs sit beside
the pipeline task and returned-result import. Local execution and command details
are expandable. Missing-data, verification and submission states remain explicit.

The interface uses the one-click report palette: ink #1e2630, blue #2863a5,
cyan #2d7e9b and white. Browse buttons select report folders, project manifests,
local read/reference resources and returned results through a local file picker.
The picker lists Windows drive letters without probing disconnected drives.
Report import automatically creates a separate workspace under the user's
Katana/workspaces data directory. Reopening the same report in the same runtime
reuses that workspace after checking the original and imported evidence hashes.
Returned result ZIPs use the existing verified restore contract; directories
and result manifests remain supported. Neither choice starts an analysis.

Click-through acceptance on 8 November 2026 selected an extracted synthetic
report, one exact event (SE:e5), four BAM/index pairs and the annotation without
typing paths. A local Python target run generated seven displayed plots;
the downloaded ZIP was selected, restored and verified in the interface.
Nextflow task preparation was also exercised. This local UI acceptance did not
submit a Nextflow/Slurm job or validate biological conclusions.
