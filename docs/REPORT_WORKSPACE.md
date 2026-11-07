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

Native Nextflow execution on Linux/Slurm remains to be verified. The local
machine has no usable Linux executor, and the configured SSH host was
unreachable during this acceptance. No formal dataset or cluster job was
launched. This is a separately launched development candidate; it does not
replace the sealed port-18773 default or claim nf-core acceptance/certification.
