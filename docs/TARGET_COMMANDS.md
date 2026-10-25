# Selected Gene / Event reports

Select genes or exact events in the browser, add them to **Gene / Event targets**,
and export a request ZIP. The ZIP contains selected completed evidence and
biological identities. BAMs, indexes, large references and execution-machine
paths stay outside it. No source project database is needed on the next machine.

This workflow creates independent target reports. It does not change the
seven-section study report or rerun differential analysis, network discovery or
MotifTrace. Native effects, P values, FDR and original testing populations remain
unchanged. A gene request includes all available measured events for that gene,
including non-significant events.

## Run locally or on a server

Install the development wheel with `pip install 'katana_splicing_tool-<version>-py3-none-any.whl[targets]'`.
Add the `next` extra for the browser service. Use the Python environment that
will perform the analysis; Windows, Linux and other compatible Python platforms
use the same commands. No SSH configuration, Dardel path or Slurm account is
part of a target request.

```sh
katana-targets inspect --request targets.zip
katana-targets plan --request targets.zip --bindings bindings.json
katana-targets run --request targets.zip --bindings bindings.json --outdir target-reports
katana-targets results --outdir target-reports
```

`python -m katana_next.targets` is equivalent to `katana-targets`. Standard
output is JSON; diagnostics go to standard error. Exit code 0 means successful
execution of the command, 1 means invalid inputs or an execution error, and 2
means the plan is not ready or requested outputs are partial/failed. Inspect
the returned per-target and per-module statuses before treating every requested
output as complete.

The separate binding file maps exported resource IDs to files accessible to
the execution environment. Use the IDs returned by `inspect`, not these example
IDs:

```json
{
  "schema": "katana.target_bindings/v1",
  "resources": {
    "resource-id-for-sample-bam": "reads/sample.bam",
    "resource-id-for-sample-index": "reads/sample.bam.bai",
    "resource-id-for-annotation": "reference/annotation.gtf"
  }
}
```

Relative paths resolve from the binding file. Change only these paths when
moving to another machine. Resource identities are checked against the exported
verification records. A sampled file fingerprint is explicitly distinguished
from a full-file SHA256. An absent or mismatched BAM is reported; it never
becomes an empty successful Sashimi track. Without BAM access, saved evidence
can still be rendered in cards and interactive tables.

Keep BAM indexes adjacent where possible. The built-in BAM reader also supports
explicit index bindings; other installed readers may require conventional
adjacent index names. Files are read where they reside rather than copied into
the report package.

## Browser workflow

1. Open an existing project and browse a completed result table.
2. Tick rows and explicitly choose **Add genes** or **Add events**. The basket
   persists across pages and retains the exact source table and browsing state.
3. In **Gene / Event targets**, export the request. Configure resource bindings
   for the execution machine and inspect the plan before starting the run.
4. Open the batch index, each target report or the complete ZIP. A returned
   result can be verified and opened through **Import results** on another
   Katana service that can access the unpacked result directory.

The default comparison scope is the union of comparisons represented by the
selected source rows; a source without comparison rows uses the declared study
comparisons. Choose **All project comparisons** or **Custom** to set an explicit
scope. An empty custom selection cannot mean all. The exported summary lists
the actual comparison IDs for every target.

`katana-open` starts the browser and Python service. On a server use
`katana-open --no-open` and access its loopback service through your normal
secure connection. An offline HTML report can display all included results;
running Python or reading BAMs requires the service or CLI.

## Restore returned results on another machine

To browse a result ZIP, extract it and open `INDEX.html`, or use **Import
results** in the matching project's browser workspace. Importing results opens
their verified reports; it does not submit new work.

To continue CLI work from the returned products, restore them into a new empty
directory. You do not need the original project database or exported request
ZIP:

```sh
katana-targets restore --results completed-targets.zip --outdir restored-targets
katana-targets inspect --request restored-targets/TARGET_REQUEST.zip
katana-targets plan --request restored-targets/TARGET_REQUEST.zip --bindings bindings.json --outdir restored-targets
katana-targets run --request restored-targets/TARGET_REQUEST.zip --bindings bindings.json --outdir restored-targets
katana-targets results --outdir restored-targets
```

`--results` also accepts an unpacked result directory or its
`result-manifest.json`. Optional `--project-id` requires an exact expected
project identity. Restore verifies the complete member hashes, target and
request identities, original comparison scopes and existing completion
records before creating a new local owner. It preserves the original result
manifest, native statistics and per-target provenance. `TARGET_REQUEST.zip`
contains the original request and evidence needed by the next plan/run;
`.target-restore.json` records the source manifest/archive SHA256 and the new
local restoration separately from the original scientific completion.

Restore performs no analysis, Sashimi read counting or figure rendering.
Rebind BAMs/indexes and references using the execution machine's paths, then
check the plan. Matching methods, evidence and resources reuse the completed
targets; restoring a result alone does not certify a different software or
resource identity as a cache hit. An occupied output directory, unknown writer
claims, unsafe archive, altered product or mismatched project is rejected;
existing products and claims are never removed to make a restore succeed.

## Batch and scheduler use

```sh
katana-targets tasks --request targets.zip --format tsv
katana-targets run --request targets.zip --bindings bindings.json --outdir one-target --target-id gene-EXACT_ID
```

Use actual exported target IDs. Each scheduler task can run one target or an
explicit subset in its own output directory. The same command works inside a
Slurm, PBS, Kubernetes, Nextflow or other scheduler-managed environment; the
allocation, resource limits and scheduling policy belong to that environment.
Katana does not submit to a cluster automatically or impose a four-task limit.
For execution within one allocation, `--workers N` sets the local worker budget;
the default derives from the available allocation. Plotting operations respect
the renderer's concurrency constraints.

Browser jobs currently run within the Python service process. Keep it running
until they finish. After an abnormal service exit, persisted active records and
unknown writer claims require inspection; the browser does not automatically
declare them stale or resubmit them. The CLI also preserves unknown claims.

Rerunning the same request with the same verified products reuses completed
targets. Cache identities include source evidence, resources and rendering
code; result files are checked before reuse. Unknown live writer claims and
unrelated output directories are not overwritten. A failed target retains its
diagnostic record and can be retried independently.

## Python interface and direct export

```python
from katana_next.targets import export_targets, plan_targets, run_targets, inspect_results, restore_results

# context is an existing Katana ProjectContext, opened by the caller.
exported = export_targets(context, [{
    "source_artifact_id": source_artifact_id,
    "row_ids": selected_row_ids,
    "target_kind": "event",  # or "gene"
    "filters": {"event_set": "stable_significant", "fdr": 0.05}
}], "new-target-request")

plan = plan_targets(exported["archive_path"], resources="bindings.json")
result = run_targets(exported["archive_path"], "target-reports", resources="bindings.json")
verified = inspect_results("target-reports")

# A returned result package can become a new resumable local workspace.
restored = restore_results("completed-targets.zip", "restored-targets")
```

For scripted export, write the selection array above to `selections.json`, then
run `katana-targets export --project katana.project.json --selections selections.json --outdir new-target-request`.
Use `--module gene_card --module event_card` when read analysis is not requested.

## Outputs

The batch includes `INDEX.html`, `REPORT.zip` and a SHA-covered result manifest.
Each target has its own `REPORT.html`, machine-readable card and provenance,
complete selected evidence in Parquet/TSV/JSON, searchable paginated offline
tables, saved sample-level plots and native Sashimi figures/counts where
resources are available. Saved network and MotifTrace evidence remains tied to
its original analysis; selecting a target does not create a new statistical fit.

The internal ledger and content-derived request ID detect accidental changes.
They are not a publisher signature. Distribute the ZIP SHA256 through a trusted
channel when external provenance authentication is required.
