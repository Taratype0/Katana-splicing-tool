# Native workflow and output directories

New studies use `nextflow run pipeline/main.nf`. The exported release uses
`nextflow run tools/katana-controller/pipeline/main.nf`. Parameters are identical.
Use `-profile`, `-c`, `-params-file`, `-work-dir` and `-resume` as Nextflow options;
sample, contrast, reference and output parameters use `--`.

One Nextflow DAG owns configuration, preflight, read processing, native models,
import, per-comparison downstream tasks, collection, report, project retention
and final delivery verification. Python and R programs run **inside tasks**;
no Python program polls one Nextflow session to start another.

## Output layout

| Directory | Contents |
| --- | --- |
| `rnaseq/qc` | FastQC and trimming reports, published per completed sample |
| `rnaseq/star` | Alignment diagnostics and splice junctions |
| `rnaseq/bam` | Sorted BAM/index when `storage_mode=retain` (default) |
| `rnaseq/salmon` | Per-sample quantification and native Salmon metadata |
| `rnaseq/matrices` | Transcript/gene expression matrices |
| `rnasplicing/rmats` | Native event tables and statistics by comparison |
| `rnasplicing/suppa` | Native PSI, differential and auxiliary SUPPA results |
| `differential/deg`, `differential/dtu`, `differential/dexseq_exon` | Native tables and fitted model files |
| `katana/tasks` | Individual downstream task products (`retain` mode) |
| `katana/project` | Independent project, inputs, completed artifact catalogue and collection |
| `report` | Approved seven-section HTML, ZIP, figures and full tables |
| `pipeline_info` | Trace, execution report, timeline, configuration and `delivery.receipt.json` |

An upstream-only run publishes `katana_bundle` as an interchange product;
complete runs retain its inputs in `katana/project` instead of publishing another
full bundle copy. Auxiliary tables go to `rnaseq/auxiliary` and reference
provenance to `pipeline_info/reference`. The report
format, algorithms, thresholds and native statistical tables are unchanged.
Successful intermediate tasks do not imply that the entire study is complete;
the final Nextflow delivery task requires both the accepted report and retained
client project. A failed required task prevents that final receipt.

## Existing results and recovery

```bash
# Existing upstream output: no alignment or upstream models are scheduled.
nextflow run pipeline/main.nf --mode downstream-only \
  --bundle /data/completed/katana_bundle -params-file analysis-resources.yaml \
  --outdir /data/study/results -work-dir /data/study/work -profile apptainer

# Imported project: runs the selected downstream graph.
nextflow run pipeline/main.nf --mode downstream-only \
  --project /data/project/project.json -params-file analysis-resources.yaml \
  --outdir /data/study/results -work-dir /data/study/work -profile apptainer
```

These commands still require the matching Katana container and resources for
enabled analyses, as described in the run guide. `design_factors` defaults to
`[condition]`; explicit multi-factor designs can be supplied in the parameter
file. Paired contrasts continue to use the declared subjects and comparison
design. Condition labels do not automatically create interaction models.

Use `-resume` for completed tasks from the same native run. Existing complete
downstream plans can be supplied through `--flow_plan` with `--project`. A new
project or new session is not automatically equivalent to a historical cache.
Keep the old software/path identity for an active historical managed run. Never
restart completed science merely to change directory names.

## Client use and work retention

After copying the result directory to another machine, rebind its local paths:

```bash
katana-flow verify-retained --directory results/katana/project
katana-open
```

Open `results/katana/project/project.json` in the browser. The retained project
includes target-analysis inputs; the report ZIP contains report assets, not BAMs.
Native Nextflow retains work for recovery. After successful delivery and retained
project readback you can clean it manually as below. Alternatively add
`--clean_work` when launching: Nextflow deletes its task work only after the
entire DAG, delivery verification and publication succeed. Failed work is kept.
Nextflow preserves task directories still referenced by earlier resumed runs;
clean those completed run names explicitly when they are no longer needed.
This replaces stage-by-stage Python cleanup; it does not delete a live producer
while downstream tasks still need it. After successful delivery and retained
project readback, `nextflow clean <run-name> -n` previews Nextflow-owned cache
cleanup; `nextflow clean <run-name> -f` applies it. Cleaning removes that run's
resume cache. Do not delete work while consumers or recovery jobs are active.
The older staged-cleanup managed command remains documented for existing runs.
