# MotifTrace

MotifTrace connects splice-event sequences, motif sites and compatible RBP
evidence. Its Python, command-line and interactive interfaces use the same
Katana workflow service, parameter definitions, input checks and result cache.
This integration retains the approved `1.0.0-frozen.20261017` scoring and
statistical methods. The package's interface and execution-policy revisions
are recorded separately in [the integration record](motiftrace/integration.v1.json).

## Install

Use Python 3.11 or newer in a virtual environment. From an accepted `main`
checkout:

```sh
python -m pip install ".[motiftrace]"
katana-motiftrace modes
```

For the browser application, install `".[motiftrace,next]"` instead.
The installation includes Katana's shared network dependencies; Qt is an
optional `desktop` extra. An independently minimized MotifTrace distribution is not yet
provided. There is no separate published `pip install motiftrace` package.

## Select an analysis

| Mode | Inputs and operation | Main outputs |
| --- | --- | --- |
| `boundary` | Saved splice-event selection, compatible genome FASTA/index and motif library/scanner; reference CLIP evidence follows its explicit configured policy. Computes event-neighborhood motifs and matched-background evidence. | Event/RBP evidence tables, motif and binding summaries, figures and methods. Motif and CLIP statistics remain separate. |
| `se_impact` | Selected SE events, compatible reference sequence and native FIMO environment. Runs the existing fixed-site sequence analysis. | Motif changes, matched controls, RNA sequence/site coordinates and result report. |
| `completed_maps` | A saved selection matching a verified completed SE evidence bundle; an additional matched-window registration is optional. Imports original results without refitting them. | Gene/event/motif browser, short RNA and motif-site tables, complete original test families and optional window browser. |
| `consequences` | Selected SE, A3SS, A5SS, MXE or RI events, compatible reference sequence, configured motif resource and matched background. Runs the existing local RNA reconstruction and consequence analysis. | Reference RNA forms and FASTA, gained/lost/retained sites, split genomic coordinates, paired-window/region results and browser. |
| `usage` | An existing supported comparison and the native expression, transcript/exon-usage, abundance and sequence inputs required by the selected operation; it can also receive a compatible completed case. | Usage-to-motif tables, sequence evidence, figures and report. Original DEG, DTU, DEXSeq and DAS endpoints retain their own statistics. |
| `parallel` | A verified completed MotifTrace/MAPP parallel-results manifest and its matching comparison. Imports fitted results. | Parallel-model browser and original method, motif, region, activity and native-statistics tables. |

Project inputs must declare their sample/comparison identities and compatible
references. Register resources and completed bundles in the project manifest;
use the interface's input panel or the project's existing configuration format.
Selection IDs refer to saved rows from a source artifact in that same project.
Do not reuse an ID from another project or change comparison direction by
renaming labels. See [project usage](USAGE.md) for project setup.

`boundary`, `se_impact`, `completed_maps` and `consequences` use saved event
selections. `usage` and `parallel` use the original comparison and their own
validated input contracts. There is no implicit quantification, differential
refit or general new-project MAPP fitting in these import/integration modes.

## Python

```python
from katana_motiftrace import MotifTrace, mode_catalog

print(mode_catalog())  # no project, server or GUI initialization

with MotifTrace("project/katana.project.json") as mt:
    # Inspect every parameter's current default and accepted choices.
    inventory = mt.inventory()
    print(inventory["modes"])

    # Use exact source row IDs from an existing eligible event table.
    selection = mt.select(
        source_artifact_id="SOURCE_ARTIFACT_ID",
        row_ids=["EXACT_EVENT_ROW_ID"],
        name="Events to inspect",
    )
    plan = mt.plan("consequences", selection_id=selection.selection_id)
    if not plan["ready"]:
        raise ValueError(plan["errors"])
    job = mt.run("consequences", selection_id=selection.selection_id)
    outputs = mt.artifacts("consequences", run_fingerprint=job.run_fingerprint)
```

An already-open `KatanaProject` exposes the same interface as
`project.motiftrace`; closing this borrowed client leaves the parent project
open. `inventory()`, `plan()` and `artifacts()` do not submit analysis jobs.
`select()` saves a selection, and `run()` explicitly executes the selected
existing workflow. `run()` returns the normal Katana job record and raises on
failure or timeout. Repeating the same validated request reuses its exact cache
when available; `force=True` is an explicit request to recompute.

## Command line

The module form `python -m katana_motiftrace` and the installed
`katana-motiftrace` command are equivalent. Place `--project` before the action.

```sh
katana-motiftrace modes
katana-motiftrace --project project/katana.project.json artifacts --mode completed_maps
katana-motiftrace --project project/katana.project.json select SOURCE_ARTIFACT_ID --rows selected-rows.json
katana-motiftrace --project project/katana.project.json plan consequences --selection-id SAVED_SELECTION_ID
katana-motiftrace --project project/katana.project.json run consequences --selection-id SAVED_SELECTION_ID
```

`selected-rows.json` is an array of exact source row IDs. Optional scientific
parameters are a JSON object passed with `--parameters parameters.json`; their
names/defaults are the same as the selected workflow's inventory and interactive
controls. A comparison-scoped mode accepts `--contrast EXACT_COMPARISON_ID`.
Use only comparisons supported by that mode; repeated flags do not turn a
single-comparison method into a pooled analysis.

For example, browse a registered completed parallel result:

```sh
katana-motiftrace --project project/katana.project.json plan parallel --contrast EXACT_COMPARISON_ID
katana-motiftrace --project project/katana.project.json run parallel --contrast EXACT_COMPARISON_ID
```

Commands emit JSON. Exit status is `0` on success, `2` for a plan that is not
ready (or a command-line syntax error), and `1` for an execution/input error.
An unready plan lists its errors and never submits work. Missing sequence,
scanner, database, original fit or matched-background evidence is reported;
it is not replaced with an invented result or silently downloaded resource.

## Interactive use

```sh
katana-open --mode browser --runtime-root ../katana-runtime
# Or, on a machine with a desktop session:
katana-open --mode desktop --runtime-root ../katana-runtime
```

Open the project, choose **Prioritize**, then select a MotifTrace mode in the
MotifTrace selector. Save events from an eligible result table when the mode
requires a selection. Inspect the input checks and parameter controls, review
the run plan, then explicitly run. Selecting a mode does not launch work.
Completed tables, figures and HTML browsers use the existing result viewers.
The selector, Python API and CLI share one mode catalog; advanced controls
continue to expose the underlying established workflows.

## Computational resources

MotifTrace does not reject valid sequence calculations because a fixed motif
count, RNA size, custom-library file size/row count, test count, paired-matrix
size or permutation-by-gene-pair count has been exceeded. All selected events,
requested motif models and eligible windows remain in their complete
statistical family. The whole-pair permutation default remains `999`, with an
integer minimum of `99` and no fixed upper count.

The older `max_tests` and `max_matrix_cells` fields are accepted for
configuration compatibility but no longer restrict a calculation. Omit them
from new parameter files. Changing their values does not allocate memory or
CPU resources.

On a cluster, request CPU, memory and wall time through the Nextflow executor
and site configuration. For local use, the operating system supplies the
available resources. Large requests can still exceed those resources; removing
an application limit does not guarantee that every calculation fits in memory.
Recorded dimension estimates describe the complete workload, rather than
serving as rejection thresholds. No events, models or windows are discarded to
make a calculation fit.

Motif-format validation, supported motif widths of 1–25 nucleotides, matching
rules and the minimum of 20 independent matched gene pairs for consequence
inference remain unchanged. These determine whether the input and statistical
analysis are valid; they are separate from computational capacity. Report
Top 10 presentation and full-table attachments also remain unchanged.

This execution policy applies to the updated source and packages or containers
built from that source. Previously built wheels, container images and running
jobs retain their original behavior until their runtime is replaced.

## Methods and interpretation

- Numerical defaults, scan scoring, sequence windows, matched controls and
  statistical procedures remain those of the approved scientific workflows.
  Effective parameters and original method identities are retained in run and
  artifact records. The API does not maintain a second set of thresholds.
- RBP identity/sets, neighborhood motif or reference CLIP evidence, and RNA
  sequence consequences are distinct evidence levels. Predicted sites and
  reference CLIP overlap do not establish binding in the current experiment.
- The frozen descriptive collection contains 5,712 models. Its 706 calibrated
  common literal motifs across 14 windows are a narrower scope; calibration
  must not be attributed to every descriptive model.
- Native P values, FDR families, missing/non-estimable results and method-specific
  comparison directions are retained. Viewing a smaller event selection does
  not redefine the original adjustment family or refit imported models.
- Historical completed results retain their original input/QC policy. Importing
  them does not certify compliance with a newer event-count filter.
- Sequence reconstruction of five event types is distinct from acceptance of
  a complete inferential analysis for each type. Current private regression
  evidence must be read at its stated scope; it is not a claim that every
  possible project, resource or external executable has been validated.

The package provides shared interfaces and distribution checks. It
does not change the frozen original application, recompute completed studies or
publish their private input files. The recorded release identity and adapter
boundary are in [the integration record](motiftrace/integration.v1.json).
