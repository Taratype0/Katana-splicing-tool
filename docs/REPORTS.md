# Study reports

Katana assembles completed results into the approved seven-section report:
design and QC; AS landscape; differential molecular results; cross-comparison
results; GO/KEGG and pathways; co-splicing networks; RBP and MotifTrace.

## Command line and Python

Install the development package with report dependencies:

```sh
python -m pip install ".[report]"
katana-report --project study/katana.project.json --outdir exports/study-report
```

The output directory must be new. Report generation reads the manifest and its
completed artifact export, or opens the artifact database in read-only mode.
It does not migrate the original project, run upstream tools, or fit missing
analyses. An explicit artifact export can be supplied with `--catalog artifacts.json`.

```python
from katana_next.reports import build_report

manifest = build_report("study/katana.project.json", "exports/study-report")
```

The same renderer is registered as `reports.study` in the application's report
workspace and shared workflow API. Its outputs are the HTML report, portable
ZIP, provenance and machine-readable report manifest. The existing
`reports.bundle` and completed-report import routes remain available.

### Selected-comparison view and read-only refresh

The main report opens on the first declared comparison. Selecting another
comparison changes the visible differential, enrichment, network and motif
results together; shared study QC remains visible. Cross-comparison panels
retain their declared axes and appear only for the selected comparison's group.
The seven chapters and blue/red palette are retained.

The network chapter embeds the module browser directly. Module genes are
complete, searchable and downloadable, with symbols, gene IDs and all event
identities. The dendrogram remains in outputs rather than the main report.
Where original fitted parameters and numeric inputs are authenticated, a graph
view can re-export edges among the 100 highest recorded absolute-kME events per
assigned module using the unchanged original adjacency kernel. The graph is an
explicit display subset (up to 1,000 edges); complete module membership is not
filtered. This does not rerun module discovery, TOM, association or enrichment.

MotifTrace sequence cards show at most ten distinct models with a recorded
native test FDR at or below 0.05. If none qualify, the main report says so instead
of filling ten nonsignificant cards. All models and original statistics remain
in outputs. The selected test-family FDR is not a new model-level FDR; MAPP's
native Bonferroni results remain a separate method.

Refresh an existing portable report without rerunning its analyses:

```sh
katana-report --from-report exports/previous-report --outdir exports/updated-report --compact
```

Older reports without embedded display-project metadata additionally require
`--project study/katana.project.json`. Optional `--network-display-parameters`
accepts a JSON mapping of comparison IDs to their original fitted parameter
files. File paths in this mapping resolve relative to the mapping file.
Unsupported or unavailable network parameters never trigger a replacement fit.

`--compact` consolidates byte-identical result copies and losslessly compresses
complete table-browser payloads, retaining every row and the native source
tables. It writes a storage audit and updated portable-file hashes. The report
must be refreshed into a new directory; prior reports and scientific files are
preserved. Current browsers decode these self-contained tables offline.

### Relative file paths

A report descriptor can live separately from the original study. Its
`project_manifest`, catalog fields, `extra_assets` and `supplemental_sources`
are relative to the descriptor file's directory. Inline descriptors supplied
as Python dictionaries use the caller's current working directory. Native
artifact paths inside a catalog remain relative to the native project
manifest, regardless of where the descriptor or exported catalog is stored.
The same native directory resolves `input_root`, `output_root`, annotation
GTF/FASTA paths, project asset paths and sample `quant_sf` paths. Resolution
uses a report-only copy; parent hashes always describe the original manifest.
Remote URIs are preserved as declared and are not fetched by this local report
loader; use locally available completed inputs for report assembly.

An explicitly supplied `--catalog` or Python `catalog=` argument is relative
to the caller's current working directory. This overrides a descriptor's
catalog path. For `extra_views`, a referenced JSON file resolves its own panel,
table and Top 10 complete-result paths relative to that views file; inline
view definitions resolve those paths relative to the report descriptor.
Absolute paths retain their original targets. Report loading does not rewrite
any descriptor, catalog or source result.

### Select comparisons for the report

Select the experimental questions to show instead of displaying every available
pairwise comparison. A descriptor can name existing contrasts and group related
questions for the cross-comparison figures:

```json
{
  "project_manifest": "study/katana.project.json",
  "report_contrast_ids": ["KO_vs_control", "WT_reexpression_vs_KO", "MUT_reexpression_vs_KO"],
  "report_contrast_groups": [
    {
      "group_id": "untreated",
      "label": "Untreated comparisons",
      "contrast_ids": ["KO_vs_control", "WT_reexpression_vs_KO", "MUT_reexpression_vs_KO"]
    }
  ]
}
```

IDs must refer to the original manifest. Each selected comparison belongs to
exactly one group, with up to four comparisons per group. The report records
both the selected and original contrast lists. The original manifest and all
completed source results are retained; selecting a subset does not rerun tests
or change their original P values or correction families.

Each group's intersection chart calculates exclusive memberships within that
group. It shows every nonempty pattern, together with set sizes and complete
membership tables, rather than truncating patterns from a larger chart. Event
identity, exon-bin identity and gene-level intersections remain separate.

## DAS selection and Top 10

RS-DAS is the default: native caller FDR and effect criteria plus IJC + SJC of
at least 10 in every declared biological sample in both contrast arms. No
direction-consistency or uncertainty cutoff is added. Historical RS flags are
retained in source tables; they do not establish this newer read-support rule.

The report-wide selector, differential section and DAS Top 10 share the same
selection. Ordinary DAS retains native significance without the count gate.
The initial selection can also be set with `--event-set significant`.

Top 10 ranks native FDR ascending, absolute native effect descending when
available, then stable object identity. Gene, event, transcript and exon units
remain separate. Missing usage effect sizes are not invented. Every entry has
expandable source details and a complete-result link.

Matched downstream views follow the selected event set. Fixed completed
imports retain their original measurement and test families; a missing matched
analysis is not silently replaced with a differently selected result. Networks
use their declared measurable PSI population, rather than significant DAS alone.

## Figures and complete results

The report includes summaries and per-comparison views. Complete tables,
searchable table pages, the all-event PSI matrix browser and figure galleries
are stored under `outputs/`. The PSI browser preserves missing measurements,
shows exact values on hover and offers native or row-standardized display.
Similarity matrices remain in the output gallery.

Registered network and MotifTrace figures, tables and local browsers retain
their comparison, module, method and model identities. Local dependencies are
copied into the portable package; original files remain unchanged. Unrun CLIP,
MAPP, external preservation or other optional analyses are identified explicitly.

For a verified completed MotifTrace family, the output files include every
member of its scientific completion manifest, including supporting RNA
sequences, matching records and per-replicate opportunity tables. Files that
are not listed in that manifest are not collected as scientific results.

QC can be rendered from completed FastQC, Trim Galore and STAR outputs. Absent
5-prime/3-prime or coverage measurements are recorded as unavailable; no new
read processing is started to fill a report panel.

Each main figure can carry explicit `data_links` supplied by its producer. A
link records `role`, `label`, `path`, `contrast_id`, `policy` and `unit`;
supported roles include `plot_data`, `full_results`, `foreground`, `background`
and `methods`. Paths in an external views JSON resolve relative to that file.
The renderer validates file contents against their declared comparison and
measurement scope, copies the data and records source and portable hashes.
`report_require_plot_data: true` rejects an unbound main figure. A whole-study
master table may be linked as complete results, but cannot stand in for the
actual plotted points.

`report_native_figures_in_main: false` keeps historical registered figures in
the output gallery while displaying the newly bound report figures in the
body. Complete tables stay outside the main report. Fixed RS-required
MotifTrace outputs remain visible with their original scope explicitly stated
when the ordinary-DAS selector changes; they are never relabelled as an
ordinary-DAS calculation.

## Provenance

`provenance.json` records source hashes, actual workflow parameters, filtering
stage, references and measurement scopes. Source P/FDR and effects are
preserved. An existing study filtered only after model fitting is described as
post-fit filtering; report assembly does not claim that upstream tests used a
new exclusion policy. New downstream computations have their own source and
completion records.

Reports contain measured results and stated methods. They do not infer an
experimental purpose, rescue mechanism, or validation programme. Selected
gene/event dossiers remain separate interactive analyses.

The complete content and visual requirements are in
[the approved report specification](REPORT_SPECIFICATION.md).
