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

QC can be rendered from completed FastQC, Trim Galore and STAR outputs. Absent
5-prime/3-prime or coverage measurements are recorded as unavailable; no new
read processing is started to fill a report panel.

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
