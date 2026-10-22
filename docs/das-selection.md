# DAS and RS-DAS selection

The default differential-splicing foreground is **RS-DAS**. It applies the existing caller significance/effect criteria and additionally requires IJC + SJC >= 10 in **every declared sample of both comparison arms**. Ordinary **DAS** retains the significance/effect criteria without that read-support filter. Caller effects, P values and FDR are preserved. Missing read counts do not pass support.

Use the workflow inventory to check each module's `das_selection` capability. A `switchable` module exposes `event_set`; a module marked `rs_required` rejects ordinary DAS. `measurement_population` modules use their own measurable PSI population, including non-significant events, and do not use either significant foreground for network construction. Imported completed MotifTrace/MAPP/usage collections preserve their own source scope; this switch does not refit or relabel them.

## Python

```python
from katana_next.client import KatanaProject

with KatanaProject("project.json") as project:
    # Read-only: normalized parameters, required parents and missing inputs.
    plan = project.plan("differential.as", event_set="significant")
    # Explicit downstream execution using the same selection.
    result = project.run("differential.as", event_set="significant")
    # Reuse only matching saved downstream artifacts; no hidden analysis.
    report = project.run("reports.bundle", event_set="significant")
```

Omit `event_set` for the RS-DAS default. Supplying conflicting values through a parameter dictionary and keyword argument raises an error.

## Command line and HTTP

```bash
katana-next --project project.json workflows
katana-next --project project.json plan differential.as --event-set significant
katana-next --project project.json run differential.as --event-set significant
katana-next --project project.json run reports.bundle --event-set stable_significant
```

The same setting is accepted in a JSON parameter file as `{"event_set":"significant"}`. HTTP plan/run requests use `{"parameters":{"event_set":"significant"},"force":false}`. The public `/api/workflows` inventory includes the capability and choices used by the desktop/web application.

## Desktop and reports

The desktop/web workflow control labels the two choices as ordinary DAS and RS-DAS. Changing the selection keeps it consistent across the other modules explicitly declaring support. A change only checks the matching cache until the user runs the module. Existing jobs, source event flags, result files and parameters are not rewritten. Restoring a saved result's settings restores its recorded selection.

The evidence-package report (`reports.bundle`) packages only completed switchable results with the requested, explicitly recorded policy. A missing matching result is **not computed for this view**; an older or RS-only result is not shown as an ordinary-DAS result. RS-required analyses remain separate and are labelled. Measurable-population networks keep their original input population. The report records the selection in HTML, PDF, inventory and JSON provenance, and ranks DAS Top 10 separately within each comparison using FDR, absolute effect and event ID. Complete source tables are retained.

The older completed-study-report import (`reports.interpreted_bundle`) packages an existing report with its original scope. It does not provide a foreground switch over already rendered figures or tables. The seven-section dummy report is a separate report-design preview; it is not a claim that the full new production report generator has been released.

No FASTQ/BAM analysis is scheduled by selecting an event set or exporting a report. A supported downstream module may need an explicitly requested new downstream result when that selection has not previously been calculated.
