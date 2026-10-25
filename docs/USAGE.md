# Using Katana

Install the accepted version from `main` using the [installation guide](INSTALL.md).
Start with the command matching your input:

| Starting point | Entry |
| --- | --- |
| Bulk FASTQ or genome-aligned BAM | Nextflow `--mode complete`: upstream, selected Katana calculations and the seven-section report. |
| A native bundle or project with downstream calculations still needed | `katana-flow import` if needed, then Nextflow `--mode downstream-only`. |
| Completed analysis results | `katana-report`: assemble the fixed report without refitting analyses. |
| Interactive result inspection and target selection | `katana-open`: start the browser and its Python service. |
| Selected genes or exact events | `katana-targets`: generate independent cards, requested Sashimi and portable target reports. |

The [entry-point guide](ENTRY_POINTS.md) gives the complete FASTQ/BAM,
existing-result, reporting and target commands. Use only the comparisons you
want to analyze. Raw-read execution is bulk RNA-seq; processed single-cell and
single-nucleus inputs retain their own supported downstream scope.

## Projects and inputs

Open an existing v2 Katana project manifest in the interface. It records samples, contrasts, input assets, reference annotations and output locations. Supported legacy projects can be migrated into a separate output location:

```sh
python -m katana_next.migration legacy-project.yaml --output-root ../analysis-results --manifest ../analysis-project/katana.project.json
```

Replace the example file names with your own paths. Migration discovers compatible files; it does not make missing references or analyses available. Check the resulting manifest and readiness report before running analyses.

Use explicit sample IDs and condition labels. Preserve donor/biological replicate identities and distinguish them from technical measurements. Inspect numerator and denominator before interpreting an effect. Select the contrasts relevant to your question rather than treating every possible comparison as a required analysis.

## Workflow execution

The [threshold reference](THRESHOLDS.md) records the active `bulk-all-samples-counts-v1` policy: AS read support requires `IJC + SJC >= 10` in every expected sample of both arms of the selected contrast. The count threshold is fixed; there is no selectable replicate-QC preset, direction-consistency gate or approximate-SE option. Missing counts or ambiguous sample mappings do not pass. Significance selection is separate.

The browser and Python interfaces use the same workflow inventory and artifact services. A workflow's readiness plan describes its inputs, parameters and limitations. Opening a result or a report is not an instruction to recompute its upstream analyses. Compatible completed artifacts can be reused; changing an input, reference or scientific parameter requires the corresponding new analysis identity.

Inspect the inventory, plan an analysis and list available artifacts:

```sh
python -m katana_next.cli --project ../analysis-project/katana.project.json workflows
python -m katana_next.cli --project ../analysis-project/katana.project.json plan reports.bundle
python -m katana_next.cli --project ../analysis-project/katana.project.json artifacts
```

For workflows with parameters, supply a JSON object using the same field names as the interface:

```sh
python -m katana_next.cli --project ../analysis-project/katana.project.json plan WORKFLOW_ID --parameters parameters.json
python -m katana_next.cli --project ../analysis-project/katana.project.json run WORKFLOW_ID --parameters parameters.json
```

Planning checks readiness. `run` executes the selected workflow. Neither command implies a complete raw-data pipeline.

### DEG thresholds at runtime

The DEG analysis defaults to an adjusted P-value cutoff of `0.05` and a minimum absolute log2 fold change of `1.0`. Both can be changed for an individual run:

| Parameter | Meaning | Default |
| --- | --- | --- |
| `fdr` | Maximum imported adjusted P-value, inclusive | `0.05` |
| `min_abs_log2fc` | Minimum absolute log2 fold change, inclusive | `1.0` |

In the desktop or browser interface, open **Differential gene expression** and edit **DEG FDR** and **Minimum \|log2FC\|** before running. For the command line, copy and edit [the DEG parameter example](../examples/parameters/deg.json), then supply it:

```sh
python -m katana_next.cli --project path/to/katana.project.json run differential.deg --parameters examples/parameters/deg.json
```

The JSON example contains:

```json
{
  "fdr": 0.05,
  "min_abs_log2fc": 1.0
}
```

The same keys are accepted by Python's `project.run("differential.deg", parameters)` and the API's `parameters` object. FDR must be between 0 and 1; the effect-size cutoff must be nonnegative. An effect-size cutoff of `0` removes the additional fold-change requirement for rows with a finite effect estimate.

These settings classify imported DEG statistics; they do not refit DESeq2 or overwrite its P-values. The normalized parameters are retained with the run and its results and participate in the cache identity. The candidate-evidence analysis has its own declared thresholds, and interpreted reports read their configured completed delivery; changing a DEG run does not silently update those analyses or an existing report.

## AS read support

Select the contrast to check, then use **AS read support**. The workflow checks all participating samples individually; it does not accept a passing group mean or a subset of supported samples. Each declared arm must be present, count arrays must map unambiguously to the expected sample IDs, and every sample needs finite, nonnegative IJC and SJC with a sum of at least ten. PSI must be finite and within 0 to 1 as an input-validity requirement; this is not a replicate-stability threshold. Invalid or unavailable inputs are blocked or carry an explicit reason; empty or missing inputs never pass.

The compatibility output names `stable_event` and `stable_significant` refer to the current count gate and its intersection with native significance. They do not mean that replicate-direction consistency, uncertainty calibration or biological stability was tested. Read-support results and significance results should be inspected separately. Historical outputs keep their original policy identity.

For a multi-contrast Landscape, each selected sample must have a valid PSI and meet the count rule, and every included contrast row for an event must pass read support. An event passing one contrast is not sufficient to import its failing measurements from another contrast. A whole-cohort network applies its own versioned all-sample policy to the entire fitted cohort, rather than only one contrast. Upstream rMATS statistics, gene-level differential models and completed results are not rewritten.

## MotifTrace and networks

MotifTrace exposes event-level sequence consequences, motif and RBP evidence, and links back to their source inputs. Choose the appropriate evidence layer: RBP identity annotations, boundary motif/reference binding evidence, and mature-RNA sequence consequences answer different questions. Native method results and valid imported results retain their source identities. The presence of a resource does not establish that every optional author-native fit has run.

For count-bearing network discovery, declare every biological sample in the fitted network and use the `network-all-samples-counts-v1` input policy. In the integrated pipeline, each selected comparison gets an independent network using only its own samples. Every participating sample must have valid PSI and IJC+SJC >=10 for an event to enter; events with a constant PSI vector are excluded. Do not use significant DAS events as the default input universe. Count and sample mapping must be available; a fraction of supported samples, a group average or a passing subset of contrasts cannot replace per-sample support.

The older missingness/count-support fractions, per-condition replicate/coverage thresholds, positive variance floor, repeatability and strong condition-response gates are no longer active input filters for this policy. They may remain diagnostics. Choose or import the correct count mode and retain JC/JCEC or other measurement provenance. External PSI-only input remains explicitly `read_support_not_assessed`; it must not be represented as having passed the ten-count rule.

For imported rMATS network matrices, inspect the source-selection table: the adapter retains its canonical source choice per condition/event type (largest file, then lexical path/side order). The new QC checks every sample in that resolved matrix. It does not compare all alternative tables for agreement, so it cannot certify that repeated sample/event observations are identical across comparisons. Counts are not assembled or median/max-aggregated across comparison tables. By default, `max_events=0` uses all eligible events; a positive cap deliberately limits a run to a recorded, variance-ranked exploratory subset.

Inspect module membership, activities, Hub evidence and the actual fitted sample count. The integrated pipeline keeps each comparison's discovery, associations, network figures and enrichment separate. A previously completed cohort-wide network can be imported with its original scope; importing it does not turn it into independent per-comparison fits. Module association, preservation and functional annotation have separate statistical meanings. Reopening old network results preserves their original policy and does not recompute them. See [network scope](NETWORK_SCOPE.md) and [execution and ranking](NETWORK_EXECUTION.md).

Integrated runs rank Hub events by `|kME|`, retaining signed kME and complete
rankings. Bootstrap is off by default. Setting `network.bootstrap.enabled` in
the flow configuration enables the optional resampling stage; ordinary kME is
not described as bootstrap confidence. Standalone network configurations retain
their explicit scope and sample minimum.

New-policy runs record `network-all-samples-counts-v1` in their resolved configuration and output provenance. Legacy QC settings may be retained as diagnostics but cannot reactivate retired gates. Do not point the new standalone run at a historical output directory: an output missing the new policy is protected against overwrite, including with force. The application uses the independently versioned engine identity for new-policy output and cache locations.

## Network module settings

The standard defaults are 20 splice events per discovery module and WGCNA merge
cut height 0.25. Preserve custom project YAML settings by leaving the optional
`min_module_size` and `merge_cut_height` fields blank in the interface, or omitting
them / supplying JSON `null` through the API. An explicit runtime value overrides
that run only; it does not edit the input YAML. The merge height must be strictly
between zero and one; module size must be an integer of at least two.

For an explicit override, save a parameters file such as `network-parameters.json`:

```json
{"scope": "discovery", "min_module_size": 20, "merge_cut_height": 0.25}
```

Inspect the plan before running:

```sh
python -m katana_next.cli --project ../analysis-project/katana.project.json plan programs.cosplice --parameters network-parameters.json
python -m katana_next.cli --project ../analysis-project/katana.project.json run programs.cosplice --parameters network-parameters.json
```

The same values work through `project.plan("programs.cosplice", parameters)` and
`project.run("programs.cosplice", parameters)` in the Python interface. These
commands run an analysis; rendering a report does not run it. Changed effective
module settings select a separate engine identity; omitted or equal-valued
overrides preserve the compatible engine cache. An engine stage must still pass
its own cache validation before reuse.

Resolved configuration, output metadata and report methods record the actual
values used. The evidence report lists values per source run. The interpreted
one-click report follows its actual module-card input chain to the network run;
if that chain or its recorded values are unavailable, it says **Not recorded**.
It never substitutes today's defaults or an unrelated newer run. New interpreted
reports use a separate run directory, preserving previous deliveries.

## Optional network labels

`programs.cosplice`, `programs.psi_network` and `programs.hub_events` default to
`apply_threshold_labels=false`. Keep this default to export continuous evidence
and complete rankings without core/peripheral/stability or hub-candidate labels.
To add the documented descriptive presets, explicitly supply:

```json
{"apply_threshold_labels": true}
```

This adds labels only; it does not delete members or refit module assignments.
`not_applied` means labeling was disabled. It is different from an estimable
candidate failing a selected preset, an analysis not run, or a metric that could
not be estimated. An optional preset-selection frequency depends on that preset;
it is not the probability of a true biological hub. The legacy composite
membership score is also identified as threshold-dependent evidence. Raw
measurements and complete per-metric ranks remain available.

Report methods record the source run's actual label policy and switch, while
preserving the default FDR reporting level of 0.05 and the native correction
families. Historical results are not relabeled by opening or exporting them.

The standalone stage/CLI route uses the same opt-in setting in project YAML:

```yaml
consensus:
  clustering_method: backbone_confidence_only
  apply_threshold_labels: false
```

Its current tables and report companion are written separately under
`label_views/`; raw scientific tables and the frozen report remain internal or
historical references. Use `CoSplicePipeline.run_stage(...)` or the corresponding
CLI stage command for these versioned views. Direct legacy numerical stage
methods retain their internal raw schema. Historical refinement modes cannot
be presented as an unfiltered partition merely by disabling their labels.

## Seven-section study reports

Generate the fixed study report from completed results:

```sh
katana-report --project results.json --outdir exports/study-report --compact --no-archive
```

Open `exports/study-report/REPORT.html`. Use the report's comparison selector
to show one selected comparison at a time. The seven sections are design and
QC; AS landscape; DEG/DAS/DTU/DEXSeq; cross-comparison results; GO/KEGG and
pathways; co-splicing networks; and RBP/MotifTrace. Per-comparison Top 10 tables,
network graphs and full member/gene tables are retained. MotifTrace displays
significant motifs with their native P/FDR and Top 10 details; complete model
and diagnostic tables remain downloadable. The
[report specification](REPORT_SPECIFICATION.md) defines the fixed layout.

RS-DAS is the default (`--event-set stable_significant`). To choose ordinary DAS
initially, use `--event-set significant`. The report-wide selector synchronizes
supported downstream views and Top 10 tables; modules requiring RS retain that
explicit requirement. Switching a displayed set never refits source statistics.

For an existing portable report, use a new output directory:

```sh
katana-report --from-report portable-report --outdir exports/refreshed-report --compact --no-archive
```

Use `--archive downloads/study-report.zip` instead of `--no-archive` for a
download package. For inputs moved between machines, use explicit `--path-map`
bindings. See [completed sources and file identities](EXISTING_RESULTS.md) for
result descriptors, multiple completed collections and matched nf-core output
mapping. An arbitrary nf-core output folder is not automatically a complete
Katana project; nf-core/rnaseq alone does not supply all splicing endpoints.

`katana-report` assembles completed evidence. `katana-flow report` is the
integrated workflow's strict completion gate after all required calculations
have passed. Missing, failed, not-estimable and genuinely zero-finding analyses
remain distinct. Neither renderer invents missing analyses or performs new
statistical fits. Exclusions are applied to the relevant analysis inputs and
recorded; an unfiltered historical result is never silently relabeled.

The following compatibility export workflows also remain available. They do
not replace the fixed seven-section study-report entry above:

| Workflow | Purpose | Prerequisites |
| --- | --- | --- |
| `reports.bundle` | Reproducibility evidence package: HTML/PDF summary, inventory, provenance and ZIP. | Compatible completed artifacts in the project's catalog. |
| `reports.interpreted_bundle` | Import/export of a compatible historical study delivery and its attachments. | A delivery configuration referencing completed standardized delivery folders, or an approved reader report with its tables and figures. |

Run the general evidence package after confirming its plan:

```sh
python -m katana_next.cli --project ../analysis-project/katana.project.json run reports.bundle
```

For the interpreted route, provide a parameters JSON containing `config_path`. Its delivery configuration declares `study`, `result_root` and `manifest_path`; the optional approved-reader route also declares `approved_reader_root`. This route is not an unrestricted report generator for any arbitrary folder.

## Gene/Event deep dive

Start `katana-open`, select rows in completed result tables and add genes or
exact events to the target basket. Export the request ZIP; BAMs and large
references remain outside it. On the machine holding the read resources:

```sh
katana-targets inspect --request targets.zip
katana-targets plan --request targets.zip --bindings bindings.json
katana-targets run --request targets.zip --bindings bindings.json --outdir target-results
katana-targets results --outdir target-results
```

Use the exported resource IDs in `bindings.json` to bind that machine's BAMs,
indexes and references. The independent reports include Gene/Event cards,
available expression/splicing/usage evidence, per-sample data, requested
Sashimi, figures and complete interactive tables. A gene request retains all
available measured events for that gene, including non-significant events.
It does not change the main seven-section report or repeat upstream/network/
MotifTrace fits.

Open `INDEX.html` or import returned results into the Katana browser. To
restore a returned ZIP into a new resumable workspace:

```sh
katana-targets restore --results completed-targets.zip --outdir restored-targets
```

The same CLI runs locally or within a server allocation. It does not submit
Slurm/PBS jobs automatically. See [target commands, scheduling and reuse](TARGET_COMMANDS.md).

## Python interface

```python
from katana_next.client import KatanaProject

with KatanaProject("../analysis-project/katana.project.json") as project:
    inventory = project.inventory()
    plan = project.plan("reports.bundle")
    if plan["ready"]:
        completed_job = project.run("reports.bundle")
```

The same package provides completed-report assembly directly from Python:

```python
from katana_next.reports import build_report

manifest = build_report("results.json", "exports/python-study-report", compact=True, create_archive=False)
```

[MotifTrace](MOTIFTRACE.md) and [target analysis](TARGET_COMMANDS.md#python-interface-and-direct-export)
also expose public Python interfaces. Use the locked standalone Linux runtime
for new native calculations; Python-only installation does not install R or
external tools. Qt remains optional and is not required for CLI/browser use.
