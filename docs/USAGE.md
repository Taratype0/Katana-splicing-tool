# Using Katana

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

For count-bearing network discovery, declare the complete fitted cohort and use the `network-all-samples-counts-v1` input policy. Every participating sample must have valid PSI and IJC+SJC >=10 for an event to enter; events with a constant PSI vector are excluded. Do not use significant DAS events as the default input universe. Count and sample mapping must be available; a fraction of supported samples, a group average or a passing subset of contrasts cannot replace per-sample support.

The older missingness/count-support fractions, per-condition replicate/coverage thresholds, positive variance floor, repeatability and strong condition-response gates are no longer active input filters for this policy. They may remain diagnostics. Choose or import the correct count mode and retain JC/JCEC or other measurement provenance. External PSI-only input remains explicitly `read_support_not_assessed`; it must not be represented as having passed the ten-count rule.

For imported rMATS network matrices, inspect the source-selection table: the adapter retains its canonical source choice per condition/event type (largest file, then lexical path/side order). The new QC checks every sample in that resolved matrix. It does not compare all alternative tables for agreement, so it cannot certify that repeated sample/event observations are identical across comparisons. Counts are not assembled or median/max-aggregated across comparison tables. By default, `max_events=0` uses all eligible events; a positive cap deliberately limits a run to a recorded, variance-ranked exploratory subset.

Inspect module membership, scores, hub evidence and stability before comparing conditions. A network fitted over a cohort is not refitted independently inside each small contrast by default. Module association, preservation and functional annotation have separate statistical meanings. This input approval does not change raw BAM-junction measurement settings or the existing module-size/merging algorithm. Reopening old network results preserves their original policy and does not recompute them.

New-policy runs record `network-all-samples-counts-v1` in their resolved configuration and output provenance. Legacy QC settings may be retained as diagnostics but cannot reactivate retired gates. Do not point the new standalone run at a historical output directory: an output missing the new policy is protected against overwrite, including with force. The application uses the independently versioned engine identity for new-policy output and cache locations.

## English reports

The two existing report entries remain available:

| Workflow | Purpose | Prerequisites |
| --- | --- | --- |
| `reports.bundle` | Reproducibility evidence package: HTML/PDF summary, inventory, provenance and ZIP. | Compatible completed artifacts in the project's catalog. |
| `reports.interpreted_bundle` | Interpreted study report: findings, design, evidence, limitations and complete attachments. | A delivery configuration referencing completed standardized delivery folders, or an approved English reader report with its tables and figures. |

The interpreted route preserves the existing report structure: main findings; experimental design; abundance and splicing evidence; replicate patterns; qualified interactions; functional interpretation; network modules; usage/candidates; supporting regulatory evidence; limitations and reproducibility attachments. Content depends on the actual completed delivery. It does not fabricate a section's analysis when prerequisites are missing.

Run the general evidence package after confirming its plan:

```sh
python -m katana_next.cli --project ../analysis-project/katana.project.json run reports.bundle
```

For the interpreted route, provide a parameters JSON containing `config_path`. Its delivery configuration declares `study`, `result_root` and `manifest_path`; the optional approved-reader route also declares `approved_reader_root`. This route is not an unrestricted report generator for any arbitrary folder.

Both report routes read completed evidence. Missing, failed, unassessable and genuinely zero-finding analyses remain distinct. Gene exclusion policies and filters must be applied to the relevant results and documented; an unfiltered delivery cannot simply be relabeled as filtered. Report output is English.

## Python interface

```python
from katana_next.client import KatanaProject

with KatanaProject("../analysis-project/katana.project.json") as project:
    inventory = project.inventory()
    plan = project.plan("reports.bundle")
    if plan["ready"]:
        completed_job = project.run("reports.bundle")
```

This interface currently shares the installed application dependencies. It does not yet constitute a separately packaged, minimal headless analysis core.
