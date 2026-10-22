# Using Katana

## Projects and inputs

Open an existing v2 Katana project manifest in the interface. It records samples, contrasts, input assets, reference annotations and output locations. Supported legacy projects can be migrated into a separate output location:

```sh
python -m katana_next.migration legacy-project.yaml --output-root ../analysis-results --manifest ../analysis-project/katana.project.json
```

Replace the example file names with your own paths. Migration discovers compatible files; it does not make missing references or analyses available. Check the resulting manifest and readiness report before running analyses.

Use explicit sample IDs and condition labels. Preserve donor/biological replicate identities and distinguish them from technical measurements. Inspect numerator and denominator before interpreting an effect. Select the contrasts relevant to your question rather than treating every possible comparison as a required analysis.

## Workflow execution

The [threshold reference](THRESHOLDS.md) records the current bulk values, exact units, literature support and limitations. In particular, the two 75% replicate-consistency checks are Katana heuristics, and approximate PSI standard error is diagnostic by default.

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

## MotifTrace and networks

MotifTrace exposes event-level sequence consequences, motif and RBP evidence, and links back to their source inputs. Choose the appropriate evidence layer: RBP identity annotations, boundary motif/reference binding evidence, and mature-RNA sequence consequences answer different questions. Native method results and valid imported results retain their source identities. The presence of a resource does not establish that every optional author-native fit has run.

Network discovery uses the qualified cohort and its declared event universe. Inspect module membership, scores, hub evidence and stability before comparing conditions. A network fitted over a cohort is not refitted independently inside each small contrast by default. Module association, preservation and functional annotation have separate statistical meanings.

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
