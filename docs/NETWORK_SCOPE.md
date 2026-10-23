# Networks for the selected experimental comparisons

New Katana flow plans use `network.scope: "contrast"`. Each selected comparison
has its own network discovery, bootstrap, condition association and module
annotation tasks. Only the samples in that comparison and its matching native
rMATS result files contribute to its PSI matrix. A third experimental group is
not pooled into that discovery.

```json
{
  "schema": "katana.flow.config/v1",
  "contrast_ids": ["KO_vs_CONTROL", "WT_rescue_vs_KO"],
  "network": {
    "scope": "contrast",
    "parameters": {"rmats_mode": "JCEC", "bootstrap_iterations": 50}
  }
}
```

The two discoveries are independent tasks in the same downstream pipeline.
Their bootstrap, association and annotation tasks depend on their own completed
fit. Each task records the comparison, ordered sample and biological sample IDs,
sample count and scope. These fields participate in plan and task identities;
an older cohort fit cannot be reused as an independent comparison fit.

At least six independent biological samples are required for this flow mode.
This is a supported formal execution mode. New task identities record the
neutral eligibility value `minimum_met` together with the actual sample count;
they do not classify the analysis as exploratory based on sample count.
Actual bootstrap diagnostics, continuous membership metrics and native P/FDR
remain in the outputs. Six samples meeting the execution minimum does not itself
establish module stability. Technical copies cannot increase biological n.
The standalone CoSplice default remains 15 samples. No frozen plan is rewritten.
Older recorded scope identities remain readable for exact result verification;
new plans and report labels use the current neutral contract.

Every module ID is local to its discovery. For example, `M001` from two contrasts
does not establish that these are the same module. Original condition-association
P values and the existing binary-trait multiple-testing family are retained
within each comparison. Module enrichment continues to use all measured genes
in that network, including unassigned genes, and the full tested term family.
The report's significance filter never recalculates FDR on displayed terms.

For a deliberate joint cohort analysis, set `network.scope: "cohort"`. This
retains the all-sample discovery, requires at least 15 independent samples, and
computes only the requested contrasts against its fixed module scores. It is a
different analysis scope and is never substituted for the contrast default.

An external YAML may be supplied as `network.config` for one discovery, or
through `network.configs` mapping individual comparison IDs to YAML files. Its
sample metadata must match the task exactly and include `biological_sample_id`
matching the project for every sample. A missing biological-unit column is
rejected before fitting. In contrast mode, native rMATS
file hashes, event types and sample-vector order must also match the declared
comparison. A broad cohort configuration is rejected instead of silently being
subsetted or ignored.

Existing finished cohort networks remain valid cohort results. Producing an
independent network for each comparison requires new discoveries with these
new scopes; changing a report heading cannot turn the old result into them.
The implementation's local tests validate planning, matrix transport and scope
isolation without refitting study networks or rerunning upstream analyses.
