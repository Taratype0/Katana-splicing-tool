# Running independent comparison networks

The flow creates one observed network for each selected contrast by default.
Each network uses only the biological samples in that contrast. Its bootstrap,
condition associations, enrichment and report views must reference that same fit.
A cohort-wide network is a separate, explicitly selected analysis.

## Full event sets and execution limits

The network includes all events passing the configured measurable-PSI checks.
The execution limit is a guard against an unbudgeted dense matrix allocation;
it is not a filter that retains the first or most variable events. When the
event count exceeds the guard, execution fails rather than reducing the input.

For an existing-results flow, explicit execution settings belong in
`network.parameters` in the flow configuration, for example:

```json
{
  "network": {
    "scope": "contrast",
    "parameters": {
      "tom_dense_max_events": 65000,
      "consensus_dense_max_events": 65000,
      "wgcna_timeout_seconds": 86400,
      "bootstrap_iterations": 50
    }
  }
}
```

These example limits require a separately sized scheduler memory/time request.
They do not allocate memory themselves. Both event guards and the R timeout are
carried into the observed fit and its bootstrap continuation and recorded in the
resolved configuration. Defaults remain conservative; do not raise the guards
without measuring a representative fit's peak memory and runtime.

The current `tom_hierarchical` / `wgcna_blockwise` implementation uses a full
event-by-event TOM. The WGCNA function name does **not** mean that memory is
bounded to small blocks: its block size equals the event universe. Memory scales
quadratically and matrix multiplication can dominate runtime. Reducing the
WGCNA block size would change the discovery procedure and is not a transparent
execution fix.

## Export and reproducibility

Large-network edge export selects the same largest adjacency values with
bounded temporary index storage. The existing edge export cap affects saved
edges only, never module discovery, event inclusion or connectivity. For large
networks, ties at the export cutoff use ascending event-matrix row/column order;
the rule is recorded in `network_diagnostics.tsv`. Small-network exports retain
their previous ordering.

Unweighted WGCNA execution no longer calculates an unused correlation matrix
before `blockwiseModules` calculates its own. TOM fitting without a gene prior
also avoids unused prior matrices and index arrays. These changes preserve the
declared correlation, signedness, power, module size and merge height.

## Acceptance before scaling out

Validate a single real comparison before submitting the remaining comparisons:

1. Confirm the exact sample IDs, biological units and complete event universe.
2. Require successful native completion and verified observed modules, activity,
   numerical tree, edge outputs and resolved configuration.
3. Complete every declared biological-unit bootstrap and validate its input,
   sample draw, kernel and output identities before reduction.
4. Validate matched condition effects, enrichment foreground/background and hit
   genes, network figures, and the comparison-specific report links.
5. Inspect peak memory and runtime, then schedule independent remaining
   comparisons through the cluster scheduler. Preserve completed products.

Passing container installation or small-fixture tests is not full-scale network
acceptance. Failed attempts and incomplete bootstrap results remain visible.
