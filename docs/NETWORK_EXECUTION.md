# Running independent comparison networks

The flow creates one observed network for each selected contrast by default.
Each network uses only the biological samples in that contrast. Its condition
associations, enrichment and report views must reference that same fit.
A cohort-wide network is a separate, explicitly selected analysis.

## Sample count

The integrated defaults are six independent samples in total per comparison,
or 15 for a cohort-wide network. Set `network.parameters.minimum_network_samples`
in `flow.json` to choose a different minimum, such as `4` for two-plus-two.
The engine accepts integers >=3. The upstream differential route still needs
at least two samples per group. Actual sample count and the requested minimum
are recorded; a computational minimum does not establish statistical power.

## Membership and optional resampling

New flow plans do not run bootstrap by default. They create an observed fit,
matched condition associations and module annotation. Events are ranked within
each module by descending absolute eigengene membership `|kME|`, retaining the
signed Pearson correlation and its P/FDR columns. Gene scores are the maximum
absolute kME among that gene's events within the module. Exact ties use ascending
event IDs or gene IDs. Connectivity and historical confidence columns remain
available as separate metrics; they do not contribute to the `abs_kME` score.
Neither kME nor the number of high-ranking events is a resampling stability test.

The sealed flow configuration explicitly records
`network.bootstrap.enabled: false` and `network.parameters.hub_metric: abs_kME`.
The observed fit disables both edge and graph bootstrap iterations. Disabled
bootstrap creates no bootstrap task, no synthetic success receipt and no claim
that stability was calculated. Cached fit products suffice to export kME tables;
this does not rerun discovery or change module membership.

To request resampling explicitly, set `network.bootstrap.enabled: true` and
`network.parameters.bootstrap_iterations` to an integer of at least two. That
plan adds a real bootstrap stage, which must finish before the optional stability
analysis can be accepted. Old sealed plans and standalone CoSplice configurations
retain their recorded behavior and composite ranking unless explicitly migrated.

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
    "bootstrap": {"enabled": false},
    "parameters": {
      "tom_dense_max_events": 65000,
      "consensus_dense_max_events": 65000,
      "wgcna_timeout_seconds": 86400,
      "hub_metric": "abs_kME"
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

Dynamic Tree Cut receives the existing square distance matrix directly, rather
than converting it to a condensed vector and reconstructing submatrices with
Python loops. The sensitivity grid prepares that distance once for all recorded
depths. This changes distance access and avoids repeated allocations; it does
not remove a sensitivity depth or change the tree, cut, PAM options or event
set. Regression checks compare complete label arrays against the previous
condensed-distance path, including tied distances.

The Python network diagnostics and the official R WGCNA module fit currently
use distinct correlation/TOM implementations. They are retained separately:
substituting one for the other would change scientific outputs. The execution
optimization does not claim to collapse these distinct products into one TOM.

Nextflow exports its actual `task.cpus` as the worker thread budget. Python BLAS
and the R child inherit the allocation; R WGCNA threading is enabled explicitly
within that budget. Standalone execution also respects the available CPU
affinity and an explicit thread budget. Thread counts are execution settings,
not changes to the correlation method or module parameters. Serial tree steps
are not guaranteed to occupy every allocated core. Timestamped stage messages
allow operators to distinguish matrix work, module fitting and tree cuts.

## Checking a new deployment

A representative comparison is useful when qualifying a new machine, image or
resource configuration. This is a deployment check, not a manual gate between
comparisons in ordinary runs. Nextflow schedules each ready comparison automatically.
For deployment qualification:

1. Confirm the exact sample IDs, biological units and complete event universe.
2. Require successful native completion and verified observed modules, activity,
   numerical tree, edge outputs and resolved configuration.
3. Verify the signed kME, absolute membership ranks and the exact recorded hub
   policy. Only when bootstrap was explicitly enabled, complete every declared
   biological-unit resample and validate its input, sample draw, kernel and output
   identities before reduction.
4. Validate matched condition effects, enrichment foreground/background and hit
   genes, network figures, and the comparison-specific report links.
5. Inspect peak memory and runtime, then schedule independent remaining
   comparisons through the cluster scheduler. Preserve completed products.

Passing container installation or small-fixture tests is not full-scale network
acceptance. Failed attempts and incomplete bootstrap results remain visible.
