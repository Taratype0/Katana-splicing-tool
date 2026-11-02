# Whole-study CoSplice sample design

Decision confirmed by the user on 2026-11-02.

- No minimum number of biological units is imposed as a study eligibility gate.
- At least 12 independent biological units are recommended. This is the project recommendation, not a literature-derived significance threshold or guarantee of power.
- Paired analysis is optional. It is enabled by explicit comparison design or the standalone `--paired` option, never by silently treating repeated measurements as independent samples.
- Preserve the same `biological_sample_id` for observations of one individual. Six observations from three pairs remain three independent biological units.
- Each paired comparison must have exactly one observation per biological unit in each arm. Technical replicates must be merged before this analysis.
- A shared network uses all declared comparable study samples. Each requested contrast retains its own direction and paired/unpaired association model. Combining conditions does not increase the independent-unit count when the individuals are repeated.
- Statistical estimability is separate from study eligibility. Missing residual degrees of freedom, singular designs or unestimable variances produce explicit unavailable tests; they do not become significant results.

This change applies only to the independent integration branch. The original sealed software and completed results remain unchanged.
