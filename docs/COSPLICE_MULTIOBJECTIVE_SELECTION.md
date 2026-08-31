# CoSplice Multi-objective Module Selection

## Purpose

The publication-oriented selector chooses one consensus partition that is simultaneously:

1. reproducible under discovery-cohort sample bootstrap;
2. structurally preserved in predeclared tuning cohorts;
3. functionally coherent against a versioned tuning gene-set collection.

This replaces the weaker practice of selecting modules only for internal stability and describing
external preservation or enrichment after the fact.

## Data separation

- **Discovery cohort**: builds the co-assignment matrix and candidate partitions.
- **Tuning cohorts**: may influence parameter selection through preservation scores.
- **Holdout cohorts**: are never loaded during selection and remain confirmatory.
- **Tuning GMT**: may influence module selection; enrichment against the same GMT is not independent
  validation.

Tuning and holdout cohort names must be disjoint. The stage fails if configured tuning inputs are
missing. Functional evidence is adaptive by default: if fewer than the predeclared number of
modules pass the functional FDR cutoff, the functional objective is disabled and its weight is
renormalized across internal stability and tuning preservation. This prevents a weak raw enrichment
score from changing the selected partition.

## Candidate partitions

The existing Leiden grid search supplies candidate partitions across graph edge thresholds and
resolutions. Grey-fraction, largest-module-fraction and module-count constraints are applied before
objective optimization.

### Internal stability

For candidate partition \(c\):

\[
S_c = 0.50\,\overline{ARI}_c
    + 0.30\,median(ARI_c)
    + 0.20\,min(ARI_c)
\]

The lower-tail term prevents a partition with a poor bootstrap failure mode from being selected only
because its average ARI is high. Values are clipped to \([0,1]\) for the joint objective. Candidate
selection additionally applies a partition-validity factor combining event coverage, module-size
balance, membership purity and program separation. A one-module partition receives a 0.25 penalty,
because a broad single cluster can have high ARI without resolving distinct biological programs.

### Tuning-cohort preservation

For each module and tuning cohort, exact shared events are used to calculate:

- absolute-correlation density similarity;
- intramodular-connectivity correlation;
- module-membership-profile correlation;
- first-program variance similarity;
- exact shared-event coverage.

The module/cohort preservation contribution is a geometric composite. Candidate preservation
\(P_c\) is the median contribution across eligible module/cohort tests. This is a continuous,
permutation-free tuning surrogate. It is not the final WGCNA-like Zsummary or Bayesian holdout
preservation statistic.

### Functional coherence

Event annotations are collapsed to unique genes within each module. Terms are tested by one-sided
hypergeometric enrichment against the candidate's annotated gene universe. BH-FDR is controlled
within each module so candidates with more modules are not penalized solely by a larger global test
count.

For module \(m\), formal FDR evidence contributes 75% and raw enrichment strength contributes 25%
after annotation-coverage weighting. The raw component can rank candidates only after the formal
FDR activation gate has passed; it cannot activate the functional objective by itself.

For the formal component of module \(m\):

\[
F_m = annotation\_coverage_m \times
\min\left(1, \frac{-\log_{10}(bestFDR_m)}
{-\log_{10}(FDR_{cutoff})}\right)
\]

Candidate functional coherence \(F_c\) is the median across eligible modules.

## Joint objective

The configured weights are:

- internal stability: 0.40;
- tuning preservation: 0.35;
- functional coherence: 0.25.

After marking the Pareto frontier, the selected feasible candidate maximizes:

\[
J_c =
\exp\left(
w_S\log(\max(S_c,\epsilon)) +
w_P\log(\max(P_c,\epsilon)) +
w_F\log(\max(F_c,\epsilon))
\right)
\]

with default \(\epsilon=0.01\). The geometric objective makes compensation deliberately limited: a
near-zero preservation or active functional score cannot be hidden by excellent internal stability.
Under the default `auto_significance` mode, a failed functional gate renormalizes 0.40/0.35 to
0.533/0.467 and sets the effective functional weight to zero. `always` requires evaluable formal
functional signal; `disabled` always excludes it.

## Outputs

`04c_consensus_modules/` contains:

- `consensus_parameter_search.tsv`: all candidates, three objectives, Pareto status and final score;
- `multiobjective_selection_audit.tsv`: selected-partition module/cohort and functional contributions;
- `consensus_diagnostics.tsv`: selected objective values and method;
- `parameters.json`: configured and effective weights, functional-gate state, tuning/holdout roles,
  input fingerprint and interpretation limits;
- `plots/multiobjective_candidate_tradeoff.png`: preservation versus functional coherence, with
  point size representing internal stability.

## Publication boundary

The selector provides a scientifically motivated parameter-selection rule. It does not itself prove
generalization or biological mechanism. A methods claim still requires:

1. a locked partition;
2. untouched external holdout cohorts;
3. formal permutation or Bayesian preservation;
4. functional replication using an independent source where possible;
5. ablation against internal-stability-only CoSplice and event-level WGCNA.
