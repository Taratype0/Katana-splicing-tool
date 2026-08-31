# CoSplice Multi-objective Validation Audit

> Production status: since v1.1.0, `backbone_confidence_only` preserves the
> official WGCNA partition. The objectives below audit a fixed partition and
> cannot alter module membership. Candidate selection remains available only
> for explicit reproduction of historical research runs.

## Purpose

The publication-oriented audit asks whether the fixed WGCNA partition is:

1. reproducible under discovery-cohort sample bootstrap;
2. structurally preserved in predeclared tuning cohorts;
3. functionally coherent against a versioned tuning gene-set collection.

These metrics support validation and reporting; they are not allowed to tune the
production module definition.

## Data separation

- **Discovery cohort**: builds the co-assignment matrix and candidate partitions.
- **Tuning cohorts**: contribute preservation evidence to the fixed-partition audit.
- **Holdout cohorts**: are never loaded during this audit and remain confirmatory.
- **Tuning GMT**: contributes functional evidence; enrichment against the same GMT is not independent
  validation and cannot change production membership.

Tuning and holdout cohort names must be disjoint. The stage fails if configured tuning inputs are
missing. Functional evidence is adaptive by default: if fewer than the predeclared number of
modules pass the functional FDR cutoff, the functional objective is disabled and its weight is
renormalized across internal stability and tuning preservation. This prevents a weak raw enrichment
score from distorting the audit; explicit historical selection modes retain their original behavior.

## Historical candidate partitions

Historical candidate generation had two frozen levels. First, a discovery-only pilot reuses one
TOM hierarchy and compares predeclared cut heights using lower-tail bootstrap ARI plus event
coverage. The selected TOM height is then frozen before the outer bootstrap and before holdout
access. Second, hierarchical cuts of the outer-bootstrap co-assignment matrix supply the final
candidate partitions. Grey-fraction, largest-module-fraction and module-count constraints are
applied before objective optimization. Sparse Leiden candidates remain available only for legacy
cache compatibility.

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
preservation statistic. If no candidate has an eligible tuning preservation test, the preservation
objective is disabled and the remaining active weights are renormalized; an unevaluable objective
cannot force every joint score to zero.

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

## Joint audit score

The configured weights are:

- internal stability: 0.40;
- tuning preservation: 0.35;
- functional coherence: 0.25.

For audit and historical reproduction, the joint score is:

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
Under `auto_significance`, a failed functional gate sets the effective functional
weight to zero. An unevaluable tuning-preservation objective is handled in the same conservative
way. All remaining active weights are renormalized to sum to one and exported. `always` requires
evaluable formal functional signal; `disabled` always excludes it.

## Outputs

`04c_consensus_modules/` contains:

- `consensus_parameter_search.tsv`: TOM pilot cuts and final consensus candidates, active objectives,
  Pareto status and final score;
- `multiobjective_selection_audit.tsv`: selected-partition module/cohort and functional contributions;
- `consensus_diagnostics.tsv`: selected objective values and method;
- `parameters.json`: configured and effective weights, functional-gate state, tuning/holdout roles,
  input fingerprint and interpretation limits;
- `plots/multiobjective_candidate_tradeoff.png`: preservation versus functional coherence, with
  point size representing internal stability.

## Publication boundary

The production audit does not prove generalization or biological mechanism and
does not manipulate module membership. A methods claim still requires:

1. a locked partition;
2. untouched external holdout cohorts;
3. formal permutation or Bayesian preservation;
4. functional replication using an independent source where possible;
5. ablation against internal-stability-only CoSplice and event-level WGCNA.
