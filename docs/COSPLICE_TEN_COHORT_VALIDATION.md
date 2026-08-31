# CoSplice locked 10-cohort validation

## Objective

This benchmark asks three separate questions:

1. Does CoSplice produce a reproducible, non-trivial event partition within a cohort?
2. Are discovery modules preserved when projected into an external cohort with a compatible event catalogue?
3. Do module-level genes and functional themes replicate independently of module discovery?

The benchmark does not assume that CoSplice is superior. Claims are emitted only when
the locked metrics support them.

## Public panel

The locked panel contains ten independent TCGA disease cohorts from TCGA SpliceSeq:
KIRC, KIRP, KICH, LUAD, LUSC, COAD, READ, BRCA, PRAD and LIHC. All use exon-skipping
events downloaded with the same filters and exact `gene_symbol|event_type|AS_ID`
matching.

These are independent patient cohorts, but not ten independent studies: they share the
TCGA sequencing programme and SpliceSeq processing ecosystem. Kidney, lung and
colorectal projections are designated organ-matched primary tests. Other projections
are cross-organ stress tests.

## Locked design

- Discovery event universe: top 500 variable eligible events selected within the
  discovery cohort only.
- External event universe: every eligible target event; no independent top-500
  truncation.
- Primary comparator: official R WGCNA on the identical discovery event universe.
- Internal stability: 20 matched sample bootstraps.
- External preservation: 500 permutations per module.
- Event matching: exact TCGA SpliceSeq event key; no fuzzy coordinate matching.
- Biological annotation is evaluated after module discovery and never tunes modules.

## Claim rules

- Internal stability is compared only where both methods yield valid multi-module
  partitions.
- A claim about the strongest preserved program requires CoSplice to have a higher
  best-module Zsummary in at least 70% of primary projections, a positive median
  paired difference and a one-sided paired Wilcoxon P value below 0.05. The report
  also retains the exact sign-test P value and a bootstrap confidence interval.
- A broader claim about the complete module set additionally requires a positive
  median-Zsummary difference and positive program-score preservation difference.
  This prevents a method with more modules from winning only because it has more
  chances to produce one high maximum Zsummary.
- Stronger functional replication requires a higher cross-cohort functional-theme
  replication metric.
- A negative or inconclusive result is retained as a result, not hidden.

## Scientific basis

- TCGA SpliceSeq provides PSI values across TCGA disease cohorts:
  <https://pubmed.ncbi.nlm.nih.gov/26602693/>
- WGCNA defines the primary network comparator:
  <https://doi.org/10.1186/1471-2105-9-559>
- Module preservation statistics motivate density/connectivity-based external
  validation:
  <https://doi.org/10.1371/journal.pcbi.1001057>
- Published cosplicing work demonstrates cross-species module preservation as a
  relevant validation concept:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC4429622/>

## Locked result

The formal run used 4,490 tumour samples, 500 discovery events per cohort,
20 matched bootstraps, and 500 preservation permutations per projected module.
All 90 directed cohort projections completed.

### Internal stability and partition usability

- CoSplice produced a valid non-trivial multi-module partition in 10/10 cohorts.
- Exact WGCNA produced a valid non-trivial multi-module partition in 4/10 cohorts.
- In the four cohorts where both partitions were valid, exact WGCNA had the
  higher bootstrap stability in all four. Therefore the statement that CoSplice
  is more internally stable than exact WGCNA is not supported.
- The supported claim is narrower: CoSplice produced a scientifically usable
  multi-module partition more consistently on this panel.

### Organ-matched external preservation

The primary analysis contains ten directed kidney, lung, and colorectal
projections.

- Every projection contained at least one CoSplice module with a higher
  best-module Zsummary than exact WGCNA.
- Median paired best-module Zsummary difference was +7.066
  (95% bootstrap CI +5.178 to +7.726).
- One-sided exact sign-test and paired Wilcoxon P values were both 0.0009766.
- This maximum statistic has a module-number opportunity effect. CoSplice had
  more valid modules and therefore more chances to contain one strong program.
- Across the complete module set, median module Zsummary difference was -0.635
  and median program-score preservation difference was -0.035. A claim that the
  complete CoSplice module set is better preserved is therefore not supported.

The defensible preservation statement is:

> In each organ-matched projection, CoSplice recovered at least one strongly
> preserved splicing program, but preservation was not uniformly stronger
> across its complete module set.

### Functional replication

- CoSplice generated 35 modules, of which 28.6% had significant GO/KEGG
  annotation; exact WGCNA generated 15 modules, of which 20.0% were annotated.
- Organ-matched term-level replication was evaluable for only two CoSplice
  projections and had median best Jaccard 0. Exact WGCNA did not provide enough
  eligible organ-matched functional projections for a paired test.
- Stronger cross-cohort functional-theme replication is not supported by this
  run. More independent studies with harmonized event and gene annotation are
  required.

## Cross-study robustness boundary

A separate challenge used three heterogeneous public datasets: GSE62669 mouse
brain exon usage, GSE141633 prostate bulk PSI, and GSE85908 single-cell PSI.
Their event catalogues are incompatible, so they test internal robustness and
partition usability rather than exact module preservation.

- CoSplice produced valid multi-module partitions in 3/3 datasets.
- Exact WGCNA produced a valid multi-module partition in 1/3 datasets.
- In GSE62669, where both methods were valid, exact WGCNA was more stable
  (minimum ARI 0.347 versus 0.194; median ARI 0.522 versus 0.290).
- In GSE141633 and GSE85908, exact WGCNA assigned every tested event to grey,
  whereas CoSplice retained two non-grey modules.

This challenge supports broader partition usability, not universal stability
superiority. It is not a substitute for additional independent cohorts with
compatible event catalogues and phenotype-aware functional validation.
