# Katana CoSplice v1.2.2

## Continuous membership confidence

- Kept the posterior-mean signed WGCNA/TOM partition fixed.
- Added an event-level continuous confidence score from co-assignment, bootstrap
  kME, split-half confirmation and optional IJC/SJC read reliability.
- Added per-module count-matched continuous-confidence and top-kME arms so a
  smaller event set cannot receive automatic preservation credit.
- Exported every score component, within-module rank and percentile for audit.

## Formal external preservation

- Added official `WGCNA::modulePreservation` output.
- Added NetRep seven-statistic preservation with per-statistic FDR control.
- Retained the existing WGCNA-like Zsummary only as a secondary diagnostic.
- Added parallel external-cohort and within-cohort arm execution with isolated
  temporary directories; no NetRep cache is written into the source tree.
- Added a formal preservation plot and machine-readable R package/version audit.
- Updated `scripts/install_wgcna.R` to install and verify both WGCNA and NetRep.

## Real-count Setd2 ablation

Using 3,000 Stable-AS events, 36 samples, seed 20260901 and 20 matched
condition-stratified bootstraps:

- posterior-mean WGCNA retained six modules and 1,624 assigned events;
- posterior-mean sample-bootstrap ARI was 0.556 versus 0.557 for raw PSI;
- the stricter posterior-draw membership stress-test ARI was 0.356 and is not
  interpreted as a clustering improvement;
- 926 events were split-half core, 464 peripheral and 234 uncertain;
- among 926 count-matched events, continuous-confidence selection increased
  median read reliability from 0.698 to 0.718 and co-assignment component from
  0.427 to 0.495 relative to top-kME, while median kME component decreased from
  0.890 to 0.882.

This internal trade-off does not establish external superiority. The formal
public-cohort result used BRCA discovery and KIRC, KIRP, KICH, LUAD and LUSC
holdouts with 1,000 permutations per arm/cohort. Across 15 matched
module-cohort pairs, continuous confidence versus top-kME had an official
WGCNA median Zsummary difference of -0.052 (FDR 0.083; 5 wins and 10 losses),
while the NetRep significant-statistic count was identical (FDR 1.0). The new
score is therefore retained for audit and multi-evidence prioritization, not as
a replacement for kME or as an external-superiority claim.
