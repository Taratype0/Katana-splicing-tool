# Katana CoSplice v1.2.0

## Method upgrade

- Added beta-posterior draw ensembles for event-correlation uncertainty diagnostics.
- Retained the official signed WGCNA soft-adjacency, TOM, Dynamic Tree Cut and module-merging backbone.
- Added cached posterior correlation mean, SD and sign-probability matrices; the production partition remains posterior-mean WGCNA because the matched Setd2 ablation found over-merging in the ensemble partition arm.
- Extended condition-stratified co-assignment bootstrap with one PSI posterior draw per run without changing the full-data WGCNA partition.
- Added bootstrap kME and hub-rank stability plus module program-score composition confidence intervals.
- Added bootstrap confidence intervals to fixed-module gained/lost/rewired tests.
- Added explicit observed-versus-residual program classes while retaining continuous overlap statistics.
- Tightened cross-annotation event mapping to exact canonical keys, one-to-one validated mappings or unambiguous same-build splice-site signatures. No arbitrary coordinate window is used.

## Reproducibility

- Posterior alpha/beta matrices and ensemble correlation caches are tracked artifacts.
- Cache fingerprints include uncertainty draw count, seed, model settings and input signatures.
- The real-count ablation compares raw WGCNA, posterior-mean WGCNA, posterior-ensemble WGCNA and full CoSplice confidence annotation.

## Frozen Setd2 ablation result

Using 3,000 Stable-AS events, 36 samples, seed 37 and 10 condition-stratified resamples:

- raw-PSI WGCNA: mean bootstrap ARI 0.546; five modules; largest-module fraction 0.379;
- posterior-mean WGCNA: mean bootstrap ARI 0.559; six modules; largest-module fraction 0.327;
- posterior-ensemble partition: mean bootstrap ARI 0.914, but largest-module fraction 0.799 and 99.6% assignment, so it was rejected as over-merged;
- full CoSplice posterior-draw confidence stress test: mean ARI 0.367; 138 core, 785 peripheral, 701 low-confidence assigned and 1,363 grey events.

The production claim is therefore improved uncertainty visibility and a small posterior-mean stability gain, not universal clustering superiority over WGCNA.

## Interpretation limits

- Posterior sign probability measures edge-direction stability, not differential-splicing significance.
- Observed/residual preservation is descriptive decomposition, not causal attribution.
- With three replicates per condition, fixed-module differential-network results remain prioritization evidence requiring external or orthogonal validation.
