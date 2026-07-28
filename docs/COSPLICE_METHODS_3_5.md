# CoSplice methods: interpretation, replicate-aware statistics, and meta-analysis

This document defines the scientific scope of CoSplice stages 6E-6G. The stages are independent: they read cached upstream artifacts and never run another stage automatically.

## 6E Module functional interpretation

Inputs include cached module assignments, sample-level PSI, canonical rMATS event metadata, optional DEG/DEXSeq/DTU support, optional RBP expression (or Salmon `quant.sf` plus GTF), genome FASTA, RBP motifs, and GO/KEGG providers. Bootstrap-consensus modules are preferred. If they do not exist, the stage may explicitly use the cached observed network and records `module_source=observed_network_fallback`; it never starts consensus computation itself.

Outputs include event module membership (kME), hub events and genes, joint evidence, RBP-expression correlation, rMAPS2-inspired positional motif enrichment, GO/KEGG enrichment, dominant-isoform switches, switch enrichment, and candidate regulator ranking. Official ATtRACT ZIP bundles, ATtRACT text files, CisBP-RNA-style long PWMs, and custom consensus tables can be scanned automatically against an indexed genome FASTA. Overlapping splice-site windows are merged before strand-oriented scanning. Local GMT is preferred for an immutable database snapshot; optional g:Profiler supplies GO and KEGG with provider/version audit when no GMT is configured. Descriptive dominant-transcript switches are reported separately from the stricter `DTU_supported_switch` enrichment scope. Module eigengene/kME follows the WGCNA framework. Motif, expression, and switch evidence prioritize mechanisms but do not establish causality.

Motif scanning and isoform-switch quantification have independent input fingerprints. A change to PWM/FASTA or `quant.sf`/GTF invalidates only the affected support cache rather than forcing both expensive computations to rerun.

## 6F Replicate-aware module statistics

For every stable event assigned to the selected cached module result, Katana fits beta-binomial likelihoods to sample-level rMATS inclusion/skipping junction counts:

- null: intercept plus optional user-configured batch;
- full: intercept plus optional batch plus condition.

The likelihood-ratio output includes global Benjamini-Hochberg FDR, overdispersion `rho`, read support, convergence diagnostics, and condition mean PSI. This complements, but never replaces, the hierarchical rMATS result.

Module program scores use OLS with HC3 robust covariance by default. If `random_effect_column` identifies repeated donor/block observations, event counts use a beta-binomial random-intercept likelihood integrated by Gauss-Hermite quadrature, and program scores use a linear mixed model random intercept. Model convergence and failures are retained in diagnostics. Only condition, an explicitly configured optional batch, and an explicitly configured repeated block are modeled. Sex, age, and cell type are not inserted implicitly. Without a valid repeated donor/block identifier the output is not described as mixed effects.

The stage performs a declared head-to-head comparison against simple maximum condition-mean delta PSI and reports rank correlation and top-N Jaccard overlap.

## 6G Cross-cohort module meta-analysis

Event harmonization first uses exact canonical event keys. Within the same genome build it can additionally match exact event type, strand, and splice-site coordinate signatures when gene/annotation labels changed. Across genome builds, a user-validated explicit mapping table is required. Katana does not merge nearby coordinates with an arbitrary base-pair window.

Discovery module eigengene loadings are fitted in the reference cohort and projected into independent target PSI matrices. Condition effects are represented as Hedges g and combined with DerSimonian-Laird random-effects meta-analysis (Q, I2, and tau2 are reported).

If independently discovered target module assignments are supplied, modules are matched one-to-one using the validated shared-event universe. Match evidence is tested with a hypergeometric test and BH-FDR. Across multiple targets, each reference module is labelled `preserved_across_cohorts`, `partially_preserved`, `not_preserved`, or `not_assessed` from the configured supported-cohort fraction. A module may instead be labelled `rewired_across_cohorts` or `partially_rewired` only when independent preservation diagnostics show negative connectivity with permutation FDR support (or empirical p-value when no adjusted value is supplied). Rewiring is never inferred from Jaccard or low overlap alone. Without target assignments Katana reports projected effects and `not_assessed`; it does not invent structural states.

The structural outputs are `cross_cohort_module_state.tsv` and `multi_cohort_module_consensus.tsv`. Their labels describe recoverability in configured independent cohorts and are not causal claims.

## Cached validation report and browser

Stage 7 reads cached artifacts only and writes a reproducible HTML summary, a PDF validation summary, and a self-contained module browser. The browser exposes assignments, hubs, regulators, enrichment, isoform switches, and preservation tables without rerunning discovery. Reports include parameters, source audit, input fingerprints, and explicit unavailable/not-assessed states.

## Head-to-head boundaries

1. Stage 6C compares CoSplice, exact WGCNA (when R is available), hierarchical correlation, and graphical lasso on one predeclared stable-AS event universe.
2. Stage 6F compares beta-binomial replicate/read-count evidence with simple delta PSI on the same events.

Internal bootstrap stability is not external validation. External preservation (6D) and cross-cohort meta-analysis (6G) are required before claiming generalization.

## Primary method references

- rMATS replicate-aware model: https://pmc.ncbi.nlm.nih.gov/articles/PMC4280593/
- LeafCutter Dirichlet-multinomial model: https://pmc.ncbi.nlm.nih.gov/articles/PMC5742080/
- DRIMSeq DTU/DEU framework: https://pmc.ncbi.nlm.nih.gov/articles/PMC5200948/
- Beta-binomial splice-count modeling: https://pmc.ncbi.nlm.nih.gov/articles/PMC8144374/
- rMAPS2 positional RBP motif enrichment: https://pmc.ncbi.nlm.nih.gov/articles/PMC7319468/
- ATtRACT RBP motifs: https://pmc.ncbi.nlm.nih.gov/articles/PMC4823821/
- CisBP-RNA/EuPRI: https://pmc.ncbi.nlm.nih.gov/articles/PMC12807593/
- WGCNA module eigengene and kME: https://pmc.ncbi.nlm.nih.gov/articles/PMC2631488/
- WGCNA module preservation: https://pmc.ncbi.nlm.nih.gov/articles/PMC3024255/
- Consensus module hubs: https://pmc.ncbi.nlm.nih.gov/articles/PMC3629234/
