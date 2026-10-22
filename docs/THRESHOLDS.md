# Thresholds: current policies and evidence

**Active policy: `bulk-all-samples-counts-v1`.** For each selected contrast, every expected sample in both arms must have **finite informative counts `IJC + SJC >= 10`**. This is a fixed per-sample rule. Counts are not pooled across samples, arms or contrasts.

The [machine-readable register](thresholds/bulk-defaults.v1.json) records this policy separately from native significance selection. The former replicate-consistency presets and approximate-uncertainty options are superseded; their exact documentation is retained only in the [historical archive](thresholds/history/8c368b7/README.md). They are not current choices or additional gates.

## AS read support

The Quality workflow is **AS read support**. It evaluates one declared numerator/denominator contrast at a time:

1. Resolve the expected sample IDs and count-slot mapping for both arms from the project's sample and contrast metadata. Both arms must be nonempty. An unknown, duplicated or ambiguous mapping cannot be assumed correct.
2. Require actual IJC, SJC and PSI observations for **every expected sample**, with the same vector length. Counts must be finite and nonnegative, and PSI must be finite and within its defined range of 0 to 1. Missing fields, unequal/incorrect array lengths, unmapped samples and malformed values cannot pass. These are input-validity checks, not extra stability thresholds. A missing arm or empty list never passes by vacuous truth.
3. Evaluate `IJC + SJC >= 10` separately for every expected sample. A valid zero count fails the threshold; an unavailable count remains unavailable rather than being silently converted to zero.
4. Pass the event only when every participating sample in both arms passes. A group mean, total over replicates, subset of supported samples or evidence borrowed from another contrast cannot satisfy the rule.

Only samples participating in the selected contrast are evaluated for that contrast. Selecting several contrasts does not pool their sample counts. Lane consolidation and technical-measurement handling belong to upstream sample preparation; lanes are not promoted to additional biological replicates.

Insufficient or ambiguous input is blocked or recorded as unavailable with a reason, as appropriate to the input error. It is never reported as passing support. A computable below-threshold result is a read-support failure, distinct from missing data. Complete result tables retain the relevant state and exclusion reasons.

This policy adds **no replicate-consistency or approximate-SE filter**. It does not establish that a study has adequate replication, sufficient power or a valid statistical design. Native model eligibility and experimental-unit requirements remain separate.

### Literature basis and its limits

Kuijper et al. [4] filtered **junction-count results** by requiring at least ten informative reads (`IJC + SJC`) in **each sample**. This is a published precedent for the active per-sample JC policy, rather than a claim that ten is optimal for every assay or experiment.

For JCEC inputs, the same numerical rule is an explicit Katana extension to the caller's additional exon-count evidence. That paper does **not** directly validate a JCEC threshold of ten. Record JC/JCEC mode with every result; do not present the two count definitions as identical.

The read-support gate does not refit rMATS, change effective lengths, recalculate native P-values or replace the model's treatment of replicate variation [2]. A support-filtered subset does not automatically inherit a new, independently established error-rate guarantee.

## Differential-result selection

| Analysis | Current default selection | Basis and limits |
| --- | --- | --- |
| DAS | Native rMATS FDR **<=0.05** and observed **absolute delta PSI >=0.10** | Published rMATS applications use this general combination [3]. Inclusive boundary operators are explicit implementation choices. This is downstream selection of imported statistics, not a replacement for the rMATS model [2]. |
| DEG | Native adjusted P **<=0.05** and observed **absolute log2FC >=1.0** | Runtime-configurable significance and observed-effect choices. DESeq2 does not prescribe this pair as its universal default [5,6]. |
| DTU | Corresponding imported adjusted value **<=0.05** | Preserve the original gene/transcript endpoint and adjustment. For stageR, the target must match the original stage-wise adjustment target [8,9]. |
| DEXSeq | Corresponding imported adjusted value **<=0.05** | A project selection level on exon/bin evidence, not a fixed DEXSeq requirement or automatically a new gene-level FDR [7]. |

An absolute delta PSI of 0.10 is ten percentage points; an absolute log2FC of 1 is a twofold ratio in either direction. These are observed-effect filters, not confidence bounds on the true effects. BH provides a multiple-testing procedure [1]; choosing 0.05 is an analysis policy, not a universal biological boundary.

DAS read support and significance are separate states. The read gate is fixed by `bulk-all-samples-counts-v1`; the default significance settings remain independently identified. DEG's existing runtime parameters `fdr` and `min_abs_log2fc` remain available as described in [usage](USAGE.md#deg-thresholds-at-runtime).

### Native statistics remain native

- **rMATS:** retain upstream test settings, including `cstat` [2,10]. Filtering an observed delta PSI at 0.10 does not mean the original null hypothesis tested a difference of 0.10.
- **DESeq2:** the documented `results()` defaults include `alpha=0.1` and `lfcThreshold=0` [6]. Changing a Katana DEG selection reclassifies imported rows; it does not rerun independent filtering, shrinkage or a thresholded-effect test. An observed fold-change filter must not be called a native test proving the true absolute log2FC exceeds 1 [5,6].
- **stageR:** the author's DTU workflow uses an adjustment target of 0.05 and explains that a different target requires repeating stage-wise adjustment [9], usually without requantifying reads. Existing stageR values are not arbitrary-target q-values. **The current generic importer does not enforce target compatibility**; recording and validating the original target remains a packaging follow-up.
- **DEXSeq:** an exon/bin adjusted P-value does not become a new gene-level FDR when a gene is selected for having a significant bin [7].

Katana retains complete source results and original effects, P-values and adjusted statistics. It does not pool gene-, transcript- and exon-level adjustments or recalculate them after the read gate.

## Output compatibility and reproducibility

Read support, native significance and their intersection must remain distinguishable. The legacy column names `stable_event` and `stable_significant` are retained for compatibility only: under this policy they mean **the all-sample count gate** and **that gate intersected with native significance/effect selection**, respectively. They no longer imply replicate-direction consistency, uncertainty calibration or biological stability.

The result records these distinct fields:

| Field | Meaning |
| --- | --- |
| `read_support_policy` | Active policy identity. |
| `passes_data_usability` | Both arms have valid, correctly sized input vectors and valid PSI/count observations. |
| `passes_read_counts` | Every expected count slot in both arms meets the fixed ten-count rule. |
| `passes_read_support` | Both input usability and the all-sample read-count rule pass. |
| `measurement_eligible` / `stable_event` | Compatibility aliases for `passes_read_support`. |
| `stable_significant` | Read support intersected with native significance and observed-effect selection. |
| `read_support_failure_reason` | Explicit invalid-input or low-count reasons. |
| `numerator_expected_samples` / `denominator_expected_samples` | Expected arm sizes resolved for the selected contrast. |

Record the policy ID, per-sample count threshold, expected sample IDs/counts, resolved arm/slot mapping, JC/JCEC mode, support/availability reasons, significance settings, native test settings and method/workflow versions. A new policy requires a new analysis/cache identity; historical completed Quality results must not be silently relabeled as current-policy results.

Accepted legacy parameter keys, when retained for call compatibility, are inert and recorded as ignored; they cannot turn the former gates back on. Retired non-default settings are not an alternative current policy.

### Selected-contrast Landscape matrices

When Landscape assembles an event-by-sample matrix for several selected contrasts, a passing row in one contrast cannot rescue a failing row in another. Every selected sample used in that matrix must have valid PSI and read support of at least ten, and all included rows for that event from the selected contrasts must be `measurement_eligible`. The matrix is not built by an OR across contrast-level support decisions. Missing or unassessable sample values are not zero counts or passing measurements.

The selected-contrast rule and the whole-cohort network rule have different sample scopes. The network policy below evaluates every sample participating in the fitted network; a passing subset of contrasts is not sufficient. Opening an old report does not apply either new policy or recompute its upstream analysis.

## Whole-cohort network input

**Approved policy: `network-all-samples-counts-v1`.** For count-bearing network inputs, use the measured, nonconstant event universe of the declared network cohort. Do **not** require that an event is a significant DAS hit or has a particular condition effect before unsupervised network discovery.

| Current requirement | Purpose | Basis and limits |
| --- | --- | --- |
| Every declared participating sample has finite PSI within 0 to 1 and actual finite, nonnegative IJC/SJC | Reject incomplete, invalid or unmapped measurements | Input validity, rather than an additional biological effect or significance cutoff. |
| `IJC + SJC >= 10` in **every** participating sample | Establish per-sample count support across the entire fitted cohort | Same per-sample JC precedent [4] as AS read support; applying the number to JCEC or a different count-producing assay is an explicit Katana extension. |
| At least two distinct PSI values in an event's complete sample vector | Exclude constant features for which correlation cannot be defined | Mathematical eligibility; no positive variance floor such as 0.0025 remains. It does not ensure reliable correlations or sufficient sample size. |
| No DAS-significance prefilter | Retain unsupervised co-variation that need not follow a prespecified condition effect | WGCNA authors advise against selecting differential-expression hits before network analysis [11]. Extending this reasoning to splicing events is declared, rather than a claim that this entire PSI policy was validated by WGCNA. |

All samples supplied to a given network fit must pass individually for an event to enter that fit. Count support cannot be established with a cohort average, pooled counts, at least two supported replicates per condition, or a fraction of supported samples. The declared cohort, sample IDs, input mode and count evidence must be recorded.

The existing rMATS matrix adapter chooses one canonical source table for each `(condition, event type)` before this gate: largest table by file size, then lexical path/side order. It evaluates the complete PSI/IJC/SJC observations in that chosen source, rather than assembling partial counts from different comparisons or taking median/maximum counts. Its source table records the chosen file, selection rule and candidate count. **This adapter does not certify agreement among the unselected candidate tables.** If the same sample/event has conflicting observations in another comparison, that discrepancy is outside this gate's evidence; a passing network event is not proof that every alternative source observation was audited. This source-selection limitation is preserved and disclosed, rather than silently changing the measurement adapter in this policy update.

The former 80% valid-PSI, 70% count-supported-sample, two-supported-replicates-per-condition and 80% condition-coverage hard filters are retired for this count-bearing policy. The former variance >=0.0025, repeatability >=0.30, condition-mean PSI span >=0.25 and strong dynamic-range gates are also retired as network-input eligibility filters. Such quantities may remain diagnostic columns, but they cannot exclude an otherwise complete, count-supported, nonconstant event under the active policy. This change does not remove a statistical model's own experimental-design requirements or change a native P-value.

### PSI-only input and count evidence

External PSI-only matrices remain a distinct optional input. Without actual counts, they must be labeled **`read_support_not_assessed`**; they cannot be declared to have passed the ten-count policy, have counts fabricated, or silently borrow a different assay's count definition. A network fitted from such a matrix must retain that limitation in its provenance and report. This approval does not authorize silently converting PSI-only input into the count-verified route.

The informative-count rule also does not alter raw BAM-junction discovery settings such as mapping quality, anchor length, event reconstruction or count definition. Those settings affect what was measured and must remain declared independently. Existing network methods, bootstrap procedures, association models and completed fits are not refitted or relabeled by this input-policy change.

### Approved module size and merging settings

The approved standard defaults are **minimum module size 20 splice events** and **merge cut height 0.25**. These are module-discovery settings, not read-support, P-value or preservation thresholds. Twenty refers to event nodes, not genes or samples; it is a clustering setting rather than a promise that every final table group meets that size. The merge cut applies to eigengene dissimilarity in the retained WGCNA method. The WGCNA authors' [tutorial](https://edo98811.github.io/WGCNA_official_documentation/FemaleLiver-02-networkConstr-man.pdf) uses 0.25; the [function documentation](https://search.r-project.org/CRAN/refmans/WGCNA/html/blockwiseModules.html) lists a default minimum of `min(20, ncol(datExpr)/2)` and a different default merge height (0.15). Thus 20/0.25 is Katana's explicit, configurable choice with methodological precedent, not a universal optimum for splicing.

Runtime `min_module_size` and `merge_cut_height` override the corresponding project YAML values for that run. Blank/null preserves the project configuration, including existing custom values. The standard engine fallback remains 20/0.25. Actual effective values are recorded in resolved configuration and artifact metadata, and displayed in the evidence report and the one-click interpreted report when its source network can be identified. Historical values without that provenance are **not recorded**, never backfilled from today's defaults. Module algorithms and completed fits are unchanged. See [module policy](thresholds/network-modules.v1.json) and [runtime examples](USAGE.md#network-module-settings).

### Continuous network evidence and optional labels

**Approved policy: `network-descriptive-labels-opt-in-v1`.** Continuous co-assignment,
kME, connectivity, margins, resampling diagnostics and complete rankings are the
default output. The Katana-specific core/peripheral, stability and hub-candidate
cutoffs are optional descriptive labels, disabled by default. Their old values
are documented in the [preset review](thresholds/network-confidence-review.md),
not presented as literature-calibrated boundaries or mandatory inclusion rules.
Enabling `apply_threshold_labels` adds labels without removing events or
changing their fitted module assignments. Disabling it produces `not_applied`
label states rather than false candidate failures or a claim that no hubs exist.

The old composite confidence score and preset bootstrap selection frequencies
depend on heuristic cutoffs. They remain traceable as legacy/preset sensitivity
evidence, not probabilities or required primary ranking scores. Numerical
measurements remain available independently of these optional labels. Invalid
inputs, undefined correlations, insufficient independent samples and unavailable
annotations retain explicit reasons; disabling labels does not manufacture
estimable measurements or override the separately approved count-support rule.

The default statistical reporting level stays **FDR 0.05**. BH correction retains
its native test families; full tested tables remain available. A significant
module association does not imply significant network edges or causal hubs.
Reports record the actual labeling switch and configured presets for their own
source run. Unknown historical settings are not backfilled. See the
[policy register](thresholds/network-labels.v1.json) for scope and source evidence.

The [network policy register](thresholds/network-input.v1.json) records this approved change separately from the [bulk contrast register](thresholds/bulk-defaults.v1.json). The default `max_events=0` includes all eligible events; a user-specified positive cap is an explicitly recorded exploratory subset ranked by cohort PSI variance, independent of condition significance. Implementation, testing and scientific validation are separate statuses; an implemented gate is not a new validation experiment.

The implemented network-input gate is versioned as **CoSplice engine `1.3.0` / workflow `programs.cosplice` `3.2.0`** within software `2.0.0b3.dev20261022`. Its register binds the current QC, stage, execution and workflow sources by SHA-256. Local input and execution-contract tests have passed; this does not claim a new native R network fit or biological validation. Resolved configuration and result provenance record the active policy and mark old QC settings as diagnostic-only. New engine/stage identities separate these runs from historical caches. Existing output lacking the current policy is rejected for overwrite, even when force is requested; completed fits retain their original identity.

## Implementation identity

This policy is implemented in software **`2.0.0b3.dev20261022`**. Its workflow versions are:

| Workflow | Version |
| --- | --- |
| `quality.stable_as` — AS read support | `3.0.0` |
| `differential.as` | `3.0.0` |
| `differential.cross_comparison` | `3.0.0` |
| `landscape.bulk` | `2.1.0` |

The [register](thresholds/bulk-defaults.v1.json) binds these identities to the current read-support core, workflow and source-contract guard hashes. These versions distinguish new-policy output and caches from historical Quality results. This is an implementation record, not a new biological validation claim or a refit of upstream statistical models. The stageR target-compatibility guard discussed above remains unimplemented and is recorded separately from the implemented read-support policy.

## Historical policy

Commit `8c368b7` documented the former presets and descriptive consistency/uncertainty calculations. Its exact Markdown and JSON bytes are [archived with hashes](thresholds/history/8c368b7/README.md). That archive explains old results only. It is neither an alternative current policy nor a source of active parameter defaults.

## References and what they support

1. Benjamini and Hochberg (1995), *Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing*. [DOI](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x). Supports the adjustment procedure, not a compulsory 0.05 choice or arbitrary-subset guarantees.
2. Shen et al. (2014), *rMATS: Robust and flexible detection of differential alternative splicing from replicate RNA-Seq data*. [DOI](https://doi.org/10.1073/pnas.1419161111). Supports native sampling/replicate modeling; not a replacement read-support model.
3. *STRAP regulates alternative splicing fidelity during lineage commitment of mouse embryonic stem cells* (2020). [DOI](https://doi.org/10.1038/s41467-020-19698-6). A primary application of FDR/effect selection, not a universal threshold validation or sample-size recommendation.
4. Kuijper et al. (2025), *Determining off-target effects of splice-switching antisense oligonucleotides using short read RNAseq in neuronally differentiated human induced pluripotent stem cells*. [DOI](https://doi.org/10.1093/hmg/ddaf153). Methods: junction-count events retained only when each sample had at least ten informative reads. This directly supports a numerical/per-sample JC precedent; the JCEC extension is Katana's declared choice.
5. Love, Huber and Anders (2014), *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2*. [DOI](https://doi.org/10.1186/s13059-014-0550-8). Native modeling and effect-threshold tests; a later observed-effect filter is different.
6. DESeq2 authors, [official vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html). Independent filtering, `alpha` and thresholded testing; preserve the exact upstream version and configuration.
7. Anders, Reyes and Huber (2012), *Detecting differential usage of exons from RNA-seq data*. [DOI](https://doi.org/10.1101/gr.133744.111); [official DEXSeq guide](https://bioconductor.org/packages/release/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html). Exon/bin modeling and interpretation.
8. Van den Berge et al. (2017), *stageR: a general stage-wise method for controlling the gene-level false discovery rate in differential expression and differential transcript usage*. [DOI](https://doi.org/10.1186/s13059-017-1277-0). Gene-level overall error-rate target.
9. Love, Soneson and Patro, [author-maintained rnaseqDTU workflow](https://www.bioconductor.org/packages/release/workflows/vignettes/rnaseqDTU/inst/doc/rnaseqDTU.html), stageR sections. Adjusted outputs are tied to their specified OFDR target.
10. Xing laboratory, [rMATS-turbo documentation](https://github.com/Xinglab/rmats-turbo#readme). Native options and count definitions; preserve upstream provenance.
11. WGCNA authors, [official FAQ](https://edo98811.github.io/WGCNA_official_documentation/faq.html), gene-filtering guidance. Advises against selecting genes by differential-expression significance before network construction; low-quality and low-variation filtering and its suitability remain data dependent. This is not validation of Katana's specific ten-count PSI rule or module parameters.

Policy and evidence reviewed on 2026-10-22. Documentation of this approved policy is not a new validation experiment. Previous installations, native statistical results and completed research deliveries remain preserved.
