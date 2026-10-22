# Bulk thresholds: values, rationale and evidence

This document covers the existing bulk differential-result and replicate-QC defaults. It distinguishes published methods, published application choices and Katana-specific heuristics. A reference for a statistical method does not validate every numerical cutoff layered on its output.

The companion [machine-readable register](thresholds/bulk-defaults.v1.json) records each value, operator, unit, evidence classification and source. Defaults remain configurable where the corresponding analysis exposes them. This documentation does not change an algorithm, an existing run or an upstream nf-core parameter.

## Differential-result selection

| Analysis | Existing selection rule | Basis and limits |
| --- | --- | --- |
| DAS | Native rMATS FDR **≤0.05** and observed **absolute ΔPSI ≥0.10** | Published rMATS applications use this general FDR/effect-size combination [3]. Katana's inclusive boundary operators are explicit implementation choices. This is downstream selection of imported statistics, not a replacement for the rMATS model [2]. |
| DEG | Native adjusted P **≤0.05** and observed **absolute log2FC ≥1.0** | Configurable significance and effect-size choices. DESeq2 provides the statistical method [5], but does not prescribe this pair as its universal default. Its documented `results()` defaults include `alpha=0.1`; its thresholded testing is separate from observed-effect selection [6]. |
| DTU | Corresponding imported adjusted value **≤0.05** | Preserve the original gene/transcript endpoint and adjustment. For stageR, the selection target must match the target used for stage-wise adjustment [8,9]; see below. |
| DEXSeq | Corresponding imported adjusted value **≤0.05** | A project selection level applied to exon/bin evidence, not a fixed DEXSeq requirement. The original exon-usage method and its testing unit are retained [7]. |

BH provides a multiple-testing procedure [1]; the choice of 0.05 is an analysis policy, not a universal biological boundary. Katana does not turn every imported adjusted statistic into a new BH analysis or pool gene-, transcript- and exon-level values.

An absolute ΔPSI of 0.10 means ten percentage points. An absolute log2FC of 1 corresponds to a twofold ratio in either direction. These are observed-effect filters: passing them does not establish a confidence bound on the true effect.

### Native tests and later filtering

- **rMATS:** retain the native test settings, including `cstat`, in upstream provenance [2,10]. Filtering an observed ΔPSI at 0.10 does not mean the original null hypothesis tested an effect of 0.10.
- **DESeq2:** changing a Katana DEG cutoff reclassifies imported rows. It does not rerun independent filtering, shrinkage or a thresholded-effect test. A post hoc fold-change filter must not be described as a native test proving that the true absolute log2FC exceeds 1 [5,6].
- **stageR:** the author's DTU workflow uses an adjustment target of 0.05 and explains that another target requires repeating the stage-wise adjustment [9]. This usually does not require requantifying reads. Existing stageR values must not be treated as arbitrary-target q-values. The current generic importer does not enforce this compatibility check; recording and validating the original adjustment target is a packaging follow-up.
- **DEXSeq:** an exon/bin-level adjusted P-value is not automatically a new gene-level FDR when a gene is selected because it has a significant bin [7].

Additional read or consistency filters do not automatically confer the original model's FDR guarantee on the resulting subset. Complete source results and their native statistics remain available.

## Read support and replicate consistency

| Quantity | Existing balanced rule | Evidence classification |
| --- | --- | --- |
| Per-sample informative count | **IJC + SJC ≥10** | **Published application precedent.** Kuijper et al. require at least ten informative reads in every sample [4]. This supports the numerical read-depth choice, not every Katana aggregation rule. |
| Finite-PSI and read-supported replicates | In **each comparison arm**, at least **2** finite-PSI samples and at least **2** samples individually meeting the read-count rule | **Katana aggregation policy with related published practice.** Published rMATS studies include two-replicate designs [3], but this does not validate two as a sufficient experimental design or validate the complete Katana `10 × 2` rule. |
| Cross-arm direction consistency | **≥0.75** | **Katana descriptive heuristic.** No directly matching primary validation for this exact statistic and cutoff was identified in this review. |
| Leave-one-out direction consistency | **≥0.75** | **Katana descriptive heuristic.** Not an independently replicated result or a refitted rMATS test. No direct validation of this exact cutoff was identified. |
| Approximate PSI-SE gate | **Disabled by default**; diagnostic reference **0.12** | **Uncalibrated diagnostic.** Neither 0.12 nor the current approximation is claimed to be a validated rMATS uncertainty cutoff. |

### What the read-count rule actually does

The ten-read criterion is evaluated for each sample using its paired `IJC_SAMPLE_*` and `SJC_SAMPLE_*` entries, not a group mean or a sum over replicates. For JC inputs these are junction-count quantities; JCEC includes the caller's additional exon-count evidence. The count mode is retained and must be reported. The cited application [4] filters junction-count results; applying the same ten-count cutoff to JCEC is an explicit Katana extension, not directly validated by that example.

Kuijper et al.'s published rule applies to **every sample** [4]. Katana's balanced rule instead requires **at least two supported samples in each arm**. Other samples may have fewer than ten counts. These policies are different, so the full Katana rule is documented as a local policy rather than attributed verbatim to that publication.

Two observed replicates allow this rule to be evaluated; they do not establish adequate power, biological replication quality or eligibility for every statistical design. Repeated observations from one biological unit must not be promoted to independent biological replicates.

### Exact definitions of the two consistency metrics

Both use the finite PSI observations, **not only samples passing the ten-count rule**. Let the expected direction be the sign of the imported, direction-mapped ΔPSI.

1. **Cross-arm metric:** compare every finite PSI in the numerator arm with every finite PSI in the denominator arm. Remove equal-value pairs. The metric is the fraction of remaining differences with the expected sign. These pairwise differences share samples and are not independent replicates.
2. **Leave-one-out metric:** remove one finite PSI observation from one arm at a time while keeping the other arm intact. Compute the new arithmetic mean difference, omit zero differences, and count the fraction retaining the expected sign. This does not refit the native model and does not implement paired-donor leave-one-out inference.

If the imported effect is zero/non-finite, or there are no eligible nonzero differences, the metric is unavailable and its consistency check does not pass. The current code does not enforce a separate finite-PSI range validation at this loader stage; documentation must not imply that it does.

The 0.75 values are interpretable screening choices, not P-values, posterior probabilities or experimentally established error rates. Their exact numerical calibration remains unestablished; sensitivity results must be distinguished from independent validation. The reviewed primary sources support investigating depth and replicate variability [2], but do not validate these two exact 0.75 checks.

### Approximate uncertainty

For each finite PSI with positive finite count `n = IJC + SJC`, the implementation calculates:

```text
approximate_SE = sqrt(max(PSI * (1 - PSI), 0) / (n + 1))
event_diagnostic = maximum approximate_SE across both arms
```

The extra `+1`, the maximum aggregation and the reference value 0.12 are Katana diagnostic choices. This is not the rMATS hierarchical-model standard error: the native model includes effective-length transformations and replicate variability [2]. It is not a calibrated confidence interval or a biological replicate standard error.

`apply_measurement_se_gate` is **false** by default. Explicitly enabling it applies the chosen preset's ceiling; that remains an experimental filter, not a literature-validated stability criterion.

## Presets, result states and effective parameters

The existing presets are project policies:

| Preset | Counts per supported sample | Supported samples per arm | Both consistency fractions | SE reference, only if explicitly enabled |
| --- | --- | --- | --- | --- |
| exploratory | ≥8 | ≥2 | ≥0.67 | ≤0.18 |
| **balanced** | **≥10** | **≥2** | **≥0.75** | **≤0.12** |
| stringent | ≥20 | ≥3 | ≥0.80 | ≤0.10 |

The exact preset combinations are not represented as published, universally validated standards. `0.67` is the numerical boundary in the code, not exactly two-thirds. Select `stability_profile="custom"` to apply user-entered QC values; other presets use their stored values even if custom-value fields are supplied.

The output distinguishes:

- `measurement_eligible`: sufficient finite observations, read support and the optional SE check; no direction/effect requirement.
- `differential_repeat_supported`: both direction-consistency checks.
- `stable_event`: both preceding states.
- `stable_significant`: `stable_event` plus the imported FDR/observed-effect selection.

Complete event tables retain failures and exclusion reasons. Descriptive consistency filtering must not define the entire measurable feature universe for a co-splicing network.

DAS can inherit a selected or newest compatible completed Quality result. In that case its effective policy comes from that result, rather than reapplying the visible balanced defaults. Record the selected source, preset, effective `criterion_*` values, parameters, JC/JCEC mode, sample/contrast identities and method versions with the output. An old report does not silently change when another run is created.

## Source and documentation status

Implementation anchors: [replicate metrics and SE](../katana_next/science/bulk.py#L483), [presets and QC](../katana_next/science/bulk.py#L557), [native significance and stable subsets](../katana_next/science/bulk.py#L654), [DEG runtime selection](../katana_next/workflows/bulk.py#L594), and [Quality-result inheritance](../katana_next/workflows/bulk.py#L128).

The related default values were retained during review on 2026-10-22. This is an explicit parameter policy and evidence audit, **not a new validation experiment**. Thresholds for other modules will be documented separately as they are reviewed. No source statistics, scientific code, frozen original installation or completed research results were changed for this documentation.

## References and what they support

1. Benjamini and Hochberg (1995), *Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing*. [DOI](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x). Supports the multiple-testing procedure, not a mandatory 0.05 choice or a guarantee for arbitrary filtered subsets.
2. Shen et al. (2014), *rMATS: Robust and flexible detection of differential alternative splicing from replicate RNA-Seq data*. [DOI](https://doi.org/10.1073/pnas.1419161111). Supports native modeling of sampling uncertainty and between-replicate variability; does not validate Katana's 75% or approximate-SE cutoffs.
3. *STRAP regulates alternative splicing fidelity during lineage commitment of mouse embryonic stem cells* (2020). [DOI](https://doi.org/10.1038/s41467-020-19698-6). A primary rMATS application using FDR/effect selection and two-replicate comparisons; not validation of every Katana rule or a sample-size recommendation.
4. Kuijper et al. (2025), *Determining off-target effects of splice-switching antisense oligonucleotides using short read RNAseq in neuronally differentiated human induced pluripotent stem cells*. [DOI](https://doi.org/10.1093/hmg/ddaf153). Methods: differential splicing; per-sample informative-read filtering. Its all-sample rule differs from Katana's per-arm supported-replicate count.
5. Love, Huber and Anders (2014), *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2*. [DOI](https://doi.org/10.1186/s13059-014-0550-8). Native modeling and effect-threshold hypothesis tests; a later effect-size filter is not that test.
6. DESeq2 authors, [official vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html). Independent filtering, `alpha` and thresholded-effect testing. Consult the exact upstream version when reproducing a result.
7. Anders, Reyes and Huber (2012), *Detecting differential usage of exons from RNA-seq data*. [DOI](https://doi.org/10.1101/gr.133744.111); [official DEXSeq guide](https://bioconductor.org/packages/release/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html). Supports exon/bin-level modeling and interpretation.
8. Van den Berge et al. (2017), *stageR: a general stage-wise method for controlling the gene-level false discovery rate in differential expression and differential transcript usage*. [DOI](https://doi.org/10.1186/s13059-017-1277-0). Supports the method's gene-level overall error-rate target.
9. Love, Soneson and Patro, [author-maintained rnaseqDTU workflow](https://www.bioconductor.org/packages/release/workflows/vignettes/rnaseqDTU/inst/doc/rnaseqDTU.html), stageR sections. The adjustment target is fixed for downstream interpretation; changing it requires repeating the adjustment.
10. Xing laboratory, [rMATS-turbo documentation](https://github.com/Xinglab/rmats-turbo#readme). Native options and count definitions; retain the exact release/configuration used upstream.

Primary papers and author documentation were checked on 2026-10-22. A failure to find a directly matching validation is recorded as a limitation, not proof that no such paper can exist.
