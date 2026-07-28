# Native 5.6 vs bundled Jutils: strict same-input validation

## Scope

This validation used the same Setd2 project, the same six selected rMATS JC
comparisons, and the same practical filters:

- `rMATS FDR <= 0.05`
- `abs(dPSI) >= 0.20`

Legacy Jutils and Katana Native 5.6 were run independently. This is a
same-input functional comparison, not a claim that the figures must be
pixel-identical or that either method is universally statistically superior.

## Filtered-event agreement

| Comparison | Bundled Jutils | Native 5.6 | Agreement |
| --- | ---: | ---: | --- |
| CLK3 WT DP | 5,603 | 5,603 | Exact |
| CLK3 WT noDP | 5,706 | 5,706 | Exact |
| MUT DP | 484 | 484 | Exact |
| MUT noDP | 318 | 318 | Exact |
| Pla WT DP | 1,688 | 1,688 | Exact |
| Pla WT noDP | 1,444 PCA events / 1,443 heatmap rows | 1,444 | Exact before the legacy heatmap duplicate-feature collapse |

The one-row Pla heatmap difference is expected: bundled Jutils collapses
duplicate `GeneName + FeatureLabel + FeatureType` rows for the supervised
heatmap. Its PCA retains all 1,444 filtered events, matching Native 5.6.

## Runtime and artifacts

| Metric | Bundled Jutils | Native 5.6 |
| --- | ---: | ---: |
| Wall time | 346.371 s | 56.736 s |
| Relative runtime | 1.0x | about 6.1x faster |
| Generated files | 76 | 47 |
| Plot formats | 55 PDF, 1 PNG | 28 PNG |
| Reusable tables | 20 TSV | 17 TSV |
| Provenance | converted inputs and command outputs | `parameters.json`, tables, plots and benchmark report |

## Output semantics

Both paths provide sample-level PSI heatmaps and PCA. The complete Katana
workflow also provides:

- reusable absolute-PSI and row-Z-score matrices;
- per-comparison and cross-comparison sample PSI PCA;
- replicate PSI QC and within-group versus between-group separation;
- exact-event overlap, signed direction agreement, Jaccard and
  denominator-aware confidence weighting;
- low-overlap warnings that prevent a high correlation based on very few
  shared events from being overinterpreted;
- direction-aware Venn analysis in 5.4.5;
- selected-event Native Sashimi in 5.5.

Bundled Jutils remains useful as a compatibility/reference implementation. It
retains classic PDF layouts, multi-tool converted input handling and arbitrary
MDS metric support. Native 5.6 deliberately applies `Top N` per comparison for
readable previews and uses PCA/SVD rather than reproducing every Jutils MDS
option.

## Conclusion

For the tested rMATS pattern-exploration scope, Native 5.6 matched the filtered
event population and supplied equivalent or stronger auditable outputs while
running substantially faster. The conclusion is limited to this validated
scope; it is not a blanket claim of superiority over every Jutils input mode or
command-line option.

Reference: [Jutils publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC9502236/).
