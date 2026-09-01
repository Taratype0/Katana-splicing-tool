# CoSplice v1.2.2 continuous-confidence validation

## Question

Does a continuous, splicing-specific membership score select externally more
preserved module members than standard WGCNA eigengene membership (`kME`)?

The continuous score combines co-assignment confidence, bootstrap kME,
split-half confirmation and optional read reliability. It never changes the
posterior-mean WGCNA partition. The fair control is the same number of top-|kME|
events selected separately within every module.

## Frozen design

- Discovery: BRCA TCGA SpliceSeq exon-skipping PSI.
- Holdouts: KIRC, KIRP, KICH, LUAD and LUSC.
- Frozen discovery universe: 500 events, seed 20260901 and 20 stratified sample
  bootstraps.
- Four arms: all WGCNA members, count-matched top-kME, split-half core and
  count-matched continuous confidence.
- Formal endpoints: official `WGCNA::modulePreservation` and NetRep's seven
  preservation statistics, each with 1,000 permutations.
- Secondary diagnostics: the historical WGCNA-like score, projected program
  score, independently refitted holdout modules, GO/KEGG overlap and
  tumor-normal association.

Public SpliceSeq has PSI but no replicate IJC/SJC. This panel validates the
sample-bootstrap components of membership confidence. The read-count component
is evaluated separately in the Setd2 real-count ablation.

## Result

Fifteen module-cohort pairs were evaluable.

| Comparison | Result |
|---|---:|
| Official WGCNA median Zsummary, continuous confidence | 4.585 |
| Official WGCNA median Zsummary, count-matched top-kME | 4.560 |
| Paired median continuous-minus-kME Zsummary | -0.052 |
| Paired wins / losses | 5 / 10 |
| Paired Wilcoxon p-value / FDR | 0.041 / 0.083 |
| Median significant NetRep statistics, both arms | 7 / 7 |
| NetRep paired FDR | 1.0 |
| Replicated GO/KEGG themes | 0 |

The difference between the two medians above is not the paired median because
pairwise differences are calculated within the same module and cohort.

## Decision

Continuous confidence did not establish external superiority over count-matched
top-kME. It improved internal read reliability and co-assignment in Setd2 but
slightly reduced kME, and that trade-off did not improve formal cross-cohort
preservation. Katana therefore exports the score and all components for
prioritization and sensitivity analysis while keeping:

1. posterior-mean signed WGCNA/TOM as the production partition;
2. kME as the standard hub/member measure;
3. split-half role and continuous confidence as complementary audit fields.

This negative result prevents tuning a new weighted score on the same holdouts.
A future change to the score requires predeclared tuning cohorts followed by a
new untouched disease-matched validation cohort.

## Reproduce

```powershell
$env:PYTHONPATH = 'E:\Project\katana'
python scripts\validate_membership_confidence_panel.py `
  --panel-root 'E:\Project\RNASEQ\Public_validation\TCGA_SpliceSeq_ES_10_cohort_panel_v1' `
  --output 'E:\Project\RNASEQ\Public_validation\CoSplice_v122_BRCA_membership_external_validation' `
  --discovery BRCA `
  --holdouts KIRC KIRP KICH LUAD LUSC `
  --event-type ES `
  --gmt 'E:\Project\RNASEQ\Public_validation\resources\GO_Biological_Process_2023.gmt' `
  --gmt 'E:\Project\RNASEQ\Public_validation\resources\KEGG_2021_Human.gmt' `
  --rscript 'D:\anaconda3\envs\katana-splicing-tool\Scripts\Rscript.exe' `
  --events 500 --bootstraps 20 --permutations 200 `
  --formal-permutations 1000 --formal-workers 4 `
  --formal-cohort-workers 5 --formal-threads-per-arm 1 `
  --min-module-size 15 --random-seed 20260901
```

Primary files are `formal_preservation_summary.tsv`,
`formal_paired_confidence_vs_top_kme.tsv`,
`02_external_evaluation/official_wgcna_module_preservation.tsv`,
`02_external_evaluation/netrep_seven_statistic_preservation.tsv` and
`plots/official_wgcna_netrep_preservation.png`.
