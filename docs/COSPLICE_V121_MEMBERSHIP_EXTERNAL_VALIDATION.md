# CoSplice v1.2.1 membership-confidence validation

## Question

Does the v1.2.1 split-half core label select events that transport better than
a simple WGCNA hub-event filter?

The fair primary comparison is not core events versus every WGCNA member.
Smaller event sets can look cleaner by construction. The validation therefore
uses three arms on one frozen WGCNA partition:

1. all non-grey WGCNA members;
2. the top absolute-kME members, selected per module and count-matched to the
   CoSplice core arm;
3. CoSplice events confirmed as core in both bootstrap split halves.

## Public panel

The runner consumes prepared TCGA SpliceSeq PSI matrices. A discovery cohort is
frozen before the runner opens the external cohorts. Public SpliceSeq does not
provide replicate IJC/SJC, so this benchmark validates sample-bootstrap
membership confidence but not the separate read-count uncertainty layer.

Primary endpoints are external WGCNA-like module preservation and program-score
preservation. Secondary endpoints are recovery of GO/KEGG themes in an
independently fitted holdout partition and projection of discovery program
weights into tumor/normal samples.

## Reproduce

```powershell
$env:PYTHONPATH = 'E:\Project\katana'
python scripts\validate_membership_confidence_panel.py `
  --panel-root 'E:\Project\RNASEQ\Public_validation\TCGA_SpliceSeq_ES_10_cohort_panel_v1' `
  --output 'E:\Project\RNASEQ\Public_validation\CoSplice_v121_BRCA_membership_external_validation' `
  --discovery BRCA `
  --holdouts KIRC KIRP KICH LUAD LUSC `
  --event-type ES `
  --gmt 'E:\Project\RNASEQ\Public_validation\resources\GO_Biological_Process_2023.gmt' `
  --gmt 'E:\Project\RNASEQ\Public_validation\resources\KEGG_2021_Human.gmt' `
  --rscript 'D:\anaconda3\envs\katana-splicing-tool\Scripts\Rscript.exe' `
  --events 500 --bootstraps 20 --permutations 200 `
  --min-module-size 15 --random-seed 20260901
```

## Current result

The BRCA discovery partition contained three non-grey modules and was evaluated
in five external cohorts. Against count-matched top-kME members, CoSplice core
members had a small positive paired median difference for Zsummary and
program-score preservation, but neither endpoint was significant after FDR
control. No GO/KEGG theme met the independent replication definition. The
confidence layer is therefore operational and scientifically controlled, but
external superiority is not established by this panel.

The KIRC discovery sensitivity run produced only one non-grey module under the
same default settings. It is retained as a negative control and must not be used
to claim module-level functional replication.

## Interpretation

This validation closes the implementation gap, not the methods-paper evidence
gap. A methods claim requires a predeclared disease-matched cohort with multiple
recoverable modules, read-count data for the uncertainty arm, and significant
improvement over the matched top-kME baseline. Cross-cancer transportability is
supporting evidence only.
