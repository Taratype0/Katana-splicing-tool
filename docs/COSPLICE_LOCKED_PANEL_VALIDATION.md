# CoSplice locked validation panel

## Question and design

This panel asks whether a WGCNA/TOM-first CoSplice partition improves reproducibility or external
preservation over exact R WGCNA when both methods receive the identical discovery-selected event
universe. It does not treat a larger module count, a coarse one-module solution, or a high ARI alone
as proof of biological superiority.

The frozen panel contains 11 runs spanning AA, AD, AP, AT, ES, ME and RI events; ES event-universe
sensitivity at 250, 500 and 1,000 events; a 68-event ME stress test; and KICH, LUAD and LUSC
holdouts. KIRC is discovery, KIRP is tuning, and each holdout is inaccessible until TOM cut height,
consensus threshold and method parameters are frozen. KICH is a same-organ kidney holdout; LUAD and
LUSC are transportability challenges rather than same-disease replication.

Three arms are fitted:

1. stability-only CoSplice using the TOM backbone and bootstrap co-assignment consensus;
2. adaptive CoSplice, which may add predeclared tuning-cohort preservation and tuning-GMT
   functional coherence only when those objectives are evaluable;
3. exact R `WGCNA::blockwiseModules` on the same events and samples.

This is a frozen historical validation of the boundary-rejecting consensus
arm. The real-count ablation later showed that rejection reduced stability, so
v1.1.0 production now preserves the WGCNA partition and uses co-assignment only
as confidence annotation. The frozen results below are retained for provenance.

## Frozen result

Across the 11 runs, the final paired rerun used one seed and identical bootstrap
sample indices for every arm. TOM-first CoSplice reached WGCNA-level stability
but did not improve the panel mean:

| Metric | CoSplice | Exact WGCNA | Interpretation |
| --- | ---: | ---: | --- |
| Mean bootstrap ARI | 0.867 | 0.871 | 3 wins, 7 ties and 1 loss for CoSplice |
| Mean holdout program-score preservation | 0.497 | 0.472 | 3 wins, 7 ties and 1 loss for CoSplice |
| Non-empty module recovery | 11/11 | 11/11 | Equal recovery |
| Valid multi-module partition after strict coverage gate | 2/11 | 2/11 | Equal strict recovery |
| Functional objective activated | 0/11 | not applicable | No tuning-GMT module passed the predeclared gate |

The primary 500-event ES run was exactly tied: both methods returned one module,
mean bootstrap ARI 0.883, KICH Zsummary-like preservation 0.281 and
program-score preservation 0.582. This is parity, not improvement.

AT was the clearest preservation tradeoff: CoSplice recovered three modules with
bootstrap ARI 0.641 and median holdout program-score preservation 0.607, versus
four WGCNA modules with ARI 0.683 and preservation 0.332. Both passed the strict
partition gate. This is improved external score preservation at the cost of
lower internal stability, not an across-metric win.

## Claim boundary

Zero-module/all-grey partitions are excluded from ARI means and paired ARI wins because a constant
all-grey partition can otherwise receive a misleading perfect ARI. The frozen evidence supports
the statement that **a WGCNA/TOM backbone plus CoSplice bootstrap consensus
matches WGCNA-level resampling stability and modestly improves average external
program-score preservation on this panel**. It does not support universal superiority over WGCNA, better results
for every AS class, or stronger biological function. Zero GO/KEGG replication is retained as a
negative result and blocks a functional-superiority claim. The adaptive arm equals stability-only
whenever preservation cannot be evaluated or functional enrichment fails the predeclared FDR gate;
this is intentional rather than a missing analysis.

## Reproduction

```powershell
katana-cosplice-locked-ablation --help
katana-cosplice-locked-panel-report --root E:\path\to\completed_locked_runs
```

Every child run records the pre-holdout freeze manifest, assignment and input SHA256 fingerprints,
holdout-access audit, method-specific assignments, preservation statistics, local enrichment tables,
plots and a standalone report. The aggregator writes `LOCKED_ABLATION_PANEL_SUMMARY.tsv`,
`locked_ablation_panel.png`, `LOCKED_ABLATION_PANEL_REPORT.md` and
`locked_panel_parameters.json`.
