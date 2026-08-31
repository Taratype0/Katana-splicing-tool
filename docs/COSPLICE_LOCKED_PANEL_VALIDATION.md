# CoSplice locked validation panel

## Question

This validation asks whether CoSplice recovers stable, non-trivial event programs that transfer to
an untouched cohort better than a same-input exact R WGCNA baseline. It does not assume that more
modules, a higher ARI, or one enriched term alone proves superiority.

## Locked design

- Discovery: TCGA SpliceSeq KIRC.
- Tuning: TCGA SpliceSeq KIRP.
- Primary same-organ holdout: TCGA SpliceSeq KICH.
- Transportability challenges: LUAD and LUSC.
- Event classes: ES, RI, AA, AD, AP, AT and ME.
- Scale sensitivity: 250, 500 and 1,000 ES events where available.
- Comparators: stability-only CoSplice, adaptive multi-objective CoSplice and exact R WGCNA.
- Selection: discovery and tuning only; holdout files are read after assignments and SHA256
  fingerprints have been frozen.
- Functional source: a predeclared local GO/KEGG GMT. Functional evidence affects selection only
  after its formal FDR activation gate passes.

Each child output contains the analysis plan, pre-holdout freeze manifest, input/assignment hashes,
holdout access audit, assignments, preservation statistics, enrichment tables, report and figure.
The aggregate command produces one TSV, PNG, Markdown report and parameter manifest.

## Current result

Across 11 frozen runs, CoSplice produced a valid multi-module partition in 9/11 settings, compared
with 2/11 for exact WGCNA. Mean holdout program-score preservation was 0.534 for CoSplice and 0.483
for WGCNA. WGCNA had higher mean bootstrap ARI (0.883 versus 0.705), and mean traditional
Zsummary-like preservation did not favor CoSplice (0.247 versus 0.139).

No run contained enough formally significant tuning-GMT modules to activate the functional
objective. Adaptive multi-objective CoSplice therefore correctly reduced to the stability and
preservation objective and selected the same partitions as stability-only CoSplice in 11/11 runs.
No independent GO/KEGG term replication was observed. These are retained negative results.

## Defensible claim

The current evidence supports a narrower architectural advantage: CoSplice more often resolves
non-trivial event-level program structure and has modestly higher program-score preservation on this
panel. It does not support universal superiority over WGCNA, superior internal stability, or
superior functional replication. TCGA cohorts are independent patient cohorts but share one data
generation and processing ecosystem; a cross-study, cross-platform holdout remains necessary for a
strong methods-paper generalization claim.

## Commands

```powershell
katana-cosplice-locked-ablation --help
katana-cosplice-locked-panel-report --root E:\path\to\completed_locked_runs
```

The TCGA SpliceSeq downloader accepts repeated `--splice-type` arguments for `ES`, `RI`, `AA`,
`AD`, `AP`, `AT` and `ME`, and merges incremental downloads into one panel manifest.
