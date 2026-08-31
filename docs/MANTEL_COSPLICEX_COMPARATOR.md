# Strict Mantel-CoSplicEx Comparator

## Scientific question

The published CoSplicEx method is the closest formal peer for a gene-level co-splicing network.
It is not an event-level replacement for CoSplice:

- CoSplicEx represents one gene by its multi-exon count profile.
- CoSplice represents one AS/exon-usage event by its sample-level usage profile.

The comparator therefore uses the same public exon-count source and the same samples, but reports
the two native biological units separately. Assignment ARI between the methods is prohibited
because a gene node is not interchangeable with an event node.

## Published algorithm reproduced

The implementation in `katana_cosplice/mantel_cosplicex.py` follows Iancu et al. and the official
`iancuo/cosplicingNetworks` source:

1. Build one sample-by-exon matrix per multi-exon gene.
2. Calculate the Canberra distance between every sample pair for each gene.
3. Correlate lower-triangle distance vectors between genes. This is the Mantel statistic with zero
   permutations used by the original network-construction code.
4. Convert the correlation to an unsigned weighted adjacency with `abs(r)^6`.
5. Cluster `1 - adjacency` with average linkage.
6. Call `dynamicTreeCut::cutreeHybrid` with the published settings: `cutHeight=0.99995`,
   `maxCoreScatter=0.995`, `deepSplit=4`, `minClusterSize=20`, PAM enabled.

The public benchmark constrains both methods to a 500-feature budget for computational and
granularity comparability. Gene selection is phenotype-blind and based on multi-exon coverage;
CoSplice event selection is the predeclared variance-ranked exon-usage selection. This is a strict
algorithmic reproduction under a matched budget, not a claim that the feature selectors are
identical.

## Public result

Dataset: GSE62669, 60 mouse brain samples, 114,691 normalized exon rows from 9,066 genes.

| Endpoint | CoSplice | Mantel-CoSplicEx |
|---|---:|---:|
| Native feature unit | AS/exon-usage event | multi-exon gene |
| Features | 500 | 500 |
| Non-grey modules | 2 | 3 |
| Grey fraction | 0.486 | 0.000 |
| Median bootstrap ARI | 0.272 | 0.813 |
| Minimum bootstrap ARI | 0.177 | 0.304 |
| Functionally annotated module fraction | 0.500 | 0.667 |
| Significant GO/KEGG terms | 3 | 29 |
| RBP-enriched modules at FDR 0.05 | 0 | 0 |

Mantel-CoSplicEx recovered a 93-gene vesicle/synapse program and a 67-gene myelin program. These
are biologically coherent for brain tissue. The large 340-gene module had no comparable
FDR-significant functional theme and should not be over-interpreted.

CoSplice retained event identity but was less stable and its three significant annotations were
broad intracellular/organelle terms. On this dataset, event resolution did not compensate for the
weaker partition stability and weaker functional aggregation.

## Defensible conclusion

The strict comparator changes the publication claim:

- CoSplice is not universally superior to published Mantel-CoSplicEx.
- Mantel-CoSplicEx is stronger on this mouse-brain cohort for stable, gene-level biological
  programs.
- CoSplice provides a different advantage: specific event identity, read/uncertainty weighting,
  replicate-aware validation and external preservation.
- Those additions are scientifically useful only when the resulting event modules are stable and
  independently validated.

The publication should present Mantel-CoSplicEx as a formal complementary comparator, not omit it
and not force it into an incompatible event-level leaderboard.

## Reproducible command

```powershell
python scripts/benchmark_mantel_cosplicex.py `
  --exon-counts E:\Project\RNASEQ\Benchmark\public_sources\GSE62669\GSE62669_RNA110714RH_normalizedExonCounts.txt.gz `
  --gtf E:\Project\RNASEQ\gencode.vM34.annotation.gtf `
  --output E:\Project\RNASEQ\Public_validation\GSE62669_strict_Mantel_CoSplicEx_20260729 `
  --rscript "C:\Program Files\R\R-4.5.0\bin\Rscript.exe" `
  --genes 500 `
  --bootstraps 20 `
  --cosplice-summary E:\Project\RNASEQ\Public_validation\GEO_three_dataset_CoSplice_challenge_20260729\GSE62669_mouse_brain_exon_usage\internal_matched_bootstrap_benchmark.tsv `
  --cosplice-biology E:\Project\RNASEQ\Public_validation\GEO_three_dataset_CoSplice_challenge_20260729\three_dataset_biological_audit.tsv
```

Primary references:

- Iancu et al., *Frontiers in Genetics* (2015):
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC4429622/>
- Official CoSplicEx source:
  <https://github.com/iancuo/cosplicingNetworks>
- WGCNA:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC2631488/>
