# CoSplice Publication Benchmark Protocol

## Purpose

This protocol tests whether Katana CoSplice adds reproducible value beyond established network and
splicing methods. It does not seek one universal winner. The primary questions are:

1. Does a method recover planted modules when ground truth is known?
2. Is its event partition stable under matched sample resampling?
3. Are discovered modules preserved in an untouched external cohort?
4. Does a module recover known biological programs without using validation data during discovery?
5. What granularity, grey-event fraction and computational cost accompany that performance?

## Fair direct comparison

The final bootstrap-consensus CoSplice partition, exact `WGCNA::blockwiseModules`, Pearson
hierarchical clustering and graphical lasso receive the same discovery-selected event-by-sample matrix. Event selection
occurs once and is not repeated per method. Zero-module and one-module partitions are degenerate
and cannot earn a stability win.

Production CoSplice exposes one final bootstrap-consensus partition. The former tight/diffuse
hybrid remains benchmark-only legacy code and is not part of the UI or default pipeline.

The direct primary endpoints are matched-bootstrap adjusted Rand index, independent-cohort module
preservation and known-truth recovery where available. Silhouette, graph modularity,
within/between-module separation, module number, grey fraction, runtime and memory are secondary
trade-off metrics. Biological interpretability is a separate validation layer: unique-gene
diversity, protection against single-gene event domination, RBP/splicing-factor enrichment,
GO/KEGG enrichment, and cross-cohort module-gene/function replication. These annotations never
define modules and are not used for parameter tuning. No equal-weight composite is a scientific
endpoint.

## Published-method scope

Not every published tool is an interchangeable module partitioner:

- **Exact WGCNA** is the primary direct comparator on the same event-by-sample matrix.
- **WGCNA/Mantel co-splicing** compares gene-level multi-exon profiles and is reproduced as a
  separate strict gene-profile track using the published Canberra-distance, Mantel-correlation,
  unsigned power-6 adjacency and dynamic hybrid tree-cut algorithm.
- **CoSpliceNet** is an RBP/motif-centred time-course regulatory-network reference.
- **Spycone** clusters longitudinal transcript/isoform usage and belongs in a time-course track.
- **MAJIQ HET** evaluates event heterogeneity and robustness; it is not a module competitor.

This separation prevents a misleading leaderboard built from incompatible biological units.

## Current public evidence

- **GSE62669**: 60 mouse brain samples from the published WGCNA/Mantel co-splicing study. The
  transparent Katana adapter converts normalized exon counts into within-gene exon-usage fractions.
  This is an exon-usage proxy, not rMATS PSI and not the original Mantel statistic.
- **GSE141633 / TCGA-PRAD**: bulk public PSI scale-sensitivity track.
- **GSE85908**: single-cell PSI boundary test; not treated as bulk-RNA-seq replication.
- **TCGA SpliceSeq KIRC, KIRP and KICH**: rotating discovery with all six directed
  external-preservation projections.
- **Known-truth simulation**: positive control with three planted modules.

The formal GSE62669 1,000-feature, 20-bootstrap run produced two valid modules for each branch.
Exact WGCNA had higher mean bootstrap ARI (`0.568` versus `0.444`), cosine silhouette (`0.590`
versus `0.517`) and within/between separation (`0.862` versus `0.730`) than the tight branch.
The tight branch retained more events (grey fraction `0.449` versus `0.533`) and had slightly higher
neutral graph modularity (`0.336` versus `0.333`). The fixed Pearson hierarchical baseline was the
most stable partition on this dataset (mean bootstrap ARI `0.626`).

A separate strict Mantel-CoSplicEx comparator was then run on the same GSE62669 exon-count source
and the same 60 samples at a matched 500-feature budget. The biological units remain intentionally
different: CoSplice uses exon-usage events and Mantel-CoSplicEx uses multi-exon genes. Direct
assignment ARI between methods is therefore prohibited. Mantel-CoSplicEx produced three valid
modules with median/minimum bootstrap ARI `0.813/0.304`, zero grey fraction and 29 FDR-significant
GO/KEGG terms. CoSplice produced two modules with median/minimum ARI `0.272/0.177`, grey fraction
`0.486` and three broad FDR-significant terms. The Mantel modules included coherent
vesicle/synapse and myelin programs. This comparator demonstrates that CoSplice is not universally
superior: on this brain cohort the published gene-level method is more stable and more strongly
functionally aggregated, while CoSplice retains the distinct advantage of event-level resolution
and splicing-specific validation. See `docs/MANTEL_COSPLICEX_COMPARATOR.md`.

The formal KIRC/KIRP/KICH run used 500 discovery-selected ES events per cohort, 20 matched
bootstraps and 500 preservation permutations per module. CoSplice produced a valid multi-module
partition in all three discovery cohorts, whereas exact WGCNA produced one broad non-grey module in
each cohort and therefore had a valid-partition fraction of zero under the predeclared
`>=2 modules` and `largest module <=75%` rule. CoSplice's median and worst minimum-bootstrap ARI
were `0.350` and `0.171`. Its median and worst best-module external-preservation Zsummary were
`6.893` and `3.774`; all six directed projections contained at least one partially preserved
module. Exact WGCNA's corresponding preservation values were `-0.173` and `-0.992`, but its higher
bootstrap ARI describes a stable degenerate partition and is not counted as a scientific win.

The biological audit found ten CoSplice modules with a median of 94 unique genes per module,
median effective-gene fraction `0.992`, worst top-gene event fraction `0.05`, and no module
dominated by one multi-event gene at the predeclared 20% cutoff. Twenty-five GO/KEGG terms were
significant across the modules, with KIRC modules C02 and C04 supporting RNA-processing/nuclear
speck and membrane/cytoskeletal-adhesion themes, respectively. However, only 20% of modules had a
significant GO/KEGG annotation, no RBP-set enrichment survived FDR 0.05, and no significant
functional theme replicated across all three cohorts. Median best cross-cohort module-gene
Jaccard was `0.133`. The current evidence therefore supports non-trivial, externally preservable
event programs, but not a claim that their biological functions are already universally
replicated.

KIRC, KIRP and KICH are three independent disease cohorts from one public resource, TCGA
SpliceSeq. They satisfy a three-cohort rotation test, not a strict three-database validation.
Independent GEO and other event resources remain necessary before making a broad cross-platform
superiority claim.

### Untouched fourth-cohort confirmation

TCGA SpliceSeq LUAD was reserved as an untouched fourth disease cohort. Discovery module
assignments, the 912-event discovery union, network settings, consensus settings and preservation
settings were fixed before LUAD was read. The official LUAD ES matrix contained 12,440 events and
573 samples (514 primary tumors and 59 adjacent normals). Exact
`gene_symbol|event_type|TCGA SpliceSeq AS ID` matching recovered 716 discovery-selected events in
LUAD; all 716 passed the predeclared validation-only missingness and non-zero-variance QC. No LUAD
variance ranking or parameter tuning was performed.

Ten discovery modules were projected from KIRC, KIRP and KICH into LUAD. None met the strict
`Zsummary >= 10` preserved threshold, four were partially preserved (`2 <= Zsummary < 10`) and six
were not preserved. The four partial results were KICH C02 (`Zsummary 6.962`), KIRC C01 (`4.675`),
KIRC C02 (`2.060`) and KIRP C01 (`4.541`). This is evidence for limited cross-cancer program
preservation, not complete module replication.

Locked discovery PC1 loadings were then used to project module scores into LUAD. Five of ten
projected modules were associated with tumor status in both the complete 514-tumor/59-normal
analysis and the stricter 57-patient matched tumor-normal sensitivity analysis. Four additional
modules were significant only in the unpaired complete-cohort analysis and are labelled as
all-sample-only rather than robust. The paired analysis is reported separately to prevent the
unbalanced complete cohort from creating an overconfident phenotype claim.

Independent LUAD module discovery used only the locked discovery-event union and produced four
non-grey modules. KIRC C02 matched LUAD C03 at gene Jaccard `0.447` and independently recovered four
FDR-significant biological themes: RNA processing, nuclear ribonucleoprotein granule, nuclear
speck and ribonucleoprotein granule. LUAD C03 also showed significant RNA-splicing and KEGG
Spliceosome enrichment. This is the clearest independent functional replication in the current
benchmark. It supports reproducible biological interpretation, but does not establish a causal
regulator or universal module identity.

The fourth-cohort result therefore narrows the defensible claim: CoSplice can recover partially
preserved and phenotype-relevant splicing programs with independently replicated RNA-processing
biology, while strict whole-module cross-cancer preservation remains incomplete. LUAD is an
independent disease cohort but still uses the TCGA SpliceSeq platform; a cross-platform cohort is
still required for a broad platform-independent claim.

An additional heterogeneous challenge set used 500 identically preselected events and ten matched
bootstraps in each of GSE62669, GSE141633/TCGA-PRAD and GSE85908. CoSplice produced two valid
modules in all three datasets. Exact WGCNA produced two valid modules in GSE62669, where it had
higher minimum bootstrap ARI (`0.347` versus `0.177`) and stronger within/between separation
(`0.527` versus `0.340`). In PRAD bulk PSI and the single-cell PSI boundary dataset, exact WGCNA
returned all 500 events to grey, while CoSplice retained valid two-module partitions with minimum
bootstrap ARI `0.418` and `0.151`. This is evidence that the methods have different failure modes,
not proof that CoSplice is universally superior.

Biological interpretation in the heterogeneous challenge was deliberately conservative.
GSE62669 CoSplice modules passed the single-gene domination audit and one of two modules had three
significant, broad cellular-component GO terms. No significant GO/KEGG term was found for PRAD,
and GSE85908 has coordinate-only event identifiers, so gene/RBP/GO interpretation was disabled
rather than fabricated. Cross-dataset module preservation was also not calculated because the
three event catalogues are not directly interchangeable.

## Reproducible commands

```powershell
python scripts/prepare_gse62669_exon_usage.py `
  --source GSE62669_RNA110714RH_normalizedExonCounts.txt.gz `
  --output-matrix GSE62669_exon_usage_proxy.tsv `
  --output-metadata GSE62669_sample_metadata.tsv `
  --output-audit GSE62669_gene_audit.tsv `
  --min-feature-sd 1e-8

python scripts/benchmark_public_event_matrix.py `
  --matrix GSE62669_exon_usage_proxy.tsv `
  --metadata GSE62669_sample_metadata.tsv `
  --dataset-id GSE62669_exon_usage_proxy `
  --output benchmark_results/GSE62669_exon_usage_proxy `
  --events 1000 `
  --bootstraps 20

python scripts/run_publication_benchmark.py `
  --manifest publication_benchmark_manifest.json `
  --output publication_head_to_head

python scripts/prepare_tcga_spliceseq.py `
  --source PSI_download_LUAD.txt `
  --output standardized_LUAD `
  --cohort LUAD

python scripts/validate_locked_external_cosplice.py `
  --discovery KIRC KIRC_ES_event_by_sample_PSI.tsv KIRC_module_assignments.tsv `
  --discovery KIRP KIRP_ES_event_by_sample_PSI.tsv KIRP_module_assignments.tsv `
  --discovery KICH KICH_ES_event_by_sample_PSI.tsv KICH_module_assignments.tsv `
  --training-parameters three_cohort_benchmark_parameters.json `
  --validation-name LUAD `
  --validation-matrix standardized_LUAD/LUAD_ES_event_by_sample_PSI.tsv `
  --validation-metadata standardized_LUAD/LUAD_sample_metadata.tsv `
  --output locked_LUAD_validation

python scripts/benchmark_three_independent_public_datasets.py `
  --gse62669-matrix GSE62669_exon_usage_proxy.tsv `
  --gse141633-psi GSE141633_PSI.txt.gz `
  --gse85908-splicing GSE85908_splicing.csv.gz `
  --gse85908-metadata GSE85908_metadata.csv.gz `
  --output public_challenge `
  --events 500 `
  --bootstraps 10
```

## Primary references

- [WGCNA](https://pmc.ncbi.nlm.nih.gov/articles/PMC2631488/)
- [WGCNA/Mantel mammalian brain co-splicing](https://pmc.ncbi.nlm.nih.gov/articles/PMC4429622/)
- [GSE62669](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62669)
- [TCGA SpliceSeq PSI downloads](https://bioinformatics.mdanderson.org/TCGASpliceSeq/PSIdownload.jsp)
- [CoSpliceNet](https://pmc.ncbi.nlm.nih.gov/articles/PMC5086072/)
- [Spycone](https://pmc.ncbi.nlm.nih.gov/articles/PMC9831059/)
- [MAJIQ HET](https://www.nature.com/articles/s41467-023-36585-y)
