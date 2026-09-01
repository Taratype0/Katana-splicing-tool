# CoSplice Method Extensions

This document defines four CoSplice method extensions, their scientific scope,
and the ablation evidence used to choose production defaults.

## 1. Read-count uncertainty

`estimate_beta_binomial_uncertainty` uses replicate-level rMATS IJC/SJC counts
to estimate condition-adjusted beta-binomial overdispersion, effective reads,
posterior PSI, posterior standard deviation and event reliability.

It is a measurement-uncertainty model, not a differential-splicing test.
Available network modes are:

- `diagnostic_only`
- `raw_reliability_weighted`
- `posterior`
- `posterior_reliability_weighted`
- `posterior_ensemble`

The production mode is `posterior`: reliability first filters low-confidence
events and posterior-mean PSI supplies the official WGCNA partition. Repeated
beta-posterior draws separately estimate edge correlation SD, sign probability
and bootstrap co-assignment confidence. These are uncertainty diagnostics, not
p-values, and do not replace the production partition. The two
reliability-weighted modes are kept
for research reproduction only. They multiply reliability into network edges
and are not production defaults because a real-count ablation showed that this
operation can destroy otherwise coherent TOM structure and return an all-grey
solution. `diagnostic_only` is the explicit raw-PSI sensitivity arm.

The v1.2 four-arm ablation compares raw-PSI WGCNA, posterior-mean WGCNA,
posterior-ensemble WGCNA, and posterior-mean WGCNA plus annotation-only
posterior-draw bootstrap confidence. All arms use identical events, seeds and stratified
resamples. Posterior-ensemble partitioning was rejected as the default because
its high ARI coincided with an over-dominant module. Full CoSplice therefore
keeps posterior-mean WGCNA and uses posterior draws only for event, module, kME,
hub and edge-confidence reporting.

From v1.2.1, confidence annotation combines label-invariant co-assignment with
fixed-module bootstrap kME. Final `core` and `peripheral` labels must reproduce
in two disjoint bootstrap halves. This corrects the large-module penalty of
pairwise co-assignment without changing the WGCNA partition or lowering the
predeclared evidence thresholds.

Scientific basis: beta-binomial read-count uncertainty is established in
replicate-aware splicing models such as FRASER and Bisbee. CoSplice uses it for
measurement weighting rather than claiming a new differential-AS p-value.

## 2. Differential co-splicing network

This stage keeps the global consensus modules fixed, builds separate
shrinkage-correlation networks within each metadata condition, and tests:

- module density gain or loss;
- topology rewiring;
- event connectivity change.

The module identity is fixed because independent community detection from
three replicates per condition is unstable. Label permutations provide density
and topology p-values. `density_FDR` and `topology_FDR` are corrected globally
across all condition-pair/module tests; pair-local FDR columns are exploratory.

Known-truth simulation recovered the gained module and produced no false
positive for a strict stable module. In Setd2, 264 tests across 66 condition
pairs and four modules yielded no globally significant density/topology change.
This is a valid negative result, not evidence that mean dPSI is unchanged.

Scientific basis: the question follows differential-network methods such as
DiffCoEx and multiWGCNA, while preserving CoSplice's event-level AS identity.

## 3. AS-primary multi-view fusion

The AS network remains primary. RBP motif, RBP expression, DTU, DEG or isoform
features may reweight an existing AS edge, but cannot create an edge unsupported
by sample-level splice usage. The AS contribution is constrained to at least
50%; the tested default uses 80%.

Known-truth edge AUROC changed from 0.991 for AS only to 0.998 with an
informative 20% auxiliary view. An adversarial random view produced 0.987,
showing limited degradation under the AS-primary guardrail.

The Setd2 audit evaluated 100,000 AS edges. RBP motif evidence covered 97.7% of
edges; mapped isoform/DTU evidence covered 14.6%. Median absolute edge
reweighting was 0.012, so auxiliary layers did not dominate the AS signal.

This implementation is an AS-primary multiplex edge fusion inspired by
Similarity Network Fusion. It is not the canonical iterative SNF algorithm and
must not be described as causal regulator inference.

## 4. Bayesian cross-cohort preservation

The hierarchical preservation model combines:

- permutation-calibrated density evidence;
- direct connectivity correlation;
- membership-profile correlation;
- program-score preservation.

Correlation evidence is Fisher transformed and pooled across cohorts with
between-cohort heterogeneity. Outputs are posterior probabilities and credible
intervals. A negative historical `Zconnectivity` is not misinterpreted as a
negative observed connectivity correlation.

On reciprocal TCGA-KIRC/KIRP/KICH projections:

| Method | Discovery modules | Fixed partial-or-better | Bayesian partial-or-better |
|---|---:|---:|---:|
| CoSplice | 10 | 30% | 30% |
| Exact WGCNA | 3 | 0% | 0% |

CoSplice resolved 3–4 modules per cohort, whereas exact WGCNA returned one
broad module per cohort. Exact WGCNA had higher bootstrap ARI, but that number
must be interpreted with its degenerate single-module partition. CoSplice's
advantage in these cohorts is program resolution and external projection, not
runtime or universal stability superiority.

## 5. Optional gene-level Mantel scaffold

CoSplice can use a gene-level exon-profile Mantel network as a bounded prior for
event-level module discovery. This is a hierarchical extension, not a replacement
for the AS network:

1. map each network event to its gene;
2. construct one multi-exon profile per eligible gene from a sample-level exon
   count matrix;
3. calculate gene-gene Mantel similarity from Canberra sample-distance vectors;
4. derive a gene-program scaffold;
5. reinforce only AS edges that already exist and whose genes share scaffold
   support.

The reinforcement is bounded:

`reinforced = AS + weight * AS * (1 - AS) * gene_prior`

Therefore a gene prior cannot create an edge where the AS network has no support,
and its configured weight is restricted to at most 0.50. The same cached prior is
used inside every sample bootstrap so that the observed partition and confidence audit use
the same method. Event identity, PSI values, and differential-splicing statistics
remain event-level.

Example:

```yaml
gene_prior:
  enabled: true
  exon_count_matrix: /path/to/sample_level_exon_counts.tsv
  feature_gene_map: null
  encoded_gene_rule: prefix_before_last_underscore
  min_exons: 2
  soft_power: 6
  prior_mode: scaffold_reinforcement
  scaffold_knn_neighbors: 25
  scaffold_min_module_size: 10
  weight: 0.10
```

The exon-count table uses exon or feature IDs as rows and sample IDs as columns.
When feature IDs do not encode a gene ID, provide `feature_gene_map` with the
configured `feature_id_column` and `gene_id_column`.

Generated audit artifacts are:

- `gene_mantel_prior_matrix.npz`;
- `gene_mantel_prior_gene_audit.tsv`;
- `gene_mantel_prior_event_audit.tsv`;
- `base_AS_adjacency`, `gene_mantel_prior`, and
  `gene_prior_adjacency_boost` in sparse network-edge output.

The feature is disabled by default. In a locked 500-event GSE62669 ablation with
20 sample bootstraps, scaffold weight 0.10 increased median bootstrap ARI from
0.333 to 0.347 and mean consensus confidence from 0.465 to 0.470 while reducing
grey fraction from 0.502 to 0.498. Weight 0.30 reduced ARI to 0.328. This is a
small cohort-specific improvement, not evidence of universal superiority.
Production use should compare weight 0 against predeclared candidate weights on
internal stability and structural metrics without using GO terms or phenotype
labels to tune the network.

## Reproducible benchmark

```powershell
python scripts\benchmark_cosplice_method_extensions.py `
  --cosplice-output <completed-cosplice-output> `
  --public-preservation <multi-cohort-preservation.tsv> `
  --output <benchmark-output-directory>
```

The gene-prior ablation is run separately because it requires a sample-level exon
count matrix:

```powershell
python scripts\benchmark_gene_prior_ablation.py `
  --event-matrix <locked-event-by-sample-matrix.tsv> `
  --sample-metadata <sample-metadata.tsv> `
  --exon-counts <exon-feature-by-sample-counts.tsv> `
  --weights 0,0.1,0.2,0.3 `
  --output <benchmark-output-directory>
```

The benchmark writes known-truth ablation, real IJC/SJC stability, real
differential-network and multi-view audits, public Bayesian preservation,
recommendations and a PNG summary. Benchmark outputs belong outside the source
tree and are never required at runtime.

## References

- Bisbee: <https://www.nature.com/articles/s41598-021-89938-2>
- FRASER: <https://www.nature.com/articles/s41467-020-20573-7>
- Similarity Network Fusion: <https://www.nature.com/articles/nmeth.2810>
- DiffCoEx: <https://pubmed.ncbi.nlm.nih.gov/20925918/>
- multiWGCNA: <https://pmc.ncbi.nlm.nih.gov/articles/PMC10039544/>
- WGCNA module preservation: <https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1001057>
- Mantel co-splicing network: <https://pmc.ncbi.nlm.nih.gov/articles/PMC4429622/>
