# Katana Next Scientific Methods

## 1. Principles

Katana Next separates four questions that the legacy product sometimes mixes:

1. Is a splice event measurable and reproducible?
2. Is its condition effect statistically supported?
3. Does it participate in a reproducible co-splicing program?
4. Is there independent evidence that helps interpret that event or program?

Read depth, stability, significance, and biological interpretation are reported
separately. They are not collapsed into one unexplained score. Unassigned or
insufficiently supported events remain unassigned rather than being forced into
a module or significance class.

## 2. Canonical bulk splicing layer

The primary data object is a sample-level competing-junction count matrix with
annotation-aware event mappings. rMATS, LeafCutter, MAJIQ, SUPPA, BAM-derived
junctions, and future callers are adapters into this layer; no single caller is
treated as the complete biological event universe.

For every event, Katana stores:

- genomic build, strand, splice sites, event type, gene and transcript links;
- inclusion/exclusion or competing-junction counts by sample;
- PSI or intron-excision estimates and uncertainty;
- annotation status and caller provenance;
- exact-event identity plus explicitly separate same-locus and same-gene links.

Exact overlap requires the same normalized event structure after coordinate
conversion. Opposite direction is still the same event but is a discordant
effect. Same-gene or nearby-locus overlap is never silently called the same
event.

## 3. Stable-event QC

Stability is estimated before differential testing and program discovery using
predeclared, configurable gates:

- minimum informative read support and sample coverage;
- minimum number or fraction of supported biological replicates per group;
- replicate PSI dispersion and leave-one-replicate-out effect consistency;
- posterior PSI uncertainty from inclusion/exclusion counts;
- mapping/annotation confidence and outlier diagnostics.

The UI offers validated presets plus a custom profile, but always shows the
actual numeric criteria and attrition at each gate. Stability is not inferred
from a small FDR alone.

## 4. Differential splicing

The default event-level test uses replicate-aware count likelihoods appropriate
to the event representation:

- beta-binomial regression for binary inclusion/exclusion evidence;
- Dirichlet-multinomial regression for multi-junction intron clusters;
- imported rMATS statistics as a caller-specific result, not a substitute for
  sample-level diagnostics.

Models report effect size, confidence/credible interval, nominal p-value, BH
FDR, read support, replicate consistency, and model diagnostics. Optional batch
is included only when identifiable. Donor or paired designs are modeled when
metadata support them. Arbitrary fixed direction inferred from filenames is not
part of the new model.

Cross-comparison analysis reports separately:

- exact-event overlap and containment;
- direction concordance among shared exact events;
- effect-size correlation with uncertainty;
- shared, condition-specific, gained, lost, and reversed effects;
- overlap enrichment against a declared tested-event universe.

## 5. CoSplice 2.0

### 5.1 Production backbone

Katana Next does not use the legacy Leiden-first or competing multi-engine path
as its production result. The workflow locks module discovery to official
signed blockwise WGCNA: data-driven soft-power diagnostics, signed adjacency,
TOM, hierarchical clustering, Dynamic Tree Cut, eigengene-based merging, and
kME membership. Dense Python TOM and exploratory sparse-network results remain
benchmark code and are not presented as a second production answer.

### 5.2 Locked baseline

Official signed blockwise WGCNA is the mandatory comparator, not something to
discard because it is old. It provides soft adjacency, TOM, hierarchical
clustering, Dynamic Tree Cut, module eigengenes and membership. Soft power is
selected from predeclared scale-free topology and connectivity diagnostics, not
hard-coded to six without evidence.

Reference: Langfelder and Horvath, 2008,
https://doi.org/10.1186/1471-2105-9-559; Dynamic Tree Cut,
https://doi.org/10.1093/bioinformatics/btm563.

### 5.3 Production stabilization

CoSplice adds an uncertainty-aware stabilization layer to the locked WGCNA
backbone:

1. Learn posterior splice usage from replicate junction counts.
2. Generate paired sample bootstraps and posterior imputations.
3. Fit the same signed blockwise network protocol in every resample.
4. Estimate event-event co-assignment, module recovery, kME uncertainty, and
   event coverage.
5. Retain stable WGCNA assignments and send weak or ambiguous boundary events
   to a grey/unassigned set. Stable modules are not repartitioned.
6. Refit module scores on the discovery data and project them unchanged into
   validation cohorts.

This is `backbone_refined_consensus`, not a de novo consensus clustering
algorithm. Its value must be demonstrated against exact WGCNA under locked,
same-input benchmarks; rejecting unstable events is not automatically an
improvement if it destroys coverage or external preservation.

### 5.4 Module statistics

Module discovery is unsupervised. Phenotype labels, GO terms, RBP sets, and
holdout data do not define or tune membership.

For each module Katana reports:

- bootstrap recovery and co-assignment confidence;
- split-half stability and module separation;
- event kME with uncertainty and grey-event rate;
- module-score association with genotype, treatment, interaction, and optional
  batch using the biological sample as the unit;
- effect size, confidence interval, permutation-calibrated p-value where
  needed, and BH FDR across all tested module-trait pairs;
- external density, connectivity, membership and score preservation;
- event and annotation coverage in every validation cohort.

Condition-specific gained/lost/rewired networks require enough independent
samples per condition. Three replicates per condition do not support formal
condition-specific correlation-network inference; in that setting Katana may
test module-score effects but must disable or label rewiring as exploratory.

### 5.5 Biological interpretation

Interpretation is applied after module discovery:

- GO/KEGG or user gene-set enrichment with an explicit tested background;
- RBP motif/eCLIP enrichment and RBP-expression/module-score association;
- DEG, DTU, isoform-switch and transcript consequence support;
- hub events, hub genes, uncertainty and Native Sashimi read evidence;
- cross-cohort functional replication.

These prioritize mechanisms; they do not establish an upstream causal RBP.
Functional coherence is an external validation endpoint and must not be used to
optimize the same modules later claimed to be functionally coherent.

### 5.6 Benchmark decision rule

Compare exact WGCNA and candidate CoSplice on identical features, samples,
seeds, and bootstrap draws. Freeze all parameters before untouched holdout
evaluation. Primary endpoints are:

- module recovery and adjusted Rand/co-assignment stability;
- external module preservation and event coverage;
- calibrated module-trait false-positive rate and power in simulation;
- independent functional replication;
- resistance to one giant module, excessive singleton modules and forced
  assignment;
- runtime and peak memory as secondary endpoints.

No weighted composite score hides tradeoffs. Results are shown as a Pareto
profile. WGCNA remains production default unless the extension demonstrates a
consistent advantage without worse calibration.

## 6. Single-cell splicing

### 6.1 Assay-aware entry gate

Single-cell splicing is a separate workflow, not another bulk page. Input
requires aligned reads with cell barcode/UMI tags, cell metadata, sample/donor
IDs, GTF, genome build, and assay chemistry. Katana distinguishes:

- droplet 3-prime/5-prime protocols;
- full-length plate protocols such as Smart-seq;
- polyA endpoint/DTU use cases;
- optional long-read confirmation.

If donor/sample IDs are absent, Katana can explore cell-state patterns but must
not report population-level differential significance by treating cells as
replicates.

### 6.2 Junction calling and confidence

UMI-aware junction evidence is aggregated without double counting. Junctions
receive annotation and statistical confidence before PSI-like quantities are
computed. This follows the problem addressed by SICILIAN, which improves
single-cell splice-junction precision and replicate concordance:
https://doi.org/10.1186/s13059-021-02434-8.

### 6.3 Quantification tracks

Katana chooses a compatible track rather than promising one universal PSI:

- **Droplet track:** annotation-free intron usage and a SpliZ-like gene-level
  splicing score for sparse data. SpliZ was introduced in Nature Methods as an
  annotation-free statistic for regulated splicing in droplet scRNA-seq:
  https://doi.org/10.1038/s41592-022-01400-x.
- **Event track:** event-level posterior PSI for full-length or sufficiently
  covered data, with BRIE2-style uncertainty-aware regression:
  https://doi.org/10.1186/s13059-021-02461-5.
- **Multi-junction track:** annotation-free intron-cluster usage and
  Dirichlet-multinomial testing, aligned with the robust scQuint task model:
  https://doi.org/10.7554/eLife.73520.
- **PolyA/DTU track:** transcript-end usage is treated separately, following
  the distinct problem addressed by Sierra:
  https://doi.org/10.1186/s13059-020-02071-7.

Newer imputation or flow-conservation methods such as SCSES and ELLIPSIS are
benchmark candidates, not automatic defaults. Imputed values cannot create
independent biological replication.

### 6.4 Statistical unit and tests

The primary confirmatory unit is the biological sample/donor. For each cell
type, Katana aggregates junction counts by donor and fits beta-binomial or
Dirichlet-multinomial models with condition and optional batch. It reports cell
count, donor count, coverage, effect size, interval, FDR, and leave-one-donor-out
consistency.

Cell-level models are reserved for supported exploratory questions and include
sample structure. P-values from thousands of cells are never presented as if
they came from thousands of independent animals or patients.

### 6.5 Cell-state and co-splicing programs

Katana supports two distinct questions:

- SpliZ/Psix-style association of splicing with a cell-state manifold or
  trajectory. Psix demonstrated probabilistic state-associated splicing and
  co-regulated exon modules without requiring predefined clusters:
  https://doi.org/10.1101/gr.276109.121.
- CoSplice programs built from donor-by-cell-type pseudobulk profiles when the
  number of independent units passes a feasibility gate.

Cells are not used as pseudo-samples for WGCNA. Program discovery is disabled
when independent-unit count is inadequate; descriptive event patterns remain
available.

### 6.6 Single-cell outputs

- assay and junction-confidence QC;
- cell-type by donor coverage matrix;
- gene-level splicing landscape and cell-state score plots;
- cell-type differential AS volcano/forest plots;
- effect-size and uncertainty heatmaps;
- sample-aware similarity/PCA;
- trajectory-associated event modules when supported;
- RBP-expression/motif support and DEG/DTU integration;
- pseudobulk Native Sashimi by cell type and donor group;
- AnnData-compatible annotations, Parquet/TSV tables, SVG/PNG figures, and a
  provenance-complete report.

MARVEL demonstrates the value of an integrated single-cell AS interface across
plate and droplet data (https://doi.org/10.1093/nar/gkac1260), while a 2026
cross-method assessment found scQuint robust across tested full-length datasets
and tasks (https://doi.org/10.1016/j.isci.2026.116090). Katana must therefore
benchmark by assay and task rather than claim one method is globally best.

## 7. Claim boundaries

Katana may claim reproducible association, overlap, preservation, functional
enrichment, or regulator prioritization when the corresponding test passes. It
must not claim causal RBP control, protein-domain disruption, condition-specific
network rewiring, or population-level single-cell significance without the
required evidence and experimental units.
