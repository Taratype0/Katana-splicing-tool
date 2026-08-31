# Katana-CoSplice

Katana-CoSplice is a standalone Python workflow and an embedded Katana module for discovering **replicate-supported co-splicing programs**. Its matrix source is explicit: either sample-level rMATS events or annotation-independent competing junction usage extracted directly from BAM files.

Implementation boundaries, stage contracts, artifact provenance, and annotation extension points are documented in [COSPLICE_ARCHITECTURE.md](COSPLICE_ARCHITECTURE.md).

It addresses a different question from differential splicing:

- Differential AS asks which events change in one contrast.
- CoSplice asks which AS events have coordinated PSI trajectories across a multi-condition cohort.

## Position relative to published approaches

Katana-CoSplice is not presented as an already proven replacement for every published co-splicing
method. WGCNA provides the established weighted-network/module-eigengene framework, and published
co-splicing studies have applied WGCNA with Mantel correlations or graphical-model selection to
exon/splicing profiles. Katana's current practical advantages are rMATS-native sample PSI/IJC/SJC
ingestion, explicit read/replicate Stable-AS QC, an official R WGCNA signed adjacency/TOM backbone,
condition-stratified bootstrap consensus, observed/residual network separation, factorial metadata association, and
parameter-fingerprinted standalone/UI stages. Katana now provides an internal same-input recovery
benchmark, but general statistical superiority still requires independent-cohort preservation;
feature completeness or one discovery cohort alone is not that proof.

Primary references: [WGCNA](https://pubmed.ncbi.nlm.nih.gov/19114008/),
[WGCNA/Mantel cosplicing](https://pmc.ncbi.nlm.nih.gov/articles/PMC4429622/),
[graphical-model splicing networks](https://pubmed.ncbi.nlm.nih.gov/26449793/), and
[cross-species coregulated splicing networks](https://genome.cshlp.org/content/27/10/1759).
The outer confidence implementation additionally follows the resampling/co-association rationale of
[consensus clustering](https://link.springer.com/article/10.1023/A:1023949509487). Production
modules use official Dynamic Tree Cut on exact TOM. Bootstrap co-assignment annotates unstable event boundaries without changing the WGCNA backbone. Leiden and
the sparse multiscale-network implementation remain available only for historical cache and
comparator compatibility. This combination is a Katana workflow; those references do not by
themselves validate it as a new universally superior method.

The network is therefore fit across all selected samples, not independently inside an `n=3` comparison. A separate residual network can remove explicitly configured covariates before calculating event-event similarity.

## Scientific workflow

1. **Input audit** resolves sample identity, conditions, rMATS bamlists, BAM/BAI coverage, quant coverage and references.
2. **Sample-level matrix** either extracts exact rMATS event PSI/IJC/SJC or builds competing-junction usage directly from BAM CIGAR `N` operations.
3. **Stable-AS QC** filters events using condition coverage, junction-read support, informative PSI range and replicate repeatability.
4. **Observed network** uses beta-posterior-mean PSI with official R WGCNA signed adjacency, exact TOM, Dynamic Tree Cut and eigengene merging across the full cohort; posterior-draw correlations are retained as edge-confidence diagnostics.
5. **Residual network** optionally removes configured covariates and repeats network discovery.
6. **Bootstrap co-assignment confidence** resamples samples within conditions and annotates confidence in the full-data WGCNA assignment without changing that assignment.
7. **Legacy compatibility** can reproduce former sparse/Leiden, tight/diffuse, and boundary-rejecting consensus outputs, but production emits one WGCNA partition with confidence annotations.
8. **Program scores** summarize each module by the first principal component of standardized event PSI.
9. **Association** relates program scores to configured metadata factors and can fit a balanced factorial model with interaction terms and an optional batch covariate.
10. **Module validation** evaluates baseline-module bootstrap recovery, TOM/consensus-parameter sensitivity and observed/residual preservation.
11. **External preservation** projects final consensus modules into explicitly configured independent cohorts and measures density, connectivity, membership and program preservation.
12. **Module evidence annotation** overlays cached DEG, rMATS, DEXSeq, DTU and RBP gene-set support without changing PSI-defined modules.
13. **Method benchmark** compares methods on the same predeclared stable-event universe and records unavailable or degenerate methods explicitly.
14. **Configured factor contrasts** optionally test a predeclared target-reference effect, reference and target perturbation effects, and their difference-in-differences on stable PSI and cached program scores.
15. **Differential network** compares fixed consensus-module density and edge topology across declared sample conditions with permutation tests and global FDR.
16. **AS-primary multi-view fusion** optionally reweights existing AS edges with mapped RBP/DTU/isoform evidence; auxiliary evidence cannot create an edge.
17. **Bayesian preservation** combines density, connectivity, membership and program preservation across external cohorts into conservative module posterior probabilities.
18. **Report** summarizes cached outputs, warnings and provenance.

The read-count uncertainty stage estimates beta-binomial overdispersion, posterior PSI and
measurement reliability from replicate IJC/SJC. It is enabled and parameter-fingerprinted by default: low-reliability events are removed before networking, posterior-mean PSI supplies the WGCNA partition, and repeated posterior draws supply edge and bootstrap confidence. Reliability is retained as event-level evidence and is not multiplied directly into TOM adjacency, because that can collapse otherwise valid modules. PSI-only public cohorts use unmodified `WGCNA::blockwiseModules`; the software never invents read uncertainty. See
[COSPLICE_METHOD_EXTENSIONS.md](COSPLICE_METHOD_EXTENSIONS.md).

The default engine uses posterior-mean PSI and official R WGCNA signed bicor adjacency, exact TOM, Dynamic Tree Cut and eigengene merging. CoSplice's method layer is posterior-draw uncertainty diagnostics, condition-stratified bootstrap confidence, differential topology, cross-cohort validation and event-level interpretation, not a replacement for WGCNA.

## Stable-AS rationale

An rMATS FDR or an observed BAM junction alone does not establish that an event is suitable for a sample-level network. CoSplice also requires that an event:

- is measured in a sufficient fraction of samples;
- has sufficient IJC+SJC support;
- has at least two read-supported replicates in at least 80% of conditions under the default preset;
- has a condition-mean PSI range of at least 0.25 and either an informative intermediate PSI fraction or a sample PSI range of at least 0.50;
- is more consistent between biological replicates than within-condition noise would predict.

The condition-coverage and PSI-range criteria adapt the unsupervised event selection used in the published tissue co-splicing atlas ([Genome Research 2017](https://genome.cshlp.org/content/27/10/1759)). Katana adds rMATS IJC+SJC and replicate-repeatability safeguards because its input is replicate-level rMATS rather than VastDB. Presets (`lenient`, `default`, `strict`) remain user-adjustable and are not universal biological laws.

## Network stability

- Every event passing stable-AS QC enters the default network; there is no fixed Top-N truncation.
- `max_events` is optional and only applies when the user explicitly requests an engineering cap.
- A robust biweight-style correlation is used by default.
- Production materializes the exact dense TOM for every Stable-AS event. `tom_dense_max_events`
  is an explicit memory-safety limit, not a ranking cutoff: larger inputs stop with an actionable
  error instead of silently dropping low-ranked events.
- WGCNA parameters are declared before bootstrap and holdout access; production does not select module boundaries using functional enrichment.
- Stability-first mode uses the complete declared `consensus_bootstrap_iterations`; there is no event-count shortcut. Both configured and effective iteration counts are exported, and downstream validation bootstrap remains a separate module-preservation test.
- Module programs correlated at or above `module_merge_correlation` are merged before final assignment.
- Adaptive membership uses a robust lower bound estimated inside each current module. The lower bound is constrained by `max(adaptive_membership_floor, 1/sqrt(n_samples))` and the fixed `min_module_membership` ceiling. This is a sample-size-aware heuristic, not a significance threshold.
- Adaptive relaxation can retain an event in its current module. Assignment from `grey` or reassignment across modules still requires the original fixed `min_module_membership`, preventing weak bridge events from merging programs.
- Events below the active membership cutoff, or below the optional membership margin, become `grey` instead of being forced into a community.
- `event_module_assignments.tsv` records the original community, final module, membership values, `membership_cutoff_used`, margin and refinement action for every selected event.
- `network_edges.tsv` records raw correlation, signed adjacency and TOM-derived retained edges;
  `network_diagnostics.tsv` records the TOM backbone, cut height, power, memory cap and membership settings.
- Bootstrap resampling is performed **within condition**, preserving group structure while resampling replicates.
- `bootstrap_support` reports how often each retained relationship survives the configured bootstrap procedure.
- `module_resolution_sensitivity.tsv` reports module count, grey events and largest-module fraction over configured TOM cut heights.
- A largest module containing at least 60% of selected events triggers a warning rather than being silently presented as a definitive program.

## Installation

From the repository root:

```powershell
python -m pip install -e .
katana-cosplice --help
```

## CLI

Initialize a configuration by discovering inputs:

```powershell
katana-cosplice init `
  --project-root E:\path\to\project `
  --output E:\path\to\katana_output\09_cosplice_network `
  --config E:\path\to\cosplice_config.yaml
```

Run stages independently:

```powershell
katana-cosplice audit --config cosplice_config.yaml
katana-cosplice build-matrix --config cosplice_config.yaml
katana-cosplice qc --config cosplice_config.yaml
katana-cosplice fit-network --config cosplice_config.yaml --network observed
katana-cosplice consensus-modules --config cosplice_config.yaml
katana-cosplice hybrid-programs --config cosplice_config.yaml
katana-cosplice associate --config cosplice_config.yaml --network observed
katana-cosplice validate-network --config cosplice_config.yaml
katana-cosplice annotate-modules --config cosplice_config.yaml
katana-cosplice benchmark-methods --config cosplice_config.yaml
katana-cosplice preserve-modules --config cosplice_config.yaml
katana-cosplice interpret-modules --config cosplice_config.yaml
katana-cosplice fit-module-models --config cosplice_config.yaml
katana-cosplice meta-analyze-modules --config cosplice_config.yaml
katana-cosplice contrast-factors --config cosplice_config.yaml
katana-cosplice report --config cosplice_config.yaml
```

`katana-cosplice run` is the only explicit chained command. Individual stage commands never run missing upstream stages automatically.

## Required inputs

- Sample metadata containing `sample` or `sample_id`, plus `condition`.
- `matrix_input.input_mode: rmats_events`: rMATS event tables (`SE`, `RI`, `A3SS`, `A5SS`, `MXE`) with replicate-level `IncLevel1/2`, IJC and SJC fields, plus bamlist files or equivalent sample ordering.
- `matrix_input.input_mode: bam_junctions`: coordinate-sorted BAM+BAI for every metadata sample and the current project's GTF. rMATS remains required as an event annotation and differential-evidence layer.

In BAM mode, junctions sharing a splice donor or acceptor are connected into competing clusters. For each sample:

```text
usage = reads supporting one junction / reads supporting all junctions in its competitor cluster
IJC = reads supporting that junction
SJC = competitor-cluster reads excluding that junction
```

This is a LeafCutter-like annotation-independent junction-usage representation, not a reimplementation of LeafCutter differential testing. One BAM junction and one compound rMATS event are different statistical units; their raw event counts must not be compared directly. BAM mode adds observed junctions outside the tested rMATS catalogue, while stable-AS QC controls read support, missingness, dynamic range and replicate consistency before network fitting.

Example configuration:

```yaml
matrix_input:
  input_mode: bam_junctions
  backend: auto
  min_mapping_quality: 1
  min_anchor_bp: 8
  max_intron_length: 500000
  min_total_junction_reads: 20
  min_samples_with_reads: 2
  min_cluster_fraction: 0.001
  exclude_duplicates: true
  max_parallel_samples: 3
```

`backend: auto` prefers local `pysam`, then Ubuntu WSL with `python3-pysam`, then the bundled pure-Python `bamnostic` fallback. Linux/HPC installations should normally use `pysam`; Windows users without WSL can use `bamnostic`, but large BAM scans are slower. Per-sample count caches make unchanged rebuilds fast.

In rMATS mode, `quant.sf`, GTF, genome FASTA and transcript FASTA are optional mechanism-annotation inputs. In BAM mode, GTF is required for project-specific junction annotation; `quant.sf` and the FASTA files remain optional for downstream interpretation.

### BAM-junction versus rMATS functional benchmark

The input-mode benchmark reads completed CoSplice caches and does not rescan BAMs or run a
network. It compares stable-feature QC, sample structure, within-condition replicate
similarity, and leave-one-replicate-out condition effects. For significant SE events it also
compares rMATS dPSI with role-adjusted inclusion/skipping junction usage re-extracted from the
same BAM alignments. This is orthogonal computational support, not an independent biological
validation cohort. Compound rMATS events and individual BAM junctions remain different units,
so their raw event counts are never treated as a direct performance score.

```powershell
katana-cosplice-input-benchmark `
  --bam-root E:\results\katana_output\09_cosplice_bam `
  --rmats-root E:\results\katana_output\09_cosplice_rmats `
  --output-dir E:\results\katana_output\09_cosplice_input_benchmark `
  --setd2-method-summary E:\setd2\06c_method_benchmark\method_benchmark_summary.tsv
```

The optional method summary keeps the input-layer question separate from the network-method
question: the Ythdc1 track tests BAM versus rMATS representations on the same samples, while a
sample-rich cohort such as Setd2 tests CoSplice versus exact WGCNA on the same event universe.

For module evidence annotation, `support_root` may point to cached per-comparison
`gene_level_integrated_candidates.tsv` files and `rbp_gene_set` may point to a custom TSV/CSV. The
stage audits every loaded source. If support files are absent, the module network remains valid and
the missing annotation is reported rather than replaced with zeros or guessed evidence.

Additional metadata columns can be configured for association and residualization. For example, a balanced design can use:

```yaml
association_factors: [condition, genotype, DP_status, treatment, batch]
factorial_factors: [genotype, DP_status, treatment]
batch_column: batch
residual_covariates: [treatment, batch]
```

This fits `genotype * DP_status * treatment + batch` using sum contrasts and Type III ANOVA. Outputs report FDR, partial eta-squared, coefficients and confidence intervals. A batch column must represent a defensible experimental batch; a replicate number used as a provisional batch is recorded as an explicit assumption rather than inferred by CoSplice. Unknown sex is excluded unless a valid sex metadata column is supplied.

## Main outputs

```text
01_audit/input_audit.tsv
02_matrix/sample_level_psi_matrix.tsv
02_matrix/sample_level_ijc_matrix.tsv
02_matrix/sample_level_sjc_matrix.tsv
02_matrix/event_metadata.tsv
02_matrix/sample_metadata_resolved.tsv
02_matrix/matrix_source_selection.tsv
02_matrix/matrix_integrity_summary.tsv
03_stable_as_qc/stable_as_event_qc.tsv
03_stable_as_qc/stable_as_condition_qc.tsv
03_stable_as_qc/stable_as_qc_summary.tsv
03_stable_as_qc/stable_psi_matrix.tsv
04_network/<observed|residual>/event_module_assignments.tsv
04_network/<observed|residual>/network_edges.tsv
04_network/<observed|residual>/module_program_scores.tsv
04_network/<observed|residual>/network_diagnostics.tsv
04_network/<observed|residual>/module_resolution_sensitivity.tsv
04c_consensus_modules/consensus_event_module_assignments.tsv
04c_consensus_modules/consensus_module_summary.tsv
04c_consensus_modules/consensus_module_program_scores.tsv
04c_consensus_modules/coassignment_bootstrap_runs.tsv
04c_consensus_modules/bootstrap_coassignment_edges.tsv
04c_consensus_modules/bootstrap_coassignment_matrix.npz
04c_consensus_modules/consensus_diagnostics.tsv
04c_consensus_modules/consensus_parameter_search.tsv
04d_hybrid_programs/hybrid_event_program_assignments.tsv
04d_hybrid_programs/hybrid_program_summary.tsv
04d_hybrid_programs/hybrid_program_scores.tsv
04d_hybrid_programs/diffuse_bootstrap_coassignment_runs.tsv
04d_hybrid_programs/tight_diffuse_program_overlap.tsv
04d_hybrid_programs/diffuse_engine_audit.tsv
04d_hybrid_programs/external_program_preservation.tsv
04d_hybrid_programs/plots/*.png
05_association/<observed|residual>/module_trait_association.tsv
05_association/<observed|residual>/factorial_module_association.tsv
05_association/<observed|residual>/factorial_module_coefficients.tsv
05_association/<observed|residual>/module_program_score_group_summary.tsv
05_association/<observed|residual>/association_design.tsv
05_association/<observed|residual>/factorial_module_association.png
06_module_validation/event_membership_stability.tsv
06_module_validation/module_bootstrap_stability.tsv
06_module_validation/bootstrap_partition_stability.tsv
06_module_validation/network_parameter_sensitivity.tsv
06_module_validation/observed_residual_preservation_pairs.tsv
06_module_validation/observed_residual_preservation_summary.tsv
06_module_validation/module_validation_report.html
06b_module_annotation/event_module_annotation.tsv
06b_module_annotation/module_gene_annotation.tsv
06b_module_annotation/module_evidence_summary.tsv
06b_module_annotation/module_rbp_enrichment.tsv
06b_module_annotation/support_source_audit.tsv
06b_module_annotation/plots/module_evidence_support_heatmap.png
06c_method_benchmark/common_event_universe.tsv
06c_method_benchmark/event_method_assignments.tsv
06c_method_benchmark/method_benchmark_summary.tsv
06c_method_benchmark/bootstrap_partition_stability.tsv
06c_method_benchmark/module_matching.tsv
06c_method_benchmark/method_availability.tsv
06c_method_benchmark/method_benchmark_report.html
06c_method_benchmark/plots/*.png
06d_external_preservation/external_module_preservation.tsv
06d_external_preservation/external_event_mapping.tsv
06d_external_preservation/external_preservation_null_summary.tsv
06d_external_preservation/external_preservation_target_audit.tsv
06d_external_preservation/plots/*.png
06e_module_interpretation/module_hub_events.tsv
06e_module_interpretation/module_hub_genes.tsv
06e_module_interpretation/event_rbp_motif_hits.tsv
06e_module_interpretation/module_rbp_motif_enrichment.tsv
06e_module_interpretation/module_go_kegg_enrichment.tsv
06e_module_interpretation/dominant_isoform_switches.tsv
06e_module_interpretation/module_isoform_switch_enrichment.tsv
06e_module_interpretation/module_regulator_prioritization.tsv
06e_module_interpretation/interpretation_source_audit.tsv
06e_module_interpretation/motif_scan_parameters.json
06e_module_interpretation/isoform_switch_parameters.json
06f_module_statistics/event_beta_binomial_condition_effect.tsv
06f_module_statistics/module_program_condition_effect.tsv
06f_module_statistics/beta_binomial_vs_simple_delta_dpsi.tsv
06f_module_statistics/beta_binomial_model_diagnostics.tsv
06g_cross_cohort_meta/exact_event_alignment.tsv
06g_cross_cohort_meta/cross_cohort_module_matching.tsv
06g_cross_cohort_meta/cross_cohort_module_state.tsv
06g_cross_cohort_meta/multi_cohort_module_consensus.tsv
06g_cross_cohort_meta/random_effects_module_meta_analysis.tsv
06g_cross_cohort_meta/meta_analysis_cohort_audit.tsv
06h_factor_contrasts/event_factor_contrasts.tsv
06h_factor_contrasts/module_factor_contrasts.tsv
06h_factor_contrasts/phenocopy_occlusion_summary.tsv
06h_factor_contrasts/factor_evidence_summary.tsv
06h_factor_contrasts/presentation_notes.md
06h_factor_contrasts/plots/*.png
07_report/katana_cosplice_report.html
07_report/katana_cosplice_validation_report.pdf
07_report/module_browser.html
state.json
```

Every computational directory includes `parameters.json`, including a parameter/input fingerprint and generated plot list.

## Configured factor contrasts and functional-state hypotheses

Stage `factor_contrasts` is optional and explicit. It never infers target/reference labels from filenames. The user predeclares a primary factor, target/reference levels, context factors, a perturbation factor, baseline and perturbation levels. Event-level effects are mean stable PSI differences with replicate-aware Welch tests. Program-level effects use the selected cached observed or residual CoSplice scores. The interaction is:

```text
(target perturbation - target baseline) - (reference perturbation - reference baseline)
```

Interactions are labelled as attenuated same-direction, enhanced same-direction, opposite response or undetermined. A statistically non-zero interaction is not automatically an occlusion result. Broad signature similarity, enriched significant-event overlap and same-direction agreement are phenocopy evidence only. They must not be reported as proof of loss of function, gain of function or causal pathway regulation without orthogonal molecular and rescue/complementation evidence.

`stable_as_event_qc.tsv` retains one audit row for every input event. To keep large cohorts practical, `stable_as_condition_qc.tsv` contains condition-level replicate/read-support detail only for events that pass the final Stable-AS rules; this scope is recorded in `03_stable_as_qc/parameters.json`.

## Module validation

- Bootstrap iterations resample biological replicates within each condition and refit the complete event partition.
- Event membership stability is the fraction of bootstrap networks in which an event returns to its optimally matched baseline module.
- Module Jaccard measures event-set recovery; whole-network adjusted Rand index measures recovery of the complete partition.
- Default labels are `robust` at validation score >=0.75, `moderate` at >=0.50 and `unstable` below 0.50. The score is the lower of mean matched-module Jaccard and mean event-membership stability.
- Parameter sensitivity compares configured TOM cut heights without changing the active cached network.
- Observed/residual preservation uses event-set Jaccard plus hypergeometric FDR. Low preservation can indicate covariate-dependent structure and is not automatically a software failure.
- These are internal robustness checks. Independent cohorts are still required for external module preservation.

## Bootstrap co-assignment confidence

- Samples are resampled with replacement within each condition, preserving the cohort design and replicate count.
- Each outer bootstrap performs one network partition; inner edge/community bootstraps are disabled to avoid nested resampling.
- Pairwise co-assignment is the fraction of bootstraps in which two events occupy the same non-grey module.
- Production keeps the full-data WGCNA partition exactly. Co-assignment is not cut into a second partition and does not send a WGCNA-assigned event to grey.
- Events with strong within-module co-assignment are labelled `core`; moderately supported events are `peripheral`; weak or ambiguous events are `uncertain` but remain in their WGCNA module. `grey` only means WGCNA itself did not assign the event.
- `partition_source=full_data_wgcna_backbone` and `confidence_action=annotation_only` make this contract machine-readable. `adjusted_rand_index_observed_vs_consensus` must equal 1 in production mode.
- Historical average-linkage, Dynamic Tree Cut, Leiden, and boundary-rejecting modes remain available only through explicit non-default configuration for reproducibility.
- The default uses 20 stratified sample bootstraps. Higher values such as 50 or 100 are optional sensitivity analyses, not a requirement for routine discovery.
- Event confidence is mean co-assignment to the selected module; the confidence margin compares that value with the next-best module. `assignment_uncertainty` is one minus confidence.
- `consensus_parameter_search.tsv` records the fixed production partition and its audit metrics. No Setd2 labels, conditions or expected biological answer enter confidence estimation.
- The complete matrix is stored as compressed NPZ (`event_keys`, `coassignment`); the TSV edge view is intentionally capped for preview.
- This is internal consensus in the discovery cohort, not independent validation.

## External module preservation

Each `external_cohorts` entry supplies an independent sample-level PSI matrix and sample metadata.
Events map by exact `event_key`; low overlap is reported rather than imputed. For each consensus module,
Katana reports target density, reference-to-target connectivity correlation, event membership-profile
correlation, PC1 variance similarity and a program-score preservation score. Random event sets of the
same size calibrate `Zdensity` and `Zconnectivity`; their mean is exported as
`Zsummary_wgcna_like`. The labels preserved (>=10), partially preserved (>=2), and not preserved (<2)
follow familiar WGCNA interpretation bands, but this statistic is explicitly not the exact
`WGCNA::modulePreservation` implementation.

## Head-to-head benchmark

The independent `benchmark_methods` stage compares Katana's **TOM-first consensus partition** with a fixed-cut
Pearson hierarchical baseline, a Ledoit-Wolf shrinkage graphical-lasso partial-correlation baseline,
and the exact R `WGCNA::blockwiseModules` implementation when Rscript and WGCNA are installed. If R
or WGCNA is absent, the method is marked unavailable; no Python approximation is presented as WGCNA.
When a completed `04c_consensus_modules` cache already exists, the benchmark includes the final
CoSplice bootstrap-consensus modules only if its event universe exactly matches the benchmark
universe. It never runs consensus automatically and never subsets full-universe modules into a
top-N benchmark, because module metrics from different discovery universes are not directly
comparable. A missing or mismatched cache is reported explicitly.

Katana searches `KATANA_RSCRIPT`, `PATH`, the active Python/Conda prefix and sibling Conda
environments for an R installation containing WGCNA. The benchmark also writes
`method_metric_ranking.tsv`, `method_strengths_limitations.tsv` and
`plots/head_to_head_scorecard.png`. The equal-weight composite in the interpretation table is a
descriptive navigation score; scientific interpretation should retain the raw stability,
separation, modularity, grey-event and runtime metrics.

The default common universe contains the 150 highest cached Stable-AS `stability_score` events. This
bounded universe is required for a fair comparison with cubic-complexity partial-correlation
estimation and does not affect the full production network. The report includes condition-stratified
bootstrap ARI, cosine silhouette, neutral KNN-graph modularity, within/between-module Pearson
separation, descriptive condition-pattern R2, module matching, runtime, and method diagnostics.

Zero-module and one-module partitions are classified as `degenerate`. Their bootstrap ARI is not
interpreted because an unchanged all-grey partition would otherwise appear spuriously perfect. These
metrics are internal; external module preservation is still required for a general superiority claim.

### Consensus-network benchmark boundary

The current frozen public panel contains 11 same-input runs across seven AS event types, three ES
event-universe sizes and three holdout cohorts. The unified-seed paired rerun
gave mean bootstrap ARI `0.867` for TOM-first CoSplice and `0.871` for exact
WGCNA, with 3 CoSplice wins, 7 ties and 1 loss. Mean holdout program-score
preservation was `0.497` versus `0.472`; both methods recovered non-empty
partitions in 11/11 runs and strict valid multi-module partitions in 2/11.
No tuning-GMT functional objective passed the predeclared FDR gate. These
results support WGCNA-level internal stability and a modest average preservation
gain; they do not establish universal or functional superiority.
See [COSPLICE_LOCKED_PANEL_VALIDATION.md](COSPLICE_LOCKED_PANEL_VALIDATION.md).

### Public reciprocal external validation

The reusable `scripts/benchmark_public_spliceseq.py` entry accepts any two event-by-sample PSI
matrices with matching event keys. Event filtering and variance ranking are performed on the
discovery cohort only. It compares bootstrap-consensus CoSplice with the exact R WGCNA package,
then projects both partitions into the untouched external cohort. The command writes matched
bootstrap tables, module assignments, permutation-calibrated preservation results, plots,
`parameters.json`, and `PUBLIC_VALIDATION_REPORT.md`.

The latest primary TCGA SpliceSeq regression uses 500 discovery-selected ES events, KIRC discovery,
KIRP tuning and untouched KICH holdout. In the paired rerun, CoSplice and exact
WGCNA were identical: one broad module, bootstrap ARI `0.883`, KICH
Zsummary-like preservation `0.281`, and program-score preservation `0.582`.
The public files are obtained from the official TCGA SpliceSeq download and the benchmark
remains a validation workflow rather than a production-network stage.

### Publication head-to-head aggregation

`scripts/run_publication_benchmark.py` combines completed public benchmark tracks without changing
their underlying module assignments. A JSON manifest identifies each accession, design, input
contract, validation role, citation and result directory. It writes dataset/method metrics, paired
CoSplice-versus-WGCNA contrasts, consensus coverage and overlap, external preservation, comparator
scope, a claim audit, a PNG scorecard, Markdown/HTML
reports and `parameters.json` with SHA-256 source fingerprints.

Only the production CoSplice partition, exact R WGCNA and fixed same-PSI baselines are pooled into
single-partition numerical comparisons. Legacy tight/diffuse outputs are historical views rather
than a second production result. WGCNA/Mantel co-splicing, CoSpliceNet, Spycone and MAJIQ
HET remain domain-specific published references because they use different biological units or
experimental designs.

```powershell
python scripts/run_publication_benchmark.py `
  --manifest configs/publication_benchmark.example.json `
  --output E:\benchmark_results\publication_head_to_head
```

The aggregation layer rejects universal-superiority wording. A zero- or one-module partition is
degenerate and cannot win a stability comparison merely by assigning every event to the same label.
Direct claims require the identical discovery-selected event universe, valid partitions from both
methods, matched resampling and, where available, untouched external-cohort preservation.

## Katana integration

Open `9. CoSplice Network` and initialize or choose a config. Page refresh reads only `state.json` and cached PNG files. Every button starts exactly one stage in a worker thread. If upstream output is missing, the current stage reports what must be run first; it does not run it automatically.

## Interpretation limits

- A co-splicing module is not automatically an RBP regulon.
- DEG, DTU and RBP enrichment are downstream support annotations; they do not redefine event membership.
- Association p-values are calculated on the same cohort used for discovery and are prioritization statistics, not independent validation.
- DEG, DTU, RBP and motif support can annotate plausible mechanisms but cannot establish causality from bulk RNA-seq alone.
- Public cohorts or an independent experiment should be used for module preservation and external validation.
