# Katana Nextflow: parameters and analysis choices

This reference describes the **current Katana implementation**, not every
option of `nf-core/rnaseq`, `rnasplice` or `differentialabundance`. For formats
and commands, use [NEXTFLOW_USAGE.md](NEXTFLOW_USAGE.md). Defaults below were
checked against `pipeline/nextflow.config`, the public validator, importer and
downstream planner for the accepted 25 October 2026 release.

The four strandedness-inference parameters below are available in the
**current `main` source workflow**. Use that checkout's
`pipeline/main.nf`, as described in [the usage guide](NEXTFLOW_USAGE.md#2-installation-and-directory-layout).
The frozen `v2.0.0-beta.5.dev20261025` wheel/export still requires explicit
strandedness and remains unchanged. The released downstream container remains
usable because the inference processes use their own pinned RSeQC/gffread images.

## Configuration layers

**Current-source interface:** [direct downstream parameters](NEXTFLOW_SIMPLE.md) now
generate internal JSON. See [all native defaults](NEXTFLOW_NATIVE_OPTIONS.md).
The advanced JSON fields below remain supported; older releases retain their old interface.

| Layer | Where it belongs | Examples |
| --- | --- | --- |
| Nextflow launcher | Single-dash Nextflow arguments | `-profile`, `-params-file`, `-c`, `-work-dir`, `-resume`. |
| Pipeline parameters | `--parameter` or a Nextflow YAML/JSON params file | `--mode`, `--input`, `--contrasts`, reference paths and rMATS settings. |
| Study/resources | Direct `--species` and resource paths; generated template, or an existing imported `project.json` | Species, exact samples/contrasts, resources, JC/JCEC and original result provenance. |
| Downstream plan | Direct switches/thresholds; generated `flow.json`, or an advanced explicit config | Selected contrasts, RS-DAS, per-module analysis parameters, independent networks, QC and MAPP. |
| Site execution | Nextflow `-c site.config` | Fully qualified process resource selectors, account policy and executor configuration. |
| Managed storage runner | `katana-flow run` arguments | `--outdir`, `--run-dir`, `--keep-work`, `--resume`, `--profile`, repeated `--site-config`. |
| Deep dive | Exported request plus `bindings.json` and target CLI arguments | Chosen source row IDs, target scope, BAM/index locations and worker budget. |

The current source workflow supports `--deg_fdr`, `--deg_log2fc`, `--das_fdr`, `--das_dpsi`,
`--dtu_fdr`, `--rmats_mode` and `--event_set`. Use `--event_set stable_significant`
for RS-DAS. There is no `--rs_das` or `--motif_top10` flag; the approved significant-only
Top 10 presentation is unchanged.

Example `inputs/params.yaml` for `nextflow ... -params-file inputs/params.yaml`:

```yaml
mode: complete
input: inputs/samples.csv
contrasts: inputs/contrasts.csv
fasta: references/genome.fa
gtf: references/annotation.gtf
reference_build: YOUR_BUILD
annotation_release: YOUR_RELEASE
species: Mus musculus
gene_sets: resources/pathways.gmt
rbp_gene_sets: resources/rbps.gmt
motif_database: resources/motifs.json
clip_peaks: resources/matched-peaks.tsv
katana_container: YOUR_MATCHING_KATANA_IMAGE_OR_SIF
dexseq_tools_container: YOUR_MATCHING_DEXSEQ_IMAGE_OR_SIF
outdir: results
```

Relative pipeline paths above resolve in the launch context. Raw-sheet read
paths resolve from the sheet; template resources from the template; flow
resource paths from the flow configuration; target bindings from the bindings
file. Preserve those bases when transporting a study. CLI `--parameter`
overrides values from the params file. Keep scientific params in the params file
or CLI and executor/resource configuration in `-c`; do not rely on a borrowed
nf-core profile to add undeclared methods to Katana.

## 1. Input/output and mode

| Parameter | Default | Requirement / effect |
| --- | --- | --- |
| `--mode` | `complete` | `complete`, `upstream-only` or `downstream-only`. No `report-only` Nextflow mode; use the report CLI. |
| `--input` | unset | Raw FASTQ/BAM samplesheet; required for complete/upstream-only. |
| `--contrasts` | unset | Selected contrast CSV/TSV; required for complete/upstream-only. |
| `--fasta` | unset | Matching genomic FASTA, required for raw-entry modes. |
| `--gtf` | unset | Matching GTF, required for raw-entry modes. |
| `--project_template` | unset | Advanced resource/metadata template; generated automatically by the direct route. |
| `--project` | unset | Katana v2 project manifest; required for downstream-only. |
| `--flow_config` | unset | Advanced downstream configuration; generated automatically when omitted. |
| `--flow_plan` | unset | Original sealed executed downstream plan; only downstream-only. Bypasses replanning, so an additional `--flow_config` is not applied. |
| `--outdir` | `results` | Published output root; keep outside input project/config roots. |
| `--publish_dir_mode` | `copy` | Upstream general publication mode. Some processes explicitly retain `copy` for portable deliverables; this is not a global storage-deduplication switch. |
| `--storage_mode` | `retain` | `retain` preserves legacy publication. `lean` avoids known large duplicate publications, retains full native statistics/models and forces portable copies. This flag alone does not delete work; the managed runner selects lean and verifies stage cleanup. |
| `--managed_preflight` | `false` | Internal managed upstream-only gate: require complete-flow resource/template preflight before upstream execution. Requires the matching project template and flow configuration; does not change scientific tests. |
| `--flow_report_revision` | `1` | Explicit presentation-only report task revision. Letters/digits/dot/underscore/hyphen; no scientific threshold change. |
| `--help` | `false` | Print concise pipeline help. This guide is the full parameter reference; an nf-core schema-generated launcher is not claimed. |

The required raw-sheet headers are `sample,condition,strandedness`, plus actual
FASTQ/BAM fields. The validator supplements these with explicit resolved
`source,assay,single_end,read_length,subject,batch` fields. See the complete
[samplesheet and contrast specification](NEXTFLOW_USAGE.md#3-samplesheet-input).

### Managed execution and cleanup commands

The optional managed wrapper below retains its explicit configuration interface.
The simplified `nextflow run` route above does not require these two JSON files.

These are **`katana-flow` CLI arguments**, not additional scientific Nextflow
parameters. They require an installation containing the updated managed
runner; old wheels/images/exported workflows are not changed automatically.
See [storage management](STORAGE_MANAGEMENT.md) for complete commands and
retention guarantees.

| Command/argument | Default | Effect |
| --- | --- | --- |
| `run --mode` | `complete` | `complete` runs native upstream then downstream; `downstream-only` starts from a prepared project. |
| Run `--outdir` | required | New final delivery root: retained project, fixed report, native supporting outputs and evidence. |
| Run `--run-dir` | sibling `OUTDIR.run` | Private ownership/state/work directory; not a shared arbitrary work root. |
| Run `--keep-work` | off | Disable automatic successful-stage work deletion; keep lean publication. |
| Run `--resume` | off | Resume this exact managed run with unchanged input/software/configuration identities. Completed verified handoffs are reused; unknown engine termination blocks blind resubmission. |
| Run `--profile` | `local` | Native Nextflow profile string, e.g. `slurm,apptainer`. |
| Repeated run `--site-config` | none | Explicit native Nextflow site configuration files; no implicit site substitution. |
| Run `--params-file` | none | Supported native scientific/resource params file. Managed routing, input/output locations and cleanup arguments are controlled by the runner. |
| Run `--input --contrasts --fasta --gtf --project_template --flow_config` | required for complete | Same bulk samplesheet, requested contrasts, references and downstream resource settings as the raw complete route. |
| Run `--project` | required for downstream-only | Prepared Katana project; source resources remain available until final retention succeeds. |
| Run `--flow_plan` | none | Exact compatible original sealed plan, downstream-only; replaces replanning with `--flow_config`. |
| Run `--katana_container --dexseq_tools_container --flow_executable` | unset; unset; `katana-flow` | Same declared runtime selection as native Nextflow; container profiles require matching usable images. |
| Run `--pipeline` | installed bundled workflow | Optional explicit `main.nf`; identity participates in resume checks. |
| Run `--nextflow-command-json` | `["nextflow"]` | Explicit executable-prefix argv, useful for a pinned native Nextflow/Java launcher; no shell evaluation. |
| `clean-run --run-dir` | required | Dry-run cleanup plan for an already completed, verified managed scope; refuses arbitrary/failed/running work. |
| Clean `--apply` | off | Apply that bounded cleanup plan after final delivery revalidation. |
| `retain-project --collection --project --outdir` | required | Retain original completed products and required resources with independent bindings; no scientific recalculation. |
| `verify-retained --directory` | required | Verify and reload a retained project at its current location. |

Automatic cleanup removes the successful upstream work after complete bundle
handoff, then successful downstream work after report/project/delivery
verification. Original inputs and shared resources are preserved. Container and
Conda caches stay outside reclaimable work and are not automatically deleted.
Deleting stage work removes its native Nextflow cache-resume capability; a
failed or unverified stage retains its work. Linux/native Nextflow is required;
Windows users may use an already configured WSL/Linux execution environment.

## 2. Reference, library and supported design

| Parameter | Default | Effect / restriction |
| --- | --- | --- |
| `--reference_build` | unset | Declared genome-build label; preferred over `--genome` when supplied. Reference files are still explicitly required. |
| `--genome` | unset | Fallback build label, not an automatic reference catalogue/download selector. |
| `--annotation_release` | unset | Declared annotation release recorded with reference provenance. |
| `--read_length` | unset | Positive nominal length fallback for a sample without its own value. All study rows must resolve to one nominal length. |
| `--allow_variable_read_length` | `true` | BAM audit permits observed positive primary lengths no greater than the declared nominal length; `false` requires exact equality. It does not allow mixed nominal study lengths. |
| `--design` | `~ condition` | Only this independent two-level model is wired across the complete bundle. Paired, blocking, batch-adjusted and interaction formulas are rejected. |
| `--analysis_purpose` | `analysis` | `analysis` or `engine_smoke`. `engine_smoke` is for explicitly non-biological engine fixtures; it cannot establish valid replication or release scientific results. Use `analysis` for research. |

`subject`/`batch` may record information without fitting it. Repeated subjects
across arms do not silently enable pairing. An extended metadata adapter or an
already fitted imported model has a separate contract; it does not expand the
raw complete route. See [design support](development/contracts/selected-contrasts-v1.md).

The upstream analysis-reference policy excludes annotated pseudogenes and
mitochondrial-encoded genes before fitting, preserving the full original
annotation and the decisions. No public top-level exclusion toggle is currently
exposed by this raw DAG. Changing a report filter is not a replacement for
that native analysis-reference policy. Existing historical results retain their
original fitting universe.

## 3. Raw-read QC, alignment and quantification

| Parameter | Default | Effect |
| --- | --- | --- |
| `--skip_trimming` | `false` | Skip Trim Galore for FASTQ reads. BAMs are not trimmed. |
| `--skip_fastqc` | `false` | Skip native FastQC only in upstream-only. Complete mode rejects `true`, because its fixed report requires measured sequencing QC. |
| `--stranded_threshold` | `0.8` | `auto` only: minimum normalized sense or antisense fraction for a forward/reverse call. Matches the nf-core classification default. |
| `--unstranded_threshold` | `0.1` | `auto` only: maximum absolute difference between normalized sense and antisense fractions for an unstranded call, after checking stranded calls. Inclusive boundary. |
| `--strandedness_sample_size` | `200000` | `auto` only: RSeQC native maximum usable gene-overlapping alignment records sampled in BAM order. These are reads/alignment records, not paired fragments or random subsamples. |
| `--strandedness_mapq` | `30` | `auto` only: native RSeQC minimum mapping quality for eligible alignments. This is independent of DEXSeq `--alignment_quality`. |
| `--star_index` | unset | Prebuilt compatible STAR index directory. If absent, FASTQ input triggers index construction; a BAM-only study does not build an unused index. |
| `--star_sjdb_overhang` | `99` | STAR genome-generation splice-junction overhang. Choose explicitly for the intended reference/library; it is not automatically `read_length-1`. |
| `--star_genome_sa_index_nbases` | unset | Optional STAR index setting, useful when a reference requires a non-default value. |
| `--salmon_seq_bias` | `true` | Enable Salmon sequence-bias correction in the retained genome-projection quantification route. |
| `--salmon_gc_bias` | `true` | Enable Salmon GC-bias correction. |
| `--salmon_archive` | pinned Salmon 2.8.0 Linux x86-64 archive | Exact executable archive for genome-BAM projection; see `pipeline/nextflow.config` for the full URL. |
| `--salmon_archive_sha256` | `10a7f7994df7ef0592ac191cc79a1abd41316463000b307647a9dcab21e1efa6` | Required identity check for that archive; changing a URL requires its actual matching hash and runtime qualification. |

The current DAG selects STAR and Salmon genome projection; it does not expose
another pipeline's `--aligner` or arbitrary quantifier selection. Retained STAR
alignment arguments include two-pass Basic, unsorted BAM output followed by
Samtools sort/index, seed 0, NH/HI/AS/NM/MD attributes and unmapped records within
BAM. Custom arguments require review of the actual process `ext.args` contract;
blindly replacing it can remove required data.

`auto` is a per-sample `strandedness` value, not a global `--strandedness` flag.
Only those rows run RSeQC 5.0.3 after the existing genome-BAM audit/sort/index;
FASTQ rows use their STAR BAM and BAM rows are not realigned. gffread 0.12.1
builds BED12 from the matched **unfiltered** GTF, independently of the filtered
native-analysis GTF. The sample-size/MAPQ defaults are RSeQC's native defaults;
the direction thresholds follow nf-core's classification convention.

Classification uses `sense/(sense+antisense)` and the corresponding reverse
fraction, retaining the failed-to-determine fraction separately. A side at or
above `stranded_threshold` wins; otherwise a difference at or below
`unstranded_threshold` calls unstranded. Missing, invalid or undetermined
evidence fails rather than falling back. Explicit declarations remain unchanged.
Mixed explicit/`auto` requests are allowed before inference, but every selected
comparison must have consistent **resolved** strandedness before native
quantification/counting/testing begins. Requested/resolved settings and native
RSeQC stdout/stderr are retained in pipeline QC and bundle provenance.
See [the inference/tool pins](../pipeline/modules/local/infer_strandedness.provenance.json)
and [nf-core's classification function](https://github.com/nf-core/rnaseq/blob/3.19.0/subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness/main.nf#L19).
Thresholds require `0.5 < stranded_threshold <= 1` and
`0 <= unstranded_threshold < 2*stranded_threshold-1`, so the classification
regions cannot overlap. The sample budget must be a positive integer and MAPQ
an integer from 0 to 255. See [native auto-inference validation and output paths](validation/auto-strandedness-20261025.md).

## 4. rMATS native testing

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--rmats_splice_diff_cutoff` | `0.0001` | Native rMATS `cstat` test setting. This is **not** Katana's observed `|delta PSI| >= 0.10` selection. |
| `--rmats_novel_splice_site` | `false` | Enable rMATS novel splice-site behavior when explicitly selected. |
| `--rmats_min_intron_len` | `50` | Native minimum intron length option. |
| `--rmats_max_exon_len` | `500` | Native maximum exon length option. |

Five event-class schemas are retained: SE, A3SS, A5SS, MXE and RI, each with JC
and JCEC native tables. JC/JCEC selection belongs in the project template and
network configuration, not a nonexistent top-level flag. Treatment is group 1;
native replicate vectors keep explicit sample order. Paired statistics are not
enabled by adding a column to the contrastsheet. Variable observed-read-length
handling is separately audited and passed to the pinned implementation.

## 5. DEXSeq DEU and DTU filtering

| Parameter | Default | Native use |
| --- | --- | --- |
| `--alignment_quality` | `10` | Alignment-quality setting for exon counting. |
| `--aggregation` | `true` | DEXSeq annotation/count aggregation setting. |
| `--n_dexseq_plot` | `10` | Native plots requested; complete native tables are retained. This is unrelated to report Top 10. |
| `--min_samps_feature_expr` | `2` | DTU filter: minimum samples meeting the feature-expression criterion. |
| `--min_samps_feature_prop` | `2` | DTU filter: minimum samples meeting the feature-proportion criterion. |
| `--min_samps_gene_expr` | `4` | DTU filter: minimum samples meeting the gene-expression criterion. |
| `--min_feature_expr` | `10` | Feature-expression cutoff used by the retained DTU filtering. |
| `--min_feature_prop` | `0.1` | Feature-proportion cutoff used by the retained DTU filtering. |
| `--min_gene_expr` | `10` | Gene-expression cutoff used by the retained DTU filtering. |

These native eligibility filters are distinct from downstream significance
selection. Preserve DEXSeq exon/bin units, DTU gene/transcript units and original
stageR adjustment targets. Imported stageR values are not arbitrary-target
q-values; a changed target may require repeating stage-wise adjustment. The
generic importer's remaining target-compatibility limitation is explicitly
recorded in [THRESHOLDS.md](THRESHOLDS.md).

## 6. SUPPA2 options

| Parameter | Default | Role |
| --- | --- | --- |
| `--suppa_per_local_event` | `true` | Generate/use local event analyses. |
| `--suppa_per_isoform` | `true` | Generate/use isoform PSI analyses. |
| `--generateevents_pool_genes` | `true` | Pool overlapping genes in the retained event generator. |
| `--generateevents_event_type` | `SE SS MX RI FL` | SUPPA event-type string; retain shell quoting when overriding. |
| `--generateevents_boundary` | `S` | Event-boundary definition. |
| `--generateevents_threshold` | `10` | Native event-generation threshold. |
| `--generateevents_exon_length` | `100` | Native exon-length setting. |
| `--psiperevent_total_filter` | `0` | Total-expression filter for event PSI. |
| `--diffsplice_local_event` | `true` | Run local-event differential analysis. |
| `--diffsplice_isoform` | `true` | Run isoform differential analysis within the enabled isoform branch. |
| `--diffsplice_method` | `empirical` | Native SUPPA differential method. |
| `--diffsplice_area` | `1000` | Native empirical-method area setting. |
| `--diffsplice_lower_bound` | `0` | Native lower-bound setting. |
| `--diffsplice_gene_correction` | `true` | Native gene-level correction option. |
| `--diffsplice_alpha` | `0.05` | Native differential-analysis alpha. |
| `--diffsplice_median` | `false` | Native median option. |
| `--diffsplice_tpm_threshold` | `0` | Native TPM eligibility threshold. |
| `--diffsplice_nan_threshold` | `0` | Native missing-value threshold. |
| `--clusterevents_local_event` | `false` | Optional local-event clustering. |
| `--clusterevents_isoform` | `false` | Optional isoform clustering. |
| `--clusterevents_sigthreshold` | unset | Optional native significance threshold for clustering. |
| `--clusterevents_dpsithreshold` | `0.05` | Clustering delta-PSI setting. |
| `--clusterevents_eps` | `0.05` | Clustering distance-neighborhood setting. |
| `--clusterevents_metric` | `euclidean` | Native clustering metric. |
| `--clusterevents_separation` | unset | Optional native clustering separation. |
| `--clusterevents_min_pts` | `20` | Native clustering minimum points. |
| `--clusterevents_method` | `DBSCAN` | Native clustering method. |

These are native SUPPA settings. They do not change the rMATS universe, the
fixed AS read-support gate or CoSplice's WGCNA module detection. The current
report's principal DAS layer uses its declared rMATS sources.

## 7. Katana `flow.json`: selected calculations

| Field | Default | Meaning |
| --- | --- | --- |
| `schema` | `katana.flow.config/v1` | Explicit schema identity recommended. |
| `contrast_ids` | All **declared** project contrasts | Nonempty, unique selected subset. Not all possible pairwise comparisons. |
| `event_set` | `stable_significant` | RS-DAS; alternative `significant` is ordinary DAS. Report switch preserves both views where supported. |
| `parameters` | `{}` | Map from exact workflow ID to that workflow's recognized runtime parameters. |
| `qc` | `{}` | Optional gene-body coverage configuration, used only with `analyses.coverage: true`. `reference_bed` can name matching BED12; otherwise generated from matching GTF. `rseqc_executable` defaults to `geneBody_coverage.py`. |
| `network` | Per-contrast network configuration | See next section. |
| `design_factors` | Not automatically inferred | Explicit factor declaration. A multifactor study requires matching fitted-interaction definitions for the strict complete report. |
| `interactions` | `[]` | Dedicated interaction entries with `id`, `design_path`, `reference_a`, `reference_b` and supported workflow parameters. Actual count inputs and model eligibility are required. |
| `mapp` | Project-derived native MAEI | `adapter: "mapp.maei"`; `inputs`, `parameters`, optional `prepare_parameters`. Empty inputs resolve the packaged author source and actual project quantification. |

The selected plan includes audit/read support, AS landscape, molecular evidence,
selected cross-comparison analysis, evidence-specific enrichments/RBP identity,
per-comparison MotifTrace/boundary/CLIP, sequencing QC, independent networks
and the native MAPP family. The extra gene-body BAM scan is opt-in. Optional
modules can be disabled through `analyses` in the flow configuration; disabled
modules are not scheduled and are recorded as not selected in the report.
See [analysis selection](ANALYSIS_AND_EXECUTION_OPTIONS.md).

### Differential thresholds and evidence selection

| Choice | Location | Accepted default / boundary |
| --- | --- | --- |
| AS read support | Fixed policy `bulk-all-samples-counts-v1` | Finite actual `IJC + SJC >= 10` for **every expected sample** in each arm; no mean/pooled/subset substitution. Not an adjustable runtime preset. |
| DAS significance | `parameters.quality.stable_as` | `fdr: 0.05`, `min_abs_dpsi: 0.10`; both runtime configurable. The `differential.as` workflow consumes this completed Quality result, rather than owning separate threshold flags. Native P/FDR remain unchanged. |
| DEG selection | `parameters.differential.deg` | `fdr: 0.05`, `min_abs_log2fc: 1.0`. Runtime configurable; classifies imported results without refitting DESeq2. |
| DTU/DEXSeq selection | `parameters.differential.dtu` | `fdr: 0.05`; retains endpoint/adjustment semantics. |
| Integrated evidence | `parameters.prioritize.evidence` | `fdr: 0.05`, `deg_min_abs_log2fc: 1.0`; keep consistent with the intended separate DEG selection. |
| Enrichment foreground | Set by the flow planner per task | `deg`, `das`, `dtu`, `dexseq` separately; their detected/testable backgrounds remain separate. |
| Enrichment display | `parameters.prioritize.enrichment` | `minimum_overlap: 2` is display-only, `top_terms: 12`, `ontology_category: "All"`, `bubble_x: "gene_ratio"`. BH family is not recalculated after display filtering. |

The report defaults to RS-DAS. Motif/CLIP selection is explicitly RS-DAS even
when `event_set` is ordinary DAS. A network uses its complete measurable,
count-supported, nonconstant event universe, **not only significant DAS**.
Complete raw results and original effects/P/FDR remain available.

### Networks

| `flow.json` field | Default | Effect |
| --- | --- | --- |
| `network.scope` | `contrast` | Separate discovery/association/annotation for each selected comparison. Explicit `cohort` exists as an alternative and changes scientific scope. |
| `network.config` | unset | Optional external engine YAML; must match the actual project/scope. |
| `network.configs` | `{}` | Optional map from selected contrast IDs to exact scoped YAMLs; incompatible with cohort scope. |
| `network.bootstrap.enabled` | `false` | Off by default; enable explicitly to schedule real resampling, not just a label. |
| `network.parameters.minimum_network_samples` | `6` for contrast, `15` for cohort | Configurable integer >=3; total independent samples in the selected network, not per group. Actual n and the requested minimum are recorded. |
| `network.parameters.rmats_mode` | Explicit request, then project metadata, then unique asset mode | `JC` or `JCEC`; ambiguous dual-mode sources fail rather than silently choose. |
| `network.parameters.min_module_size` | `20` for generated configuration | Events, not genes. Omitted/null preserves explicit external YAML settings. |
| `network.parameters.merge_cut_height` | `0.25` for generated configuration | Eigengene-dissimilarity merging; omitted/null preserves explicit external YAML settings. |
| `network.parameters.hub_metric` | `abs_kME` | Retain signed kME and full ranking; `legacy_composite` is an explicitly different optional metric. |
| `network.parameters.apply_threshold_labels` | `false` | Opt-in descriptive labels only; no threshold-based member deletion/refit. |
| `network.parameters.max_events` | `0` | All eligible events. Positive cap is an explicit variance-ranked subset; not the default. |
| `network.parameters.bootstrap_iterations` | `50` only when enabled | Integer at least 2 for an explicitly enabled bootstrap. No default bootstrap work when disabled. |
| `network.parameters.tom_dense_max_events` | `null` | Legacy compatibility field; no event-count protection is enforced. |
| `network.parameters.wgcna_timeout_seconds` | `null` | Legacy compatibility field; no internal R-WGCNA timeout. Scheduler time remains user-configurable. |
| `network.parameters.rscript` | Runtime `Rscript` / declared path | Actual executable in the selected environment. |

Generated contrast networks default to six distinct independent biological
samples in total; cohort scope defaults to 15. Both are configurable through
`network.parameters.minimum_network_samples` (integer >=3). Scope identities and actual sample count are recorded.
The production engine retains observed-PSI WGCNA/TOM, `bicor`, signed network
and soft power 6 unless a compatible explicit engine configuration changes
them. Scientific defaults and full rankings remain separate from optional
descriptive labels. Resolved engine configuration is the authoritative record.

### MotifTrace and boundary/CLIP

| Runtime field | Default | Purpose |
| --- | --- | --- |
| `parameters.prioritize.motiftrace_consequences.model_ids` | empty | All supplied models; optional predeclared native IDs/hashes. No fixed model-count cap or silent first-model selection. |
| `...custom_motif_path` | empty | Optional species-matched consensus TSV/Minimal MEME with required metadata, replacing the public library for that calculation. |
| `...cluster_permutations` | `999` | Whole-pair swaps; integer >=99, with no fixed upper count. Maximum-cluster P is FWER, not enrichment FDR. |
| `...max_tests` | `null` | Deprecated compatibility field for older parameter files; does not limit the number of tests. Omit from new configurations. |
| `...max_matrix_cells` | `null` | Deprecated compatibility field for older parameter files; does not limit paired-matrix elements. Omit from new configurations. |
| `parameters.prioritize.rbp_mechanism.sequence_scope` | `boundary_association` | Retained boundary analysis mode; completed-map/parallel/consequence modes have separate inputs and contracts. |
| `...window_bp` | `100` | Bases on each side of the splice boundary, same convention for matched foreground/background; allowed 10–1000. |
| `...site_pvalue` | `0.0001` | FIMO site screening, not enrichment FDR; allowed `1e-8`–`0.01`. |
| `...motif_evidence` | `all_evidence_separate` | Direct, homology and predicted evidence retained separately; `direct_experiment` selects the direct evidence subset. |
| `...binding_mode` | Planner forces `required` | Strict flow requires supplied matched CLIP evidence; standalone workflow default `if_configured` is not the complete-flow setting. |
| `parameters.prioritize.motiftrace_parallel.parallel_bundle_asset_id` | empty | Optional exact completed/calibrated parallel-family asset; no new calibration inferred from native MAPP output alone. |

MotifTrace has no fixed application cap on model counts, RNA size,
custom-library file size/row counts, test counts, paired-matrix elements,
or permutation-by-gene-pair elements. Workload estimates are informational.
Configure actual CPU, memory and wall-time requests through the executor/site
configuration; the legacy fields above do not allocate resources. Complete
event/model/window families and the established P/FDR/FWER calculations are
retained.

`...` means the same workflow prefix as the preceding full field, not a JSON
key. Use actual nested keys as shown in the
[configuration example](NEXTFLOW_USAGE.md#6-downstream-analysis-configuration).

The report presents **up to ten significant** sequence models with their native
statistical identities, motif sequences and recorded gene contributors. It
does not force ten discoveries when fewer pass. Full model/window/cluster
families remain in files. Native P, BH FDR and cluster FWER are distinct; some
models legitimately share an adjusted value. The presentation keeps one
recorded native window/full-sequence test per model and ranks eligible models
by that test's `q_new_window_family`, then paired raw P, then model ID. This
recorded test's FDR is **not a newly fitted whole-model FDR**. Fewer than 20 independent matched
pairs yields non-estimable new region/window inference rather than fabricated
significance. See [MOTIFTRACE.md](MOTIFTRACE.md).

## 8. Execution environment and per-job resources

| Parameter | Default | Effect |
| --- | --- | --- |
| `--katana_container` | unset | Required for complete/downstream container profiles; matching full runtime image/SIF. |
| `--dexseq_tools_container` | unset | Required for raw FASTQ/BAM container execution; separate verified DEXSeq tools image/SIF. |
| `--katana_stage_containers` | `{}` | Advanced map of execution groups to compatible images. Main control still requires the full image. Use a structured params file, not an unvalidated shell string. |
| `--flow_executable` | `katana-flow` | Container profiles permit only `katana-flow` or `/opt/conda/bin/katana-flow`; mounted host launchers are rejected. |
| `--slurm_account` | unset | Account passed to Slurm when configured. |
| `--slurm_partition` | unset | Slurm queue/partition. |
| `--max_cpus` | unset | Optional per-task CPU ceiling, not a global concurrency cap. |
| `--max_memory` | unset | Optional per-task memory ceiling, with unit e.g. `128.GB`/`128 GB` as accepted by Nextflow. |
| `--max_time` | unset | Optional per-task time ceiling, e.g. `24.h`/`24 h`. |

Recognized image execution-group keys are `control`, `report`, `audit`,
`quality`, `as`, `deg`, `dtu`, `landscape`, `evidence`, `enrichment`,
`motiftrace`, `overlap`, `reference`, `qc`, `coverage`, `network_fit`,
`network_followup`, `interaction`, `mapp_prepare`, `mapp`. Unknown keys or empty
image references fail early.

Profiles are `local`, `docker`, `singularity`, `apptainer`, `conda`, `slurm`.
Combine scheduler and engine, for example `-profile slurm,apptainer`. Do not
mix engines. The Conda profile enables module recipes; install the full native
environment explicitly for Katana standalone/native execution.

Important first-attempt resource requests:

| Process class | CPUs | Memory | Time |
| --- | --- | --- | --- |
| General upstream | 1 | 6 GB | 4 h |
| Upstream `process_low` | 2 | 12 GB | 4 h |
| Upstream `process_medium` | 6 | 36 GB | 8 h |
| Upstream `process_high` | 12 | 72 GB | 16 h |
| STAR index override | 16 | 96 GB; 128 GB on allowed retry | Inherited matching process-label time; inspect resolved task. |
| Katana control | 1 | 4 GB | 2 h |
| Katana downstream task | task `spec.resources.cpus`, else 1 | task setting, else 8 GB | task setting, else 12 h |
| Katana report | 2 | 12 GB | 8 h |

Actual process labels, selectors and task specs can override general defaults.
Check the trace and scheduler allocation. Never assume a theoretical override
was matched. Upstream eligible resource failures have at most one configured
retry; downstream tasks have `errorStrategy = 'finish'` without automatic retry,
allowing already submitted independent work to finish for a later `-resume`.
The Slurm profile has unbounded executor submission queue (`queueSize = 0`),
with admission controlled by Slurm/site limits. Native tool threads should
match each task allocation.

## 9. Deep-dive options

These are `katana-targets` arguments, **not** Nextflow pipeline params:

| Command/argument | Effect |
| --- | --- |
| `export --project --selections --outdir` | Create request from exact existing source artifact/row identities. |
| Repeated export `--contrast-id` | Explicit declared comparison subset; does not enumerate all pairs. |
| Repeated export `--module` | `gene_card`, `event_card`, `sashimi`. |
| `inspect --request` | Verify/show request and required logical resources. |
| `plan --request --bindings [--outdir]` | Check target evidence/resources/cache readiness without analysis. |
| `run --request --bindings --outdir` | Execute/reuse requested targets, producing independent reports. |
| Repeated run `--target-id` | Run exact target subset; use distinct output directories for independent scheduler tasks. |
| Run `--workers` | Optional local budget; default derives from the available allocation. |
| `tasks --request --format tsv` | List exported target IDs for scheduler integration. |
| `results --outdir` | Verify completed target products and manifest. |
| `restore --results --outdir [--project-id]` | Verify returned ZIP/directory and create a new resumable workspace; no refitting. |

Requests contain completed evidence and resource identities, not large BAMs.
Binding paths are machine-specific; identities are portable. Missing resources
produce partial/blocked states, not successful empty Sashimi. Exit 0 means that
command succeeded; target exit 2 means not ready or partial/failed requested
outputs. Inspect per-target/module statuses. See
[TARGET_COMMANDS.md](TARGET_COMMANDS.md) for complete contracts and examples.

## 10. Reproducibility record and references

Keep source revision, matching wheel/image hashes, effective params, samples,
selected contrasts, resource versions/hashes, JC/JCEC and original tool
directions. New resource/settings identities require the corresponding new
analysis identity. Preserve old results, failed logs and cache records.

Primary implementation references:

- [Pipeline defaults](../pipeline/nextflow.config), [entry validation](../pipeline/main.nf), [raw metadata validator](../pipeline/bin/katana_validate_inputs.py).
- [Upstream DAG](../pipeline/workflows/katana_core.nf), [base resources](../pipeline/conf/upstream_base.config), [process arguments](../pipeline/conf/upstream_modules.config).
- [Downstream planner](../katana_next/flow/planner.py), [native importer](../katana_next/flow/importer.py), [network scope](../katana_next/flow/network_inputs.py).
- [Threshold policies and cited scientific basis](THRESHOLDS.md), [network overrides](../katana_next/workflows/programs.py), [MotifTrace definitions](../katana_next/workflows/motiftrace_consequences.py), [boundary/CLIP definitions](../katana_next/workflows/rbp_mechanism.py).

Official documentation structure and input conventions were checked against
[rnaseq](https://nf-co.re/rnaseq/usage),
[rnasplice](https://nf-co.re/rnasplice/usage),
[rnasplice parameters](https://nf-co.re/rnasplice/parameters), and
[differentialabundance](https://nf-co.re/differentialabundance/usage).
Their complete current parameter lists are not interchangeable with Katana's.

## Optional downstream analyses and resource-only resume

`flow.json` accepts `analyses` with boolean switches for `landscape`,
`cross_comparison`, `enrichment`, `rbp_identity`, `motiftrace`, `motif_binding`,
`mapp`, `coverage`, `sequencing_qc`, `network` and `interaction`. Every switch defaults to `true` except `coverage`, which defaults to `false` (extra RSeQC BAM scan is opt-in).
A `false` value omits its tasks and report content. Core design/QC and molecular
results remain. See [selection, resource and resume examples](ANALYSIS_AND_EXECUTION_OPTIONS.md).

The managed `katana-flow run` command accepts `--resource-config FILE.json`.
This execution-only file can change during `--resume`; scientific configuration,
inputs and arbitrary site configurations retain their normal identity checks.
