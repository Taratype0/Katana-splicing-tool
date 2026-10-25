# Katana Nextflow: parameters and analysis choices

This reference describes the **current Katana implementation**, not every
option of `nf-core/rnaseq`, `rnasplice` or `differentialabundance`. For formats
and commands, use [NEXTFLOW_USAGE.md](NEXTFLOW_USAGE.md). Defaults below were
checked against `pipeline/nextflow.config`, the public validator, importer and
downstream planner for the accepted 25 October 2026 release.

The four strandedness-inference parameters below are a **development-workflow
addition** on `codex/package-development`. Use that checkout's
`pipeline/main.nf`, as described in [the usage guide](NEXTFLOW_USAGE.md#2-installation-and-directory-layout).
The frozen `v2.0.0-beta.5.dev20261025` wheel/export still requires explicit
strandedness and does not gain these parameters from a documentation update.

## Configuration layers

| Layer | Where it belongs | Examples |
| --- | --- | --- |
| Nextflow launcher | Single-dash Nextflow arguments | `-profile`, `-params-file`, `-c`, `-work-dir`, `-resume`. |
| Pipeline parameters | `--parameter` or a Nextflow YAML/JSON params file | `--mode`, `--input`, `--contrasts`, reference paths and rMATS settings. |
| Study/resources | `project-template.json` or an imported v2 `project.json` | Species, exact samples/contrasts, resources, JC/JCEC and original result provenance. |
| Downstream plan | `flow.json`, schema `katana.flow.config/v1` | Selected contrasts, RS-DAS, per-module analysis parameters, independent networks, QC and MAPP. |
| Site execution | Nextflow `-c site.config` | Fully qualified process resource selectors, account policy and executor configuration. |
| Deep dive | Exported request plus `bindings.json` and target CLI arguments | Chosen source row IDs, target scope, BAM/index locations and worker budget. |

Do not invent `--deg_fdr`, `--rs_das`, `--motif_top10` or `--rmats_mode` flags.
Those are not public top-level pipeline parameters. Use the supported
configuration field or the appropriate command below.

Example `inputs/params.yaml` for `nextflow ... -params-file inputs/params.yaml`:

```yaml
mode: complete
input: inputs/samples.csv
contrasts: inputs/contrasts.csv
fasta: references/genome.fa
gtf: references/annotation.gtf
reference_build: YOUR_BUILD
annotation_release: YOUR_RELEASE
project_template: inputs/project-template.json
flow_config: inputs/flow.json
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
| `--project_template` | unset | Resource/metadata template, required for complete mode. |
| `--project` | unset | Katana v2 project manifest; required for downstream-only. |
| `--flow_config` | unset | Downstream configuration; required for complete, and downstream-only unless reusing `--flow_plan`. |
| `--flow_plan` | unset | Original sealed executed downstream plan; only downstream-only. Bypasses replanning, so an additional `--flow_config` is not applied. |
| `--outdir` | `results` | Published output root; keep outside input project/config roots. |
| `--publish_dir_mode` | `copy` | Upstream general publication mode. Some processes explicitly retain `copy` for portable deliverables; this is not a global storage-deduplication switch. |
| `--flow_report_revision` | `1` | Explicit presentation-only report task revision. Letters/digits/dot/underscore/hyphen; no scientific threshold change. |
| `--help` | `false` | Print concise pipeline help. This guide is the full parameter reference; an nf-core schema-generated launcher is not claimed. |

The required raw-sheet headers are `sample,condition,strandedness`, plus actual
FASTQ/BAM fields. The validator supplements these with explicit resolved
`source,assay,single_end,read_length,subject,batch` fields. See the complete
[samplesheet and contrast specification](NEXTFLOW_USAGE.md#3-samplesheet-input).

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
| `--stranded_threshold` | `0.8` | Development `auto` only: minimum normalized sense or antisense fraction for a forward/reverse call. Matches the nf-core classification default. |
| `--unstranded_threshold` | `0.1` | Development `auto` only: maximum absolute difference between normalized sense and antisense fractions for an unstranded call, after checking stranded calls. Inclusive boundary. |
| `--strandedness_sample_size` | `200000` | Development `auto` only: RSeQC native maximum usable gene-overlapping alignment records sampled in BAM order. These are reads/alignment records, not paired fragments or random subsamples. |
| `--strandedness_mapq` | `30` | Development `auto` only: native RSeQC minimum mapping quality for eligible alignments. This is independent of DEXSeq `--alignment_quality`. |
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
| `qc` | `{}` | Gene-body coverage configuration. `reference_bed` can name matching BED12; otherwise generated from matching GTF. `rseqc_executable` defaults to `geneBody_coverage.py`. |
| `network` | Per-contrast network configuration | See next section. |
| `design_factors` | Not automatically inferred | Explicit factor declaration. A multifactor study requires matching fitted-interaction definitions for the strict complete report. |
| `interactions` | `[]` | Dedicated interaction entries with `id`, `design_path`, `reference_a`, `reference_b` and supported workflow parameters. Actual count inputs and model eligibility are required. |
| `mapp` | Project-derived native MAEI | `adapter: "mapp.maei"`; `inputs`, `parameters`, optional `prepare_parameters`. Empty inputs resolve the packaged author source and actual project quantification. |

The complete plan includes audit/read support, AS landscape, molecular evidence,
selected cross-comparison analysis, evidence-specific enrichments/RBP identity,
per-comparison MotifTrace/boundary/CLIP, sequencing/gene-body QC, independent
networks and the native MAPP family. It is not a universal collection of
`skip_*` flags. To run a subset as a module, use its standalone CLI/API and
retain its partial scope rather than claim complete-report acceptance.

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
| `network.parameters.rmats_mode` | Explicit request, then project metadata, then unique asset mode | `JC` or `JCEC`; ambiguous dual-mode sources fail rather than silently choose. |
| `network.parameters.min_module_size` | `20` for generated configuration | Events, not genes. Omitted/null preserves explicit external YAML settings. |
| `network.parameters.merge_cut_height` | `0.25` for generated configuration | Eigengene-dissimilarity merging; omitted/null preserves explicit external YAML settings. |
| `network.parameters.hub_metric` | `abs_kME` | Retain signed kME and full ranking; `legacy_composite` is an explicitly different optional metric. |
| `network.parameters.apply_threshold_labels` | `false` | Opt-in descriptive labels only; no threshold-based member deletion/refit. |
| `network.parameters.max_events` | `0` | All eligible events. Positive cap is an explicit variance-ranked subset; not the default. |
| `network.parameters.bootstrap_iterations` | `50` only when enabled | Integer at least 2 for an explicitly enabled bootstrap. No default bootstrap work when disabled. |
| `network.parameters.tom_dense_max_events` | Engine `3000` | Positive integer computation/storage strategy limit; not a significant-event threshold. |
| `network.parameters.wgcna_timeout_seconds` | Engine `1800` | Native WGCNA subprocess timeout; distinct from scheduler wall time. |
| `network.parameters.rscript` | Runtime `Rscript` / declared path | Actual executable in the selected environment. |

Generated contrast networks require six distinct independent biological samples;
cohort scope requires 15. Scope identities and actual sample count are recorded.
The production engine retains observed-PSI WGCNA/TOM, `bicor`, signed network
and soft power 6 unless a compatible explicit engine configuration changes
them. Scientific defaults and full rankings remain separate from optional
descriptive labels. Resolved engine configuration is the authoritative record.

### MotifTrace and boundary/CLIP

| Runtime field | Default | Purpose |
| --- | --- | --- |
| `parameters.prioritize.motiftrace_consequences.model_ids` | empty | All supplied models; optional predeclared native IDs/hashes. Up to 256 models under the accepted implementation; no silent first-model selection. |
| `...custom_motif_path` | empty | Optional species-matched consensus TSV/Minimal MEME with required metadata, replacing the public library for that calculation. |
| `...cluster_permutations` | `999` | Whole-pair swaps; allowed 99–9999. Maximum-cluster P is FWER, not enrichment FDR. |
| `...max_tests` | `1000000` | Complete-family computation budget; exceeding it fails without truncation. |
| `...max_matrix_cells` | `8000000` | Paired float64 matrix budget; no removal of matched pairs to fit it. |
| `parameters.prioritize.rbp_mechanism.sequence_scope` | `boundary_association` | Retained boundary analysis mode; completed-map/parallel/consequence modes have separate inputs and contracts. |
| `...window_bp` | `100` | Bases on each side of the splice boundary, same convention for matched foreground/background; allowed 10–1000. |
| `...site_pvalue` | `0.0001` | FIMO site screening, not enrichment FDR; allowed `1e-8`–`0.01`. |
| `...motif_evidence` | `all_evidence_separate` | Direct, homology and predicted evidence retained separately; `direct_experiment` selects the direct evidence subset. |
| `...binding_mode` | Planner forces `required` | Strict flow requires supplied matched CLIP evidence; standalone workflow default `if_configured` is not the complete-flow setting. |
| `parameters.prioritize.motiftrace_parallel.parallel_bundle_asset_id` | empty | Optional exact completed/calibrated parallel-family asset; no new calibration inferred from native MAPP output alone. |

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
