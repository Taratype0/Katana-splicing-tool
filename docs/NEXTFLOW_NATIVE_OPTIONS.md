# Native Nextflow options (test branch)

All 114 declared options from `pipeline/nextflow.config`.
Use these as `--name value` or keys in one optional `-params-file`.
Single-dash launcher options such as `-profile`, `-c`, `-work-dir`, and `-resume` are separate.

See [run examples and explanations](NEXTFLOW_SIMPLE.md), [scientific parameter details](NEXTFLOW_PARAMETERS.md),
and [the complete editable template](../examples/nextflow/params.complete.yaml).
The older wheel and `main` do not yet include the new flat downstream options.

| Parameter | Declared default |
| --- | --- |
| `--mode` | `"complete"` |
| `--input` | `null` |
| `--contrasts` | `null` |
| `--fasta` | `null` |
| `--gtf` | `null` |
| `--project` | `null` |
| `--project_template` | `null` |
| `--flow_config` | `null` |
| `--flow_plan` | `null` |
| `--flow_report_revision` | `"1"` |
| `--flow_executable` | `"katana-flow"` |
| `--katana_container` | `null` |
| `--dexseq_tools_container` | `null` |
| `--katana_stage_containers` | `{}` |
| `--outdir` | `"results"` |
| `--publish_dir_mode` | `"copy"` |
| `--storage_mode` | `"retain"` |
| `--managed_preflight` | `false` |
| `--help` | `false` |
| `--study_name` | `"Katana study"` |
| `--species` | `null` |
| `--gene_sets` | `null` |
| `--rbp_gene_sets` | `null` |
| `--motif_database` | `null` |
| `--clip_peaks` | `null` |
| `--fimo_runtime` | `null` |
| `--resource_metadata` | `{}` |
| `--contrast_ids` | `null` |
| `--rmats_mode` | `"JCEC"` |
| `--event_set` | `"stable_significant"` |
| `--deg_fdr` | `0.05` |
| `--deg_log2fc` | `1.0` |
| `--das_fdr` | `0.05` |
| `--das_dpsi` | `0.1` |
| `--dtu_fdr` | `0.05` |
| `--network_min_module_size` | `20` |
| `--network_merge_cut_height` | `0.25` |
| `--skip_landscape` | `false` |
| `--skip_cross_comparison` | `false` |
| `--skip_enrichment` | `false` |
| `--skip_rbp_identity` | `false` |
| `--skip_motiftrace` | `false` |
| `--skip_motif_binding` | `false` |
| `--skip_mapp` | `false` |
| `--skip_network` | `false` |
| `--skip_interaction` | `false` |
| `--skip_sequencing_qc` | `false` |
| `--with_genebody_coverage` | `false` |
| `--with_network_bootstrap` | `false` |
| `--genome` | `null` |
| `--reference_build` | `null` |
| `--annotation_release` | `null` |
| `--read_length` | `null` |
| `--allow_variable_read_length` | `true` |
| `--stranded_threshold` | `0.8` |
| `--unstranded_threshold` | `0.1` |
| `--strandedness_sample_size` | `200000` |
| `--strandedness_mapq` | `30` |
| `--design` | `"~ condition"` |
| `--analysis_purpose` | `"analysis"` |
| `--skip_trimming` | `false` |
| `--skip_fastqc` | `false` |
| `--star_index` | `null` |
| `--star_sjdb_overhang` | `99` |
| `--star_genome_sa_index_nbases` | `null` |
| `--salmon_archive` | `"https://github.com/COMBINE-lab/salmon/releases/download/v2.8.0/salmon-cli-x86_64-unknown-linux-gnu.tar.xz"` |
| `--salmon_archive_sha256` | `"10a7f7994df7ef0592ac191cc79a1abd41316463000b307647a9dcab21e1efa6"` |
| `--salmon_seq_bias` | `true` |
| `--salmon_gc_bias` | `true` |
| `--alignment_quality` | `10` |
| `--aggregation` | `true` |
| `--n_dexseq_plot` | `10` |
| `--min_samps_feature_expr` | `2` |
| `--min_samps_feature_prop` | `2` |
| `--min_samps_gene_expr` | `4` |
| `--min_feature_expr` | `10` |
| `--min_feature_prop` | `0.1` |
| `--min_gene_expr` | `10` |
| `--rmats_splice_diff_cutoff` | `0.0001` |
| `--rmats_novel_splice_site` | `false` |
| `--rmats_min_intron_len` | `50` |
| `--rmats_max_exon_len` | `500` |
| `--suppa_per_local_event` | `true` |
| `--suppa_per_isoform` | `true` |
| `--generateevents_pool_genes` | `true` |
| `--generateevents_event_type` | `"SE SS MX RI FL"` |
| `--generateevents_boundary` | `"S"` |
| `--generateevents_threshold` | `10` |
| `--generateevents_exon_length` | `100` |
| `--psiperevent_total_filter` | `0` |
| `--diffsplice_local_event` | `true` |
| `--diffsplice_isoform` | `true` |
| `--diffsplice_method` | `"empirical"` |
| `--diffsplice_area` | `1000` |
| `--diffsplice_lower_bound` | `0` |
| `--diffsplice_gene_correction` | `true` |
| `--diffsplice_alpha` | `0.05` |
| `--diffsplice_median` | `false` |
| `--diffsplice_tpm_threshold` | `0` |
| `--diffsplice_nan_threshold` | `0` |
| `--clusterevents_local_event` | `false` |
| `--clusterevents_isoform` | `false` |
| `--clusterevents_sigthreshold` | `null` |
| `--clusterevents_dpsithreshold` | `0.05` |
| `--clusterevents_eps` | `0.05` |
| `--clusterevents_metric` | `"euclidean"` |
| `--clusterevents_separation` | `null` |
| `--clusterevents_min_pts` | `20` |
| `--clusterevents_method` | `"DBSCAN"` |
| `--slurm_account` | `null` |
| `--slurm_partition` | `null` |
| `--max_cpus` | `null` |
| `--max_memory` | `null` |
| `--max_time` | `null` |
