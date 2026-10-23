include { VALIDATE_INPUTS } from '../modules/local/validate_inputs'
include { VALIDATE_REFERENCE } from '../modules/local/validate_reference'
include { STAGE_FASTQ } from '../modules/local/stage_fastq'
include { PREPARE_REFERENCES } from '../modules/local/prepare_references'
include { FILTER_ANALYSIS_GTF; ANALYSIS_TRANSCRIPTS } from '../modules/local/analysis_reference'
include { ADAPT_SUPPA_CONTRASTS } from '../modules/local/adapt_suppa_contrasts'
include { REFERENCE_FAIDX } from '../modules/local/reference_faidx'
include { CAPTURE_BAM_AUDIT } from '../modules/local/capture_bam_audit'
include { AUDIT_GENOME_BAM } from '../modules/local/audit_genome_bam'
include { COLLATE_GENOME_BAM } from '../modules/local/collate_genome_bam'
include { SALMON_GENOME_PROJECTION } from '../modules/local/salmon_genome_projection'
include { KATANA_EXPORT } from '../modules/local/katana_export'
include { KATANA_DIFFERENTIAL } from '../modules/local/katana_differential'

include { FASTQC } from '../modules/nf-core/fastqc/main'
include { FASTQC as BAM_FASTQC } from '../modules/nf-core/fastqc/main'
include { TRIMGALORE } from '../modules/nf-core/trimgalore/main'
include { STAR_GENOMEGENERATE } from '../modules/nf-core/star/genomegenerate/main'
include { STAR_ALIGN } from '../modules/nf-core/star/align/main'
include { SAMTOOLS_SORT } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index/main'
include { DEXSEQ_ANNOTATION } from '../modules/local/dexseq_annotation'
include { DEXSEQ_DEU } from '../subworkflows/local/dexseq_deu'
include { TX2GENE_TXIMPORT } from '../subworkflows/local/tx2gene_tximport'
include { DRIMSEQ_DEXSEQ_DTU } from '../subworkflows/local/drimseq_dexseq_dtu'
include { RMATS } from '../subworkflows/local/rmats'
include { SUPPA } from '../subworkflows/local/suppa'

def artifact(kind, path, relative_path, extra = [:]) {
    [extra + [kind: kind, relative_path: relative_path], path]
}

def qcArtifacts(meta, reports, kind, phase) {
    [reports].flatten().collect { p -> artifact('qc', p, "qc/${phase}/${meta.id}/${p.name}",
        [sample_id: meta.id, qc_type: kind, qc_phase: phase]) }
}

def tableArtifacts(files, kind, directory, pattern, extra = [:]) {
    [files].flatten().collect { p ->
        def match = p.name =~ pattern
        if (!match.matches()) error("Unexpected ${kind} result name: ${p.name}")
        artifact(kind, p, "${directory}/${p.name}", extra + [comparison_id: match[0][1]])
    }
}

def nativeDtuArtifacts(files, kind, directory, pattern, contrasts) {
    [files].flatten().collect { p ->
        def match = p.name =~ pattern
        if (!match.matches()) error("Unexpected ${kind} result name: ${p.name}")
        def pair_name = match[0][2]
        def cmp = contrasts.find { "${it.treatment}-${it.control}" == pair_name }
        if (!cmp) error("DTU result has no declared contrast: ${pair_name}")
        def direction = match[0][1] == 'DEXSeqResults' ? 'control_over_treatment_native_header_authoritative' : 'not_applicable'
        artifact(kind, p, "${directory}/${match[0][1]}.${cmp.contrast}.tsv", [comparison_id: cmp.contrast,
            native_fold_change_direction: direction, contrast_level_order: [cmp.treatment, cmp.control]])
    }
}

workflow KATANA_PIPELINE {
    take:
    ch_samplesheet
    ch_contrastsheet
    ch_fasta
    ch_gtf
    main:

    // Stage source parents for container-visible validation. CSV-relative paths
    // resolve against the original samplesheet, not an ephemeral task directory.
    ch_raw_files = ch_samplesheet.splitCsv(header: true).combine(ch_samplesheet).flatMap { row, sheet ->
        def input_parent = sheet.parent
        ['fastq_1', 'fastq_2', 'genome_bam', 'transcriptome_bam'].findAll { row[it] }.collect { key ->
            def candidate = row[key].toString()
            if (candidate ==~ /^[A-Za-z][A-Za-z0-9+.-]*:\/\/.*/) {
                error("Sample ${row.sample}: ${key} must be a local or mounted file path in v1; download ${candidate} first")
            }
            file(candidate.startsWith('/') || candidate =~ /^[A-Za-z]:[\\\/]/ ? candidate : input_parent.resolve(candidate).toString(), checkIfExists: true)
        }
    }.unique().collect()

    VALIDATE_INPUTS(ch_samplesheet, ch_contrastsheet, ch_raw_files, ch_fasta, ch_gtf,
        Channel.value(file("${projectDir}/bin/katana_validate_inputs.py", checkIfExists: true)),
        ch_samplesheet.map { it.parent.toString() },
        params.read_length ?: '', params.design, params.analysis_purpose)

    VALIDATE_REFERENCE(VALIDATE_INPUTS.out.audit, ch_fasta, ch_gtf,
        Channel.value(file("${projectDir}/bin/kanata_validate_reference.py", checkIfExists: true)),
        params.reference_build ?: params.genome ?: '', params.annotation_release ?: '')
    // Keep gated study singletons reusable for every sample/contrast.
    ch_checked_fasta = ch_fasta.combine(VALIDATE_REFERENCE.out.audit).map { fasta, audit -> fasta }.first()
    ch_checked_gtf = ch_gtf.combine(VALIDATE_REFERENCE.out.audit).map { gtf, audit -> gtf }.first()

    ch_rows = VALIDATE_INPUTS.out.samples.splitCsv(header: true).collect().map { rows ->
        rows.withIndex().collect { row, i -> row + [order_index: i] }
    }.map { [records: it] }.combine(VALIDATE_REFERENCE.out.audit).map { container, audit -> container.records }.first()
    ch_contrasts = VALIDATE_INPUTS.out.contrasts.splitCsv(header: true).collect().map { rows ->
        def repeated = rows.groupBy { "${it.treatment}-${it.control}" }.findAll { k, v -> v.size() > 1 }
        if (repeated) error('Repeated treatment/control pairs under different contrast names are not supported in complete mode')
        rows
    }
    ch_rows_container = ch_rows.map { [records: it] }
    ch_contrasts_container = ch_contrasts.map { [records: it] }
    ch_rmats_read_length = ch_rows.map { rows ->
        def values = rows.collect { it.read_length as int }.unique()
        if (values.size() != 1) error('Complete mode currently requires a single declared read_length across the study')
        values[0]
    }
    ch_sources = ch_rows.flatMap { rows -> rows }.map { row ->
        def meta = [id: row.sample, condition: row.condition, source: row.source, assay: row.assay,
            single_end: row.single_end.toBoolean(), strandedness: row.strandedness,
            read_length: row.read_length as int, order_index: row.order_index,
            subject: row.subject ?: '', batch: row.batch ?: '']
        if (row.source == 'fastq') {
            tuple(meta, row.fastq_2 ? [file(row.fastq_1), file(row.fastq_2)] : [file(row.fastq_1)])
        } else {
            if (!row.genome_bam) error("Sample ${row.sample} lacks a genome BAM; transcriptome BAM alone cannot support a complete bundle")
            tuple(meta, file(row.genome_bam))
        }
    }.branch { fastq: it[0].source == 'fastq'; bam: it[0].source != 'fastq' }

    PREPARE_REFERENCES(ch_checked_fasta, ch_checked_gtf)
    FILTER_ANALYSIS_GTF(PREPARE_REFERENCES.out.gtf,
        Channel.value(file("${projectDir}/bin/kanata_filter_analysis_gtf.py", checkIfExists: true)))
    ANALYSIS_TRANSCRIPTS(PREPARE_REFERENCES.out.fasta, FILTER_ANALYSIS_GTF.out.gtf)
    ch_versions = PREPARE_REFERENCES.out.versions
    ch_qc_files = Channel.empty()
    ch_qc_artifacts = Channel.empty()
    STAGE_FASTQ(ch_sources.fastq, Channel.value(file("${projectDir}/bin/kanata_stage_fastq.py", checkIfExists: true)))
    ch_fastq = STAGE_FASTQ.out.reads
    if (!params.skip_fastqc) {
        FASTQC(ch_fastq)
        // BAM-only entry measures the records supplied by the user. These
        // results are not raw-read, trimming or original-aligner QC.
        BAM_FASTQC(ch_sources.bam)
        ch_versions = ch_versions.mix(FASTQC.out.versions)
        ch_versions = ch_versions.mix(BAM_FASTQC.out.versions)
        ch_qc_files = ch_qc_files.mix(FASTQC.out.zip.map { meta, reports -> reports })
        ch_qc_files = ch_qc_files.mix(BAM_FASTQC.out.zip.map { meta, reports -> reports })
        ch_qc_artifacts = ch_qc_artifacts.mix(
            FASTQC.out.zip.flatMap { meta, reports -> qcArtifacts(meta, reports, 'fastqc', 'raw') },
            FASTQC.out.html.flatMap { meta, reports -> qcArtifacts(meta, reports, 'fastqc', 'raw') })
        ch_qc_artifacts = ch_qc_artifacts.mix(
            BAM_FASTQC.out.zip.flatMap { meta, reports -> qcArtifacts(meta, reports, 'fastqc', 'bam_records') },
            BAM_FASTQC.out.html.flatMap { meta, reports -> qcArtifacts(meta, reports, 'fastqc', 'bam_records') })
    }
    if (!params.skip_trimming) {
        TRIMGALORE(ch_fastq)
        ch_fastq = TRIMGALORE.out.reads
        ch_versions = ch_versions.mix(TRIMGALORE.out.versions)
        ch_qc_files = ch_qc_files.mix(TRIMGALORE.out.log.map { meta, reports -> reports },
            TRIMGALORE.out.zip.map { meta, reports -> reports })
        ch_qc_artifacts = ch_qc_artifacts.mix(
            TRIMGALORE.out.log.flatMap { meta, reports -> qcArtifacts(meta, reports, 'trim_galore', 'trimming') },
            TRIMGALORE.out.zip.flatMap { meta, reports -> qcArtifacts(meta, reports, 'fastqc', 'trimmed') })
    }

    if (params.star_index) {
        ch_star_index = Channel.value(file(params.star_index, checkIfExists: true))
    } else {
        // Do not build an unused STAR index for a BAM-only study.
        ch_fastq_reference = PREPARE_REFERENCES.out.fasta.combine(ch_rows_container)
            .filter { fasta, container -> container.records.any { it.source == 'fastq' } }.map { fasta, container -> fasta }
        STAR_GENOMEGENERATE(ch_fastq_reference, PREPARE_REFERENCES.out.gtf)
        ch_star_index = STAR_GENOMEGENERATE.out.index.first()
        ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
    }
    STAR_ALIGN(ch_fastq, ch_star_index, PREPARE_REFERENCES.out.gtf, false, 'ILLUMINA', '')
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)
    ch_qc_files = ch_qc_files.mix(STAR_ALIGN.out.log_final.map { meta, report -> report })
    ch_qc_artifacts = ch_qc_artifacts.mix(STAR_ALIGN.out.log_final.flatMap { meta, report -> qcArtifacts(meta, report, 'star', 'alignment') })
    ch_genome_input = STAR_ALIGN.out.bam.mix(ch_sources.bam)

    // Header/reference identity, layout and primary-record length are checked
    // before a BAM reaches any shared quantification or statistical module.
    REFERENCE_FAIDX(PREPARE_REFERENCES.out.fasta)
    CAPTURE_BAM_AUDIT(ch_genome_input, Channel.value(file("${projectDir}/bin/katana_bam_primary_metrics.awk", checkIfExists: true)))
    AUDIT_GENOME_BAM(CAPTURE_BAM_AUDIT.out.captured, REFERENCE_FAIDX.out.fai, REFERENCE_FAIDX.out.dict,
        Channel.value(file("${projectDir}/bin/katana_audit_bam.py", checkIfExists: true)), params.allow_variable_read_length)
    // The captured BAM retains a fixed staged basename as a symlink. Arbitrary
    // user filenames must never collide with native SAMTOOLS_SORT's output name.
    ch_audited_genome_input = CAPTURE_BAM_AUDIT.out.staged_bam.map { meta, bam -> tuple(meta.id, meta, bam) }
        .join(AUDIT_GENOME_BAM.out.audit.map { meta, audit -> tuple(meta.id, audit) })
        .map { id, meta, bam, audit -> tuple(meta, bam) }
    SAMTOOLS_SORT(ch_audited_genome_input)
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)
    COLLATE_GENOME_BAM(SAMTOOLS_SORT.out.bam)
    ch_versions = ch_versions.mix(CAPTURE_BAM_AUDIT.out.versions, SAMTOOLS_SORT.out.versions, SAMTOOLS_INDEX.out.versions, COLLATE_GENOME_BAM.out.versions)

    SALMON_GENOME_PROJECTION(COLLATE_GENOME_BAM.out.bam, PREPARE_REFERENCES.out.fasta,
        FILTER_ANALYSIS_GTF.out.gtf, Channel.value(file(params.salmon_archive)), params.salmon_archive_sha256)
    ch_versions = ch_versions.mix(SALMON_GENOME_PROJECTION.out.versions)

    // Count each sample once, then fit only the declared contrasts.
    DEXSEQ_ANNOTATION(FILTER_ANALYSIS_GTF.out.gtf, params.aggregation)
    DEXSEQ_DEU(FILTER_ANALYSIS_GTF.out.gtf, SAMTOOLS_SORT.out.bam, DEXSEQ_ANNOTATION.out.gff,
        VALIDATE_INPUTS.out.samples, VALIDATE_INPUTS.out.contrasts,
        params.n_dexseq_plot, params.aggregation, params.alignment_quality)
    TX2GENE_TXIMPORT(SALMON_GENOME_PROJECTION.out.results.map { meta, directory -> tuple(meta, [directory]) }, FILTER_ANALYSIS_GTF.out.gtf)
    DRIMSEQ_DEXSEQ_DTU(TX2GENE_TXIMPORT.out.txi_dtu, TX2GENE_TXIMPORT.out.tximport_tx2gene,
        VALIDATE_INPUTS.out.samples, VALIDATE_INPUTS.out.contrasts, params.n_dexseq_plot,
        params.min_samps_gene_expr, params.min_samps_feature_expr, params.min_samps_feature_prop,
        params.min_feature_expr, params.min_feature_prop, params.min_gene_expr)
    ch_versions = ch_versions.mix(DEXSEQ_ANNOTATION.out.versions, DEXSEQ_DEU.out.versions,
        TX2GENE_TXIMPORT.out.versions, DRIMSEQ_DEXSEQ_DTU.out.versions)

    // Asynchronous task completion must never determine rMATS vector ordering.
    ch_genome_conditions = SAMTOOLS_SORT.out.bam.map { meta, bam -> tuple(meta.condition, [meta, bam]) }
        .groupTuple().map { condition, pairs ->
            def ordered = pairs.sort { a, b -> a[0].order_index <=> b[0].order_index }
            tuple(condition, ordered.collect { it[0] }, ordered.collect { it[1] })
        }
    RMATS(VALIDATE_INPUTS.out.contrasts, ch_genome_conditions, FILTER_ANALYSIS_GTF.out.gtf,
        false, ch_rmats_read_length, params.rmats_splice_diff_cutoff,
        params.rmats_novel_splice_site, params.rmats_min_intron_len, params.rmats_max_exon_len, false)
    ch_versions = ch_versions.mix(RMATS.out.versions)

    ADAPT_SUPPA_CONTRASTS(VALIDATE_INPUTS.out.contrasts)
    SUPPA(FILTER_ANALYSIS_GTF.out.gtf, TX2GENE_TXIMPORT.out.suppa_tpm,
        VALIDATE_INPUTS.out.samples, ADAPT_SUPPA_CONTRASTS.out.contrasts,
        params.suppa_per_local_event, params.generateevents_boundary, params.generateevents_threshold,
        params.generateevents_exon_length, params.generateevents_event_type, params.generateevents_pool_genes,
        params.psiperevent_total_filter, params.diffsplice_local_event, params.diffsplice_method,
        params.diffsplice_area, params.diffsplice_lower_bound, params.diffsplice_alpha,
        params.diffsplice_tpm_threshold, params.diffsplice_nan_threshold, params.diffsplice_gene_correction,
        false, params.diffsplice_median, params.clusterevents_local_event, params.clusterevents_dpsithreshold,
        params.clusterevents_eps, params.clusterevents_metric, params.clusterevents_min_pts,
        params.clusterevents_method, params.clusterevents_sigthreshold ?: false,
        params.clusterevents_separation ?: false, params.suppa_per_isoform)
    ch_versions = ch_versions.mix(SUPPA.out.versions)

    KATANA_DIFFERENTIAL([id: 'study'], SALMON_GENOME_PROJECTION.out.results.map { it[1] }.collect(),
        VALIDATE_INPUTS.out.samples, VALIDATE_INPUTS.out.contrasts, FILTER_ANALYSIS_GTF.out.gtf, params.design, '')
    ch_versions = ch_versions.mix(KATANA_DIFFERENTIAL.out.versions)

    // Every artifact has explicit provenance and destination. The export task
    // stages channel outputs directly; it never races asynchronous publishDir.
    ch_artifacts = PREPARE_REFERENCES.out.fasta.map { artifact('genome_fasta', it, 'refs/genome.fasta') }
        .mix(Channel.value(file("${projectDir}/assets/software-provenance.lock.json", checkIfExists: true)).map { artifact('provenance', it, 'provenance/software.lock.json') })
        .mix(Channel.value(file("${projectDir}/assets/migration-source-manifest.json", checkIfExists: true)).map { artifact('provenance', it, 'provenance/implementation-migration.json') })
        .mix(FILTER_ANALYSIS_GTF.out.gtf.map { artifact('gtf', it, 'refs/annotation.gtf') })
        .mix(PREPARE_REFERENCES.out.gtf.map { artifact('provenance', it, 'provenance/reference/annotation.unfiltered.gtf') })
        .mix(ANALYSIS_TRANSCRIPTS.out.transcripts.map { artifact('transcript_fasta', it, 'refs/transcripts.fasta') })
        .mix(FILTER_ANALYSIS_GTF.out.decisions.map { artifact('provenance', it, 'provenance/gene_exclusions/annotation_decisions.tsv') })
        .mix(FILTER_ANALYSIS_GTF.out.receipt.map { artifact('provenance', it, 'provenance/gene_exclusions/filtering.receipt.json') })
        .mix(PREPARE_REFERENCES.out.hashes.map { artifact('provenance', it, 'provenance/reference/unfiltered_reference.sha256') })
        .mix(VALIDATE_INPUTS.out.audit.map { artifact('provenance', it, 'pipeline_info/input_audit.json') })
        .mix(VALIDATE_REFERENCE.out.audit.map { artifact('provenance', it, 'provenance/reference_audit.json') })
        .mix(STAGE_FASTQ.out.audit.map { meta, p -> artifact('provenance', p, "provenance/fastq_staging/${meta.id}.json", [sample_id: meta.id]) })
        .mix(ADAPT_SUPPA_CONTRASTS.out.contrasts.map { artifact('provenance', it, 'provenance/contrastsheet.suppa.csv') })
        .mix(AUDIT_GENOME_BAM.out.audit.map { meta, p -> artifact('provenance', p, "provenance/bam_audit/${meta.id}.json", [sample_id: meta.id]) })
        .mix(CAPTURE_BAM_AUDIT.out.captured.flatMap { meta, header, metrics -> [
            artifact('provenance', header, "provenance/bam_audit/${meta.id}.header.sam", [sample_id: meta.id]),
            artifact('provenance', metrics, "provenance/bam_audit/${meta.id}.primary_metrics.tsv", [sample_id: meta.id])] })
        .mix(REFERENCE_FAIDX.out.fai.map { artifact('auxiliary', it, 'refs/genome.fasta.fai') })
        .mix(REFERENCE_FAIDX.out.dict.map { artifact('provenance', it, 'provenance/reference.genome.dict') })
        .mix(SALMON_GENOME_PROJECTION.out.results.map { meta, directory -> artifact('quant', directory.resolve('quant.sf'), "quant/${meta.id}/quant.sf", [sample_id: meta.id]) })
        .mix(SAMTOOLS_SORT.out.bam.map { meta, bam -> artifact('bam', bam, "bam/${meta.id}.bam", [sample_id: meta.id]) })
        .mix(SAMTOOLS_INDEX.out.bai.map { meta, bai -> artifact('bai', bai, "bam/${meta.id}.bam.bai", [sample_id: meta.id]) })
        .mix(KATANA_DIFFERENTIAL.out.results.flatMap { meta, files -> tableArtifacts(files, 'deg', 'DEG', /^(.+)\.deseq2\.results\.tsv$/) })
        .mix(KATANA_DIFFERENTIAL.out.normalised_counts.map { meta, p -> artifact('normalized_counts', p, 'DEG/all.normalised_counts.tsv') })
        .mix(KATANA_DIFFERENTIAL.out.annotation.map { meta, p -> artifact('annotation', p, 'refs/annotation.tsv') })
        .mix(KATANA_DIFFERENTIAL.out.tx2gene.map { meta, p -> artifact('matrix', p, 'refs/tx2gene.tsv') })
        .mix(KATANA_DIFFERENTIAL.out.run_metadata.map { meta, p -> artifact('provenance', p, "provenance/${p.name}") })
        .mix(KATANA_DIFFERENTIAL.out.session_info.map { meta, p -> artifact('provenance', p, "provenance/${p.name}") })
        .mix(KATANA_DIFFERENTIAL.out.coverage.map { meta, p -> artifact('provenance', p, "provenance/${p.name}") })
        .mix(KATANA_DIFFERENTIAL.out.gene_matrices.flatMap { meta, files -> [files].flatten().collect { p -> artifact('matrix', p, "quant/matrices/${p.name}") } })
        .mix(KATANA_DIFFERENTIAL.out.normalization_factors.map { meta, p -> artifact('matrix', p, "DEG/${p.name}") })
        .mix(KATANA_DIFFERENTIAL.out.models.map { meta, p -> artifact('provenance', p, "provenance/${p.name}") })
        .mix(DEXSEQ_ANNOTATION.out.gff.map { artifact('auxiliary', it, 'dexseq_exon/annotation/DEXSeq.gff') })
        .mix(DEXSEQ_DEU.out.dexseq_clean_txt.flatten().map { artifact('matrix', it, "dexseq_exon/counts/${it.name}") })
        .mix(DEXSEQ_DEU.out.dexseq_gene_results_csv.flatMap { files -> tableArtifacts(files, 'dexseq', 'dexseq_exon', /^perGeneQValue\.(.+)\.csv$/, [direction: 'not_applicable']) })
        .mix(DEXSEQ_DEU.out.dexseq_exon_results_csv.flatMap { files -> tableArtifacts(files, 'auxiliary', 'dexseq_exon', /^DEXSeqResults\.(.+)\.csv$/, [native_fold_change_direction: 'control_over_treatment_native_header_authoritative']) })
        .mix(TX2GENE_TXIMPORT.out.tpm_transcript.map { artifact('matrix', it, "quant/matrices/${it.name}") })
        .mix(TX2GENE_TXIMPORT.out.counts_transcript.map { artifact('matrix', it, "quant/matrices/${it.name}") })
        .mix(TX2GENE_TXIMPORT.out.tpm_gene.map { artifact('matrix', it, "quant/matrices/${it.name}") })
        .mix(TX2GENE_TXIMPORT.out.counts_gene.map { artifact('matrix', it, "quant/matrices/${it.name}") })
        .mix(SUPPA.out.suppa_local_psi.map { artifact('matrix', it, "suppa/matrices/${it.name}") })
        .mix(SUPPA.out.suppa_isoform_psi.map { artifact('matrix', it, "suppa/matrices/${it.name}") })
        .mix(ch_qc_artifacts)
        .mix(ch_versions.flatten().filter { it != null }.unique().map { p -> artifact('provenance', p, "provenance/versions/${p.parent.name}.versions.yml") })

    // Native DTU scripts use treatment-control filenames. Translate their names
    // to the declared ID without altering any statistical values or directions.
    ch_dtu_artifacts = DRIMSEQ_DEXSEQ_DTU.out.dexseq_exon_results_tsv.map { [files: it] }
        .combine(ch_contrasts_container).flatMap { result, cmp ->
            nativeDtuArtifacts(result.files, 'dtu', 'dtu/results/dexseq', /^(DEXSeqResults)\.(.+)\.tsv$/, cmp.records)
        }
        .mix(DRIMSEQ_DEXSEQ_DTU.out.dexseq_gene_results_tsv.map { [files: it] }
            .combine(ch_contrasts_container).flatMap { result, cmp ->
                nativeDtuArtifacts(result.files, 'auxiliary', 'dtu/results/dexseq', /^(perGeneQValue)\.(.+)\.tsv$/, cmp.records)
            })
        .mix(DRIMSEQ_DEXSEQ_DTU.out.stager_padj_tsv.map { [files: it] }
            .combine(ch_contrasts_container).flatMap { result, cmp ->
                nativeDtuArtifacts(result.files, 'stager', 'dtu/results/stager', /^(getAdjustedPValues)\.(.+)\.tsv$/, cmp.records)
            })

    ch_rmats_artifacts = RMATS.out.rmats_post.map { [files: it] }.combine(ch_contrasts_container).flatMap { file_container, contrast_container ->
        // The native glob also emits its tmp directory; bundle artifacts are
        // regular result files. Missing required event/mode tables remain a
        // fatal error in the full exporter.
        [file_container.files].flatten().findAll { java.nio.file.Files.isRegularFile(it) }.collect { p ->
            def pair_name = p.parent.parent.name
            def cmp = contrast_container.records.find { "${it.treatment}-${it.control}" == pair_name }
            if (!cmp) error("rMATS directory has no declared contrast: ${pair_name}")
            def match = p.name =~ /^(SE|RI|A3SS|A5SS|MXE)\.MATS\.(JC|JCEC)\.txt$/
            def spec = [comparison_id: cmp.contrast]
            if (match.matches()) spec += [event_type: match[0][1], mode: match[0][2]]
            artifact(match.matches() ? 'rmats' : 'auxiliary', p, "rmats/${cmp.contrast}/rmats_post/${p.name}", spec)
        }
    }
    ch_suppa_artifacts = SUPPA.out.dpsi_local.map { treatment, control, files -> [treatment: treatment, control: control, files: files] }
        .combine(ch_contrasts_container).flatMap { result, contrast_container ->
        def cmp = contrast_container.records.find { it.control == result.treatment && it.treatment == result.control }
        if (!cmp) error('SUPPA output has no declared contrast')
        [result.files].flatten().collect { p -> artifact('suppa', p, "suppa/${cmp.contrast}/diffsplice/per_local_event/${p.name}", [comparison_id: cmp.contrast,
            direction: 'treatment_minus_control', source_input1: cmp.control, source_input2: cmp.treatment]) }
    }
    ch_suppa_isoform_artifacts = SUPPA.out.dpsi_isoform.map { treatment, control, files -> [treatment: treatment, control: control, files: files] }
        .combine(ch_contrasts_container).flatMap { result, contrast_container ->
        def cmp = contrast_container.records.find { it.control == result.treatment && it.treatment == result.control }
        if (!cmp) error('SUPPA isoform output has no declared contrast')
        [result.files].flatten().collect { p -> artifact('auxiliary', p, "suppa/${cmp.contrast}/diffsplice/per_isoform/${p.name}", [comparison_id: cmp.contrast,
            direction: 'treatment_minus_control', source_input1: cmp.control, source_input2: cmp.treatment]) }
    }
    ch_export_inputs = ch_artifacts.mix(ch_dtu_artifacts, ch_rmats_artifacts, ch_suppa_artifacts, ch_suppa_isoform_artifacts).collect(flat: false)
        .map { artifacts -> tuple(artifacts.collect { it[0] }, artifacts.collect { it[1] }) }
    ch_study_metadata = ch_rows_container.combine(ch_contrasts_container).map { row_container, contrast_container ->
        def rows = row_container.records
        def contrasts = contrast_container.records
        [analysis_purpose: params.analysis_purpose, design: params.design,
            assay: 'bulk_rna_seq', assay_validation: 'user_declared_or_defaulted_not_inferred_from_files',
            reference_declarations: [genome_build: params.reference_build ?: params.genome ?: null,
                annotation_release: params.annotation_release ?: null,
                validation: 'user_declared_labels_not_independently_verified'],
            quantification_method: 'salmon_2.8.0_genome_projection',
            gene_exclusions: [exclude_pseudogenes: true, exclude_mitochondrial: true,
                application_stage: 'before_quantification_counts_and_differential_models',
                receipt: 'provenance/gene_exclusions/filtering.receipt.json',
                raw_annotation: 'provenance/reference/annotation.unfiltered.gtf',
                alignment_annotation: 'unfiltered_reference_preserved',
                already_completed_projects_refitted: false],
            rmats_read_length_semantics: 'declared nominal length for effective event lengths; upstream variable-read-length and clipping enabled',
            allow_variable_read_length: params.allow_variable_read_length,
            rmats_effective_settings: [version: '4.1.2', cstat: params.rmats_splice_diff_cutoff,
                novelSS: params.rmats_novel_splice_site,
                mil: params.rmats_novel_splice_site ? params.rmats_min_intron_len : null,
                mel: params.rmats_novel_splice_site ? params.rmats_max_exon_len : null,
                readLength: rows[0].read_length as int, variable_read_length: true,
                allow_clipping: true, paired_stats: false],
            dexseq_native_fold_change_direction: 'control_over_treatment_native_header_authoritative; gene/stageR q-values have no direction',
            rmats_sample_order: rows.groupBy { it.condition }.collectEntries { condition, values -> [(condition): values.collect { it.sample }] },
            comparisons: contrasts.collect { [comparison_id: it.contrast, treatment: it.treatment, control: it.control,
                rmats_group1: it.treatment, rmats_group2: it.control, deg_target: it.treatment, deg_reference: it.control,
                suppa_input1: it.control, suppa_input2: it.treatment, suppa_direction: 'treatment_minus_control'] }]
    }
    KATANA_EXPORT(ch_export_inputs, VALIDATE_INPUTS.out.samples, VALIDATE_INPUTS.out.contrasts,
        Channel.value(file("${projectDir}/bin/katana_bundle.py", checkIfExists: true)), ch_study_metadata)

    emit:
    bundle = KATANA_EXPORT.out.bundle
    quant = SALMON_GENOME_PROJECTION.out.results
    bam = SAMTOOLS_SORT.out.bam
    bai = SAMTOOLS_INDEX.out.bai
    samplesheet = VALIDATE_INPUTS.out.samples
    contrastsheet = VALIDATE_INPUTS.out.contrasts
    versions = ch_versions
    qc_files = ch_qc_files.flatten()
}
