include { GFFREAD_STRAND_MODELS; RSEQC_INFEREXPERIMENT } from '../../modules/local/infer_strandedness'
include { RESOLVE_SAMPLE_STRANDEDNESS } from '../../modules/local/resolve_sample_strandedness'

workflow RESOLVE_STRANDEDNESS {
    take:
    ch_bam
    ch_bai
    ch_samplesheet
    ch_contrastsheet
    ch_gtf
    ch_resolver
    ch_validator

    main:
    ch_indexed = ch_bam.map { meta, bam -> tuple(meta.id, meta, bam) }
        .join(ch_bai.map { meta, bai -> tuple(meta.id, bai) }, failOnDuplicate: true, failOnMismatch: true)
        .map { id, meta, bam, bai -> tuple(meta, bam, bai) }
    ch_auto = ch_indexed.filter { meta, bam, bai -> meta.strandedness == 'auto' }
    // Build a matched model only when at least one library requests inference.
    ch_auto_model = ch_samplesheet.splitCsv(header: true).filter { it.strandedness == 'auto' }
        .map { true }.first().combine(ch_gtf).map { needed, gtf -> gtf }
    GFFREAD_STRAND_MODELS(ch_auto_model)
    RSEQC_INFEREXPERIMENT(ch_auto, GFFREAD_STRAND_MODELS.out.bed.first(),
        params.strandedness_sample_size, params.strandedness_mapq)
    ch_evidence = RSEQC_INFEREXPERIMENT.out.txt.map { meta, report -> tuple(meta.id, report) }
        .join(RSEQC_INFEREXPERIMENT.out.log.map { meta, log -> tuple(meta.id, log) }, failOnDuplicate: true, failOnMismatch: true)
        .toList().map { records ->
            def ordered = records.sort(false) { a, b -> a[0] <=> b[0] }
            tuple(ordered.collect { it[0] }, ordered.collect { it[1] }, ordered.collect { it[2] })
        }
    RESOLVE_SAMPLE_STRANDEDNESS(ch_samplesheet, ch_contrastsheet, ch_evidence, ch_resolver, ch_validator,
        params.stranded_threshold, params.unstranded_threshold, params.design, params.analysis_purpose)
    ch_resolved_meta = RESOLVE_SAMPLE_STRANDEDNESS.out.samplesheet.splitCsv(header: true)
        .map { row -> tuple(row.sample, row.strandedness) }
    // The whole-study resolver validates actual libraries within every contrast
    // before any BAM can reach Salmon, DEXSeq counting or rMATS preparation.
    ch_resolved_bam = ch_bam.map { meta, bam -> tuple(meta.id, meta, bam) }
        .join(ch_resolved_meta, failOnDuplicate: true, failOnMismatch: true)
        .map { id, meta, bam, strand ->
            if (!(strand in ['forward', 'reverse', 'unstranded'])) error("Unresolved strand for ${id}")
            tuple(meta + [strandedness: strand], bam)
        }
    ch_resolved_bai = ch_indexed.map { meta, bam, bai -> tuple(meta.id, meta, bai) }
        .join(ch_resolved_meta, failOnDuplicate: true, failOnMismatch: true)
        .map { id, meta, bai, strand -> tuple(meta + [strandedness: strand], bai) }

    emit:
    bam = ch_resolved_bam
    bai = ch_resolved_bai
    samplesheet = RESOLVE_SAMPLE_STRANDEDNESS.out.samplesheet
    audit = RESOLVE_SAMPLE_STRANDEDNESS.out.audit
    summary = RESOLVE_SAMPLE_STRANDEDNESS.out.summary
    native_reports = RSEQC_INFEREXPERIMENT.out.txt
    native_logs = RSEQC_INFEREXPERIMENT.out.log
    bed = GFFREAD_STRAND_MODELS.out.bed
    versions = GFFREAD_STRAND_MODELS.out.versions.mix(RSEQC_INFEREXPERIMENT.out.versions)
}
