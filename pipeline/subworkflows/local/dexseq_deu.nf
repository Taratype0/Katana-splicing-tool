//
// DEXSeq DEU subworkflow
//

include { DEXSEQ_ANNOTATION   } from '../../modules/local/dexseq_annotation'
include { DEXSEQ_COUNT        } from '../../modules/local/dexseq_count'
include { DEXSEQ_EXON         } from '../../modules/local/dexseq_exon'
include { CONTRAST_TASK_SHEETS } from '../../modules/local/split_contrasts'

def stableCountVersion(counts, versions) {
    // A cached task can arrive before a task that originally finished first.
    // Select one declared version output by sample identity, never arrival or
    // its work-directory hash. Keep counts and versions paired by producer.
    def sample_directories = counts.map { meta, files ->
        def paths = [files].flatten()
        if (paths.size() != 1) error('DEXSeq count version selection requires one count file per sample')
        tuple(paths[0].parent.toString(), meta.id.toString())
    }
    def version_directories = versions.map { path -> tuple(path.parent.toString(), path) }
    sample_directories.join(version_directories, failOnDuplicate: true, failOnMismatch: true)
        .toSortedList { a, b -> a[1] <=> b[1] }
        .map { rows ->
            if (!rows || rows.collect { it[1] }.unique(false).size() != rows.size()) {
                error('DEXSeq count version selection requires unique sample identities')
            }
            rows[0][2]
        }
}

workflow DEXSEQ_DEU {

    take:

    gtf                // path gtf
    ch_genome_bam      // bam channel
    ch_dexseq_gff      // path dexseq gff
    ch_samplesheet     // Channel.fromPath(params.input)
    ch_contrastsheet   // Channel.fromPath()
    n_dexseq_plot      // val: numeric
    aggregation        // params.aggregation
    alignment_quality  // params.alignment_quality

    main:

    ch_versions = Channel.empty()

    if (!ch_dexseq_gff) {

        //
        // MODULE: DEXSeq Annotation
        //

        DEXSEQ_ANNOTATION (
            gtf,
            aggregation
        )

        ch_versions = ch_versions.mix(DEXSEQ_ANNOTATION.out.versions.first())

        ch_dexseq_gff = DEXSEQ_ANNOTATION.out.gff

    }

    ch_genome_bam_gff = ch_genome_bam.combine(ch_dexseq_gff)

    //
    // MODULE: DEXSeq Count
    //

    DEXSEQ_COUNT (
        ch_genome_bam_gff,
        alignment_quality
    )

    ch_versions = ch_versions.mix(stableCountVersion(DEXSEQ_COUNT.out.dexseq_clean_txt, DEXSEQ_COUNT.out.versions))

    //
    // MODULE: DEXSeq DEU
    //

    CONTRAST_TASK_SHEETS(ch_contrastsheet.first())
    DEXSEQ_EXON (
        DEXSEQ_COUNT.out.dexseq_clean_txt.map{ it[1] }.collect(),
        ch_dexseq_gff.first(),
        ch_samplesheet.first(),
        CONTRAST_TASK_SHEETS.out.sheets,
        n_dexseq_plot
    )

    ch_versions = ch_versions.mix(CONTRAST_TASK_SHEETS.out.versions, DEXSEQ_EXON.out.versions)

    emit:

    dexseq_clean_txt        = DEXSEQ_COUNT.out.dexseq_clean_txt.map{ it[1] }.collect()

    dexseq_exon_dataset_rds = DEXSEQ_EXON.out.dexseq_exon_dataset_rds
    dexseq_exon_results_rds = DEXSEQ_EXON.out.dexseq_exon_results_rds
    dexseq_gene_results_rds = DEXSEQ_EXON.out.dexseq_gene_results_rds
    dexseq_exon_results_csv = DEXSEQ_EXON.out.dexseq_exon_results_csv
    dexseq_gene_results_csv = DEXSEQ_EXON.out.dexseq_gene_results_csv
    dexseq_plot_results_pdf = DEXSEQ_EXON.out.dexseq_plot_results_pdf

    versions                = ch_versions                                  // channel: [ versions.yml ]

}
