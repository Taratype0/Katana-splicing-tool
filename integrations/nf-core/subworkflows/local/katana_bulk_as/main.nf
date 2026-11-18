include { KATANA_SPLICENET   } from '../../../modules/local/katana_splicenet/main'
include { KATANA_SPLICETRACE } from '../../../modules/local/katana_splicetrace/main'

workflow KATANA_BULK_AS {
    take:
    ch_study_usage     // channel: [ val(meta), path(usage), path(samples), path(support), path(annotation) ]
    ch_gmt             // channel: path(gmt) or []
    ch_rmats_contrasts // channel: [ val(meta), path(rmats_dir_or_tar_or_tsv) ]
    ch_fasta           // channel: path(fasta)
    ch_fai             // channel: path(fai)
    ch_gtf             // channel: path(gtf) or []

    main:
    ch_versions = Channel.empty()

    // 1. Study-wide ExcessBioVar + Adaptive Signed TOM Co-Splicing Modules
    KATANA_SPLICENET(ch_study_usage, ch_gmt)
    ch_versions = ch_versions.mix(KATANA_SPLICENET.out.versions.first())

    // 2. Per-contrast RBP Discovery + Mature mRNA + GTF/UniProt Proteome Consequences + Cross-Module Bridge
    KATANA_SPLICETRACE(
        ch_rmats_contrasts,
        ch_fasta,
        ch_fai,
        ch_gtf,
        KATANA_SPLICENET.out.module_table.first()
    )
    ch_versions = ch_versions.mix(KATANA_SPLICETRACE.out.versions.first())

    emit:
    splicenet_results   = KATANA_SPLICENET.out.results
    module_table        = KATANA_SPLICENET.out.module_table
    splicetrace_results = KATANA_SPLICETRACE.out.results
    versions            = ch_versions
}
