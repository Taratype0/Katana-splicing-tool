include { KATANA_SCSPLICENET   } from '../../../modules/local/katana_scsplicenet/main'
include { KATANA_SCSPLICETRACE } from '../../../modules/local/katana_scsplicetrace/main'

workflow KATANA_SMARTSEQ3_AS {
    take:
    ch_sc_net_inputs   // channel: [ val(meta), path(usage), path(metadata), path(event_genes), path(support) ]
    ch_gmt             // channel: path(gmt) or []
    ch_sc_trace_inputs // channel: [ val(meta), path(events), path(psi), path(metadata), path(support) ]
    ch_fasta           // channel: path(fasta)
    ch_fai             // channel: path(fai)
    ch_gtf             // channel: path(gtf) or []
    val_numerator      // val: string
    val_denominator    // val: string

    main:
    ch_versions = Channel.empty()

    // 1. Single-Cell ExcessBioVar + Fisher-z Shrunk Signed TOM Co-Splicing Modules
    KATANA_SCSPLICENET(ch_sc_net_inputs, ch_gmt)
    ch_versions = ch_versions.mix(KATANA_SCSPLICENET.out.versions.first())

    // 2. Single-Cell RBP Discovery + Mature mRNA + GTF/UniProt Proteome Consequences + Module Bridge
    KATANA_SCSPLICETRACE(
        ch_sc_trace_inputs,
        ch_fasta,
        ch_fai,
        ch_gtf,
        KATANA_SCSPLICENET.out.module_table.first(),
        val_numerator,
        val_denominator
    )
    ch_versions = ch_versions.mix(KATANA_SCSPLICETRACE.out.versions.first())

    emit:
    sc_splicenet_results   = KATANA_SCSPLICENET.out.results
    sc_module_table        = KATANA_SCSPLICENET.out.module_table
    sc_splicetrace_results = KATANA_SCSPLICETRACE.out.results
    versions               = ch_versions
}
