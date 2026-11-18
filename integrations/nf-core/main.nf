#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { KATANA_BULK_AS      } from './subworkflows/local/katana_bulk_as/main'
include { KATANA_SMARTSEQ3_AS } from './subworkflows/local/katana_smartseq3_as/main'

workflow {
    ch_fasta = Channel.value(file(params.fasta, checkIfExists: true))
    ch_fai   = Channel.value(file(params.fai ?: "${params.fasta}.fai", checkIfExists: true))
    ch_gtf   = params.gtf ? Channel.value(file(params.gtf, checkIfExists: true)) : Channel.value([])
    ch_gmt   = params.gmt ? Channel.value(file(params.gmt, checkIfExists: true)) : Channel.value([])

    if (params.mode == 'bulk') {
        ch_study = Channel.of([
            [id: 'study'],
            file(params.usage, checkIfExists: true),
            file(params.samples, checkIfExists: true),
            params.support ? file(params.support, checkIfExists: true) : [],
            params.annotation ? file(params.annotation, checkIfExists: true) : []
        ])
        ch_rmats = Channel.of([[id: 'contrast'], file(params.rmats, checkIfExists: true)])
        KATANA_BULK_AS(ch_study, ch_gmt, ch_rmats, ch_fasta, ch_fai, ch_gtf)
    } else {
        ch_sc_net = Channel.of([
            [id: 'sc_cohort'],
            file(params.psi, checkIfExists: true),
            file(params.metadata, checkIfExists: true),
            file(params.events, checkIfExists: true),
            params.support ? file(params.support, checkIfExists: true) : []
        ])
        ch_sc_trace = Channel.of([
            [id: 'sc_contrast'],
            file(params.events, checkIfExists: true),
            file(params.psi, checkIfExists: true),
            file(params.metadata, checkIfExists: true),
            params.support ? file(params.support, checkIfExists: true) : []
        ])
        KATANA_SMARTSEQ3_AS(ch_sc_net, ch_gmt, ch_sc_trace, ch_fasta, ch_fai, ch_gtf, params.numerator, params.denominator)
    }
}
