include { FILTER_ANALYSIS_GTF } from '../modules/local/analysis_reference'

workflow {
    FILTER_ANALYSIS_GTF(Channel.value(file(params.fixture_gtf, checkIfExists: true)),
        Channel.value(file(params.filter_script, checkIfExists: true)))
    FILTER_ANALYSIS_GTF.out.gtf.view { "filtered_gtf=${it}" }
    FILTER_ANALYSIS_GTF.out.receipt.view { "filter_receipt=${it}" }
}
