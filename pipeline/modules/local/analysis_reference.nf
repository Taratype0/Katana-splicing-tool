process FILTER_ANALYSIS_GTF {
    label 'process_single'
    container 'docker.io/library/python:3.12.11-bookworm'
    input:
    path gtf
    path filter_script
    output:
    path 'analysis_reference/annotation.analysis.gtf', emit: gtf
    path 'analysis_reference/annotation_decisions.tsv', emit: decisions
    path 'analysis_reference/filtering.receipt.json', emit: receipt
    script:
    def quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    """
    python3 ${quote(filter_script)} --gtf ${quote(gtf)} --outdir analysis_reference
    """
}

process ANALYSIS_TRANSCRIPTS {
    label 'process_low'
    container "${workflow.containerEngine == 'singularity' ? 'https://depot.galaxyproject.org/singularity/gffread:0.12.1--h8b12597_0' : 'biocontainers/gffread:0.12.1--h8b12597_0'}"
    input:
    path fasta
    path gtf
    output:
    path 'transcripts.analysis.fasta', emit: transcripts
    script:
    def quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    """
    gffread ${quote(gtf)} -g ${quote(fasta)} -w transcripts.analysis.fasta
    test -s transcripts.analysis.fasta
    """
}
