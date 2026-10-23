process STAGE_FASTQ {
    tag "${meta.id}"
    label 'process_low'
    container 'docker.io/library/python:3.12.11-bookworm'

    input:
    tuple val(meta), path(reads, stageAs: 'source_reads/??????/*')
    path helper, stageAs: 'kanata_stage_fastq.py'

    output:
    tuple val(meta), path("${meta.id}_R*.fastq.gz"), emit: reads
    tuple val(meta), path('fastq_staging_audit.json'), emit: audit

    script:
    if (!(meta.id.toString() ==~ /^[A-Za-z][A-Za-z0-9_.-]*$/)) {
        error('Invalid sample ID for FASTQ staging')
    }
    def ordered_reads = reads instanceof List ? reads : [reads]
    def request = groovy.json.JsonOutput.toJson([reads: ordered_reads.collect { it.toString() }, single_end: meta.single_end])
    """
    cat > reads.json <<'KANATA_FASTQ_REQUEST'
    ${request}
    KANATA_FASTQ_REQUEST
    python3 kanata_stage_fastq.py --reads-json reads.json --sample "${meta.id}" --output-dir .
    """
}
