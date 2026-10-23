process CAPTURE_BAM_AUDIT {
    tag "${meta.id}"
    label 'process_medium'
    container "${workflow.containerEngine == 'singularity' ? 'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' : 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'}"
    input:
    tuple val(meta), path(bam, stageAs: 'katana_source.bam')
    path metrics_helper
    output:
    tuple val(meta), path('audit.header.sam'), path('audit.primary_metrics.tsv'), emit: captured
    tuple val(meta), path('katana_source.bam', includeInputs: true), emit: staged_bam
    path 'versions.yml', emit: versions
    script:
    """
    samtools quickcheck -v "$bam"
    samtools view -H "$bam" > audit.header.sam
    samtools view -@ ${task.cpus} -F 2304 "$bam" | awk -f "$metrics_helper" > audit.primary_metrics.tsv
    printf 'samtools: %s\\n' "\$(samtools --version | head -n 1)" > versions.yml
    """
}

