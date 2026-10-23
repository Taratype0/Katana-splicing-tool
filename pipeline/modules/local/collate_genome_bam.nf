process COLLATE_GENOME_BAM {
    tag "${meta.id}"
    label 'process_medium'
    container "${workflow.containerEngine == 'singularity' ? 'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' : 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'}"
    input:
    tuple val(meta), path(bam)
    output:
    tuple val(meta), path("${meta.id}.name_grouped.bam"), emit: bam
    path 'versions.yml', emit: versions
    script:
    """
    samtools quickcheck -v "$bam"
    samtools collate -@ ${task.cpus} -o ${meta.id}.name_grouped.bam "$bam"
    samtools quickcheck -v ${meta.id}.name_grouped.bam
    printf 'samtools: %s\\n' "\$(samtools --version | head -n 1)" > versions.yml
    """
}

