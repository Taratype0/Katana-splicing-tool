process REFERENCE_FAIDX {
    label 'process_low'
    container "${workflow.containerEngine == 'singularity' ? 'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' : 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'}"
    input:
    path fasta
    output:
    path '*.fai', emit: fai
    path 'genome.dict', emit: dict
    script:
    """
    samtools faidx "$fasta"
    samtools dict "$fasta" > genome.dict
    """
}

