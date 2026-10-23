process PREPARE_REFERENCES {
    label 'process_low'
    container "${workflow.containerEngine == 'singularity' ? 'https://depot.galaxyproject.org/singularity/gffread:0.12.1--h8b12597_0' : 'quay.io/biocontainers/gffread:0.12.1--h8b12597_0'}"
    input:
    path fasta
    path gtf
    output:
    path 'genome.fasta', emit: fasta
    path 'annotation.gtf', emit: gtf
    path 'transcripts.fasta', emit: transcripts
    path 'reference.sha256', emit: hashes
    path 'versions.yml', emit: versions
    script:
    def shellQuote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    def fasta_cmd = fasta.name.endsWith('.gz') ? "gzip -cd -- ${shellQuote(fasta)}" : "cat -- ${shellQuote(fasta)}"
    def gtf_cmd = gtf.name.endsWith('.gz') ? "gzip -cd -- ${shellQuote(gtf)}" : "cat -- ${shellQuote(gtf)}"
    """
    $fasta_cmd > genome.fasta
    $gtf_cmd > annotation.gtf
    gffread annotation.gtf -g genome.fasta -w transcripts.fasta
    test -s genome.fasta && test -s annotation.gtf && test -s transcripts.fasta
    sha256sum genome.fasta annotation.gtf transcripts.fasta > reference.sha256
    printf 'gffread: %s\\n' "\$(gffread --version 2>&1)" > versions.yml
    """
}

