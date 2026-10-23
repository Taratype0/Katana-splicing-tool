process VALIDATE_INPUTS {
    label 'process_low'
    container 'docker.io/library/python:3.12.11-bookworm'
    input:
    path samplesheet
    path contrastsheet
    path raw_files, stageAs: 'validation_inputs/??????/*'
    path fasta
    path gtf
    path validator
    val samplesheet_base
    val read_length
    val design
    val analysis_purpose
    output:
    path 'validated/samplesheet.valid.csv', emit: samples
    path 'validated/contrastsheet.valid.csv', emit: contrasts
    path 'validated/input_audit.json', emit: audit
    script:
    def shellQuote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    def read_args = read_length ? "--read-length ${read_length}" : ''
    """
    python3 ${shellQuote(validator)} --samplesheet ${shellQuote(samplesheet)} --samplesheet-base ${shellQuote(samplesheet_base)} \\
      --contrasts ${shellQuote(contrastsheet)} --outdir validated --fasta ${shellQuote(fasta)} --gtf ${shellQuote(gtf)} \\
      --design ${shellQuote(design)} --analysis-purpose ${shellQuote(analysis_purpose)} $read_args
    """
}

