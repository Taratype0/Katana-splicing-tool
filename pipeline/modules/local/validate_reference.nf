process VALIDATE_REFERENCE {
    label 'process_low'
    container 'docker.io/library/python:3.12.11-bookworm'

    input:
    path metadata_audit
    path fasta
    path gtf
    path validator
    val reference_build
    val annotation_release

    output:
    path 'reference_audit.json', emit: audit

    script:
    def shellQuote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    def build_arg = reference_build ? "--reference-build ${shellQuote(reference_build)}" : ''
    def release_arg = annotation_release ? "--annotation-release ${shellQuote(annotation_release)}" : ''
    """
    test -s ${shellQuote(metadata_audit)}
    python3 ${shellQuote(validator)} --fasta ${shellQuote(fasta)} --gtf ${shellQuote(gtf)} \\
      --output reference_audit.json $build_arg $release_arg
    """
}
