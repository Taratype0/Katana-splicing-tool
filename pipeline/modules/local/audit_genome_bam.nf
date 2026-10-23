process AUDIT_GENOME_BAM {
    tag "${meta.id}"
    label 'process_low'
    container 'docker.io/library/python:3.12.11-bookworm'
    input:
    tuple val(meta), path(header), path(metrics)
    path fai
    path reference_dict
    path helper
    val allow_variable_read_length
    output:
    tuple val(meta), path('bam_audit.json'), emit: audit
    script:
    def metadata = groovy.json.JsonOutput.toJson(meta)
    def variable_args = allow_variable_read_length ? '--allow-variable-read-length' : ''
    """
    cat > meta.json <<'KATANA_BAM_META'
    ${metadata}
    KATANA_BAM_META
    python3 "$helper" --header "$header" --primary-metrics "$metrics" --fai "$fai" --reference-dict "$reference_dict" \\
      --meta-json meta.json --output bam_audit.json $variable_args
    """
}

