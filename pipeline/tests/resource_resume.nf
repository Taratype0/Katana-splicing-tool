nextflow.enable.dsl=2

process RESOURCE_CHECK {
    label 'katana_downstream'
    tag "${spec.id}"
    input:
    val spec
    output:
    path 'complete.txt'
    script:
    """
    if [ '${spec.id}' = 'retry' ] && [ ! -f '${params.fixture_marker}' ]; then
        touch '${params.fixture_marker}'
        exit 137
    fi
    echo '${spec.id}' > complete.txt
    """
}

workflow {
    RESOURCE_CHECK(Channel.of(
        [id:'plain', execution_group:'network_fit', resources:[cpus:1, memory:'16 MB', time:'1m']],
        [id:'retry', execution_group:'motiftrace', resources:[cpus:1, memory:'16 MB', time:'1m']]
    ))
}
