process KATANA_SCSPLICENET {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/katana-sc-splicenet:2.2.2' :
        'ghcr.io/taratype0/katana-sc-splicenet:2.2.2' }"

    input:
    tuple val(meta), path(usage), path(metadata), path(event_genes), path(support)
    path gmt

    output:
    tuple val(meta), path("${prefix}/*")                        , emit: results
    path "${prefix}/sc_cosplice_modules_v21.tsv"                , emit: module_table
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def sup_arg = (support && support.toString() != '') ? "--support ${support}" : ""
    def gmt_arg = (gmt && gmt.toString() != '') ? "--gmt ${gmt}" : ""
    """
    katana-sc-splicenet \
        --usage ${usage} \
        --metadata ${metadata} \
        --event-genes ${event_genes} \
        ${sup_arg} \
        ${gmt_arg} \
        --threads ${task.cpus} \
        --output-dir ${prefix} \
        --emit-versions-yml versions.yml \
        --process-name "${task.process}" \
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix} && touch ${prefix}/sc_cosplice_edges_v21.tsv ${prefix}/sc_cosplice_modules_v21.tsv
    katana-sc-splicenet --emit-versions-yml versions.yml --process-name "${task.process}"
    """
}
