process KATANA_SPLICENET {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/katana-splicenet:2.3.2' :
        'ghcr.io/taratype0/katana-splicenet:2.3.2' }"

    input:
    tuple val(meta), path(usage), path(samples), path(support), path(annotation)
    path gmt

    output:
    tuple val(meta), path("${prefix}/*")                        , emit: results
    path "${prefix}/module_assignments_v2.tsv"                  , emit: module_table
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def sup_arg = (support && support.toString() != '') ? "--support ${support}" : ""
    def ann_arg = (annotation && annotation.toString() != '') ? "--annotation ${annotation}" : ""
    def gmt_arg = (gmt && gmt.toString() != '') ? "--gmt ${gmt}" : ""
    """
    katana-splicenet \
        --usage ${usage} \
        --samples ${samples} \
        ${sup_arg} \
        ${ann_arg} \
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
    mkdir -p ${prefix} && touch ${prefix}/module_assignments_v2.tsv ${prefix}/run_receipt.json
    katana-splicenet --emit-versions-yml versions.yml --process-name "${task.process}"
    """
}
