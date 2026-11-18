process KATANA_SPLICETRACE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/katana-splicetrace:2.4.3' :
        'ghcr.io/taratype0/katana-splicetrace:2.4.3' }"

    input:
    tuple val(meta), path(rmats_input)
    path fasta
    path fai
    path gtf
    path splicenet_modules

    output:
    tuple val(meta), path("${prefix}/*"), emit: results
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def gtf_arg = (gtf && gtf.toString() != '') ? "--gtf ${gtf}" : ""
    def mod_arg = (splicenet_modules && splicenet_modules.toString() != '') ? "--splicenet-modules ${splicenet_modules}" : ""
    """
    katana-splicetrace \
        --rmats ${rmats_input} \
        --fasta ${fasta} \
        ${gtf_arg} \
        ${mod_arg} \
        --threads ${task.cpus} \
        --output-dir ${prefix} \
        --emit-versions-yml versions.yml \
        --process-name "${task.process}" \
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix} && touch ${prefix}/rbp_ranking_v22.tsv ${prefix}/run_receipt.json
    katana-splicetrace --emit-versions-yml versions.yml --process-name "${task.process}"
    """
}
