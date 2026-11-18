process KATANA_SCSPLICETRACE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/katana-sc-splicetrace:2.4.1' :
        'ghcr.io/taratype0/katana-sc-splicetrace:2.4.1' }"

    input:
    tuple val(meta), path(events), path(psi), path(metadata), path(support)
    path fasta
    path fai
    path gtf
    path splicenet_modules
    val numerator
    val denominator

    output:
    tuple val(meta), path("${prefix}/*"), emit: results
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def sup_arg = (support && support.toString() != '') ? "--support ${support}" : ""
    def gtf_arg = (gtf && gtf.toString() != '') ? "--gtf ${gtf}" : ""
    def mod_arg = (splicenet_modules && splicenet_modules.toString() != '') ? "--splicenet-modules ${splicenet_modules}" : ""
    """
    katana-sc-splicetrace \
        --events ${events} \
        --psi ${psi} \
        --metadata ${metadata} \
        ${sup_arg} \
        --fasta ${fasta} \
        ${gtf_arg} \
        ${mod_arg} \
        --threads ${task.cpus} \
        --numerator "${numerator}" \
        --denominator "${denominator}" \
        --output-dir ${prefix} \
        --emit-versions-yml versions.yml \
        --process-name "${task.process}" \
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix} && touch ${prefix}/sc_splicetrace_rbp_ranking_v24.tsv ${prefix}/sc_splicetrace_summary_v24.json
    katana-sc-splicetrace --emit-versions-yml versions.yml --process-name "${task.process}"
    """
}
