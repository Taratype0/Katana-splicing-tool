// nf-core/rnasplice 2.0.0 task granularity, with the existing rMATS 4.1.2
// runtime and analysis arguments retained for numerical compatibility.
process RMATS_PREP_BAM {
    tag "${meta.id}"
    label 'process_single'

    conda 'bioconda::r-pairadise=1.0.0 bioconda::rmats=4.1.2'
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8ea76ff0a6a4c7e5c818fd4281abf918f92eeeae:121e48ab4817ec619c157a346458efca1ccf3c0a-0' :
        'quay.io/biocontainers/mulled-v2-8ea76ff0a6a4c7e5c818fd4281abf918f92eeeae:121e48ab4817ec619c157a346458efca1ccf3c0a-0'}"

    input:
    tuple val(meta), path(bam)
    path gtf
    val read_length
    val splice_diff_cutoff
    val novel_splice_site
    val min_intron_len
    val max_exon_len

    output:
    tuple val(meta), path("${meta.id}.rmats"), emit: rmats_temp
    path "${meta.id}.rmats_prep.log", emit: log
    path 'versions.yml', emit: versions

    script:
    def read_type = meta.single_end ? 'single' : 'paired'
    def strand = [forward:'fr-secondstrand', reverse:'fr-firststrand', unstranded:'fr-unstranded'][meta.strandedness]
    if (!strand) error("Unresolved rMATS strandedness for ${meta.id}")
    def novel = novel_splice_site ? "--novelSS --mil ${min_intron_len ?: 50} --mel ${max_exon_len ?: 500}" : ''
    """
    printf '%s\\n' '${bam}' > bamlist.txt
    rmats.py --task prep --gtf ${gtf} --b1 bamlist.txt \\
        --od rmats_prep --tmp rmats_temp -t ${read_type} --libType ${strand} \\
        --readLength ${read_length} --variable-read-length --allow-clipping \\
        --nthread ${task.cpus} --tstat ${task.cpus} --cstat ${splice_diff_cutoff} \\
        ${novel} > ${meta.id}.rmats_prep.log
    # A single BAM produces one record. Renaming avoids cross-task timestamp
    # collisions when independent prep jobs finish at the same time.
    cp rmats_temp/*.rmats ${meta.id}.rmats
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rmats: \$(rmats.py --version | sed 's/v//g')
    END_VERSIONS
    """
}
