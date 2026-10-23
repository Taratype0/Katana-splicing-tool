nextflow.enable.dsl = 2

// Installation/transport acceptance only: no study matrices or model fits.
process IMAGE_RUNTIME {
    tag { probe_kind }
    container params.image
    cpus 1
    memory '2 GB'
    time '10 min'
    publishDir params.outdir, mode: 'copy'

    input:
    val probe_kind
    path checker

    output:
    path "${probe_kind}"

    script:
    """
    export KATANA_CONTAINER=1
    unset PYTHONPATH PYTHONHOME PYTHONUSERBASE R_HOME R_LIBS R_LIBS_USER R_LIBS_SITE
    export PYTHONNOUSERSITE=1 PYTHONSAFEPATH=1
    python ${checker} ${probe_kind}
    """
}

workflow {
    if (!params.image || params.uid == null || params.gid == null || !params.outdir) {
        error('Use check_container_runtime.py with an existing local image and a new output directory')
    }
    IMAGE_RUNTIME(Channel.of('network', 'mapp', 'motif'), file("${projectDir}/inside.py"))
}
