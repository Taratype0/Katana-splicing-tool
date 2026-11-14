// Aggregate completed logs only. No read alignment, BAM scan or QC rerun.
process KATANA_MULTIQC {
    label 'process_single'
    cpus 1
    memory { 4.GB * task.attempt }
    time { 1.h * task.attempt }
    conda 'bioconda::multiqc=1.35'
    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c8/c8e346f4f6080eadf1253505e6ff09ef004454fc18e8d672006fd7b222cc412e/data'
        : 'community.wave.seqera.io/library/multiqc:1.35--c17fb751507e9dfc'}"
    input:
    path qc_files, stageAs: 'inputs/??????/*'
    output:
    path 'multiqc_report.html', emit: report
    path 'multiqc_data/multiqc_data.json', emit: data
    script:
    """
    multiqc inputs --filename multiqc_report.html --outdir . --data-format json --force
    test -s multiqc_report.html
    test -s multiqc_data/multiqc_data.json
    """
}
