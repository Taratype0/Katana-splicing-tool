// Split import and differential models to reuse two pinned nf-core container
// specifications; neither image is assumed to contain the other's package.
// Import pin: nf-core/rnasplice 1.0.4 modules/local/tximport.nf.
// DESeq2 pin: nf-core/differentialabundance 1.5.0 DESEQ2_DIFFERENTIAL.
process KATANA_TXIMPORT {
    tag "${meta.id}"
    label 'process_medium'
    conda 'bioconda::bioconductor-tximeta=1.8.0'
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-tximeta:1.8.0--r40_0' :
        'biocontainers/bioconductor-tximeta:1.8.0--r40_0'}"
    publishDir "${params.outdir}/DEG", mode: 'copy'

    input:
    val meta
    path quant_dirs
    path samplesheet
    path gtf
    val extra_args

    output:
    tuple val(meta), path('gene_txi.rds'), emit: gene_txi
    tuple val(meta), path('annotation.tsv'), emit: annotation
    tuple val(meta), path('tx2gene.tsv'), emit: tx2gene
    tuple val(meta), path('gene.*.tsv'), emit: gene_matrices
    tuple val(meta), path('tx2gene_coverage.tsv'), emit: coverage
    tuple val(meta), path('import_metadata.tsv'), emit: run_metadata
    tuple val(meta), path('import.R_sessionInfo.log'), emit: session_info
    path 'statistics.import.versions.yml', emit: versions

    script:
    def args = extra_args ?: ''
    """
    katana_differential.R \\
        --stage import \\
        --samples '${samplesheet}' \\
        --gtf '${gtf}' \\
        --quant-root . \\
        --outdir . \\
        ${args}
    """
}

process KATANA_DESEQ2_CONTRAST {
    tag "${contrast_id}"
    label 'process_medium'
    conda 'bioconda::bioconductor-deseq2=1.34.0'
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-deseq2:1.34.0--r41hc247a5b_3' :
        'biocontainers/bioconductor-deseq2:1.34.0--r41hc247a5b_3'}"
    publishDir "${params.outdir}/DEG", mode: 'copy'

    input:
    tuple val(meta), path(gene_txi)
    path samplesheet
    tuple val(contrast_id), path(contrastsheet)
    val design
    val extra_args

    output:
    tuple val(meta), path('*.deseq2.results.tsv'), emit: results
    tuple val(meta), path('*.normalised_counts.tsv'), emit: contrast_normalised_counts
    tuple val(meta), path('*.dds.rds'), emit: models_rds
    tuple val(meta), path('*.deseq2.model.txt'), emit: model_formulas
    path 'models.*.tsv', emit: model_summary
    path 'differential_metadata.*.tsv', emit: run_metadata
    tuple val(meta), path('differential.*.R_sessionInfo.log'), emit: session_info
    path 'statistics.differential.*.versions.yml', emit: versions

    script:
    def args = extra_args ?: ''
    def formula = design ?: '~ condition'
    """
    katana_differential.R \\
        --stage contrast \\
        --samples '${samplesheet}' \\
        --contrasts '${contrastsheet}' \\
        --gene-rds '${gene_txi}' \\
        --design '${formula}' \\
        --outdir . \\
        ${args}
    """
}

process KATANA_NORMALIZE {
    tag "${meta.id}"
    label 'process_medium'
    conda 'bioconda::bioconductor-deseq2=1.34.0'
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-deseq2:1.34.0--r41hc247a5b_3' :
        'biocontainers/bioconductor-deseq2:1.34.0--r41hc247a5b_3'}"
    publishDir "${params.outdir}/DEG", mode: 'copy'
    input:
    tuple val(meta), path(gene_txi)
    path samplesheet
    val extra_args
    output:
    tuple val(meta), path('all.normalised_counts.tsv'), emit: normalised_counts
    tuple val(meta), path('all.normalization_factors.tsv'), emit: normalization_factors
    tuple val(meta), path('gene.deseq2_integer_counts.tsv'), emit: integer_counts
    path 'statistics.normalization.versions.yml', emit: versions
    script:
    def args = extra_args ?: ''
    """
    katana_differential.R --stage normalize --samples '${samplesheet}' --gene-rds '${gene_txi}' --outdir . ${args}
    """

    stub:
    """
    touch all.normalised_counts.tsv all.normalization_factors.tsv gene.deseq2_integer_counts.tsv
    printf 'KATANA_NORMALIZE: {mode: stub}\\n' > statistics.normalization.versions.yml
    """
}

process KATANA_DESEQ2_METADATA {
    label 'process_single'
    container 'docker.io/library/python:3.12.11-bookworm'
    publishDir "${params.outdir}/DEG", mode: 'copy'
    input:
    val meta
    path contrastsheet
    path model_files, stageAs: 'models/*'
    path metadata_files, stageAs: 'metadata/*'
    path merger
    output:
    tuple val(meta), path('models.tsv'), emit: models
    tuple val(meta), path('differential_metadata.tsv'), emit: metadata
    script:
    """
    python3 '${merger}' --contrasts '${contrastsheet}' --models models/*.tsv --metadata metadata/*.tsv
    """
}

include { CONTRAST_TASK_SHEETS } from './split_contrasts'

workflow KATANA_DIFFERENTIAL {
    take:
    meta
    quant_dirs
    samplesheet
    contrastsheet
    gtf
    design
    extra_args

    main:
    KATANA_TXIMPORT(meta, quant_dirs, samplesheet, gtf, extra_args)
    CONTRAST_TASK_SHEETS(contrastsheet.first())
    ch_gene_txi = KATANA_TXIMPORT.out.gene_txi.first()
    ch_samples = samplesheet.first()
    KATANA_NORMALIZE(ch_gene_txi, ch_samples, extra_args)
    KATANA_DESEQ2_CONTRAST(ch_gene_txi, ch_samples, CONTRAST_TASK_SHEETS.out.sheets, design, extra_args)
    KATANA_DESEQ2_METADATA(meta, contrastsheet.first(),
        KATANA_DESEQ2_CONTRAST.out.model_summary.flatten().collect(),
        KATANA_DESEQ2_CONTRAST.out.run_metadata.flatten().collect(),
        Channel.value(file("${projectDir}/bin/kanata_merge_deseq_metadata.py", checkIfExists: true)))

    emit:
    results = KATANA_DESEQ2_CONTRAST.out.results
    normalised_counts = KATANA_NORMALIZE.out.normalised_counts
    annotation = KATANA_TXIMPORT.out.annotation
    tx2gene = KATANA_TXIMPORT.out.tx2gene
    run_metadata = KATANA_TXIMPORT.out.run_metadata.mix(KATANA_DESEQ2_METADATA.out.metadata)
    session_info = KATANA_TXIMPORT.out.session_info.mix(KATANA_DESEQ2_CONTRAST.out.session_info)
    coverage = KATANA_TXIMPORT.out.coverage
    gene_matrices = KATANA_TXIMPORT.out.gene_matrices
    normalization_factors = KATANA_NORMALIZE.out.normalization_factors
    models = KATANA_DESEQ2_METADATA.out.models
    versions = KATANA_TXIMPORT.out.versions.mix(KATANA_NORMALIZE.out.versions, KATANA_DESEQ2_CONTRAST.out.versions, CONTRAST_TASK_SHEETS.out.versions)
}
