process DEXSEQ_DTU {
    tag "${contrast_id}"
    label 'process_high'

    conda "bioconda::bioconductor-dexseq=1.36.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dexseq:1.36.0--r40_0' :
        'quay.io/biocontainers/bioconductor-dexseq:1.36.0--r40_0' }"

    input:
    path drimseq_sample_data
    path drimseq_d_counts
    tuple val(contrast_id), path(drimseq_contrast_data)
    val ntop

    output:
    path "DEXSeqDataSet.*.rds"  , emit: dexseq_exon_dataset_rds
    path "DEXSeqResults.*.rds"  , emit: dexseq_exon_results_rds
    path "DEXSeqResults.*.tsv"  , emit: dexseq_exon_results_tsv
    path "perGeneQValue.*.rds"  , emit: dexseq_gene_results_rds
    path "perGeneQValue.*.tsv"  , emit: dexseq_gene_results_tsv
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export R_BIOC_WORKERS=${task.cpus}
    export KANATA_SINGLE_CONTRAST=1
    run_dexseq_dtu.R $drimseq_sample_data $drimseq_contrast_data $drimseq_d_counts $ntop

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-dexseq:  \$(Rscript -e "library(DEXSeq); cat(as.character(packageVersion('DEXSeq')))")
    END_VERSIONS
    """

    stub:
    """
    python3 - '${drimseq_contrast_data}' <<'PY'
    import csv
    from pathlib import Path
    row = next(csv.DictReader(open(__import__('sys').argv[1])))
    name = row['treatment'] + '-' + row['control']
    for prefix, extension in [('DEXSeqDataSet','rds'),('DEXSeqResults','rds'),('DEXSeqResults','tsv'),('perGeneQValue','rds'),('perGeneQValue','tsv')]:
        Path(prefix+'.'+name+'.'+extension).touch()
    PY
    printf 'DEXSEQ_DTU: {mode: stub}\\n' > versions.yml
    """

}
