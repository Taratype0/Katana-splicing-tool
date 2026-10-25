// RSeQC interface adapted from nf-core/modules commit
// 05954dab2ff481bcb999f24455da29a5828af08d. See the accompanying provenance JSON.
// The caller supplies only samples explicitly requesting strandedness=auto.

process GFFREAD_STRAND_MODELS {
    label 'process_low'

    conda "${moduleDir}/infer_strandedness.gffread.environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.1--h8b12597_0' :
        'quay.io/biocontainers/gffread:0.12.1--h8b12597_0'}"

    input:
    path gtf

    output:
    path 'strand_models.bed', emit: bed
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gffread "$gtf" --bed -o strand_models.raw.bed
    # Keep exactly BED12, including exon blocks and transcript strand. gffread
    # may append attributes after column 12; they are not part of this contract.
    awk 'BEGIN { FS = OFS = "\\t"; models = 0 }
         /^#/ || /^track/ || /^browser/ || NF == 0 { next }
         {
             if (NF < 12 || \$2 < 0 || \$3 <= \$2 || (\$6 != "+" && \$6 != "-") || \$10 < 1) {
                 print "Invalid stranded BED12 model at line " NR > "/dev/stderr"
                 exit 1
             }
             print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12
             models++
         }
         END {
             if (models == 0) {
                 print "No stranded transcript models were produced from the annotation" > "/dev/stderr"
                 exit 1
             }
         }' strand_models.raw.bed > strand_models.bed
    test -s strand_models.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """

    stub:
    """
    touch strand_models.bed
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: 0.12.1
    END_VERSIONS
    """
}

process RSEQC_INFEREXPERIMENT {
    tag "$meta.id"
    // Native inference streams one alignment at a time on a single CPU.
    label 'process_single'

    conda "${moduleDir}/infer_strandedness.rseqc.environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rseqc:5.0.3--py39hf95cd2a_0' :
        'quay.io/biocontainers/rseqc:5.0.3--py39hf95cd2a_0'}"

    input:
    tuple val(meta), path(bam), path(bai)
    path bed
    val sample_size
    val mapq

    output:
    tuple val(meta), path('*.infer_experiment.txt'), emit: txt
    tuple val(meta), path('*.infer_experiment.log'), emit: log
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (!(sample_size.toString() ==~ /^[0-9]+$/) || !(mapq.toString() ==~ /^[0-9]+$/)) {
        error('strandedness_sample_size and strandedness_mapq must be integers')
    }
    def requested_samples = sample_size as Integer
    def minimum_mapq = mapq as Integer
    if (requested_samples < 1) error('strandedness_sample_size must be a positive integer')
    if (minimum_mapq < 0 || minimum_mapq > 255) error('strandedness_mapq must be between 0 and 255')
    """
    # The shared BAM audit/sort/index route supplies this tuple. RSeQC itself
    # streams alignment records; its usable-read count is retained in stderr.
    test -s "$bam"
    test -s "$bai"
    test -s "$bed"
    infer_experiment.py \\
        -i "$bam" \\
        -r "$bed" \\
        -s $requested_samples \\
        -q $minimum_mapq \\
        > "${prefix}.infer_experiment.txt" \\
        2> "${prefix}.infer_experiment.log"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(infer_experiment.py --version | sed -e "s/infer_experiment.py //g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Empty stub files cannot be mistaken for an observed strandedness call.
    touch "${prefix}.infer_experiment.txt" "${prefix}.infer_experiment.log"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: 5.0.3
    END_VERSIONS
    """
}
