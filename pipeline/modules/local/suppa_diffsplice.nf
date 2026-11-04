process DIFFSPLICE {
    tag "${cond1}-${cond2}"
    label 'process_high'

    conda "bioconda::suppa=2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/suppa:2.3--py36_0' :
        'quay.io/biocontainers/suppa:2.3--py36_0' }"

    input:
    path events
    tuple val(cond1), val(cond2), path(tpm1), path(tpm2)
    tuple val(cond1), val(cond2), path(psi1), path(psi2)
    val prefix
    val diffsplice_method             // val params.diffsplice_method
    val diffsplice_area               // val params.diffsplice_area
    val diffsplice_lower_bound        // val params.diffsplice_lower_bound
    val diffsplice_alpha              // val params.diffsplice_alpha
    val diffsplice_tpm_threshold      // val params.diffsplice_tpm_threshold
    val diffsplice_nan_threshold      // val params.diffsplice_nan_threshold
    val diffsplice_gene_correction    // val params.diffsplice_gene_correction
    val diffsplice_paired             // val params.diffsplice_paired
    val diffsplice_median             // val params.diffsplice_median

    output:
    tuple val(cond1), val(cond2), path("*.dpsi")       , emit: dpsi
    tuple val(cond1), val(cond2), path("*.psivec")     , emit: psivec
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: //  Calculate differential analysis between conditions

    def gc     = diffsplice_gene_correction ? "-gc" : ''
    def is_paired = diffsplice_paired instanceof Map ? diffsplice_paired[cond1 + '\t' + cond2] : diffsplice_paired
    if (is_paired == null) error("Missing paired-model declaration for ${cond1}-${cond2}")
    def excluded = is_paired instanceof Map ? (is_paired.exclude_samples ?: '') : ''
    is_paired = is_paired instanceof Map ? is_paired.paired : is_paired
    def selected_tpm1 = excluded ? 'selected-1.tpm' : tpm1
    def selected_tpm2 = excluded ? 'selected-2.tpm' : tpm2
    def selected_psi1 = excluded ? 'selected-1.psi' : psi1
    def selected_psi2 = excluded ? 'selected-2.psi' : psi2
    def select_columns = excluded ? "python3 \"\$(command -v katana_subset_suppa.py)\" --exclude '${excluded}' '${tpm1}' '${selected_tpm1}' '${tpm2}' '${selected_tpm2}' '${psi1}' '${selected_psi1}' '${psi2}' '${selected_psi2}'" : ''
    def pa     = is_paired ? "-pa" : ''
    def median = diffsplice_median ? "-me" : ''

    """
    ${select_columns}
    suppa.py \\
        diffSplice \\
        -m $diffsplice_method \\
        $gc $pa -s -c $median \\
        -a $diffsplice_area \\
        -l $diffsplice_lower_bound \\
        -al $diffsplice_alpha \\
        -th $diffsplice_tpm_threshold \\
        -nan $diffsplice_nan_threshold \\
        -i $events \\
        -p $selected_psi1 $selected_psi2 \\
        -e $selected_tpm1 $selected_tpm2 \\
        -o ${cond1}-${cond2}_${prefix}_diffsplice

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        suppa: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('suppa').version)")
    END_VERSIONS
    """
}
