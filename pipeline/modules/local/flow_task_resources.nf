// Explicit external motif inputs participate in staging and Nextflow cache identity.
// They remain read-only data; configured executables are checked separately.
// Results enter later levels through declared channels, never publishDir scans.
process KATANA_FLOW_RESOURCE_TASK {
    tag "${spec.id}"
    label 'katana_downstream'

    input:
    tuple val(spec), val(plan_text), path(dependency_results, stageAs: 'dependencies/*'), path(runtime_data_files, stageAs: 'runtime-data/??????/*')
    path project_root, stageAs: 'project'
    val project_filename

    output:
    path "task-results/${spec.id}", emit: result
    path "runtime-${spec.id}.json", emit: runtime_probe

    script:
    def quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    def executable = quote(params.flow_executable)
    def encoded_plan = plan_text.bytes.encodeBase64().toString()
    def isolation = workflow.containerEngine ? '''export KATANA_CONTAINER=1
    unset PYTHONPATH PYTHONHOME PYTHONUSERBASE R_HOME R_LIBS R_LIBS_USER R_LIBS_SITE R_ENVIRON R_ENVIRON_USER R_PROFILE R_PROFILE_USER VIRTUAL_ENV
    export PYTHONNOUSERSITE=1 PYTHONSAFEPATH=1''' : ''
    def previous = dependency_results instanceof List ? dependency_results : [dependency_results]
    def dependency_args = previous ? '--dependency-receipts ' + previous.collect { quote(it.resolve('task.receipt.json')) }.join(' ') : ''
    """
    ${isolation}
    printf '%s' '${encoded_plan}' | base64 -d > task-plan.json
    export KATANA_CPUS=${task.cpus}
    export OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} MKL_NUM_THREADS=${task.cpus} BLIS_NUM_THREADS=${task.cpus} NUMEXPR_NUM_THREADS=${task.cpus}
    ${executable} runtime-check --plan task-plan.json --project ${quote(project_root.resolve(project_filename))} \\
        --task-id ${quote(spec.id)} --runtime-label ${quote(task.container ?: 'configured-host')} --out ${quote("runtime-${spec.id}.json")}
    ${executable} task --plan task-plan.json --project ${quote(project_root.resolve(project_filename))} --task-id ${quote(spec.id)} \\
        --outdir ${quote("task-results/${spec.id}")} ${dependency_args}
    test -s ${quote("task-results/${spec.id}/task.receipt.json")}
    """
}
