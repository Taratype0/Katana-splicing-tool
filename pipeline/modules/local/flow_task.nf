// Every invocation is one independently scheduled, isolated Katana task.
// Results enter later levels through declared channels, never publishDir scans.
process KATANA_FLOW_TASK {
    tag "${spec.id}"
    label 'katana_downstream'

    input:
    tuple val(spec), path(plan), path(dependency_results, stageAs: 'dependencies/*')
    path project_root, stageAs: 'project'
    val project_filename

    output:
    path "task-results/${spec.id}", emit: result
    path "runtime-${spec.id}.json", emit: runtime_probe

    script:
    def quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    def executable = quote(params.flow_executable)
    def isolation = workflow.containerEngine ? '''export KATANA_CONTAINER=1
    unset PYTHONPATH PYTHONHOME PYTHONUSERBASE R_HOME R_LIBS R_LIBS_USER R_LIBS_SITE R_ENVIRON R_ENVIRON_USER R_PROFILE R_PROFILE_USER VIRTUAL_ENV
    export PYTHONNOUSERSITE=1 PYTHONSAFEPATH=1''' : ''
    def previous = dependency_results instanceof List ? dependency_results : [dependency_results]
    def dependency_args = previous ? '--dependency-receipts ' + previous.collect { quote(it.resolve('task.receipt.json')) }.join(' ') : ''
    """
    ${isolation}
    ${executable} runtime-check --plan ${quote(plan)} --project ${quote(project_root.resolve(project_filename))} \\
        --task-id ${quote(spec.id)} --runtime-label ${quote(task.container ?: 'configured-host')} --out ${quote("runtime-${spec.id}.json")}
    ${executable} task --plan ${quote(plan)} --project ${quote(project_root.resolve(project_filename))} --task-id ${quote(spec.id)} \\
        --outdir ${quote("task-results/${spec.id}")} ${dependency_args}
    test -s ${quote("task-results/${spec.id}/task.receipt.json")}
    """
}
