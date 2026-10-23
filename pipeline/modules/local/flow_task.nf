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

    script:
    def quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    def executable = quote(params.flow_executable)
    def previous = dependency_results instanceof List ? dependency_results : [dependency_results]
    def dependency_args = previous ? '--dependency-receipts ' + previous.collect { quote(it.resolve('task.receipt.json')) }.join(' ') : ''
    """
    ${executable} task --plan ${quote(plan)} --project ${quote(project_root.resolve(project_filename))} --task-id ${quote(spec.id)} \\
        --outdir ${quote("task-results/${spec.id}")} ${dependency_args}
    test -s ${quote("task-results/${spec.id}/task.receipt.json")}
    """
}
