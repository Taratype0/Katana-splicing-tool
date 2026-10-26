process POLICY_CHECK {
    label 'katana_downstream'
    tag "${spec.id}"
    input:
    val spec
    output:
    path 'allocation.tsv'
    script:
    """
    printf '%s\t%s\t%s\t%s\t%s\n' '${spec.id}' '${task.cpus}' '${task.memory}' '${task.time}' '${task.queue ?: "unset"}' > allocation.tsv
    """
}
