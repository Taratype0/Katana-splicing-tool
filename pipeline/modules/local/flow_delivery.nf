include { containerRuntimePrelude } from './flow_control'

process KATANA_FLOW_RETAIN {
    label 'katana_control'
    input:
    path collection
    path project_root, stageAs: 'source-project'
    val project_filename
    output:
    path 'project', emit: project
    script:
    def quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    """
    ${containerRuntimePrelude(workflow.containerEngine)}
    ${quote(params.flow_executable)} retain-project --collection ${quote(collection)} \\
        --project ${quote(project_root.resolve(project_filename))} --outdir project
    test -s project/retention.receipt.json
    """
}

process KATANA_FLOW_DELIVERY {
    label 'katana_control'
    input:
    path retained_project, stageAs: 'project'
    path report
    path verifier, stageAs: 'verify_delivery.py'
    output:
    path 'delivery.receipt.json', emit: receipt
    script:
    """
    ${containerRuntimePrelude(workflow.containerEngine)}
    python3 "${verifier}" --project "${retained_project}" --report "${report}" --out delivery.receipt.json
    """
}
