def containerRuntimePrelude(engine) {
    return engine ? '''export KATANA_CONTAINER=1
    unset PYTHONPATH PYTHONHOME PYTHONUSERBASE R_HOME R_LIBS R_LIBS_USER R_LIBS_SITE R_ENVIRON R_ENVIRON_USER R_PROFILE R_PROFILE_USER VIRTUAL_ENV
    export PYTHONNOUSERSITE=1 PYTHONSAFEPATH=1''' : ''
}

process KATANA_FLOW_PREFLIGHT {
    label 'katana_control'
    input:
    path project_template_root, stageAs: 'project-template'
    val project_template_filename
    path flow_config, stageAs: 'flow-config'
    val flow_config_filename
    output:
    path 'preflight.json', emit: receipt
    script:
    def quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    """
    ${containerRuntimePrelude(workflow.containerEngine)}
    ${quote(params.flow_executable)} preflight --config ${quote(flow_config.resolve(flow_config_filename))} \\
        --project-template ${quote(project_template_root.resolve(project_template_filename))} --out preflight.json
    test -s preflight.json
    """
}

process KATANA_FLOW_IMPORT {
    label 'katana_control'

    input:
    path bundle
    path project_template_root, stageAs: 'project-template'
    val project_template_filename

    output:
    path 'workspace', emit: workspace

    script:
    def quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    """
    ${containerRuntimePrelude(workflow.containerEngine)}
    ${quote(params.flow_executable)} import --bundle ${quote(bundle)} \\
        --project-template ${quote(project_template_root.resolve(project_template_filename))} --outdir workspace
    test -s workspace/project.json
    """
}

process KATANA_FLOW_PLAN {
    label 'katana_control'

    input:
    path project_root, stageAs: 'project'
    val project_filename
    path flow_config, stageAs: 'flow-config'
    val flow_config_filename

    output:
    path 'plan.json', emit: plan

    script:
    def quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    """
    ${containerRuntimePrelude(workflow.containerEngine)}
    ${quote(params.flow_executable)} plan --require-ready --project ${quote(project_root.resolve(project_filename))} \\
        --config ${quote(flow_config.resolve(flow_config_filename))} --out plan.json
    test -s plan.json
    """
}

process KATANA_FLOW_COLLECT {
    label 'katana_control'

    input:
    path plan
    path task_results, stageAs: 'task-results/*'

    output:
    path 'collection', emit: collection

    script:
    def quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    def results = task_results instanceof List ? task_results : [task_results]
    if (!results) error('Cannot collect an empty task plan')
    def result_args = results.collect { quote(it) }.join(' ')
    """
    ${containerRuntimePrelude(workflow.containerEngine)}
    ${quote(params.flow_executable)} collect --plan ${quote(plan)} \\
        --task-results ${result_args} --outdir collection
    test -s collection/collection.receipt.json
    """
}

process KATANA_FLOW_REPORT {
    label 'katana_report'

    input:
    path collection

    output:
    path 'report', emit: report

    script:
    def quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    """
    ${containerRuntimePrelude(workflow.containerEngine)}
    ${quote(params.flow_executable)} report --collection ${quote(collection)} --outdir report
    test -s report/REPORT.html
    test -s report/report-manifest.json
    """
}
