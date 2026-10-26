import groovy.json.JsonSlurper
import groovy.json.JsonOutput

include { KATANA_FLOW_PLAN; KATANA_FLOW_VALIDATE_PLAN; KATANA_FLOW_COLLECT; KATANA_FLOW_REPORT } from '../modules/local/flow_control'
include { KATANA_FLOW_TASK as AUDIT } from '../modules/local/flow_task'
include { KATANA_FLOW_TASK as QUALITY } from '../modules/local/flow_task'
include { KATANA_FLOW_TASK as AS } from '../modules/local/flow_task'
include { KATANA_FLOW_TASK as DEG } from '../modules/local/flow_task'
include { KATANA_FLOW_TASK as DTU } from '../modules/local/flow_task'
include { KATANA_FLOW_TASK as LANDSCAPE } from '../modules/local/flow_task'
include { KATANA_FLOW_TASK as EVIDENCE } from '../modules/local/flow_task'
include { KATANA_FLOW_TASK as ENRICHMENT } from '../modules/local/flow_task'
include { KATANA_FLOW_TASK as MOTIFTRACE } from '../modules/local/flow_task'
include { KATANA_FLOW_RESOURCE_TASK as MOTIFTRACE_RESOURCE } from '../modules/local/flow_task_resources'
include { KATANA_FLOW_TASK as OVERLAP } from '../modules/local/flow_task'
include { KATANA_FLOW_TASK as QC } from '../modules/local/flow_task'
include { KATANA_FLOW_TASK as REFERENCE } from '../modules/local/flow_task'
include { KATANA_FLOW_TASK as COVERAGE } from '../modules/local/flow_task'
include { KATANA_FLOW_TASK as NETWORK_FIT } from '../modules/local/flow_task'
include { KATANA_FLOW_TASK as NETWORK_FOLLOWUP } from '../modules/local/flow_task'
include { KATANA_FLOW_TASK as INTERACTION } from '../modules/local/flow_task'
include { KATANA_FLOW_TASK as MAPP } from '../modules/local/flow_task'
include { KATANA_FLOW_TASK as MAPP_PREPARE } from '../modules/local/flow_task'

def tasksForGroup(planFile, group) {
    def plan = new JsonSlurper().parse(planFile.toFile())
    def parents = [audit: [], quality: ['audit'], as: ['quality'], deg: ['quality'],
        dtu: ['quality'], landscape: ['quality'], evidence: ['as', 'deg', 'dtu'],
        enrichment: ['evidence'], motiftrace: ['as'], overlap: ['as', 'deg', 'dtu'],
        reference: ['audit'], qc: ['audit', 'reference'], coverage: ['qc'], network_fit: ['audit'],
        network_followup: ['network_fit'], interaction: ['audit'],
        mapp_prepare: ['as'], mapp: ['audit', 'mapp_prepare']]
    if (!(plan.tasks instanceof List) || plan.tasks.isEmpty()) error('The Katana task plan is empty')
    def ids = plan.tasks.collect { it.id }
    if (ids.unique(false).size() != ids.size()) error('The Katana task plan has duplicate IDs')
    plan.tasks.each { task ->
        if (!(task.id ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) error("Unsafe task ID: ${task.id}")
        if (!parents.containsKey(task.execution_group)) error("Unsupported execution group: ${task.execution_group}")
        (task.dependencies ?: []).each { dependency ->
            def parent = plan.tasks.find { it.id == dependency }
            if (!parent || !(parent.execution_group in parents[task.execution_group])) {
                error("Invalid task dependency group: ${task.id} -> ${dependency}")
            }
        }
    }
    plan.tasks.findAll { it.execution_group == group }.collect { tuple(it, JsonOutput.toJson(plan.execution_plans ? plan.execution_plans[it.id] : plan)) }
}

def readyTasks(plan, group, completed) {
    def tasks = plan.flatMap { tasksForGroup(it, group) }
    def direct = tasks.filter { spec, planFile -> !(spec.dependencies ?: []) }
        .map { spec, planFile -> tuple(spec, planFile, []) }
    def waiting = tasks.filter { spec, planFile -> (spec.dependencies ?: []) }
        .flatMap { spec, planFile -> spec.dependencies.collect { dependency -> tuple(dependency, spec, planFile) } }
    def finished = completed.map { result -> tuple(result.name, result) }
    def dependent = waiting.combine(finished, by: 0)
        .map { dependency, spec, planFile, result -> tuple(groupKey(spec.id, spec.dependencies.size()), spec, planFile, result) }
        .groupTuple()
        .map { key, specs, plans, results -> tuple(specs[0], plans[0], results.sort(false) { it.name }) }
    direct.mix(dependent)
}

def motifRuntimeInputs(spec) {
    // Custom models are task data, not installed software. Declare them as path
    // inputs so Docker/Apptainer can see them and changed bytes invalidate cache.
    def model = spec.parameters?.custom_motif_path
    model ? [file(model, checkIfExists: true)] : []
}

workflow KATANA_DOWNSTREAM {
    take:
    project_root
    project_filename
    flow_config
    flow_config_filename

    main:
    if (params.flow_plan) {
        // A validation receipt is a gate, never a new copy of the task plan.
        // Preserve the original Path so Nextflow can reuse completed tasks.
        original_plan = Channel.value(file(params.flow_plan, checkIfExists: true))
        KATANA_FLOW_VALIDATE_PLAN(original_plan, project_root, project_filename)
        plan = original_plan.combine(KATANA_FLOW_VALIDATE_PLAN.out.receipt)
            .map { original, checked -> original }
    } else {
        // Directory metadata alone does not bind edits inside a configuration directory.
        // These small content inputs invalidate planning, while workers consume only their ancestry.
        configuration_content = flow_config.map { it.resolve(flow_config_filename).text }
        project_manifest_content = project_root.map { it.resolve(project_filename).text }
        KATANA_FLOW_PLAN(project_root, project_filename, flow_config, flow_config_filename,
            configuration_content, project_manifest_content)
        plan = KATANA_FLOW_PLAN.out.plan
    }

    // Each task starts when its own named dependency directories arrive.
    // There are no whole-level barriers, maxForks or application queue limits.
    AUDIT(readyTasks(plan, 'audit', Channel.empty()), project_root, project_filename)
    QUALITY(readyTasks(plan, 'quality', AUDIT.out.result), project_root, project_filename)
    AS(readyTasks(plan, 'as', QUALITY.out.result), project_root, project_filename)
    DEG(readyTasks(plan, 'deg', QUALITY.out.result), project_root, project_filename)
    DTU(readyTasks(plan, 'dtu', QUALITY.out.result), project_root, project_filename)
    LANDSCAPE(readyTasks(plan, 'landscape', QUALITY.out.result), project_root, project_filename)
    molecular = AS.out.result.mix(DEG.out.result, DTU.out.result)
    EVIDENCE(readyTasks(plan, 'evidence', molecular), project_root, project_filename)
    ENRICHMENT(readyTasks(plan, 'enrichment', EVIDENCE.out.result), project_root, project_filename)
    // Classify each ready task once. A single branch makes the two routes
    // mutually exclusive, including when tasks arrive concurrently.
    motif_routes = readyTasks(plan, 'motiftrace', AS.out.result).branch {
        external: !!it[0].parameters?.custom_motif_path
        bundled: true
    }
    MOTIFTRACE(motif_routes.bundled, project_root, project_filename)
    motif_resource_tasks = motif_routes.external
        .map { spec, planFile, results -> tuple(spec, planFile, results, motifRuntimeInputs(spec)) }
    MOTIFTRACE_RESOURCE(motif_resource_tasks, project_root, project_filename)
    OVERLAP(readyTasks(plan, 'overlap', molecular), project_root, project_filename)
    REFERENCE(readyTasks(plan, 'reference', AUDIT.out.result), project_root, project_filename)
    QC(readyTasks(plan, 'qc', AUDIT.out.result.mix(REFERENCE.out.result)), project_root, project_filename)
    COVERAGE(readyTasks(plan, 'coverage', QC.out.result), project_root, project_filename)
    NETWORK_FIT(readyTasks(plan, 'network_fit', AUDIT.out.result), project_root, project_filename)
    NETWORK_FOLLOWUP(readyTasks(plan, 'network_followup', NETWORK_FIT.out.result), project_root, project_filename)
    INTERACTION(readyTasks(plan, 'interaction', AUDIT.out.result), project_root, project_filename)
    MAPP_PREPARE(readyTasks(plan, 'mapp_prepare', AS.out.result), project_root, project_filename)
    MAPP(readyTasks(plan, 'mapp', AUDIT.out.result.mix(MAPP_PREPARE.out.result)), project_root, project_filename)

    all_results = AUDIT.out.result.mix(QUALITY.out.result, AS.out.result, DEG.out.result,
        DTU.out.result, LANDSCAPE.out.result, EVIDENCE.out.result, ENRICHMENT.out.result,
        MOTIFTRACE.out.result, MOTIFTRACE_RESOURCE.out.result, OVERLAP.out.result, REFERENCE.out.result, QC.out.result, COVERAGE.out.result,
        NETWORK_FIT.out.result, NETWORK_FOLLOWUP.out.result, INTERACTION.out.result,
        MAPP_PREPARE.out.result, MAPP.out.result).toSortedList { a, b -> a.name <=> b.name }
    // Collection/report verify producer input hashes too; the same external
    // model files must remain visible there without copying scientific outputs.
    runtime_data = plan.map { planFile ->
        new JsonSlurper().parse(planFile.toFile()).tasks.collectMany { motifRuntimeInputs(it) }
            .unique().sort { it.toString() }
    }
    KATANA_FLOW_COLLECT(plan, all_results, runtime_data)
    KATANA_FLOW_REPORT(KATANA_FLOW_COLLECT.out.collection, runtime_data)

    emit:
    report = KATANA_FLOW_REPORT.out.report
    collection = KATANA_FLOW_COLLECT.out.collection
    task_results = all_results
    task_plan = plan
}
