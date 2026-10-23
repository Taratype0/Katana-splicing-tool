#!/usr/bin/env nextflow
// Unified orchestration. Scientific upstream modules retain their copied source
// identities; Katana execution and reporting use the installed Python package.
include { KATANA_PIPELINE } from './workflows/katana_core'
include { KATANA_DOWNSTREAM } from './workflows/katana_downstream'
include { KATANA_FLOW_PREFLIGHT; KATANA_FLOW_IMPORT } from './modules/local/flow_control'

def requireParameters(names) {
    names.each { name -> if (!params[name]) error("Required parameter missing: --${name}") }
}

def inputDirectory(path) {
    def directory = path.toAbsolutePath().normalize()
    def outputs = file(params.outdir).toAbsolutePath().normalize()
    def work = workflow.workDir.toAbsolutePath().normalize()
    if (outputs.startsWith(directory) || work.startsWith(directory)) {
        error("Input directory ${directory} contains the output/work directory. Keep configuration/project inputs in a separate directory.")
    }
    directory
}

workflow {
    if (params.help) {
        log.info '''Katana bulk analysis and report
  --mode complete         --input samples.csv --contrasts contrasts.csv --fasta genome.fa --gtf annotation.gtf
                          --project_template project-template.json --flow_config flow.json
  --mode upstream-only    --input samples.csv --contrasts contrasts.csv --fasta genome.fa --gtf annotation.gtf
  --mode downstream-only  --project project.json --flow_config flow.json
  --outdir results        -profile local|slurm,docker|singularity|apptainer
Existing completed analysis: use downstream-only. Use -resume for this run's Nextflow cache.
Required calculations are verified before the fixed seven-section report is accepted.'''
    } else {
        if (!(params.mode in ['complete', 'upstream-only', 'downstream-only'])) error('Invalid --mode')
        if (!(params.katana_stage_containers instanceof Map)) error('--katana_stage_containers must be a map of execution groups to images')
        def groups = ['control', 'report', 'audit', 'quality', 'as', 'deg', 'dtu', 'landscape', 'evidence',
                      'enrichment', 'motiftrace', 'overlap', 'reference', 'qc', 'coverage', 'network_fit',
                      'network_followup', 'interaction', 'mapp_prepare', 'mapp']
        if (params.katana_stage_containers.any { key, value -> !(key in groups) || !(value instanceof String) || !value.trim() }) {
            error('Katana stage containers require known execution-group keys and nonempty image references')
        }
        if (params.mode != 'upstream-only' && (workflow.containerEngine in ['docker', 'singularity', 'apptainer']) && !params.katana_container) {
            error('A Katana container profile requires --katana_container pointing to your built full runtime image; host fallback is disabled')
        }
        if (params.mode != 'upstream-only' && (workflow.containerEngine in ['docker', 'singularity', 'apptainer']) && !(params.flow_executable in ['katana-flow', '/opt/conda/bin/katana-flow'])) {
            error('Container profiles require the image launcher katana-flow or /opt/conda/bin/katana-flow; mounted host launchers are not allowed')
        }
        if (params.mode == 'complete' && params.skip_fastqc) {
            error('Complete fixed-report mode requires native sequencing FastQC. Use upstream-only for an explicitly skipped QC stage.')
        }
        if (params.mode == 'downstream-only') {
            requireParameters(['project', 'flow_config'])
            project = file(params.project, checkIfExists: true)
            config_file = file(params.flow_config, checkIfExists: true)
            KATANA_DOWNSTREAM(Channel.value(inputDirectory(project.parent)), project.name,
                Channel.value(inputDirectory(config_file.parent)), config_file.name)
        } else {
            requireParameters(['input', 'contrasts', 'fasta', 'gtf'])
            samples = Channel.value(file(params.input, checkIfExists: true))
            contrasts = Channel.value(file(params.contrasts, checkIfExists: true))
            fasta = Channel.value(file(params.fasta, checkIfExists: true))
            gtf = Channel.value(file(params.gtf, checkIfExists: true))
            if (params.mode == 'complete') {
                requireParameters(['project_template', 'flow_config'])
                template_file = file(params.project_template, checkIfExists: true)
                template = Channel.value(inputDirectory(template_file.parent))
                config_file = file(params.flow_config, checkIfExists: true)
                config = Channel.value(inputDirectory(config_file.parent))
                KATANA_FLOW_PREFLIGHT(template, template_file.name, config, config_file.name)
                samples = samples.combine(KATANA_FLOW_PREFLIGHT.out.receipt).map { sheet, ready -> sheet }.first()
                contrasts = contrasts.combine(KATANA_FLOW_PREFLIGHT.out.receipt).map { sheet, ready -> sheet }.first()
                fasta = fasta.combine(KATANA_FLOW_PREFLIGHT.out.receipt).map { ref, ready -> ref }.first()
                gtf = gtf.combine(KATANA_FLOW_PREFLIGHT.out.receipt).map { ref, ready -> ref }.first()
            }
            KATANA_PIPELINE(samples, contrasts, fasta, gtf)
            if (params.mode == 'complete') {
                KATANA_FLOW_IMPORT(KATANA_PIPELINE.out.bundle, template, template_file.name)
                KATANA_DOWNSTREAM(KATANA_FLOW_IMPORT.out.workspace.first(), 'project.json', config, config_file.name)
            }
        }
    }
}
