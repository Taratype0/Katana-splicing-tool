#!/usr/bin/env nextflow
// Unified orchestration. Scientific upstream modules retain their copied source
// identities; Katana execution and reporting use the installed Python package.
include { KATANA_PIPELINE } from './workflows/katana_core'
include { KATANA_DOWNSTREAM } from './workflows/katana_downstream'
include { KATANA_FLOW_CONFIGURE; KATANA_FLOW_PREFLIGHT; KATANA_FLOW_IMPORT } from './modules/local/flow_control'

// Flat public parameters become small, cacheable task inputs. Resource files are
// native path inputs so container mounts and Nextflow invalidation see them.
def automaticConfiguration() {
    def keys = ['mode', 'study_name', 'species', 'reference_build', 'annotation_release',
                'event_set', 'rmats_mode', 'contrast_ids', 'resource_metadata',
                'deg_fdr', 'deg_log2fc', 'das_fdr', 'das_dpsi', 'dtu_fdr',
                'network_min_module_size', 'network_merge_cut_height',
                'with_genebody_coverage', 'with_network_bootstrap']
    keys += ['landscape', 'cross_comparison', 'enrichment', 'rbp_identity', 'motiftrace',
             'motif_binding', 'mapp', 'network', 'interaction', 'sequencing_qc'].collect { "skip_${it}".toString() }
    def request = keys.collectEntries { key -> [(key): params[key]] }
    def names = ['gene_sets', 'rbp_gene_sets', 'motif_database', 'clip_peaks', 'fimo_runtime'].findAll { params[it] }
    request.resource_names = names
    def resources = names.collect { file(params[it], checkIfExists: true) }
    [request, resources]
}

def requireParameters(names) {
    names.each { name -> if (!params[name]) error("Required parameter missing: --${name}") }
}

def validateStrandednessParameters() {
    def stranded
    def unstranded
    try {
        stranded = Double.parseDouble(params.stranded_threshold.toString())
        unstranded = Double.parseDouble(params.unstranded_threshold.toString())
    } catch (Exception ignored) {
        error('Strandedness thresholds must be finite numeric values')
    }
    if (!Double.isFinite(stranded) || stranded <= 0.5 || stranded > 1 ||
        !Double.isFinite(unstranded) || unstranded < 0 || unstranded >= 2 * stranded - 1) {
        error('Require 0.5 < stranded_threshold <= 1 and 0 <= unstranded_threshold < 2*stranded_threshold-1')
    }
    if (!(params.strandedness_sample_size.toString() ==~ /^[0-9]+$/) ||
        !(params.strandedness_mapq.toString() ==~ /^[0-9]+$/)) {
        error('strandedness_sample_size and strandedness_mapq must be integers')
    }
    if (new BigInteger(params.strandedness_sample_size.toString()) < 1 ||
        new BigInteger(params.strandedness_sample_size.toString()) > Integer.MAX_VALUE ||
        new BigInteger(params.strandedness_mapq.toString()) > 255) {
        error('Require positive strandedness_sample_size <= 2147483647 and 0 <= strandedness_mapq <= 255')
    }
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
                          --species 'Mus musculus' --gene_sets pathways.gmt --rbp_gene_sets rbps.gmt
                          --motif_database motifs.json --clip_peaks matched-peaks.tsv
                          --dexseq_tools_container <built-image-or-SIF>
  --mode upstream-only    --input samples.csv --contrasts contrasts.csv --fasta genome.fa --gtf annotation.gtf
  --mode downstream-only  --project project.json
  --mode downstream-only  --project project.json --flow_plan /absolute/original/plan.json -resume
  --outdir results        -profile local|slurm,docker|singularity|apptainer
Existing completed analysis: use downstream-only. Use -resume for this run's Nextflow cache.
Use -c examples/nextflow/kth.config for the supplied Dardel shared-queue example.
Optional analyses: --skip_network --skip_motiftrace --skip_motif_binding --skip_mapp
  --skip_enrichment --skip_rbp_identity --skip_landscape --skip_cross_comparison
  --skip_interaction --skip_sequencing_qc
Extra gene-body BAM scanning: --with_genebody_coverage (off by default).
Internal JSON is generated automatically. Advanced --project_template/--flow_config still work.
See docs/NEXTFLOW_SIMPLE.md and examples/nextflow/params.complete.yaml for all settings.
Samplesheet strandedness: forward|reverse|unstranded|auto. Auto uses matched genome-BAM RSeQC evidence.
Required calculations are verified before the fixed seven-section report is accepted.'''
    } else {
        if (!(params.mode in ['complete', 'upstream-only', 'downstream-only'])) error('Invalid --mode')
        if (!(params.storage_mode in ['retain', 'lean'])) error('Invalid --storage_mode; use retain or lean')
        if (params.flow_plan && params.mode != 'downstream-only') error('--flow_plan reuses completed downstream work and requires --mode downstream-only')
        if (!(params.katana_stage_containers instanceof Map)) error('--katana_stage_containers must be a map of execution groups to images')
        def groups = ['control', 'report', 'audit', 'quality', 'as', 'deg', 'dtu', 'landscape', 'evidence',
                      'enrichment', 'motiftrace', 'overlap', 'reference', 'qc', 'coverage', 'network_fit',
                      'network_followup', 'interaction', 'mapp_prepare', 'mapp']
        if (params.katana_stage_containers.any { key, value -> !(key in groups) || !(value instanceof String) || !value.trim() }) {
            error('Katana stage containers require known execution-group keys and nonempty image references')
        }
        if ((params.mode != 'upstream-only' || params.managed_preflight) && (workflow.containerEngine in ['docker', 'singularity', 'apptainer']) && !params.katana_container) {
            error('A Katana container profile requires --katana_container pointing to your built full runtime image; host fallback is disabled')
        }
        if ((params.mode != 'upstream-only' || params.managed_preflight) && (workflow.containerEngine in ['docker', 'singularity', 'apptainer']) && !(params.flow_executable in ['katana-flow', '/opt/conda/bin/katana-flow'])) {
            error('Container profiles require the image launcher katana-flow or /opt/conda/bin/katana-flow; mounted host launchers are not allowed')
        }
        if ((params.mode == 'complete' || params.managed_preflight) && params.skip_fastqc && !params.skip_sequencing_qc) {
            error('Complete fixed-report mode requires native sequencing FastQC. Use upstream-only for an explicitly skipped QC stage.')
        }
        if (params.mode == 'downstream-only') {
            requireParameters(['project'])
            project = file(params.project, checkIfExists: true)
            if (params.flow_plan) {
                if (params.flow_config) log.info 'Using the original sealed --flow_plan; --flow_config is not replanned or reapplied.'
                KATANA_DOWNSTREAM(Channel.value(inputDirectory(project.parent)), project.name,
                    Channel.empty(), '')
            } else {
                if (params.flow_config) {
                    config_file = file(params.flow_config, checkIfExists: true)
                    KATANA_DOWNSTREAM(Channel.value(inputDirectory(project.parent)), project.name,
                        Channel.value(inputDirectory(config_file.parent)), config_file.name)
                } else {
                    def automatic = automaticConfiguration()
                    KATANA_FLOW_CONFIGURE(automatic[0], automatic[1], file("${projectDir}/bin/katana_configure.py"))
                    KATANA_DOWNSTREAM(Channel.value(inputDirectory(project.parent)), project.name,
                        KATANA_FLOW_CONFIGURE.out.configuration.first(), 'flow.json')
                }
            }
        } else {
            validateStrandednessParameters()
            requireParameters(['input', 'contrasts', 'fasta', 'gtf'])
            if ((workflow.containerEngine in ['docker', 'singularity', 'apptainer']) && !params.dexseq_tools_container) {
                error('Raw FASTQ/BAM execution requires --dexseq_tools_container pointing to the built DEXSeq tools image or SIF. See containers/dexseq-tools/README.md; no unpublished registry image is assumed.')
            }
            samples = Channel.value(file(params.input, checkIfExists: true))
            contrasts = Channel.value(file(params.contrasts, checkIfExists: true))
            fasta = Channel.value(file(params.fasta, checkIfExists: true))
            gtf = Channel.value(file(params.gtf, checkIfExists: true))
            if (params.mode == 'complete' || params.managed_preflight) {
                runtime_resources = []
                if (params.project_template || params.flow_config) {
                    requireParameters(['project_template', 'flow_config'])
                    template_file = file(params.project_template, checkIfExists: true)
                    template = Channel.value(inputDirectory(template_file.parent))
                    config_file = file(params.flow_config, checkIfExists: true)
                    config = Channel.value(inputDirectory(config_file.parent))
                    template_name = template_file.name
                    config_name = config_file.name
                    log.info 'Using advanced project_template/flow_config; their downstream selections take precedence over flat downstream options.'
                } else {
                    def automatic = automaticConfiguration()
                    runtime_resources = automatic[1]
                    KATANA_FLOW_CONFIGURE(automatic[0], runtime_resources, file("${projectDir}/bin/katana_configure.py"))
                    template = KATANA_FLOW_CONFIGURE.out.configuration.first()
                    config = template
                    template_name = 'project-template.json'
                    config_name = 'flow.json'
                }
                KATANA_FLOW_PREFLIGHT(template, template_name, config, config_name, runtime_resources)
                samples = samples.combine(KATANA_FLOW_PREFLIGHT.out.receipt).map { sheet, ready -> sheet }.first()
                contrasts = contrasts.combine(KATANA_FLOW_PREFLIGHT.out.receipt).map { sheet, ready -> sheet }.first()
                fasta = fasta.combine(KATANA_FLOW_PREFLIGHT.out.receipt).map { ref, ready -> ref }.first()
                gtf = gtf.combine(KATANA_FLOW_PREFLIGHT.out.receipt).map { ref, ready -> ref }.first()
            }
            KATANA_PIPELINE(samples, contrasts, fasta, gtf)
            if (params.mode == 'complete') {
                KATANA_FLOW_IMPORT(KATANA_PIPELINE.out.bundle, template, template_name, config, config_name, runtime_resources)
                KATANA_DOWNSTREAM(KATANA_FLOW_IMPORT.out.workspace.first(), 'project.json', config, config_name)
            }
        }
    }
}
