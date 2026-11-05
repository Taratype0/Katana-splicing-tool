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
                'event_set', 'rmats_mode', 'contrast_ids', 'resource_metadata', 'design_factors',
                'deg_fdr', 'deg_log2fc', 'das_fdr', 'das_dpsi', 'dtu_fdr',
                'network_min_module_size', 'network_merge_cut_height',
                'with_genebody_coverage', 'with_network_bootstrap']
    keys += ['landscape', 'cross_comparison', 'enrichment', 'rbp_identity', 'motiftrace',
             'motif_binding', 'mapp', 'network', 'interaction', 'sequencing_qc'].collect { "skip_${it}".toString() }
    def request = keys.collectEntries { key -> [(key): params[key]] }
    def kinds = ['gene_sets', 'rbp_gene_sets', 'motif_database', 'clip_peaks', 'fimo_runtime']
    def resources = []
    def names = []
    def entries = []
    kinds.findAll { params[it] }.each { name ->
        def supplied = params[name] instanceof List ? params[name] : params[name].toString().split(',').toList()
        def paths = supplied.collectMany { pattern ->
            def matched = file(pattern.toString().trim(), checkIfExists: true)
            matched instanceof List ? matched : [matched]
        }.sort { it.toString() }
        paths.eachWithIndex { path, index ->
            resources << path
            names << name
            entries << [asset_id: paths.size() == 1 ? name : "${name}_${index + 1}".toString(), metadata: [:]]
        }
    }
    if (params.bundle) request.mode = 'complete'
    request.resource_names = names
    request.resource_entries = entries
    request.resource_cache = file(params.resource_cache).toAbsolutePath().normalize().toString()
    request.resource_catalog = []
    def raw = params.mode != 'downstream-only' || params.bundle
    if (raw && params.resource_bundle != 'none') {
        def bundle = params.resource_bundle
        def mouse = params.species?.toString()?.replace(' ', '_')?.toLowerCase() in ['mus_musculus', 'mouse']
        if (bundle == 'auto') bundle = mouse && params.reference_build == 'GRCm39' ? 'mouse-grcm39-v1' : 'none'
        if (bundle != 'none') {
            if (bundle != 'mouse-grcm39-v1' || !mouse || params.reference_build != 'GRCm39')
                error('The mouse-grcm39-v1 resource bundle requires Mus_musculus and GRCm39. Use --resource_bundle none with your own databases for another species/build.')
            def catalog = new groovy.json.JsonSlurper().parse(file("${projectDir}/assets/resources/mouse-grcm39-v1.json").toFile())
            def wanted = []
            if (!params.skip_enrichment || !params.skip_network) wanted << 'gene_sets'
            if (!params.skip_rbp_identity) wanted << 'rbp_gene_sets'
            if (!params.skip_motiftrace || !params.skip_motif_binding) wanted << 'motif_database'
            if (!params.skip_motif_binding) wanted << 'clip_peaks'
            request.resource_catalog = catalog.assets.findAll { it.parameter in wanted && !(it.parameter in names) }
        }
    }

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
                          --resource_bundle auto  (mouse / GRCm39 defaults; optional local resource overrides)
  --mode upstream-only    --input samples.csv --contrasts contrasts.csv --fasta genome.fa --gtf annotation.gtf
  --mode downstream-only  --project project.json
  --mode downstream-only  --bundle /data/completed/katana_bundle
  --mode downstream-only  --project project.json --flow_plan /absolute/original/plan.json -resume
  --outdir results        -profile local|slurm,docker|singularity|apptainer
Existing completed analysis: use downstream-only. Use -resume for this run's Nextflow cache.
Use -c examples/nextflow/kth.config for the supplied Dardel shared-queue example.
Optional analyses: --skip_network --skip_motiftrace --skip_motif_binding --skip_mapp
  --skip_enrichment --skip_rbp_identity --skip_landscape --skip_cross_comparison
  --skip_interaction --skip_sequencing_qc
Extra gene-body BAM scanning: --with_genebody_coverage (off by default).
Internal JSON is generated automatically. Advanced --project_template/--flow_config still work.
Use nextflow -bg run to launch in the background. See docs/NEXTFLOW_SIMPLE.md for direct CLI examples.
Containers download automatically for Docker, Singularity and Apptainer. No hand-written JSON is required.
Samplesheet strandedness: forward|reverse|unstranded|auto. Auto uses matched genome-BAM RSeQC evidence.
Required calculations are verified before the fixed seven-section report is accepted.'''
    } else {
        if (!(params.mode in ['complete', 'upstream-only', 'downstream-only'])) error('Invalid --mode')
        if (!(params.storage_mode in ['retain', 'lean'])) error('Invalid --storage_mode; use retain or lean')
        if (!(params.output_layout in ['standard', 'legacy'])) error('Invalid --output_layout')
        if (params.clean_work && (params.output_layout != 'standard' || !params.native_delivery || params.publish_dir_mode != 'copy')) {
            error('--clean_work requires standard copied outputs and native delivery verification')
        }
        if (params.bundle && (params.mode != 'downstream-only' || params.project || params.flow_plan)) {
            error('--bundle selects completed upstream results; use downstream-only without --project or --flow_plan')
        }
        if (params.flow_plan && params.mode != 'downstream-only') error('--flow_plan reuses completed downstream work and requires --mode downstream-only')
        if (!(params.katana_stage_containers instanceof Map)) error('--katana_stage_containers must map execution groups or exact task IDs to images')
        // Exact IDs are checked against the real task plan before workers start.
        if (params.katana_stage_containers.any { key, value -> !(key instanceof String) ||
                !(key ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/) || !(value instanceof String) || !value.trim() }) {
            error('Katana container overrides require safe execution-group/task keys and nonempty image references')
        }
        if ((params.mode != 'upstream-only' || params.managed_preflight) && (workflow.containerEngine in ['docker', 'singularity', 'apptainer']) && !params.katana_container) {
            error('No Katana container is configured; host fallback is disabled. Restore the release defaults or supply --katana_container.')
        }
        if ((params.mode != 'upstream-only' || params.managed_preflight) && (workflow.containerEngine in ['docker', 'singularity', 'apptainer']) && !(params.flow_executable in ['katana-flow', '/opt/conda/bin/katana-flow'])) {
            error('Container profiles require the image launcher katana-flow or /opt/conda/bin/katana-flow; mounted host launchers are not allowed')
        }
        if ((params.mode == 'complete' || params.managed_preflight) && params.skip_fastqc && !params.skip_sequencing_qc) {
            error('Complete fixed-report mode requires native sequencing FastQC. Use upstream-only for an explicitly skipped QC stage.')
        }
        if (params.mode == 'downstream-only' && params.bundle) {
            completed_bundle = Channel.value(inputDirectory(file(params.bundle, checkIfExists: true)))
            runtime_resources = []
            if (params.project_template || params.flow_config) {
                requireParameters(['project_template', 'flow_config'])
                template_file = file(params.project_template, checkIfExists: true)
                config_file = file(params.flow_config, checkIfExists: true)
                template = Channel.value(inputDirectory(template_file.parent))
                config = Channel.value(inputDirectory(config_file.parent))
                template_name = template_file.name
                config_name = config_file.name
            } else {
                def automatic = automaticConfiguration()
                runtime_resources = automatic[1]
                KATANA_FLOW_CONFIGURE(automatic[0], runtime_resources, file("${projectDir}/bin/katana_configure.py"))
                template = KATANA_FLOW_CONFIGURE.out.configuration.first()
                config = template
                template_name = 'project-template.json'
                config_name = 'flow.json'
            }
            KATANA_FLOW_IMPORT(completed_bundle, template, template_name, config, config_name, runtime_resources)
            KATANA_DOWNSTREAM(KATANA_FLOW_IMPORT.out.workspace.first(), 'project.json', config, config_name)
        } else if (params.mode == 'downstream-only') {
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
                error('No DEXSeq container is configured; restore the release defaults or supply --dexseq_tools_container.')
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
