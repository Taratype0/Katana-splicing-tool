import groovy.json.JsonOutput

process RESOLVE_SAMPLE_STRANDEDNESS {
    label 'process_single'
    container 'docker.io/library/python:3.12.11-bookworm'
    input:
    path samplesheet
    path contrastsheet
    tuple val(sample_ids), path(reports, stageAs: 'strand_reports/*'), path(logs, stageAs: 'strand_logs/*')
    path resolver
    path validator
    val stranded_threshold
    val unstranded_threshold
    val design
    val analysis_purpose
    output:
    path 'strandedness/samplesheet.resolved.csv', emit: samplesheet
    path 'strandedness/strandedness.audit.json', emit: audit
    path 'strandedness/strandedness.summary.tsv', emit: summary
    script:
    def evidence = [:]
    def report_list = reports instanceof List ? reports : [reports]
    def log_list = logs instanceof List ? logs : [logs]
    if (sample_ids.size() != report_list.size() || sample_ids.size() != log_list.size()) {
        error('Strandedness evidence IDs and native files have different lengths')
    }
    sample_ids.eachWithIndex { id, i ->
        if (evidence.containsKey(id)) error("Duplicate native strandedness evidence for ${id}")
        evidence[id] = [report: report_list[i].toString(), log: log_list[i].toString()]
    }
    def payload = JsonOutput.toJson(evidence).getBytes('UTF-8').encodeBase64().toString()
    def quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    """
    printf '%s' '${payload}' | python3 -c 'import base64,sys; from pathlib import Path; Path("evidence.json").write_bytes(base64.b64decode(sys.stdin.buffer.read(),validate=True))'
    python3 ${quote(resolver)} --samplesheet ${quote(samplesheet)} --contrasts ${quote(contrastsheet)} \\
      --evidence-map evidence.json --outdir strandedness \\
      --stranded-threshold ${stranded_threshold} --unstranded-threshold ${unstranded_threshold} \\
      --design ${quote(design)} --analysis-purpose ${quote(analysis_purpose)}
    """
}
