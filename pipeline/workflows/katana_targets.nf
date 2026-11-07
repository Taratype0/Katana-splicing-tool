// Independent selected-target DAG. No FASTQ, alignment or statistical refits.
process TARGET_REQUEST_INSPECT {
    tag 'verify target request'
    cpus 1
    container { params.katana_container ?: null }
    input:
    path request
    output:
    path 'targets.tsv'
    script:
    """
    python -m katana_next.targets tasks --request '${request}' --format tsv > targets.tsv
    """
}

process TARGET_DETAIL {
    tag { target_id }
    cpus 1
    container { params.katana_container ?: null }
    publishDir "${params.outdir}/targets", mode: 'copy'
    input:
    val target_id
    path request
    val resource_ids
    path resources, stageAs: 'resources/resource??/*'
    output:
    path "${target_id}"
    script:
    def staged = resources instanceof List ? resources : [resources]
    def bindings = groovy.json.JsonOutput.toJson([schema:'katana.target_bindings/v1', resources:
        [resource_ids, staged.collect { it.toString() }].transpose().collectEntries { [(it[0]): it[1]] }])
    def encoded = bindings.bytes.encodeBase64().toString()
    """
    python - <<'PY'
    import base64
    from pathlib import Path
    Path('bindings.json').write_bytes(base64.b64decode('${encoded}'))
    PY
    set +e
    python -m katana_next.targets run --request '${request}' --bindings bindings.json --outdir '${target_id}' --target-id '${target_id}' --workers ${task.cpus} > execution.json
    code=\$?
    set -e
    if [ "\$code" -ne 0 ] && [ "\$code" -ne 2 ]; then cat execution.json; exit "\$code"; fi
    cp execution.json '${target_id}/pipeline-execution.json'
    """
}

workflow KATANA_TARGETS {
    take:
    request
    resource_ids
    resources
    main:
    TARGET_REQUEST_INSPECT(request)
    ids = TARGET_REQUEST_INSPECT.out.splitCsv(header:true, sep:'\t').map { row ->
        if (!(row.target_id ==~ /(?:gene|event)-[a-f0-9]+/)) error('Invalid exact target ID')
        row.target_id
    }
    TARGET_DETAIL(ids, request, resource_ids, resources)
    emit:
    reports = TARGET_DETAIL.out
}
