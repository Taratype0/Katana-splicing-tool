import groovy.json.JsonOutput

process KATANA_EXPORT {
    label 'process_low'
    container 'docker.io/library/python:3.12.11-bookworm'
    input:
    tuple val(artifact_metadata), path(artifact_files, stageAs: 'artifacts/??????/*')
    path samplesheet
    path contrastsheet
    path exporter
    val study_metadata
    output:
    path 'katana_bundle', emit: bundle
    path 'source_map.json', emit: source_map
    script:
    def paths = artifact_files instanceof List ? artifact_files : [artifact_files]
    if (paths.size() != artifact_metadata.size()) error('Exporter artifact paths and metadata differ in length')
    def entries = []
    for (int i = 0; i < paths.size(); i++) {
        entries << (artifact_metadata[i] + [source: paths[i].toString()])
    }
    def source_map = study_metadata + [schema_version: 1, artifacts: entries]
    // Encode the UTF-8 payload before embedding it in the shell script. A
    // multiline interpolated JSON string can prevent Nextflow's indentation
    // stripping and leave a heredoc terminator indented in .command.sh.
    def source_json_b64 = JsonOutput.prettyPrint(JsonOutput.toJson(source_map))
        .getBytes('UTF-8').encodeBase64().toString()
    """
    printf '%s' '${source_json_b64}' | python3 -c 'import base64, sys; from pathlib import Path; Path("source_map.json").write_bytes(base64.b64decode(sys.stdin.buffer.read(), validate=True))'
    python3 "$exporter" --source-map source_map.json --output-dir katana_bundle \\
      --samples "$samplesheet" --contrasts "$contrastsheet" --profile full
    """
}
