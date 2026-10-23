process SALMON_GENOME_PROJECTION {
    tag "${meta.id}"
    label 'process_medium'
    // Official prebuilt Linux binary is hash checked, not installed at runtime.
    container 'docker.io/library/python:3.12.11-bookworm'
    input:
    tuple val(meta), path(bam)
    path fasta
    path gtf
    path salmon_archive
    val expected_sha256
    output:
    tuple val(meta), path("${meta.id}"), emit: results
    path 'versions.yml', emit: versions
    script:
    def libtype = meta.single_end ? (meta.strandedness == 'forward' ? 'SF' : meta.strandedness == 'reverse' ? 'SR' : 'U') : (meta.strandedness == 'forward' ? 'ISF' : meta.strandedness == 'reverse' ? 'ISR' : 'IU')
    def bias_args = [params.salmon_seq_bias ? '--seqBias' : '', params.salmon_gc_bias ? '--gcBias' : ''].join(' ')
    """
    echo '${expected_sha256}  ${salmon_archive}' | sha256sum --check --strict
    python3 - <<'PY'
    import pathlib, tarfile
    with tarfile.open('${salmon_archive}') as archive:
        for member in archive.getmembers():
            if pathlib.PurePosixPath(member.name).name == 'salmon' and member.isfile():
                pathlib.Path('salmon').write_bytes(archive.extractfile(member).read())
                pathlib.Path('salmon').chmod(0o755)
                break
        else:
            raise SystemExit('Pinned Salmon archive has no salmon executable')
    PY
    ./salmon --version | grep -E '2\\.8\\.0'
    ./salmon quant -a "$bam" --annotation "$gtf" --genome "$fasta" -l $libtype \\
      --geneMap "$gtf" -p ${task.cpus} $bias_args -o ${meta.id}
    test -s ${meta.id}/quant.sf
    printf 'salmon: %s\\n' "\$(./salmon --version)" > versions.yml
    """
}
