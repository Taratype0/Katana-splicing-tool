process SPLIT_CONTRASTS {
    label 'process_single'
    container 'docker.io/library/python:3.12.11-bookworm'

    input:
    path contrastsheet
    path splitter

    output:
    path 'contrasts/contrast.*.csv', emit: sheets
    path 'versions.yml', emit: versions

    script:
    """
    python3 '${splitter}' '${contrastsheet}' contrasts
    printf 'SPLIT_CONTRASTS:\\n  python: %s\\n' "\$(python3 --version | cut -d' ' -f2)" > versions.yml
    """
}

workflow CONTRAST_TASK_SHEETS {
    take:
    contrastsheet

    main:
    SPLIT_CONTRASTS(contrastsheet, Channel.value(file("${projectDir}/bin/kanata_split_contrasts.py", checkIfExists: true)))
    ch_sheets = SPLIT_CONTRASTS.out.sheets.flatten().map { sheet ->
        tuple(sheet.baseName.toString().replaceFirst(/^contrast\./, ''), sheet)
    }

    emit:
    sheets = ch_sheets
    versions = SPLIT_CONTRASTS.out.versions
}
