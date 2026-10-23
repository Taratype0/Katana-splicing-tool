process ADAPT_SUPPA_CONTRASTS {
    label 'process_low'
    container 'docker.io/library/python:3.12.11-bookworm'
    input:
    path contrastsheet
    output:
    path 'contrastsheet.suppa.csv', emit: contrasts
    script:
    """
    python3 - <<'PY'
    import csv
    with open('${contrastsheet}', newline='') as src:
        reader = csv.DictReader(src)
        fields = reader.fieldnames
        rows = list(reader)
    # SUPPA2 reports PSI(second input) - PSI(first input). Place the declared
    # control first and treatment second; preserve the comparison identifier.
    for row in rows:
        row['treatment'], row['control'] = row['control'], row['treatment']
    with open('contrastsheet.suppa.csv', 'w', newline='') as dst:
        writer = csv.DictWriter(dst, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    PY
    """
}

