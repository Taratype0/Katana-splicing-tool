"""Resolve display labels from the project's own reference, without inference."""
from pathlib import Path


def report_gene_names(project, artifacts, destination):
    import pandas as pd
    from .sources import read_table
    tables = []
    for record in artifacts:
        if record.get('name') == 'Project GTF gene-name resolution':
            path = Path(record['path'])
            if path.is_file():
                frame = read_table(path)
                if {'gene_id', 'gene_symbol'} <= set(frame):
                    tables.append(frame[['gene_id', 'gene_symbol']])
    gtf = (project.get('annotation') or {}).get('gtf')
    if gtf and Path(gtf).is_file():
        from katana_next.science.isoform import transcript_gene_map
        cache = Path(destination) / 'outputs' / 'annotation' / 'reference-gene-names.parquet'
        cache.parent.mkdir(parents=True, exist_ok=True)
        tables.append(transcript_gene_map(Path(gtf), cache)[['gene_id', 'gene_symbol']])
    return (pd.concat(tables, ignore_index=True).drop_duplicates()
            if tables else pd.DataFrame(columns=['gene_id', 'gene_symbol']))


def annotate_report_genes(frame, annotation):
    from katana_next.science.enrichment import fill_gene_symbols
    result = fill_gene_symbols(frame, annotation)
    if 'gene_symbol' in frame:
        result['source_gene_symbol'] = frame['gene_symbol']
    return result
