"""Local ORA over a complete, scoped result universe; never refit source tests."""
import json
from pathlib import Path
from uuid import uuid4

import pandas as pd

from katana_next.science.enrichment import gene_token, read_gene_sets, overrepresentation
from katana_next.services.report_import import sha
from katana_next.services.report_query import query_table, clauses_for, quote, json_rows


def enrich(context, artifact, filters, *, library_path, gene_column, test_column,
           library_species, minimum_term_size=1, maximum_term_size=2000):
    library = Path(library_path).expanduser().resolve()
    if not library.is_file() or library.suffix.lower() not in {'.gmt', '.tsv', '.csv'}:
        raise ValueError('Choose a local GMT or term/gene TSV/CSV pathway library')
    build = context.manifest.annotation.genome_build.lower()
    if (library_species == 'human' and build in {'mm10', 'mm39', 'grcm38', 'grcm39'} or
            library_species == 'mouse' and build in {'hg19', 'hg38', 'grch37', 'grch38'}):
        raise ValueError('Pathway library species does not match the project reference')
    metadata = artifact.metadata.get('source_metadata', artifact.metadata)
    scope = str(metadata.get('result_scope', metadata.get('source_result_scope', 'not_declared'))).lower()
    if (metadata.get('significant_only') or metadata.get('filtered_to_significant') or metadata.get('full_results') is False
            or scope in {'significant_only', 'hits_only', 'significant', 'selected', 'subset', 'top_n', 'foreground_only'}):
        raise ValueError('A significant-only/subset table cannot define the background; select the complete DEG/DAS result table')
    name = (artifact.name + ' ' + artifact.metadata.get('source_relative_path', '')).lower()
    if any(token in name for token in ('top10', 'top 10', 'top-10', '-significant-', '_significant_', 'significant genes', 'significant events')):
        raise ValueError('Select the complete result table, not a significant-only or Top 10 view')
    if artifact.analysis_id == 'imported.report' and sha(artifact.path) != artifact.metadata['imported_table_sha256']:
        raise ValueError('Imported evidence changed after import')
    with query_table(artifact) as (db, columns):
        gene = quote(gene_column, columns)
        if gene_column not in {'gene_id', 'gene_symbol'}:
            raise ValueError('Use gene_id or gene_symbol to map a pathway library')
        tested = f'TRY_CAST({quote(test_column, columns)} AS DOUBLE)'
        background_where, background_values = clauses_for(filters, columns, background=True)
        db.execute('CREATE TEMP TABLE scoped AS SELECT * FROM evidence' + background_where, background_values)
        for column in ('comparison_id', 'cell_type', 'celltype', 'cluster', 'cluster_id'):
            if column in columns and db.execute(f'SELECT count(DISTINCT {quote(column, columns)}) FROM scoped').fetchone()[0] > 1:
                raise ValueError(f'Select one {column} before enrichment; statistical populations must remain separate')
        background_frame = db.execute(f'SELECT {gene} AS gene FROM scoped WHERE isfinite({tested}) AND {tested} BETWEEN 0 AND 1').fetchdf()
        foreground_where, foreground_values = clauses_for(filters, columns)
        foreground = db.execute(f'SELECT {gene} AS gene FROM evidence' + foreground_where, foreground_values).fetchdf()
    tested_genes = {gene_token(g) for g in background_frame.gene} - {''}
    selected = {gene_token(g) for g in foreground.gene} - {''}
    if not selected:
        raise ValueError('No selected genes match this identifier column')
    if not tested_genes:
        raise ValueError('No testable background genes; choose the original P/FDR column from complete results')
    library_sha = sha(library)
    gene_sets = {term: {gene_token(g) for g in genes} - {''} for term, genes in read_gene_sets(library).items()}
    if sha(library) != library_sha:
        raise ValueError('Library changed while being read')
    annotated = set().union(*gene_sets.values()) if gene_sets else set()
    background = tested_genes & annotated
    mapped_selected = selected & background
    if not mapped_selected:
        raise ValueError('No selected, testable genes map to this library. Check gene IDs/symbols and species')
    result = overrepresentation(mapped_selected, background, gene_sets, minimum_overlap=1,
                                minimum_term_size=minimum_term_size, maximum_term_size=maximum_term_size)
    if result.empty:
        raise ValueError('No pathways meet the background term-size range')
    run_id = uuid4().hex
    destination = context.manifest.output_root / 'local-enrichment' / run_id
    destination.mkdir(parents=True, exist_ok=False)
    audit = {'schema': 'katana.local_report_enrichment/v1', 'run_id': run_id,
        'project_id': context.manifest.project_id, 'assay_mode': context.manifest.metadata.get('assay_mode', 'unknown'),
        'source_artifact_id': artifact.artifact_id, 'source_sha256': sha(artifact.path),
        'library': library.name, 'library_sha256': library_sha, 'library_species': library_species,
        'gene_column': gene_column, 'test_column': test_column, 'filters': filters.model_dump(),
        'selected_unique_genes': len(selected), 'tested_unique_genes': len(tested_genes),
        'annotated_background_genes': len(background), 'mapped_selected_genes': len(mapped_selected),
        'unmapped_or_untestable_selected': sorted(selected - mapped_selected),
        'minimum_term_size': minimum_term_size, 'maximum_term_size': maximum_term_size,
        'tested_terms': len(result), 'background_rule': 'all genes with valid original P/FDR in the same comparison and biological scope, intersected with this library; event rows deduplicated to genes',
        'completeness': 'all rows delivered in the selected table; undisclosed upstream truncation cannot be certified',
        'method': 'one-sided hypergeometric ORA; BH across all eligible terms including zero overlaps',
        'native_statistics_recomputed': False}
    result.to_csv(destination / 'enrichment.tsv', sep='\t', index=False)
    pd.DataFrame({'gene': sorted(background)}).to_csv(destination / 'background.tsv', sep='\t', index=False)
    pd.DataFrame({'gene': sorted(mapped_selected)}).to_csv(destination / 'foreground.tsv', sep='\t', index=False)
    (destination / 'audit.json').write_text(json.dumps(audit, indent=2), encoding='utf-8')
    return {'audit': audit, 'rows': json_rows(result), 'columns': list(result.columns)}
