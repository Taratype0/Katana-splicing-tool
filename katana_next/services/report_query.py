"""One full-population query contract for tables, plots, targets and ORA."""
from contextlib import contextmanager
import json

import duckdb

SCOPES = ('cell_type', 'celltype', 'cluster', 'cluster_id', 'donor_id', 'event_type', 'condition', 'passes_read_support', 'passes_precision_support', 'passes_legacy_all10', 'nmd_asymmetry_rescued', 'nmd_switch_class', 'frame_status', 'splicenet_module')


def result_table(context, artifact_id):
    artifact = context.catalog.get_artifact(artifact_id)
    if (not artifact or artifact.project_id != context.manifest.project_id
            or artifact.kind != 'table' or artifact.path.suffix != '.parquet'):
        raise ValueError('A complete local result table is required')
    if artifact.analysis_id != 'imported.report':
        from katana_next.workflows.sources import validate_source
        validate_source(context, artifact)
    return artifact


@contextmanager
def query_table(artifact):
    with duckdb.connect() as db:
        db.read_parquet(str(artifact.path)).create_view('original_evidence')
        original_columns = [r[0] for r in db.execute('DESCRIBE original_evidence').fetchall()]
        if 'gene_id' in original_columns and 'gene_name' not in original_columns:
            from katana_next.services.report_names import display_names
            names = display_names(artifact)
            db.register('display_names', names)
            symbol = 'CAST(gene_symbol AS VARCHAR)' if 'gene_symbol' in original_columns else "''"
            native = f"({symbol} IS NOT NULL AND trim({symbol}) NOT IN ('', 'NA', 'nan') AND NOT starts_with({symbol}, 'ENS'))"
            db.execute(f"""CREATE TEMP VIEW evidence AS SELECT original_evidence.*,
                CASE WHEN {native} THEN {symbol} ELSE _mapped_name END AS gene_name,
                CASE WHEN {native} THEN 'source gene_symbol' WHEN _mapped_name IS NOT NULL THEN 'report identifier mapping' ELSE 'unresolved' END AS gene_name_source
                FROM original_evidence LEFT JOIN display_names ON regexp_replace(CAST(gene_id AS VARCHAR), '\\.[0-9]+$', '') = _name_id""")
        else:
            db.execute('CREATE TEMP VIEW evidence AS SELECT * FROM original_evidence')
        columns = [r[0] for r in db.execute('DESCRIBE evidence').fetchall()]
        yield db, columns


def quote(column, columns):
    if column not in columns:
        raise ValueError(f'Unknown column: {column}')
    return '"' + column.replace('"', '""') + '"'


def clauses_for(filters, columns, *, background=False):
    clauses, values = [], []
    def numeric(name):
        return f'TRY_CAST({quote(name, columns)} AS DOUBLE)'
    if filters.comparison:
        clauses.append(quote('comparison_id', columns) + ' = ?')
        values.append(filters.comparison)
    for column, value in filters.scopes.items():
        if column not in SCOPES:
            raise ValueError('Unsupported biological scope column')
        if value:
            clauses.append(f'CAST({quote(column, columns)} AS VARCHAR) = ?')
            values.append(value)
    if not background:
        if filters.identifiers:
            names = [c for c in ('gene_id', 'gene_symbol', 'gene_name', 'event_id') if c in columns]
            if not names:
                raise ValueError('This table has no gene/event identifiers')
            clauses.append('(' + ' OR '.join(f'CAST({quote(c, columns)} AS VARCHAR) IN (SELECT unnest(?))' for c in names) + ')')
            values.extend([filters.identifiers] * len(names))
        if filters.search:
            clauses.append('(' + ' OR '.join(f'contains(lower(CAST({quote(c, columns)} AS VARCHAR)), lower(?))' for c in columns) + ')')
            values.extend([filters.search] * len(columns))
        if filters.row_ids is not None:
            clauses.append(quote('row_id', columns) + ' IN (SELECT unnest(?))')
            values.append(filters.row_ids)
        if filters.fdr_max is not None or filters.fdr_min is not None:
            fdr = numeric(filters.fdr_column)
            clauses.append(f'isfinite({fdr}) AND {fdr} BETWEEN 0 AND 1')
            for bound, operator in ((filters.fdr_max, '<='), (filters.fdr_min, '>=')):
                if bound is not None:
                    clauses.append(f'{fdr} {operator} ?'); values.append(bound)
        effect = filters.effect_column
        if filters.effect_min is not None or filters.effect_low is not None or filters.effect_high is not None or filters.direction != 'both':
            num = numeric(effect)
            clauses.append(f'isfinite({num})')
            for bound, operator, expr in ((filters.effect_min, '>=', f'abs({num})'),
                                           (filters.effect_low, '>=', num), (filters.effect_high, '<=', num)):
                if bound is not None:
                    clauses.append(f'{expr} {operator} ?'); values.append(bound)
            if filters.direction != 'both':
                clauses.append(f'{num} {">" if filters.direction == "positive" else "<"} 0')
    return (' WHERE ' + ' AND '.join(clauses) if clauses else ''), values


def json_rows(frame):
    # pandas preserves nonfinite/missing as JSON null; source Parquet is unchanged.
    return json.loads(frame.to_json(orient='records', double_precision=15))


def browse(artifact, filters):
    with query_table(artifact) as (db, columns):
        where, values = clauses_for(filters, columns)
        total = db.execute('SELECT count(*) FROM evidence' + where, values).fetchone()[0]
        order = quote(filters.sort_column or 'row_id', columns)
        if filters.sort_column in {filters.fdr_column, filters.effect_column}:
            order = f'TRY_CAST({order} AS DOUBLE)'
        order += ' DESC' if filters.sort_desc else ' ASC'
        frame = db.execute('SELECT * FROM evidence' + where + f' ORDER BY {order} NULLS LAST, row_id LIMIT ? OFFSET ?',
                           values + [filters.limit, filters.offset]).fetchdf()
        scopes = {}
        for column in ('comparison_id', *SCOPES):
            if column in columns:
                scopes[column] = [r[0] for r in db.execute(f'SELECT DISTINCT CAST({quote(column, columns)} AS VARCHAR) FROM evidence WHERE {quote(column, columns)} IS NOT NULL ORDER BY 1').fetchall()]
        return {'columns': columns, 'rows': json_rows(frame), 'total': total,
                'source_rows': artifact.row_count, 'offset': filters.offset, 'limit': filters.limit, 'scopes': scopes}


def select_rows(artifact, filters):
    with query_table(artifact) as (db, columns):
        where, values = clauses_for(filters, columns)
        names = [c for c in ('row_id', 'gene_id', 'gene_symbol', 'gene_name', 'event_id', 'event_type', 'comparison_id', 'cell_type') if c in columns]
        frame = db.execute('SELECT ' + ','.join(quote(c, columns) for c in names) + ' FROM evidence' + where + ' ORDER BY row_id', values).fetchdf()
        return {'rows': json_rows(frame), 'total': len(frame), 'scope': 'all matching rows; no page limit'}


def plot(artifact, filters):
    with query_table(artifact) as (db, columns):
        where, values = clauses_for(filters, columns)
        x, p = [f'TRY_CAST({quote(c, columns)} AS DOUBLE)' for c in (filters.effect_column, filters.fdr_column)]
        label = "coalesce(" + ','.join(f"NULLIF(CAST({quote(c, columns)} AS VARCHAR), '')" for c in ('gene_name', 'gene_symbol', 'gene_id', 'event_id', 'row_id') if c in columns) + ")"
        support = 'coalesce(TRY_CAST(passes_read_support AS BOOLEAN), false)' if 'passes_read_support' in columns else 'true'
        db.execute(f'CREATE TEMP TABLE plotted AS SELECT row_id, {label} AS label, {x} AS x, {p} AS p, {support} AS read_supported FROM evidence' + where, values)
        total = db.execute('SELECT count(*) FROM plotted').fetchone()[0]
        db.execute('CREATE TEMP VIEW valid AS SELECT *, -log10(greatest(p, 1e-300)) AS y FROM plotted WHERE isfinite(x) AND isfinite(p) AND p BETWEEN 0 AND 1')
        count, lo, hi, ymax = db.execute('SELECT count(*), min(x), max(x), max(y) FROM valid').fetchone()
        extent = max(1., abs(lo or 0.), abs(hi or 0.)) * 1.05
        lo, hi = -extent, extent
        ymax = max(1., ymax or 0.)
        result = {'total': total, 'plotted': count, 'unplottable': total-count, 'x_min': lo, 'x_max': hi,
                  'y_max': ymax, 'zero_fdr_display_floor': 1e-300, 'scope': 'all matching rows'}
        if count <= 150000:
            result.update(mode='points', points=db.execute('SELECT x, y, CAST(row_id AS VARCHAR), CAST(label AS VARCHAR), read_supported, p FROM valid ORDER BY row_id').fetchall())
        else:
            # Full-population density, never sampling or truncating to the current page.
            bins = db.execute('SELECT least(239, floor((x-?)/?*240))::INTEGER AS bucket_x, least(119, floor(y/?*120))::INTEGER AS bucket_y, count(*) FROM valid GROUP BY bucket_x, bucket_y ORDER BY bucket_x, bucket_y', [lo, hi-lo, ymax]).fetchall()
            result.update(mode='density', bins=bins, grid=[240, 120])
        return result
