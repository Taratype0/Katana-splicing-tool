"""Read project manifests and artifact exports without opening their database."""
from pathlib import Path
import copy
import hashlib, json


def digest(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def _absolute_path(value, base):
    path = Path(value).expanduser()
    return (path if path.is_absolute() else base / path).resolve()


def _is_remote_uri(value):
    return '://' in str(value)


def _native_paths(raw, base):
    """Resolve local manifest path fields on a copy, never arbitrary labels."""
    native = copy.deepcopy(raw)

    def normalize(record, keys):
        for key in keys:
            if record.get(key) and not _is_remote_uri(record[key]):
                record[key] = str(_absolute_path(record[key], base))

    normalize(native, ('input_root', 'output_root'))
    normalize(native.get('annotation') or {}, ('gtf', 'genome_fasta', 'transcript_fasta'))
    for asset in native.get('assets', []):
        normalize(asset, ('path',))
    for sample in native.get('samples', []):
        normalize(sample.get('assets') or {}, ('quant_sf',))
    return native


def _view_paths(raw, base):
    """Resolve only declared view file references, without changing metadata."""
    views = copy.deepcopy(raw)
    for section, keys in (('panels', ('path', 'image', 'svg')),
                          ('tables', ('path', 'source', 'view')),
                          ('top10', ('full_path', 'full_table'))):
        for record in views.get(section, []):
            for key in keys:
                if record.get(key):
                    record[key] = str(_absolute_path(record[key], base))
    return views


def load_source(project, catalog=None):
    """Load a descriptor using its own directory, or cwd for inline input.

    An explicit catalog argument follows command-line cwd semantics. Paths in
    a descriptor follow that descriptor, while native artifact paths follow the
    native project manifest. A file of extra views owns its own relative paths.
    """
    return _load_source(project, catalog, inline_base=Path.cwd())


def _load_source(project, catalog=None, *, inline_base):
    source_files = []
    if isinstance(project, dict):
        value = copy.deepcopy(project)
        descriptor_base = inline_base
    else:
        location = Path(project).expanduser().resolve()
        value = json.loads(location.read_text('utf-8-sig'))
        descriptor_base = location.parent
        source_files.append(location)
    native_base = descriptor_base
    if value.get('project_manifest'):
        manifest = _absolute_path(value['project_manifest'], descriptor_base)
        native = json.loads(manifest.read_text('utf-8-sig'))
        source_files.append(manifest)
        native_base = manifest.parent
    else:
        native = value.get('project', value)
    # Identity belongs to the exact source JSON, before paths are normalized in
    # a report-only copy. Existing supplementary evidence binds that identity.
    native_location = manifest if value.get('project_manifest') else (location if not isinstance(project, dict) and 'project' not in value else None)
    manifest_identity = digest(native_location) if native_location else hashlib.sha256(json.dumps(native, sort_keys=True, separators=(',', ':'), default=str).encode()).hexdigest()
    native = _native_paths(native, native_base)
    records = value.get('artifacts')
    explicit_catalog = catalog is not None
    catalog_base = Path.cwd() if catalog is not None else descriptor_base
    catalog = catalog if catalog is not None else (value.get('artifact_catalog') or value.get('catalog_export') or value.get('native_artifact_catalog'))
    if isinstance(records, str):
        if not explicit_catalog:
            catalog, catalog_base = records, descriptor_base
        records = None
    output_value = native.get('output_root', native_base)
    output_root = None if _is_remote_uri(output_value) else _absolute_path(output_value, native_base)
    if records is None:
        candidates = [_absolute_path(catalog, catalog_base)] if catalog else [native_base / 'native-artifact-catalog.json']
        if not catalog and output_root is not None:
            candidates.append(output_root / '.katana/artifact_manifest.json')
        for path in candidates:
            if path.is_file():
                records = json.loads(path.read_text('utf-8-sig'))
                records = records.get('artifacts', records) if isinstance(records, dict) else records
                source_files.append(path)
                break
    if records is None:
        if output_root is None:
            raise ValueError('Remote result URIs are not read by the report loader; supply a local completed artifact export.')
        db = output_root / '.katana/catalog.duckdb'
        if not db.is_file():
            raise ValueError('No completed artifact export found; supply catalog= with a JSON artifact export.')
        import duckdb
        with duckdb.connect(str(db), read_only=True) as connection:
            rows = connection.execute('SELECT record_json FROM artifacts WHERE project_id = ? ORDER BY created_at DESC', [native['project_id']]).fetchall()
            records = [json.loads(row[0]) for row in rows]
        source_files.append(db)
    if not isinstance(records, list):
        raise ValueError('Artifact export must contain a list of completed artifact records.')
    artifacts, seen = [], set()
    for raw in records:
        record = copy.deepcopy(raw)
        if record.get('status', 'complete') not in ('complete', 'completed'):
            continue
        if record.get('project_id') and record['project_id'] != native.get('project_id'):
            raise ValueError('Artifact belongs to a different project; use an explicitly verified supplemental source.')
        path = _absolute_path(record.get('path', ''), native_base)
        if not path.is_file():
            record['missing_reason'] = 'Recorded completed artifact file is unavailable.'
        record['path'] = str(path.resolve())
        if path.is_file():
            from .rich_results import _completion
            complete, evidence = _completion(record, path)
            if not complete:
                continue
            record.setdefault('metadata', {})['report_completion_evidence'] = evidence
        key = (record.get('analysis_id'), record.get('name'), tuple(record.get('contrast_ids', [])), record.get('parameters', {}).get('sequence_scope'), record.get('parameters', {}).get('event_set'))
        if key in seen:
            continue
        seen.add(key)
        artifacts.append(record)
    extras = copy.deepcopy(value.get('extra_assets', []))
    views = []
    for asset in extras:
        meta = asset.get('metadata', {})
        if asset.get('project_id') != native.get('project_id') or meta.get('source_project_manifest_sha256', meta.get('source_project_manifest_SHA256')) != manifest_identity:
            raise ValueError('Extra report assets require matching project identity and parent manifest hash.')
        if asset.get('status', 'complete') not in ('complete', 'completed'):
            raise ValueError('Extra report assets must be explicitly completed.')
        if asset.get('path'):
            asset['path'] = str(_absolute_path(asset['path'], descriptor_base))
        for key, path in meta.get('exports', {}).items():
            if key in ('png', 'svg') and path:
                meta['exports'][key] = str(_absolute_path(path, descriptor_base))
    for raw in value.get('extra_views', []):
        view_base = descriptor_base
        if isinstance(raw, str):
            path = _absolute_path(raw, descriptor_base)
            raw = json.loads(path.read_text('utf-8-sig'))
            source_files.append(path)
            view_base = path.parent
        if raw.get('project_id') != native.get('project_id') or raw.get('source_project_manifest_sha256') != manifest_identity:
            raise ValueError('Extra report views require matching project identity and parent manifest hash.')
        if raw.get('status', 'complete') not in ('complete', 'completed'):
            raise ValueError('Extra report views must be completed; running or failed collections cannot be presented as completed results.')
        views.append(_view_paths(raw, view_base))
    supplemental = []
    for item in value.get('supplemental_sources', []):
        if isinstance(item, str):
            item = str(_absolute_path(item, descriptor_base))
        child = _load_source(item, inline_base=descriptor_base)
        parent = child['project'].get('metadata', {})
        parent_id = parent.get('source_project_id')
        expected = parent.get('source_project_manifest_SHA256')
        actual = manifest_identity
        same_samples = [(s['sample_id'], s.get('biological_sample_id'), s.get('condition')) for s in native.get('samples', [])] == [(s['sample_id'], s.get('biological_sample_id'), s.get('condition')) for s in child['project'].get('samples', [])]
        same_contrasts = native.get('contrasts') == child['project'].get('contrasts')
        if not (same_samples and same_contrasts and expected and expected == actual):
            raise ValueError('Supplemental results require the exact parent manifest hash, samples and comparisons.')
        if parent_id and parent_id != native.get('project_id'):
            raise ValueError('Supplemental result parent identity does not match the source study.')
        supplemental.append(child)
    synthetic = bool(value.get('dummy', value.get('synthetic', native.get('metadata', {}).get('dummy', native.get('metadata', {}).get('synthetic', False)))))
    return {'project': native, 'artifacts': artifacts, 'supplemental': supplemental, 'extra_assets': extras, 'extra_views': views, 'dummy': synthetic,
            'search_roots': value.get('search_roots', []), 'source_files': source_files,
            'mode': value.get('mode', 'single-cell' if str(native.get('metadata', {}).get('assay', '')).startswith(('sc', 'sn', 'single')) else 'bulk')}


def read_table(path, *, columns=None, filters=None):
    path = Path(path)
    import pandas as pd
    if path.suffix == '.parquet':
        return pd.read_parquet(path, columns=columns, filters=filters)
    return pd.read_csv(path, sep='\t' if any(s.lower() in ('.tsv', '.txt') for s in path.suffixes) else ',', usecols=columns)


def table_schema(path):
    path = Path(path)
    if path.suffix == '.parquet':
        import pyarrow.parquet as pq
        meta = pq.ParquetFile(path)
        return {'columns': meta.schema_arrow.names, 'rows': meta.metadata.num_rows}
    import gzip
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8-sig') as handle:
        header = handle.readline().rstrip('\r\n').split('\t' if any(s.lower() in ('.tsv', '.txt') for s in path.suffixes) else ',')
        count = sum(1 for _ in handle)
    return {'columns': header, 'rows': count}
